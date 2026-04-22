"""Main annotation script for TEProf2 v3.

Annotates TACO-merged GTF transcripts with TE and gene information.
"""
from __future__ import annotations

import re
import sys
import pickle
from pathlib import Path
from typing import Optional

import pandas as pd
import pysam
import typer

from teprof3.annotation.gene_annotator import GeneAnnotator
from teprof3.annotation.te_query import query_te_at_tss
from teprof3.cli.output import split_outputs

app = typer.Typer(add_completion=False)

# ---------------------------------------------------------------------------
# Column schema
# ---------------------------------------------------------------------------

# TE classes where TSS is expected inside the TE body (LTR/LINE have internal Pol II promoters)
_TSS_INSIDE_TE_CLASSES = frozenset([
    "LTR", "LINE",
    "Unknown__ClassI_LTR_Gypsy", "Unknown__ClassI_LTR_Copia", "Unknown__ClassI_LTR",
    "Unknown__ClassI_nLTR_LINE_L1", "Unknown__ClassI_nLTR_LINE_RTE",
    "Unknown__ClassI_nLTR_LINE_Jockey", "Unknown__ClassI_nLTR_LINE_CR1",
    "Unknown__ClassI_nLTR_LINE_I", "Unknown__ClassI_nLTR_LINE_R1",
    "Unknown__ClassI_nLTR_LINE_R2", "Unknown__ClassI_nLTR_LINE",
    "Unknown__ClassI_nLTR_PLE",
])


def _te_driven_mode(strict_hit: bool, strict_te_class: str,
                    window_hit: bool, window_te_class: str) -> str:
    """Classify TE-driven mode based on hit type and TE class.

    - 'strict': TSS inside LTR/LINE body (TE provides the promoter directly)
    - 'window': MITE/DNA transposon upstream of TSS (TE acts as upstream regulator)
    - 'none': no TE hit in either mode
    """
    if strict_hit and strict_te_class in _TSS_INSIDE_TE_CLASSES:
        return "strict"
    if window_hit:
        return "window"
    if strict_hit:
        # strict hit but not LTR/LINE — treat as window-like (MITE inside TSS region)
        return "window"
    return "none"


EXPECTED_COLUMNS = [
    # base
    "transcript_id", "chrom", "start", "end", "strand",
    "n_exons", "exon_coords", "splicing",
    # tss gene annotation (12)
    "tss_gene_type", "tss_gene_name", "tss_chr", "tss_start", "tss_end",
    "tss_element", "tss_element_num", "tss_start_codon",
    "tss_transcript_start", "tss_transcript_end",
    "tss_element_list", "tss_ref_transcript_id",
    # tend gene annotation (12)
    "tend_gene_type", "tend_gene_name", "tend_chr", "tend_start", "tend_end",
    "tend_element", "tend_element_num", "tend_start_codon",
    "tend_transcript_start", "tend_transcript_end",
    "tend_element_list", "tend_ref_transcript_id",
    # strict TE (11)
    "strict_hit", "strict_te_name", "strict_te_class", "strict_te_family",
    "strict_te_chrom", "strict_te_start", "strict_te_end",
    "strict_te_strand", "strict_te_score",
    "strict_all_te_names", "strict_n_te_overlaps",
    # window TE (5)
    "window_hit", "window_te_name", "window_te_class",
    "window_te_family", "window_n_te_overlaps",
    # driven mode
    "te_driven_mode",
]

_TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------

def parse_gtf(source) -> pd.DataFrame:
    """Parse TACO-merged GTF and return one row per transcript.

    Args:
        source: file path (str/Path) or file-like object (for testing)

    Returns:
        DataFrame with columns: transcript_id, chrom, start, end, strand,
        tss, tend, n_exons, splicing, exon_coords
    """
    col_names = [
        "chrom", "source", "feature", "start", "end",
        "score", "strand", "frame", "attributes",
    ]
    raw = pd.read_csv(
        source,
        sep="\t",
        comment="#",
        header=None,
        names=col_names,
        dtype={"start": int, "end": int},
    )

    def _extract_tid(attr: str) -> str:
        m = _TRANSCRIPT_ID_RE.search(attr)
        return m.group(1) if m else ""

    raw["transcript_id"] = raw["attributes"].apply(_extract_tid)

    transcripts = raw[raw["feature"] == "transcript"].copy()
    exons = raw[raw["feature"] == "exon"].copy()

    # Aggregate exon info per transcript
    def _agg_exons(grp: pd.DataFrame) -> pd.Series:
        n = len(grp)
        # Sort by start for consistent ordering
        grp_sorted = grp.sort_values("start")
        strand = grp_sorted.iloc[0]["strand"]
        chrom = grp_sorted.iloc[0]["chrom"]
        # Build 0-based exon coords string (only first 2 exons needed per spec)
        parts = [strand, chrom]
        for _, row in grp_sorted.head(2).iterrows():
            parts.append(str(row["start"] - 1))  # convert to 0-based
            parts.append(str(row["end"]))
        coords = ",".join(parts)
        return pd.Series({"n_exons": n, "exon_coords": coords})

    exon_agg = exons.groupby("transcript_id").apply(_agg_exons).reset_index()

    # Merge exon info onto transcripts
    merged = transcripts.merge(exon_agg, on="transcript_id", how="left")
    merged["n_exons"] = merged["n_exons"].fillna(0).astype(int)
    merged["exon_coords"] = merged["exon_coords"].fillna("")

    # Compute TSS and tend (0-based)
    plus = merged["strand"] == "+"
    merged["tss"] = 0
    merged["tend"] = 0
    merged.loc[plus, "tss"] = merged.loc[plus, "start"] - 1
    merged.loc[~plus, "tss"] = merged.loc[~plus, "end"]
    merged.loc[plus, "tend"] = merged.loc[plus, "end"]
    merged.loc[~plus, "tend"] = merged.loc[~plus, "start"] - 1

    # Splicing
    merged["splicing"] = merged["n_exons"].apply(lambda n: "Yes" if n >= 2 else "No")

    return merged[
        ["transcript_id", "chrom", "start", "end", "strand",
         "tss", "tend", "n_exons", "splicing", "exon_coords"]
    ].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Tabix rmsk handler
# ---------------------------------------------------------------------------

class _TabixHandler:
    """Wraps pysam.TabixFile to provide .query(chrom, start, end)."""

    def __init__(self, tbx: pysam.TabixFile) -> None:
        self._tbx = tbx

    def query(self, chrom: str, start: int, end: int) -> list:
        records = []
        try:
            for row in self._tbx.fetch(chrom, max(0, start), end):
                parts = row.split("\t")
                records.append((
                    parts[0],
                    int(parts[1]),
                    int(parts[2]),
                    parts[3],
                    parts[4],
                    parts[5] if len(parts) > 5 else ".",
                ))
        except ValueError:
            pass  # contig not in index
        return records


# ---------------------------------------------------------------------------
# Core annotation loop
# ---------------------------------------------------------------------------

def build_records(
    df: pd.DataFrame,
    gene_annotator: GeneAnnotator,
    rmsk_handler,
    te_classification: dict,
    strict_bp: int,
    window_bp: int,
) -> list[dict]:
    """Annotate each transcript row and return list of record dicts.

    Each record contains all EXPECTED_COLUMNS plus _strict_hit and _window_hit.
    """
    records = []
    for _, row in df.iterrows():
        tid = row["transcript_id"]
        chrom = row["chrom"]
        strand = row["strand"]
        tss = int(row["tss"])
        tend = int(row["tend"])

        # TE annotation
        te_result = query_te_at_tss(
            chrom, tss, strand, strict_bp, window_bp, rmsk_handler, te_classification
        )
        s = te_result["strict"]
        w = te_result["window"]

        # Gene annotation
        tss_ann = gene_annotator.annotate_tss(chrom, tss, strand)
        tend_ann = gene_annotator.annotate_tend(chrom, tend, strand)

        rec = {
            "transcript_id": tid,
            "chrom": chrom,
            "start": int(row["start"]),
            "end": int(row["end"]),
            "strand": strand,
            "n_exons": int(row["n_exons"]),
            "exon_coords": row["exon_coords"],
            "splicing": row["splicing"],
            # tss gene annotation
            **tss_ann,
            # tend gene annotation
            **tend_ann,
            # strict TE
            "strict_hit": s["hit"],
            "strict_te_name": s["primary_te_name"],
            "strict_te_class": s["te_class"],
            "strict_te_family": s["te_family"],
            "strict_te_chrom": s["te_chrom"],
            "strict_te_start": s["te_start"],
            "strict_te_end": s["te_end"],
            "strict_te_strand": s["te_strand"],
            "strict_te_score": s["te_score"],
            "strict_all_te_names": s["all_te_names"],
            "strict_n_te_overlaps": s["n_te_overlaps"],
            # window TE
            "window_hit": w["hit"],
            "window_te_name": w["primary_te_name"],
            "window_te_class": w["te_class"],
            "window_te_family": w["te_family"],
            "window_n_te_overlaps": w["n_te_overlaps"],
            # driven mode
            "te_driven_mode": _te_driven_mode(
                s["hit"], s["te_class"], w["hit"], w["te_class"]
            ),
            # internal flags for split_outputs
            "_strict_hit": s["hit"],
            "_window_hit": w["hit"],
        }
        records.append(rec)
    return records


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@app.command()
def annotate(
    gtf_file: Path = typer.Argument(..., help="TACO-merged GTF file"),
    rmsk_bed: Path = typer.Argument(..., help="bgzipped+tabix RepeatMasker BED"),
    rmsk_annotation: Path = typer.Argument(..., help="TSV: name, class, family (no header)"),
    gencode_plus: Path = typer.Argument(..., help="Pickle: + strand gencode dict"),
    gencode_minus: Path = typer.Argument(..., help="Pickle: - strand gencode dict"),
    strict_bp: int = typer.Option(1, "--strict-bp", help="Strict window half-width"),
    window_bp: int = typer.Option(500, "--window-bp", help="Upstream window size"),
    include_single_exon: bool = typer.Option(False, "--include-single-exon"),
    output_dir: Path = typer.Option(Path("."), "--output-dir"),
    sample: str = typer.Option("", "--sample", help="Output file prefix (default: GTF stem)"),
) -> None:
    """Annotate TACO transcripts with TE and gene information."""
    # Validate inputs
    for p in (gtf_file, rmsk_bed, rmsk_annotation, gencode_plus, gencode_minus):
        if not p.exists():
            print(f"ERROR: file not found: {p}", file=sys.stderr)
            raise typer.Exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    sample = sample or gtf_file.stem

    print(f"[teprof3] Loading gencode dicts...", file=sys.stderr)
    with open(gencode_plus, "rb") as f:
        plus_dic = pickle.load(f)
    with open(gencode_minus, "rb") as f:
        minus_dic = pickle.load(f)
    gene_annotator = GeneAnnotator(plus_dic, minus_dic)

    print(f"[teprof3] Loading TE classification...", file=sys.stderr)
    te_ann_df = pd.read_csv(rmsk_annotation, sep="\t", header=None, names=["name", "class", "family"])
    te_classification = {
        row["name"]: (row["class"], row["family"])
        for _, row in te_ann_df.iterrows()
    }

    print(f"[teprof3] Parsing GTF: {gtf_file}", file=sys.stderr)
    transcript_df = parse_gtf(gtf_file)
    print(f"[teprof3] {len(transcript_df)} transcripts loaded", file=sys.stderr)

    print(f"[teprof3] Opening tabix file: {rmsk_bed}", file=sys.stderr)
    tbx = pysam.TabixFile(str(rmsk_bed))
    rmsk_handler = _TabixHandler(tbx)

    print(f"[teprof3] Annotating transcripts...", file=sys.stderr)
    records = build_records(
        transcript_df, gene_annotator, rmsk_handler, te_classification,
        strict_bp=strict_bp, window_bp=window_bp,
    )
    result_df = pd.DataFrame(records)

    strict_df, window_df, all_df = split_outputs(result_df, include_single_exon)

    all_path = output_dir / f"{sample}_all.tsv"
    strict_path = output_dir / f"{sample}_strict.tsv"
    window_path = output_dir / f"{sample}_window.tsv"

    all_df.to_csv(all_path, sep="\t", index=False)
    strict_df.to_csv(strict_path, sep="\t", index=False)
    window_df.to_csv(window_path, sep="\t", index=False)

    print(f"[teprof3] Written: {all_path} ({len(all_df)} rows)", file=sys.stderr)
    print(f"[teprof3] Written: {strict_path} ({len(strict_df)} rows)", file=sys.stderr)
    print(f"[teprof3] Written: {window_path} ({len(window_df)} rows)", file=sys.stderr)
