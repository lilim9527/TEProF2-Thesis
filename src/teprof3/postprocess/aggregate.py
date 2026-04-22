"""
Post-processing aggregation and filtering for TEProf2 v3.

Replaces v1's R scripts with pure Python/pandas.
"""
from __future__ import annotations

import pandas as pd
from pathlib import Path

_KEEP_CLASSES = {"LINE", "SINE", "SINE?", "LTR", "LTR?", "DNA", "DNA?", "Retroposon"}


def _exon1_len(coords: str, _sentinel: int = 999_999) -> int:
    """Parse exon1 length from exon_coords string: strand,chr,e1_start,e1_end,...
    Returns a large sentinel on parse failure so the row is rejected by the filter."""
    try:
        p = coords.split(",")
        return int(p[3]) - int(p[2])
    except (IndexError, ValueError):
        return _sentinel


def filter_by_exon1_length(df: pd.DataFrame, max_bp: int = 2588) -> pd.DataFrame:
    """Keep rows where first exon length <= max_bp."""
    return df[df["exon_coords"].apply(_exon1_len) <= max_bp].reset_index(drop=True)


def filter_by_exon_skip(df: pd.DataFrame, max_skip: int = 2) -> pd.DataFrame:
    """Keep rows where n_exons - 2 <= max_skip (i.e. n_exons <= max_skip + 2)."""
    return df[df["n_exons"].astype(int) - 2 <= max_skip].reset_index(drop=True)


def filter_by_min_samples(df: pd.DataFrame, min_samples: int = 1) -> pd.DataFrame:
    """Keep transcripts present in at least min_samples distinct samples."""
    counts = df.groupby("transcript_id")["sample"].nunique()
    keep = counts[counts >= min_samples].index
    return df[df["transcript_id"].isin(keep)].reset_index(drop=True)


def filter_by_te_class(df: pd.DataFrame) -> pd.DataFrame:
    """Remove rows with Unknown or unrecognized TE classes."""
    return df[df["te_class"].isin(_KEEP_CLASSES)].reset_index(drop=True)


def filter_by_tpm(df: pd.DataFrame, min_tpm: float = 1.0) -> pd.DataFrame:
    """Keep rows with TPM >= min_tpm."""
    return df[df["tpm"].astype(float) >= min_tpm].reset_index(drop=True)


def filter_by_intron_reads(df: pd.DataFrame, min_reads: int = 1) -> pd.DataFrame:
    """Keep rows with intron_reads >= min_reads (skipped if column absent)."""
    if "intron_reads" not in df.columns:
        return df
    return df[df["intron_reads"].astype(int) >= min_reads].reset_index(drop=True)


def aggregate(
    input_dir: str | Path,
    mode: str,
    exon1_max_bp: int = 2588,
    exon_skip_max: int = 2,
    min_samples: int = 1,
    filter_te_class: bool = False,
    min_tpm: float = 1.0,
    min_intron_reads: int = 1,
    treatment_samples: str = "",
    treatment_exclusive: bool = False,
) -> pd.DataFrame:
    """
    Load all *_{mode}.tsv files from input_dir, apply filters, return combined DataFrame.

    Args:
        input_dir: Directory containing per-sample TSV files
        mode: File suffix to match ("strict", "window", or "all")
        exon1_max_bp: Max first-exon length in bp
        exon_skip_max: Max number of exon skips (n_exons - 2)
        min_samples: Min distinct samples a transcript must appear in
        filter_te_class: If True, remove Unknown/unrecognized TE classes
        min_tpm: Minimum TPM threshold
        min_intron_reads: Minimum intron-spanning reads
        treatment_samples: Comma-separated sample name substrings for treatment group
        treatment_exclusive: If True, keep only transcripts exclusive to treatment group

    Returns:
        Filtered DataFrame with all samples combined
    """
    files = sorted(Path(input_dir).glob(f"*_{mode}.tsv"))
    if not files:
        raise FileNotFoundError(f"No *_{mode}.tsv files found in {input_dir}")

    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        df["sample"] = f.stem.replace(f"_{mode}", "")
        frames.append(df)

    df = pd.concat(frames, ignore_index=True)
    df = filter_by_exon1_length(df, exon1_max_bp)
    df = filter_by_exon_skip(df, exon_skip_max)
    df = filter_by_tpm(df, min_tpm)
    df = filter_by_intron_reads(df, min_intron_reads)
    if filter_te_class:
        df = filter_by_te_class(df)
    df = filter_by_min_samples(df, min_samples)

    if treatment_samples:
        t_list = [s.strip() for s in treatment_samples.split(",")]
        df["_is_treatment"] = df["sample"].apply(
            lambda s: any(t in s for t in t_list)
        )
        if treatment_exclusive:
            t_ids = set(df[df["_is_treatment"]]["transcript_id"])
            c_ids = set(df[~df["_is_treatment"]]["transcript_id"])
            df = df[df["transcript_id"].isin(t_ids - c_ids)]
        df = df.drop(columns=["_is_treatment"])

    return df
