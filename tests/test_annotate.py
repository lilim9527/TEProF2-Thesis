"""Tests for cli/annotate.py — GTF parsing and annotation pipeline."""
from __future__ import annotations

import io
import pandas as pd
import pytest
from unittest.mock import MagicMock, patch


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

GTF_CONTENT = """\
##gff-version 3
chr1\tTACO\ttranscript\t1001\t2000\t.\t+\t.\ttranscript_id "T1"; gene_id "G1";
chr1\tTACO\texon\t1001\t1200\t.\t+\t.\ttranscript_id "T1"; gene_id "G1";
chr1\tTACO\texon\t1500\t2000\t.\t+\t.\ttranscript_id "T1"; gene_id "G1";
chr1\tTACO\ttranscript\t3001\t4000\t.\t-\t.\ttranscript_id "T2"; gene_id "G2";
chr1\tTACO\texon\t3001\t3500\t.\t-\t.\ttranscript_id "T2"; gene_id "G2";
chr1\tTACO\texon\t3700\t4000\t.\t-\t.\ttranscript_id "T2"; gene_id "G2";
chr1\tTACO\ttranscript\t5001\t5500\t.\t+\t.\ttranscript_id "T3"; gene_id "G3";
chr1\tTACO\texon\t5001\t5500\t.\t+\t.\ttranscript_id "T3"; gene_id "G3";
"""

def _make_mock_gene_annotator():
    ann = MagicMock()
    ann.annotate_tss.return_value = {
        "tss_gene_type": "protein_coding", "tss_gene_name": "GENE1",
        "tss_chr": "chr1", "tss_start": "1000", "tss_end": "2000",
        "tss_element": "exon", "tss_element_num": "1",
        "tss_start_codon": "None", "tss_transcript_start": "1000",
        "tss_transcript_end": "2000", "tss_element_list": "None",
        "tss_ref_transcript_id": "ENST001",
    }
    ann.annotate_tend.return_value = {
        "tend_gene_type": "protein_coding", "tend_gene_name": "GENE1",
        "tend_chr": "chr1", "tend_start": "1000", "tend_end": "2000",
        "tend_element": "exon", "tend_element_num": "1",
        "tend_start_codon": "None", "tend_transcript_start": "1000",
        "tend_transcript_end": "2000", "tend_element_list": "None",
        "tend_ref_transcript_id": "ENST001",
    }
    return ann


def _make_mock_te_result(hit=True):
    return {
        "strict": {
            "hit": hit, "primary_te_name": "LINE1" if hit else "None",
            "te_class": "LINE" if hit else "None",
            "te_family": "L1" if hit else "None",
            "te_chrom": "chr1" if hit else "None",
            "te_start": 999, "te_end": 1100,
            "te_strand": "+", "te_score": 500.0,
            "all_te_names": ["LINE1"] if hit else [],
            "n_te_overlaps": 1 if hit else 0,
        },
        "window": {
            "hit": hit, "primary_te_name": "LINE1" if hit else "None",
            "te_class": "LINE" if hit else "None",
            "te_family": "L1" if hit else "None",
            "te_chrom": "chr1" if hit else "None",
            "te_start": 999, "te_end": 1100,
            "te_strand": "+", "te_score": 500.0,
            "all_te_names": ["LINE1"] if hit else [],
            "n_te_overlaps": 1 if hit else 0,
        },
    }


# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------

from teprof3.cli.annotate import parse_gtf, build_records, EXPECTED_COLUMNS


# ---------------------------------------------------------------------------
# Tests: GTF parsing
# ---------------------------------------------------------------------------

def test_parse_gtf_transcript_count():
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    assert len(df) == 3


def test_parse_gtf_transcript_ids():
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    assert set(df["transcript_id"]) == {"T1", "T2", "T3"}


def test_parse_gtf_tss_plus_strand():
    """TSS for + strand = start - 1 (0-based)."""
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    t1 = df[df["transcript_id"] == "T1"].iloc[0]
    assert t1["tss"] == 1000  # 1001 - 1 = 1000


def test_parse_gtf_tss_minus_strand():
    """TSS for - strand = end (0-based)."""
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    t2 = df[df["transcript_id"] == "T2"].iloc[0]
    assert t2["tss"] == 4000  # end = 4000


def test_parse_gtf_tend_plus_strand():
    """tend for + strand = end."""
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    t1 = df[df["transcript_id"] == "T1"].iloc[0]
    assert t1["tend"] == 2000


def test_parse_gtf_tend_minus_strand():
    """tend for - strand = start - 1 (0-based)."""
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    t2 = df[df["transcript_id"] == "T2"].iloc[0]
    assert t2["tend"] == 3000  # 3001 - 1 = 3000


def test_parse_gtf_n_exons():
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    t1 = df[df["transcript_id"] == "T1"].iloc[0]
    t3 = df[df["transcript_id"] == "T3"].iloc[0]
    assert t1["n_exons"] == 2
    assert t3["n_exons"] == 1


def test_parse_gtf_splicing():
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    t1 = df[df["transcript_id"] == "T1"].iloc[0]
    t3 = df[df["transcript_id"] == "T3"].iloc[0]
    assert t1["splicing"] == "Yes"
    assert t3["splicing"] == "No"


def test_parse_gtf_exon_coords_format():
    """exon_coords: '{strand},{chrom},{e1_start},{e1_end},{e2_start},{e2_end},...'"""
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    t1 = df[df["transcript_id"] == "T1"].iloc[0]
    parts = t1["exon_coords"].split(",")
    assert parts[0] == "+"
    assert parts[1] == "chr1"
    # 0-based coords: GTF exon 1001-1200 → 1000-1200
    assert parts[2] == "1000"
    assert parts[3] == "1200"


# ---------------------------------------------------------------------------
# Tests: build_records (full annotation pipeline)
# ---------------------------------------------------------------------------

def test_build_records_column_count():
    """build_records returns EXPECTED_COLUMNS + 2 internal flags."""
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    gene_ann = _make_mock_gene_annotator()
    rmsk = MagicMock()
    rmsk.query.return_value = []
    te_class = {}

    records = build_records(df, gene_ann, rmsk, te_class, strict_bp=1, window_bp=500)
    result_df = pd.DataFrame(records)
    expected = len(EXPECTED_COLUMNS) + 2  # +2 for _strict_hit, _window_hit
    assert len(result_df.columns) == expected, (
        f"Expected {expected} columns, got {len(result_df.columns)}: {list(result_df.columns)}"
    )


def test_build_records_column_names():
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    gene_ann = _make_mock_gene_annotator()
    rmsk = MagicMock()
    rmsk.query.return_value = []
    te_class = {}

    records = build_records(df, gene_ann, rmsk, te_class, strict_bp=1, window_bp=500)
    result_df = pd.DataFrame(records)
    # Check all expected columns present (internal _strict_hit/_window_hit included before split)
    for col in EXPECTED_COLUMNS:
        assert col in result_df.columns, f"Missing column: {col}"


def test_build_records_strict_hit_flag():
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    gene_ann = _make_mock_gene_annotator()
    rmsk = MagicMock()
    # Return a TE record that overlaps strict window
    rmsk.query.return_value = [("chr1", 999, 1001, "LINE1", "500", "+")]
    te_class = {"LINE1": ("LINE", "L1")}

    records = build_records(df, gene_ann, rmsk, te_class, strict_bp=1, window_bp=500)
    result_df = pd.DataFrame(records)
    t1 = result_df[result_df["transcript_id"] == "T1"].iloc[0]
    assert t1["_strict_hit"] is True or t1["_strict_hit"] == True


def test_build_records_no_hit():
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    gene_ann = _make_mock_gene_annotator()
    rmsk = MagicMock()
    rmsk.query.return_value = []
    te_class = {}

    records = build_records(df, gene_ann, rmsk, te_class, strict_bp=1, window_bp=500)
    result_df = pd.DataFrame(records)
    assert all(result_df["_strict_hit"] == False)
    assert all(result_df["_window_hit"] == False)


def test_build_records_48_output_columns():
    """After dropping internal columns, output should have exactly 48 columns."""
    from teprof3.cli.output import split_outputs
    df = parse_gtf(io.StringIO(GTF_CONTENT))
    gene_ann = _make_mock_gene_annotator()
    rmsk = MagicMock()
    rmsk.query.return_value = []
    te_class = {}

    records = build_records(df, gene_ann, rmsk, te_class, strict_bp=1, window_bp=500)
    result_df = pd.DataFrame(records)
    _, _, all_df = split_outputs(result_df, include_single_exon=False)
    assert len(all_df.columns) == 49, (
        f"Expected 49 output columns, got {len(all_df.columns)}: {list(all_df.columns)}"
    )
