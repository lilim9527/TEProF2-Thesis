"""Tests for postprocess/aggregate.py filter functions."""
import pandas as pd
from teprof3.postprocess.aggregate import (
    filter_by_exon1_length,
    filter_by_exon_skip,
    filter_by_min_samples,
    filter_by_te_class,
    filter_by_tpm,
)


def _df():
    return pd.DataFrame({
        "transcript_id": ["T1", "T2", "T3"],
        "exon_coords": [
            "+,chr1,100,200,200,500",       # exon1 len = 100
            "+,chr1,100,2700,2700,3000",    # exon1 len = 2600
            "+,chr1,100,150,150,200",       # exon1 len = 50
        ],
        "n_exons": [2, 3, 2],
        "tpm": [5.0, 0.5, 2.0],
        "te_class": ["LINE", "SINE", "Unknown"],
        "sample": ["S1", "S1", "S2"],
    })


def test_exon1_length_filters_long():
    """T2 has exon1 len 2600 which is <= 2588? No — 2600 > 2588, should be filtered."""
    result = filter_by_exon1_length(_df(), 2588)
    assert "T2" not in result["transcript_id"].values
    assert len(result) == 2


def test_exon1_length_keeps_short():
    result = filter_by_exon1_length(_df(), 2588)
    assert "T1" in result["transcript_id"].values
    assert "T3" in result["transcript_id"].values


def test_exon_skip_filters_excess():
    df = _df()
    df["n_exons"] = [2, 5, 3]
    # max_skip=2 means n_exons - 2 <= 2, so n_exons <= 4; T2 (5 exons) filtered
    result = filter_by_exon_skip(df, 2)
    assert "T2" not in result["transcript_id"].values


def test_exon_skip_keeps_within_limit():
    df = _df()
    df["n_exons"] = [2, 4, 3]
    result = filter_by_exon_skip(df, 2)
    assert len(result) == 3


def test_min_samples_filters_insufficient():
    df = _df()
    df["sample"] = ["S1", "S1", "S1"]  # all same sample
    result = filter_by_min_samples(df, 2)
    assert len(result) == 0


def test_min_samples_keeps_sufficient():
    df = _df()
    # T1 appears in S1 only, T3 in S2 only — neither meets min_samples=2
    # But if we duplicate T1 across samples:
    df2 = pd.concat([df, df.assign(sample="S2")], ignore_index=True)
    result = filter_by_min_samples(df2, 2)
    # All transcripts now appear in both S1 and S2
    assert len(result) > 0


def test_te_class_removes_unknown():
    result = filter_by_te_class(_df())
    assert "Unknown" not in result["te_class"].values


def test_te_class_keeps_known():
    result = filter_by_te_class(_df())
    assert "LINE" in result["te_class"].values
    assert "SINE" in result["te_class"].values


def test_tpm_filters_low():
    result = filter_by_tpm(_df(), 1.0)
    assert "T2" not in result["transcript_id"].values  # tpm=0.5


def test_tpm_keeps_high():
    result = filter_by_tpm(_df(), 1.0)
    assert "T1" in result["transcript_id"].values  # tpm=5.0
