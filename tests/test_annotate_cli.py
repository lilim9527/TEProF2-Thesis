"""Tests for cli/output.py — 3-file split logic."""
import pandas as pd
from teprof3.cli.output import split_outputs as _split_outputs


def _sample_df():
    return pd.DataFrame({
        "transcript_id": ["T1", "T2", "T3", "T4"],
        "splicing": ["Yes", "Yes", "No", "No"],
        "_strict_hit": [True, False, True, False],
        "_window_hit": [True, True, True, False],
        "tpm": [5.0, 2.0, 1.0, 0.5],
    })


def test_split_outputs_strict_requires_splicing():
    strict, window, all_ = _split_outputs(_sample_df(), include_single_exon=False)
    assert list(strict["transcript_id"]) == ["T1"]


def test_split_outputs_window_requires_splicing():
    strict, window, all_ = _split_outputs(_sample_df(), include_single_exon=False)
    assert set(window["transcript_id"]) == {"T1", "T2"}


def test_split_outputs_all_no_filter():
    strict, window, all_ = _split_outputs(_sample_df(), include_single_exon=False)
    assert len(all_) == 4


def test_split_outputs_include_single_exon():
    strict, window, all_ = _split_outputs(_sample_df(), include_single_exon=True)
    assert "T3" in strict["transcript_id"].values


def test_split_outputs_drops_internal_columns():
    strict, window, all_ = _split_outputs(_sample_df(), include_single_exon=False)
    for df in (strict, window, all_):
        assert "_strict_hit" not in df.columns
        assert "_window_hit" not in df.columns
