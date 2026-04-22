"""CLI output splitting for TEProf2 v3 — 3-file output logic."""
from __future__ import annotations

import pandas as pd


def split_outputs(
    df: pd.DataFrame, include_single_exon: bool
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Split annotated DataFrame into three output files.

    Args:
        df: Annotated DataFrame with columns: splicing, _strict_hit, _window_hit
        include_single_exon: If True, single-exon transcripts pass the splicing filter

    Returns:
        (strict_df, window_df, all_df) — internal columns dropped from all three.
    """
    splicing_ok = (df["splicing"] == "Yes") | include_single_exon
    drop_cols = ["_strict_hit", "_window_hit"]

    strict = df[splicing_ok & df["_strict_hit"]].drop(columns=drop_cols, errors="ignore").reset_index(drop=True)
    window = df[splicing_ok & df["_window_hit"]].drop(columns=drop_cols, errors="ignore").reset_index(drop=True)
    all_ = df.drop(columns=drop_cols, errors="ignore").reset_index(drop=True)
    return strict, window, all_
