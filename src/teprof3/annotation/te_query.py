"""
TE query functions for TEProf2 v3.

Standalone dual-mode TE detection with Unknown priority fix.
Extracted from te_annotator.py (v2) for use in v3.
"""
from __future__ import annotations

import re

_DEEPTE_PATTERN = re.compile(r"rnd-\d+_family-\d+")

_NONE_TE: dict = {
    "hit": False, "primary_te_name": "None", "te_class": "None",
    "te_family": "None", "te_chrom": "None", "te_start": 0,
    "te_end": 0, "te_strand": "None", "te_score": 0.0,
    "all_te_names": [], "n_te_overlaps": 0,
}


def _select_primary(records: list) -> tuple | None:
    """Pick best TE record: named > DeepTE-reclassified > Unknown, then highest score."""
    if not records:
        return None
    named = [r for r in records if "Unknown" not in r[3]]
    deepte = [r for r in records if _DEEPTE_PATTERN.search(r[3])]
    pool = named if named else (deepte if deepte else records)
    return max(pool, key=lambda r: float(r[4]) if r[4] else 0.0)


def _build_te_result(records: list, te_classification: dict) -> dict:
    """Build result dict from a list of tabix records.

    Simple repeats (names not in te_classification) are filtered out before
    selecting the primary TE, so they never appear as hits.
    """
    if not records:
        return dict(_NONE_TE)
    # Filter to real TEs only (names present in te_classification)
    te_records = [r for r in records if r[3] in te_classification]
    if not te_records:
        return dict(_NONE_TE)
    primary = _select_primary(te_records)
    name = primary[3]
    te_class, te_family = te_classification.get(name, ("Unknown", "Unknown"))
    all_names = sorted(
        set(r[3] for r in te_records),
        key=lambda n: (1 if "Unknown" in n else 0, n),
    )
    return {
        "hit": True,
        "primary_te_name": name,
        "te_class": te_class,
        "te_family": te_family,
        "te_chrom": primary[0],
        "te_start": int(primary[1]),
        "te_end": int(primary[2]),
        "te_strand": primary[5] if len(primary) > 5 else "None",
        "te_score": float(primary[4]) if primary[4] else 0.0,
        "all_te_names": all_names,
        "n_te_overlaps": len(te_records),
    }


def query_te_at_tss(
    chrom: str,
    tss: int,
    strand: str,
    strict_bp: int,
    window_bp: int,
    rmsk_handler,
    te_classification: dict,
) -> dict:
    """
    Query TE overlaps at TSS in two modes.

    Args:
        chrom: Chromosome name
        tss: Transcription start site (0-based)
        strand: "+" or "-"
        strict_bp: Half-width for strict window (TSS ± strict_bp)
        window_bp: Upstream window size for window mode
        rmsk_handler: Object with .query(chrom, start, end) returning records
        te_classification: {te_name: (te_class, te_family)}

    Returns:
        {"strict": result_dict, "window": result_dict}
        Each result_dict has keys: hit, primary_te_name, te_class, te_family,
        te_chrom, te_start, te_end, te_strand, te_score, all_te_names, n_te_overlaps
    """
    s_start, s_end = tss - strict_bp, tss + strict_bp
    if strand == "+":
        w_start, w_end = tss - window_bp, tss
    else:
        w_start, w_end = tss, tss + window_bp

    strict_recs = [
        r for r in rmsk_handler.query(chrom, s_start, s_end)
        if int(r[1]) < s_end and int(r[2]) > s_start
    ]
    window_recs = [
        r for r in rmsk_handler.query(chrom, w_start, w_end)
        if int(r[1]) < w_end and int(r[2]) > w_start
    ]

    return {
        "strict": _build_te_result(strict_recs, te_classification),
        "window": _build_te_result(window_recs, te_classification),
    }
