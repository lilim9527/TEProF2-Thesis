"""Tests for query_te_at_tss — dual-mode TE detection with Unknown priority fix."""
from unittest.mock import MagicMock
from teprof3.annotation.te_query import query_te_at_tss


def _make_rmsk(records):
    handler = MagicMock()
    handler.query.return_value = records
    return handler


TE_CLASS = {
    "LINE1": ("LINE", "L1"),
    "AluSx": ("SINE", "Alu"),
    "Unknown_1": ("Unknown", "Unknown"),
}


def test_strict_hit():
    rmsk = _make_rmsk([("chr1", 999, 1100, "LINE1", 500.0, "+")])
    result = query_te_at_tss("chr1", 1000, "+", 1, 500, rmsk, TE_CLASS)
    assert result["strict"]["hit"] is True
    assert result["strict"]["primary_te_name"] == "LINE1"
    assert result["strict"]["te_class"] == "LINE"


def test_strict_miss_window_hit():
    rmsk = _make_rmsk([("chr1", 800, 900, "AluSx", 300.0, "+")])
    result = query_te_at_tss("chr1", 1000, "+", 1, 500, rmsk, TE_CLASS)
    assert result["strict"]["hit"] is False
    assert result["window"]["hit"] is True
    assert result["window"]["primary_te_name"] == "AluSx"


def test_unknown_fallback():
    """When only Unknown TEs overlap, still report them."""
    rmsk = _make_rmsk([("chr1", 999, 1001, "Unknown_1", 100.0, "+")])
    result = query_te_at_tss("chr1", 1000, "+", 1, 500, rmsk, TE_CLASS)
    assert result["strict"]["hit"] is True
    assert result["strict"]["primary_te_name"] == "Unknown_1"


def test_named_preferred_over_unknown():
    """Named TE wins over Unknown when both overlap strict window."""
    rmsk = _make_rmsk([
        ("chr1", 999, 1001, "Unknown_1", 100.0, "+"),
        ("chr1", 998, 1002, "LINE1", 200.0, "+"),
    ])
    result = query_te_at_tss("chr1", 1000, "+", 1, 500, rmsk, TE_CLASS)
    assert result["strict"]["primary_te_name"] == "LINE1"


def test_no_hit():
    rmsk = _make_rmsk([])
    result = query_te_at_tss("chr1", 1000, "+", 1, 500, rmsk, TE_CLASS)
    assert result["strict"]["hit"] is False
    assert result["strict"]["primary_te_name"] == "None"
    assert result["strict"]["te_start"] == 0
    assert result["strict"]["all_te_names"] == []


def test_minus_strand_window_direction():
    """Minus strand window extends downstream (higher coords)."""
    rmsk = _make_rmsk([("chr1", 1001, 1200, "LINE1", 500.0, "-")])
    result = query_te_at_tss("chr1", 1000, "-", 1, 500, rmsk, TE_CLASS)
    assert result["window"]["hit"] is True


def test_all_te_names_sorted_unknown_last():
    """all_te_names: named TEs before Unknown."""
    rmsk = _make_rmsk([
        ("chr1", 999, 1001, "Unknown_1", 100.0, "+"),
        ("chr1", 998, 1002, "LINE1", 200.0, "+"),
    ])
    result = query_te_at_tss("chr1", 1000, "+", 1, 500, rmsk, TE_CLASS)
    names = result["strict"]["all_te_names"]
    assert names[0] == "LINE1"
    assert "Unknown_1" in names


def test_n_te_overlaps():
    rmsk = _make_rmsk([
        ("chr1", 999, 1001, "LINE1", 200.0, "+"),
        ("chr1", 998, 1002, "AluSx", 150.0, "+"),
    ])
    result = query_te_at_tss("chr1", 1000, "+", 1, 500, rmsk, TE_CLASS)
    assert result["strict"]["n_te_overlaps"] == 2


def test_deepte_preferred_over_plain_unknown():
    """DeepTE-reclassified (rnd-N_family-N) wins over plain Unknown."""
    te_class_ext = {
        **TE_CLASS,
        "rnd-4_family-123": ("LINE", "L1"),
    }
    rmsk = _make_rmsk([
        ("chr1", 999, 1001, "Unknown_1", 50.0, "+"),
        ("chr1", 998, 1002, "rnd-4_family-123", 80.0, "+"),
    ])
    result = query_te_at_tss("chr1", 1000, "+", 1, 500, rmsk, te_class_ext)
    assert result["strict"]["primary_te_name"] == "rnd-4_family-123"
