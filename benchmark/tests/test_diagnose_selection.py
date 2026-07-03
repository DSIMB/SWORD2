"""Unit tests for the selection diagnostic's pure logic (no sklearn / no big CSVs)."""
from __future__ import annotations

import csv

import pytest

from benchmark.diagnose_selection import (
    FEATURE_COLS,
    Candidate,
    baseline_selectors,
    count_bias,
    gap_attribution,
    load_candidates,
    per_entry_pick_ndo,
    rank1,
    signature,
)


# --------------------------------------------------------------------------- #
# signature() — the numbering-invariant join key
# --------------------------------------------------------------------------- #
def test_signature_matches_across_dump_and_scores_formats():
    # dump: space-separated domains, ';'-separated segments.
    # scores: comma-separated domains, '_'-separated segments.
    assert signature("0-9 10-19") == signature("0-9,10-19")
    assert signature("0-44;120-165 45-119") == signature("0-44_120-165,45-119")


def test_signature_distinguishes_asymmetric_splits():
    assert signature("0-59 60-199") != signature("0-139 140-199")


def test_signature_is_invariant_to_a_numbering_offset():
    # Same relative structure shifted by a constant offset -> identical signature,
    # because it keys on segment *lengths*, not absolute indices.
    assert signature("0-9 10-29") == signature("100-109 110-129")


def test_signature_orders_domains_by_position_not_input_order():
    assert signature("60-99 0-59") == signature("0-59 60-99")


def test_signature_empty_is_none():
    assert signature("") is None
    assert signature('   "" ') is None


# --------------------------------------------------------------------------- #
# helpers to build synthetic candidates
# --------------------------------------------------------------------------- #
def mk(entry, variant, n_pred, n_true, ndo, feats=None):
    return Candidate(
        entry_id=entry,
        variant=variant,
        n_pred=n_pred,
        n_true=n_true,
        n_residues=100,
        ndo=ndo,
        iou=ndo,
        boundary_f1_10=ndo,
        features=feats,
    )


@pytest.fixture
def toy():
    # entry A: rank-1 (2 dom) NDO .70; a better 2-dom (.75); a 3-dom that is oracle
    #          and matches the true count (.90).
    return {
        "A": [
            mk("A", "optimal", 2, 3, 0.70),
            mk("A", "alternative", 2, 3, 0.75),
            mk("A", "alternative", 3, 3, 0.90),
        ]
    }


# --------------------------------------------------------------------------- #
# gap attribution
# --------------------------------------------------------------------------- #
def test_gap_attribution_components_sum_to_total(toy):
    g = gap_attribution(toy)
    assert g["within_count_component"] == pytest.approx(0.05)  # .75 - .70
    assert g["count_component"] == pytest.approx(0.15)  # .90 - .75
    assert g["total_gap"] == pytest.approx(0.20)  # oracle .90 - rank1 .70
    assert g["within_count_component"] + g["count_component"] == pytest.approx(g["total_gap"])


def test_gap_attribution_components_nonnegative_on_random_like_data():
    # rank-1 is always one candidate at its own count, oracle >= best-at-count,
    # so both components must be >= 0 by construction.
    by_entry = {
        "A": [mk("A", "optimal", 2, 2, 0.6), mk("A", "alternative", 4, 2, 0.4)],
        "B": [mk("B", "optimal", 1, 2, 0.5), mk("B", "alternative", 2, 2, 0.95)],
    }
    g = gap_attribution(by_entry)
    assert g["within_count_component"] >= 0
    assert g["count_component"] >= 0


def test_gap_attribution_reachability_flags(toy):
    g = gap_attribution(toy)
    assert g["frac_rank1_count_eq_true"] == pytest.approx(0.0)  # rank1=2, true=3
    assert g["frac_true_count_reachable"] == pytest.approx(1.0)  # a 3-dom candidate exists
    assert g["frac_rank1_count_is_best"] == pytest.approx(0.0)  # best-NDO count is 3, not 2


# --------------------------------------------------------------------------- #
# count bias
# --------------------------------------------------------------------------- #
def test_count_bias_detects_under_segmentation(toy):
    cb = count_bias(toy)
    assert cb["mean_delta"] == pytest.approx(-1.0)  # rank1 predicts 2, truth 3
    assert cb["frac_under"] == pytest.approx(1.0)
    assert cb["histogram"] == {-1: 1}


# --------------------------------------------------------------------------- #
# selectors
# --------------------------------------------------------------------------- #
def test_rank1_prefers_optimal_variant(toy):
    assert rank1(toy["A"]).ndo == pytest.approx(0.70)


def test_baseline_selectors_key_values(toy):
    b = baseline_selectors(toy)
    assert b["rank1 (current)"]["mean"] == pytest.approx(0.70)
    assert b["oracle (best NDO)"]["mean"] == pytest.approx(0.90)
    # a perfect count selector restricts to n_pred == n_true (3) -> only the .90 cand
    assert b["true count (perfect count)"]["mean"] == pytest.approx(0.90)


def test_per_entry_pick_is_argmax():
    by_entry = {"A": [mk("A", "optimal", 2, 2, 0.3), mk("A", "alternative", 3, 2, 0.8)]}
    picked = per_entry_pick_ndo(by_entry, lambda c: c.ndo)
    assert picked["A"] == pytest.approx(0.8)


# --------------------------------------------------------------------------- #
# join / loading from CSVs
# --------------------------------------------------------------------------- #
def _write_csv(path, fieldnames, rows):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def test_load_candidates_joins_features_by_signature(tmp_path):
    scores = tmp_path / "scores.csv"
    feats = tmp_path / "feats.csv"
    score_fields = [
        "tool", "variant", "entry_id", "pred_chopping",
        "n_pred_domains", "n_true_domains", "n_residues", "ndo", "iou", "boundary_f1_10",
    ]
    _write_csv(
        scores,
        score_fields,
        [
            # scores use comma/underscore formatting
            {"tool": "sword2-rust", "variant": "optimal", "entry_id": "X",
             "pred_chopping": "0-49,50-99", "n_pred_domains": 2, "n_true_domains": 2,
             "n_residues": 100, "ndo": 0.8, "iou": 0.8, "boundary_f1_10": 0.7},
            # a non-sword2 row that must be ignored
            {"tool": "merizo", "variant": "cuda", "entry_id": "X",
             "pred_chopping": "0-99", "n_pred_domains": 1, "n_true_domains": 2,
             "n_residues": 100, "ndo": 0.5, "iou": 0.5, "boundary_f1_10": 0.1},
        ],
    )
    feat_fields = ["entry_id", "delineation", *FEATURE_COLS]
    _write_csv(
        feats,
        feat_fields,
        [  # dump uses space/semicolon formatting for the same 2-domain partition
            {"entry_id": "X", "delineation": "0-49 50-99",
             **{c: 1.0 for c in FEATURE_COLS}, "max_cr": 0.42},
        ],
    )
    by_entry = load_candidates(scores, feats)
    assert set(by_entry) == {"X"}
    assert len(by_entry["X"]) == 1  # merizo row ignored
    cand = by_entry["X"][0]
    assert cand.features is not None
    assert cand.features["max_cr"] == pytest.approx(0.42)
    assert cand.derived["n_residues"] == pytest.approx(100.0)
