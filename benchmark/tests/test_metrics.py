from pathlib import Path

import pandas as pd
import pytest

from benchmark.metrics import chopping_from_named_bounds, score_choppings


CHAINSAW_BENCHMARK = Path(
    "/home/chili/cretin/PROJECTS/chainsaw/data_and_benchmarking/"
    "chainsaw_model_v3_on_cath1363_test.csv"
)


def test_chainsaw_published_metrics_are_reproduced_for_sample_rows():
    rows = pd.read_csv(CHAINSAW_BENCHMARK).head(25)

    for row in rows.itertuples(index=False):
        true = chopping_from_named_bounds(row.true_dbounds, row.true_dnames)
        pred = chopping_from_named_bounds(row.pred_dbounds, row.pred_dnames)
        metrics = score_choppings(true, pred)

        assert metrics.ndo == pytest.approx(row.ndo, abs=1e-12)
        assert metrics.boundary_dist_score == pytest.approx(row.boundary_dist_score, abs=1e-12)
        assert metrics.d_count_acc == row.d_count_acc
        assert metrics.d_count_dev == row.d_count_dev
        if pd.isna(row.multi_ndo):
            assert pd.isna(metrics.multi_ndo)
        else:
            assert metrics.multi_ndo == pytest.approx(row.multi_ndo, abs=1e-12)


def test_metrics_support_pipe_and_comma_domain_delimiters():
    pipe_metrics = score_choppings("0-9|10-19", "0-8|9-19", 20)
    comma_metrics = score_choppings("0-9,10-19", "0-8,9-19", 20)

    assert pipe_metrics == comma_metrics
    assert pipe_metrics.n_true_domains == 2
    assert pipe_metrics.n_pred_domains == 2


def test_extended_metrics_are_perfect_for_identical_choppings():
    metrics = score_choppings("0-4,5-9", "0-4,5-9", 10)

    assert metrics.domain_count_bias == 0.0
    assert metrics.over_split == 0.0
    assert metrics.merge == 0.0
    assert metrics.boundary_f1_10 == 1.0
    assert metrics.median_boundary_error == 0.0
    assert metrics.worst_boundary_error == 0.0
    assert metrics.pred_coverage == 1.0
    assert metrics.pred_linker_fraction == 0.0
    assert metrics.pairwise_f1 == 1.0
    assert metrics.adjusted_rand == 1.0
    assert metrics.normalized_mutual_info == 1.0
    assert metrics.variation_of_information == 0.0
    assert metrics.matched_dice == 1.0
    assert metrics.matched_jaccard == 1.0
    assert metrics.exact_match == 1.0


def test_extended_metrics_flag_merge_and_linker_coverage():
    metrics = score_choppings("0-4,5-9", "0-3", 10)

    assert metrics.domain_count_bias == -1.0
    assert metrics.over_split == 0.0
    assert metrics.merge == 1.0
    assert metrics.pred_coverage == 0.4
    assert metrics.pred_linker_fraction == 0.6
    assert metrics.boundary_recall_10 == 0.0
    assert metrics.exact_match == 0.0


def test_score_choppings_tolerates_shared_boundary_residue_once():
    metrics = score_choppings("0-4,4-9", "0-4,5-9", 10)

    assert metrics.n_true_domains == 2
    assert metrics.n_pred_domains == 2
    assert metrics.pred_coverage == 1.0
