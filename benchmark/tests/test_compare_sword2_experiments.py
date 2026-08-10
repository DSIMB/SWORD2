import pandas as pd
import pytest

from benchmark.compare_sword2_experiments import (
    build_difference_report,
    compare_experiment_coverage,
    compare_experiment_scores,
)


def test_compare_experiment_scores_reports_paired_deltas():
    baseline = pd.DataFrame(
        [
            {"entry_id": "a", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.5},
            {"entry_id": "b", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.8},
        ]
    )
    experiment = pd.DataFrame(
        [
            {"entry_id": "a", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.7},
            {"entry_id": "b", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.7},
        ]
    )

    rows = compare_experiment_scores(baseline, {"semantic": experiment}, metrics=["ndo"])

    assert rows == [
        {
            "experiment": "semantic",
            "metric": "ndo",
            "better_when": "higher",
            "n_pairs": 2,
            "baseline_mean": pytest.approx(0.65),
            "experiment_mean": pytest.approx(0.7),
            "mean_delta": pytest.approx(0.05),
            "median_delta": pytest.approx(0.05),
            "experiment_better": 1,
            "ties": 0,
            "baseline_better": 1,
            "experiment_better_rate": pytest.approx(0.5),
        }
    ]


def test_compare_experiment_scores_defaults_to_all_numeric_metrics_and_direction():
    baseline = pd.DataFrame(
        [
            {
                "entry_id": "a",
                "tool": "sword2-rust",
                "variant": "optimal",
                "ndo": 0.5,
                "d_count_dev": 2.0,
                "runtime_s": 1.0,
            },
            {
                "entry_id": "b",
                "tool": "sword2-rust",
                "variant": "optimal",
                "ndo": 0.8,
                "d_count_dev": 0.0,
                "runtime_s": 1.0,
            },
        ]
    )
    experiment = pd.DataFrame(
        [
            {
                "entry_id": "a",
                "tool": "sword2-rust",
                "variant": "optimal",
                "ndo": 0.7,
                "d_count_dev": 1.0,
                "runtime_s": 2.0,
            },
            {
                "entry_id": "b",
                "tool": "sword2-rust",
                "variant": "optimal",
                "ndo": 0.7,
                "d_count_dev": 1.0,
                "runtime_s": 1.0,
            },
        ]
    )

    rows = compare_experiment_scores(baseline, {"current": experiment})
    by_metric = {row["metric"]: row for row in rows}

    assert set(by_metric) == {"ndo", "d_count_dev", "runtime_s"}
    assert by_metric["ndo"]["better_when"] == "higher"
    assert by_metric["ndo"]["experiment_better"] == 1
    assert by_metric["ndo"]["baseline_better"] == 1
    assert by_metric["d_count_dev"]["better_when"] == "lower"
    assert by_metric["d_count_dev"]["experiment_better"] == 1
    assert by_metric["d_count_dev"]["baseline_better"] == 1
    assert by_metric["runtime_s"]["better_when"] == "lower"
    assert by_metric["runtime_s"]["experiment_better"] == 0
    assert by_metric["runtime_s"]["baseline_better"] == 1


def test_difference_report_includes_coverage_and_metric_tables():
    baseline = pd.DataFrame(
        [
            {"entry_id": "a", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.5, "pred_chopping": "0-9"},
            {"entry_id": "b", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.8, "pred_chopping": "0-4,5-9"},
        ]
    )
    experiment = pd.DataFrame(
        [
            {"entry_id": "a", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.7, "pred_chopping": "0-4,5-9"},
            {"entry_id": "c", "tool": "sword2-rust", "variant": "optimal", "ndo": 0.9, "pred_chopping": "0-9"},
        ]
    )

    coverage = compare_experiment_coverage(baseline, {"current": experiment})
    report = build_difference_report(baseline, {"current": experiment}, top_n=1)

    assert coverage == [
        {
            "experiment": "current",
            "baseline_rows": 2,
            "experiment_rows": 2,
            "paired_rows": 1,
            "experiment_only_rows": 1,
            "baseline_only_rows": 1,
            "changed_partitions": 1,
        }
    ]
    assert "Coverage" in report
    assert "Metric Deltas" in report
    assert "Top ndo Changes" in report
    assert "current" in report
