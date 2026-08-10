"""Compare paired SWORD2 benchmark score CSVs across experiment runs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_BASELINE_PATH = Path("benchmark/results_rust_original/scores.csv")
DEFAULT_EXPERIMENT_PATH = Path("benchmark/results_no_experiments/scores.csv")
EPSILON = 1e-12

HIGHER_IS_BETTER = {
    "ndo",
    "boundary_dist_score",
    "d_count_acc",
    "iou",
    "multi_ndo",
    "boundary_precision_5",
    "boundary_recall_5",
    "boundary_f1_5",
    "boundary_precision_10",
    "boundary_recall_10",
    "boundary_f1_10",
    "boundary_precision_20",
    "boundary_recall_20",
    "boundary_f1_20",
    "pred_coverage",
    "pairwise_precision",
    "pairwise_recall",
    "pairwise_f1",
    "adjusted_rand",
    "normalized_mutual_info",
    "matched_dice",
    "matched_jaccard",
    "exact_match",
}

LOWER_IS_BETTER = {
    "d_count_dev",
    "over_split",
    "merge",
    "median_boundary_error",
    "worst_boundary_error",
    "pred_linker_fraction",
    "variation_of_information",
    "runtime_s",
    "peak_rss_mb",
}

NEUTRAL_METRICS = {
    "n_residues",
    "n_true_domains",
    "n_pred_domains",
    "domain_count_bias",
}

NON_METRIC_COLUMNS = {
    "dataset",
    "entry_id",
    "pdb_id",
    "chain_id",
    "tool",
    "variant",
    "partition",
    "true_chopping",
    "pred_chopping",
}


def _select_variant(scores: pd.DataFrame, tool: str, variant: str) -> pd.DataFrame:
    selected = scores[(scores["tool"] == tool) & (scores["variant"] == variant)]
    return selected.drop_duplicates(subset=["entry_id"], keep="last").set_index("entry_id")


def metric_direction(metric: str) -> str:
    """Return whether higher, lower, or neither direction is better for a metric."""
    if metric in HIGHER_IS_BETTER:
        return "higher"
    if metric in LOWER_IS_BETTER:
        return "lower"
    if metric in NEUTRAL_METRICS:
        return "neutral"
    if any(token in metric for token in ["error", "dev", "runtime", "rss", "variation"]):
        return "lower"
    if any(
        token in metric
        for token in [
            "precision",
            "recall",
            "f1",
            "acc",
            "ndo",
            "iou",
            "dice",
            "jaccard",
            "coverage",
            "match",
            "rand",
            "info",
        ]
    ):
        return "higher"
    return "neutral"


def discover_numeric_metrics(
    baseline: pd.DataFrame,
    experiments: dict[str, pd.DataFrame],
) -> list[str]:
    """Return numeric score columns shared by all compared score tables."""
    common_columns = set(baseline.columns)
    for experiment_scores in experiments.values():
        common_columns &= set(experiment_scores.columns)

    metrics: list[str] = []
    for column in baseline.columns:
        if column not in common_columns or column in NON_METRIC_COLUMNS:
            continue
        frames = [baseline, *experiments.values()]
        if all(pd.api.types.is_numeric_dtype(frame[column]) for frame in frames):
            metrics.append(column)
    return metrics


def compare_experiment_scores(
    baseline: pd.DataFrame,
    experiments: dict[str, pd.DataFrame],
    metrics: list[str] | None = None,
    *,
    tool: str = "sword2-rust",
    variant: str = "optimal",
) -> list[dict[str, object]]:
    """Return paired score deltas for experiment score tables against a baseline."""
    metric_names = metrics or discover_numeric_metrics(baseline, experiments)
    baseline_rows = _select_variant(baseline, tool, variant)
    rows: list[dict[str, object]] = []

    for experiment_name, experiment_scores in experiments.items():
        experiment_rows = _select_variant(experiment_scores, tool, variant)
        for metric in metric_names:
            if metric not in baseline_rows.columns or metric not in experiment_rows.columns:
                continue

            paired = baseline_rows[[metric]].join(
                experiment_rows[[metric]],
                how="inner",
                lsuffix="_baseline",
                rsuffix="_experiment",
            )
            paired = paired.dropna()
            if paired.empty:
                continue

            baseline_col = f"{metric}_baseline"
            experiment_col = f"{metric}_experiment"
            baseline_values = pd.to_numeric(paired[baseline_col], errors="coerce")
            experiment_values = pd.to_numeric(paired[experiment_col], errors="coerce")
            deltas = experiment_values - baseline_values
            ties = deltas.abs() <= EPSILON
            direction = metric_direction(metric)
            if direction == "higher":
                experiment_better = deltas > EPSILON
                baseline_better = deltas < -EPSILON
            elif direction == "lower":
                experiment_better = deltas < -EPSILON
                baseline_better = deltas > EPSILON
            else:
                experiment_better = None
                baseline_better = None

            rows.append(
                {
                    "experiment": experiment_name,
                    "metric": metric,
                    "better_when": direction,
                    "n_pairs": int(len(paired)),
                    "baseline_mean": float(baseline_values.mean()),
                    "experiment_mean": float(experiment_values.mean()),
                    "mean_delta": float(deltas.mean()),
                    "median_delta": float(deltas.median()),
                    "experiment_better": None if experiment_better is None else int(experiment_better.sum()),
                    "ties": int(ties.sum()),
                    "baseline_better": None if baseline_better is None else int(baseline_better.sum()),
                    "experiment_better_rate": None
                    if experiment_better is None
                    else float(experiment_better.mean()),
                }
            )

    return rows


def compare_experiment_coverage(
    baseline: pd.DataFrame,
    experiments: dict[str, pd.DataFrame],
    *,
    tool: str = "sword2-rust",
    variant: str = "optimal",
) -> list[dict[str, object]]:
    """Return row coverage and partition-change counts for each experiment."""
    baseline_rows = _select_variant(baseline, tool, variant)
    rows: list[dict[str, object]] = []

    for experiment_name, experiment_scores in experiments.items():
        experiment_rows = _select_variant(experiment_scores, tool, variant)
        paired_index = baseline_rows.index.intersection(experiment_rows.index)
        changed_partitions = None
        if "pred_chopping" in baseline_rows.columns and "pred_chopping" in experiment_rows.columns:
            paired = baseline_rows.loc[paired_index, ["pred_chopping"]].join(
                experiment_rows.loc[paired_index, ["pred_chopping"]],
                lsuffix="_baseline",
                rsuffix="_experiment",
            )
            changed_partitions = int(
                (paired["pred_chopping_baseline"] != paired["pred_chopping_experiment"]).sum()
            )

        rows.append(
            {
                "experiment": experiment_name,
                "baseline_rows": int(len(baseline_rows)),
                "experiment_rows": int(len(experiment_rows)),
                "paired_rows": int(len(paired_index)),
                "experiment_only_rows": int(len(experiment_rows.index.difference(baseline_rows.index))),
                "baseline_only_rows": int(len(baseline_rows.index.difference(experiment_rows.index))),
                "changed_partitions": changed_partitions,
            }
        )

    return rows


def top_metric_changes(
    baseline: pd.DataFrame,
    experiment: pd.DataFrame,
    *,
    metric: str,
    tool: str = "sword2-rust",
    variant: str = "optimal",
    top_n: int = 10,
) -> list[dict[str, object]]:
    baseline_rows = _select_variant(baseline, tool, variant)
    experiment_rows = _select_variant(experiment, tool, variant)
    if metric not in baseline_rows.columns or metric not in experiment_rows.columns:
        return []

    extra_columns = [metric]
    for optional_column in ["n_pred_domains", "pred_chopping"]:
        if optional_column in baseline_rows.columns and optional_column in experiment_rows.columns:
            extra_columns.append(optional_column)

    paired = baseline_rows[extra_columns].join(
        experiment_rows[extra_columns],
        how="inner",
        lsuffix="_baseline",
        rsuffix="_experiment",
    )
    paired = paired.dropna(subset=[f"{metric}_baseline", f"{metric}_experiment"])
    if paired.empty:
        return []

    paired["delta"] = paired[f"{metric}_experiment"] - paired[f"{metric}_baseline"]
    paired["abs_delta"] = paired["delta"].abs()
    paired = paired.sort_values("abs_delta", ascending=False).head(top_n)

    rows: list[dict[str, object]] = []
    for entry_id, row in paired.iterrows():
        result = {
            "entry_id": entry_id,
            "baseline": float(row[f"{metric}_baseline"]),
            "experiment": float(row[f"{metric}_experiment"]),
            "delta": float(row["delta"]),
        }
        if "n_pred_domains_baseline" in row and "n_pred_domains_experiment" in row:
            result["baseline_domains"] = int(row["n_pred_domains_baseline"])
            result["experiment_domains"] = int(row["n_pred_domains_experiment"])
        if "pred_chopping_baseline" in row and "pred_chopping_experiment" in row:
            result["partition_changed"] = bool(row["pred_chopping_baseline"] != row["pred_chopping_experiment"])
        rows.append(result)
    return rows


def _format_float(value: object) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.6f}"


def _format_rate(value: object) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.3f}"


def _format_count(value: object) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return str(int(value))


def _table(rows: list[dict[str, object]], *, float_columns: set[str] | None = None) -> str:
    if not rows:
        return "No rows."
    frame = pd.DataFrame(rows)
    if "experiment_better_rate" in frame.columns:
        frame["experiment_better_rate"] = frame["experiment_better_rate"].map(_format_rate)
    for count_column in ["experiment_better", "baseline_better"]:
        if count_column in frame.columns:
            frame[count_column] = frame[count_column].map(_format_count)
    formatters = {column: _format_float for column in float_columns or set() if column in frame.columns}
    return frame.to_string(index=False, formatters=formatters)


def build_difference_report(
    baseline: pd.DataFrame,
    experiments: dict[str, pd.DataFrame],
    metrics: list[str] | None = None,
    *,
    tool: str = "sword2-rust",
    variant: str = "optimal",
    top_metric: str = "ndo",
    top_n: int = 10,
) -> str:
    metric_rows = compare_experiment_scores(baseline, experiments, metrics, tool=tool, variant=variant)
    coverage_rows = compare_experiment_coverage(baseline, experiments, tool=tool, variant=variant)

    lines = [
        "SWORD2 Difference Analysis",
        f"Tool: {tool}",
        f"Variant: {variant}",
        "",
        "Convention: mean_delta = experiment_mean - baseline_mean.",
        "For higher/lower metrics, experiment_better counts paired entries where the experiment improved.",
        "Neutral metrics are reported for context; their better counts are NA.",
        "",
        "Coverage",
        _table(coverage_rows),
        "",
        "Metric Deltas",
        _table(
            metric_rows,
            float_columns={"baseline_mean", "experiment_mean", "mean_delta", "median_delta"},
        ),
    ]

    if top_n > 0:
        for experiment_name, experiment_scores in experiments.items():
            rows = top_metric_changes(
                baseline,
                experiment_scores,
                metric=top_metric,
                tool=tool,
                variant=variant,
                top_n=top_n,
            )
            lines.extend(
                [
                    "",
                    f"Top {top_metric} Changes: {experiment_name}",
                    _table(rows, float_columns={"baseline", "experiment", "delta"}),
                ]
            )

    return "\n".join(lines)


def _parse_experiment(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("experiments must use NAME=PATH")
    name, path = value.split("=", 1)
    if not name:
        raise argparse.ArgumentTypeError("experiment name cannot be empty")
    return name, Path(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        type=Path,
        default=None,
        help=f"Baseline scores.csv (default: {DEFAULT_BASELINE_PATH})",
    )
    parser.add_argument(
        "--experiment",
        type=_parse_experiment,
        action="append",
        default=None,
        help="Experiment score CSV as NAME=PATH; may be repeated",
    )
    parser.add_argument("--out", type=Path, default=None, help="Optional output CSV")
    parser.add_argument(
        "--metrics",
        default=None,
        help="Comma-separated metric columns to compare (default: all shared numeric metrics)",
    )
    parser.add_argument("--top-metric", default="ndo", help="Metric used for the top-changes section")
    parser.add_argument("--top-n", type=int, default=10, help="Number of largest absolute changes to show")
    parser.add_argument("--tool", default="sword2-rust")
    parser.add_argument("--variant", default="optimal")
    return parser.parse_args()


def _resolve_inputs(args: argparse.Namespace) -> tuple[Path, list[tuple[str, Path]]]:
    baseline_path = args.baseline or DEFAULT_BASELINE_PATH
    if not baseline_path.exists():
        raise SystemExit(f"Baseline scores CSV does not exist: {baseline_path}")

    experiments = args.experiment
    if experiments is None:
        if not DEFAULT_EXPERIMENT_PATH.exists():
            raise SystemExit(
                "No --experiment was provided and the default experiment CSV does not exist: "
                f"{DEFAULT_EXPERIMENT_PATH}"
            )
        experiments = [(DEFAULT_EXPERIMENT_PATH.parent.name, DEFAULT_EXPERIMENT_PATH)]
    for _, path in experiments:
        if not path.exists():
            raise SystemExit(f"Experiment scores CSV does not exist: {path}")
    return baseline_path, experiments


def main() -> int:
    args = parse_args()
    baseline_path, experiment_paths = _resolve_inputs(args)
    baseline = pd.read_csv(baseline_path)
    experiments = {name: pd.read_csv(path) for name, path in experiment_paths}
    metrics = None if args.metrics is None else [metric.strip() for metric in args.metrics.split(",") if metric.strip()]
    rows = compare_experiment_scores(
        baseline,
        experiments,
        metrics=metrics,
        tool=args.tool,
        variant=args.variant,
    )
    results = pd.DataFrame(rows)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.out, index=False)
    print(
        build_difference_report(
            baseline,
            experiments,
            metrics=metrics,
            tool=args.tool,
            variant=args.variant,
            top_metric=args.top_metric,
            top_n=args.top_n,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
