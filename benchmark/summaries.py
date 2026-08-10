from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


HEADLINE_METRICS = [
    "ndo",
    "boundary_dist_score",
    "boundary_f1_10",
    "matched_dice",
    "pairwise_f1",
    "d_count_acc",
    "iou",
]


def headline_scores(scores: pd.DataFrame) -> pd.DataFrame:
    sword_optimal = scores["tool"].str.startswith("sword2") & (scores["variant"] == "optimal")
    one_shot = ~scores["tool"].str.startswith("sword2")
    return scores[sword_optimal | one_shot].copy()


def _bootstrap_ci(values: pd.Series, seed: int = 0, n_bootstrap: int = 1000) -> tuple[float, float]:
    clean = values.dropna().to_numpy(dtype=float)
    if len(clean) == 0:
        return float("nan"), float("nan")
    if len(clean) == 1:
        return float(clean[0]), float(clean[0])
    rng = np.random.default_rng(seed)
    samples = rng.choice(clean, size=(n_bootstrap, len(clean)), replace=True).mean(axis=1)
    low, high = np.quantile(samples, [0.025, 0.975])
    return float(low), float(high)


def _summary_rows(scores: pd.DataFrame, metrics: list[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    available_metrics = [metric for metric in metrics if metric in scores.columns]
    for (tool, variant), group in scores.groupby(["tool", "variant"], dropna=False):
        for metric in available_metrics:
            ci_low, ci_high = _bootstrap_ci(group[metric])
            rows.append(
                {
                    "tool": tool,
                    "variant": variant,
                    "metric": metric,
                    "n": int(group[metric].dropna().shape[0]),
                    "mean": float(group[metric].mean()),
                    "median": float(group[metric].median()),
                    "ci95_low": ci_low,
                    "ci95_high": ci_high,
                }
            )
    return rows


def write_metric_summaries(scores: pd.DataFrame, output_dir: Path) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [
        output_dir / "headline_metric_summary.csv",
        output_dir / "all_partition_metric_summary.csv",
    ]
    pd.DataFrame(_summary_rows(headline_scores(scores), HEADLINE_METRICS)).to_csv(paths[0], index=False)
    pd.DataFrame(_summary_rows(scores, HEADLINE_METRICS)).to_csv(paths[1], index=False)
    return paths


def write_win_rates(scores: pd.DataFrame, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    data = headline_scores(scores)
    rows: list[dict[str, object]] = []
    tools = sorted(data["tool"].dropna().unique())
    metrics = [metric for metric in HEADLINE_METRICS if metric in data.columns]
    for metric in metrics:
        pivot = data.pivot_table(index="entry_id", columns="tool", values=metric, aggfunc="first")
        for left in tools:
            for right in tools:
                if left >= right or left not in pivot or right not in pivot:
                    continue
                paired = pivot[[left, right]].dropna()
                if paired.empty:
                    continue
                delta = paired[left] - paired[right]
                rows.append(
                    {
                        "metric": metric,
                        "tool_a": left,
                        "tool_b": right,
                        "n_pairs": int(len(paired)),
                        "tool_a_win_rate": float((delta > 0).mean()),
                        "tool_b_win_rate": float((delta < 0).mean()),
                        "tie_rate": float((delta == 0).mean()),
                        "median_delta_a_minus_b": float(delta.median()),
                    }
                )
    path = output_dir / "headline_win_rates.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path
