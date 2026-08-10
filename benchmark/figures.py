from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/sword2-benchmark-matplotlib")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from benchmark.summaries import HEADLINE_METRICS, headline_scores


TOOL_ORDER = ["sword2-rust", "sword2-original", "merizo", "chainsaw"]
TOOL_PALETTE = {
    "sword2-rust": "#1f77b4",
    "sword2-original": "#d62728",
    "merizo": "#2ca02c",
    "chainsaw": "#9467bd",
}
METRIC_LABELS = {
    "ndo": "NDO",
    "boundary_dist_score": "CASP boundary score",
    "boundary_f1_10": "Boundary F1 (+/-10 residues)",
    "matched_dice": "Matched domain Dice",
    "pairwise_f1": "Pairwise same-domain F1",
    "d_count_acc": "Exact domain-count accuracy",
    "iou": "Weighted best-domain IoU",
}


def _setup_theme() -> None:
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 220,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _ordered_tools(data: pd.DataFrame) -> list[str]:
    present = set(data["tool"].dropna().unique())
    ordered = [tool for tool in TOOL_ORDER if tool in present]
    ordered.extend(sorted(present.difference(ordered)))
    return ordered


def _savefig(path: Path) -> list[Path]:
    path.parent.mkdir(parents=True, exist_ok=True)
    png_path = path.with_suffix(".png")
    plt.tight_layout()
    plt.savefig(png_path)
    plt.close()
    return [png_path]


def _violinplot_with_points(
    *,
    data: pd.DataFrame,
    x: str,
    y: str,
    order: list[str] | None = None,
    hue: str | None = None,
    hue_order: list[str] | None = None,
    palette: dict[str, str] | None = None,
    legend: str | bool = "auto",
    dodge_points: bool = False,
    point_alpha: float = 0.16,
    point_size: float = 2.3,
) -> None:
    violin_kwargs: dict[str, object] = {
        "data": data,
        "x": x,
        "y": y,
        "order": order,
        "cut": 0,
        "inner": "quartile",
        "density_norm": "width",
        "linewidth": 1,
        "saturation": 0.8,
    }
    if hue is not None:
        violin_kwargs.update(
            {
                "hue": hue,
                "hue_order": hue_order,
                "palette": palette,
                "legend": legend,
            }
        )
    sns.violinplot(**violin_kwargs)

    strip_kwargs: dict[str, object] = {
        "data": data,
        "x": x,
        "y": y,
        "order": order,
        "color": "#222222",
        "alpha": point_alpha,
        "size": point_size,
        "jitter": 0.25,
    }
    if hue is not None and dodge_points:
        strip_kwargs.update(
            {
                "hue": hue,
                "hue_order": hue_order,
                "dodge": True,
                "palette": {tool: "#222222" for tool in data[hue].dropna().unique()},
                "legend": False,
            }
        )
        strip_kwargs.pop("color")
    sns.stripplot(**strip_kwargs)


def plot_metric_distributions(scores: pd.DataFrame, output_dir: Path) -> list[Path]:
    _setup_theme()
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    data = headline_scores(scores)
    for metric in [metric for metric in HEADLINE_METRICS if metric in data.columns]:
        path = output_dir / f"headline_{metric}_by_tool.png"
        plt.figure(figsize=(10, 5.5))
        order = _ordered_tools(data)
        _violinplot_with_points(
            data=data,
            x="tool",
            y=metric,
            hue="tool",
            order=order,
            hue_order=order,
            palette=TOOL_PALETTE,
            legend=False,
        )
        plt.xlabel("")
        plt.ylabel(METRIC_LABELS.get(metric, metric))
        plt.title(f"{METRIC_LABELS.get(metric, metric)} by tool (headline predictions)")
        plt.xticks(rotation=20, ha="right")
        paths.extend(_savefig(path))
    return paths


def _run_level_data(scores: pd.DataFrame, runs: pd.DataFrame | None = None) -> pd.DataFrame:
    lengths = scores[["entry_id", "tool", "n_residues"]].drop_duplicates()
    score_data = headline_scores(scores).drop_duplicates(["entry_id", "tool"]).copy()
    if runs is not None and not runs.empty:
        data = runs.merge(lengths, on=["entry_id", "tool"], how="left")
        fallback_cols = [
            col
            for col in ["entry_id", "tool", "n_residues", "runtime_s", "peak_rss_mb"]
            if col in score_data.columns
        ]
        fallback = score_data[fallback_cols]
        data = data.merge(fallback, on=["entry_id", "tool"], how="left", suffixes=("", "_score"))
        for column in ["n_residues", "runtime_s", "peak_rss_mb"]:
            score_column = f"{column}_score"
            if score_column in data.columns:
                if column in data.columns:
                    data[column] = data[column].combine_first(data[score_column])
                else:
                    data[column] = data[score_column]
                data = data.drop(columns=[score_column])
        return data.copy()
    return score_data


def plot_runtime_memory(scores: pd.DataFrame, output_dir: Path, runs: pd.DataFrame | None = None) -> list[Path]:
    _setup_theme()
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    data = _run_level_data(scores, runs=runs)
    for metric in ["runtime_s", "peak_rss_mb"]:
        if metric not in data.columns or data[metric].dropna().empty:
            continue
        plot_data = data.dropna(subset=[metric, "n_residues"]).copy()
        path = output_dir / f"{metric}_vs_length.png"
        plt.figure(figsize=(9, 5.5))
        sns.scatterplot(
            data=plot_data,
            x="n_residues",
            y=metric,
            hue="tool",
            hue_order=_ordered_tools(plot_data),
            palette=TOOL_PALETTE,
            alpha=0.75,
            s=36,
        )
        plt.yscale("log")
        plt.xlabel("Residues")
        plt.ylabel("Runtime (s, log scale)" if metric == "runtime_s" else "Peak RSS (MB, log scale)")
        plt.title("Run-level resource usage by protein length")
        plt.legend(title="Tool", frameon=True)
        paths.extend(_savefig(path))
    return paths


def plot_sword2_oracle_gap(scores: pd.DataFrame, output_dir: Path) -> list[Path]:
    _setup_theme()
    data = scores[scores["tool"].str.startswith("sword2") & scores["variant"].isin(["optimal", "oracle"])]
    metrics = [metric for metric in ["ndo", "boundary_f1_10", "matched_dice", "pairwise_f1"] if metric in data.columns]
    rows: list[dict[str, object]] = []
    for metric in metrics:
        pivot = data.pivot_table(index=["entry_id", "tool"], columns="variant", values=metric, aggfunc="first")
        if "optimal" not in pivot or "oracle" not in pivot:
            continue
        gap = (pivot["oracle"] - pivot["optimal"]).dropna()
        for (entry_id, tool), value in gap.items():
            rows.append({"entry_id": entry_id, "tool": tool, "metric": METRIC_LABELS.get(metric, metric), "gap": value})
    if not rows:
        return []
    plot_data = pd.DataFrame(rows)
    plt.figure(figsize=(10, 5.5))
    _violinplot_with_points(
        data=plot_data,
        x="metric",
        y="gap",
        hue="tool",
        hue_order=_ordered_tools(plot_data),
        palette=TOOL_PALETTE,
        dodge_points=True,
        point_alpha=0.12,
        point_size=2,
    )
    plt.legend(title="Tool", frameon=True, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)
    plt.xlabel("")
    plt.ylabel("Oracle minus optimal")
    plt.title("Potential gain from SWORD2 alternatives")
    plt.xticks(rotation=20, ha="right")
    return _savefig(output_dir / "sword2_oracle_gap.png")


def _partition_rank(partition: str) -> int | None:
    if partition == "Optimal partition":
        return 1
    prefix = "Alternative partition "
    if partition.startswith(prefix):
        return int(partition.removeprefix(prefix)) + 1
    return None


def plot_sword2_topk_curve(scores: pd.DataFrame, output_dir: Path) -> list[Path]:
    _setup_theme()
    data = scores[scores["tool"].str.startswith("sword2") & scores["variant"].isin(["optimal", "alternative"])].copy()
    if data.empty or "ndo" not in data.columns:
        return []
    data["rank"] = data["partition"].map(_partition_rank)
    data = data.dropna(subset=["rank"])
    rows: list[dict[str, object]] = []
    for tool, tool_data in data.groupby("tool"):
        max_rank = int(tool_data["rank"].max())
        for entry_id, group in tool_data.groupby("entry_id"):
            ranked = group.groupby("rank", as_index=False)["ndo"].max()
            best_so_far: float | None = None
            for k in range(1, max_rank + 1):
                at_rank = ranked[ranked["rank"] == k]["ndo"]
                if not at_rank.empty:
                    rank_best = float(at_rank.max())
                    best_so_far = rank_best if best_so_far is None else max(best_so_far, rank_best)
                if best_so_far is not None:
                    rows.append({"tool": tool, "entry_id": entry_id, "k": k, "best_ndo": best_so_far})
    if not rows:
        return []
    curve = pd.DataFrame(rows).groupby(["tool", "k"], as_index=False)["best_ndo"].mean()
    plt.figure(figsize=(9, 5.5))
    sns.lineplot(data=curve, x="k", y="best_ndo", hue="tool", marker="o", palette=TOOL_PALETTE)
    plt.xlabel("Top-k SWORD2 partitionings considered (entries carried forward)")
    plt.ylabel("Mean best NDO")
    plt.title("SWORD2 top-k alternative value (fixed entry cohort)")
    plt.legend(title="Tool", frameon=True)
    return _savefig(output_dir / "sword2_topk_ndo.png")


def plot_alternatives_per_entry(scores: pd.DataFrame, output_dir: Path) -> list[Path]:
    _setup_theme()
    data = scores[scores["tool"].str.startswith("sword2") & scores["variant"].isin(["optimal", "alternative"])]
    if data.empty:
        return []
    counts = data.groupby(["entry_id", "tool"], as_index=False).size()
    counts["n_alternatives"] = counts["size"] - 1
    plt.figure(figsize=(8, 5))
    order = _ordered_tools(counts)
    _violinplot_with_points(
        data=counts,
        x="tool",
        y="n_alternatives",
        hue="tool",
        order=order,
        hue_order=order,
        palette=TOOL_PALETTE,
        legend=False,
    )
    plt.xlabel("")
    plt.ylabel("Alternative partitionings per entry")
    plt.title("SWORD2 alternative count distribution")
    plt.xticks(rotation=20, ha="right")
    return _savefig(output_dir / "sword2_alternatives_per_entry.png")


def plot_paired_deltas(scores: pd.DataFrame, output_dir: Path, baseline: str = "merizo") -> list[Path]:
    _setup_theme()
    data = headline_scores(scores)
    if baseline not in set(data["tool"]):
        return []
    metrics = [metric for metric in ["ndo", "boundary_f1_10", "matched_dice", "pairwise_f1"] if metric in data.columns]
    rows: list[dict[str, object]] = []
    for metric in metrics:
        pivot = data.pivot_table(index="entry_id", columns="tool", values=metric, aggfunc="first")
        if baseline not in pivot:
            continue
        for tool in [tool for tool in pivot.columns if tool != baseline]:
            paired = pivot[[baseline, tool]].dropna()
            for entry_id, row in paired.iterrows():
                rows.append(
                    {
                        "entry_id": entry_id,
                        "tool": tool,
                        "metric": METRIC_LABELS.get(metric, metric),
                        "delta": row[tool] - row[baseline],
                    }
                )
    if not rows:
        return []
    plot_data = pd.DataFrame(rows)
    plt.figure(figsize=(10, 5.5))
    _violinplot_with_points(
        data=plot_data,
        x="metric",
        y="delta",
        hue="tool",
        hue_order=_ordered_tools(plot_data),
        palette=TOOL_PALETTE,
        dodge_points=True,
        point_alpha=0.12,
        point_size=2,
    )
    plt.axhline(0, color="#333333", linewidth=1)
    plt.xlabel("")
    plt.ylabel(f"Delta vs {baseline}")
    plt.title(f"Paired headline metric deltas vs {baseline}")
    plt.xticks(rotation=20, ha="right")
    plt.legend(title="Tool", frameon=True, loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0)
    return _savefig(output_dir / f"paired_deltas_vs_{baseline}.png")
