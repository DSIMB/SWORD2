from __future__ import annotations

import csv
import math
import statistics
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from benchmark.datasets import CathEntry, strip_cath_labels
from benchmark.metrics import score_choppings
from benchmark.numbering import StructureNumbering, map_author_chopping
from benchmark.runners.base import PartitionPrediction


@dataclass(frozen=True)
class ScoreRow:
    dataset: str
    entry_id: str
    pdb_id: str
    chain_id: str
    tool: str
    variant: str
    partition: str
    true_chopping: str
    pred_chopping: str
    n_residues: int
    n_true_domains: int
    n_pred_domains: int
    ndo: float
    boundary_dist_score: float
    d_count_acc: float
    d_count_dev: float
    iou: float
    multi_ndo: float
    domain_count_bias: float
    over_split: float
    merge: float
    boundary_precision_5: float
    boundary_recall_5: float
    boundary_f1_5: float
    boundary_precision_10: float
    boundary_recall_10: float
    boundary_f1_10: float
    boundary_precision_20: float
    boundary_recall_20: float
    boundary_f1_20: float
    median_boundary_error: float
    worst_boundary_error: float
    pred_coverage: float
    pred_linker_fraction: float
    pairwise_precision: float
    pairwise_recall: float
    pairwise_f1: float
    adjusted_rand: float
    normalized_mutual_info: float
    variation_of_information: float
    matched_dice: float
    matched_jaccard: float
    exact_match: float
    runtime_s: float | None = None
    peak_rss_mb: float | None = None


@dataclass(frozen=True)
class RunRow:
    dataset: str
    entry_id: str
    pdb_id: str
    chain_id: str
    tool: str
    variant: str
    reused_output: bool
    returncode: int | None
    runtime_s: float | None
    peak_rss_mb: float | None
    command: str
    cwd: str
    raw_path: str
    batch_runtime_s: float | None = None
    batch_peak_rss_mb: float | None = None
    batch_n_entries: int | None = None


@dataclass(frozen=True)
class FailureRow:
    dataset: str
    entry_id: str
    pdb_id: str
    chain_id: str
    tool: str
    stage: str
    message: str
    returncode: int | None = None
    command: str | None = None
    stderr: str | None = None


def score_prediction(
    entry: CathEntry,
    numbering: StructureNumbering,
    prediction: PartitionPrediction,
    runtime_s: float | None = None,
    peak_rss_mb: float | None = None,
) -> ScoreRow:
    true_chopping = map_author_chopping(strip_cath_labels(entry.chopping), numbering, entry.chain_id)
    metrics = score_choppings(true_chopping, prediction.chopping, n_res=numbering.n_residues)
    return ScoreRow(
        dataset=entry.dataset,
        entry_id=entry.entry_id,
        pdb_id=entry.pdb_id,
        chain_id=entry.chain_id,
        tool=prediction.tool,
        variant=prediction.variant,
        partition=prediction.name,
        true_chopping=true_chopping,
        pred_chopping=prediction.chopping,
        n_residues=numbering.n_residues,
        n_true_domains=metrics.n_true_domains,
        n_pred_domains=metrics.n_pred_domains,
        ndo=metrics.ndo,
        boundary_dist_score=metrics.boundary_dist_score,
        d_count_acc=metrics.d_count_acc,
        d_count_dev=metrics.d_count_dev,
        iou=metrics.iou,
        multi_ndo=metrics.multi_ndo,
        domain_count_bias=metrics.domain_count_bias,
        over_split=metrics.over_split,
        merge=metrics.merge,
        boundary_precision_5=metrics.boundary_precision_5,
        boundary_recall_5=metrics.boundary_recall_5,
        boundary_f1_5=metrics.boundary_f1_5,
        boundary_precision_10=metrics.boundary_precision_10,
        boundary_recall_10=metrics.boundary_recall_10,
        boundary_f1_10=metrics.boundary_f1_10,
        boundary_precision_20=metrics.boundary_precision_20,
        boundary_recall_20=metrics.boundary_recall_20,
        boundary_f1_20=metrics.boundary_f1_20,
        median_boundary_error=metrics.median_boundary_error,
        worst_boundary_error=metrics.worst_boundary_error,
        pred_coverage=metrics.pred_coverage,
        pred_linker_fraction=metrics.pred_linker_fraction,
        pairwise_precision=metrics.pairwise_precision,
        pairwise_recall=metrics.pairwise_recall,
        pairwise_f1=metrics.pairwise_f1,
        adjusted_rand=metrics.adjusted_rand,
        normalized_mutual_info=metrics.normalized_mutual_info,
        variation_of_information=metrics.variation_of_information,
        matched_dice=metrics.matched_dice,
        matched_jaccard=metrics.matched_jaccard,
        exact_match=metrics.exact_match,
        runtime_s=runtime_s,
        peak_rss_mb=peak_rss_mb,
    )


def add_sword2_oracle(rows: list[ScoreRow]) -> list[ScoreRow]:
    by_entry: dict[tuple[str, str], list[ScoreRow]] = {}
    for row in rows:
        if row.tool.startswith("sword2") and row.variant in {"optimal", "alternative"}:
            by_entry.setdefault((row.entry_id, row.tool), []).append(row)

    oracle_rows: list[ScoreRow] = []
    for candidates in by_entry.values():
        best = max(candidates, key=lambda row: (row.ndo, row.boundary_dist_score, row.iou))
        oracle_rows.append(
            ScoreRow(
                **{
                    **asdict(best),
                    "variant": "oracle",
                    "partition": f"Best alternative by NDO ({best.partition})",
                }
            )
        )
    return [*rows, *oracle_rows]


def _composite_s(row: ScoreRow) -> float:
    """Composite selection score S = mean(NDO, IoU, boundary_f1_10, matched_dice, count_term).

    count_term = 1 / (1 + |n_pred - n_true|) so it is in [0,1] like the other components.
    NaN components are skipped; returns 0.0 when all components are NaN.
    """
    count_term = 1.0 / (1.0 + abs(row.n_pred_domains - row.n_true_domains))
    vals = [row.ndo, row.iou, row.boundary_f1_10, row.matched_dice, count_term]
    finite = [v for v in vals if not math.isnan(v)]
    return statistics.mean(finite) if finite else 0.0


def add_sword2_oracle_s(rows: list[ScoreRow]) -> list[ScoreRow]:
    """Add an oracle_s row for each (entry_id, tool) — the alternative with highest S score."""
    by_entry: dict[tuple[str, str], list[ScoreRow]] = {}
    for row in rows:
        if row.tool.startswith("sword2") and row.variant in {"optimal", "alternative"}:
            by_entry.setdefault((row.entry_id, row.tool), []).append(row)

    oracle_rows: list[ScoreRow] = []
    for candidates in by_entry.values():
        best = max(candidates, key=_composite_s)
        oracle_rows.append(
            ScoreRow(
                **{
                    **asdict(best),
                    "variant": "oracle_s",
                    "partition": f"Best alternative by S ({best.partition})",
                }
            )
        )
    return [*rows, *oracle_rows]


def write_scores_csv(rows: list[ScoreRow], path: Path) -> None:
    if not rows:
        raise ValueError("Cannot write an empty scores CSV")
    _write_dataclass_csv(rows, path, ScoreRow)


def write_runs_csv(rows: list[RunRow], path: Path) -> None:
    _write_dataclass_csv(rows, path, RunRow)


def write_failures_csv(rows: list[FailureRow], path: Path) -> None:
    _write_dataclass_csv(rows, path, FailureRow)


def _write_dataclass_csv(rows: list, path: Path, row_type: type) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        fieldnames = [field.name for field in fields(row_type)]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
