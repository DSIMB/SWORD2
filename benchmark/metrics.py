from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.optimize import linear_sum_assignment

from benchmark.numbering import _split_segment, infer_n_res, map_zero_based_chopping, split_domains
from benchmark.scorers.domain_boundary_distance_score import (
    boundary_distance_score,
    get_true_boundary_res,
)
from benchmark.scorers.ndo_score import ndo_score


@dataclass(frozen=True)
class MetricResult:
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
    n_true_domains: int
    n_pred_domains: int


def _n_domains(domain_dict: dict[str, list[int]]) -> int:
    return len([name for name in domain_dict if name != "linker"])


def _make_domain_dict(domain_str: str, n_res: int) -> dict[str, list[int]]:
    domain_dict: dict[str, list[int]] = {}
    assigned: set[int] = set()
    for i, domain in enumerate(split_domains(domain_str), start=1):
        residues: list[int] = []
        for segment in domain:
            start, end = _split_segment(segment)
            lo, hi = sorted((int(start), int(end)))
            for residue in range(lo, hi + 1):
                if n_res > 0 and not 0 <= residue < n_res:
                    continue
                if residue in assigned:
                    continue
                residues.append(residue)
                assigned.add(residue)
        if residues:
            domain_dict[f"D{i}"] = residues
    domain_dict["linker"] = list(set(range(n_res)).difference(assigned)) if n_res > 0 else []
    return domain_dict


def _domain_iou(true: dict[str, list[int]], pred: dict[str, list[int]]) -> float:
    true_domains = [(name, set(res)) for name, res in true.items() if name != "linker"]
    pred_domains = [(name, set(res)) for name, res in pred.items() if name != "linker"]
    total = sum(len(res) for _, res in true_domains)
    if total == 0:
        return 0.0

    weighted = 0.0
    for _, true_res in true_domains:
        best = 0.0
        for _, pred_res in pred_domains:
            union = len(true_res | pred_res)
            if union:
                best = max(best, len(true_res & pred_res) / union)
        weighted += len(true_res) * best
    return weighted / total


def _comb2(value: int) -> float:
    return float(value * (value - 1) / 2)


def _domain_sets(domain_dict: dict[str, list[int]], include_linker: bool = False) -> list[set[int]]:
    return [
        set(residues)
        for name, residues in sorted(domain_dict.items())
        if include_linker or name != "linker"
    ]


def _assignment_labels(domain_dict: dict[str, list[int]], n_res: int) -> list[int]:
    labels = [-1] * n_res
    for label, (name, residues) in enumerate(
        sorted((name, residues) for name, residues in domain_dict.items() if name != "linker")
    ):
        for residue in residues:
            if 0 <= residue < n_res:
                labels[residue] = label
    return labels


def _boundary_positions(domain_dict: dict[str, list[int]], n_res: int) -> set[int]:
    labels = _assignment_labels(domain_dict, n_res)
    return {
        idx
        for idx in range(n_res - 1)
        if labels[idx] >= 0 and labels[idx + 1] >= 0 and labels[idx] != labels[idx + 1]
    }


def _match_boundaries(
    true_boundaries: set[int],
    pred_boundaries: set[int],
    tolerance: int,
) -> tuple[int, int, int]:
    unmatched_true = set(true_boundaries)
    matches = 0
    for pred_boundary in sorted(pred_boundaries):
        candidates = [
            true_boundary
            for true_boundary in unmatched_true
            if abs(true_boundary - pred_boundary) <= tolerance
        ]
        if not candidates:
            continue
        nearest = min(candidates, key=lambda true_boundary: abs(true_boundary - pred_boundary))
        unmatched_true.remove(nearest)
        matches += 1
    return matches, len(pred_boundaries), len(true_boundaries)


def _precision_recall_f1(matches: int, predicted: int, actual: int) -> tuple[float, float, float]:
    precision = 1.0 if predicted == 0 and actual == 0 else (matches / predicted if predicted else 0.0)
    recall = 1.0 if predicted == 0 and actual == 0 else (matches / actual if actual else 0.0)
    f1 = 0.0 if precision + recall == 0.0 else 2 * precision * recall / (precision + recall)
    return float(precision), float(recall), float(f1)


def _boundary_metrics(
    true: dict[str, list[int]],
    pred: dict[str, list[int]],
    n_res: int,
) -> dict[str, float]:
    true_boundaries = _boundary_positions(true, n_res)
    pred_boundaries = _boundary_positions(pred, n_res)
    metrics: dict[str, float] = {}
    for tolerance in (5, 10, 20):
        matches, predicted, actual = _match_boundaries(true_boundaries, pred_boundaries, tolerance)
        precision, recall, f1 = _precision_recall_f1(matches, predicted, actual)
        metrics[f"boundary_precision_{tolerance}"] = precision
        metrics[f"boundary_recall_{tolerance}"] = recall
        metrics[f"boundary_f1_{tolerance}"] = f1

    if not true_boundaries and not pred_boundaries:
        errors = [0.0]
    elif not true_boundaries or not pred_boundaries:
        errors = [float(n_res)]
    else:
        errors = [
            float(min(abs(true_boundary - pred_boundary) for pred_boundary in pred_boundaries))
            for true_boundary in true_boundaries
        ]
    metrics["median_boundary_error"] = float(np.median(errors))
    metrics["worst_boundary_error"] = float(max(errors))
    return metrics


def _coverage(domain_dict: dict[str, list[int]], n_res: int) -> tuple[float, float]:
    if n_res <= 0:
        return 0.0, 0.0
    covered = set().union(*_domain_sets(domain_dict)) if _domain_sets(domain_dict) else set()
    coverage = len(covered) / n_res
    return float(coverage), float(1.0 - coverage)


def _pairwise_f1(true: dict[str, list[int]], pred: dict[str, list[int]]) -> tuple[float, float, float]:
    true_domains = _domain_sets(true)
    pred_domains = _domain_sets(pred)
    true_pairs = sum(_comb2(len(domain)) for domain in true_domains)
    pred_pairs = sum(_comb2(len(domain)) for domain in pred_domains)
    tp = 0.0
    for true_domain in true_domains:
        for pred_domain in pred_domains:
            tp += _comb2(len(true_domain & pred_domain))
    return _precision_recall_f1(int(tp), int(pred_pairs), int(true_pairs))


def _information_metrics(
    true: dict[str, list[int]],
    pred: dict[str, list[int]],
    n_res: int,
) -> tuple[float, float, float]:
    if n_res <= 0:
        return 0.0, 0.0, 0.0
    true_sets = _domain_sets(true, include_linker=True)
    pred_sets = _domain_sets(pred, include_linker=True)
    true_sizes = [len(domain) for domain in true_sets]
    pred_sizes = [len(domain) for domain in pred_sets]

    h_true = -sum((size / n_res) * math.log(size / n_res) for size in true_sizes if size)
    h_pred = -sum((size / n_res) * math.log(size / n_res) for size in pred_sizes if size)
    mutual_info = 0.0
    for true_domain in true_sets:
        for pred_domain in pred_sets:
            overlap = len(true_domain & pred_domain)
            if overlap:
                mutual_info += (overlap / n_res) * math.log((overlap * n_res) / (len(true_domain) * len(pred_domain)))
    normalized_mi = 1.0 if h_true == 0.0 and h_pred == 0.0 else (
        mutual_info / math.sqrt(h_true * h_pred) if h_true and h_pred else 0.0
    )
    variation = h_true + h_pred - 2 * mutual_info
    return float(mutual_info), float(normalized_mi), float(max(0.0, variation))


def _adjusted_rand(true: dict[str, list[int]], pred: dict[str, list[int]], n_res: int) -> float:
    if n_res <= 1:
        return 1.0
    true_sets = _domain_sets(true, include_linker=True)
    pred_sets = _domain_sets(pred, include_linker=True)
    sum_comb = sum(_comb2(len(true_domain & pred_domain)) for true_domain in true_sets for pred_domain in pred_sets)
    row_comb = sum(_comb2(len(domain)) for domain in true_sets)
    col_comb = sum(_comb2(len(domain)) for domain in pred_sets)
    total_comb = _comb2(n_res)
    expected = row_comb * col_comb / total_comb if total_comb else 0.0
    max_index = 0.5 * (row_comb + col_comb)
    denominator = max_index - expected
    return 1.0 if denominator == 0.0 else float((sum_comb - expected) / denominator)


def _matched_domain_scores(true: dict[str, list[int]], pred: dict[str, list[int]]) -> tuple[float, float]:
    true_domains = _domain_sets(true)
    pred_domains = _domain_sets(pred)
    total_true = sum(len(domain) for domain in true_domains)
    if not true_domains or not pred_domains or total_true == 0:
        return 0.0, 0.0

    dice = np.zeros((len(true_domains), len(pred_domains)))
    jaccard = np.zeros_like(dice)
    for i, true_domain in enumerate(true_domains):
        for j, pred_domain in enumerate(pred_domains):
            intersection = len(true_domain & pred_domain)
            dice[i, j] = 2 * intersection / (len(true_domain) + len(pred_domain))
            union = len(true_domain | pred_domain)
            jaccard[i, j] = intersection / union if union else 0.0

    rows, cols = linear_sum_assignment(-dice)
    weighted_dice = 0.0
    weighted_jaccard = 0.0
    for row, col in zip(rows, cols):
        weight = len(true_domains[row]) / total_true
        weighted_dice += weight * dice[row, col]
        weighted_jaccard += weight * jaccard[row, col]
    return float(weighted_dice), float(weighted_jaccard)


def _exact_match(true: dict[str, list[int]], pred: dict[str, list[int]]) -> float:
    true_signature = {frozenset(domain) for domain in _domain_sets(true)}
    pred_signature = {frozenset(domain) for domain in _domain_sets(pred)}
    return 1.0 if true_signature == pred_signature and set(true.get("linker", [])) == set(pred.get("linker", [])) else 0.0


def score_domains(true: dict[str, list[int]], pred: dict[str, list[int]]) -> MetricResult:
    n_true = _n_domains(true)
    n_pred = _n_domains(pred)
    ndo = float(ndo_score(true, pred))
    boundaries = get_true_boundary_res(true)
    boundaries["boundary_res"] = np.asarray(boundaries["boundary_res"])
    boundary = float(boundary_distance_score(pred, boundaries))
    n_res = len(set().union(*_domain_sets(true, include_linker=True)))
    boundary_extra = _boundary_metrics(true, pred, n_res)
    pred_coverage, pred_linker_fraction = _coverage(pred, n_res)
    pairwise_precision, pairwise_recall, pairwise_f1 = _pairwise_f1(true, pred)
    _, normalized_mutual_info, variation_of_information = _information_metrics(true, pred, n_res)
    matched_dice, matched_jaccard = _matched_domain_scores(true, pred)
    return MetricResult(
        ndo=ndo,
        boundary_dist_score=boundary,
        d_count_acc=1.0 if n_true == n_pred else 0.0,
        d_count_dev=float(abs(n_true - n_pred)),
        iou=float(_domain_iou(true, pred)),
        multi_ndo=ndo if n_true > 1 else float("nan"),
        domain_count_bias=float(n_pred - n_true),
        over_split=1.0 if n_pred > n_true else 0.0,
        merge=1.0 if n_pred < n_true else 0.0,
        boundary_precision_5=boundary_extra["boundary_precision_5"],
        boundary_recall_5=boundary_extra["boundary_recall_5"],
        boundary_f1_5=boundary_extra["boundary_f1_5"],
        boundary_precision_10=boundary_extra["boundary_precision_10"],
        boundary_recall_10=boundary_extra["boundary_recall_10"],
        boundary_f1_10=boundary_extra["boundary_f1_10"],
        boundary_precision_20=boundary_extra["boundary_precision_20"],
        boundary_recall_20=boundary_extra["boundary_recall_20"],
        boundary_f1_20=boundary_extra["boundary_f1_20"],
        median_boundary_error=boundary_extra["median_boundary_error"],
        worst_boundary_error=boundary_extra["worst_boundary_error"],
        pred_coverage=pred_coverage,
        pred_linker_fraction=pred_linker_fraction,
        pairwise_precision=pairwise_precision,
        pairwise_recall=pairwise_recall,
        pairwise_f1=pairwise_f1,
        adjusted_rand=_adjusted_rand(true, pred, n_res),
        normalized_mutual_info=normalized_mutual_info,
        variation_of_information=variation_of_information,
        matched_dice=matched_dice,
        matched_jaccard=matched_jaccard,
        exact_match=_exact_match(true, pred),
        n_true_domains=n_true,
        n_pred_domains=n_pred,
    )


def score_choppings(
    true_chopping: str,
    pred_chopping: str,
    n_res: int | None = None,
) -> MetricResult:
    true_norm = map_zero_based_chopping(true_chopping)
    pred_norm = map_zero_based_chopping(pred_chopping)
    if n_res is None:
        # Chainsaw's published benchmark scorer includes unassigned truth residues
        # as linker/NDRs, but does not synthesize linker residues in predictions.
        true = _make_domain_dict(true_norm, infer_n_res([true_norm, pred_norm]))
        pred = _make_domain_dict(pred_norm, 0)
        return score_domains(true, pred)
    true = _make_domain_dict(true_norm, n_res)
    pred = _make_domain_dict(pred_norm, n_res)
    return score_domains(true, pred)


def domain_count(chopping: str) -> int:
    return len(split_domains(chopping))


def chopping_from_named_bounds(bounds: str, names: str) -> str:
    """Group published Chainsaw boundary segments by their domain names."""
    grouped: dict[str, list[str]] = {}
    order: list[str] = []
    bound_parts = [part for part in bounds.split("|") if part]
    name_parts = [part for part in names.split("|") if part]
    if len(bound_parts) != len(name_parts):
        raise ValueError(f"Bounds/name length mismatch: {bounds!r} vs {names!r}")
    for bound, name in zip(bound_parts, name_parts):
        if name not in grouped:
            grouped[name] = []
            order.append(name)
        grouped[name].append(bound)
    return ",".join("_".join(grouped[name]) for name in order)
