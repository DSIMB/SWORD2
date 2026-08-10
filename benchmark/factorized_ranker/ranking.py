"""Reference symmetrized Borda ranking and exact tie rules."""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from numbers import Integral, Real

import numpy as np

from benchmark.factorized_ranker.pairs import pair_vector


def _probability(value: object, description: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{description} must be a real probability")
    probability = float(value)
    if not math.isfinite(probability) or not 0.0 <= probability <= 1.0:
        raise ValueError(f"{description} must be finite and in [0,1]")
    return probability


def symmetrized_probability(
    predict_probability: Callable[[np.ndarray], float],
    shared: Sequence[float],
    left: Sequence[float],
    right: Sequence[float],
) -> float:
    """Average both directional predictions into P(left beats right)."""
    if not callable(predict_probability):
        raise ValueError("predict_probability must be callable")
    left_right = pair_vector(shared, left, right)
    right_left = pair_vector(shared, right, left)
    p_left_right = _probability(
        predict_probability(left_right), "left-to-right prediction"
    )
    p_right_left = _probability(
        predict_probability(right_left), "right-to-left prediction"
    )
    result = 0.5 * (p_left_right + 1.0 - p_right_left)
    return _probability(result, "symmetrized prediction")


def _canonical_items(
    item_ids: Sequence[int] | Sequence[str],
) -> tuple[int, ...] | tuple[str, ...]:
    values = list(item_ids)
    if not values:
        raise ValueError("Borda item IDs are empty")
    if all(isinstance(value, Integral) and not isinstance(value, (bool, np.bool_)) for value in values):
        normalized: tuple[int, ...] | tuple[str, ...] = tuple(int(value) for value in values)
        if any(value <= 0 for value in normalized):
            raise ValueError("count item IDs must be positive")
    elif all(isinstance(value, str) for value in values):
        strings = tuple(values)
        if any(not value or "\0" in value for value in strings):
            raise ValueError("candidate item IDs must be nonempty strings without NUL")
        normalized = strings
    else:
        raise ValueError("Borda item IDs must have one homogeneous supported type")
    if len(set(normalized)) != len(normalized):
        raise ValueError("Borda item IDs contain duplicates")
    return tuple(sorted(normalized))


def normalized_borda(
    item_ids: Sequence[int] | Sequence[str],
    probability: Callable[[int | str, int | str], float],
) -> dict[int | str, float]:
    """Compute canonical arithmetic-mean pairwise win scores."""
    items = _canonical_items(item_ids)
    if not callable(probability):
        raise ValueError("probability must be callable")
    if len(items) == 1:
        return {items[0]: 1.0}

    cached: dict[tuple[int | str, int | str], float] = {}
    for left_index, left in enumerate(items):
        for right in items[left_index + 1 :]:
            p_left_right = _probability(
                probability(left, right), "left-to-right Borda probability"
            )
            p_right_left = _probability(
                probability(right, left), "right-to-left Borda probability"
            )
            cached[(left, right)] = _probability(
                0.5 * (p_left_right + 1.0 - p_right_left),
                "symmetrized Borda probability",
            )

    scores: dict[int | str, float] = {}
    denominator = len(items) - 1
    for item in items:
        total = 0.0
        for opponent in items:
            if opponent == item:
                continue
            if item < opponent:
                win_probability = cached[(item, opponent)]
            else:
                win_probability = 1.0 - cached[(opponent, item)]
            total += win_probability
        scores[item] = _probability(total / denominator, "normalized Borda score")
    return scores


def _mapping_items(scores: Mapping[object, object], description: str) -> list[tuple[object, object]]:
    if not isinstance(scores, Mapping):
        raise ValueError(f"{description} scores must be a mapping")
    try:
        items = list(scores.items())
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError(f"{description} scores cannot be materialized") from error
    if not items:
        raise ValueError(f"{description} scores are empty")
    return items


def _finite_score(value: object, description: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{description} score must be real")
    score = float(value)
    if not math.isfinite(score):
        raise ValueError(f"{description} score must be finite")
    return score


def select_count(scores: Mapping[int, float], legacy_count: int) -> int:
    normalized: dict[int, float] = {}
    for raw_count, raw_score in _mapping_items(scores, "count"):
        if (
            not isinstance(raw_count, Integral)
            or isinstance(raw_count, (bool, np.bool_))
            or int(raw_count) <= 0
        ):
            raise ValueError("count score keys must be positive integers")
        count = int(raw_count)
        if count in normalized:
            raise ValueError("count scores contain duplicate keys")
        normalized[count] = _finite_score(raw_score, "count")
    maximum = max(normalized.values())
    tied = sorted(count for count, score in normalized.items() if score == maximum)
    legacy_valid = (
        isinstance(legacy_count, Integral)
        and not isinstance(legacy_count, (bool, np.bool_))
        and int(legacy_count) > 0
    )
    if legacy_valid and int(legacy_count) in tied:
        return int(legacy_count)
    return tied[0]


def select_candidate(scores: Mapping[str, float]) -> str:
    normalized: dict[str, float] = {}
    for raw_candidate, raw_score in _mapping_items(scores, "candidate"):
        if not isinstance(raw_candidate, str) or not raw_candidate or "\0" in raw_candidate:
            raise ValueError("candidate score keys must be nonempty strings without NUL")
        if raw_candidate in normalized:
            raise ValueError("candidate scores contain duplicate keys")
        normalized[raw_candidate] = _finite_score(raw_score, "candidate")
    maximum = max(normalized.values())
    return min(candidate for candidate, score in normalized.items() if score == maximum)
