from __future__ import annotations

import numpy as np
import pytest

from benchmark.factorized_ranker.ranking import (
    normalized_borda,
    select_candidate,
    select_count,
    symmetrized_probability,
)


def test_symmetrization_calls_both_orientations_once() -> None:
    calls: list[np.ndarray] = []

    def predictor(vector: np.ndarray) -> float:
        calls.append(vector.copy())
        return 0.8 if len(calls) == 1 else 0.3

    probability = symmetrized_probability(predictor, [10.0], [3.0], [1.0])
    assert len(calls) == 2
    assert probability == pytest.approx(0.5 * (0.8 + 1.0 - 0.3))
    np.testing.assert_array_equal(calls[0], [10.0, 2.0, 2.0])
    np.testing.assert_array_equal(calls[1], [10.0, -2.0, 2.0])


def test_ties_follow_approved_rules() -> None:
    assert select_count({2: 0.5, 3: 0.5}, legacy_count=3) == 3
    assert select_count({2: 0.5, 3: 0.5}, legacy_count=4) == 2
    assert select_candidate({"0-4 5-9": 0.5, "0-5 6-9": 0.5}) == "0-4 5-9"


def test_normalized_borda_calls_two_orientations_once_and_is_canonical() -> None:
    calls: list[tuple[str, str]] = []
    directional = {
        ("a", "b"): 0.8,
        ("b", "a"): 0.4,
        ("a", "c"): 0.2,
        ("c", "a"): 0.6,
        ("b", "c"): 0.7,
        ("c", "b"): 0.1,
    }

    def probability(left: str, right: str) -> float:
        calls.append((left, right))
        return directional[(left, right)]

    scores = normalized_borda(["c", "a", "b"], probability)
    assert list(scores) == ["a", "b", "c"]
    assert calls == [
        ("a", "b"),
        ("b", "a"),
        ("a", "c"),
        ("c", "a"),
        ("b", "c"),
        ("c", "b"),
    ]
    assert scores == {
        "a": pytest.approx((0.7 + 0.3) / 2),
        "b": pytest.approx((0.3 + 0.8) / 2),
        "c": pytest.approx((0.7 + 0.2) / 2),
    }


def test_borda_input_permutation_is_bit_identical() -> None:
    def probability(left: int, right: int) -> float:
        return 0.75 if left < right else 0.2

    first = normalized_borda([3, 1, 2], probability)
    second = normalized_borda([2, 3, 1], probability)
    assert first == second
    assert list(first) == [1, 2, 3]


def test_one_item_avoids_callback_and_empty_fails() -> None:
    def forbidden(*_args) -> float:
        raise AssertionError("callback must not run")

    assert normalized_borda([2], forbidden) == {2: 1.0}
    with pytest.raises(ValueError):
        normalized_borda([], forbidden)


@pytest.mark.parametrize(
    "items",
    [[1, "2"], [True, 2], [0, 2], ["", "a"], ["a", "a"], [1, 1]],
)
def test_borda_rejects_mixed_invalid_or_duplicate_ids(items) -> None:
    with pytest.raises(ValueError):
        normalized_borda(items, lambda _left, _right: 0.5)


@pytest.mark.parametrize("value", [np.nan, np.inf, -0.1, 1.1, "0.5", True])
def test_directional_probability_validation_is_strict(value: object) -> None:
    with pytest.raises(ValueError):
        symmetrized_probability(lambda _vector: value, [1.0], [2.0], [3.0])
    with pytest.raises(ValueError):
        normalized_borda([1, 2], lambda _left, _right: value)


def test_hand_calculated_two_item_borda() -> None:
    calls = {(1, 2): 0.9, (2, 1): 0.3}
    assert normalized_borda([2, 1], lambda left, right: calls[(left, right)]) == {
        1: pytest.approx(0.8),
        2: pytest.approx(0.2),
    }


def test_selectors_validate_maps_and_do_not_round_near_ties() -> None:
    assert select_count({2: 0.5000000000000001, 3: 0.5}, legacy_count=3) == 2
    assert select_candidate({"b": 0.5, "a": 0.5000000000000001}) == "a"
    for scores in ({}, {0: 0.5}, {True: 0.5}, {2: np.nan}, {2: "0.5"}):
        with pytest.raises(ValueError):
            select_count(scores, legacy_count=2)
    for scores in ({}, {"": 0.5}, {"a": np.inf}, {"a": False}):
        with pytest.raises(ValueError):
            select_candidate(scores)


def test_invalid_legacy_count_is_not_preferred_in_a_tie() -> None:
    assert select_count({2: 0.5, 3: 0.5}, legacy_count=True) == 2
    assert select_count({2: 0.5, 3: 0.5}, legacy_count=0) == 2
