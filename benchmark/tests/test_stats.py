from __future__ import annotations

from collections.abc import Iterator, Mapping

import numpy as np
import pandas as pd
import pytest

import benchmark.stats as stats
from benchmark.stats import paired_chain_bootstrap, paired_wilcoxon


def test_paired_bootstrap_uses_only_common_chain_ids() -> None:
    result = paired_chain_bootstrap(
        {"a": 0.7, "b": 0.9, "only_base": 1.0},
        {"a": 0.8, "b": 1.0, "only_new": 0.0},
        n_resamples=1000,
        seed=37,
    )
    assert result.n == 2
    assert result.mean_delta == pytest.approx(0.1)
    assert all(type(value) is float and np.isfinite(value) for value in (
        result.mean_delta,
        result.ci_low,
        result.ci_high,
    ))


def test_paired_bootstrap_is_deterministic_and_order_independent() -> None:
    first = paired_chain_bootstrap({"b": 0.2, "a": 0.1}, {"a": 0.3, "b": 0.4})
    second = paired_chain_bootstrap({"a": 0.1, "b": 0.2}, {"b": 0.4, "a": 0.3})
    assert first == second


def test_lower_is_better_reverses_delta_exactly() -> None:
    higher = paired_chain_bootstrap({"a": 1.0}, {"a": 3.0}, higher_is_better=True)
    lower = paired_chain_bootstrap({"a": 1.0}, {"a": 3.0}, higher_is_better=False)
    assert lower.mean_delta == -higher.mean_delta
    assert lower.ci_low == -higher.ci_high
    assert lower.ci_high == -higher.ci_low


def test_single_chain_has_point_interval() -> None:
    result = paired_chain_bootstrap({"a": 0.1}, {"a": 0.4}, n_resamples=17)
    assert result.n == 1
    assert result.mean_delta == pytest.approx(0.3)
    assert result.ci_low == result.ci_high == result.mean_delta


def test_bootstrap_chunk_sizes_consume_one_identical_rng_stream() -> None:
    delta = np.asarray([0.1, -0.3, 0.7], dtype=np.float64)
    expected = stats._bootstrap_means(delta, 1001, np.random.default_rng(37), chunk_size=1)
    for chunk_size in (17, 256, 2000):
        actual = stats._bootstrap_means(
            delta,
            1001,
            np.random.default_rng(37),
            chunk_size=chunk_size,
        )
        assert np.array_equal(actual, expected)
    low, high = np.quantile(expected, [0.025, 0.975], method="linear")
    result = paired_chain_bootstrap(
        {"a": 0.0, "b": 0.0, "c": 0.0},
        {"a": 0.1, "b": -0.3, "c": 0.7},
        n_resamples=1001,
        seed=37,
    )
    assert (result.ci_low, result.ci_high) == (float(low), float(high))


def test_private_bootstrap_helper_never_requests_more_than_chunk_size() -> None:
    class RecordingRng:
        def __init__(self) -> None:
            self.requests: list[tuple[int, int]] = []

        def integers(self, _low, high, *, size, dtype):
            assert high == 2 and dtype == np.int64
            self.requests.append(size)
            return np.zeros(size, dtype=dtype)

    rng = RecordingRng()
    means = stats._bootstrap_means(
        np.asarray([1.0, 2.0]), 600, rng, chunk_size=stats.BOOTSTRAP_CHUNK_SIZE
    )
    assert means.shape == (600,)
    assert [rows for rows, _columns in rng.requests] == [256, 256, 88]
    assert all(columns == 2 for _rows, columns in rng.requests)


class DuplicateItemsMapping(Mapping[str, float]):
    def __getitem__(self, key: str) -> float:
        return 1.0

    def __iter__(self) -> Iterator[str]:
        return iter(("a",))

    def __len__(self) -> int:
        return 1

    def items(self):
        return [("a", 1.0), ("a", 2.0)]


@pytest.mark.parametrize(
    "baseline,experiment,kwargs",
    [
        ({"a": 1.0}, {"b": 2.0}, {}),
        ({"": 1.0}, {"": 2.0}, {}),
        ({1: 1.0}, {1: 2.0}, {}),
        ({"a": True}, {"a": 2.0}, {}),
        ({"a": float("nan")}, {"a": 2.0}, {}),
        ({"a": 1.0}, {"a": float("inf")}, {}),
        ({"a": 1.0}, {"a": 2.0}, {"n_resamples": 0}),
        ({"a": 1.0}, {"a": 2.0}, {"n_resamples": True}),
        ({"a": 1.0}, {"a": 2.0}, {"n_resamples": 1.5}),
        ({"a": 1.0}, {"a": 2.0}, {"seed": -1}),
        ({"a": 1.0}, {"a": 2.0}, {"seed": True}),
        ({"a": 1.0}, {"a": 2.0}, {"seed": 2**64}),
        ({"a": 1.0}, {"a": 2.0}, {"higher_is_better": 1}),
    ],
)
def test_paired_bootstrap_rejects_invalid_inputs(baseline, experiment, kwargs) -> None:
    with pytest.raises(ValueError):
        paired_chain_bootstrap(baseline, experiment, **kwargs)


def test_paired_bootstrap_rejects_duplicate_mapping_items() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        paired_chain_bootstrap(DuplicateItemsMapping(), {"a": 2.0})


def test_existing_paired_wilcoxon_contract_is_unchanged() -> None:
    scores = pd.DataFrame(
        [
            {"entry_id": "a", "tool": "base", "metric": 0.5},
            {"entry_id": "a", "tool": "new", "metric": 0.7},
            {"entry_id": "b", "tool": "base", "metric": 0.6},
            {"entry_id": "b", "tool": "new", "metric": 0.8},
        ]
    )
    result = paired_wilcoxon(scores, "base", "new", "metric")
    assert result.n_pairs == 2
    assert result.median_delta == pytest.approx(-0.2)
