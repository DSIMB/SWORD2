from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Integral, Real
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


BOOTSTRAP_CHUNK_SIZE = 256


@dataclass(frozen=True)
class PairedComparison:
    metric: str
    baseline_tool: str
    comparison_tool: str
    n_pairs: int
    baseline_median: float
    comparison_median: float
    median_delta: float
    p_value: float


@dataclass(frozen=True)
class PairedBootstrap:
    n: int
    mean_delta: float
    ci_low: float
    ci_high: float


def _materialize_metric_mapping(
    values: Mapping[str, float], description: str
) -> dict[str, float]:
    if not isinstance(values, Mapping):
        raise ValueError(f"{description} must be a mapping")
    try:
        items = list(values.items())
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError(f"{description} must expose finite keyed values") from error
    materialized: dict[str, float] = {}
    for key, raw_value in items:
        if not isinstance(key, str) or not key:
            raise ValueError(f"{description} contains a non-string or empty chain ID")
        if key in materialized:
            raise ValueError(f"{description} contains a duplicate chain ID")
        if isinstance(raw_value, (bool, np.bool_)) or not isinstance(raw_value, Real):
            raise ValueError(f"{description} contains a non-real value")
        value = float(raw_value)
        if not np.isfinite(value):
            raise ValueError(f"{description} contains a non-finite value")
        materialized[key] = value
    return materialized


def _bootstrap_means(
    delta: np.ndarray,
    n_resamples: int,
    rng: Any,
    *,
    chunk_size: int = BOOTSTRAP_CHUNK_SIZE,
) -> np.ndarray:
    """Draw bootstrap means in bounded batches from one unchanged RNG stream."""
    if not isinstance(chunk_size, Integral) or isinstance(chunk_size, bool) or chunk_size <= 0:
        raise ValueError("bootstrap chunk size must be a positive integer")
    means = np.empty(n_resamples, dtype=np.float64)
    n_chains = len(delta)
    offset = 0
    while offset < n_resamples:
        batch_size = min(int(chunk_size), n_resamples - offset)
        draws = rng.integers(
            0,
            n_chains,
            size=(batch_size, n_chains),
            dtype=np.int64,
        )
        means[offset : offset + batch_size] = delta[draws].mean(
            axis=1, dtype=np.float64
        )
        offset += batch_size
    if not np.isfinite(means).all():
        raise ValueError("bootstrap produced a non-finite result")
    return means


def paired_chain_bootstrap(
    baseline_by_chain: Mapping[str, float],
    experiment_by_chain: Mapping[str, float],
    *,
    n_resamples: int = 10_000,
    seed: int = 37,
    higher_is_better: bool = True,
) -> PairedBootstrap:
    """Return a deterministic paired chain-level mean-delta interval."""
    baseline = _materialize_metric_mapping(baseline_by_chain, "baseline")
    experiment = _materialize_metric_mapping(experiment_by_chain, "experiment")
    if (
        not isinstance(n_resamples, Integral)
        or isinstance(n_resamples, bool)
        or n_resamples <= 0
    ):
        raise ValueError("n_resamples must be a positive integer")
    if (
        not isinstance(seed, Integral)
        or isinstance(seed, bool)
        or seed < 0
        or seed > 2**64 - 1
    ):
        raise ValueError("seed must be an integer in 0..=2**64-1")
    if type(higher_is_better) is not bool:
        raise ValueError("higher_is_better must be a bool")

    chain_ids = sorted(set(baseline) & set(experiment))
    if not chain_ids:
        raise ValueError("paired bootstrap has no common chain IDs")
    baseline_values = np.asarray([baseline[chain_id] for chain_id in chain_ids], dtype=np.float64)
    experiment_values = np.asarray(
        [experiment[chain_id] for chain_id in chain_ids], dtype=np.float64
    )
    delta = experiment_values - baseline_values
    if not higher_is_better:
        delta = -delta
    if not np.isfinite(delta).all():
        raise ValueError("paired bootstrap deltas are non-finite")

    rng = np.random.default_rng(int(seed))
    means = _bootstrap_means(delta, int(n_resamples), rng)
    low, high = np.quantile(means, [0.025, 0.975], method="linear")
    mean_delta = float(delta.mean(dtype=np.float64))
    result = PairedBootstrap(
        n=int(len(chain_ids)),
        mean_delta=mean_delta,
        ci_low=float(low),
        ci_high=float(high),
    )
    if not all(np.isfinite(value) for value in (result.mean_delta, result.ci_low, result.ci_high)):
        raise ValueError("paired bootstrap summary is non-finite")
    return result


def paired_wilcoxon(
    scores: pd.DataFrame,
    baseline_tool: str,
    comparison_tool: str,
    metric: str,
    variant: str | None = None,
) -> PairedComparison:
    data = scores
    if variant is not None:
        data = data[data["variant"] == variant]
    pivot = data.pivot_table(index="entry_id", columns="tool", values=metric, aggfunc="first")
    paired = pivot[[baseline_tool, comparison_tool]].dropna()
    if paired.empty:
        raise ValueError(f"No paired rows for {baseline_tool} vs {comparison_tool}")
    deltas = paired[baseline_tool] - paired[comparison_tool]
    p_value = float(wilcoxon(deltas).pvalue) if len(deltas) > 1 else 1.0
    return PairedComparison(
        metric=metric,
        baseline_tool=baseline_tool,
        comparison_tool=comparison_tool,
        n_pairs=len(paired),
        baseline_median=float(paired[baseline_tool].median()),
        comparison_median=float(paired[comparison_tool].median()),
        median_delta=float(deltas.median()),
        p_value=p_value,
    )
