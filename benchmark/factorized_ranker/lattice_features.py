"""Deterministic lattice-relative and per-count feature formulas."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence

from benchmark.factorized_ranker.schema import (
    COUNT_ITEM_FEATURES,
    COUNT_SUMMARIES,
    COUNT_SUMMARY_SOURCES,
    RELATIVE_CORE_FEATURES,
    RELATIVE_HIERARCHY_FEATURES,
)


_HIERARCHY_INTEGER_FIELDS = (
    "hierarchy_first_appearance_level",
    "hierarchy_persistence_levels",
    "hierarchy_path_count",
)
_HIERARCHY_FLOAT_FIELDS = (
    "hierarchy_parent_merge_margin",
    "hierarchy_child_merge_margin",
)


def _finite(row: Mapping[str, object], name: str) -> float:
    try:
        value = float(row[name])
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"invalid {name}") from exc
    if not math.isfinite(value):
        raise ValueError(f"non-finite {name}")
    return value


def _integral(row: Mapping[str, object], name: str, *, positive: bool) -> int:
    if isinstance(row.get(name), bool):
        raise ValueError(f"invalid {name}")
    value = _finite(row, name)
    if not value.is_integer() or (value <= 0.0 if positive else value < 0.0):
        raise ValueError(f"invalid {name}")
    return int(value)


def _identity(row: Mapping[str, object]) -> tuple[int, int, str]:
    source_index = _integral(row, "source_index", positive=False)
    num_domains = _integral(row, "num_domains", positive=True)
    canonical = row.get("canonical_delineation")
    if not isinstance(canonical, str) or not canonical:
        raise ValueError("invalid canonical_delineation")
    return source_index, num_domains, canonical


def _unique_identities(
    rows: Sequence[Mapping[str, object]],
) -> list[tuple[int, int, str]]:
    identities = [_identity(row) for row in rows]
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate candidate identity")
    count_canonicals = [(num_domains, canonical) for _, num_domains, canonical in identities]
    if len(set(count_canonicals)) != len(count_canonicals):
        raise ValueError("duplicate candidate canonical delineation for count")
    return identities


def _boundary_set(row: Mapping[str, object]) -> tuple[int, ...]:
    raw = row.get("sequential_boundaries")
    if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
        raise ValueError("invalid sequential_boundaries")
    values: list[int] = []
    for value in raw:
        if isinstance(value, bool):
            raise ValueError("invalid sequential_boundaries")
        try:
            numeric = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("invalid sequential_boundaries") from exc
        if not math.isfinite(numeric) or not numeric.is_integer() or numeric < 0.0:
            raise ValueError("invalid sequential_boundaries")
        values.append(int(numeric))
    if values != sorted(set(values)):
        raise ValueError("sequential_boundaries must be sorted and unique")
    return tuple(values)


def _symmetric_boundary_distance(left: Sequence[int], right: Sequence[int]) -> float:
    """Return the bidirectional nearest-cut distance in residue-index units."""
    if not left and not right:
        return 0.0
    if not left or not right:
        raise ValueError("inconsistent empty sibling boundary evidence")

    def directed(source: Sequence[int], target: Sequence[int]) -> float:
        return sum(min(abs(value - other) for other in target) for value in source) / len(
            source
        )

    result = 0.5 * (directed(left, right) + directed(right, left))
    if not math.isfinite(result):
        raise ValueError("non-finite sibling boundary distance")
    return result


def add_sibling_relative_features(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Return copied rows with the frozen 47 lattice-relative fields populated."""
    identities = _unique_identities(rows)
    validated: list[dict[str, object]] = []
    for row, identity in zip(rows, identities):
        source_values = {
            name: _finite(row, name) for name in RELATIVE_CORE_FEATURES
        }
        hierarchy = {
            name: float(_integral(row, name, positive=False))
            for name in _HIERARCHY_INTEGER_FIELDS
        }
        hierarchy.update({name: _finite(row, name) for name in _HIERARCHY_FLOAT_FIELDS})
        validated.append(
            {
                "identity": identity,
                "row": row,
                "source_values": source_values,
                "boundaries": _boundary_set(row),
                "hierarchy": hierarchy,
                "max_cr": _finite(row, "max_cr"),
                "density_min": _finite(row, "density_min"),
            }
        )

    groups: dict[int, list[dict[str, object]]] = defaultdict(list)
    for item in validated:
        groups[item["identity"][1]].append(item)  # type: ignore[index]

    populated: list[dict[str, object]] = []
    for item in validated:
        identity = item["identity"]
        group = groups[identity[1]]  # type: ignore[index]
        output = dict(item["row"])  # type: ignore[arg-type]
        for name in RELATIVE_CORE_FEATURES:
            value = item["source_values"][name]  # type: ignore[index]
            group_values = [other["source_values"][name] for other in group]  # type: ignore[index]
            less = sum(other < value for other in group_values)
            equal = sum(other == value for other in group_values)
            output[f"sibling_percentile_{name}"] = (
                less + 0.5 * equal
            ) / len(group_values)
        for name in RELATIVE_CORE_FEATURES:
            value = item["source_values"][name]  # type: ignore[index]
            ordered = sorted(other["source_values"][name] for other in group)  # type: ignore[index]
            middle = len(ordered) // 2
            median = (
                ordered[middle]
                if len(ordered) % 2
                else (ordered[middle - 1] + ordered[middle]) / 2.0
            )
            output[f"sibling_median_delta_{name}"] = value - median

        output.update(item["hierarchy"])  # type: ignore[arg-type]
        siblings = [other for other in group if other is not item]
        if siblings:
            nearest = min(
                siblings,
                key=lambda sibling: (
                    _symmetric_boundary_distance(
                        item["boundaries"], sibling["boundaries"]  # type: ignore[arg-type]
                    ),
                    sibling["identity"][2],  # type: ignore[index]
                ),
            )
            output["sibling_nearest_cr_delta"] = (
                item["max_cr"] - nearest["max_cr"]  # type: ignore[operator]
            )
            output["sibling_nearest_density_delta"] = (
                item["density_min"] - nearest["density_min"]  # type: ignore[operator]
            )
        else:
            output["sibling_nearest_cr_delta"] = 0.0
            output["sibling_nearest_density_delta"] = 0.0

        if not all(
            math.isfinite(float(output[name])) for name in RELATIVE_HIERARCHY_FEATURES
        ):
            raise ValueError("non-finite relative feature")
        populated.append(output)

    return sorted(
        populated,
        key=lambda row: (
            _integral(row, "num_domains", positive=True),
            row["canonical_delineation"],
        ),
    )


def build_count_rows(
    chain_features: Mapping[str, object],
    candidate_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Build the frozen 98-value per-count rows from one complete population."""
    if not candidate_rows:
        raise ValueError("candidate population is empty")
    identities = _unique_identities(candidate_rows)
    counts = [identity[1] for identity in identities]
    total = len(candidate_rows)
    unique_counts = sorted(set(counts))

    if _integral(chain_features, "chain_candidate_total", positive=False) != total:
        raise ValueError("candidate population mismatch")
    if (
        _integral(chain_features, "chain_available_count_total", positive=False)
        != len(unique_counts)
    ):
        raise ValueError("count population mismatch")

    frequencies = {count: counts.count(count) for count in unique_counts}
    modal = min(unique_counts, key=lambda count: (-frequencies[count], count))
    if _finite(chain_features, "chain_modal_count") != float(modal):
        raise ValueError("modal count mismatch")
    expected_histogram_counts = [0] * 21
    for count, frequency in frequencies.items():
        expected_histogram_counts[count - 1 if count <= 20 else 20] += frequency
    for index, frequency in enumerate(expected_histogram_counts):
        name = (
            f"chain_count_hist_{index + 1}"
            if index < 20
            else "chain_count_hist_21_plus"
        )
        if _finite(chain_features, name) != frequency / total:
            raise ValueError("count histogram mismatch")

    grouped: dict[int, list[Mapping[str, object]]] = defaultdict(list)
    for row, count in zip(candidate_rows, counts):
        grouped[count].append(row)
    summaries: dict[int, dict[str, float]] = {}
    for count in unique_counts:
        values: dict[str, float] = {}
        for source in COUNT_SUMMARY_SOURCES:
            source_values = [_finite(row, source) for row in grouped[count]]
            for stat, value in zip(
                ("min", "mean", "max"),
                (
                    min(source_values),
                    sum(source_values) / len(source_values),
                    max(source_values),
                ),
            ):
                values[f"count_{source}_{stat}"] = value
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError("non-finite count summary")
        summaries[count] = values

    result: list[dict[str, object]] = []
    for index, count in enumerate(unique_counts):
        lower = unique_counts[index - 1] if index > 0 else None
        higher = unique_counts[index + 1] if index + 1 < len(unique_counts) else None
        values: dict[str, float] = {
            "count_num_domains": float(count),
            "count_n_candidates": float(frequencies[count]),
            "count_candidate_fraction": frequencies[count] / total,
            "count_modal_distance": abs(float(count) - modal),
            "count_has_lower": float(lower is not None),
            "count_has_higher": float(higher is not None),
            "count_lower_gap": float(count - lower) if lower is not None else 0.0,
            "count_higher_gap": float(higher - count) if higher is not None else 0.0,
            **summaries[count],
        }
        for name in COUNT_SUMMARIES:
            values[f"{name}_delta_lower"] = (
                summaries[count][name] - summaries[lower][name]
                if lower is not None
                else 0.0
            )
        for name in COUNT_SUMMARIES:
            values[f"{name}_delta_higher"] = (
                summaries[count][name] - summaries[higher][name]
                if higher is not None
                else 0.0
            )
        if set(values) != set(COUNT_ITEM_FEATURES):
            raise ValueError("count feature schema mismatch")
        ordered = {name: float(values[name]) for name in COUNT_ITEM_FEATURES}
        if not all(math.isfinite(value) for value in ordered.values()):
            raise ValueError("non-finite count feature")
        result.append(ordered)
    return result
