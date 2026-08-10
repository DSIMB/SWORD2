import json
import math
from pathlib import Path

import pytest

from benchmark.factorized_ranker.lattice_features import (
    _symmetric_boundary_distance,
    add_sibling_relative_features,
    build_count_rows,
)
from benchmark.factorized_ranker.schema import (
    BASE_CANDIDATE_FEATURES,
    BOUNDARY_LOCAL_FEATURES,
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    DISCONTINUITY_FEATURES,
    DOMAIN_CONDITIONED_FEATURES,
    GLOBAL_FEATURES,
    RELATIVE_CORE_FEATURES,
    RELATIVE_HIERARCHY_FEATURES,
)


FIXTURE_PATH = (
    Path(__file__).parents[1] / "fixtures" / "factorized_features_synthetic.json"
)


def canonical_json(value):
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def candidate_row(
    source_index,
    num_domains,
    canonical,
    *,
    min_size,
    max_cr,
    density_min,
    boundaries,
    legacy_distance=None,
):
    seed = source_index + 1.0
    row = {
        name: seed + index / 10.0
        for index, name in enumerate(BASE_CANDIDATE_FEATURES)
    }
    row.update(
        {
            "source_index": source_index,
            "num_domains": num_domains,
            "canonical_delineation": canonical,
            "sequential_boundaries": list(boundaries),
            "legacy_distance": (
                seed / 10.0 if legacy_distance is None else legacy_distance
            ),
            "hierarchy_first_appearance_level": source_index % 3,
            "hierarchy_persistence_levels": source_index % 4 + 1,
            "hierarchy_parent_merge_margin": seed / 20.0,
            "hierarchy_child_merge_margin": seed / 25.0,
            "hierarchy_path_count": source_index + 2,
            "num_domains": float(num_domains),
            "min_size": float(min_size),
            "max_cr": float(max_cr),
            "density_min": float(density_min),
            "modal_count_distance": float(abs(num_domains - 2)),
        }
    )
    for offset, name in enumerate(
        (*DOMAIN_CONDITIONED_FEATURES, *BOUNDARY_LOCAL_FEATURES)
    ):
        row[name] = seed * 10.0 + offset / 100.0
    for offset, name in enumerate(DISCONTINUITY_FEATURES):
        row[name] = seed * 100.0 + offset / 100.0
    return row


def sibling_rows(values):
    return [
        candidate_row(
            index,
            2,
            f"candidate-{index}",
            min_size=value,
            max_cr=10.0 + index,
            density_min=20.0 + index,
            boundaries=[index + 2],
        )
        for index, value in enumerate(values)
    ]


def global_fixture(rows):
    counts = [int(row["num_domains"]) for row in rows]
    frequencies = {count: counts.count(count) for count in set(counts)}
    modal = min(frequencies, key=lambda count: (-frequencies[count], count))
    values = {name: 0.0 for name in GLOBAL_FEATURES}
    values["chain_candidate_total"] = float(len(rows))
    values["chain_available_count_total"] = float(len(frequencies))
    values["chain_modal_count"] = float(modal)
    for count, frequency in frequencies.items():
        suffix = str(count) if count <= 20 else "21_plus"
        values[f"chain_count_hist_{suffix}"] += frequency / len(rows)
    return values


def count_row(rows, num_domains):
    return next(row for row in rows if row["count_num_domains"] == num_domains)


def test_sibling_percentile_uses_midrank_and_is_permutation_invariant():
    rows = sibling_rows(values=[1.0, 1.0, 3.0])
    first = add_sibling_relative_features(rows)
    second = add_sibling_relative_features(list(reversed(rows)))
    assert canonical_json(first) == canonical_json(second)
    assert [r["sibling_percentile_min_size"] for r in first] == [
        1 / 3,
        1 / 3,
        5 / 6,
    ]


def test_sibling_median_delta_handles_odd_and_even_groups_exactly():
    odd = add_sibling_relative_features(sibling_rows([1.0, 2.0, 100.0]))
    assert [row["sibling_median_delta_min_size"] for row in odd] == [
        -1.0,
        0.0,
        98.0,
    ]

    even = add_sibling_relative_features(sibling_rows([1.0, 3.0]))
    assert [row["sibling_median_delta_min_size"] for row in even] == [
        -1.0,
        1.0,
    ]


def test_asymmetric_boundary_set_distance_is_bidirectional_mean():
    assert _symmetric_boundary_distance([1], [0, 2, 10]) == 0.5 * (
        1.0 + 11.0 / 3.0
    )


def test_nearest_sibling_ties_choose_lexicographically_smallest_canonical():
    rows = [
        candidate_row(
            0,
            2,
            "middle",
            min_size=1.0,
            max_cr=10.0,
            density_min=3.0,
            boundaries=[3],
        ),
        candidate_row(
            1,
            2,
            "alpha",
            min_size=2.0,
            max_cr=4.0,
            density_min=1.0,
            boundaries=[2],
        ),
        candidate_row(
            2,
            2,
            "zulu",
            min_size=3.0,
            max_cr=7.0,
            density_min=2.0,
            boundaries=[4],
        ),
    ]

    middle = next(
        row
        for row in add_sibling_relative_features(rows)
        if row["canonical_delineation"] == "middle"
    )
    assert middle["sibling_nearest_cr_delta"] == 6.0
    assert middle["sibling_nearest_density_delta"] == 2.0


def test_one_candidate_group_has_zero_nearest_sibling_deltas():
    row = candidate_row(
        0,
        7,
        "only",
        min_size=2.0,
        max_cr=8.0,
        density_min=5.0,
        boundaries=[],
    )
    [result] = add_sibling_relative_features([row])
    assert result["sibling_percentile_min_size"] == 0.5
    assert result["sibling_nearest_cr_delta"] == 0.0
    assert result["sibling_nearest_density_delta"] == 0.0


def test_sibling_population_adds_only_frozen_relative_fields():
    rows = sibling_rows([1.0, 2.0, 3.0])
    before = [dict(row) for row in rows]
    result = add_sibling_relative_features(rows)

    assert rows == before
    assert tuple(name for name in result[0] if name not in before[0]) == (
        *RELATIVE_HIERARCHY_FEATURES[:40],
        "sibling_nearest_cr_delta",
        "sibling_nearest_density_delta",
    )
    for original, populated in zip(before, result):
        assert all(populated[name] == value for name, value in original.items())
        assert all(math.isfinite(populated[name]) for name in RELATIVE_HIERARCHY_FEATURES)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda row: row.pop("hierarchy_path_count"),
        lambda row: row.__setitem__(RELATIVE_CORE_FEATURES[0], math.nan),
        lambda row: row.__setitem__("sequential_boundaries", [3, 2]),
        lambda row: row.__setitem__("num_domains", 0),
    ],
)
def test_sibling_population_rejects_missing_or_invalid_evidence(mutation):
    rows = sibling_rows([1.0, 2.0])
    mutation(rows[1])
    with pytest.raises(ValueError):
        add_sibling_relative_features(rows)


def test_nonconsecutive_neighbor_deltas_include_gaps():
    candidates = [
        candidate_row(
            count,
            count,
            f"count-{count}",
            min_size=count,
            max_cr=count * 2.0,
            density_min=count * 3.0,
            boundaries=[],
        )
        for count in (2, 4, 7)
    ]
    counts = build_count_rows(global_fixture(candidates), candidates)
    four = next(row for row in counts if row["count_num_domains"] == 4)
    assert four["count_lower_gap"] == 2.0
    assert four["count_higher_gap"] == 3.0
    assert four["count_max_cr_mean_delta_lower"] == (
        four["count_max_cr_mean"] - count_row(counts, 2)["count_max_cr_mean"]
    )


def test_count_rows_use_exact_population_modal_summary_and_missing_neighbor_values():
    rows = [
        candidate_row(
            0,
            2,
            "a",
            min_size=1.0,
            max_cr=2.0,
            density_min=8.0,
            boundaries=[2],
            legacy_distance=0.5,
        ),
        candidate_row(
            1,
            2,
            "b",
            min_size=3.0,
            max_cr=6.0,
            density_min=4.0,
            boundaries=[3],
            legacy_distance=1.5,
        ),
        candidate_row(
            2,
            4,
            "c",
            min_size=2.0,
            max_cr=10.0,
            density_min=2.0,
            boundaries=[1, 3, 5],
            legacy_distance=2.5,
        ),
    ]
    counts = build_count_rows(global_fixture(rows), rows)
    two = count_row(counts, 2)
    four = count_row(counts, 4)

    assert tuple(two) == COUNT_ITEM_FEATURES
    assert len(two) == 98
    assert two["count_n_candidates"] == 2.0
    assert two["count_candidate_fraction"] == 2.0 / 3.0
    assert two["count_modal_distance"] == 0.0
    assert two["count_legacy_distance_min"] == 0.5
    assert two["count_legacy_distance_mean"] == 1.0
    assert two["count_legacy_distance_max"] == 1.5
    assert two["count_max_cr_mean"] == 4.0
    assert two["count_has_lower"] == 0.0
    assert two["count_lower_gap"] == 0.0
    assert all(
        value == 0.0
        for name, value in two.items()
        if name.endswith("_delta_lower")
    )
    assert four["count_has_higher"] == 0.0
    assert four["count_higher_gap"] == 0.0
    assert all(math.isfinite(value) for row in counts for value in row.values())


@pytest.mark.parametrize(
    "mutate_global",
    [
        lambda global_: global_.__setitem__("chain_candidate_total", 99.0),
        lambda global_: global_.__setitem__("chain_available_count_total", 99.0),
        lambda global_: global_.__setitem__("chain_modal_count", 4.0),
        lambda global_: global_.__setitem__("chain_count_hist_2", 0.0),
    ],
)
def test_count_rows_fail_closed_for_mixed_global_population(mutate_global):
    rows = sibling_rows([1.0, 2.0, 3.0])
    global_ = global_fixture(rows)
    mutate_global(global_)
    with pytest.raises(ValueError):
        build_count_rows(global_, rows)


def test_count_rows_validate_every_histogram_bin():
    rows = sibling_rows([1.0, 2.0, 3.0])
    for count in range(1, 21):
        global_ = global_fixture(rows)
        global_[f"chain_count_hist_{count}"] += 0.125
        with pytest.raises(ValueError):
            build_count_rows(global_, rows)
    global_ = global_fixture(rows)
    global_["chain_count_hist_21_plus"] += 0.125
    with pytest.raises(ValueError):
        build_count_rows(global_, rows)


@pytest.mark.parametrize("value", [-1.0, 0.5, math.nan, math.inf])
def test_count_rows_require_integral_finite_global_population_fields(value):
    rows = sibling_rows([1.0, 2.0, 3.0])
    for name in ("chain_candidate_total", "chain_available_count_total"):
        global_ = global_fixture(rows)
        global_[name] = value
        with pytest.raises(ValueError):
            build_count_rows(global_, rows)
    global_ = global_fixture(rows)
    global_["chain_modal_count"] = value
    with pytest.raises(ValueError):
        build_count_rows(global_, rows)


def test_count_rows_reject_duplicate_identity_and_nonfinite_source():
    rows = sibling_rows([1.0, 2.0])
    duplicate = [rows[0], dict(rows[0])]
    with pytest.raises(ValueError):
        build_count_rows(global_fixture(duplicate), duplicate)

    rows[1]["legacy_distance"] = math.inf
    with pytest.raises(ValueError):
        build_count_rows(global_fixture(rows), rows)


def test_count_summary_order_is_source_major_and_boundary_coil_fraction_is_last():
    row = candidate_row(
        0,
        2,
        "a",
        min_size=11.0,
        max_cr=12.0,
        density_min=13.0,
        boundaries=[2],
        legacy_distance=10.0,
    )
    for name, value in {
        "mean_density": 14.0,
        "contact_q_mean": 15.0,
        "contact_q_max": 16.0,
        "n_segments": 17.0,
        "n_discontinuous": 18.0,
        "boundary_coil_fraction": 19.0,
    }.items():
        row[name] = value
    [count] = build_count_rows(global_fixture([row]), [row])
    assert list(count.values())[8:38] == [
        value
        for value in range(10, 20)
        for _ in ("min", "mean", "max")
    ]
    assert COUNT_ITEM_FEATURES[35:38] == (
        "count_boundary_coil_fraction_min",
        "count_boundary_coil_fraction_mean",
        "count_boundary_coil_fraction_max",
    )


def synthetic_rows():
    return [
        candidate_row(0, 2, "0-2 3-11", min_size=1, max_cr=2, density_min=10, boundaries=[2]),
        candidate_row(1, 2, "0-3 4-11", min_size=1, max_cr=4, density_min=9, boundaries=[3]),
        candidate_row(2, 2, "0-4 5-11", min_size=3, max_cr=8, density_min=7, boundaries=[4]),
        candidate_row(3, 4, "0-1 2-4 5-8 9-11", min_size=2, max_cr=5, density_min=6, boundaries=[1, 4, 8]),
        candidate_row(4, 4, "0-2 3-5 6-8 9-11", min_size=3, max_cr=9, density_min=3, boundaries=[2, 5, 8]),
        candidate_row(5, 7, "0 1 2 3-4 5-6 7-8 9-11", min_size=1, max_cr=12, density_min=1, boundaries=[0, 1, 2, 4, 6, 8]),
    ]


def synthetic_fixture(rows=None):
    inputs = sorted(
        synthetic_rows() if rows is None else rows,
        key=lambda row: (row["num_domains"], row["canonical_delineation"]),
    )
    populated = add_sibling_relative_features(inputs)
    global_ = global_fixture(inputs)
    count_rows = build_count_rows(global_, populated)
    relative = set(RELATIVE_HIERARCHY_FEATURES)
    return {
        "candidate_feature_names": list(CANDIDATE_FEATURES),
        "count_feature_names": list(COUNT_ITEM_FEATURES),
        "global_feature_names": list(GLOBAL_FEATURES),
        "global_values": [global_[name] for name in GLOBAL_FEATURES],
        "candidate_inputs": [
            {
                "source_index": row["source_index"],
                "num_domains": int(row["num_domains"]),
                "delineation": row["canonical_delineation"],
                "canonical_delineation": row["canonical_delineation"],
                "legacy_distance": row["legacy_distance"],
                "sequential_boundaries": row["sequential_boundaries"],
                "hierarchy": {
                    "first_appearance_level": row["hierarchy_first_appearance_level"],
                    "persistence_levels": row["hierarchy_persistence_levels"],
                    "parent_merge_margin": row["hierarchy_parent_merge_margin"],
                    "child_merge_margin": row["hierarchy_child_merge_margin"],
                    "hierarchy_path_count": row["hierarchy_path_count"],
                },
                "values_before": [
                    0.0 if name in relative else row[name]
                    for name in CANDIDATE_FEATURES
                ],
            }
            for row in inputs
        ],
        "candidate_rows": [
            {
                "source_index": row["source_index"],
                "num_domains": int(row["num_domains"]),
                "canonical_delineation": row["canonical_delineation"],
                "legacy_distance": row["legacy_distance"],
                "values": [row[name] for name in CANDIDATE_FEATURES],
            }
            for row in populated
        ],
        "count_rows": [
            {
                "num_domains": int(row["count_num_domains"]),
                "values": [row[name] for name in COUNT_ITEM_FEATURES],
            }
            for row in count_rows
        ],
    }


def test_synthetic_fixture_has_exact_committed_canonical_json_bytes():
    payload = (canonical_json(synthetic_fixture()) + "\n").encode("utf-8")
    assert FIXTURE_PATH.read_bytes() == payload


def test_synthetic_fixture_is_byte_identical_for_reversed_inputs():
    forward = (canonical_json(synthetic_fixture(synthetic_rows())) + "\n").encode()
    reverse = (
        canonical_json(synthetic_fixture(list(reversed(synthetic_rows())))) + "\n"
    ).encode()
    assert reverse == forward
