import math

import numpy as np
import pytest

from benchmark.factorized_ranker.integrity import validate_partition
from benchmark.factorized_ranker.schema import (
    BOUNDARY_LOCAL_FEATURES,
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    DISCONTINUITY_FEATURES,
    DOMAIN_CONDITIONED_FEATURES,
    GLOBAL_FEATURES,
    pair_feature_names,
)
from benchmark.factorized_ranker.structure_features import (
    DsspHydrogenBond,
    DsspResidue,
    PeelingSummary,
    aggregate_measurements,
    compute_boundary_local_features,
    compute_discontinuity_features,
    compute_domain_conditioned_features,
    compute_global_features,
)


def _dssp(sequence: str) -> tuple[DsspResidue, ...]:
    return tuple(
        DsspResidue(
            index=index,
            ss=code,
            bridge_partners=(),
            sheet_label="",
            hydrogen_bonds=(),
            kappa=0.0,
            alpha=0.0,
        )
        for index, code in enumerate(sequence)
    )


def _continuous_fixture():
    return validate_partition("0-3 4-7", n_residues=8, declared_domains=2)


def _octahedron() -> np.ndarray:
    return np.asarray(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ]
    )


def test_schema_is_unique_finite_and_has_no_competitor_fields():
    all_names = (*GLOBAL_FEATURES, *COUNT_ITEM_FEATURES, *CANDIDATE_FEATURES)
    assert len(all_names) == len(set(all_names))
    assert not any(
        "merizo" in name.lower() or "chainsaw" in name.lower()
        for name in all_names
    )

    pair_names = pair_feature_names(GLOBAL_FEATURES, COUNT_ITEM_FEATURES)
    assert pair_names[: len(GLOBAL_FEATURES)] == GLOBAL_FEATURES
    assert pair_names[len(GLOBAL_FEATURES)] == "diff__count_num_domains"
    assert pair_names[-1] == "abs_diff__count_boundary_coil_fraction_max_delta_higher"


def test_boundary_aggregation_is_exact():
    values = aggregate_measurements(
        [{"hbond_count": 2.0}, {"hbond_count": 4.0}]
    )
    assert values["boundary_hbond_count_min"] == 2.0
    assert values["boundary_hbond_count_mean"] == 3.0
    assert values["boundary_hbond_count_max"] == 4.0


def test_global_features_use_exact_geometry_dssp_and_count_formulas():
    coordinates = _octahedron()
    contacts = np.zeros((6, 6), dtype=float)
    contacts[0, 5] = contacts[5, 0] = 0.75
    values = compute_global_features(
        coordinates,
        contacts,
        _dssp("HHHCEE"),
        PeelingSummary(n_levels=3, finest_pu_count=5),
        [2, 2, 4, 21],
    )

    assert tuple(values) == GLOBAL_FEATURES
    assert values["chain_rg_normalized"] == pytest.approx(1.0 / 6 ** (1.0 / 3.0))
    assert values["chain_inertia_ratio_21"] == pytest.approx(1.0)
    assert values["chain_inertia_ratio_31"] == pytest.approx(1.0)
    assert values["chain_nonlocal_contact_density"] == 0.0
    assert values["chain_contact_order"] == pytest.approx(5.0 / 6.0)
    assert values["chain_helix_fraction"] == 0.5
    assert values["chain_strand_fraction"] == pytest.approx(2.0 / 6.0)
    assert values["chain_coil_fraction"] == pytest.approx(1.0 / 6.0)
    assert values["chain_helix_blocks"] == 1.0
    assert values["chain_strand_blocks"] == 1.0
    assert values["chain_count_hist_2"] == 2.0
    assert values["chain_count_hist_4"] == 1.0
    assert values["chain_count_hist_21_plus"] == 1.0
    assert values["chain_modal_count"] == 2.0
    assert all(math.isfinite(value) for value in values.values())


def test_domain_conditioned_features_use_exact_contact_mass_definitions():
    coordinates = np.column_stack(
        (np.arange(14, dtype=float), np.zeros(14), np.zeros(14))
    )
    contacts = np.zeros((14, 14), dtype=float)
    contacts[:4, :4] = 0.5
    contacts[4:, 4:] = 0.75
    contacts[:4, 4:] = contacts[4:, :4] = 0.25
    np.fill_diagonal(contacts, 0.0)
    partition = validate_partition("0-3 4-13", 14, 2)

    values = compute_domain_conditioned_features(partition, coordinates, contacts)

    assert tuple(values) == DOMAIN_CONDITIONED_FEATURES
    assert values["smallest_size_fraction"] == pytest.approx(4.0 / 14.0)
    assert values["largest_size_fraction"] == pytest.approx(10.0 / 14.0)
    assert values["smallest_internal_contact_density"] == 0.0
    assert values["largest_internal_contact_density"] == pytest.approx(0.75)
    assert values["smallest_contact_order"] == pytest.approx(5.0 / 42.0)
    assert values["largest_contact_order"] == pytest.approx(11.0 / 42.0)
    assert values["domain_internal_contact_fraction_min"] == pytest.approx(3.0 / 13.0)
    assert values["domain_internal_contact_fraction_max"] == pytest.approx(33.75 / 43.75)
    assert values["domain_conductance_max"] == pytest.approx(10.0 / 16.0)
    assert values["domain_conductance_min"] == pytest.approx(10.0 / 77.5)
    assert all(math.isfinite(value) for value in values.values())


def test_boundary_features_deduplicate_structural_evidence_and_wrap_angles():
    coordinates = np.column_stack(
        (np.arange(10, dtype=float), np.zeros(10), np.zeros(10))
    )
    contacts = np.zeros((10, 10), dtype=float)
    contacts[0, 9] = contacts[9, 0] = 0.9
    dssp = list(_dssp("HHHCCEEECC"))
    dssp[3] = DsspResidue(
        3, "coil", (), "", (DsspHydrogenBond(6, -1.0),), 0.0, 0.0
    )
    dssp[7] = DsspResidue(
        7, "strand", (), "", (DsspHydrogenBond(2, -2.0),), 0.0, 0.0
    )
    dssp[4] = DsspResidue(
        4, "coil", (5,), "A", (), 20.0, 170.0
    )
    dssp[5] = DsspResidue(
        5, "strand", (4,), "A", (), 38.0, -170.0
    )
    partition = validate_partition("0-4 5-9", 10, 2)

    values = compute_boundary_local_features(
        partition, coordinates, contacts, tuple(dssp)
    )

    assert tuple(values) == BOUNDARY_LOCAL_FEATURES
    assert values["boundary_inside_coil_mean"] == 1.0
    assert values["boundary_sse_terminus_distance_mean"] == 0.0
    assert values["boundary_hbond_count_mean"] == 2.0
    assert values["boundary_hbond_energy_kcal_mean"] == -3.0
    assert values["boundary_bridge_count_mean"] == 1.0
    assert values["boundary_sheet_link_count_mean"] == 1.0
    assert values["boundary_insulation_w8_mean"] == pytest.approx(1.0 - 0.9 / 25.0)
    assert values["boundary_long_range_contact_density_mean"] == pytest.approx(0.3)
    assert values["boundary_bend_change_mean"] == pytest.approx(0.1)
    assert values["boundary_virtual_dihedral_change_mean"] == pytest.approx(1.0 / 9.0)
    assert all(math.isfinite(value) for value in values.values())


def test_boundary_inside_sse_distance_and_dihedral_sentinel_are_exact():
    coordinates = np.zeros((10, 3), dtype=float)
    contacts = np.zeros((10, 10), dtype=float)
    dssp = list(_dssp("HHHCCEEECC"))
    dssp[1] = DsspResidue(1, "H", (), "", (), 10.0, 360.0)
    dssp[2] = DsspResidue(2, "H", (), "", (), 10.0, 20.0)
    partition = validate_partition("0-1 2-9", 10, 2)

    values = compute_boundary_local_features(
        partition, coordinates, contacts, tuple(dssp)
    )

    assert values["boundary_inside_helix_mean"] == 1.0
    assert values["boundary_inside_strand_mean"] == 0.0
    assert values["boundary_inside_coil_mean"] == 0.0
    assert values["boundary_sse_terminus_distance_mean"] == 1.0
    assert values["boundary_virtual_dihedral_change_mean"] == 0.0


def test_continuous_candidate_has_explicit_zero_discontinuity_vector():
    values = compute_discontinuity_features(
        _continuous_fixture(), np.zeros((8, 8)), _dssp("CCCCCCCC")
    )
    assert tuple(values) == DISCONTINUITY_FEATURES
    assert values["has_discontinuity"] == 0.0
    assert all(values[name] == 0.0 for name in DISCONTINUITY_FEATURES[1:])


def test_discontinuous_partition_aggregates_segment_affinity_and_sheet_links():
    partition = validate_partition("0-1;6-7 2-5", 8, 2)
    contacts = np.zeros((8, 8), dtype=float)
    contacts[np.ix_([0, 1], [6, 7])] = 0.8
    contacts[np.ix_([6, 7], [0, 1])] = 0.8
    contacts[np.ix_([0, 1, 6, 7], [2, 3, 4, 5])] = 0.2
    contacts[np.ix_([2, 3, 4, 5], [0, 1, 6, 7])] = 0.2
    dssp = list(_dssp("EEEEEEEE"))
    dssp[0] = DsspResidue(0, "E", (6,), "A", (), 0.0, 0.0)
    dssp[6] = DsspResidue(6, "E", (0,), "A", (), 0.0, 0.0)

    values = compute_discontinuity_features(partition, contacts, tuple(dssp))

    assert tuple(values) == DISCONTINUITY_FEATURES
    assert values["has_discontinuity"] == 1.0
    for stat in ("min", "mean", "max"):
        assert values[f"segment_same_domain_affinity_{stat}"] == pytest.approx(0.8)
        assert values[f"segment_affinity_margin_{stat}"] == pytest.approx(0.6)
        assert values[f"segment_long_range_internal_capture_{stat}"] == 0.0
        assert values[f"segment_interface_span_entropy_{stat}"] == 0.0
        assert values[f"segment_same_domain_sheet_links_{stat}"] == 1.0
    assert all(math.isfinite(value) for value in values.values())


def test_discontinuity_capture_excludes_contacts_within_the_segment():
    partition = validate_partition("0-9;16-17 10-15", 18, 2)
    contacts = np.zeros((18, 18), dtype=float)
    contacts[0, 8] = contacts[8, 0] = 1.0
    contacts[0, 16] = contacts[16, 0] = 1.0

    values = compute_discontinuity_features(
        partition, contacts, _dssp("C" * 18)
    )

    assert values["segment_long_range_internal_capture_min"] == 1.0
    assert values["segment_long_range_internal_capture_mean"] == 1.0
    assert values["segment_long_range_internal_capture_max"] == 1.0
