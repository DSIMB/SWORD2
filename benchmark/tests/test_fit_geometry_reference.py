import numpy as np
import pytest

from benchmark.fit_geometry_reference import (
    ca_number_density,
    effective_sphere_radius,
    fit_gamma,
    fit_power_law,
    gyration_tensor,
    interdomain_contacts,
    parse_zero_based_domains,
    relative_shape_anisotropy,
    sym3x3_eigenvalues,
)


def test_octahedron_is_a_perfect_sphere():
    coords = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ]
    )
    eig = sym3x3_eigenvalues(gyration_tensor(coords, list(range(6))))
    assert eig[0] == pytest.approx(eig[1], abs=1e-9)
    assert eig[1] == pytest.approx(eig[2], abs=1e-9)
    assert relative_shape_anisotropy(eig) == pytest.approx(0.0, abs=1e-9)


def test_collinear_points_have_unit_anisotropy():
    coords = np.array([[x, 0.0, 0.0] for x in range(-2, 3)])
    eig = sym3x3_eigenvalues(gyration_tensor(coords, list(range(5))))
    assert relative_shape_anisotropy(eig) == pytest.approx(1.0, abs=1e-9)


def test_eigenvalues_match_known_symmetric_matrix():
    # [[2,1,0],[1,2,0],[0,0,3]] has eigenvalues {3, 3, 1}.
    m = np.array([[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    eig = sym3x3_eigenvalues(m)
    assert eig[0] == pytest.approx(3.0, abs=1e-6)
    assert eig[1] == pytest.approx(3.0, abs=1e-6)
    assert eig[2] == pytest.approx(1.0, abs=1e-6)


def test_density_matches_sphere_volume_formula():
    r_eff = effective_sphere_radius(1.0)
    assert r_eff == pytest.approx(np.sqrt(5.0 / 3.0))
    density = ca_number_density(100, r_eff)
    expected = 100 / ((4.0 / 3.0) * np.pi * r_eff**3)
    assert density == pytest.approx(expected)


def test_zero_radius_gives_zero_density():
    assert ca_number_density(10, 0.0) == 0.0


def test_interdomain_contacts_counts_pairs_within_cutoff():
    coords = np.array(
        [
            [0.0, 0.0, 0.0],  # 0: domain a
            [10.0, 0.0, 0.0],  # 1: domain a
            [3.0, 0.0, 0.0],  # 2: domain b, 3A from 0 -> contact
            [3.0, 0.0, 0.0],  # 3: domain b, 3A from 0, 7A from 1 -> contact both
            [50.0, 0.0, 0.0],  # 4: domain b, far from both -> no contact
        ]
    )
    count = interdomain_contacts(coords, [0, 1], [2, 3, 4])
    assert count == 4


def test_parse_zero_based_domains_simple():
    domains = parse_zero_based_domains("0-2,3-5")
    assert domains == [[0, 1, 2], [3, 4, 5]]


def test_parse_zero_based_domains_discontinuous_segments():
    domains = parse_zero_based_domains("0-2_7-8,3-6")
    assert domains == [[0, 1, 2, 7, 8], [3, 4, 5, 6]]


def test_fit_power_law_recovers_known_law():
    ns = np.array([10.0, 100.0, 1000.0])
    a_true, b_true = 2.0, 0.4
    rs = a_true * ns**b_true
    a, b = fit_power_law(ns, rs)
    assert a == pytest.approx(a_true, rel=1e-6)
    assert b == pytest.approx(b_true, rel=1e-6)


def test_fit_gamma_recovers_known_ratio():
    r_a = np.array([5.0, 8.0, 10.0])
    r_b = np.array([5.0, 6.0, 12.0])
    gamma_true = 0.07
    denom = np.pi * np.minimum(r_a, r_b) ** 2
    contacts = gamma_true * denom
    gamma = fit_gamma(contacts, r_a, r_b)
    assert gamma == pytest.approx(gamma_true, rel=1e-6)


def test_fit_gamma_zero_when_no_pairs():
    assert fit_gamma(np.array([]), np.array([]), np.array([])) == 0.0
