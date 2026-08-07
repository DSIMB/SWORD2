import numpy as np
import pytest

from benchmark.candidate_geometry import (
    CANDIDATE_GEOMETRY_FIELDS,
    candidate_geometry_features,
    parse_delineation,
)


def test_parse_delineation_preserves_discontinuous_domains():
    assert parse_delineation("0-4;10-14 5-9") == [[(0, 4), (10, 14)], [(5, 9)]]


def test_candidate_geometry_features_are_finite_and_capture_topology():
    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [3.8, 0.0, 0.0],
            [0.0, 3.8, 0.0],
            [0.0, 0.0, 3.8],
            [20.0, 0.0, 0.0],
            [23.8, 0.0, 0.0],
            [20.0, 3.8, 0.0],
            [20.0, 0.0, 3.8],
        ]
    )

    features = candidate_geometry_features(coordinates, "0-3 4-7")

    assert set(features) == set(CANDIDATE_GEOMETRY_FIELDS)
    assert all(np.isfinite(value) for value in features.values())
    assert features["n_segments"] == 2.0
    assert features["n_discontinuous"] == 0.0
    assert features["size_balance"] == 1.0
    assert features["largest_domain_fraction"] == 0.5


def test_candidate_geometry_features_reject_out_of_range_delineation():
    coordinates = np.zeros((4, 3))
    with pytest.raises(ValueError, match="invalid candidate partition"):
        candidate_geometry_features(coordinates, "0-3 4-5")


@pytest.mark.parametrize("delineation", ["", "3-1", "not-a-segment"])
def test_candidate_geometry_features_reject_empty_or_malformed_delineation(delineation):
    coordinates = np.zeros((4, 3))
    with pytest.raises(ValueError, match="invalid candidate partition"):
        candidate_geometry_features(coordinates, delineation)


def test_candidate_geometry_features_reject_empty_coordinates():
    with pytest.raises(ValueError, match="invalid candidate partition"):
        candidate_geometry_features(np.empty((0, 3)), "0")


@pytest.mark.parametrize("kind", ["negative", "above_one", "asymmetric"])
def test_candidate_geometry_features_reject_invalid_contact_probabilities(kind):
    coordinates = np.zeros((4, 3))
    contacts = np.eye(4)
    if kind == "negative":
        contacts[0, 1] = contacts[1, 0] = -0.1
    elif kind == "above_one":
        contacts[0, 1] = contacts[1, 0] = 1.1
    else:
        contacts[0, 1] = 0.5

    with pytest.raises(ValueError, match="invalid candidate partition"):
        candidate_geometry_features(coordinates, "0-3", contacts)
