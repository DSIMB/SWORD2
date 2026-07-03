import numpy as np
import pytest

from benchmark.fit_count_calibration import fit_linear_count, round_accuracy


def test_fit_recovers_known_line():
    lengths = np.array([100.0, 200.0, 300.0, 400.0])
    counts = 1.0 + 0.01 * lengths  # exact line: intercept 1, slope 0.01
    intercept, slope = fit_linear_count(lengths, counts)
    assert intercept == pytest.approx(1.0, abs=1e-6)
    assert slope == pytest.approx(0.01, abs=1e-6)


def test_round_accuracy_perfect_when_counts_on_line():
    lengths = np.array([150.0, 300.0])
    counts = np.array([2.0, 4.0])  # round(1.0 + 0.01*L) = 2, 4
    acc, bias = round_accuracy(lengths, counts, 1.0, 0.01)
    assert acc == pytest.approx(1.0)
    assert bias == pytest.approx(0.0)


def test_round_accuracy_clamps_to_min_one_domain():
    # a very short chain must never predict < 1 domain
    lengths = np.array([10.0])
    counts = np.array([1.0])
    acc, bias = round_accuracy(lengths, counts, -5.0, 0.001)
    assert acc == pytest.approx(1.0)  # clamped prediction 1 == true 1
