"""
Unit tests for benchmark/train_reranker.py helper functions.

Run with:
    python3 -m pytest benchmark/tests/test_train_reranker.py -v
"""

import pytest
import pandas as pd

from benchmark.train_reranker import (
    parse_domain_sizes,
    parse_junctions,
    compute_junction_support_for_group,
    mean_junction_support,
    compute_features,
    FEATURES,
)


# ---------------------------------------------------------------------------
# parse_domain_sizes
# ---------------------------------------------------------------------------

class TestParseDomainSizes:
    def test_simple_three_domains(self):
        # "0-29" = 30 res, "30-70" = 41, "71-112" = 42
        assert parse_domain_sizes("0-29 30-70 71-112") == [30, 41, 42]

    def test_discontinuous_domain(self):
        # "71-112;147-182" is one domain: 42 + 36 = 78 residues
        assert parse_domain_sizes("0-29 71-112;147-182 183-213") == [30, 78, 31]

    def test_full_training_example(self):
        # Taken from the first rows of training_table.csv
        sizes = parse_domain_sizes("0-29 30-70 71-112;147-182 113-146 183-213")
        assert sizes == [30, 41, 78, 34, 31]

    def test_single_domain(self):
        assert parse_domain_sizes("0-99") == [100]

    def test_empty_string(self):
        assert parse_domain_sizes("") == []


# ---------------------------------------------------------------------------
# parse_junctions
# ---------------------------------------------------------------------------

class TestParseJunctions:
    def test_simple_two_domains(self):
        # "0-29" -> (0, 29); "30-70" -> (30, 70)
        assert parse_junctions("0-29 30-70") == [0, 29, 30, 70]

    def test_discontinuous_token(self):
        # "71-112;147-182": dash=2 → start=71; rfind('-') at index 10 → end=182
        assert parse_junctions("71-112;147-182") == [71, 182]

    def test_mixed_continuous_and_discontinuous(self):
        junctions = parse_junctions("0-29 71-112;147-182 183-213")
        assert junctions == [0, 29, 71, 182, 183, 213]

    def test_single_domain(self):
        assert parse_junctions("0-99") == [0, 99]

    def test_empty_string(self):
        assert parse_junctions("") == []


# ---------------------------------------------------------------------------
# compute_junction_support_for_group
# ---------------------------------------------------------------------------

class TestComputeJunctionSupportForGroup:
    def test_three_candidates(self):
        delineations = [
            "0-50 51-100",   # junctions: 0, 50, 51, 100
            "0-50 51-100",   # junctions: 0, 50, 51, 100
            "0-60 61-100",   # junctions: 0, 60, 61, 100
        ]
        support = compute_junction_support_for_group(delineations)
        # 0 and 100 appear in all 3 candidates
        assert support[0] == pytest.approx(1.0)
        assert support[100] == pytest.approx(1.0)
        # 50 and 51 appear in 2 of 3
        assert support[50] == pytest.approx(2 / 3)
        assert support[51] == pytest.approx(2 / 3)
        # 60 and 61 appear in 1 of 3
        assert support[60] == pytest.approx(1 / 3)
        assert support[61] == pytest.approx(1 / 3)
        # Positions not present are absent from dict (not 0.0)
        assert 49 not in support

    def test_single_candidate(self):
        support = compute_junction_support_for_group(["0-99"])
        assert support[0] == pytest.approx(1.0)
        assert support[99] == pytest.approx(1.0)

    def test_support_values_bounded(self):
        """All support values must be in (0, 1]."""
        delineations = [
            "0-50 51-100",
            "0-60 61-100",
            "0-40 41-100",
        ]
        support = compute_junction_support_for_group(delineations)
        for val in support.values():
            assert 0.0 < val <= 1.0


# ---------------------------------------------------------------------------
# mean_junction_support
# ---------------------------------------------------------------------------

class TestMeanJunctionSupport:
    def _support(self):
        return {0: 1.0, 50: 2 / 3, 51: 2 / 3, 100: 1.0}

    def test_all_positions_present(self):
        support = self._support()
        # junctions("0-50 51-100") = [0, 50, 51, 100]
        # mean = (1.0 + 2/3 + 2/3 + 1.0) / 4
        expected = (1.0 + 2 / 3 + 2 / 3 + 1.0) / 4
        assert mean_junction_support("0-50 51-100", support) == pytest.approx(expected)

    def test_missing_positions_default_to_zero(self):
        # Support only knows about 0 and 100; 60 and 61 are missing
        support = {0: 1.0, 100: 1.0}
        # junctions("0-60 61-100") = [0, 60, 61, 100]
        # mean = (1.0 + 0 + 0 + 1.0) / 4 = 0.5
        assert mean_junction_support("0-60 61-100", support) == pytest.approx(0.5)

    def test_empty_delineation_returns_zero(self):
        assert mean_junction_support("", {0: 1.0}) == pytest.approx(0.0)

    def test_empty_support_returns_zero(self):
        # Positions are parsed but nothing is in the support dict
        assert mean_junction_support("0-99", {}) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# compute_features (integration test on tiny synthetic DataFrame)
# ---------------------------------------------------------------------------

class TestComputeFeatures:
    """Feature computation on a three-candidate synthetic chain."""

    @pytest.fixture
    def tiny_df(self):
        rows = [
            # candidate 0: 2 domains, continuous — oracle
            dict(
                chain_id='testA', num_domains=2, min_size=50, max_cr=0.8,
                density_min=2.0, mean_density=2.5,
                delineation='0-49 50-99',
                n_true_domains=2, n_pred_domains=2,
                ndo=0.9, iou=0.85, boundary_f1_10=0.8,
                matched_dice=0.85, d_count_acc=1.0, S=0.88, is_oracle_s=1,
            ),
            # candidate 1: 3 domains, continuous
            dict(
                chain_id='testA', num_domains=3, min_size=25, max_cr=0.7,
                density_min=1.5, mean_density=2.0,
                delineation='0-49 50-74 75-99',
                n_true_domains=2, n_pred_domains=3,
                ndo=0.7, iou=0.65, boundary_f1_10=0.6,
                matched_dice=0.65, d_count_acc=0.0, S=0.52, is_oracle_s=0,
            ),
            # candidate 2: 1 domain
            dict(
                chain_id='testA', num_domains=1, min_size=100, max_cr=0.9,
                density_min=3.0, mean_density=3.0,
                delineation='0-99',
                n_true_domains=2, n_pred_domains=1,
                ndo=0.4, iou=0.4, boundary_f1_10=0.0,
                matched_dice=0.4, d_count_acc=0.0, S=0.3, is_oracle_s=0,
            ),
        ]
        return pd.DataFrame(rows)

    def test_delineation_column_dropped(self, tiny_df):
        result = compute_features(tiny_df)
        assert 'delineation' not in result.columns

    def test_new_feature_columns_exist(self, tiny_df):
        result = compute_features(tiny_df)
        for col in ('n_discontinuous', 'size_balance', 'largest_domain_frac',
                    'mean_junction_support'):
            assert col in result.columns, f"Missing column: {col}"

    def test_n_discontinuous(self, tiny_df):
        result = compute_features(tiny_df)
        # None of the three delineations have ';'
        assert list(result['n_discontinuous']) == [0, 0, 0]

    def test_size_balance_row0(self, tiny_df):
        result = compute_features(tiny_df)
        # "0-49 50-99": both domains are 50 res → balance = 50/50 = 1.0
        assert result.iloc[0]['size_balance'] == pytest.approx(1.0)

    def test_size_balance_row1(self, tiny_df):
        result = compute_features(tiny_df)
        # "0-49 50-74 75-99": sizes [50, 25, 25], mean≈33.33, min=25
        # balance = 25 / (100/3) = 75/100 = 0.75
        assert result.iloc[1]['size_balance'] == pytest.approx(0.75)

    def test_size_balance_single_domain(self, tiny_df):
        result = compute_features(tiny_df)
        # "0-99": single domain → min == mean → balance = 1.0
        assert result.iloc[2]['size_balance'] == pytest.approx(1.0)

    def test_largest_domain_frac(self, tiny_df):
        result = compute_features(tiny_df)
        assert result.iloc[0]['largest_domain_frac'] == pytest.approx(0.5)   # 50/100
        assert result.iloc[1]['largest_domain_frac'] == pytest.approx(0.5)   # 50/100
        assert result.iloc[2]['largest_domain_frac'] == pytest.approx(1.0)   # 100/100

    def test_mean_junction_support_row0(self, tiny_df):
        result = compute_features(tiny_df)
        # All candidates: junctions for each:
        #   row0 "0-49 50-99"       -> [0, 49, 50, 99]
        #   row1 "0-49 50-74 75-99" -> [0, 49, 50, 74, 75, 99]
        #   row2 "0-99"             -> [0, 99]
        # n=3; support: {0:1, 49:2/3, 50:2/3, 99:1, 74:1/3, 75:1/3}
        # row0 mjs = (1 + 2/3 + 2/3 + 1) / 4 = (10/3) / 4 = 5/6
        assert result.iloc[0]['mean_junction_support'] == pytest.approx(5 / 6)

    def test_mean_junction_support_row2(self, tiny_df):
        result = compute_features(tiny_df)
        # row2 junctions = [0, 99]; both have support 1.0 → mjs = 1.0
        assert result.iloc[2]['mean_junction_support'] == pytest.approx(1.0)

    def test_features_list_completeness(self, tiny_df):
        """All nine features from FEATURES constant must be present after compute_features."""
        result = compute_features(tiny_df)
        for feat in FEATURES:
            assert feat in result.columns, f"Missing FEATURES entry: {feat}"
