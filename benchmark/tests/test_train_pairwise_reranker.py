import json

import numpy as np
import pandas as pd

from benchmark.train_pairwise_reranker import (
    FEATURES,
    build_pairs,
    pairwise_logistic_loss_and_grad,
    train,
)


def _toy_table() -> pd.DataFrame:
    # Two chains. In each, candidate "a" is strictly better (higher ndo) and
    # has a higher energy_z and coil_fraction than candidate "b" — a
    # separable synthetic case the trainer should learn perfectly.
    rows = []
    for chain in ["c1", "c2"]:
        rows.append(dict(chain_id=chain, num_domains=2, min_size=30, max_cr=0.3,
                          density_min=1.0, mean_density=2.0, boundary_coil_fraction=0.9,
                          energy_z=-3.0, modal_count_distance=0.0, ndo=0.9))
        rows.append(dict(chain_id=chain, num_domains=4, min_size=10, max_cr=0.6,
                          density_min=0.5, mean_density=1.0, boundary_coil_fraction=0.1,
                          energy_z=1.0, modal_count_distance=2.0, ndo=0.3))
    return pd.DataFrame(rows)


def test_build_pairs_only_compares_within_chain():
    df = _toy_table()
    pairs = build_pairs(df)
    assert len(pairs) == 2  # one pair per chain, not cross-chain
    for i, j, label in pairs:
        assert df.loc[i, "chain_id"] == df.loc[j, "chain_id"]
        assert label in (1, -1)


def test_pairwise_logistic_loss_decreases_with_correct_sign():
    # A large positive margin in the "correct" direction should have lower
    # loss than the same margin in the "wrong" direction.
    w = np.array([1.0])
    x_i = np.array([[2.0]])
    x_j = np.array([[0.0]])
    labels = np.array([1])
    loss_correct, _ = pairwise_logistic_loss_and_grad(w, 0.0, x_i, x_j, labels)
    loss_wrong, _ = pairwise_logistic_loss_and_grad(-w, 0.0, x_i, x_j, labels)
    assert loss_correct < loss_wrong


def test_train_converges_on_separable_toy_data(tmp_path):
    df = _toy_table()
    weights_path = tmp_path / "weights.json"
    result = train(df, features=FEATURES, epochs=500, lr=0.5, out_path=weights_path)

    with open(weights_path) as f:
        saved = json.load(f)
    assert saved["features"] == FEATURES
    assert len(saved["weights"]) == len(FEATURES)

    # The trained model should rank candidate "a" above "b" in both chains.
    assert result["train_accuracy"] == 1.0
