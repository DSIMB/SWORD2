"""Pairwise ranking-loss trainer for SWORD2 candidate reranking.

Trains on within-chain candidate pairs, predicting sign(ndo_i - ndo_j), so
the model can never learn a corpus-wide domain-count bias the way a
pointwise classifier can (see project memory: Phase B's pointwise logistic
learned -1.33 weight on num_domains and systematically under-segmented).

Usage:
    python -m benchmark.train_pairwise_reranker \\
        --training-table benchmark/data/training_table.csv \\
        --out benchmark/data/pairwise_reranker_weights.json
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

# Fixed feature order — must match sword2-lib/src/sword/reranker.rs::FEATURES exactly.
FEATURES = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "energy_z",
    "modal_count_distance",
]


def build_pairs(df: pd.DataFrame) -> list[tuple[int, int, int]]:
    """Build (i, j, label) triples from within-chain candidate pairs.

    label = 1 means df.loc[i] should rank above df.loc[j] (higher ndo);
    label = -1 means the reverse. Pairs with equal ndo are skipped (no
    signal). Returns integer-position index pairs into df.reset_index(drop=True).
    """
    df = df.reset_index(drop=True)
    pairs: list[tuple[int, int, int]] = []
    for _, group in df.groupby("chain_id"):
        idx = group.index.tolist()
        for a in range(len(idx)):
            for b in range(a + 1, len(idx)):
                i, j = idx[a], idx[b]
                ndo_i, ndo_j = df.loc[i, "ndo"], df.loc[j, "ndo"]
                if ndo_i == ndo_j:
                    continue
                label = 1 if ndo_i > ndo_j else -1
                pairs.append((i, j, label))
    return pairs


def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-x))


def pairwise_logistic_loss_and_grad(
    w: np.ndarray,
    b: float,
    x_i: np.ndarray,
    x_j: np.ndarray,
    labels: np.ndarray,
) -> tuple[float, np.ndarray]:
    """RankNet-style pairwise logistic loss and its gradient w.r.t. w.

    score(x) = w . x + b. margin = label * (score(x_i) - score(x_j)).
    loss = mean(log(1 + exp(-margin))), a smooth surrogate for "score_i
    should exceed score_j whenever label == 1".
    """
    diff = x_i - x_j  # (n_pairs, n_features)
    margin = labels * (diff @ w + b - b)  # bias cancels in the difference; kept for API symmetry
    # Numerically stable log(1 + exp(-margin))
    loss_terms = np.logaddexp(0.0, -margin)
    loss = float(np.mean(loss_terms))

    # d/dw of mean(log(1+exp(-margin))) where margin = label * (diff @ w)
    sig = _sigmoid(-margin)  # = 1 - sigmoid(margin)
    grad = -np.mean((sig * labels)[:, None] * diff, axis=0)
    return loss, grad


def train(
    df: pd.DataFrame,
    features: list[str],
    epochs: int = 2000,
    lr: float = 0.1,
    l2: float = 1e-3,
    out_path: Path | None = None,
) -> dict:
    """Train a linear pairwise reranker with full-batch gradient descent.

    Features are chain-local z-score normalized before pairing (matching the
    normalization the Rust inference side applies at run time — see
    sword2-lib/src/sword/reranker.rs).
    """
    df = df.reset_index(drop=True).copy()
    for col in features:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=features + ["ndo"]).reset_index(drop=True)

    # Chain-local z-score normalization (per chain_id group), matching inference.
    def _zscore(group: pd.DataFrame) -> pd.DataFrame:
        for col in features:
            std = group[col].std(ddof=0)
            mean = group[col].mean()
            group[col] = (group[col] - mean) / (std + 1e-8)
        return group

    normed = df.groupby("chain_id", group_keys=False).apply(_zscore)

    pairs = build_pairs(pd.concat([df["chain_id"], df["ndo"]], axis=1))
    if not pairs:
        raise ValueError("No trainable pairs — need at least 2 candidates with differing ndo per chain")

    x = normed[features].to_numpy(dtype=float)
    idx_i = np.array([p[0] for p in pairs])
    idx_j = np.array([p[1] for p in pairs])
    labels = np.array([p[2] for p in pairs], dtype=float)
    x_i, x_j = x[idx_i], x[idx_j]

    w = np.zeros(len(features))
    b = 0.0
    for epoch in range(epochs):
        loss, grad = pairwise_logistic_loss_and_grad(w, b, x_i, x_j, labels)
        grad = grad + l2 * w
        w -= lr * grad
        if epoch % max(1, epochs // 10) == 0:
            log.info("epoch %d: loss=%.4f", epoch, loss)

    scores_i = x_i @ w
    scores_j = x_j @ w
    pred = np.sign(scores_i - scores_j)
    train_accuracy = float(np.mean(pred == labels))

    result = {
        "features": features,
        "weights": w.tolist(),
        "bias": 0.0,  # bias cancels in pairwise scoring; kept for schema stability
        "feature_means": [],  # per-chain normalization is applied at inference time, not global
        "feature_stds": [],
        "train_accuracy": train_accuracy,
        "n_pairs": len(pairs),
    }

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            json.dump({k: result[k] for k in ("features", "weights", "bias")}, f, indent=2)
        log.info("Wrote weights to %s (train_accuracy=%.3f, n_pairs=%d)", out_path, train_accuracy, len(pairs))

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Train the pairwise candidate reranker")
    parser.add_argument("--training-table", type=Path, default=Path("benchmark/data/training_table.csv"))
    parser.add_argument("--out", type=Path, default=Path("benchmark/data/pairwise_reranker_weights.json"))
    parser.add_argument("--epochs", type=int, default=2000)
    parser.add_argument("--lr", type=float, default=0.1)
    args = parser.parse_args()

    df = pd.read_csv(args.training_table)
    train(df, features=FEATURES, epochs=args.epochs, lr=args.lr, out_path=args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
