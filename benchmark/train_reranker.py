#!/usr/bin/env python3
"""
Train a logistic re-ranker for SWORD2 domain candidate selection (Phase B.2).

Reads the training table produced by SWORD2_DUMP_CANDIDATES, computes
delineation-based features, trains a binary logistic regression with
chain-local z-score normalization (L-BFGS-B via scipy), and writes
the fitted weights to JSON.

Usage:
    python benchmark/train_reranker.py \
        [--training benchmark/data/training_table.csv] \
        [--output benchmark/data/reranker_weights.json]
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize


# Feature names in fixed order expected by the re-ranker.
FEATURES = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "n_discontinuous",
    "size_balance",
    "largest_domain_frac",
    "mean_junction_support",
]


# ---------------------------------------------------------------------------
# Delineation parsing helpers (exact spec from task brief)
# ---------------------------------------------------------------------------

def parse_domain_sizes(delineation: str) -> list:
    """Return total residue count for each domain.

    Each whitespace-separated token is one domain; ';'-separated segments
    within a token form a discontinuous domain.
    """
    sizes = []
    for dom_token in delineation.split():
        size = 0
        for seg in dom_token.split(';'):
            parts = seg.split('-')
            if len(parts) == 2:
                try:
                    size += int(parts[1]) - int(parts[0]) + 1
                except ValueError:
                    pass
        sizes.append(size)
    return sizes


def parse_junctions(delineation: str) -> list:
    """Return junction positions: start and end of each domain token.

    Mirrors the logic in Rust junctions.rs.
    """
    junctions = []
    for dom_token in delineation.split():
        # start = number before first '-'
        dash = dom_token.find('-')
        if dash > 0:
            try:
                junctions.append(int(dom_token[:dash]))
            except ValueError:
                pass
        # end = number after last '-'
        last_dash = dom_token.rfind('-')
        if last_dash >= 0:
            end_str = dom_token[last_dash + 1:]
            try:
                junctions.append(
                    int(
                        end_str.split(';')[0]
                        if ';' not in dom_token[last_dash:]
                        else end_str
                    )
                )
            except ValueError:
                pass
    return junctions


def compute_junction_support_for_group(delineations: list) -> dict:
    """Compute per-position support across all candidates in one chain.

    Returns {pos: count / n_candidates}.
    """
    n = len(delineations)
    counts: dict = {}
    for d in delineations:
        for pos in parse_junctions(d):
            counts[pos] = counts.get(pos, 0) + 1
    return {pos: cnt / n for pos, cnt in counts.items()}


def mean_junction_support(delineation: str, support: dict) -> float:
    """Mean junction support for a single candidate given its chain's support dict."""
    positions = parse_junctions(delineation)
    if not positions:
        return 0.0
    return sum(support.get(p, 0.0) for p in positions) / len(positions)


# ---------------------------------------------------------------------------
# Feature computation over full DataFrame
# ---------------------------------------------------------------------------

def compute_features(df: pd.DataFrame) -> pd.DataFrame:
    """Compute the four delineation-derived features and add them to *df*.

    Also drops the ``delineation`` column to free memory.  Works in two
    passes to avoid storing per-row junction position lists (memory-saving).
    """
    chain_ids = df['chain_id'].tolist()
    delineations = df['delineation'].tolist()
    n_rows = len(df)

    n_disc_arr = [0] * n_rows
    size_bal_arr = [0.0] * n_rows
    ldf_arr = [0.0] * n_rows

    # Per-chain accumulators for junction support
    chain_counts: dict = {}   # chain_id -> n_candidates
    chain_jcounts: dict = {}  # chain_id -> {pos: count}

    print(f"  Pass 1 ({n_rows:,} rows): per-row features + junction counts...", flush=True)
    t0 = time.time()

    for i, (chain_id, d) in enumerate(zip(chain_ids, delineations)):
        if i > 0 and i % 2_000_000 == 0:
            elapsed = time.time() - t0
            pct = 100 * i / n_rows
            print(f"    {pct:.0f}%  ({elapsed:.0f}s elapsed)", flush=True)

        # n_discontinuous: domains containing ';'
        tokens = d.split()
        n_disc_arr[i] = sum(1 for t in tokens if ';' in t)

        # size_balance and largest_domain_frac
        sizes = parse_domain_sizes(d)
        if sizes:
            total = sum(sizes)
            mean_s = total / len(sizes)
            sb = min(sizes) / mean_s if mean_s > 0 else 0.0
            size_bal_arr[i] = max(0.0, min(1.0, sb))
            ldf_arr[i] = max(sizes) / total if total > 0 else 0.0

        # Accumulate junction counts
        positions = parse_junctions(d)
        if chain_id not in chain_counts:
            chain_counts[chain_id] = 0
            chain_jcounts[chain_id] = {}
        chain_counts[chain_id] += 1
        jc = chain_jcounts[chain_id]
        for p in positions:
            jc[p] = jc.get(p, 0) + 1

    print(
        f"  Pass 1 done in {time.time()-t0:.1f}s "
        f"({len(chain_counts):,} unique chains)",
        flush=True,
    )

    # Build support dicts from accumulated counts
    chain_support: dict = {
        cid: {pos: cnt / chain_counts[cid] for pos, cnt in jc.items()}
        for cid, jc in chain_jcounts.items()
    }

    print("  Pass 2: mean junction support...", flush=True)
    t1 = time.time()

    mjs_arr = [0.0] * n_rows
    for i, (chain_id, d) in enumerate(zip(chain_ids, delineations)):
        positions = parse_junctions(d)
        if positions:
            support = chain_support.get(chain_id, {})
            mjs_arr[i] = sum(support.get(p, 0.0) for p in positions) / len(positions)

    print(f"  Pass 2 done in {time.time()-t1:.1f}s", flush=True)

    df = df.copy()
    df['n_discontinuous'] = n_disc_arr
    df['size_balance'] = size_bal_arr
    df['largest_domain_frac'] = ldf_arr
    df['mean_junction_support'] = mjs_arr
    df.drop(columns=['delineation'], inplace=True)

    return df


# ---------------------------------------------------------------------------
# Chain-local z-score normalization
# ---------------------------------------------------------------------------

def chain_local_normalize(df: pd.DataFrame) -> np.ndarray:
    """Return (N, 9) float64 array of chain-local z-score normalized features.

    For each chain, subtract per-feature mean and divide by per-feature std
    (with epsilon=1e-8).  Single-candidate chains get z=0 for all features.
    """
    print("  Chain-local z-score normalization...", flush=True)
    t0 = time.time()

    X = df[FEATURES].astype(np.float64)
    chain_mean = df.groupby('chain_id')[FEATURES].transform('mean')
    chain_std = df.groupby('chain_id')[FEATURES].transform('std', ddof=1).fillna(0.0)  # ddof=1: sample std; single-candidate chains give NaN → filled with 0 below

    X_norm = (X - chain_mean) / (chain_std + 1e-8)

    print(f"  Normalization done in {time.time()-t0:.1f}s", flush=True)
    return X_norm.values.astype(np.float64)


# ---------------------------------------------------------------------------
# Logistic regression (binary cross-entropy, L-BFGS-B)
# ---------------------------------------------------------------------------

def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -500.0, 500.0)))


def _loss_and_grad(params: np.ndarray, X: np.ndarray, y: np.ndarray):
    """Binary cross-entropy loss + gradient (returned jointly for L-BFGS-B jac=True)."""
    w = params[:-1]
    b = params[-1]
    logits = X @ w + b
    probs = _sigmoid(logits)
    probs_c = np.clip(probs, 1e-10, 1.0 - 1e-10)
    loss = -np.mean(y * np.log(probs_c) + (1.0 - y) * np.log(1.0 - probs_c))
    error = probs - y
    n = float(len(y))
    grad_w = (X.T @ error) / n
    grad_b = error.sum() / n
    return float(loss), np.append(grad_w, grad_b)


def train_logistic(X: np.ndarray, y: np.ndarray):
    """Fit logistic regression. Returns (weights array of shape (9,), bias float)."""
    print("  L-BFGS-B optimizer...", flush=True)
    t0 = time.time()

    params0 = np.zeros(X.shape[1] + 1)
    result = minimize(
        fun=_loss_and_grad,
        x0=params0,
        args=(X, y),
        method='L-BFGS-B',
        jac=True,
        options={'maxiter': 1000, 'ftol': 1e-9, 'gtol': 1e-6},
    )

    elapsed = time.time() - t0
    print(f"  Optimizer done in {elapsed:.1f}s | {result.message}", flush=True)
    if not result.success:
        print(f"  WARNING: optimizer did not converge: {result.message}", file=sys.stderr)

    weights = result.x[:-1]
    bias = float(result.x[-1])
    return weights, bias


# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------

def compute_oracle_hit_rate(df: pd.DataFrame, scores: np.ndarray):
    """Fraction of chains where the top-scored candidate is the oracle_s.

    Returns (hit_rate, n_chains).
    """
    tmp = pd.DataFrame({
        'chain_id': df['chain_id'].values,
        'is_oracle_s': df['is_oracle_s'].values,
        'score': scores,
    })
    idx_max = tmp.groupby('chain_id')['score'].idxmax()
    n_chains = len(idx_max)
    n_hits = int(tmp.loc[idx_max, 'is_oracle_s'].sum())
    return n_hits / n_chains, n_chains


def compute_baseline_hit_rate(df: pd.DataFrame) -> float:
    """Fraction of chains where the candidate with n_pred==n_true is oracle_s.

    Chains with no exact match are counted as misses.
    """
    tmp = df[['chain_id', 'n_true_domains', 'n_pred_domains', 'is_oracle_s']]
    exact = tmp[tmp['n_pred_domains'] == tmp['n_true_domains']]
    n_total_chains = tmp['chain_id'].nunique()
    # For each chain with an exact-count candidate, was any of them oracle_s?
    hits_per_chain = exact.groupby('chain_id')['is_oracle_s'].max()
    return int(hits_per_chain.sum()) / n_total_chains


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Train SWORD2 logistic re-ranker and write weights JSON."
    )
    parser.add_argument(
        '--training',
        type=Path,
        default=Path('benchmark/data/training_table.csv'),
        help="Path to training_table.csv (default: benchmark/data/training_table.csv)",
    )
    parser.add_argument(
        '--output',
        type=Path,
        default=Path('benchmark/data/reranker_weights.json'),
        help="Output path for reranker_weights.json (default: benchmark/data/reranker_weights.json)",
    )
    args = parser.parse_args()

    if not args.training.exists():
        print(f"ERROR: training table not found: {args.training}", file=sys.stderr)
        sys.exit(1)

    # Fast schema validation before loading 2.8 GB
    REQUIRED_COLUMNS = [
        "chain_id", "num_domains", "min_size", "max_cr", "density_min",
        "mean_density", "delineation", "n_true_domains", "n_pred_domains",
        "is_oracle_s", "S",
    ]
    header = pd.read_csv(args.training, nrows=0)
    missing = [c for c in REQUIRED_COLUMNS if c not in header.columns]
    if missing:
        sys.exit(f"Error: training CSV is missing required columns: {missing}")

    # ------------------------------------------------------------------
    # Load
    # ------------------------------------------------------------------
    print(f"[1/5] Loading {args.training} ...", flush=True)
    t_start = time.time()
    dtypes = {
        'is_oracle_s': 'int8',
        'num_domains': 'int16',
        'min_size': 'int16',
        'n_true_domains': 'int16',
        'n_pred_domains': 'int16',
    }
    df = pd.read_csv(args.training, dtype=dtypes)
    print(
        f"    Loaded {len(df):,} rows in {time.time()-t_start:.1f}s",
        flush=True,
    )

    # ------------------------------------------------------------------
    # Compute delineation-based features
    # ------------------------------------------------------------------
    print("[2/5] Computing delineation features...", flush=True)
    df = compute_features(df)

    # ------------------------------------------------------------------
    # Chain-local normalization
    # ------------------------------------------------------------------
    print("[3/5] Normalizing features...", flush=True)
    X = chain_local_normalize(df)
    y = df['is_oracle_s'].values.astype(np.float64)

    # ------------------------------------------------------------------
    # Train
    # ------------------------------------------------------------------
    print("[4/5] Training logistic regression...", flush=True)
    weights, bias = train_logistic(X, y)

    # ------------------------------------------------------------------
    # Evaluate
    # ------------------------------------------------------------------
    print("[5/5] Evaluating...", flush=True)
    scores = X @ weights + bias
    train_hit_rate, n_chains = compute_oracle_hit_rate(df, scores)
    baseline_hit_rate = compute_baseline_hit_rate(df)

    print()
    print(f"Training oracle hit rate: {train_hit_rate*100:.1f}%  ({n_chains:,} chains)")
    print(f"Baseline oracle hit rate: {baseline_hit_rate*100:.1f}%")
    print(f"Total elapsed:            {time.time()-t_start:.1f}s")

    if train_hit_rate < 0.50:
        print(
            "WARNING: training hit rate < 50% — something may be wrong!",
            file=sys.stderr,
        )

    # ------------------------------------------------------------------
    # Write output
    # ------------------------------------------------------------------
    output_data = {
        "features": FEATURES,
        "weights": weights.tolist(),
        "bias": bias,
        "train_oracle_hit_rate": float(train_hit_rate),
        "baseline_oracle_hit_rate": float(baseline_hit_rate),
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as fh:
        json.dump(output_data, fh, indent=2)
    print(f"\nWeights written to {args.output}", flush=True)


if __name__ == '__main__':
    main()
