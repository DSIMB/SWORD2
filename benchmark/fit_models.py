"""Refit prediction_model() and distance_model() constants from a training table.

Usage:
    python -m benchmark.fit_models \\
        --training benchmark/data/training_table.csv \\
        --output benchmark/data/fitted_constants.json

Inputs:
  training_table.csv  -- produced by benchmark/build_training_table.py.
  Columns used:
    chain_id, num_domains, max_cr, mean_density, is_oracle_s

Outputs:
  fitted_constants.json  -- new constants for both models, ready to paste into Rust.
  Also prints Rust snippet to stdout.

Algorithm:
  prediction_model():
    For each chain we have a sequence of (max_cr, mean_density) values at each num_domains
    level, plus a flag for which level is the S-oracle.  We simulate the prediction_model
    reverse-scan in Python and minimise the fraction of chains where the selected
    num_domains != the oracle num_domains (0/1 loss, Nelder-Mead).

  distance_model():
    Within each (chain_id, num_domains) group we rank candidates by the model's signed
    distance score and minimise the mean rank of the S-oracle candidate.
"""
from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Current (baseline) constants
# ---------------------------------------------------------------------------

BASELINE_PREDICTION = {
    "diag_intercept": 2.818831,
    "diag_slope": 3.582524,
    "diag_inter_v": 0.09434462,
    "horizontal_lim": 3.166823,
    "vertical_lim": 0.231845,
}

BASELINE_DISTANCE = {
    "min_x": 0.008501,
    "max_x": 1.099581,
    "min_y": 1.797981,
    "max_y": 4.124218,
    "horizontal": 0.5884363,
    "vertical": 0.2046999,
    "diagonal_intercept": 0.4474396,
    "diagonal_slope": 1.680319,
    "diagonal_vertical": 0.2075470,
    "diagonal_horizontal": 0.7867766,
}


# ---------------------------------------------------------------------------
# prediction_model simulation
# ---------------------------------------------------------------------------

def _predict_one(cr: float, cpd: float, params: list[float]) -> int:
    diag_intercept, diag_slope, diag_inter_v, horizontal_lim, vertical_lim = params
    if cr <= diag_inter_v:
        theo_cpd = horizontal_lim
    elif cr >= vertical_lim:
        theo_cpd = 10_000.0
    else:
        theo_cpd = cr * diag_slope + diag_intercept
    return 0 if theo_cpd - cpd > 0.0 else 1


def _select_n_dom(group: pd.DataFrame, params: list[float]) -> int:
    """Simulate prediction_model() + n_dom reverse scan. Returns selected n_dom.

    Mirrors Rust exactly:
      - One prediction per domain count level (first candidate at each level, skip nd=1)
      - predictions_rev scan: n_dom = i+1 where i is 0-based index of first pred=0
      - Fallback: n_dom = len(predictions) + 1 = max_nd
    """
    rows = group.sort_values("num_domains", ascending=False).reset_index(drop=True)
    if len(rows) <= 1:
        return 1  # Rust: predictions_rev.len()=0, n_dom defaults to 0+1=1
    # Skip the last element (nd=1), matching Rust's relevant_measure[..last]
    working = rows.iloc[:-1]
    # One prediction per domain count level — first occurrence (matches prediction_model)
    predictions: list[int] = []
    seen: set[int] = set()
    for _, r in working.iterrows():
        nd = int(r["num_domains"])
        if nd not in seen:
            seen.add(nd)
            predictions.append(_predict_one(float(r["max_cr"]), float(r["mean_density"]), params))
    # predictions is descending nd order; reversed → ascending nd order
    for i, pred in enumerate(reversed(predictions)):
        if pred == 0:
            return i + 1  # Rust: n_dom = i + 1
    return len(predictions) + 1  # Rust fallback: predictions_rev.len() + 1


def _make_prediction_arrays(
    training_df: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pre-process training data for vectorized loss evaluation.

    Returns (cr, cpd, chain_boundaries, oracle_s_per_chain, best_s_per_level) where:
      - cr[i], cpd[i]: features for the i-th (chain, nd-level) pair, ascending nd, excl. nd=1
      - chain_boundaries: indices into cr/cpd/best_s_per_level where each chain starts/ends
      - oracle_s_per_chain[j]: best S achievable for chain j (0 if chain not in training)
      - best_s_per_level[i]: best S at the i-th (chain, nd-level) pair
    Loss = mean(oracle_S - selected_S) — minimise quality gap.
    """
    # Per (chain, nd): best S across all candidates at that nd level
    best_s_per_chain_nd = (
        training_df.groupby(["chain_id", "num_domains"])["S"].max()
    )
    # Oracle S per chain
    oracle_s_per_chain_map = training_df.groupby("chain_id")["S"].max()

    # Dedup: one representative row per (chain_id, num_domains) for cr/cpd features
    sorted_df = training_df.sort_values(
        ["chain_id", "num_domains"], ascending=[True, False]
    )
    deduped = (
        sorted_df
        .groupby(["chain_id", "num_domains"], sort=False)
        .first()
        .reset_index()
    )

    # Remove the minimum nd per chain (nd=1, excluded in prediction_model)
    min_nd = deduped.groupby("chain_id")["num_domains"].transform("min")
    deduped = deduped[deduped["num_domains"] > min_nd].copy()

    # Sort: chain_id asc, num_domains asc
    deduped = deduped.sort_values(["chain_id", "num_domains"]).reset_index(drop=True)

    # Attach best_S per (chain, nd) level via merge
    best_s_df = best_s_per_chain_nd.reset_index()
    best_s_df.columns = ["chain_id", "num_domains", "best_S"]
    deduped = deduped.merge(best_s_df, on=["chain_id", "num_domains"], how="left")

    # Build chain boundary slices
    chain_ids_arr = deduped["chain_id"].values
    chain_boundaries = np.where(
        np.concatenate(([True], chain_ids_arr[1:] != chain_ids_arr[:-1], [True]))
    )[0]
    unique_chains = [chain_ids_arr[i] for i in chain_boundaries[:-1]]

    oracle_s_arr = np.array(
        [oracle_s_per_chain_map.get(cid, 0.0) for cid in unique_chains], dtype=np.float64
    )
    best_s_per_level = deduped["best_S"].values.astype(np.float64)
    cr = deduped["max_cr"].values.astype(np.float64)
    cpd = deduped["mean_density"].values.astype(np.float64)

    return cr, cpd, chain_boundaries, oracle_s_arr, best_s_per_level


def _prediction_loss_vec(
    params: np.ndarray,
    cr: np.ndarray,
    cpd: np.ndarray,
    chain_boundaries: np.ndarray,
    oracle_s_arr: np.ndarray,
    best_s_per_level: np.ndarray,
) -> float:
    """Fully vectorized S-quality loss: mean(oracle_S - S_at_selected_nd) across chains.

    Selection logic mirrors Rust prediction_model reverse scan exactly:
      - pred0[i]=True → scan stops at ascending-nd position i → n_dom = i+1
      - n_dom=1 (i=0) → single-domain partition (S approximated as 0)
      - n_dom=k>1 → best_S at nd=k = best_s_per_level[start + i - 1]
      - all pred0 False → fallback n_dom=max_nd = best_s_per_level[end-1]
    """
    di, ds, dv, hl, vl = params
    theo = np.where(cr <= dv, hl, np.where(cr >= vl, 1e4, cr * ds + di))
    pred0 = theo > cpd  # True = Rust pred 0 (good zone) → scan stops

    starts = chain_boundaries[:-1]
    ends = chain_boundaries[1:]
    n_per_chain = ends - starts

    # For each row, assign a score for "first True detection":
    # score[i] = global position i if pred0[i] else large value
    score = np.where(pred0, np.arange(len(pred0), dtype=np.int64), len(pred0))

    # Per chain: global index of first True (or len(pred0) if none)
    first_true_global = np.minimum.reduceat(score, starts)
    first_true_local = first_true_global - starts  # local index within chain

    has_true = first_true_local < n_per_chain

    # Build selected_s per chain:
    #   ~has_true                → fallback → best_s_per_level[end - 1]
    #   has_true & local == 0   → n_dom=1  → 0.0
    #   has_true & local > 0    → n_dom>1  → best_s_per_level[first_true_global - 1]
    fallback_s = best_s_per_level[np.clip(ends - 1, 0, len(best_s_per_level) - 1)]
    ndgt1_idx = np.clip(first_true_global - 1, 0, len(best_s_per_level) - 1)
    ndgt1_s = best_s_per_level[ndgt1_idx]

    selected_s = np.where(
        ~has_true,
        fallback_s,
        np.where(first_true_local == 0, 0.0, ndgt1_s),
    )

    gap = np.maximum(0.0, oracle_s_arr - selected_s)
    return float(gap.mean())


def fit_prediction_model(training_df: pd.DataFrame) -> dict[str, float]:
    log.info("Fitting prediction_model() ...")

    cr, cpd, chain_boundaries, oracle_s_arr, best_s_per_level = _make_prediction_arrays(training_df)
    n_chains = len(chain_boundaries) - 1
    log.info("  %d chains, %d (chain, nd-level) pairs after dedup", n_chains, len(cr))

    args = (cr, cpd, chain_boundaries, oracle_s_arr, best_s_per_level)

    x0 = np.array([
        BASELINE_PREDICTION["diag_intercept"],
        BASELINE_PREDICTION["diag_slope"],
        BASELINE_PREDICTION["diag_inter_v"],
        BASELINE_PREDICTION["horizontal_lim"],
        BASELINE_PREDICTION["vertical_lim"],
    ])
    baseline_loss = _prediction_loss_vec(x0, *args)
    log.info("  Baseline S-gap loss: %.4f", baseline_loss)

    result = minimize(
        lambda p: _prediction_loss_vec(p, *args),
        x0,
        method="Nelder-Mead",
        options={"maxiter": 20_000, "xatol": 1e-7, "fatol": 1e-7, "disp": True},
    )

    fitted = result.x
    fitted_loss = _prediction_loss_vec(fitted, *args)
    log.info(
        "  Fitted S-gap loss: %.4f  delta=%.4f",
        fitted_loss,
        baseline_loss - fitted_loss,
    )

    keys = ["diag_intercept", "diag_slope", "diag_inter_v", "horizontal_lim", "vertical_lim"]
    return dict(zip(keys, fitted.tolist()))


# ---------------------------------------------------------------------------
# distance_model simulation
# ---------------------------------------------------------------------------

# Cross-product reference point hardcoded in distance_model.rs (must match Rust exactly)
_DIAG_BX: float = 0.3288426


def _distance_score(cr: float, cpd: float, params: list[float]) -> float:
    """Signed distance from quality boundary — exact mirror of distance_model.rs.

    Returns positive in the good zone (ZONE 2/3/4), negative in bad zones.
    """
    (
        min_x, max_x, min_y, max_y,
        horizontal, vertical,
        diagonal_intercept, diagonal_slope, diagonal_vertical, diagonal_horizontal,
    ) = params

    mx = (cr - min_x) / (max_x - min_x + 1e-12)
    my = (cpd - min_y) / (max_y - min_y + 1e-12)

    dist_diag = (
        abs(diagonal_slope * mx - my + diagonal_intercept)
        / math.sqrt(1.0 + diagonal_slope * diagonal_slope)
    )
    dist_horiz = abs(my - horizontal)
    dist_vert = abs(mx - vertical)

    # Side of diagonal line through A=(0, diagonal_intercept) and B=(_DIAG_BX, 1.0)
    # d < 0 → below/right of diagonal (bad side)
    d = _DIAG_BX * (my - diagonal_intercept) - (1.0 - diagonal_intercept) * mx

    if my > horizontal and mx < vertical and d < 0.0:
        # ZONE 1: above horizontal, left of vertical, below diagonal
        return -min(dist_diag, dist_horiz, dist_vert)

    if (
        (my < diagonal_horizontal and my > horizontal and d < 0.0)
        or (mx > diagonal_vertical and mx < vertical and d < 0.0)
        or (my < diagonal_horizontal and mx > diagonal_vertical and d < 0.0)
    ):
        # ZONE 7, 8, 9
        return -dist_diag

    if my < horizontal:
        # ZONE 6
        return -dist_horiz

    if mx > vertical:
        # ZONE 5
        return -dist_vert

    # ZONE 2, 3, 4 — good zone
    return min(dist_diag, dist_horiz, dist_vert)


def _distance_score_vec(
    mx: np.ndarray,
    my: np.ndarray,
    horizontal: float,
    vertical: float,
    diagonal_intercept: float,
    diagonal_slope: float,
    diagonal_vertical: float,
    diagonal_horizontal: float,
) -> np.ndarray:
    """Vectorized zone classification on already-normalized (mx, my) arrays."""
    denom = math.sqrt(1.0 + diagonal_slope * diagonal_slope)
    dist_diag = np.abs(diagonal_slope * mx - my + diagonal_intercept) / denom
    dist_horiz = np.abs(my - horizontal)
    dist_vert = np.abs(mx - vertical)

    d = _DIAG_BX * (my - diagonal_intercept) - (1.0 - diagonal_intercept) * mx

    z1 = (my > horizontal) & (mx < vertical) & (d < 0.0)
    z789 = (d < 0.0) & (
        ((my < diagonal_horizontal) & (my > horizontal))
        | ((mx > diagonal_vertical) & (mx < vertical))
        | ((my < diagonal_horizontal) & (mx > diagonal_vertical))
    )
    z6 = my < horizontal
    z5 = mx > vertical

    min3 = np.minimum(np.minimum(dist_diag, dist_horiz), dist_vert)
    return np.where(z1, -min3,
           np.where(z789, -dist_diag,
           np.where(z6, -dist_horiz,
           np.where(z5, -dist_vert,
           min3))))


def _make_distance_arrays(
    training_df: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """Pre-process distance model data into flat arrays for vectorized loss.

    Returns (cr_all, cpd_all, is_oracle_all, group_starts, group_labels) or None if
    no valid groups exist.
    """
    cr_parts: list[np.ndarray] = []
    cpd_parts: list[np.ndarray] = []
    oracle_parts: list[np.ndarray] = []
    group_sizes: list[int] = []

    for _, grp in training_df.groupby(["chain_id", "num_domains"]):
        if len(grp) > 1 and grp["is_oracle_s"].any():
            cr_parts.append(grp["max_cr"].values.astype(np.float64))
            cpd_parts.append(grp["mean_density"].values.astype(np.float64))
            oracle_parts.append(grp["is_oracle_s"].values.astype(bool))
            group_sizes.append(len(grp))

    if not group_sizes:
        return None

    cr_all = np.concatenate(cr_parts)
    cpd_all = np.concatenate(cpd_parts)
    is_oracle_all = np.concatenate(oracle_parts)
    group_starts = np.concatenate([[0], np.cumsum(group_sizes[:-1])]).astype(np.int64)
    group_labels = np.repeat(np.arange(len(group_sizes), dtype=np.int32), group_sizes)

    return cr_all, cpd_all, is_oracle_all, group_starts, group_labels


def _distance_loss_vec(
    params: np.ndarray,
    cr_all: np.ndarray,
    cpd_all: np.ndarray,
    is_oracle_all: np.ndarray,
    group_starts: np.ndarray,
    group_labels: np.ndarray,
) -> float:
    """Fully vectorized distance loss: mean rank of oracle candidate in each group."""
    min_x, max_x, min_y, max_y, horizontal, vertical, diag_int, diag_slope, diag_vert, diag_horiz = params

    mx = (cr_all - min_x) / (max_x - min_x + 1e-12)
    my = (cpd_all - min_y) / (max_y - min_y + 1e-12)

    scores = _distance_score_vec(mx, my, horizontal, vertical, diag_int, diag_slope, diag_vert, diag_horiz)

    # Per-group oracle score (max score among oracle rows in each group)
    oracle_scores_masked = np.where(is_oracle_all, scores, -np.inf)
    oracle_per_group = np.maximum.reduceat(oracle_scores_masked, group_starts)

    # Rank = number of rows in the same group with score > oracle score
    oracle_broadcast = oracle_per_group[group_labels]
    total_rank = float(np.sum(scores > oracle_broadcast))
    n_groups = len(group_starts)
    return total_rank / max(n_groups, 1)


def fit_distance_model(training_df: pd.DataFrame) -> dict[str, float]:
    log.info("Fitting distance_model() ...")

    dist_arrays = _make_distance_arrays(training_df)
    if dist_arrays is None:
        log.warning("  No groups with multiple candidates — distance model fitting skipped")
        return dict(BASELINE_DISTANCE)

    cr_all, cpd_all, is_oracle_all, group_starts, group_labels = dist_arrays
    n_groups = len(group_starts)
    log.info("  %d groups, %d total rows", n_groups, len(cr_all))

    x0 = np.array([
        BASELINE_DISTANCE["min_x"],
        BASELINE_DISTANCE["max_x"],
        BASELINE_DISTANCE["min_y"],
        BASELINE_DISTANCE["max_y"],
        BASELINE_DISTANCE["horizontal"],
        BASELINE_DISTANCE["vertical"],
        BASELINE_DISTANCE["diagonal_intercept"],
        BASELINE_DISTANCE["diagonal_slope"],
        BASELINE_DISTANCE["diagonal_vertical"],
        BASELINE_DISTANCE["diagonal_horizontal"],
    ])

    baseline_loss = _distance_loss_vec(x0, cr_all, cpd_all, is_oracle_all, group_starts, group_labels)
    log.info("  Baseline mean oracle rank: %.4f", baseline_loss)

    result = minimize(
        lambda p: _distance_loss_vec(p, cr_all, cpd_all, is_oracle_all, group_starts, group_labels),
        x0,
        method="Nelder-Mead",
        options={"maxiter": 50_000, "xatol": 1e-6, "fatol": 1e-6, "disp": True},
    )

    fitted = result.x
    fitted_loss = _distance_loss_vec(fitted, cr_all, cpd_all, is_oracle_all, group_starts, group_labels)
    log.info("  Fitted mean oracle rank: %.4f  delta=%.4f", fitted_loss, baseline_loss - fitted_loss)

    keys = [
        "min_x", "max_x", "min_y", "max_y",
        "horizontal", "vertical",
        "diagonal_intercept", "diagonal_slope", "diagonal_vertical", "diagonal_horizontal",
    ]
    return dict(zip(keys, fitted.tolist()))


# ---------------------------------------------------------------------------
# Rust snippet generation
# ---------------------------------------------------------------------------

def _rust_prediction_snippet(params: dict[str, float]) -> str:
    return (
        "// prediction_model() constants — paste into sword2-lib/src/sword/mod.rs:392-396\n"
        f"let diag_intercept: f64 = {params['diag_intercept']};\n"
        f"let diag_slope: f64 = {params['diag_slope']};\n"
        f"let diag_inter_v: f64 = {params['diag_inter_v']};\n"
        f"let horizontal_lim: f64 = {params['horizontal_lim']};\n"
        f"let vertical_lim: f64 = {params['vertical_lim']};\n"
    )


def _rust_distance_snippet(params: dict[str, float]) -> str:
    return (
        "// distance_model() constants — paste into sword2-lib/src/sword/distance_model.rs:19-34\n"
        f"let min_x: f64 = {params['min_x']};\n"
        f"let max_x: f64 = {params['max_x']};\n"
        f"let min_y: f64 = {params['min_y']};\n"
        f"let max_y: f64 = {params['max_y']};\n"
        f"let horizontal: f64 = {params['horizontal']};\n"
        f"let vertical: f64 = {params['vertical']};\n"
        f"let diagonal_intercept: f64 = {params['diagonal_intercept']};\n"
        f"let diagonal_slope: f64 = {params['diagonal_slope']};\n"
        f"let diagonal_vertical: f64 = {params['diagonal_vertical']};\n"
        f"let diagonal_horizontal: f64 = {params['diagonal_horizontal']};\n"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Refit SWORD2 selection model constants")
    parser.add_argument("--training", type=Path, default=Path("benchmark/data/training_table.csv"))
    parser.add_argument("--output", type=Path, default=Path("benchmark/data/fitted_constants.json"))
    parser.add_argument("--skip-distance", action="store_true",
                        help="Only refit prediction_model (faster)")
    args = parser.parse_args()

    log.info("Loading training table from %s", args.training)
    df = pd.read_csv(args.training)
    log.info("  %d rows, %d chains", len(df), df["chain_id"].nunique())

    fitted_pred = fit_prediction_model(df)

    if args.skip_distance:
        fitted_dist = dict(BASELINE_DISTANCE)
    else:
        fitted_dist = fit_distance_model(df)

    output = {
        "prediction_model": {
            "baseline": BASELINE_PREDICTION,
            "fitted": fitted_pred,
        },
        "distance_model": {
            "baseline": BASELINE_DISTANCE,
            "fitted": fitted_dist,
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2))
    log.info("Saved fitted constants to %s", args.output)

    print("\n" + "=" * 60)
    print(_rust_prediction_snippet(fitted_pred))
    print(_rust_distance_snippet(fitted_dist))
    print("=" * 60)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
