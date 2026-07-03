#!/usr/bin/env python3
"""Diagnose *why* SWORD2's rank-1 domain selection loses ~0.12 NDO to its own
oracle, before anyone builds a third reranker.

Everything the benchmark needs already exists:

* ``results_bypass/scores.csv`` — every SWORD2 candidate (variant ``optimal`` =
  rank-1, ``alternative`` = the rest) scored against CATH ground truth, plus the
  Merizo / Chainsaw rows. Tier-1 analyses use only this.
* ``data/cath663_candidate_features.csv`` — per-candidate features emitted by the
  ``SWORD2_DUMP_CANDIDATES`` path (see ``dump_cath663_features.py``). Joined to
  ``scores.csv`` by a numbering-invariant chopping *signature* it gives, per
  candidate, the features a selector could use *and* the true NDO it would score.

The module answers five questions (see ``main`` for the printed report):

  Tier 1 (scores only)
    1. Gap attribution — how much of oracle−rank1 NDO is a wrong *domain count*
       vs a wrong *partition at the chosen count*.
    2. Cheap baseline selectors — does any trivial rule already beat rank-1?
  Tier 2 (features ⋈ scores)
    3. Feature-blind ceiling — best top-1 NDO a strong model (gradient-boosted
       trees, cross-validated by chain) can reach with the current features.
       This is the verdict on "is it the features or the model?".
    4. Per-feature signal — which features separate the best candidate; does
       ``energy_z`` earn its ~2x runtime cost.
    5. Pairwise-reranker autopsy — reproduce the trained reranker's top-1 NDO and
       explain how 75% pairwise accuracy still regressed top-1.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
DEFAULT_SCORES = REPO / "benchmark/results_bypass/scores.csv"
DEFAULT_FEATURES = REPO / "benchmark/data/cath663_candidate_features.csv"
DEFAULT_WEIGHTS = REPO / "benchmark/data/pairwise_reranker_weights.json"

# The corrected candidate feature set (after the field-offset fix in
# sword2-lib/src/sword/mod.rs: density_min is fields[5], mean_density fields[6]).
FEATURE_COLS = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "energy_z",
    "modal_count_distance",
]


# --------------------------------------------------------------------------- #
# Join: numbering-invariant chopping signature
# --------------------------------------------------------------------------- #
def signature(chopping: str) -> tuple | None:
    """Canonical key shared by a dump ``delineation`` (0-based, space-separated
    domains / ';'-separated segments) and a ``scores.csv`` ``pred_chopping``
    (0-based, comma-separated domains / '_'-separated segments).

    Returns a tuple over domains (ordered by first residue) of the tuple of
    segment *lengths* (ordered by start). Lengths are invariant to the two 0-based
    index spaces, so it matches candidates without reconciling numbering.
    """
    chopping = chopping.strip().strip('"')
    if not chopping:
        return None
    domains = []
    for dom in re.split(r"[ ,]+", chopping):
        dom = dom.strip()
        if not dom:
            continue
        segs = []
        for seg in re.split(r"[;_]", dom):
            seg = seg.strip()
            if not seg:
                continue
            match = re.match(r"^(-?\d+)(?:-(-?\d+))?$", seg)
            if not match:
                raise ValueError(f"bad segment {seg!r} in {chopping!r}")
            lo = int(match.group(1))
            hi = int(match.group(2)) if match.group(2) else lo
            lo, hi = min(lo, hi), max(lo, hi)
            segs.append((lo, hi - lo + 1))
        if not segs:
            continue
        segs.sort()
        domains.append((segs[0][0], tuple(length for _, length in segs)))
    domains.sort()
    return tuple(seg_lengths for _, seg_lengths in domains)


# --------------------------------------------------------------------------- #
# Data model
# --------------------------------------------------------------------------- #
@dataclass
class Candidate:
    entry_id: str
    variant: str  # "optimal" (rank-1) or "alternative"
    n_pred: int
    n_true: int
    n_residues: int
    ndo: float
    iou: float
    boundary_f1_10: float
    features: dict[str, float] | None = None  # corrected dump features (or None)
    derived: dict[str, float] = field(default_factory=dict)


def _load_feature_index(features_path: Path) -> dict[tuple[str, tuple], dict[str, float]]:
    index: dict[tuple[str, tuple], dict[str, float]] = {}
    for row in csv.DictReader(features_path.open()):
        key = (row["entry_id"], signature(row["delineation"]))
        # First write wins (the <0.1% signature collisions are near-identical).
        index.setdefault(key, {col: float(row[col]) for col in FEATURE_COLS})
    return index


def load_candidates(
    scores_path: Path = DEFAULT_SCORES,
    features_path: Path | None = DEFAULT_FEATURES,
) -> dict[str, list[Candidate]]:
    """Return SWORD2 candidates grouped by entry, features joined when available."""
    feat_index = _load_feature_index(features_path) if features_path else {}
    by_entry: dict[str, list[Candidate]] = defaultdict(list)
    for row in csv.DictReader(scores_path.open()):
        if row["tool"] != "sword2-rust" or row["variant"] not in ("optimal", "alternative"):
            continue
        feats = feat_index.get((row["entry_id"], signature(row["pred_chopping"])))
        by_entry[row["entry_id"]].append(
            Candidate(
                entry_id=row["entry_id"],
                variant=row["variant"],
                n_pred=int(row["n_pred_domains"]),
                n_true=int(row["n_true_domains"]),
                n_residues=int(row["n_residues"]),
                ndo=float(row["ndo"]),
                iou=float(row["iou"]),
                boundary_f1_10=float(row["boundary_f1_10"]),
                features=dict(feats) if feats else None,
            )
        )
    _add_derived_features(by_entry)
    return dict(by_entry)


def _add_derived_features(by_entry: dict[str, list[Candidate]]) -> None:
    """Cheap, inference-available features a count selector might want."""
    for cands in by_entry.values():
        crs = [c.features["max_cr"] for c in cands if c.features]
        max_cr_chain = max(crs) if crs else 1.0
        for c in cands:
            if not c.features:
                continue
            c.derived = {
                "n_residues": float(c.n_residues),
                "log_n_residues": float(np.log(max(c.n_residues, 1))),
                "rel_max_cr": c.features["max_cr"] / max_cr_chain if max_cr_chain else 0.0,
            }


def feature_coverage(by_entry: dict[str, list[Candidate]]) -> tuple[int, int]:
    total = sum(len(v) for v in by_entry.values())
    matched = sum(1 for v in by_entry.values() for c in v if c.features)
    return matched, total


# --------------------------------------------------------------------------- #
# Selector evaluation (the common currency: mean top-1 NDO)
# --------------------------------------------------------------------------- #
def rank1(cands: list[Candidate]) -> Candidate:
    for c in cands:
        if c.variant == "optimal":
            return c
    return cands[0]


def per_entry_pick_ndo(by_entry, score_fn, restrict=None) -> dict[str, float]:
    """argmax ``score_fn`` per entry -> {entry: true NDO of the pick}.

    ``restrict(cands)`` may narrow the candidate pool first (returns a subset;
    empty subset falls back to the full pool).
    """
    out = {}
    for eid, cands in by_entry.items():
        pool = cands
        if restrict is not None:
            sub = restrict(cands)
            pool = sub if sub else cands
        out[eid] = max(pool, key=score_fn).ndo
    return out


def bootstrap_ci(values: list[float], n_boot: int = 2000, seed: int = 0) -> tuple[float, float]:
    if not values:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    arr = np.asarray(values, dtype=float)
    means = arr[rng.integers(0, len(arr), size=(n_boot, len(arr)))].mean(axis=1)
    return (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)))


def _mean_ci(per_entry: dict[str, float]) -> dict:
    vals = list(per_entry.values())
    lo, hi = bootstrap_ci(vals)
    return {"mean": float(np.mean(vals)), "ci": (lo, hi), "n": len(vals)}


# --------------------------------------------------------------------------- #
# Tier 1 — gap attribution
# --------------------------------------------------------------------------- #
def gap_attribution(by_entry: dict[str, list[Candidate]]) -> dict:
    """Split oracle−rank1 NDO into a within-count and a count component.

    within_count = best NDO at rank-1's own domain count − rank-1 NDO
    count        = oracle NDO − best NDO at rank-1's own domain count
    (both >= 0, and they sum to oracle − rank1)
    """
    rows = []
    for eid, cands in by_entry.items():
        r1 = rank1(cands)
        oracle = max(c.ndo for c in cands)
        best_at_pred = max(c.ndo for c in cands if c.n_pred == r1.n_pred)
        at_true = [c.ndo for c in cands if c.n_pred == r1.n_true]
        best_at_true = max(at_true) if at_true else float("nan")
        best_count = max(cands, key=lambda c: c.ndo).n_pred
        rows.append(
            {
                "entry_id": eid,
                "rank1": r1.ndo,
                "oracle": oracle,
                "best_at_pred_count": best_at_pred,
                "best_at_true_count": best_at_true,
                "within_count": best_at_pred - r1.ndo,
                "count": oracle - best_at_pred,
                "rank1_count": r1.n_pred,
                "true_count": r1.n_true,
                "best_ndo_count": best_count,
                "has_true_count_candidate": bool(at_true),
            }
        )
    mean = lambda k: float(np.mean([r[k] for r in rows]))
    n = len(rows)
    return {
        "n": n,
        "rank1_ndo": mean("rank1"),
        "oracle_ndo": mean("oracle"),
        "total_gap": mean("oracle") - mean("rank1"),
        "within_count_component": mean("within_count"),
        "count_component": mean("count"),
        "within_count_share": mean("within_count") / (mean("oracle") - mean("rank1")),
        "count_share": mean("count") / (mean("oracle") - mean("rank1")),
        "frac_rank1_count_is_best": float(np.mean([r["rank1_count"] == r["best_ndo_count"] for r in rows])),
        "frac_rank1_count_eq_true": float(np.mean([r["rank1_count"] == r["true_count"] for r in rows])),
        "frac_true_count_reachable": float(np.mean([r["has_true_count_candidate"] for r in rows])),
        "rows": rows,
    }


def count_bias(by_entry: dict[str, list[Candidate]]) -> dict:
    """Distribution of rank-1 (predicted − true) domain count."""
    deltas = [rank1(c).n_pred - rank1(c).n_true for c in by_entry.values()]
    hist = Counter(deltas)
    return {
        "mean_delta": float(np.mean(deltas)),
        "frac_under": float(np.mean([d < 0 for d in deltas])),
        "frac_exact": float(np.mean([d == 0 for d in deltas])),
        "frac_over": float(np.mean([d > 0 for d in deltas])),
        "histogram": dict(sorted(hist.items())),
    }


# --------------------------------------------------------------------------- #
# Tier 1 — cheap baseline selectors
# --------------------------------------------------------------------------- #
def _modal_count(cands: list[Candidate]) -> int:
    counts = Counter(c.n_pred for c in cands)
    top = max(counts.values())
    return min(k for k, v in counts.items() if v == top)


def baseline_selectors(by_entry: dict[str, list[Candidate]]) -> dict[str, dict]:
    have_feats = all(c.features for v in by_entry.values() for c in v)

    def at_modal(cands):
        m = _modal_count(cands)
        return [c for c in cands if c.n_pred == m]

    def at_true(cands):
        return [c for c in cands if c.n_pred == c.n_true]

    selectors: dict[str, dict[str, float]] = {}
    selectors["rank1 (current)"] = {eid: rank1(c).ndo for eid, c in by_entry.items()}
    selectors["oracle (best NDO)"] = per_entry_pick_ndo(by_entry, lambda c: c.ndo)
    selectors["fewest domains"] = per_entry_pick_ndo(by_entry, lambda c: -c.n_pred)
    selectors["most domains"] = per_entry_pick_ndo(by_entry, lambda c: c.n_pred)
    # modal count, tie-broken toward SWORD2's own preference (highest max_cr if
    # features present, else fewest domains then first)
    tiebreak = (lambda c: c.features["max_cr"]) if have_feats else (lambda c: -c.n_pred)
    selectors["modal count"] = per_entry_pick_ndo(by_entry, tiebreak, restrict=at_modal)
    selectors["true count (perfect count)"] = per_entry_pick_ndo(by_entry, lambda c: c.ndo, restrict=at_true)
    if have_feats:
        selectors["highest max_cr"] = per_entry_pick_ndo(by_entry, lambda c: c.features["max_cr"])
        selectors["lowest density_min"] = per_entry_pick_ndo(by_entry, lambda c: -c.features["density_min"])
    return {name: _mean_ci(pe) for name, pe in selectors.items()}


def tool_comparison(scores_path: Path = DEFAULT_SCORES) -> dict:
    """Headline metrics per tool/variant straight from scores.csv."""
    agg: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in csv.DictReader(scores_path.open()):
        key = (row["tool"], row["variant"])
        for metric in ("ndo", "d_count_acc", "boundary_f1_10", "iou"):
            val = row[metric]
            if val not in ("", "nan"):
                agg[key][metric].append(float(val))
    out = {}
    for (tool, variant), metrics in agg.items():
        out[f"{tool}/{variant}"] = {m: float(np.mean(v)) for m, v in metrics.items() if v}
        out[f"{tool}/{variant}"]["n"] = len(next(iter(metrics.values())))
    return out


# --------------------------------------------------------------------------- #
# Tier 2 — feature matrix
# --------------------------------------------------------------------------- #
def build_matrix(by_entry, cols) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[Candidate]]:
    """Rows = candidates with features. Returns X, y=ndo, groups=entry, cands."""
    X, y, groups, cands = [], [], [], []
    for eid, clist in by_entry.items():
        for c in clist:
            if not c.features:
                continue
            vec = []
            for col in cols:
                if col in c.features:
                    vec.append(c.features[col])
                elif col in c.derived:
                    vec.append(c.derived[col])
                else:
                    raise KeyError(col)
            X.append(vec)
            y.append(c.ndo)
            groups.append(eid)
            cands.append(c)
    return np.asarray(X), np.asarray(y), np.asarray(groups), cands


def _grouped_topk_ndo(scores: np.ndarray, groups: np.ndarray, cands: list[Candidate]) -> dict[str, float]:
    """Per group, pick argmax score -> {entry: true NDO}."""
    best: dict[str, tuple[float, float]] = {}
    for s, g, c in zip(scores, groups, cands):
        if g not in best or s > best[g][0]:
            best[g] = (s, c.ndo)
    return {g: ndo for g, (_, ndo) in best.items()}


def _cv_ceiling(by_entry, cols, seed=0) -> dict:
    """Chain-cross-validated best top-1 NDO reachable with feature `cols`.

    Trains gradient-boosted trees (regressor on NDO and classifier on
    is-oracle-best) under GroupKFold(5) so no chain leaks train->test, then on
    held-out chains picks argmax prediction and scores its true NDO. The better
    of the two models is the ceiling for this feature set.
    """
    from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
    from sklearn.model_selection import GroupKFold

    X, y, groups, cands = build_matrix(by_entry, cols)
    best_ndo = defaultdict(float)
    for g, ndo in zip(groups, y):
        best_ndo[g] = max(best_ndo[g], ndo)
    is_best = np.array([1 if ndo >= best_ndo[g] - 1e-9 else 0 for g, ndo in zip(groups, y)])

    gkf = GroupKFold(n_splits=5)
    reg_pred = np.zeros(len(y))
    clf_pred = np.zeros(len(y))
    for tr, te in gkf.split(X, y, groups):
        reg = HistGradientBoostingRegressor(random_state=seed, max_iter=300, learning_rate=0.05)
        reg.fit(X[tr], y[tr])
        reg_pred[te] = reg.predict(X[te])
        clf = HistGradientBoostingClassifier(random_state=seed, max_iter=300, learning_rate=0.05)
        clf.fit(X[tr], is_best[tr])
        clf_pred[te] = clf.predict_proba(X[te])[:, 1]
    reg_stat = _mean_ci(_grouped_topk_ndo(reg_pred, groups, cands))
    clf_stat = _mean_ci(_grouped_topk_ndo(clf_pred, groups, cands))
    best = reg_stat if reg_stat["mean"] >= clf_stat["mean"] else clf_stat
    return {"regressor": reg_stat, "classifier": clf_stat, "ceiling": best}


def feature_blind_ceiling(by_entry, cols=None, augmented_cols=None, seed=0) -> dict:
    """Ceiling for the current and augmented feature sets + best single feature."""
    cols = cols or FEATURE_COLS
    result: dict = {"features": cols, "current": _cv_ceiling(by_entry, cols, seed)}
    if augmented_cols:
        result["augmented"] = _cv_ceiling(by_entry, augmented_cols, seed)
        result["augmented_features"] = augmented_cols

    single = {}
    for col in cols:
        pe_hi = per_entry_pick_ndo(by_entry, lambda c, k=col: _feat(c, k))
        pe_lo = per_entry_pick_ndo(by_entry, lambda c, k=col: -_feat(c, k))
        single[col] = max(float(np.mean(list(pe_hi.values()))), float(np.mean(list(pe_lo.values()))))
    result["best_single_feature"] = dict(sorted(single.items(), key=lambda kv: -kv[1]))
    return result


def ceiling_ablations(by_entry, seed=0) -> dict[str, dict]:
    """Ceiling for named feature subsets — isolates each feature group's value."""
    geometry = ["num_domains", "min_size", "max_cr", "density_min", "mean_density"]
    subsets = {
        "geometry only (5)": geometry,
        "geometry + energy_z": geometry + ["energy_z"],
        "current − energy_z": [c for c in FEATURE_COLS if c != "energy_z"],
        "current − coil": [c for c in FEATURE_COLS if c != "boundary_coil_fraction"],
        "current + chain length": FEATURE_COLS + ["n_residues", "log_n_residues"],
    }
    return {name: _cv_ceiling(by_entry, cols, seed)["ceiling"] for name, cols in subsets.items()}


def _feat(c: Candidate, key: str) -> float:
    if c.features and key in c.features:
        return c.features[key]
    return c.derived.get(key, 0.0)


# --------------------------------------------------------------------------- #
# Tier 2 — per-feature signal
# --------------------------------------------------------------------------- #
def per_feature_signal(by_entry, cols=None) -> dict:
    """AUC of each feature for 'is oracle-best candidate' + Spearman vs NDO."""
    from scipy.stats import spearmanr
    from sklearn.metrics import roc_auc_score

    cols = cols or FEATURE_COLS
    X, y, groups, cands = build_matrix(by_entry, cols)
    best_ndo = defaultdict(float)
    for g, ndo in zip(groups, y):
        best_ndo[g] = max(best_ndo[g], ndo)
    is_best = np.array([1 if ndo >= best_ndo[g] - 1e-9 else 0 for g, ndo in zip(groups, y)])

    out = {}
    for j, col in enumerate(cols):
        xj = X[:, j]
        if np.allclose(xj, xj[0]):
            out[col] = {"auc": 0.5, "abs_auc": 0.0, "spearman_ndo": 0.0, "note": "constant"}
            continue
        auc = roc_auc_score(is_best, xj)
        rho = spearmanr(xj, y).statistic
        out[col] = {"auc": float(auc), "abs_auc": abs(auc - 0.5), "spearman_ndo": float(rho)}
    return dict(sorted(out.items(), key=lambda kv: -kv[1]["abs_auc"]))


# --------------------------------------------------------------------------- #
# Tier 2 — pairwise reranker autopsy
# --------------------------------------------------------------------------- #
def _reranker_model_vectors(cands: list[Candidate]) -> np.ndarray:
    """Reconstruct the feature vector the trained reranker actually saw.

    Per reranker.rs the vector is
      [num_domains, min_size, max_cr, density_min, mean_density, coil, energy_z,
       modal_count_distance]
    but the field-offset bug meant density_min was always 0 and the 'mean_density'
    slot actually held the real density_min. The corrected dump exposes the real
    density_min, so we map: model_density_min=0, model_mean_density=real density_min.
    """
    vecs = []
    for c in cands:
        f = c.features
        vecs.append(
            [
                f["num_domains"],
                f["min_size"],
                f["max_cr"],
                0.0,  # density_min as the model saw it (dead constant)
                f["density_min"],  # 'mean_density' slot actually held real density_min
                f["boundary_coil_fraction"],
                f["energy_z"],
                f["modal_count_distance"],
            ]
        )
    return np.asarray(vecs)


def _zscore_per_chain(vecs: np.ndarray) -> np.ndarray:
    if len(vecs) < 2:
        return np.zeros_like(vecs)
    mean = vecs.mean(axis=0)
    std = vecs.std(axis=0) + 1e-8
    return (vecs - mean) / std


def pairwise_autopsy(by_entry, weights_path: Path = DEFAULT_WEIGHTS) -> dict:
    weights_json = json.loads(weights_path.read_text())
    w = np.asarray(weights_json["weights"])
    bias = float(weights_json.get("bias", 0.0))

    pick_ndo = {}
    count_deltas = []
    pair_total = pair_correct = 0
    margin_buckets = defaultdict(lambda: [0, 0])  # bucket -> [correct, total]
    decisive_correct = decisive_total = 0  # best vs 2nd-best per chain

    for eid, cands in by_entry.items():
        feats = [c for c in cands if c.features]
        if len(feats) < 2:
            pick_ndo[eid] = rank1(cands).ndo
            continue
        vecs = _zscore_per_chain(_reranker_model_vectors(feats))
        scores = vecs @ w + bias
        pick = int(np.argmax(scores))
        pick_ndo[eid] = feats[pick].ndo
        count_deltas.append(feats[pick].n_pred - feats[pick].n_true)

        ndos = np.array([c.ndo for c in feats])
        order = np.argsort(-ndos)
        # decisive pair: the true-best vs the runner-up
        if len(order) >= 2 and ndos[order[0]] > ndos[order[1]] + 1e-9:
            i, j = order[0], order[1]
            decisive_total += 1
            decisive_correct += int((scores[i] - scores[j]) > 0)
        for a in range(len(feats)):
            for b in range(a + 1, len(feats)):
                if abs(ndos[a] - ndos[b]) < 1e-9:
                    continue
                pair_total += 1
                correct = int(np.sign(scores[a] - scores[b]) == np.sign(ndos[a] - ndos[b]))
                pair_correct += correct
                bucket = min(int(abs(ndos[a] - ndos[b]) / 0.1), 4)  # 0-.1,.1-.2,.2-.3,.3-.4,.4+
                margin_buckets[bucket][0] += correct
                margin_buckets[bucket][1] += 1

    stat = _mean_ci(pick_ndo)
    labels = {0: "0.0-0.1", 1: "0.1-0.2", 2: "0.2-0.3", 3: "0.3-0.4", 4: "0.4+"}
    by_margin = {
        labels[k]: {"accuracy": c / t if t else float("nan"), "n": t}
        for k, (c, t) in sorted(margin_buckets.items())
    }
    return {
        "reranker_top1_ndo": stat,
        "pairwise_accuracy_overall": pair_correct / pair_total if pair_total else float("nan"),
        "pairwise_accuracy_decisive": decisive_correct / decisive_total if decisive_total else float("nan"),
        "pairwise_accuracy_by_margin": by_margin,
        "reranker_count_bias_mean": float(np.mean(count_deltas)) if count_deltas else float("nan"),
        "weights": dict(zip(weights_json["features"], weights_json["weights"])),
    }


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def run_all(scores_path=DEFAULT_SCORES, features_path=DEFAULT_FEATURES, weights_path=DEFAULT_WEIGHTS) -> dict:
    by_entry = load_candidates(scores_path, features_path)
    matched, total = feature_coverage(by_entry)
    augmented = FEATURE_COLS + ["n_residues", "log_n_residues", "rel_max_cr"]
    return {
        "coverage": {"matched": matched, "total": total, "pct": 100 * matched / total},
        "tools": tool_comparison(scores_path),
        "gap_attribution": gap_attribution(by_entry),
        "count_bias": count_bias(by_entry),
        "baselines": baseline_selectors(by_entry),
        "ceiling": feature_blind_ceiling(by_entry, FEATURE_COLS, augmented),
        "ceiling_ablations": ceiling_ablations(by_entry),
        "feature_signal": per_feature_signal(by_entry),
        "autopsy": pairwise_autopsy(by_entry, weights_path),
    }


def _fmt_ci(stat: dict) -> str:
    return f"{stat['mean']:.4f}  [{stat['ci'][0]:.4f}, {stat['ci'][1]:.4f}]"


def print_report(r: dict) -> None:
    p = print
    p("=" * 78)
    p("SWORD2 SELECTION DIAGNOSIS")
    p("=" * 78)
    cov = r["coverage"]
    p(f"\nfeature join coverage: {cov['matched']}/{cov['total']} ({cov['pct']:.2f}%)")

    p("\n-- Tool comparison (mean) " + "-" * 52)
    p(f"{'tool/variant':<28}{'NDO':>8}{'d_cnt':>8}{'bf1_10':>8}{'iou':>8}{'n':>6}")
    for name in ["sword2-rust/optimal", "sword2-rust/oracle", "sword2-rust/oracle_s",
                 "merizo/cuda", "chainsaw/single"]:
        if name in r["tools"]:
            t = r["tools"][name]
            p(f"{name:<28}{t.get('ndo',0):>8.4f}{t.get('d_count_acc',0):>8.4f}"
              f"{t.get('boundary_f1_10',0):>8.4f}{t.get('iou',0):>8.4f}{t.get('n',0):>6}")

    g = r["gap_attribution"]
    p("\n-- 1. Gap attribution " + "-" * 56)
    p(f"rank-1 NDO {g['rank1_ndo']:.4f}  ->  oracle NDO {g['oracle_ndo']:.4f}   (gap {g['total_gap']:.4f})")
    p(f"  within-count component: {g['within_count_component']:.4f}  ({100*g['within_count_share']:.1f}% of gap)")
    p(f"  count component:        {g['count_component']:.4f}  ({100*g['count_share']:.1f}% of gap)")
    p(f"  rank-1 count == best-NDO count: {100*g['frac_rank1_count_is_best']:.1f}% of chains")
    p(f"  rank-1 count == true count:     {100*g['frac_rank1_count_eq_true']:.1f}% of chains")
    p(f"  a true-count candidate exists:  {100*g['frac_true_count_reachable']:.1f}% of chains")

    cb = r["count_bias"]
    p("\n-- Count bias (rank-1 predicted - true domains) " + "-" * 30)
    p(f"  mean delta {cb['mean_delta']:+.3f}   under {100*cb['frac_under']:.1f}%  "
      f"exact {100*cb['frac_exact']:.1f}%  over {100*cb['frac_over']:.1f}%")
    p(f"  histogram: {cb['histogram']}")

    p("\n-- 2. Baseline selectors (mean top-1 NDO) " + "-" * 36)
    for name, stat in sorted(r["baselines"].items(), key=lambda kv: -kv[1]["mean"]):
        p(f"  {name:<28}{_fmt_ci(stat)}")

    c = r["ceiling"]
    p("\n-- 3. Feature-blind ceiling (chain-CV top-1 NDO) " + "-" * 29)
    p(f"  current {len(c['features'])} features:")
    p(f"     regressor : {_fmt_ci(c['current']['regressor'])}")
    p(f"     classifier: {_fmt_ci(c['current']['classifier'])}")
    p(f"     CEILING   : {_fmt_ci(c['current']['ceiling'])}")
    if "augmented" in c:
        p(f"  augmented (+{', '.join(c['augmented_features'][len(FEATURE_COLS):])}):")
        p(f"     CEILING   : {_fmt_ci(c['augmented']['ceiling'])}")
    p("  best single-feature selectors:")
    for k, v in list(c["best_single_feature"].items())[:4]:
        p(f"     {k:<26}{v:.4f}")
    p("  ablations (ceiling by feature subset):")
    for name, stat in r["ceiling_ablations"].items():
        p(f"     {name:<26}{_fmt_ci(stat)}")

    p("\n-- 4. Per-feature signal (AUC for is-oracle-best; Spearman vs NDO) " + "-" * 12)
    p(f"  {'feature':<26}{'AUC':>8}{'|AUC-.5|':>10}{'rho(NDO)':>10}")
    for feat, s in r["feature_signal"].items():
        note = f"  {s.get('note','')}" if s.get("note") else ""
        p(f"  {feat:<26}{s['auc']:>8.3f}{s['abs_auc']:>10.3f}{s['spearman_ndo']:>10.3f}{note}")

    a = r["autopsy"]
    p("\n-- 5. Pairwise reranker autopsy " + "-" * 46)
    p(f"  reranker top-1 NDO: {_fmt_ci(a['reranker_top1_ndo'])}   (vs rank-1 {g['rank1_ndo']:.4f})")
    p(f"  pairwise accuracy overall : {a['pairwise_accuracy_overall']:.3f}")
    p(f"  pairwise accuracy decisive: {a['pairwise_accuracy_decisive']:.3f}   (true-best vs runner-up)")
    p(f"  reranker count bias mean  : {a['reranker_count_bias_mean']:+.3f}")
    p("  pairwise accuracy by NDO margin:")
    for bucket, v in a["pairwise_accuracy_by_margin"].items():
        p(f"     margin {bucket:<8} acc {v['accuracy']:.3f}  (n={v['n']})")
    p("  trained weights:")
    for feat, wt in a["weights"].items():
        p(f"     {feat:<26}{wt:+.4f}")
    p("\n" + "=" * 78)


def make_figures(results: dict, outdir: Path) -> list[Path]:
    """Three figures that carry the diagnosis. Clean static bars, value-labeled."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    ink, accent, muted, warn = "#2f3640", "#3b6ea5", "#9aa5b1", "#c0603a"
    paths: list[Path] = []
    g = results["gap_attribution"]
    ceil = results["ceiling"]
    abl = results["ceiling_ablations"]
    base = results["baselines"]
    merizo = results["tools"]["merizo/cuda"]["ndo"]

    # Fig 1 — the NDO ladder: how far each lever can carry selection.
    rungs = [
        ("rank-1 (ships today)", g["rank1_ndo"], muted),
        ("ceiling · current features", ceil["current"]["ceiling"]["mean"], accent),
        ("ceiling · + chain length", abl["current + chain length"]["mean"], accent),
        ("ceiling · augmented", ceil["augmented"]["ceiling"]["mean"], accent),
        ("perfect count selector", base["true count (perfect count)"]["mean"], ink),
        ("oracle (best candidate)", g["oracle_ndo"], ink),
    ]
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ys = range(len(rungs))
    ax.barh([y for y in ys], [v for _, v, _ in rungs], color=[c for *_, c in rungs], height=0.62)
    ax.set_yticks(list(ys))
    ax.set_yticklabels([n for n, *_ in rungs])
    ax.invert_yaxis()
    ax.set_xlim(0.70, 0.92)
    ax.set_xlabel("mean top-1 NDO (CATH-663)")
    ax.axvline(merizo, color=warn, ls="--", lw=1.2)
    ax.text(merizo + 0.001, 5.55, f"Merizo {merizo:.3f}", color=warn, fontsize=8,
            ha="center", va="center", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none"))
    for y, (_, v, _) in zip(ys, rungs):
        ax.text(v + 0.002, y, f"{v:.3f}", va="center", fontsize=9, color=ink)
    ax.set_title("Selection NDO ladder: current features barely clear rank-1;\n"
                 "the gap to oracle needs count-calibration, not another reranker", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p1 = outdir / "fig1_ndo_ladder.png"
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    paths.append(p1)

    # Fig 2 — marginal ceiling value of each add-on feature.
    cur = abl["current − energy_z"]["mean"]
    cur_coil = abl["current − coil"]["mean"]
    full = ceil["current"]["ceiling"]["mean"]
    marg = [
        ("energy_z\n(≈2x runtime)", full - cur, warn),
        ("boundary coil", full - cur_coil, muted),
        ("chain length\n(free)", abl["current + chain length"]["mean"] - full, accent),
    ]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.bar([n for n, *_ in marg], [v for _, v, _ in marg], color=[c for *_, c in marg], width=0.6)
    ax.set_ylabel("Δ ceiling NDO (marginal)")
    for i, (_, v, _) in enumerate(marg):
        ax.text(i, v + 0.0006, f"{v:+.3f}", ha="center", fontsize=9)
    ax.axhline(0, color=ink, lw=0.8)
    ax.set_title("Marginal value of each feature to the ceiling\nenergy_z's ~2x runtime buys almost nothing", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p2 = outdir / "fig2_feature_marginal_value.png"
    fig.savefig(p2, dpi=150)
    plt.close(fig)
    paths.append(p2)

    # Fig 3 — why 75% pairwise accuracy still regresses top-1.
    bm = results["autopsy"]["pairwise_accuracy_by_margin"]
    labels = list(bm.keys())
    accs = [bm[k]["accuracy"] for k in labels]
    ns = [bm[k]["n"] for k in labels]
    fig, ax = plt.subplots(figsize=(6.8, 3.8))
    bars = ax.bar(labels, accs, color=accent, width=0.68)
    ax.axhline(0.5, color=warn, ls="--", lw=1.2)
    ax.text(len(labels) - 0.5, 0.51, "chance", color=warn, fontsize=8, ha="right")
    ax.set_ylim(0.4, 0.85)
    ax.set_ylabel("pairwise ranking accuracy")
    ax.set_xlabel("NDO margin between the two candidates")
    for b, a, n in zip(bars, accs, ns):
        ax.text(b.get_x() + b.get_width() / 2, a + 0.006, f"{a:.2f}\nn={n//1000}k", ha="center", fontsize=8)
    ax.set_title("Reranker is near-chance on the close pairs that decide top-1;\n"
                 "high average accuracy comes only from easy, far-apart pairs", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    p3 = outdir / "fig3_pairwise_accuracy_by_margin.png"
    fig.savefig(p3, dpi=150)
    plt.close(fig)
    paths.append(p3)
    return paths


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scores", type=Path, default=DEFAULT_SCORES)
    ap.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    ap.add_argument("--weights", type=Path, default=DEFAULT_WEIGHTS)
    ap.add_argument("--json-out", type=Path, default=None)
    ap.add_argument("--figures-dir", type=Path, default=None)
    args = ap.parse_args()

    results = run_all(args.scores, args.features, args.weights)
    print_report(results)
    if args.json_out:
        args.json_out.write_text(json.dumps(results, indent=2, default=str))
        print(f"wrote {args.json_out}")
    if args.figures_dir:
        for path in make_figures(results, args.figures_dir):
            print(f"wrote {path}")


if __name__ == "__main__":
    main()
