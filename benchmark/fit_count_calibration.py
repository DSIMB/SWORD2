"""Fit the analytical domain-count calibration: expected_ndom = a + b * n_residues.

Mirrors benchmark/fit_models.py: fits constants offline and prints a Rust snippet to
paste into sword2-lib/src/sword/count_calibration.rs. Fit on CATH-17287; the constants
correct the distance_model selector's systematic under-segmentation (see
benchmark/DIAGNOSIS.md).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from benchmark.datasets import load_dataset

REPO = Path(__file__).resolve().parents[1]


def fit_linear_count(lengths: np.ndarray, counts: np.ndarray) -> tuple[float, float]:
    design = np.vstack([np.ones_like(lengths, dtype=float), lengths.astype(float)]).T
    coef, *_ = np.linalg.lstsq(design, counts.astype(float), rcond=None)
    return float(coef[0]), float(coef[1])


def round_accuracy(lengths, counts, intercept, len_coef) -> tuple[float, float]:
    pred = np.clip(np.round(intercept + len_coef * np.asarray(lengths, float)), 1, None)
    counts = np.asarray(counts, float)
    return float(np.mean(pred == counts)), float(np.mean(pred - counts))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train", default="cath17287")
    ap.add_argument("--eval", default="cath663")
    ap.add_argument("--lambda-default", type=float, default=0.05)
    ap.add_argument("--out", type=Path, default=REPO / "benchmark/data/count_calibration.json")
    args = ap.parse_args()

    train = load_dataset(args.train)
    eval_entries = load_dataset(args.eval)
    eval_ids = {e.entry_id for e in eval_entries}
    # keep the fit clean of the held-out benchmark
    train = [e for e in train if e.entry_id not in eval_ids]

    lengths = np.array([e.n_residues for e in train], float)
    counts = np.array([e.n_domains for e in train], float)
    intercept, len_coef = fit_linear_count(lengths, counts)
    acc, bias = round_accuracy(lengths, counts, intercept, len_coef)
    e_len = np.array([e.n_residues for e in eval_entries], float)
    e_cnt = np.array([e.n_domains for e in eval_entries], float)
    eacc, ebias = round_accuracy(e_len, e_cnt, intercept, len_coef)

    print(f"train {args.train}: n={len(train)}  intercept={intercept:.6f}  len_coef={len_coef:.6f}")
    print(f"  train round d_count_acc={acc:.3f} bias={bias:+.3f}")
    print(f"  eval  {args.eval} round d_count_acc={eacc:.3f} bias={ebias:+.3f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"intercept": intercept, "len_coef": len_coef, "lambda": args.lambda_default}, indent=2))
    print(f"\nwrote {args.out}\n\nRust snippet for count_calibration.rs Default:")
    print(f"    Self {{ intercept: {intercept:.6}, len_coef: {len_coef:.6}, "
          f"lambda: {args.lambda_default} }}")


if __name__ == "__main__":
    main()
