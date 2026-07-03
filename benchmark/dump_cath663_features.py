#!/usr/bin/env python3
"""Generate a candidate-feature dump over a CATH dataset for selection diagnostics.

Runs the release ``sword2`` binary on each cached single-chain structure with
``SWORD2_DUMP_CANDIDATES`` enabled. That env var makes SWORD2 write, for every
post-shortlist candidate (the ``relevant_measure2`` set the reranker would see),
one row of features:

    num_domains, min_size, max_cr, density_min, mean_density,
    delineation (0-based, space-separated domains / ';'-separated segments),
    boundary_coil_fraction, energy_z, modal_count_distance

energy_z is always computed (the pipeline builds a 200-shuffle reranker energy
config unconditionally), so ``-E`` is not needed.

Each row is tagged with ``entry_id`` (taken from the cached chain filename, which
matches ``scores.csv``) and all rows are concatenated into one CSV that
``diagnose_selection.py`` joins against ``results_bypass/scores.csv``.
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark.datasets import load_dataset

REPO = Path(__file__).resolve().parents[1]
BINARY = REPO / "target/release/sword2"
CHAINS = REPO / "benchmark/cache/chains"

DUMP_COLUMNS = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "delineation",
    "boundary_coil_fraction",
    "energy_z",
    "modal_count_distance",
]


def run_one(entry_id: str, tmpdir: str, threads: int) -> tuple[str, list[dict] | None, str | None]:
    pdb = CHAINS / f"{entry_id}.pdb"
    if not pdb.exists():
        return entry_id, None, "missing structure"
    dump = Path(tmpdir) / f"{entry_id}.csv"
    out = Path(tmpdir) / f"out_{entry_id}"
    env = dict(os.environ, SWORD2_DUMP_CANDIDATES=str(dump), RUST_LOG="error")
    try:
        subprocess.run(
            [str(BINARY), "-i", str(pdb), "-o", str(out), "--threads", str(threads)],
            cwd=REPO,
            env=env,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=300,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        return entry_id, None, f"run failed: {exc}"
    if not dump.exists():
        return entry_id, None, "no dump produced"
    return entry_id, list(csv.DictReader(dump.open())), None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", default="cath663")
    ap.add_argument("--out", type=Path, default=REPO / "benchmark/data/cath663_candidate_features.csv")
    ap.add_argument("--jobs", type=int, default=8, help="parallel sword2 invocations")
    ap.add_argument("--threads", type=int, default=2, help="rayon threads per invocation")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    if not BINARY.exists():
        sys.exit(f"binary not found: {BINARY} (run `cargo build --release`)")

    entries = load_dataset(args.dataset, limit=args.limit)
    out_rows: list[dict] = []
    failures: list[tuple[str, str]] = []
    with tempfile.TemporaryDirectory() as tmp:
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_one, e.entry_id, tmp, args.threads): e.entry_id for e in entries}
            done = 0
            for fut in as_completed(futs):
                entry_id, rows, err = fut.result()
                done += 1
                if err:
                    failures.append((entry_id, err))
                else:
                    for row in rows:
                        row["entry_id"] = entry_id
                        out_rows.append(row)
                if done % 100 == 0:
                    print(f"{done}/{len(entries)} done, {len(failures)} failures", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["entry_id", *DUMP_COLUMNS], extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_rows)

    print(
        f"wrote {len(out_rows)} candidate rows for {len(entries) - len(failures)}/{len(entries)} "
        f"entries to {args.out}"
    )
    if failures:
        print(f"{len(failures)} failures (first 10): {failures[:10]}")


if __name__ == "__main__":
    main()
