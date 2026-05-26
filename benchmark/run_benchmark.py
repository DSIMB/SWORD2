#!/usr/bin/env python3
"""Benchmark harness for the mypmfs scoring backend across protein sizes.

Times the scoring engine(s) on a set of structures spanning different sizes and
reports wall-clock (median / min / stddev) plus the produced Pseudo-energy and
Z-score (as a produced-output sanity check, not an equivalence assertion).

Today only the C++ engine (`bin/mypmfs-master/scoring_omp`) exists. The Rust
scorer (stream 4 of the migration) does not exist yet; when it lands, add one
entry to ENGINES (see the commented `rust` stub) and the comparison table grows
a column automatically.

The scorer is *isolated* from the rest of the pipeline: we use `sword2` itself
(with --disable-energies) purely as the PDB-cleaning front-end to produce the
exact cleaned/renumbered `input.pdb` the scorer consumes, then time the scorer
on that file. This matches how the pipeline invokes scoring_omp:

    scoring_omp -i <cleaned.pdb> -d <potential_dir> [-q <residues>] -z -s <shuffles>

Config baked into 025_30_100_potential/parameters.log: REPRES=CA, DISTMAX=15,
DISTMIN=0, DIFFMIN=3, DIFFMAX=5100; default num_shuffles=2000.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

# Repo root = parent of this file's directory (benchmark/..).
REPO_ROOT = Path(__file__).resolve().parent.parent

DEFAULT_PDB_IDS = ["1CRN", "1UBQ", "1AKE", "1JX4", "1QCF"]

ENERGY_RE = re.compile(r"^Pseudo-energy\s*=\s*(.+)$", re.MULTILINE)
ZSCORE_RE = re.compile(r"^Z-score\s*=\s*(.+)$", re.MULTILINE)


@dataclass
class Engine:
    """A scoring engine: builds the argv to score a cleaned PDB file.

    `build` returns (argv, env_overrides). `zscore=False` requests raw energy
    only (deterministic distance/interpolation cost, no decoy shuffling).
    """

    name: str
    build: callable  # (pdb, potential_dir, shuffles, cpu, zscore) -> (list[str], dict)


def cpp_engine(scoring_bin: str) -> Engine:
    def build(pdb, potential_dir, shuffles, cpu, zscore):
        argv = [scoring_bin, "-i", pdb, "-d", potential_dir]
        if zscore:
            argv += ["-z", "-s", str(shuffles)]
        # scoring_omp parallelizes via OpenMP; threads come from OMP_NUM_THREADS,
        # not a CLI flag. cpu==0 means "let OpenMP decide" (omit override).
        env = {"OMP_NUM_THREADS": str(cpu)} if cpu > 0 else {}
        return argv, env

    return Engine("cpp", build)


# --- Rust engine (pure-Rust scorer via the `sword2 score` subcommand) ---
#
# Mirrors the C++ invocation: whole-structure single call, scoring the same
# cleaned PDB with the same potentials.
#
# Fairness note: time the SAME whole-structure single call shown here. Do NOT
# compare against pipeline-level timing — the production C++ path forks one
# process per PU/domain, so a pipeline comparison measures process-startup
# amortization, not the scoring algorithm.
def rust_engine(sword2_bin: str) -> Engine:
    def build(pdb, potential_dir, shuffles, cpu, zscore):
        argv = [sword2_bin, "score", "--pdb", pdb,
                "--potential-dir", potential_dir, "--cpu", str(cpu)]
        if zscore:
            argv += ["--shuffles", str(shuffles)]
        else:
            argv += ["--no-zscore"]
        return argv, {}

    return Engine("rust", build)


def run_once(argv: list[str], env_overrides: dict) -> tuple[float, str]:
    """Run argv once, return (wall_seconds, stdout). Raises on non-zero exit."""
    env = os.environ.copy()
    env.update(env_overrides)
    start = time.perf_counter()
    proc = subprocess.run(
        argv, capture_output=True, text=True, env=env, cwd=REPO_ROOT
    )
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        raise RuntimeError(
            f"command failed ({proc.returncode}): {' '.join(argv)}\n"
            f"stderr:\n{proc.stderr.strip()}"
        )
    return elapsed, proc.stdout


def parse_outputs(stdout: str) -> tuple[str, str]:
    e = ENERGY_RE.search(stdout)
    z = ZSCORE_RE.search(stdout)
    return (e.group(1).strip() if e else "", z.group(1).strip() if z else "")


def clean_structure(sword2_bin: str, pdb_id: str, workdir: Path) -> Path:
    """Fetch+clean a structure via sword2 --disable-energies; return cleaned PDB.

    Cached: if a cleaned input.pdb already exists under workdir/<id>/, reuse it.
    """
    out_dir = workdir / pdb_id
    existing = glob.glob(str(out_dir / "*" / "input.pdb"))
    if existing:
        return Path(existing[0])

    out_dir.mkdir(parents=True, exist_ok=True)
    argv = [sword2_bin, "-p", pdb_id, "-o", str(out_dir), "--disable-energies"]
    proc = subprocess.run(argv, capture_output=True, text=True, cwd=REPO_ROOT)
    if proc.returncode != 0:
        raise RuntimeError(
            f"sword2 cleaning failed for {pdb_id} ({proc.returncode}):\n"
            f"{proc.stderr.strip()}"
        )
    found = glob.glob(str(out_dir / "*" / "input.pdb"))
    if not found:
        raise RuntimeError(
            f"no cleaned input.pdb produced for {pdb_id} under {out_dir}"
        )
    return Path(found[0])


def count_residues(pdb_path: Path) -> int:
    """Count residues by counting CA ATOM records in the cleaned single chain."""
    n = 0
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                n += 1
    return n


@dataclass
class Row:
    pdb_id: str
    n_residues: int
    engine: str
    mode: str  # "raw" | "full"
    median_s: float
    min_s: float
    stddev_s: float
    energy: str
    zscore: str


def benchmark_one(
    engine: Engine,
    pdb_path: Path,
    potential_dir: str,
    shuffles: int,
    cpu: int,
    zscore: bool,
    repeats: int,
) -> tuple[list[float], str, str]:
    argv, env = engine.build(str(pdb_path), potential_dir, shuffles, cpu, zscore)
    # Warmup (also surfaces failures before the timed loop).
    _, stdout = run_once(argv, env)
    energy, zsc = parse_outputs(stdout)
    times = []
    for _ in range(repeats):
        t, out = run_once(argv, env)
        times.append(t)
        if zscore:  # capture the (varying) z-score from a timed run too
            _, zsc = parse_outputs(out)
    return times, energy, zsc


def fmt_table(rows: list[Row]) -> str:
    header = [
        "pdb_id", "n_residues", "engine", "mode",
        "median_s", "min_s", "stddev_s", "energy", "zscore",
    ]
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        lines.append("| " + " | ".join([
            r.pdb_id, str(r.n_residues), r.engine, r.mode,
            f"{r.median_s:.4f}", f"{r.min_s:.4f}", f"{r.stddev_s:.4f}",
            r.energy, r.zscore,
        ]) + " |")
    return "\n".join(lines)


def write_csv(rows: list[Row], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "pdb_id", "n_residues", "engine", "mode",
            "median_s", "min_s", "stddev_s", "energy", "zscore",
        ])
        for r in rows:
            w.writerow([
                r.pdb_id, r.n_residues, r.engine, r.mode,
                f"{r.median_s:.6f}", f"{r.min_s:.6f}", f"{r.stddev_s:.6f}",
                r.energy, r.zscore,
            ])


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pdb-ids", nargs="+", default=DEFAULT_PDB_IDS,
                   help=f"PDB ids to benchmark (default: {' '.join(DEFAULT_PDB_IDS)})")
    p.add_argument("--repeats", type=int, default=5,
                   help="timed runs per measurement (default: 5)")
    p.add_argument("--shuffles", type=int, default=2000,
                   help="Z-score decoys (default: 2000, matches pipeline)")
    p.add_argument("--cpu", type=int, default=0,
                   help="OMP threads; 0 = let OpenMP decide (default: 0)")
    p.add_argument("--engines", nargs="+", default=["cpp", "rust"],
                   help="engines to run (default: cpp rust)")
    p.add_argument("--sword2-bin",
                   default=str(REPO_ROOT / "target" / "release" / "sword2"))
    p.add_argument("--scoring-bin",
                   default=str(REPO_ROOT / "bin" / "mypmfs-master" / "scoring_omp"))
    p.add_argument("--potential-dir",
                   default=str(REPO_ROOT / "bin" / "mypmfs-master"
                               / "025_30_100_potential"))
    p.add_argument("--workdir", default=str(REPO_ROOT / "benchmark" / "structures"))
    p.add_argument("--out-dir", default=str(REPO_ROOT / "benchmark"))
    args = p.parse_args()

    for label, path in [("sword2", args.sword2_bin), ("scoring", args.scoring_bin)]:
        if "cpp" in args.engines and label == "scoring" and not Path(path).exists():
            print(f"error: {label} binary not found: {path}\n"
                  f"build it with: make -C bin/mypmfs-master", file=sys.stderr)
            return 1
        if label == "sword2" and not Path(path).exists():
            print(f"error: {label} binary not found: {path}\n"
                  f"build it with: cargo build --release", file=sys.stderr)
            return 1
    if not Path(args.potential_dir).is_dir():
        print(f"error: potential dir not found: {args.potential_dir}", file=sys.stderr)
        return 1

    # Build the requested engines.
    engines: list[Engine] = []
    for name in args.engines:
        if name == "cpp":
            engines.append(cpp_engine(args.scoring_bin))
        elif name == "rust":
            engines.append(rust_engine(args.sword2_bin))
        else:
            print(f"error: unknown engine '{name}'. Available: cpp, rust.",
                  file=sys.stderr)
            return 1

    workdir = Path(args.workdir)
    rows: list[Row] = []

    for pdb_id in args.pdb_ids:
        print(f"\n=== {pdb_id} ===", file=sys.stderr)
        try:
            cleaned = clean_structure(args.sword2_bin, pdb_id, workdir)
        except RuntimeError as e:
            print(f"  skipped ({e})", file=sys.stderr)
            continue
        n_res = count_residues(cleaned)
        print(f"  cleaned -> {cleaned}  ({n_res} residues)", file=sys.stderr)

        for engine in engines:
            for mode, zscore in [("raw", False), ("full", True)]:
                print(f"  {engine.name}/{mode}: timing "
                      f"{args.repeats} runs ...", file=sys.stderr)
                try:
                    times, energy, zsc = benchmark_one(
                        engine, cleaned, args.potential_dir,
                        args.shuffles, args.cpu, zscore, args.repeats,
                    )
                except RuntimeError as e:
                    print(f"    FAILED: {e}", file=sys.stderr)
                    continue
                rows.append(Row(
                    pdb_id=pdb_id, n_residues=n_res,
                    engine=engine.name, mode=mode,
                    median_s=statistics.median(times),
                    min_s=min(times),
                    stddev_s=statistics.stdev(times) if len(times) > 1 else 0.0,
                    energy=energy, zscore=(zsc if zscore else ""),
                ))

    if not rows:
        print("\nNo successful measurements.", file=sys.stderr)
        return 1

    table = fmt_table(rows)
    print("\n" + table + "\n")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(rows, out_dir / "results.csv")
    with open(out_dir / "results.md", "w") as fh:
        fh.write("# scoring_omp benchmark results\n\n")
        fh.write(f"- repeats: {args.repeats}, shuffles: {args.shuffles}, "
                 f"cpu(OMP): {args.cpu or 'auto'}\n")
        fh.write(f"- engines: {', '.join(args.engines)}\n\n")
        fh.write(table + "\n")
    print(f"wrote {out_dir/'results.csv'} and {out_dir/'results.md'}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
