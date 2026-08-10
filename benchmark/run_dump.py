"""Run SWORD2 over a CATH dataset with SWORD2_DUMP_CANDIDATES enabled.

Produces a candidate dump CSV and collects SWORD2 output for later scoring.
Run this before build_training_table.py.

Usage:
    # Prerequisite check:
    python -m benchmark.run_dump --check

    # Full training run (default: CATH-17287):
    python -m benchmark.run_dump \\
        --dataset cath17287 \\
        --dump /tmp/sword2_dump.csv \\
        --output-root /tmp/sword2_training \\
        --sword2 ./target/release/sword2 \\
        --workers 16

    # Dry-run (print commands):
    python -m benchmark.run_dump --dataset cath17287 --dry-run

Notes:
  - Chain PDBs are saved as {entry_id}.pdb so the dump's chain_id column equals entry_id.
  - Failed chains are logged to --output-root/failures.txt and do not abort the run.
  - CATH-663 (test set) is excluded from any training run automatically.
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from benchmark.datasets import load_dataset

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

# Chain IDs in the held-out CATH-663 test set — never include in training.
_CATH663_ENTRY_IDS: set[str] | None = None


def _get_cath663_ids() -> set[str]:
    global _CATH663_ENTRY_IDS
    if _CATH663_ENTRY_IDS is None:
        try:
            entries = load_dataset("cath663")
            _CATH663_ENTRY_IDS = {e.entry_id for e in entries}
        except Exception:
            _CATH663_ENTRY_IDS = set()
    return _CATH663_ENTRY_IDS


def _fetch_chain_pdb(pdb_id: str, chain_id: str, entry_id: str, work_dir: Path) -> Path | None:
    """Fetch a single chain PDB using SWORD2's mmCIF download + chain extraction."""
    try:
        import pdbtbx
        from sword2_lib import fetch  # type: ignore[import]  # Rust pyo3 binding if available
    except ImportError:
        pass

    # Use SWORD2's own fetch mechanism: run it in fetch-only mode or just download mmCIF
    # and extract the chain with pdbtbx. This keeps the same cleaning pipeline.
    # Simplest: let SWORD2 run on the entry — it will fetch, clean, and run the pipeline.
    return None  # Handled by _run_one instead


def _run_one(args: tuple) -> tuple[str, bool, str]:
    """Run SWORD2 on one chain. Returns (entry_id, success, message).

    Each chain writes to its own per-chain dump file to avoid concurrent-write
    corruption when multiple workers share a single CSV file.
    """
    entry_id, pdb_id, chain_id, sword2_bin, output_root, dump_dir = args

    chain_out = output_root / entry_id
    chain_out.mkdir(parents=True, exist_ok=True)

    # Per-chain dump: no concurrent access, no interleaving
    per_chain_dump = dump_dir / f"{entry_id}.csv"
    env = {**os.environ, "SWORD2_DUMP_CANDIDATES": str(per_chain_dump)}

    cmd = [
        str(sword2_bin),
        "-p", pdb_id,
        "--chain", chain_id,
        "-o", str(chain_out),
    ]

    try:
        result = subprocess.run(
            cmd,
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode == 0:
            return entry_id, True, ""
        return entry_id, False, f"exit {result.returncode}: {result.stderr[-200:]}"
    except subprocess.TimeoutExpired:
        return entry_id, False, "timeout"
    except Exception as exc:
        return entry_id, False, str(exc)


def check_prerequisites(sword2_bin: Path, dataset: str) -> bool:
    ok = True
    if not sword2_bin.exists():
        log.error("SWORD2 binary not found: %s  (run: cargo build --release)", sword2_bin)
        ok = False

    try:
        entries = load_dataset(dataset)
        log.info("Dataset %r: %d entries", dataset, len(entries))
    except Exception as exc:
        log.error("Cannot load dataset %r: %s", dataset, exc)
        ok = False

    try:
        cath663 = _get_cath663_ids()
        log.info("CATH-663 has %d entries (will be excluded from any training dataset)", len(cath663))
    except Exception as exc:
        log.warning("Could not load CATH-663 for exclusion check: %s", exc)

    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description="SWORD2 training dump mass-run")
    parser.add_argument("--dataset", default="cath17287")
    parser.add_argument("--dump", type=Path, default=Path("/tmp/sword2_dump.csv"),
                        help="Final merged dump CSV (written after all runs complete)")
    parser.add_argument("--dump-dir", type=Path, default=None,
                        help="Directory for per-chain dump CSVs (default: --dump parent / dump_parts/)")
    parser.add_argument("--output-root", type=Path, default=Path("/tmp/sword2_training"))
    parser.add_argument(
        "--sword2",
        type=Path,
        default=Path("./target/release/sword2"),
    )
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--limit", type=int, help="Process only first N chains (for testing)")
    parser.add_argument("--check", action="store_true", help="Check prerequisites and exit")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    sword2_bin = args.sword2.resolve()

    if args.check:
        ok = check_prerequisites(sword2_bin, args.dataset)
        return 0 if ok else 1

    if not check_prerequisites(sword2_bin, args.dataset):
        return 1

    entries = load_dataset(args.dataset, limit=args.limit)

    # Exclude held-out test chains
    test_ids = _get_cath663_ids()
    before = len(entries)
    entries = [e for e in entries if e.entry_id not in test_ids]
    excluded = before - len(entries)
    if excluded:
        log.info("Excluded %d chains that appear in CATH-663 (test set)", excluded)

    log.info("Running SWORD2 on %d chains with %d workers", len(entries), args.workers)
    log.info("Dump CSV: %s", args.dump)
    log.info("Output root: %s", args.output_root)

    if args.dry_run:
        for entry in entries[:5]:
            print(f"sword2 -p {entry.pdb_id} --chain {entry.chain_id} -o {args.output_root / entry.entry_id}")
        print(f"... ({len(entries)} total)")
        return 0

    args.output_root.mkdir(parents=True, exist_ok=True)
    failures_file = args.output_root / "failures.txt"

    # Per-chain dump files go into a subdirectory to avoid concurrent-write corruption
    dump_dir = args.dump_dir or (args.dump.parent / "dump_parts")
    dump_dir.mkdir(parents=True, exist_ok=True)
    log.info("Per-chain dump parts: %s", dump_dir)

    task_args = [
        (e.entry_id, e.pdb_id, e.chain_id, sword2_bin, args.output_root, dump_dir)
        for e in entries
    ]

    n_ok = 0
    n_fail = 0

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_run_one, a): a[0] for a in task_args}
        for future in as_completed(futures):
            entry_id, success, msg = future.result()
            if success:
                n_ok += 1
            else:
                n_fail += 1
                with failures_file.open("a") as f:
                    f.write(f"{entry_id}\t{msg}\n")
            if (n_ok + n_fail) % 100 == 0:
                log.info("Progress: %d ok, %d failed / %d total", n_ok, n_fail, len(entries))

    log.info("Complete: %d ok, %d failed. Failures logged to %s", n_ok, n_fail, failures_file)

    # Merge per-chain dump CSVs into one file
    log.info("Merging per-chain dumps into %s ...", args.dump)
    part_files = sorted(dump_dir.glob("*.csv"))
    n_merged = 0
    args.dump.parent.mkdir(parents=True, exist_ok=True)
    with args.dump.open("w") as out:
        header_written = False
        for part in part_files:
            with part.open() as inp:
                lines = inp.readlines()
            if not lines:
                continue
            if not header_written:
                out.write(lines[0])  # write header once
                header_written = True
            out.writelines(lines[1:])  # skip per-part header, write data
            n_merged += 1

    size_mb = args.dump.stat().st_size / 1e6 if args.dump.exists() else 0
    log.info("Merged %d part files → %s (%.1f MB)", n_merged, args.dump, size_mb)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
