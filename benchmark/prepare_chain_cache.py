"""Populate benchmark/cache/{raw,chains} for a dataset (parallel, resumable).

Two phases:
  1. Download each *unique* pdb_id's raw PDB from RCSB (parallel). Entries
     sharing a pdb_id (multi-chain complexes) are deduplicated so each raw
     file is fetched exactly once, avoiding redundant requests.
  2. Extract each entry's single chain from its already-downloaded raw file
     (parallel, no network).

Safe to re-run: existing raw/chain files are reused, only missing ones are
(re-)fetched/extracted.

Usage:
    python -m benchmark.prepare_chain_cache --dataset cath17287 --workers 12
    python -m benchmark.prepare_chain_cache --dataset cath17287 --limit 20 --dry-run
"""
from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import URLError, HTTPError

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark.datasets import CathEntry, load_dataset
from benchmark.structures import download_pdb, extract_single_chain

DEFAULT_CACHE_DIR = Path("benchmark/cache")
DEFAULT_FAILURES_LOG = Path("benchmark/cache/prepare_failures.txt")
MAX_ATTEMPTS = 3
RETRY_BACKOFF_SECONDS = 2.0


def _download_one(pdb_id: str, raw_dir: Path) -> tuple[str, bool, str]:
    last_error = ""
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            download_pdb(pdb_id, raw_dir)
            return pdb_id, True, ""
        except (URLError, HTTPError, TimeoutError, OSError) as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            if attempt < MAX_ATTEMPTS:
                time.sleep(RETRY_BACKOFF_SECONDS * attempt)
    return pdb_id, False, last_error


def _extract_one(entry: CathEntry, cache_dir: Path) -> tuple[str, bool, str]:
    chain_pdb = cache_dir / "chains" / f"{entry.entry_id}.pdb"
    if chain_pdb.exists():
        return entry.entry_id, True, "cached"

    raw_pdb = cache_dir / "raw" / f"{entry.pdb_id.lower()}.pdb"
    if not raw_pdb.exists():
        return entry.entry_id, False, "raw pdb missing (download failed)"

    try:
        extract_single_chain(raw_pdb, entry.chain_id, chain_pdb)
        return entry.entry_id, True, "extracted"
    except (ValueError, OSError) as exc:
        return entry.entry_id, False, f"{type(exc).__name__}: {exc}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="cath17287", help="Dataset key, e.g. cath17287")
    parser.add_argument("--limit", type=int, default=None, help="Limit entries for smoke runs")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--failures-log", type=Path, default=DEFAULT_FAILURES_LOG)
    parser.add_argument("--workers", type=int, default=8, help="Parallel workers")
    parser.add_argument("--dry-run", action="store_true", help="List entries without fetching")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    entries = load_dataset(args.dataset, limit=args.limit)
    print(f"{len(entries)} entries in dataset {args.dataset!r}")

    if args.dry_run:
        for entry in entries[:20]:
            print(f"  {entry.entry_id} ({entry.pdb_id} chain {entry.chain_id})")
        return 0

    cache_dir = args.cache_dir
    raw_dir = cache_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    (cache_dir / "chains").mkdir(parents=True, exist_ok=True)
    args.failures_log.parent.mkdir(parents=True, exist_ok=True)

    unique_pdb_ids = sorted({entry.pdb_id for entry in entries})
    print(f"{len(unique_pdb_ids)} unique pdb_ids to download")

    download_failures: dict[str, str] = {}
    n_done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_download_one, pdb_id, raw_dir): pdb_id for pdb_id in unique_pdb_ids}
        for future in as_completed(futures):
            pdb_id, ok, message = future.result()
            n_done += 1
            if not ok:
                download_failures[pdb_id] = message
            if n_done % 500 == 0 or n_done == len(unique_pdb_ids):
                print(
                    f"[download {n_done}/{len(unique_pdb_ids)}] failed={len(download_failures)}",
                    flush=True,
                )

    n_ok = 0
    n_extracted = 0
    n_failed = 0
    failures: list[tuple[str, str]] = []

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_extract_one, entry, cache_dir): entry for entry in entries}
        for i, future in enumerate(as_completed(futures), start=1):
            entry_id, ok, message = future.result()
            if ok:
                n_ok += 1
                if message == "extracted":
                    n_extracted += 1
            else:
                n_failed += 1
                failures.append((entry_id, message))

            if i % 500 == 0 or i == len(entries):
                print(f"[extract {i}/{len(entries)}] ok={n_ok} failed={n_failed}", flush=True)

    if failures:
        with args.failures_log.open("w") as handle:
            for entry_id, message in failures:
                handle.write(f"{entry_id}\t{message}\n")
        print(f"{len(failures)} failures logged to {args.failures_log}")

    print(f"Done: {n_ok}/{len(entries)} cached ({n_extracted} newly extracted), {n_failed} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
