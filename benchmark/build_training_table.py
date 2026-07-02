"""Build a training table from a SWORD2 candidate dump and a CATH reference CSV.

Usage (fast, recommended):
    python -m benchmark.build_training_table \\
        --dump-dir /tmp/dump_parts/ \\
        --reference cath17287 \\
        --output benchmark/data/training_table.csv \\
        --workers 16

Usage (single merged CSV):
    python -m benchmark.build_training_table \\
        --dump /path/to/dump.csv \\
        --reference cath17287 \\
        --output benchmark/data/training_table.csv

Coordinate mapping:
  The dump delineation uses 0-based CA-array positions (same coordinate system as
  map_author_chopping output). No offset conversion is needed: delineation_to_chopping()
  reformats domain/segment separators and the result is directly comparable to
  the 0-based sequential chopping from map_author_chopping().
  The PDB for coordinate mapping is read from output_dir/intermediate/<chain_id>.pdb.
"""
from __future__ import annotations

import argparse
import csv
import io
import logging
import math
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any

from benchmark.datasets import load_dataset, read_merizo_csv, strip_cath_labels
from benchmark.metrics import score_choppings
from benchmark.numbering import map_author_chopping, numbering_from_pdb

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

# Set by _init_worker in each subprocess
_worker_ref: dict[str, tuple[str, str, int]] | None = None

FIELDNAMES = [
    "chain_id",
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "delineation",
    "boundary_coil_fraction",
    "energy_z",
    "modal_count_distance",
    "n_true_domains",
    "n_pred_domains",
    "ndo",
    "iou",
    "boundary_f1_10",
    "matched_dice",
    "d_count_acc",
    "S",
    "is_oracle_s",
]


def delineation_to_chopping(delineation: str) -> str:
    """Convert MeasureLine delineation to score_choppings chopping format.

    MeasureLine: space=domain separator, semicolon=segment separator within domain.
    Chopping:    comma=domain separator, underscore=segment separator within domain.
    """
    domains = delineation.strip().split()
    return ",".join(
        "_".join(s.strip() for s in d.split(";") if s.strip())
        for d in domains
    )


def _load_reference(reference_arg: str) -> dict[str, tuple[str, str, int]]:
    """Load reference choppings, indexed by both CATH entry_id and dump chain_id formats."""
    if Path(reference_arg).exists():
        entries = read_merizo_csv(Path(reference_arg), dataset="custom")
    else:
        entries = load_dataset(reference_arg)
    ref: dict[str, tuple[str, str, int]] = {}
    for e in entries:
        val = (e.chain_id, e.chopping, e.n_residues)
        ref[e.entry_id] = val                           # "12e8H"
        ref[f"{e.pdb_id.upper()}_{e.chain_id}"] = val  # "12E8_H"
    return ref


def _true_chopping_for_chain(
    raw_chopping: str,
    chain_id: str,
    pdb_path: Path | None,
) -> tuple[str, int]:
    """Convert CATH author chopping to 0-based sequential chopping. Returns (chopping, n_res)."""
    stripped = strip_cath_labels(raw_chopping)
    if pdb_path is not None and pdb_path.exists():
        try:
            numbering = numbering_from_pdb(pdb_path, chain_id=chain_id)
            true_chop = map_author_chopping(stripped, numbering, chain_id=chain_id)
            return true_chop, numbering.n_residues
        except Exception as exc:
            log.debug("PDB numbering failed for %s: %s", pdb_path, exc)
    return stripped, 0


def _find_pdb(output_dir: Path, chain_id: str) -> Path | None:
    for name in [chain_id, chain_id.lower()]:
        p = output_dir / "intermediate" / f"{name}.pdb"
        if p.exists():
            return p
    return None


def _score_candidates(
    chain_id: str,
    candidates: list[dict[str, str]],
    reference: dict[str, tuple[str, str, int]],
) -> list[dict[str, Any]]:
    """Score all candidates for one chain. Returns list of output row dicts."""
    if chain_id not in reference:
        return []

    ref_chain_id, raw_chopping, ref_n_res = reference[chain_id]

    output_dir = Path(candidates[0].get("output_dir", "")) if candidates else Path()
    pdb_path = _find_pdb(output_dir, chain_id)
    true_chop, n_res = _true_chopping_for_chain(raw_chopping, ref_chain_id, pdb_path)
    if n_res == 0:
        n_res = ref_n_res

    scored: list[tuple[dict[str, str], float, Any]] = []

    for row in candidates:
        raw_del = row.get("delineation")
        if not raw_del:
            continue
        delineation = raw_del.strip().strip('"')
        if not delineation:
            continue

        pred_chop = delineation_to_chopping(delineation)
        try:
            with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
                metrics = score_choppings(true_chop, pred_chop, n_res=n_res if n_res > 0 else None)
        except Exception:
            continue

        count_term = 1.0 / (1.0 + abs(metrics.n_pred_domains - metrics.n_true_domains))
        vals = [metrics.ndo, metrics.iou, metrics.boundary_f1_10, metrics.matched_dice, count_term]
        finite = [v for v in vals if not math.isnan(v)]
        s_val = statistics.mean(finite) if finite else 0.0
        scored.append((row, s_val, metrics))

    if not scored:
        return []

    best_s = max(s for _, s, _ in scored)

    rows: list[dict[str, Any]] = []
    for row, s_val, metrics in scored:
        rows.append({
            "chain_id": chain_id,
            "num_domains": row["num_domains"],
            "min_size": row["min_size"],
            "max_cr": row["max_cr"],
            "density_min": row["density_min"],
            "mean_density": row["mean_density"],
            "delineation": row.get("delineation", "").strip().strip('"'),
            "boundary_coil_fraction": row.get("boundary_coil_fraction", ""),
            "energy_z": row.get("energy_z", ""),
            "modal_count_distance": row.get("modal_count_distance", ""),
            "n_true_domains": metrics.n_true_domains,
            "n_pred_domains": metrics.n_pred_domains,
            "ndo": metrics.ndo,
            "iou": metrics.iou,
            "boundary_f1_10": metrics.boundary_f1_10,
            "matched_dice": metrics.matched_dice,
            "d_count_acc": metrics.d_count_acc,
            "S": s_val,
            "is_oracle_s": 1 if abs(s_val - best_s) < 1e-9 else 0,
        })
    return rows


# ---------------------------------------------------------------------------
# Multiprocessing worker (per-chain file)
# ---------------------------------------------------------------------------

def _init_worker(ref: dict[str, tuple[str, str, int]]) -> None:
    global _worker_ref
    _worker_ref = ref


def _score_chain_file(part_file: Path) -> list[dict[str, Any]] | None:
    """Worker entry point: process one per-chain dump CSV file."""
    global _worker_ref
    chain_id = part_file.stem  # entry_id = filename without .csv

    try:
        with part_file.open(newline="") as f:
            candidates = list(csv.DictReader(f))
    except Exception:
        return None

    if not candidates:
        return None

    return _score_candidates(chain_id, candidates, _worker_ref)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Main processing functions
# ---------------------------------------------------------------------------

def build_from_dump_dir(
    dump_dir: Path,
    reference: dict[str, tuple[str, str, int]],
    output_path: Path,
    workers: int = 8,
) -> None:
    """Stream per-chain dump files from a directory — fast, low-memory."""
    part_files = sorted(dump_dir.glob("*.csv"))
    log.info("Processing %d per-chain dump files with %d workers", len(part_files), workers)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_chains = n_rows = n_skipped = 0

    with output_path.open("w", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=FIELDNAMES)
        writer.writeheader()

        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_init_worker,
            initargs=(reference,),
        ) as pool:
            for chain_rows in pool.map(_score_chain_file, part_files, chunksize=20):
                if not chain_rows:
                    n_skipped += 1
                else:
                    n_chains += 1
                    n_rows += len(chain_rows)
                    for row in chain_rows:
                        writer.writerow(row)

                done = n_chains + n_skipped
                if done % 500 == 0:
                    log.info(
                        "Progress: %d / %d  (%d scored, %d skipped, %d rows)",
                        done, len(part_files), n_chains, n_skipped, n_rows,
                    )

    log.info("Done: %d chains scored, %d skipped, %d rows written", n_chains, n_skipped, n_rows)


def build_from_dump_csv(
    dump_path: Path,
    reference: dict[str, tuple[str, str, int]],
    output_path: Path,
) -> None:
    """Process a merged dump CSV — streaming groupby (lower RAM than list+sort)."""
    from itertools import groupby

    log.info("Reading dump CSV: %s", dump_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_chains = n_rows = n_skipped = 0

    # Collect rows grouped by chain_id using a dict (single linear pass, no sort)
    chain_buckets: dict[str, list[dict[str, str]]] = {}
    with dump_path.open(newline="") as src:
        for row in csv.DictReader(src):
            cid = row.get("chain_id", "")
            if cid:
                chain_buckets.setdefault(cid, []).append(row)

    log.info("Loaded %d chains from dump", len(chain_buckets))

    with output_path.open("w", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=FIELDNAMES)
        writer.writeheader()

        for chain_id, candidates in chain_buckets.items():
            chain_rows = _score_candidates(chain_id, candidates, reference)
            if chain_rows:
                n_chains += 1
                n_rows += len(chain_rows)
                for row in chain_rows:
                    writer.writerow(row)
            else:
                n_skipped += 1

            done = n_chains + n_skipped
            if done % 500 == 0:
                log.info(
                    "Progress: %d / %d  (%d scored, %d skipped, %d rows)",
                    done, len(chain_buckets), n_chains, n_skipped, n_rows,
                )

    log.info("Done: %d chains scored, %d skipped, %d rows written", n_chains, n_skipped, n_rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Build SWORD2 training table")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--dump", type=Path, help="Merged dump CSV (all chains)")
    src.add_argument("--dump-dir", type=Path, help="Directory of per-chain dump CSVs (faster)")
    parser.add_argument("--reference", default="cath17287",
                        help="Named dataset or path to Merizo CSV (default: cath17287)")
    parser.add_argument("--output", type=Path, default=Path("benchmark/data/training_table.csv"))
    parser.add_argument("--workers", type=int, default=8,
                        help="Parallel workers when using --dump-dir (default: 8)")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    log.info("Loading reference: %s", args.reference)
    reference = _load_reference(args.reference)
    log.info("Reference: %d entries (dual-keyed)", len(reference))

    if args.dump_dir:
        build_from_dump_dir(args.dump_dir, reference, args.output, workers=args.workers)
    else:
        build_from_dump_csv(args.dump, reference, args.output)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
