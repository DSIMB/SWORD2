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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from benchmark.candidate_geometry import CANDIDATE_GEOMETRY_FIELDS
from benchmark.datasets import CathEntry, load_dataset, read_merizo_csv
from benchmark.factorized_ranker.integrity import (
    IntegrityError,
    RejectionCode,
    RejectionRecord,
    map_cath_reference,
    validate_partition,
)
from benchmark.metrics import score_choppings

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)
REPO = Path(__file__).resolve().parents[1]

# Set by _init_worker in each subprocess
_worker_ref: dict[str, CathEntry] | None = None
_worker_chain_cache_dir: Path | None = None

BASE_CANDIDATE_FIELDS = [
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

FIELDNAMES = [
    "chain_id",
    *BASE_CANDIDATE_FIELDS,
    *CANDIDATE_GEOMETRY_FIELDS,
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

REJECTION_FIELDNAMES = ["chain_id", "scope", "code", "detail", "delineation"]

REQUIRED_NUMERIC_CANDIDATE_FIELDS = [
    "num_domains",
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "boundary_coil_fraction",
    "modal_count_distance",
    *CANDIDATE_GEOMETRY_FIELDS,
]
REQUIRED_CANDIDATE_COLUMNS = {"delineation", *REQUIRED_NUMERIC_CANDIDATE_FIELDS}


@dataclass(frozen=True)
class ScoredChain:
    rows: list[dict[str, Any]]
    rejections: list[RejectionRecord]


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


def _load_reference(reference_arg: str) -> dict[str, CathEntry]:
    """Load reference choppings, indexed by both CATH entry_id and dump chain_id formats."""
    if Path(reference_arg).exists():
        entries = read_merizo_csv(Path(reference_arg), dataset="custom")
    else:
        entries = load_dataset(reference_arg)
    ref: dict[str, CathEntry] = {}

    def add_alias(alias: str, entry: CathEntry) -> None:
        previous = ref.get(alias)
        if previous is not None and previous != entry:
            raise ValueError(
                f"reference alias collision for {alias!r}: "
                f"{previous.entry_id!r} and {entry.entry_id!r}"
            )
        ref[alias] = entry

    for e in entries:
        add_alias(e.entry_id, e)                           # "12e8H"
        add_alias(f"{e.pdb_id.upper()}_{e.chain_id}", e)  # "12E8_H"
    return ref


def _find_pdb(
    output_dir: Path | None,
    chain_id: str,
    chain_cache_dir: Path | None = None,
    canonical_entry_id: str | None = None,
) -> Path | None:
    """Find the canonical structure used to map CATH author numbers.

    The benchmark cache is preferred because it persists after temporary SWORD
    output directories are removed and is also what the benchmark scorer uses.
    """
    names = list(dict.fromkeys([canonical_entry_id, chain_id, chain_id.lower()]))
    if chain_cache_dir is not None:
        for name in names:
            if name is None:
                continue
            cached = chain_cache_dir / f"{name}.pdb"
            if cached.exists():
                return cached
    if output_dir is None:
        return None
    for name in names:
        if name is None:
            continue
        p = output_dir / "intermediate" / f"{name}.pdb"
        if p.exists():
            return p
    return None


def _score_candidates(
    chain_id: str,
    candidates: list[dict[str, str]],
    reference: dict[str, CathEntry],
    chain_cache_dir: Path | None = None,
) -> ScoredChain:
    """Score valid candidates for one chain and explain every exclusion."""
    if chain_id not in reference:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=chain_id,
                    scope="chain",
                    code=RejectionCode.MISSING_REFERENCE,
                    detail="chain is absent from the reference dataset",
                )
            ],
        )

    entry = reference[chain_id]
    canonical_chain_id = entry.entry_id

    if not candidates:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="chain",
                    code=RejectionCode.SCHEMA_MISMATCH,
                    detail="candidate dump contains no rows",
                )
            ],
        )

    raw_output_dir = candidates[0].get("output_dir", "").strip()
    output_dir = Path(raw_output_dir) if raw_output_dir else None
    pdb_path = _find_pdb(
        output_dir,
        chain_id,
        chain_cache_dir,
        canonical_entry_id=entry.entry_id,
    )
    if pdb_path is None:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="chain",
                    code=RejectionCode.MISSING_CHAIN_PDB,
                    detail="canonical chain PDB is unavailable",
                )
            ],
        )

    try:
        canonical_reference = map_cath_reference(entry, pdb_path)
    except IntegrityError as exc:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="chain",
                    code=exc.code,
                    detail=exc.detail,
                )
            ],
        )
    except Exception as exc:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="chain",
                    code=RejectionCode.AUTHOR_MAPPING_FAILED,
                    detail=f"unexpected mapping failure ({type(exc).__name__})",
                )
            ],
        )

    scored: list[tuple[dict[str, str], str, str, float, Any]] = []
    rejections: list[RejectionRecord] = []

    for row in candidates:
        raw_del = row.get("delineation")
        delineation = raw_del.strip().strip('"') if raw_del else ""

        invalid_field: str | None = None
        for field in REQUIRED_NUMERIC_CANDIDATE_FIELDS:
            raw_value = row.get(field, "")
            try:
                value = float(raw_value)
            except (TypeError, ValueError):
                invalid_field = field
                break
            if not math.isfinite(value):
                invalid_field = field
                break
        raw_energy = row.get("energy_z", "")
        normalized_energy = raw_energy or ""
        if raw_energy not in (None, ""):
            try:
                energy = float(raw_energy)
            except (TypeError, ValueError):
                normalized_energy = ""
            else:
                if not math.isfinite(energy):
                    normalized_energy = ""
        if invalid_field is not None:
            rejections.append(
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="candidate",
                    code=RejectionCode.NONFINITE_FEATURE,
                    detail=f"feature {invalid_field!r} is missing or non-finite",
                    delineation=delineation or None,
                )
            )
            continue

        declared_value = float(row["num_domains"])
        if not declared_value.is_integer() or declared_value <= 0:
            rejections.append(
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="candidate",
                    code=RejectionCode.CANDIDATE_DOMAIN_COUNT_MISMATCH,
                    detail="declared domain count is not a positive integer",
                    delineation=delineation or None,
                )
            )
            continue
        declared_domains = int(declared_value)

        try:
            partition = validate_partition(
                delineation,
                n_residues=canonical_reference.n_residues,
                declared_domains=declared_domains,
            )
        except IntegrityError as exc:
            rejections.append(
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="candidate",
                    code=exc.code,
                    detail=exc.detail,
                    delineation=delineation or None,
                )
            )
            continue

        pred_chop = delineation_to_chopping(partition.canonical_delineation)
        try:
            with redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
                metrics = score_choppings(
                    canonical_reference.chopping,
                    pred_chop,
                    n_res=canonical_reference.n_residues,
                )
        except Exception as exc:
            rejections.append(
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="candidate",
                    code=RejectionCode.SCORING_FAILED,
                    detail=f"partition scoring failed ({type(exc).__name__})",
                    delineation=partition.canonical_delineation,
                )
            )
            continue

        count_term = 1.0 / (1.0 + abs(metrics.n_pred_domains - metrics.n_true_domains))
        vals = [
            metrics.ndo,
            metrics.iou,
            metrics.boundary_f1_10,
            metrics.matched_dice,
            metrics.d_count_acc,
            count_term,
        ]
        if not all(math.isfinite(value) for value in vals):
            rejections.append(
                RejectionRecord(
                    chain_id=canonical_chain_id,
                    scope="candidate",
                    code=RejectionCode.SCORING_FAILED,
                    detail="partition scoring produced a non-finite metric",
                    delineation=partition.canonical_delineation,
                )
            )
            continue
        s_val = statistics.mean(vals[:4] + [count_term])
        scored.append(
            (
                row,
                partition.canonical_delineation,
                normalized_energy,
                s_val,
                metrics,
            )
        )

    if not scored:
        return ScoredChain(rows=[], rejections=rejections)

    best_s = max(s for _, _, _, s, _ in scored)

    rows: list[dict[str, Any]] = []
    for row, canonical_delineation, normalized_energy, s_val, metrics in scored:
        candidate_fields = {
            field: row.get(field, "") for field in BASE_CANDIDATE_FIELDS
        }
        candidate_fields["delineation"] = canonical_delineation
        candidate_fields["energy_z"] = normalized_energy
        rows.append(
            {
                "chain_id": canonical_chain_id,
                **candidate_fields,
                **{
                    field: row.get(field, "")
                    for field in CANDIDATE_GEOMETRY_FIELDS
                },
                "n_true_domains": metrics.n_true_domains,
                "n_pred_domains": metrics.n_pred_domains,
                "ndo": metrics.ndo,
                "iou": metrics.iou,
                "boundary_f1_10": metrics.boundary_f1_10,
                "matched_dice": metrics.matched_dice,
                "d_count_acc": metrics.d_count_acc,
                "S": s_val,
                "is_oracle_s": 1 if abs(s_val - best_s) < 1e-9 else 0,
            }
        )
    return ScoredChain(rows=rows, rejections=rejections)


# ---------------------------------------------------------------------------
# Multiprocessing worker (per-chain file)
# ---------------------------------------------------------------------------

def _init_worker(ref: dict[str, CathEntry], chain_cache_dir: Path | None) -> None:
    global _worker_ref, _worker_chain_cache_dir
    _worker_ref = ref
    _worker_chain_cache_dir = chain_cache_dir


def _score_chain_file(part_file: Path) -> ScoredChain:
    """Worker entry point: process one per-chain dump CSV file."""
    global _worker_ref
    chain_id = part_file.stem  # entry_id = filename without .csv

    try:
        with part_file.open(newline="") as f:
            reader = csv.DictReader(f)
            fieldnames = set(reader.fieldnames or [])
            missing = sorted(REQUIRED_CANDIDATE_COLUMNS - fieldnames)
            candidates = list(reader) if not missing else []
    except Exception as exc:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=chain_id,
                    scope="chain",
                    code=RejectionCode.SCHEMA_MISMATCH,
                    detail=f"could not read candidate CSV ({type(exc).__name__})",
                )
            ],
        )

    if missing:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=chain_id,
                    scope="chain",
                    code=RejectionCode.SCHEMA_MISMATCH,
                    detail=f"candidate CSV is missing columns: {', '.join(missing)}",
                )
            ],
        )

    if not candidates:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=chain_id,
                    scope="chain",
                    code=RejectionCode.SCHEMA_MISMATCH,
                    detail="candidate dump contains no rows",
                )
            ],
        )

    observed_ids = {
        value.strip()
        for row in candidates
        for field in ("chain_id", "entry_id")
        if (value := row.get(field, "")) and value.strip()
    }
    if observed_ids and observed_ids != {chain_id}:
        return ScoredChain(
            rows=[],
            rejections=[
                RejectionRecord(
                    chain_id=chain_id,
                    scope="chain",
                    code=RejectionCode.SCHEMA_MISMATCH,
                    detail="candidate row identity disagrees with part filename",
                )
            ],
        )

    if _worker_ref is None:
        raise RuntimeError("training-table worker was not initialized")
    return _score_candidates(
        chain_id,
        candidates,
        _worker_ref,
        chain_cache_dir=_worker_chain_cache_dir,
    )


# ---------------------------------------------------------------------------
# Main processing functions
# ---------------------------------------------------------------------------

def write_rejections(path: Path, records: Sequence[RejectionRecord]) -> None:
    """Write a stable rejection manifest shared by both input modes."""
    ordered = sorted(
        records,
        key=lambda record: (
            record.chain_id,
            record.scope,
            record.code.value,
            record.delineation or "",
            record.detail,
        ),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=REJECTION_FIELDNAMES,
            lineterminator="\n",
        )
        writer.writeheader()
        for record in ordered:
            writer.writerow(
                {
                    "chain_id": record.chain_id,
                    "scope": record.scope,
                    "code": record.code.value,
                    "detail": record.detail,
                    "delineation": record.delineation or "",
                }
            )


def _default_rejections_path(output_path: Path) -> Path:
    return output_path.with_suffix(".rejections.csv")


def build_from_dump_dir(
    dump_dir: Path,
    reference: dict[str, CathEntry],
    output_path: Path,
    workers: int = 8,
    chain_cache_dir: Path | None = None,
    rejections_path: Path | None = None,
) -> None:
    """Stream per-chain dump files from a directory — fast, low-memory."""
    part_files = sorted(dump_dir.glob("*.csv"))
    log.info("Processing %d per-chain dump files with %d workers", len(part_files), workers)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_chains = n_rows = n_skipped = 0
    rejections: list[RejectionRecord] = []

    with output_path.open("w", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=FIELDNAMES)
        writer.writeheader()

        with ProcessPoolExecutor(
            max_workers=workers,
            initializer=_init_worker,
            initargs=(reference, chain_cache_dir),
        ) as pool:
            for result in pool.map(_score_chain_file, part_files, chunksize=20):
                rejections.extend(result.rejections)
                if not result.rows:
                    n_skipped += 1
                else:
                    n_chains += 1
                    n_rows += len(result.rows)
                    for row in result.rows:
                        writer.writerow(row)

                done = n_chains + n_skipped
                if done % 500 == 0:
                    log.info(
                        "Progress: %d / %d  (%d scored, %d skipped, %d rows)",
                        done, len(part_files), n_chains, n_skipped, n_rows,
                    )

    rejection_output = rejections_path or _default_rejections_path(output_path)
    write_rejections(rejection_output, rejections)
    log.info(
        "Done: %d chains scored, %d skipped, %d rows written, %d rejections",
        n_chains,
        n_skipped,
        n_rows,
        len(rejections),
    )


def build_from_dump_csv(
    dump_path: Path,
    reference: dict[str, CathEntry],
    output_path: Path,
    chain_cache_dir: Path | None = None,
    rejections_path: Path | None = None,
) -> None:
    """Process a merged dump CSV — streaming groupby (lower RAM than list+sort)."""
    log.info("Reading dump CSV: %s", dump_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    n_chains = n_rows = n_skipped = 0
    rejections: list[RejectionRecord] = []

    # Collect rows grouped by chain_id using a dict (single linear pass, no sort)
    chain_buckets: dict[str, list[dict[str, str]]] = {}
    with dump_path.open(newline="") as src:
        reader = csv.DictReader(src)
        fieldnames = set(reader.fieldnames or [])
        identity_fields = fieldnames.intersection({"chain_id", "entry_id"})
        missing = sorted(REQUIRED_CANDIDATE_COLUMNS - fieldnames)
        if not identity_fields:
            rejections.append(
                RejectionRecord(
                    chain_id="",
                    scope="chain",
                    code=RejectionCode.SCHEMA_MISMATCH,
                    detail="dump CSV has no chain_id or entry_id column",
                )
            )
        elif missing:
            rejections.append(
                RejectionRecord(
                    chain_id="",
                    scope="chain",
                    code=RejectionCode.SCHEMA_MISMATCH,
                    detail=f"dump CSV is missing columns: {', '.join(missing)}",
                )
            )
        else:
            for row_number, row in enumerate(reader, start=2):
                cid = (row.get("chain_id") or row.get("entry_id") or "").strip()
                if not cid:
                    rejections.append(
                        RejectionRecord(
                            chain_id="",
                            scope="chain",
                            code=RejectionCode.SCHEMA_MISMATCH,
                            detail=f"dump row {row_number} has no chain identity",
                        )
                    )
                    continue
                chain_buckets.setdefault(cid, []).append(row)

    log.info("Loaded %d chains from dump", len(chain_buckets))

    with output_path.open("w", newline="") as dst:
        writer = csv.DictWriter(dst, fieldnames=FIELDNAMES)
        writer.writeheader()

        for chain_id in sorted(chain_buckets):
            candidates = chain_buckets[chain_id]
            result = _score_candidates(
                chain_id,
                candidates,
                reference,
                chain_cache_dir=chain_cache_dir,
            )
            rejections.extend(result.rejections)
            if result.rows:
                n_chains += 1
                n_rows += len(result.rows)
                for row in result.rows:
                    writer.writerow(row)
            else:
                n_skipped += 1

            done = n_chains + n_skipped
            if done % 500 == 0:
                log.info(
                    "Progress: %d / %d  (%d scored, %d skipped, %d rows)",
                    done, len(chain_buckets), n_chains, n_skipped, n_rows,
                )

    rejection_output = rejections_path or _default_rejections_path(output_path)
    write_rejections(rejection_output, rejections)
    log.info(
        "Done: %d chains scored, %d skipped, %d rows written, %d rejections",
        n_chains,
        n_skipped,
        n_rows,
        len(rejections),
    )


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
    parser.add_argument(
        "--rejections",
        type=Path,
        default=None,
        help="Rejected chains/candidates CSV (default: <output stem>.rejections.csv)",
    )
    parser.add_argument(
        "--chain-cache-dir",
        type=Path,
        default=REPO / "benchmark/cache/chains",
        help="Canonical single-chain PDB cache used to map CATH author residue numbers",
    )
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
        build_from_dump_dir(
            args.dump_dir,
            reference,
            args.output,
            workers=args.workers,
            chain_cache_dir=args.chain_cache_dir,
            rejections_path=args.rejections,
        )
    else:
        build_from_dump_csv(
            args.dump,
            reference,
            args.output,
            chain_cache_dir=args.chain_cache_dir,
            rejections_path=args.rejections,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
