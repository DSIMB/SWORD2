"""Acquire exact Rust factorized feature lattices into a resumable raw corpus."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Mapping, Sequence

from benchmark.datasets import dataset_path, load_dataset
from benchmark.factorized_ranker.corpus import (
    ACQUISITION_REJECTION_FIELDS,
    feature_schema_hash,
    normalize_argv,
)
from benchmark.factorized_ranker.schema import (
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    GLOBAL_FEATURES,
    SEED,
)


REPO = Path(__file__).resolve().parents[1]
BINARY = REPO / "target/release/sword2"
CHAINS = REPO / "benchmark/cache/chains"
_PROHIBITED_PREDICTOR_TOKENS = ("merizo", "chainsaw")


def required_dump_fields() -> list[str]:
    return [
        "chain_id",
        "canonical_delineation",
        "source_index",
        "legacy_distance",
        *GLOBAL_FEATURES,
        *COUNT_ITEM_FEATURES,
        *CANDIDATE_FEATURES,
    ]


OUTPUT_FIELDS = required_dump_fields()
_FLOAT_FIELDS = {"legacy_distance", *GLOBAL_FEATURES, *COUNT_ITEM_FEATURES, *CANDIDATE_FEATURES}


def validate_dump_header(fields: Sequence[str] | None) -> None:
    observed = list(fields or [])
    duplicates = sorted({field for field in observed if observed.count(field) > 1})
    if duplicates:
        raise ValueError(f"duplicate dump header names: {', '.join(duplicates)}")
    prohibited = sorted(
        field
        for field in observed
        if any(token in field.casefold() for token in _PROHIBITED_PREDICTOR_TOKENS)
    )
    if prohibited:
        raise ValueError(f"external predictor columns are prohibited: {', '.join(prohibited)}")
    expected = required_dump_fields()
    if observed != expected:
        missing = [field for field in expected if field not in observed]
        unknown = [field for field in observed if field not in expected]
        detail = []
        if missing:
            detail.append(f"missing {', '.join(missing)}")
        if unknown:
            detail.append(f"unknown {', '.join(unknown)}")
        if not detail:
            detail.append("field order differs")
        raise ValueError("dump header must match exact frozen order (" + "; ".join(detail) + ")")


def _finite(raw: str, field: str) -> float:
    try:
        value = float(raw)
    except (TypeError, ValueError) as error:
        raise ValueError(f"field {field} is malformed") from error
    if not math.isfinite(value):
        raise ValueError(f"field {field} is non-finite")
    return value


def validate_dump_rows(
    rows: Sequence[Mapping[str, str]], expected_chain_id: str | None = None
) -> list[dict[str, str]]:
    if not rows:
        raise ValueError("candidate dump is header-only")
    validated: list[dict[str, str]] = []
    identities: set[tuple[int, int, str]] = set()
    count_canonicals: set[tuple[int, str]] = set()
    exact_fields = set(required_dump_fields())
    for raw in rows:
        if set(raw) != exact_fields:
            raise ValueError("candidate row columns do not match the exact dump header")
        row = {field: raw.get(field, "") for field in required_dump_fields()}
        chain_id = row["chain_id"]
        if not chain_id or (expected_chain_id is not None and chain_id != expected_chain_id):
            raise ValueError("candidate row identity disagrees with part filename")
        canonical = row["canonical_delineation"]
        if not canonical:
            raise ValueError("canonical delineation is empty")
        try:
            source_index = int(row["source_index"])
        except ValueError as error:
            raise ValueError("source_index is not an ordinary nonnegative integer") from error
        if source_index < 0 or str(source_index) != row["source_index"]:
            raise ValueError("source_index is not an ordinary nonnegative integer")
        for field in _FLOAT_FIELDS:
            _finite(row[field], field)
        num_domains_value = _finite(row["num_domains"], "num_domains")
        if num_domains_value <= 0 or not num_domains_value.is_integer():
            raise ValueError("num_domains is not a positive integer")
        num_domains = int(num_domains_value)
        count_domains = _finite(row["count_num_domains"], "count_num_domains")
        if count_domains != num_domains:
            raise ValueError("candidate/count domain identity mismatch")
        identity = (source_index, num_domains, canonical)
        if identity in identities or (num_domains, canonical) in count_canonicals:
            raise ValueError("duplicate candidate identity")
        identities.add(identity)
        count_canonicals.add((num_domains, canonical))
        validated.append(row)
    return sorted(
        validated,
        key=lambda row: (
            int(float(row["num_domains"])),
            row["canonical_delineation"],
            int(row["source_index"]),
        ),
    )


def _read_validated(path: Path, expected_chain_id: str | None = None) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        validate_dump_header(reader.fieldnames)
        return validate_dump_rows(list(reader), expected_chain_id)


def valid_resume_part(path: Path, expected_chain_id: str) -> bool:
    try:
        _read_validated(path, expected_chain_id)
    except (OSError, csv.Error, ValueError):
        return False
    return True


def _atomic_write_csv(path: Path, rows: Sequence[Mapping[str, str]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _rejection(chain_id: str, code: str, detail: str) -> dict[str, str]:
    return {
        "chain_id": chain_id,
        "scope": "chain",
        "code": code,
        "detail": detail,
        "delineation": "",
    }


def _run_one(
    entry_id: str,
    binary: Path,
    chain_dir: Path,
    parts_dir: Path,
    scratch_dir: Path,
    threads: int,
    timeout: int,
) -> tuple[str, dict[str, str] | None]:
    part_path = parts_dir / f"{entry_id}.csv"
    # A failed refresh must never leave an old invalid part eligible for merge.
    part_path.unlink(missing_ok=True)
    pdb_path = chain_dir / f"{entry_id}.pdb"
    if not pdb_path.is_file():
        return entry_id, _rejection(entry_id, "missing_chain_pdb", "canonical chain PDB is unavailable")

    work_dir = Path(tempfile.mkdtemp(prefix=f"{entry_id}-", dir=scratch_dir))
    try:
        raw_dump = work_dir / "candidates.csv"
        env = dict(os.environ, SWORD2_DUMP_CANDIDATES=str(raw_dump), RUST_LOG="error")
        try:
            result = subprocess.run(
                [
                    str(binary),
                    "-i",
                    str(pdb_path),
                    "-o",
                    str(work_dir / "out"),
                    "--threads",
                    str(threads),
                ],
                cwd=REPO,
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return entry_id, _rejection(entry_id, "schema_mismatch", f"SWORD2 timed out after {timeout} seconds")
        if result.returncode != 0:
            return entry_id, _rejection(entry_id, "schema_mismatch", f"SWORD2 exited with status {result.returncode}")
        if not raw_dump.is_file():
            return entry_id, _rejection(entry_id, "schema_mismatch", "Rust candidate dump is missing")
        try:
            rows = _read_validated(raw_dump, entry_id)
        except (OSError, csv.Error, ValueError) as error:
            return entry_id, _rejection(entry_id, "schema_mismatch", f"Rust candidate dump is invalid ({type(error).__name__})")
        _atomic_write_csv(part_path, rows, required_dump_fields())
        return entry_id, None
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _merge_parts(entries: Sequence[str], parts_dir: Path, output: Path) -> tuple[int, int]:
    merged: list[dict[str, str]] = []
    chains = 0
    for entry_id in entries:
        part = parts_dir / f"{entry_id}.csv"
        if not part.exists():
            continue
        rows = _read_validated(part, entry_id)
        merged.extend(rows)
        chains += 1
    merged.sort(
        key=lambda row: (
            row["chain_id"],
            int(float(row["num_domains"])),
            row["canonical_delineation"],
            int(row["source_index"]),
        )
    )
    _atomic_write_csv(output, merged, required_dump_fields())
    return chains, len(merged)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _atomic_json(path: Path, value: Mapping[str, object]) -> None:
    data = (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _write_rejections(path: Path, rejections: Sequence[Mapping[str, str]]) -> None:
    unique = {
        (row["chain_id"], row["scope"], row["code"], row.get("delineation", ""), row["detail"]): row
        for row in rejections
    }
    _atomic_write_csv(path, [unique[key] for key in sorted(unique)], ACQUISITION_REJECTION_FIELDS)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default="cath17287")
    parser.add_argument("--out", type=Path, default=REPO / "benchmark/data/cath17287_factorized.csv")
    parser.add_argument("--parts-dir", type=Path, default=None)
    parser.add_argument("--chain-dir", type=Path, default=CHAINS)
    parser.add_argument("--binary", type=Path, default=BINARY)
    parser.add_argument("--jobs", type=int, default=16)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)

    if args.dataset != "cath17287":
        parser.error("the factorized development corpus accepts only dataset cath17287")
    if args.seed != SEED:
        parser.error(f"the factorized corpus seed is frozen at {SEED}")
    if args.jobs <= 0 or args.threads <= 0 or args.timeout <= 0:
        parser.error("jobs, threads, and timeout must be positive")
    if not args.binary.is_file():
        parser.error("SWORD2 binary is unavailable")

    entries = load_dataset(args.dataset)
    random.Random(args.seed).shuffle(entries)
    if args.limit is not None:
        entries = entries[: args.limit]
    entry_ids = [entry.entry_id for entry in entries]
    parts_dir = args.parts_dir or args.out.with_suffix("").with_name(args.out.stem + "_parts")
    parts_dir.mkdir(parents=True, exist_ok=True)
    pending = [
        entry_id
        for entry_id in entry_ids
        if not (args.resume and valid_resume_part(parts_dir / f"{entry_id}.csv", entry_id))
    ]
    scratch_dir = Path(tempfile.mkdtemp(prefix="sword2-factorized-corpus-"))
    failures: list[dict[str, str]] = []
    try:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(
                    _run_one,
                    entry_id,
                    args.binary,
                    args.chain_dir,
                    parts_dir,
                    scratch_dir,
                    args.threads,
                    args.timeout,
                ): entry_id
                for entry_id in pending
            }
            for completed, future in enumerate(as_completed(futures), start=1):
                _, rejection = future.result()
                if rejection is not None:
                    failures.append(rejection)
                if completed % 100 == 0 or completed == len(pending):
                    print(f"processed {completed}/{len(pending)}; failures={len(failures)}", file=sys.stderr)
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)

    chains, rows = _merge_parts(entry_ids, parts_dir, args.out)
    rejection_path = args.out.with_suffix(".rejections.csv")
    _write_rejections(rejection_path, failures)
    raw_argv = [
        "benchmark.dump_candidate_corpus",
        *(sys.argv[1:] if argv is None else argv),
    ]
    provenance = {
        "dataset": args.dataset,
        "dataset_sha256": _sha256(dataset_path(args.dataset)),
        "binary_sha256": _sha256(args.binary),
        "git_commit": _git_commit(),
        "seed": args.seed,
        "jobs": args.jobs,
        "threads": args.threads,
        "feature_schema_hash": feature_schema_hash(),
        "dump_argv_normalized": normalize_argv(raw_argv),
    }
    _atomic_json(args.out.with_suffix(".provenance.json"), provenance)
    _atomic_json(parts_dir / "provenance.json", provenance)
    _write_rejections(parts_dir / "rejections.csv", failures)
    print(f"wrote {rows} candidates from {chains}/{len(entry_ids)} chains to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
