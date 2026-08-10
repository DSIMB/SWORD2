"""Normalize exact Rust factorized dumps and attach only Task 1 CATH labels."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

from benchmark.build_training_table import _load_reference, _score_candidates
from benchmark.candidate_geometry import CANDIDATE_GEOMETRY_FIELDS
from benchmark.datasets import CathEntry, dataset_path
from benchmark.dump_candidate_corpus import (
    _sha256,
    required_dump_fields,
    validate_dump_header,
    validate_dump_rows,
)
from benchmark.factorized_ranker.corpus import (
    ACQUISITION_REJECTION_FIELDS,
    CorpusPaths,
    feature_schema_hash,
    normalize_argv,
    write_corpus,
)
from benchmark.factorized_ranker.schema import (
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    GLOBAL_FEATURES,
    SEED,
)


REPO = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = REPO / "target/release/sword2"
DEFAULT_CHAIN_CACHE = REPO / "benchmark/cache/chains"

_DUMP_VALUE_OPTIONS = {
    "--dataset": "dataset",
    "--out": "out",
    "--parts-dir": "parts_dir",
    "--chain-cache-dir": "chain_cache_dir",
    "--binary": "binary",
    "--jobs": "jobs",
    "--threads": "threads",
    "--timeout": "timeout",
    "--limit": "limit",
    "--seed": "seed",
}
_DUMP_PATH_ROLES = {
    "out": "$OUT",
    "parts_dir": "$PARTS_DIR",
    "chain_cache_dir": "$CHAIN_CACHE_DIR",
    "binary": "$BINARY",
}
_DUMP_DEFAULTS: dict[str, Any] = {
    "dataset": "cath17287",
    "out": "$OUT",
    "parts_dir": None,
    "chain_cache_dir": "$CHAIN_CACHE_DIR",
    "binary": "$BINARY",
    "jobs": 16,
    "threads": 1,
    "timeout": 300,
    "limit": None,
    "seed": SEED,
    "resume": False,
}


def _ordinary_positive_integer(raw: str, option: str) -> int:
    try:
        value = int(raw)
    except ValueError as error:
        raise ValueError(f"dump argv option {option} is not an integer") from error
    if value <= 0 or str(value) != raw:
        raise ValueError(f"dump argv option {option} is not an ordinary positive integer")
    return value


def _parse_dump_argv_normalized(argv: Sequence[str]) -> dict[str, Any]:
    if not argv or argv[0] != "benchmark.dump_candidate_corpus":
        raise ValueError("dump argv program identity mismatch")
    parsed = dict(_DUMP_DEFAULTS)
    seen: set[str] = set()
    index = 1
    while index < len(argv):
        token = argv[index]
        if not isinstance(token, str) or not token.startswith("--"):
            raise ValueError("dump argv contains an unexpected positional token")
        option, separator, inline_value = token.partition("=")
        if option == "--resume":
            if separator or option in seen:
                raise ValueError("dump argv contains a duplicate or valued --resume option")
            seen.add(option)
            parsed["resume"] = True
            index += 1
            continue
        field = _DUMP_VALUE_OPTIONS.get(option)
        if field is None:
            raise ValueError(f"dump argv contains unknown option {option}")
        if option in seen:
            raise ValueError(f"dump argv contains duplicate option {option}")
        seen.add(option)
        if separator:
            raw_value = inline_value
        else:
            if index + 1 >= len(argv):
                raise ValueError(f"dump argv option {option} has no value")
            raw_value = argv[index + 1]
            index += 1
        if not isinstance(raw_value, str) or raw_value == "":
            raise ValueError(f"dump argv option {option} has an invalid value")
        if field in _DUMP_PATH_ROLES:
            expected_role = _DUMP_PATH_ROLES[field]
            if raw_value != expected_role:
                raise ValueError(f"dump argv option {option} has the wrong path role")
            parsed[field] = raw_value
        elif field in {"jobs", "threads", "timeout", "limit", "seed"}:
            parsed[field] = _ordinary_positive_integer(raw_value, option)
        else:
            parsed[field] = raw_value
        index += 1
    return parsed


def _expected_git_commit() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _format_vector(row: Mapping[str, str], fields: Sequence[str]) -> tuple[str, ...]:
    result = []
    for field in fields:
        try:
            value = float(row[field])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f"missing or malformed feature {field}") from error
        if not math.isfinite(value):
            raise ValueError(f"non-finite feature {field}")
        result.append(format(value, ".17g"))
    return tuple(result)


def _schema_rejection(chain_id: str, detail: str) -> dict[str, str]:
    return {
        "chain_id": chain_id,
        "scope": "chain",
        "code": "schema_mismatch",
        "detail": detail,
        "delineation": "",
    }


def validate_raw_population(
    chain_id: str, rows: Sequence[Mapping[str, str]]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Validate repeated globals/counts before selecting a representative."""
    if not rows:
        raise ValueError("raw chain has no candidates")
    global_vectors = {_format_vector(row, GLOBAL_FEATURES) for row in rows}
    if len(global_vectors) != 1:
        raise ValueError("conflicting repeated global fields")
    count_vectors: dict[int, set[tuple[str, ...]]] = {}
    candidate_counts: set[int] = set()
    for row in rows:
        count = int(float(row["num_domains"]))
        candidate_counts.add(count)
        count_vectors.setdefault(count, set()).add(_format_vector(row, COUNT_ITEM_FEATURES))
    if any(len(vectors) != 1 for vectors in count_vectors.values()):
        raise ValueError("conflicting repeated count fields")
    count_keys = {int(float(next(iter(vectors))[0])) for vectors in count_vectors.values()}
    if count_keys != candidate_counts or set(count_vectors) != candidate_counts:
        raise ValueError("count rows do not match candidate count groups")

    chain_row: dict[str, Any] = {
        "chain_id": chain_id,
        **dict(zip(GLOBAL_FEATURES, next(iter(global_vectors)))),
    }
    count_rows: list[dict[str, Any]] = []
    for count in sorted(count_vectors):
        count_rows.append(
            {
                "chain_id": chain_id,
                **dict(zip(COUNT_ITEM_FEATURES, next(iter(count_vectors[count])))),
            }
        )
    return chain_row, count_rows


def _adapter_row(row: Mapping[str, str]) -> dict[str, str]:
    fields = (
        "min_size",
        "max_cr",
        "density_min",
        "mean_density",
        "boundary_coil_fraction",
        "modal_count_distance",
        *CANDIDATE_GEOMETRY_FIELDS,
    )
    return {
        "num_domains": row["num_domains"],
        "delineation": row["canonical_delineation"],
        "energy_z": "",
        **{field: row[field] for field in fields},
    }


def normalize_raw_rows(
    raw_rows: Sequence[Mapping[str, str]],
    reference: Mapping[str, CathEntry],
    chain_cache_dir: Path,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[Any],
]:
    """Validate raw populations, score by canonical key, and preserve Rust features."""
    canonical_groups: dict[str, list[dict[str, str]]] = {}
    unresolved_groups: dict[str, list[dict[str, str]]] = {}
    for raw in raw_rows:
        row = dict(raw)
        raw_chain_id = row["chain_id"]
        entry = reference.get(raw_chain_id)
        if entry is None:
            unresolved_groups.setdefault(raw_chain_id, []).append(row)
        else:
            canonical_groups.setdefault(entry.entry_id, []).append(row)

    chain_rows: list[dict[str, Any]] = []
    count_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    rejections: list[Any] = []

    for chain_id in sorted(unresolved_groups):
        scored = _score_candidates(
            chain_id, unresolved_groups[chain_id], dict(reference), chain_cache_dir=chain_cache_dir
        )
        rejections.extend(scored.rejections)

    for canonical_chain_id in sorted(canonical_groups):
        rows = canonical_groups[canonical_chain_id]
        if len({row["chain_id"] for row in rows}) != 1:
            rejections.append(
                _schema_rejection(
                    canonical_chain_id,
                    "canonical chain is sourced from multiple raw chain identities",
                )
            )
            continue
        try:
            base_chain, base_counts = validate_raw_population(canonical_chain_id, rows)
        except ValueError as error:
            rejections.append(_schema_rejection(canonical_chain_id, str(error)))
            continue

        scored = _score_candidates(
            canonical_chain_id,
            [_adapter_row(row) for row in rows],
            dict(reference),
            chain_cache_dir=chain_cache_dir,
        )
        rejections.extend(scored.rejections)
        chain_rejections = [
            rejection
            for rejection in scored.rejections
            if getattr(rejection, "scope", None) == "chain"
            or (isinstance(rejection, Mapping) and rejection.get("scope") == "chain")
        ]
        if chain_rejections:
            if scored.rows:
                rejections.append(
                    _schema_rejection(
                        canonical_chain_id,
                        "Task 1 scoring identity join is not one-to-one",
                    )
                )
            continue
        scored_by_canonical: dict[str, Mapping[str, Any]] = {}
        duplicate = False
        for scored_row in scored.rows:
            canonical = str(scored_row["delineation"])
            if canonical in scored_by_canonical:
                duplicate = True
                break
            scored_by_canonical[canonical] = scored_row
        raw_by_canonical = {row["canonical_delineation"]: row for row in rows}
        if len(raw_by_canonical) != len(rows):
            duplicate = True
        rejected_canonicals: list[str] = []
        for rejection in scored.rejections:
            if isinstance(rejection, Mapping):
                scope = rejection.get("scope")
                delineation = rejection.get("delineation")
            else:
                scope = getattr(rejection, "scope", None)
                delineation = getattr(rejection, "delineation", None)
            if scope == "candidate" and delineation:
                rejected_canonicals.append(str(delineation))
        scored_identities = set(scored_by_canonical)
        rejected_identities = set(rejected_canonicals)
        raw_identities = set(raw_by_canonical)
        if (
            duplicate
            or len(rejected_identities) != len(rejected_canonicals)
            or scored_identities & rejected_identities
            or scored_identities | rejected_identities != raw_identities
        ):
            rejections.append(
                _schema_rejection(canonical_chain_id, "Task 1 scoring identity join is not one-to-one")
            )
            continue
        if not scored_by_canonical:
            continue

        n_true_values = {int(float(row["n_true_domains"])) for row in scored_by_canonical.values()}
        if len(n_true_values) != 1:
            rejections.append(_schema_rejection(canonical_chain_id, "Task 1 true-domain labels disagree"))
            continue
        n_true = next(iter(n_true_values))
        base_chain["n_true_domains"] = n_true
        chain_rows.append(base_chain)
        count_rows.extend(base_counts)

        for canonical in sorted(scored_by_canonical):
            raw = raw_by_canonical[canonical]
            labels = scored_by_canonical[canonical]
            candidate_rows.append(
                {
                    "chain_id": canonical_chain_id,
                    "canonical_delineation": canonical,
                    "source_index": raw["source_index"],
                    "legacy_distance": raw["legacy_distance"],
                    **{field: raw[field] for field in CANDIDATE_FEATURES},
                    **{
                        field: labels[field]
                        for field in (
                            "n_true_domains",
                            "n_pred_domains",
                            "ndo",
                            "iou",
                            "boundary_f1_10",
                            "matched_dice",
                            "d_count_acc",
                            "S",
                            "is_oracle_s",
                        )
                    },
                }
            )
    return chain_rows, count_rows, candidate_rows, rejections


def _read_raw_csv(
    path: Path,
    expected_chain_id: str | None = None,
    *,
    allow_empty: bool = False,
) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        validate_dump_header(reader.fieldnames)
        rows = list(reader)
    if allow_empty and not rows:
        return []
    return validate_dump_rows(rows, expected_chain_id)


def _load_raw(dump: Path | None, dump_dir: Path | None) -> list[dict[str, str]]:
    if dump is not None:
        return _read_raw_csv(dump, allow_empty=True)
    assert dump_dir is not None
    rows: list[dict[str, str]] = []
    for part in sorted(dump_dir.glob("*.csv")):
        if part.name == "rejections.csv":
            continue
        rows.extend(_read_raw_csv(part, part.stem))
    return sorted(
        rows,
        key=lambda row: (
            row["chain_id"],
            int(float(row["num_domains"])),
            row["canonical_delineation"],
            int(row["source_index"]),
        ),
    )


def _read_rejections(path: Path | None, explicit: bool) -> list[dict[str, str]]:
    if path is None or not path.exists():
        if explicit:
            raise FileNotFoundError("explicit acquisition rejection input is missing")
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != ACQUISITION_REJECTION_FIELDS:
            raise ValueError("acquisition rejection input has the wrong schema")
        return [dict(row) for row in reader]


def _load_and_verify_provenance(path: Path, dataset: str, binary: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ValueError("acquisition provenance sidecar is missing")
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError("acquisition provenance must be a JSON object")
    required_fields = {
        "dataset",
        "dataset_sha256",
        "binary_sha256",
        "git_commit",
        "seed",
        "jobs",
        "threads",
        "timeout",
        "limit",
        "resume",
        "feature_schema_hash",
        "dump_argv_normalized",
    }
    if set(value) != required_fields:
        raise ValueError("acquisition provenance schema mismatch")
    integer_fields = ("seed", "jobs", "threads", "timeout")
    if any(type(value.get(field)) is not int or value[field] <= 0 for field in integer_fields):
        raise ValueError("acquisition provenance contains an invalid integer setting")
    if value.get("limit") is not None and (
        type(value["limit"]) is not int or value["limit"] <= 0
    ):
        raise ValueError("acquisition provenance limit is invalid")
    if type(value.get("resume")) is not bool:
        raise ValueError("acquisition provenance resume is invalid")
    expected = {
        "dataset": dataset,
        "dataset_sha256": _sha256(dataset_path(dataset)),
        "binary_sha256": _sha256(binary),
        "git_commit": _expected_git_commit(),
        "feature_schema_hash": feature_schema_hash(),
    }
    for key, expected_value in expected.items():
        if value.get(key) != expected_value:
            raise ValueError(f"acquisition provenance {key} mismatch")
    dump_argv = value.get("dump_argv_normalized")
    if not isinstance(dump_argv, list) or not all(isinstance(token, str) for token in dump_argv):
        raise ValueError("acquisition provenance dump argv is invalid")
    parsed_argv = _parse_dump_argv_normalized(dump_argv)
    for field in ("dataset", "seed", "jobs", "threads", "timeout", "limit", "resume"):
        if parsed_argv[field] != value.get(field):
            raise ValueError(f"acquisition provenance {field} does not match dump argv")
    if value["seed"] != SEED:
        raise ValueError("acquisition provenance seed mismatch")
    if parsed_argv["binary"] != "$BINARY":
        raise ValueError("acquisition provenance binary path role mismatch")
    return value


def _provenance_path(dump: Path | None, dump_dir: Path | None) -> Path:
    if dump is not None:
        return dump.with_suffix(".provenance.json")
    assert dump_dir is not None
    return dump_dir / "provenance.json"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--dump", type=Path)
    source.add_argument("--dump-dir", type=Path)
    parser.add_argument("--dataset", default="cath17287")
    parser.add_argument("--chain-cache-dir", type=Path, default=DEFAULT_CHAIN_CACHE)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--rejections", type=Path, default=None)
    args = parser.parse_args(argv)

    if args.dataset != "cath17287":
        parser.error("the factorized development corpus accepts only dataset cath17287")
    if not args.binary.is_file():
        parser.error("SWORD2 binary is unavailable")

    provenance = _load_and_verify_provenance(
        _provenance_path(args.dump, args.dump_dir), args.dataset, args.binary
    )
    raw_rows = _load_raw(args.dump, args.dump_dir)
    explicit_rejections = args.rejections is not None
    if args.rejections is not None:
        rejection_input = args.rejections
    elif args.dump is not None:
        rejection_input = args.dump.with_suffix(".rejections.csv")
    else:
        rejection_input = args.dump_dir / "rejections.csv"
    acquisition_rejections = _read_rejections(rejection_input, explicit_rejections)

    reference = _load_reference(args.dataset)
    chains, counts, candidates, scoring_rejections = normalize_raw_rows(
        raw_rows, reference, args.chain_cache_dir
    )
    raw_argv = [
        "benchmark.build_factorized_corpus",
        *(sys.argv[1:] if argv is None else argv),
    ]
    normalized_provenance = {
        "dataset": provenance["dataset"],
        "dataset_sha256": provenance["dataset_sha256"],
        "binary_sha256": provenance["binary_sha256"],
        "git_commit": provenance["git_commit"],
        "dump_argv_normalized": provenance["dump_argv_normalized"],
        "build_argv_normalized": normalize_argv(raw_argv),
        "seed": provenance["seed"],
        "jobs": provenance["jobs"],
        "threads": provenance["threads"],
        "timeout": provenance["timeout"],
        "limit": provenance["limit"],
        "resume": provenance["resume"],
        "feature_schema_hash": provenance["feature_schema_hash"],
    }
    write_corpus(
        CorpusPaths.at(args.out_dir),
        chains,
        counts,
        candidates,
        [*acquisition_rejections, *scoring_rejections],
        normalized_provenance,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
