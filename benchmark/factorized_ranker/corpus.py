"""Deterministic normalized corpus schemas and atomic writers."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from benchmark.factorized_ranker.schema import (
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    GLOBAL_FEATURES,
    SCHEMA_VERSION,
)


ACQUISITION_REJECTION_FIELDS = NORMALIZED_REJECTION_FIELDS = (
    "chain_id",
    "scope",
    "code",
    "detail",
    "delineation",
)
CHAIN_FIELDS = ("chain_id", "n_true_domains", *GLOBAL_FEATURES)
COUNT_FIELDS = ("chain_id", *COUNT_ITEM_FEATURES)
CANDIDATE_FIELDS = (
    "candidate_id",
    "chain_id",
    "canonical_delineation",
    "source_index",
    "legacy_distance",
    *CANDIDATE_FEATURES,
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

_PATH_OPTIONS = {
    "--dump": "$DUMP",
    "--dump-dir": "$DUMP_DIR",
    "--chain-cache-dir": "$CHAIN_CACHE_DIR",
    "--out-dir": "$OUT_DIR",
    "--binary": "$BINARY",
    "--rejections": "$REJECTIONS",
    "--parts-dir": "$PARTS_DIR",
    "--out": "$OUT",
}
_INTEGER_FIELDS = {
    "source_index",
    "n_true_domains",
    "n_pred_domains",
    "is_oracle_s",
}
_STRING_FIELDS = {
    "candidate_id",
    "chain_id",
    "canonical_delineation",
    "scope",
    "code",
    "detail",
    "delineation",
}


@dataclass(frozen=True)
class CorpusPaths:
    chains: Path
    counts: Path
    candidates: Path
    rejections: Path
    manifest: Path

    @classmethod
    def at(cls, directory: Path) -> "CorpusPaths":
        directory = Path(directory)
        return cls(
            chains=directory / "chains.csv",
            counts=directory / "counts.csv",
            candidates=directory / "candidates.csv",
            rejections=directory / "rejections.csv",
            manifest=directory / "corpus_manifest.json",
        )


def _canonical_json(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")


def feature_schema_document() -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "global_features": list(GLOBAL_FEATURES),
        "count_item_features": list(COUNT_ITEM_FEATURES),
        "candidate_features": list(CANDIDATE_FEATURES),
    }


def feature_schema_hash() -> str:
    return hashlib.sha256(_canonical_json(feature_schema_document()).rstrip(b"\n")).hexdigest()


def manifest_hash(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def normalize_argv(argv: Sequence[str]) -> list[str]:
    """Replace path option values with stable role tokens."""
    normalized: list[str] = []
    index = 0
    while index < len(argv):
        token = str(argv[index])
        option, separator, _value = token.partition("=")
        if separator and option in _PATH_OPTIONS:
            normalized.append(f"{option}={_PATH_OPTIONS[option]}")
            index += 1
            continue
        normalized.append(token)
        if token in _PATH_OPTIONS:
            if index + 1 >= len(argv):
                raise ValueError(f"path option {token} has no value")
            normalized.append(_PATH_OPTIONS[token])
            index += 2
        else:
            index += 1
    return normalized


def _raw_value(row: Mapping[str, Any], field: str) -> Any:
    if field not in row or row[field] is None or row[field] == "":
        raise ValueError(f"missing required value {field!r}")
    return row[field]


def _format_integer(value: Any, field: str) -> str:
    if isinstance(value, bool):
        number = int(value)
    else:
        try:
            numeric = float(value)
        except (TypeError, ValueError) as error:
            raise ValueError(f"field {field!r} is not an integer") from error
        if not math.isfinite(numeric) or not numeric.is_integer():
            raise ValueError(f"field {field!r} is not a finite integer")
        number = int(numeric)
    return str(number)


def _format_float(value: Any, field: str) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"field {field!r} is not finite") from error
    if not math.isfinite(numeric):
        raise ValueError(f"field {field!r} is not finite")
    return format(numeric, ".17g")


def _format_row(row: Mapping[str, Any], fields: Sequence[str]) -> dict[str, str]:
    formatted: dict[str, str] = {}
    for field in fields:
        value = _raw_value(row, field)
        if field in _STRING_FIELDS:
            formatted[field] = str(value)
        elif field in _INTEGER_FIELDS:
            formatted[field] = _format_integer(value, field)
        else:
            formatted[field] = _format_float(value, field)
    return formatted


def _deduplicate_identical(
    rows: Iterable[Mapping[str, Any]],
    fields: Sequence[str],
    key_fields: Sequence[str],
    description: str,
) -> list[dict[str, str]]:
    unique: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        formatted = _format_row(row, fields)
        key = tuple(formatted[field] for field in key_fields)
        previous = unique.get(key)
        if previous is not None and previous != formatted:
            raise ValueError(f"conflicting repeated {description} fields")
        unique[key] = formatted
    return list(unique.values())


def _rejection_row(value: Any) -> dict[str, str]:
    if isinstance(value, Mapping):
        source = value
    else:
        source = {field: getattr(value, field, "") for field in NORMALIZED_REJECTION_FIELDS}
    row: dict[str, str] = {}
    for field in NORMALIZED_REJECTION_FIELDS:
        raw = source.get(field, "")
        if isinstance(raw, Enum):
            raw = raw.value
        row[field] = "" if raw is None else str(raw)
    if row["scope"] not in {"chain", "candidate"}:
        raise ValueError("rejection scope must be chain or candidate")
    if not row["code"]:
        raise ValueError("rejection code is required")
    return row


def _atomic_write(path: Path, data: bytes) -> None:
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


def _csv_bytes(fields: Sequence[str], rows: Sequence[Mapping[str, str]]) -> bytes:
    from io import StringIO

    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8")


def _validate_manifest_tree(value: object) -> None:
    forbidden = {"timestamp", "hostname", "username", "manifest_sha256"}
    if isinstance(value, Mapping):
        for key, child in value.items():
            if str(key) in forbidden:
                raise ValueError("manifest provenance contains a forbidden key")
            _validate_manifest_tree(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            _validate_manifest_tree(child)
    elif isinstance(value, str) and value.startswith("/"):
        raise ValueError("manifest provenance contains an absolute path")


def write_corpus(
    paths: CorpusPaths,
    chain_rows: Iterable[Mapping[str, Any]],
    count_rows: Iterable[Mapping[str, Any]],
    candidate_rows: Iterable[Mapping[str, Any]],
    rejections: Iterable[Any],
    provenance: Mapping[str, Any],
) -> dict[str, object]:
    """Validate and atomically write the four normalized tables and manifest."""
    _validate_manifest_tree(dict(provenance))
    chains = _deduplicate_identical(chain_rows, CHAIN_FIELDS, ("chain_id",), "global")
    counts = _deduplicate_identical(
        count_rows,
        COUNT_FIELDS,
        ("chain_id", "count_num_domains"),
        "count",
    )

    candidates: list[dict[str, str]] = []
    candidate_keys: set[tuple[str, str]] = set()
    true_domains: dict[str, str] = {row["chain_id"]: row["n_true_domains"] for row in chains}
    for raw in candidate_rows:
        row_input = dict(raw)
        chain_id = str(_raw_value(row_input, "chain_id"))
        canonical = str(_raw_value(row_input, "canonical_delineation"))
        candidate_id = hashlib.sha256((chain_id + "\0" + canonical).encode("utf-8")).hexdigest()
        if row_input.get("candidate_id") not in (None, "", candidate_id):
            raise ValueError("candidate_id does not match canonical identity")
        row_input["candidate_id"] = candidate_id
        row = _format_row(row_input, CANDIDATE_FIELDS)
        key = (row["chain_id"], row["canonical_delineation"])
        if key in candidate_keys:
            raise ValueError("duplicate candidate canonical identity")
        candidate_keys.add(key)
        if chain_id not in true_domains or true_domains[chain_id] != row["n_true_domains"]:
            raise ValueError("n_true_domains disagrees across accepted chain rows")
        candidates.append(row)

    chain_ids = {row["chain_id"] for row in chains}
    if {row["chain_id"] for row in counts} != chain_ids:
        raise ValueError("count rows do not cover exactly the accepted chains")
    if {row["chain_id"] for row in candidates} != chain_ids:
        raise ValueError("candidate rows do not cover exactly the accepted chains")
    count_keys = {(row["chain_id"], int(float(row["count_num_domains"]))) for row in counts}
    candidate_count_keys = {
        (row["chain_id"], int(float(row["num_domains"]))) for row in candidates
    }
    if not candidate_count_keys.issubset(count_keys):
        raise ValueError("accepted candidate count has no complete-lattice count row")

    chains.sort(key=lambda row: row["chain_id"])
    counts.sort(key=lambda row: (row["chain_id"], int(float(row["count_num_domains"]))))
    candidates.sort(
        key=lambda row: (
            row["chain_id"],
            int(float(row["num_domains"])),
            row["canonical_delineation"],
            int(row["source_index"]),
        )
    )
    rejection_map: dict[tuple[str, str, str, str, str], dict[str, str]] = {}
    for value in rejections:
        row = _rejection_row(value)
        key = (row["chain_id"], row["scope"], row["code"], row["delineation"], row["detail"])
        rejection_map[key] = row
    ordered_rejections = [rejection_map[key] for key in sorted(rejection_map)]

    table_data = {
        "chains": _csv_bytes(CHAIN_FIELDS, chains),
        "counts": _csv_bytes(COUNT_FIELDS, counts),
        "candidates": _csv_bytes(CANDIDATE_FIELDS, candidates),
        "rejections": _csv_bytes(NORMALIZED_REJECTION_FIELDS, ordered_rejections),
    }
    for name, data in table_data.items():
        _atomic_write(getattr(paths, name), data)

    rejected_chains = {row["chain_id"] for row in ordered_rejections if row["chain_id"]}
    candidate_rejection_counts: dict[str, int] = {}
    for row in ordered_rejections:
        if row["scope"] == "candidate":
            candidate_rejection_counts[row["code"]] = candidate_rejection_counts.get(row["code"], 0) + 1
    tables = {
        name: {"sha256": hashlib.sha256(data).hexdigest(), "rows": len(rows)}
        for name, data, rows in (
            ("chains", table_data["chains"], chains),
            ("counts", table_data["counts"], counts),
            ("candidates", table_data["candidates"], candidates),
            ("rejections", table_data["rejections"], ordered_rejections),
        )
    }
    manifest: dict[str, object] = {
        **dict(provenance),
        "feature_schema_version": SCHEMA_VERSION,
        "feature_schema_hash": feature_schema_hash(),
        "tables": tables,
        "accepted_unique_chains": len(chain_ids),
        "rejected_unique_chains": len(rejected_chains),
        "candidate_rejection_counts_by_code": dict(sorted(candidate_rejection_counts.items())),
    }
    _validate_manifest_tree(manifest)
    encoded_manifest = _canonical_json(manifest)
    _atomic_write(paths.manifest, encoded_manifest)
    return manifest
