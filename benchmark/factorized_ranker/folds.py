"""Leakage-resistant structural folds and canonical fold manifests."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from dataclasses import dataclass
from io import StringIO
from numbers import Integral
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

from benchmark.datasets import (
    CathEntry,
    cath_family_combination,
    cath_family_labels,
    dataset_path,
    load_dataset,
)
from benchmark.factorized_ranker.corpus import (
    CANDIDATE_FIELDS,
    CHAIN_FIELDS,
    COUNT_FIELDS,
    NORMALIZED_REJECTION_FIELDS,
    feature_schema_hash,
)
from benchmark.factorized_ranker.schema import GLOBAL_FEATURES, SCHEMA_VERSION


CountBin = Literal["1", "2", "3", "4", "5+"]
LengthBin = Literal["<250", "250-349", "350-449", "450+"]
LabelCohort = Literal["seen", "unseen", "unknown"]

COUNT_BINS: tuple[CountBin, ...] = ("1", "2", "3", "4", "5+")
LENGTH_BINS: tuple[LengthBin, ...] = ("<250", "250-349", "350-449", "450+")
FOLD_MANIFEST_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class FoldAssignment:
    chain_id: str
    pdb_id: str
    family_combination: tuple[str, ...]
    family_labels: tuple[str, ...]
    true_count_bin: CountBin
    length_bin: LengthBin
    component_id: str
    fold: int


@dataclass(frozen=True)
class FoldManifest:
    schema_version: int
    seed: int
    n_folds: int
    dataset_sha256: str
    corpus_manifest_sha256: str
    chains_sha256: str
    assignments: tuple[FoldAssignment, ...]


def _require_integer(value: object, description: str, *, minimum: int = 0) -> int:
    if not isinstance(value, Integral) or isinstance(value, bool):
        raise ValueError(f"{description} must be an integer")
    number = int(value)
    if number < minimum:
        raise ValueError(f"{description} is out of range")
    return number


def _normalize_pdb_id(value: object) -> str:
    if not isinstance(value, str):
        raise ValueError("PDB ID must be a string")
    normalized = value.strip().lower()
    if len(normalized) != 4 or not normalized.isascii() or not normalized.isalnum():
        raise ValueError("PDB ID must be exactly four ASCII alphanumeric characters")
    return normalized


def _validate_chain_id(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value.strip() != value
        or "\0" in value
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise ValueError("entry/chain ID must be nonempty canonical text without NUL")
    return value


def _true_count_bin(value: object) -> CountBin:
    count = _require_integer(value, "true domain count", minimum=1)
    if count >= 5:
        return "5+"
    return str(count)  # type: ignore[return-value]


def _length_bin(value: object) -> LengthBin:
    length = _require_integer(value, "chain length", minimum=1)
    if length < 250:
        return "<250"
    if length < 350:
        return "250-349"
    if length < 450:
        return "350-449"
    return "450+"


def _component_id(chain_ids: Sequence[str]) -> str:
    canonical = sorted(chain_ids)
    return hashlib.sha256("\0".join(canonical).encode("utf-8")).hexdigest()


def _validate_entry(entry: CathEntry) -> tuple[str, str, tuple[str, ...], tuple[str, ...]]:
    if not isinstance(entry, CathEntry):
        raise ValueError("fold entries must be CathEntry records")
    chain_id = _validate_chain_id(entry.entry_id)
    pdb_id = _normalize_pdb_id(entry.pdb_id)
    _true_count_bin(entry.n_domains)
    _length_bin(entry.n_residues)
    labels = cath_family_labels(entry.chopping)
    combination = cath_family_combination(entry.chopping)
    if combination != tuple(sorted(labels)):
        raise ValueError("CATH family combination is not canonical")
    return chain_id, pdb_id, labels, combination


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, index: int) -> int:
        while self.parent[index] != index:
            self.parent[index] = self.parent[self.parent[index]]
            index = self.parent[index]
        return index

    def union(self, left: int, right: int) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if left_root < right_root:
            self.parent[right_root] = left_root
        else:
            self.parent[left_root] = right_root


def connected_components(entries: Sequence[CathEntry]) -> list[tuple[CathEntry, ...]]:
    """Union entries through exact PDB IDs or nonempty family multisets."""
    validated_entries: list[
        tuple[CathEntry, tuple[str, str, tuple[str, ...], tuple[str, ...]]]
    ] = []
    seen_ids: set[str] = set()
    for entry in entries:
        validated = _validate_entry(entry)
        if validated[0] in seen_ids:
            raise ValueError("duplicate canonical entry ID")
        seen_ids.add(validated[0])
        validated_entries.append((entry, validated))
    validated_entries.sort(key=lambda item: item[1][0])
    ordered = [item[0] for item in validated_entries]
    metadata = [item[1] for item in validated_entries]

    union = _UnionFind(len(ordered))
    by_pdb: dict[str, int] = {}
    by_combination: dict[tuple[str, ...], int] = {}
    for index, (_chain_id, pdb_id, _labels, combination) in enumerate(metadata):
        if pdb_id in by_pdb:
            union.union(index, by_pdb[pdb_id])
        else:
            by_pdb[pdb_id] = index
        if combination:
            if combination in by_combination:
                union.union(index, by_combination[combination])
            else:
                by_combination[combination] = index

    groups: dict[int, list[CathEntry]] = {}
    for index, entry in enumerate(ordered):
        groups.setdefault(union.find(index), []).append(entry)
    components = [tuple(sorted(group, key=lambda entry: entry.entry_id)) for group in groups.values()]
    return sorted(
        components,
        key=lambda component: _component_id([entry.entry_id for entry in component]),
    )


def _component_vector(component: Sequence[CathEntry]) -> tuple[int, ...]:
    counts = {name: 0 for name in COUNT_BINS}
    lengths = {name: 0 for name in LENGTH_BINS}
    for entry in component:
        counts[_true_count_bin(entry.n_domains)] += 1
        lengths[_length_bin(entry.n_residues)] += 1
    return (
        len(component),
        *(counts[name] for name in COUNT_BINS),
        *(lengths[name] for name in LENGTH_BINS),
    )


def assign_folds(
    entries: Sequence[CathEntry],
    n_folds: int = 5,
    seed: int = 37,
) -> tuple[FoldAssignment, ...]:
    """Assign whole structural components using fixed marginal-balance costs."""
    fold_count = _require_integer(n_folds, "n_folds", minimum=2)
    seed_value = _require_integer(seed, "seed", minimum=0)
    components = connected_components(entries)
    if len(components) < fold_count:
        raise ValueError("fewer structural components than folds")
    if sum(len(component) for component in components) < fold_count:
        raise ValueError("fewer assignments than folds")

    component_rows: list[tuple[tuple[CathEntry, ...], str, tuple[int, ...], int]] = []
    for component in components:
        identifier = _component_id([entry.entry_id for entry in component])
        vector = _component_vector(component)
        largest_stratum = max(vector[1:])
        component_rows.append((component, identifier, vector, largest_stratum))
    component_rows.sort(
        key=lambda item: (
            -item[2][0],
            -item[3],
            hashlib.sha256(f"{seed_value}:{item[1]}".encode("utf-8")).hexdigest(),
            item[1],
        )
    )

    totals = [sum(item[2][index] for item in component_rows) for index in range(10)]
    targets = [total / fold_count for total in totals]
    fold_vectors = [[0] * 10 for _ in range(fold_count)]
    component_folds: dict[str, int] = {}
    for component_index, (_component, identifier, vector, _largest) in enumerate(
        component_rows
    ):
        remaining_components = len(component_rows) - component_index
        empty_folds = [
            fold for fold in range(fold_count) if fold_vectors[fold][0] == 0
        ]
        if len(empty_folds) > remaining_components:
            raise ValueError("too few remaining components to populate every fold")
        eligible_folds = (
            empty_folds
            if empty_folds and len(empty_folds) == remaining_components
            else list(range(fold_count))
        )
        choices: list[tuple[float, int, int]] = []
        for fold in eligible_folds:
            cost = sum(
                (
                    (fold_vectors[fold][index] + vector[index] - targets[index])
                    / max(targets[index], 1.0)
                )
                ** 2
                for index in range(10)
            )
            if not math.isfinite(cost):
                raise ValueError("fold balancing produced a non-finite cost")
            choices.append((cost, fold_vectors[fold][0], fold))
        chosen_fold = min(choices)[2]
        component_folds[identifier] = chosen_fold
        fold_vectors[chosen_fold] = [
            current + addition
            for current, addition in zip(fold_vectors[chosen_fold], vector, strict=True)
        ]

    assignments: list[FoldAssignment] = []
    for component, identifier, _vector, _largest in component_rows:
        for entry in component:
            chain_id, pdb_id, labels, combination = _validate_entry(entry)
            assignments.append(
                FoldAssignment(
                    chain_id=chain_id,
                    pdb_id=pdb_id,
                    family_combination=combination,
                    family_labels=labels,
                    true_count_bin=_true_count_bin(entry.n_domains),
                    length_bin=_length_bin(entry.n_residues),
                    component_id=identifier,
                    fold=component_folds[identifier],
                )
            )
    result = tuple(sorted(assignments, key=lambda assignment: assignment.chain_id))
    validate_folds(result)
    if len({assignment.fold for assignment in result}) != fold_count:
        raise ValueError("fold assignment left an empty fold")
    return result


def _assignment_components(
    assignments: Sequence[FoldAssignment],
) -> list[list[FoldAssignment]]:
    union = _UnionFind(len(assignments))
    by_pdb: dict[str, int] = {}
    by_combination: dict[tuple[str, ...], int] = {}
    for index, assignment in enumerate(assignments):
        if assignment.pdb_id in by_pdb:
            union.union(index, by_pdb[assignment.pdb_id])
        else:
            by_pdb[assignment.pdb_id] = index
        if assignment.family_combination:
            if assignment.family_combination in by_combination:
                union.union(index, by_combination[assignment.family_combination])
            else:
                by_combination[assignment.family_combination] = index
    groups: dict[int, list[FoldAssignment]] = {}
    for index, assignment in enumerate(assignments):
        groups.setdefault(union.find(index), []).append(assignment)
    return list(groups.values())


def validate_folds(assignments: Sequence[FoldAssignment]) -> None:
    """Fail closed unless assignments are a complete canonical component partition."""
    if not assignments:
        raise ValueError("fold assignments are empty")
    seen_chains: set[str] = set()
    observed_folds: set[int] = set()
    for assignment in assignments:
        if not isinstance(assignment, FoldAssignment):
            raise ValueError("fold assignment has the wrong type")
        chain_id = _validate_chain_id(assignment.chain_id)
        if chain_id in seen_chains:
            raise ValueError("duplicate fold chain ID")
        seen_chains.add(chain_id)
        if _normalize_pdb_id(assignment.pdb_id) != assignment.pdb_id:
            raise ValueError("fold PDB ID is not canonical")
        if not isinstance(assignment.family_labels, tuple) or any(
            not isinstance(label, str)
            or not label
            or label.strip() != label
            or label == "999_999"
            for label in assignment.family_labels
        ):
            raise ValueError("fold family labels are invalid")
        if (
            not isinstance(assignment.family_combination, tuple)
            or assignment.family_combination != tuple(sorted(assignment.family_labels))
        ):
            raise ValueError("fold family combination is not the exact sorted multiset")
        if assignment.true_count_bin not in COUNT_BINS:
            raise ValueError("fold true-count bin is invalid")
        if assignment.length_bin not in LENGTH_BINS:
            raise ValueError("fold length bin is invalid")
        _require_hash(assignment.component_id, "component_id")
        fold = _require_integer(assignment.fold, "fold", minimum=0)
        observed_folds.add(fold)

    if len(observed_folds) < 2 or observed_folds != set(range(max(observed_folds) + 1)):
        raise ValueError("fold numbers are noncontiguous or contain an empty fold")
    for group in _assignment_components(assignments):
        expected_component = _component_id([assignment.chain_id for assignment in group])
        if {assignment.component_id for assignment in group} != {expected_component}:
            raise ValueError("component ID does not match transitive structural membership")
        if len({assignment.fold for assignment in group}) != 1:
            raise ValueError("a structural component crosses folds")


def _validation_fold(assignments: Sequence[FoldAssignment], validation_fold: int) -> int:
    validate_folds(assignments)
    fold = _require_integer(validation_fold, "validation_fold", minimum=0)
    if fold not in {assignment.fold for assignment in assignments}:
        raise ValueError("validation_fold is absent")
    return fold


def individual_label_seen(
    assignments: Sequence[FoldAssignment],
    validation_fold: int,
) -> dict[str, bool]:
    fold = _validation_fold(assignments, validation_fold)
    training_labels = {
        label
        for assignment in assignments
        if assignment.fold != fold
        for label in assignment.family_labels
    }
    validation_labels = sorted(
        {
            label
            for assignment in assignments
            if assignment.fold == fold
            for label in assignment.family_labels
        }
    )
    return {label: label in training_labels for label in validation_labels}


def chain_label_cohorts(
    assignments: Sequence[FoldAssignment],
    validation_fold: int,
) -> dict[str, LabelCohort]:
    fold = _validation_fold(assignments, validation_fold)
    training_labels = {
        label
        for assignment in assignments
        if assignment.fold != fold
        for label in assignment.family_labels
    }
    cohorts: dict[str, LabelCohort] = {}
    for assignment in sorted(assignments, key=lambda value: value.chain_id):
        if assignment.fold != fold:
            continue
        labels = set(assignment.family_labels)
        if not labels:
            cohort: LabelCohort = "unknown"
        elif labels.issubset(training_labels):
            cohort = "seen"
        else:
            cohort = "unseen"
        cohorts[assignment.chain_id] = cohort
    return cohorts


def _require_hash(value: object, description: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{description} must be a lowercase SHA-256")
    return value


def _assignment_payload(assignment: FoldAssignment) -> dict[str, object]:
    return {
        "chain_id": assignment.chain_id,
        "pdb_id": assignment.pdb_id,
        "family_combination": list(assignment.family_combination),
        "family_labels": list(assignment.family_labels),
        "true_count_bin": assignment.true_count_bin,
        "length_bin": assignment.length_bin,
        "component_id": assignment.component_id,
        "fold": assignment.fold,
    }


def _manifest_payload(manifest: FoldManifest) -> dict[str, object]:
    return {
        "schema_version": manifest.schema_version,
        "seed": manifest.seed,
        "n_folds": manifest.n_folds,
        "dataset_sha256": manifest.dataset_sha256,
        "corpus_manifest_sha256": manifest.corpus_manifest_sha256,
        "chains_sha256": manifest.chains_sha256,
        "assignments": [
            _assignment_payload(assignment)
            for assignment in sorted(manifest.assignments, key=lambda value: value.chain_id)
        ],
    }


def _canonical_fold_bytes(manifest: FoldManifest) -> bytes:
    return (
        json.dumps(
            _manifest_payload(manifest),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def _atomic_write(path: Path, data: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def write_fold_manifest(
    path: Path,
    assignments: Sequence[FoldAssignment],
    dataset_sha256: str,
    corpus_sha256: str,
    *,
    chains_sha256: str,
    seed: int = 37,
) -> str:
    validate_folds(assignments)
    seed_value = _require_integer(seed, "seed", minimum=0)
    manifest = FoldManifest(
        schema_version=FOLD_MANIFEST_SCHEMA_VERSION,
        seed=seed_value,
        n_folds=max(assignment.fold for assignment in assignments) + 1,
        dataset_sha256=_require_hash(dataset_sha256, "dataset_sha256"),
        corpus_manifest_sha256=_require_hash(corpus_sha256, "corpus_manifest_sha256"),
        chains_sha256=_require_hash(chains_sha256, "chains_sha256"),
        assignments=tuple(sorted(assignments, key=lambda value: value.chain_id)),
    )
    data = _canonical_fold_bytes(manifest)
    _atomic_write(Path(path), data)
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value}")


def _parse_json(data: bytes) -> Any:
    try:
        return json.loads(
            data.decode("utf-8"),
            object_pairs_hook=_unique_object,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("manifest is not valid UTF-8 JSON") from error


_ASSIGNMENT_KEYS = {
    "chain_id",
    "pdb_id",
    "family_combination",
    "family_labels",
    "true_count_bin",
    "length_bin",
    "component_id",
    "fold",
}
_FOLD_MANIFEST_KEYS = {
    "schema_version",
    "seed",
    "n_folds",
    "dataset_sha256",
    "corpus_manifest_sha256",
    "chains_sha256",
    "assignments",
}


def _assignment_from_payload(value: object) -> FoldAssignment:
    if not isinstance(value, dict) or set(value) != _ASSIGNMENT_KEYS:
        raise ValueError("fold assignment schema mismatch")
    combination = value["family_combination"]
    labels = value["family_labels"]
    if not isinstance(combination, list) or not all(isinstance(label, str) for label in combination):
        raise ValueError("family_combination must be a string array")
    if not isinstance(labels, list) or not all(isinstance(label, str) for label in labels):
        raise ValueError("family_labels must be a string array")
    return FoldAssignment(
        chain_id=value["chain_id"],
        pdb_id=value["pdb_id"],
        family_combination=tuple(combination),
        family_labels=tuple(labels),
        true_count_bin=value["true_count_bin"],
        length_bin=value["length_bin"],
        component_id=value["component_id"],
        fold=value["fold"],
    )


def load_fold_manifest(
    path: Path,
    *,
    expected_sha256: str | None = None,
    expected_dataset_sha256: str | None = None,
    expected_corpus_sha256: str | None = None,
    expected_chains_sha256: str | None = None,
) -> FoldManifest:
    data = Path(path).read_bytes()
    actual_hash = hashlib.sha256(data).hexdigest()
    if expected_sha256 is not None and actual_hash != _require_hash(
        expected_sha256, "expected_sha256"
    ):
        raise ValueError("fold manifest SHA-256 mismatch")
    payload = _parse_json(data)
    if not isinstance(payload, dict) or set(payload) != _FOLD_MANIFEST_KEYS:
        raise ValueError("fold manifest schema mismatch")
    schema_version = _require_integer(payload["schema_version"], "schema_version", minimum=1)
    if schema_version != FOLD_MANIFEST_SCHEMA_VERSION:
        raise ValueError("unsupported fold manifest schema version")
    seed = _require_integer(payload["seed"], "seed", minimum=0)
    n_folds = _require_integer(payload["n_folds"], "n_folds", minimum=2)
    assignment_values = payload["assignments"]
    if not isinstance(assignment_values, list):
        raise ValueError("fold manifest assignments must be an array")
    assignments = tuple(_assignment_from_payload(value) for value in assignment_values)
    validate_folds(assignments)
    if n_folds != max(assignment.fold for assignment in assignments) + 1:
        raise ValueError("fold manifest n_folds disagrees with assignments")
    manifest = FoldManifest(
        schema_version=schema_version,
        seed=seed,
        n_folds=n_folds,
        dataset_sha256=_require_hash(payload["dataset_sha256"], "dataset_sha256"),
        corpus_manifest_sha256=_require_hash(
            payload["corpus_manifest_sha256"], "corpus_manifest_sha256"
        ),
        chains_sha256=_require_hash(payload["chains_sha256"], "chains_sha256"),
        assignments=assignments,
    )
    for observed, expected, description in (
        (manifest.dataset_sha256, expected_dataset_sha256, "dataset SHA-256"),
        (manifest.corpus_manifest_sha256, expected_corpus_sha256, "corpus SHA-256"),
        (manifest.chains_sha256, expected_chains_sha256, "chains SHA-256"),
    ):
        if expected is not None and observed != _require_hash(expected, description):
            raise ValueError(f"fold manifest {description} mismatch")
    if data != _canonical_fold_bytes(manifest):
        raise ValueError("fold manifest is not canonically serialized or ordered")
    return manifest


_TASK8_MANIFEST_KEYS = {
    "dataset",
    "dataset_sha256",
    "binary_sha256",
    "git_commit",
    "dump_argv_normalized",
    "build_argv_normalized",
    "seed",
    "jobs",
    "threads",
    "timeout",
    "limit",
    "resume",
    "feature_schema_version",
    "feature_schema_hash",
    "tables",
    "accepted_unique_chains",
    "rejected_unique_chains",
    "candidate_rejection_counts_by_code",
}
_TABLE_FIELDS = {
    "chains": CHAIN_FIELDS,
    "counts": COUNT_FIELDS,
    "candidates": CANDIDATE_FIELDS,
    "rejections": NORMALIZED_REJECTION_FIELDS,
}


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _canonical_task8_manifest(payload: object) -> bytes:
    return (
        json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")


def _require_child_file(corpus_dir: Path, name: str) -> Path:
    root = Path(corpus_dir).resolve(strict=True)
    path = (Path(corpus_dir) / name).resolve(strict=True)
    if path.parent != root or not path.is_file():
        raise ValueError(f"corpus input {name} is not a direct file beneath corpus-dir")
    return path


def _read_canonical_table(path: Path, fields: Sequence[str]) -> tuple[bytes, list[dict[str, str]]]:
    data = path.read_bytes()
    if not data.endswith(b"\n") or b"\r" in data:
        raise ValueError(f"corpus table {path.name} has noncanonical newlines")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"corpus table {path.name} is not UTF-8") from error
    reader = csv.DictReader(StringIO(text, newline=""))
    if tuple(reader.fieldnames or ()) != tuple(fields):
        raise ValueError(f"corpus table {path.name} header mismatch")
    rows = list(reader)
    expected_fields = set(fields)
    if any(set(row) != expected_fields or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"corpus table {path.name} row schema mismatch")
    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    if output.getvalue().encode("utf-8") != data:
        raise ValueError(f"corpus table {path.name} is not canonically serialized")
    return data, rows


def _validate_task8_manifest(payload: object, dataset: str) -> dict[str, Any]:
    if not isinstance(payload, dict) or set(payload) != _TASK8_MANIFEST_KEYS:
        raise ValueError("Task 8 corpus manifest schema mismatch")
    if payload["dataset"] != dataset or dataset != "cath17287":
        raise ValueError("Task 8 corpus dataset mismatch")
    _require_hash(payload["dataset_sha256"], "dataset_sha256")
    _require_hash(payload["binary_sha256"], "binary_sha256")
    if (
        not isinstance(payload["git_commit"], str)
        or len(payload["git_commit"]) not in {40, 64}
        or any(character not in "0123456789abcdef" for character in payload["git_commit"])
    ):
        raise ValueError("Task 8 corpus git_commit is invalid")
    for field in ("dump_argv_normalized", "build_argv_normalized"):
        if not isinstance(payload[field], list) or not all(
            isinstance(token, str) for token in payload[field]
        ):
            raise ValueError(f"Task 8 corpus {field} is invalid")
    if payload["seed"] != 37:
        raise ValueError("Task 8 corpus seed mismatch")
    for field in ("jobs", "threads", "timeout"):
        _require_integer(payload[field], field, minimum=1)
    if payload["limit"] is not None:
        _require_integer(payload["limit"], "limit", minimum=1)
    if type(payload["resume"]) is not bool:
        raise ValueError("Task 8 corpus resume setting is invalid")
    if payload["feature_schema_version"] != SCHEMA_VERSION:
        raise ValueError("Task 8 corpus feature schema version mismatch")
    if payload["feature_schema_hash"] != feature_schema_hash():
        raise ValueError("Task 8 corpus feature schema hash mismatch")
    for field in ("accepted_unique_chains", "rejected_unique_chains"):
        _require_integer(payload[field], field, minimum=0)
    rejection_counts = payload["candidate_rejection_counts_by_code"]
    if not isinstance(rejection_counts, dict) or any(
        not isinstance(code, str)
        or not code
        or not isinstance(count, Integral)
        or isinstance(count, bool)
        or count < 0
        for code, count in rejection_counts.items()
    ):
        raise ValueError("Task 8 corpus candidate rejection counts are invalid")
    tables = payload["tables"]
    if not isinstance(tables, dict) or set(tables) != set(_TABLE_FIELDS):
        raise ValueError("Task 8 corpus table manifest mismatch")
    for name, descriptor in tables.items():
        if not isinstance(descriptor, dict) or set(descriptor) != {"sha256", "rows"}:
            raise ValueError(f"Task 8 corpus table descriptor {name} is invalid")
        _require_hash(descriptor["sha256"], f"{name} table SHA-256")
        _require_integer(descriptor["rows"], f"{name} table row count", minimum=0)
    return payload


def _load_accepted_entries(
    corpus_dir: Path,
    dataset: str,
) -> tuple[tuple[CathEntry, ...], str, str, str]:
    """Load exact accepted Task 8 chains and their independently hashed sources."""
    if dataset != "cath17287":
        raise ValueError("fold construction accepts only cath17287")
    manifest_path = _require_child_file(Path(corpus_dir), "corpus_manifest.json")
    manifest_data = manifest_path.read_bytes()
    manifest_payload = _parse_json(manifest_data)
    manifest = _validate_task8_manifest(manifest_payload, dataset)
    if manifest_data != _canonical_task8_manifest(manifest):
        raise ValueError("Task 8 corpus manifest is not canonically serialized")
    corpus_manifest_sha256 = _sha256(manifest_data)

    table_rows: dict[str, list[dict[str, str]]] = {}
    table_bytes: dict[str, bytes] = {}
    for name, fields in _TABLE_FIELDS.items():
        path = _require_child_file(Path(corpus_dir), f"{name}.csv")
        descriptor = manifest["tables"][name]
        raw_data = path.read_bytes()
        if _sha256(raw_data) != descriptor["sha256"]:
            raise ValueError(f"Task 8 corpus {name} table hash mismatch")
        data, rows = _read_canonical_table(path, fields)
        if len(rows) != descriptor["rows"]:
            raise ValueError(f"Task 8 corpus {name} table row-count mismatch")
        table_rows[name] = rows
        table_bytes[name] = data

    chains: dict[str, int] = {}
    for row in table_rows["chains"]:
        chain_id = _validate_chain_id(row["chain_id"])
        if chain_id in chains:
            raise ValueError("duplicate accepted chain ID")
        raw_count = row["n_true_domains"]
        try:
            count = int(raw_count)
        except ValueError as error:
            raise ValueError("accepted n_true_domains is not an ordinary integer") from error
        if count <= 0 or str(count) != raw_count:
            raise ValueError("accepted n_true_domains is not a positive ordinary integer")
        for field in GLOBAL_FEATURES:
            try:
                value = float(row[field])
            except ValueError as error:
                raise ValueError(f"accepted global feature {field} is malformed") from error
            if not math.isfinite(value) or format(value, ".17g") != row[field]:
                raise ValueError(f"accepted global feature {field} is not canonical finite text")
        chains[chain_id] = count
    if not chains or len(chains) != manifest["accepted_unique_chains"]:
        raise ValueError("accepted Task 8 chain population is empty or inconsistent")
    for table_name in ("counts", "candidates"):
        observed = {row["chain_id"] for row in table_rows[table_name]}
        if observed != set(chains):
            raise ValueError(f"Task 8 {table_name} chains disagree with accepted population")
    for row in table_rows["candidates"]:
        raw_true = row["n_true_domains"]
        try:
            candidate_true = int(raw_true)
        except ValueError as error:
            raise ValueError("candidate n_true_domains is malformed") from error
        if str(candidate_true) != raw_true or chains[row["chain_id"]] != candidate_true:
            raise ValueError("candidate n_true_domains disagrees with chains.csv")

    source_path = Path(dataset_path(dataset))
    dataset_data = source_path.read_bytes()
    dataset_sha256 = _sha256(dataset_data)
    if dataset_sha256 != manifest["dataset_sha256"]:
        raise ValueError("Task 8 corpus dataset hash mismatch")
    metadata = load_dataset(dataset)
    metadata_by_id: dict[str, CathEntry] = {}
    for entry in metadata:
        if entry.dataset != dataset:
            raise ValueError("source metadata dataset mismatch")
        entry_id = _validate_chain_id(entry.entry_id)
        if entry_id in metadata_by_id:
            raise ValueError("duplicate source metadata entry_id")
        metadata_by_id[entry_id] = entry

    accepted: list[CathEntry] = []
    for chain_id in sorted(chains):
        entry = metadata_by_id.get(chain_id)
        if entry is None:
            raise ValueError("accepted chain ID is missing from source metadata")
        if _require_integer(entry.n_domains, "metadata true count", minimum=1) != chains[chain_id]:
            raise ValueError("accepted n_true_domains disagrees with source metadata")
        _validate_entry(entry)
        accepted.append(entry)
    return (
        tuple(accepted),
        dataset_sha256,
        corpus_manifest_sha256,
        _sha256(table_bytes["chains"]),
    )
