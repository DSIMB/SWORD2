"""Canonical post-tooling runtime freeze for the factorized selector."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import stat
import subprocess
import unicodedata
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from benchmark.datasets import dataset_path, load_dataset
from benchmark.factorized_ranker.corpus import (
    CANDIDATE_FIELDS,
    CHAIN_FIELDS,
    COUNT_FIELDS,
    NORMALIZED_REJECTION_FIELDS,
)
from benchmark.factorized_ranker.folds import (
    _canonical_task8_manifest,
    _validate_task8_manifest,
    load_fold_manifest,
)
from benchmark.factorized_ranker.model_artifact import verify_top_level_manifest


RUNTIME_FREEZE_SCHEMA_VERSION = 2

BUILD_PROFILE = "release"
BUILD_COMMAND = ["cargo", "build", "--locked", "--release", "--bin", "sword2"]
CACHE_CONTRACT = {
    "default_promotion_compatible": False,
    "factorized_cache_hit_behavior": "whole_chain_legacy_fallback",
    "typed_context_reloaded": False,
}

RUNTIME_FREEZE_KEYS = {
    "schema_version",
    "model_manifest_sha256",
    "model_source_git_commit",
    "model_artifact_sha256s",
    "development_population",
    "runtime_source_git_commit",
    "runtime_input_sha256s",
    "runtime_input_tree_sha256",
    "evidence_tool_sha256s",
    "evidence_tool_tree_sha256",
    "cargo_lock_sha256",
    "rustc_version_verbose",
    "cargo_version",
    "target_triple",
    "build_profile",
    "build_command",
    "build_environment",
    "binary_name",
    "binary_size",
    "binary_sha256",
    "selector_cache_contract",
}

MODEL_ARTIFACT_FIELDS = (
    "corpus_manifest_sha256",
    "fold_manifest_sha256",
    "cv_report_sha256",
    "ablation_report_sha256",
    "oof_predictions_sha256",
    "count_model_sha256",
    "candidate_model_sha256",
    "golden_sha256",
    "standalone_baseline_sha256",
    "generated_rust_sha256",
)

EVIDENCE_TOOL_PATHS_V1 = (
    "benchmark/factorized_ranker/runtime_freeze.py",
    "benchmark/freeze_factorized_runtime.py",
    "benchmark/evaluate_factorized_acceptance.py",
    "benchmark/run_benchmark.py",
    "benchmark/score.py",
    "benchmark/runners/base.py",
    "benchmark/runners/sword2_rust.py",
    "benchmark/compare_sword2_experiments.py",
    "benchmark/datasets.py",
    "benchmark/metrics.py",
    "benchmark/numbering.py",
    "benchmark/stats.py",
    "benchmark/structures.py",
)
EVIDENCE_TOOL_PATHS_V2 = tuple(
    sorted(
        {
            *EVIDENCE_TOOL_PATHS_V1,
            "benchmark/factorized_ranker/eligibility.py",
            "benchmark/freeze_factorized_eligibility.py",
        }
    )
)
EVIDENCE_TOOL_PATHS = EVIDENCE_TOOL_PATHS_V2

_TABLE_FIELDS: dict[str, Sequence[str]] = {
    "chains": CHAIN_FIELDS,
    "counts": COUNT_FIELDS,
    "candidates": CANDIDATE_FIELDS,
    "rejections": NORMALIZED_REJECTION_FIELDS,
}


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value}")


def canonical_json_bytes(payload: object) -> bytes:
    try:
        return (
            json.dumps(
                payload,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=False,
                allow_nan=False,
            ).encode("utf-8")
            + b"\n"
        )
    except (TypeError, ValueError) as error:
        raise ValueError("payload is not canonical finite JSON data") from error


def load_canonical_json(path: Path) -> dict[str, object]:
    data = Path(path).read_bytes()
    try:
        payload = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"{path} is not valid canonical JSON") from error
    if not isinstance(payload, dict) or canonical_json_bytes(payload) != data:
        raise ValueError(f"{path} is not canonically serialized")
    return payload


def _canonical_id(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value.strip() != value
        or unicodedata.normalize("NFC", value) != value
        or "\0" in value
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise ValueError("identity must be nonempty canonical Unicode text")
    return value


def canonical_id_set_hash(ids: Iterable[str]) -> str:
    materialized = [_canonical_id(value) for value in ids]
    if len(materialized) != len(set(materialized)):
        raise ValueError("identity collection contains duplicates")
    data = (
        json.dumps(
            sorted(materialized),
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode("utf-8")
        + b"\n"
    )
    return hashlib.sha256(data).hexdigest()


def _stable_file_bytes(path: Path) -> bytes:
    path = Path(path)
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError(f"input is not a regular nonsymlink file: {path}")
    data = path.read_bytes()
    after = path.lstat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after or len(data) != before.st_size:
        raise ValueError(f"input changed while hashing: {path}")
    return data


def _stable_file_hash_and_size(path: Path) -> tuple[str, int]:
    path = Path(path)
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError(f"input is not a regular nonsymlink file: {path}")
    digest = hashlib.sha256()
    byte_count = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
            byte_count += len(chunk)
    after = path.lstat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after or byte_count != before.st_size:
        raise ValueError(f"input changed while hashing: {path}")
    return digest.hexdigest(), byte_count


def _stable_symlink_hash_and_size(path: Path) -> tuple[str, int]:
    before = path.lstat()
    if not stat.S_ISLNK(before.st_mode):
        raise ValueError(f"input is not a symlink: {path}")
    target = os.fsencode(os.readlink(path))
    after = path.lstat()
    identity_before = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
    identity_after = (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    if identity_before != identity_after:
        raise ValueError(f"input changed while hashing: {path}")
    return hashlib.sha256(b"symlink\0" + target).hexdigest(), len(target)


def sha256_file(path: Path) -> str:
    return _stable_file_hash_and_size(path)[0]


def _tree_hash(mapping: Mapping[str, str]) -> str:
    digest = hashlib.sha256()
    for name, value in sorted(mapping.items()):
        name_bytes = name.encode("utf-8")
        value_bytes = value.encode("ascii")
        digest.update(len(name_bytes).to_bytes(8, "big"))
        digest.update(name_bytes)
        digest.update(len(value_bytes).to_bytes(8, "big"))
        digest.update(value_bytes)
    return digest.hexdigest()


def _canonical_relative_path(value: str) -> str:
    if (
        not value
        or value.startswith("/")
        or ".." in Path(value).parts
        or unicodedata.normalize("NFC", value) != value
        or "\0" in value
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise ValueError("tree contains a noncanonical relative path")
    return value


def hash_directory_tree(path: Path) -> str:
    raw_root = Path(path)
    root_info = raw_root.lstat()
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        raise ValueError("tree input must be a regular directory")
    root = raw_root.resolve(strict=True)
    mapping: dict[str, str] = {}
    children = sorted(root.rglob("*"))
    snapshot: dict[str, tuple[int, int, int, int, int]] = {}
    for child in children:
        relative = _canonical_relative_path(child.relative_to(root).as_posix())
        info = child.lstat()
        if stat.S_ISLNK(info.st_mode):
            raise ValueError(f"tree contains a symlink: {relative}")
        if stat.S_ISDIR(info.st_mode):
            continue
        if not stat.S_ISREG(info.st_mode):
            raise ValueError(f"tree contains a nonregular file: {relative}")
        if relative in mapping:
            raise ValueError("tree contains duplicate normalized paths")
        mapping[relative] = sha256_file(child)
        snapshot[relative] = (
            info.st_mode,
            info.st_dev,
            info.st_ino,
            info.st_size,
            info.st_mtime_ns,
        )
    if not mapping:
        raise ValueError("tree input is empty")
    after: dict[str, tuple[int, int, int, int, int]] = {}
    for child in sorted(root.rglob("*")):
        info = child.lstat()
        relative = _canonical_relative_path(child.relative_to(root).as_posix())
        if stat.S_ISDIR(info.st_mode):
            continue
        after[relative] = (
            info.st_mode,
            info.st_dev,
            info.st_ino,
            info.st_size,
            info.st_mtime_ns,
        )
    if after != snapshot:
        raise ValueError("tree changed while hashing")
    return _tree_hash(mapping)


def hash_file_or_tree(path: Path) -> dict[str, object]:
    path = Path(path)
    info = path.lstat()
    if stat.S_ISLNK(info.st_mode):
        digest, byte_count = _stable_symlink_hash_and_size(path)
        return {
            "kind": "file",
            "byte_count": byte_count,
            "file_count": 1,
            "sha256": digest,
        }
    if stat.S_ISREG(info.st_mode):
        digest, byte_count = _stable_file_hash_and_size(path)
        return {
            "kind": "file",
            "byte_count": byte_count,
            "file_count": 1,
            "sha256": digest,
        }
    if not stat.S_ISDIR(info.st_mode):
        raise ValueError(f"artifact is not a regular file or directory: {path}")
    root = path.resolve(strict=True)
    mapping: dict[str, str] = {}
    byte_count = 0
    snapshot: dict[str, tuple[int, int, int, int, int]] = {}
    for child in sorted(root.rglob("*")):
        child_info = child.lstat()
        relative = _canonical_relative_path(child.relative_to(root).as_posix())
        if stat.S_ISLNK(child_info.st_mode):
            digest, size = _stable_symlink_hash_and_size(child)
            mapping[relative] = digest
            byte_count += size
            snapshot[relative] = (
                child_info.st_mode,
                child_info.st_dev,
                child_info.st_ino,
                child_info.st_size,
                child_info.st_mtime_ns,
            )
            continue
        if stat.S_ISDIR(child_info.st_mode):
            continue
        if not stat.S_ISREG(child_info.st_mode):
            raise ValueError(f"artifact tree contains a nonregular file: {relative}")
        digest, size = _stable_file_hash_and_size(child)
        mapping[relative] = digest
        byte_count += size
        snapshot[relative] = (
            child_info.st_mode,
            child_info.st_dev,
            child_info.st_ino,
            child_info.st_size,
            child_info.st_mtime_ns,
        )
    if not mapping:
        raise ValueError("artifact tree is empty")
    after: dict[str, tuple[int, int, int, int, int]] = {}
    for child in sorted(root.rglob("*")):
        child_info = child.lstat()
        relative = _canonical_relative_path(child.relative_to(root).as_posix())
        if stat.S_ISDIR(child_info.st_mode):
            continue
        after[relative] = (
            child_info.st_mode,
            child_info.st_dev,
            child_info.st_ino,
            child_info.st_size,
            child_info.st_mtime_ns,
        )
    if after != snapshot:
        raise ValueError("artifact tree changed while hashing")
    return {
        "kind": "tree",
        "byte_count": byte_count,
        "file_count": len(mapping),
        "sha256": _tree_hash(mapping),
    }


def _require_hash(value: object, description: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{description} must be a lowercase SHA-256")
    return value


def _require_plain_int(value: object, description: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{description} must be an integer >= {minimum}")
    return value


def _git(repo_root: Path, *args: str) -> bytes:
    completed = subprocess.run(
        ["git", "-C", os.fspath(repo_root), *args],
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise ValueError(
            f"git {' '.join(args)} failed: {completed.stderr.decode('utf-8', 'replace').strip()}"
        )
    return completed.stdout


def _git_head(repo_root: Path) -> str:
    value = _git(repo_root, "rev-parse", "HEAD").decode("ascii").strip()
    if len(value) != 40 or any(character not in "0123456789abcdef" for character in value):
        raise ValueError("repository HEAD is not a full SHA-1 commit")
    return value


def _runtime_paths(repo_root: Path) -> tuple[str, ...]:
    tracked = _git(repo_root, "ls-files", "-z").decode("utf-8").split("\0")
    result: list[str] = []
    exact = {
        "Cargo.toml",
        "Cargo.lock",
        "sword2-cli/Cargo.toml",
        "sword2-lib/Cargo.toml",
        ".cargo/config.toml",
        ".cargo/config",
        "rust-toolchain",
        "rust-toolchain.toml",
        "sword2-cli/build.rs",
        "sword2-lib/build.rs",
    }
    for path in tracked:
        if not path:
            continue
        if path in exact or (
            (path.startswith("sword2-cli/src/") or path.startswith("sword2-lib/src/"))
            and path.endswith(".rs")
        ):
            result.append(path)
    required = {"Cargo.toml", "Cargo.lock", "sword2-cli/Cargo.toml", "sword2-lib/Cargo.toml"}
    if not required.issubset(result):
        raise ValueError("tracked runtime input set is incomplete")
    return tuple(sorted(result))


def _require_clean_paths(repo_root: Path, paths: Sequence[str], *, runtime_roots: bool) -> None:
    if not paths:
        raise ValueError("Git input path set is empty")
    status = _git(
        repo_root,
        "status",
        "--porcelain=v1",
        "--untracked-files=all",
        "--",
        *paths,
    )
    if status:
        raise ValueError("tracked source/tool inputs are dirty")
    if runtime_roots:
        untracked = _git(
            repo_root,
            "ls-files",
            "--others",
            "--exclude-standard",
            "-z",
            "--",
            "sword2-cli/src",
            "sword2-lib/src",
        )
        if untracked:
            raise ValueError("runtime source roots contain an untracked file")


def _hash_git_inputs(
    repo_root: Path,
    paths: Sequence[str],
    *,
    runtime_roots: bool,
) -> tuple[dict[str, str], str]:
    _require_clean_paths(repo_root, paths, runtime_roots=runtime_roots)
    mapping: dict[str, str] = {}
    for relative in paths:
        path = repo_root / relative
        if not path.exists() or path.is_symlink() or not path.is_file():
            raise ValueError(f"Git input is missing or nonregular: {relative}")
        mapping[relative] = sha256_file(path)
    return mapping, _tree_hash(mapping)


def _canonical_csv_ids(
    path: Path,
    expected_fields: Sequence[str],
    *,
    unique: bool,
) -> tuple[set[str], int, list[dict[str, str]] | None]:
    before = Path(path).lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError(f"{path.name} is not a regular nonsymlink file")
    ends_with_lf = False
    with Path(path).open("rb") as raw_handle:
        while chunk := raw_handle.read(1024 * 1024):
            if b"\r" in chunk:
                raise ValueError(f"{path.name} has noncanonical newlines")
            ends_with_lf = chunk.endswith(b"\n")
    if before.st_size == 0 or not ends_with_lf:
        raise ValueError(f"{path.name} has noncanonical newlines")
    ids: set[str] = set()
    rows = 0
    retained_rows: list[dict[str, str]] | None = [] if path.name == "rejections.csv" else None
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != tuple(expected_fields):
            raise ValueError(f"{path.name} header mismatch")
        for row in reader:
            if set(row) != set(expected_fields) or any(value is None for value in row.values()):
                raise ValueError(f"{path.name} row schema mismatch")
            entry_id = _canonical_id(row["chain_id"])
            if unique and entry_id in ids:
                raise ValueError(f"{path.name} contains duplicate accepted identities")
            ids.add(entry_id)
            rows += 1
            if retained_rows is not None:
                retained_rows.append(row)
    after = Path(path).lstat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise ValueError(f"{path.name} changed while parsing")
    return ids, rows, retained_rows


def _load_oof_ids(path: Path) -> set[str]:
    before = Path(path).lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError("OOF predictions are not a regular nonsymlink file")
    ends_with_lf = False
    with Path(path).open("rb") as raw_handle:
        while chunk := raw_handle.read(1024 * 1024):
            if b"\r" in chunk:
                raise ValueError("OOF predictions have noncanonical newlines")
            ends_with_lf = chunk.endswith(b"\n")
    if before.st_size == 0 or not ends_with_lf:
        raise ValueError("OOF predictions have noncanonical newlines")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "chain_id" not in reader.fieldnames:
            raise ValueError("OOF predictions lack chain_id")
        ids: set[str] = set()
        for row in reader:
            entry_id = _canonical_id(row.get("chain_id"))
            if entry_id in ids:
                raise ValueError("OOF predictions contain duplicate chain IDs")
            ids.add(entry_id)
    if not ids:
        raise ValueError("OOF predictions are empty")
    after = Path(path).lstat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    ):
        raise ValueError("OOF predictions changed while parsing")
    return ids


def _verify_development_population(
    *,
    model_manifest: Mapping[str, object],
    model_manifest_path: Path,
    corpus_dir: Path,
    fold_manifest_path: Path,
    oof_predictions: Path,
) -> dict[str, object]:
    corpus_info = Path(corpus_dir).lstat()
    if stat.S_ISLNK(corpus_info.st_mode) or not stat.S_ISDIR(corpus_info.st_mode):
        raise ValueError("development corpus is not a regular nonsymlink directory")
    corpus_root = Path(corpus_dir).resolve(strict=True)
    corpus_manifest_path = corpus_root / "corpus_manifest.json"
    corpus_manifest_bytes = _stable_file_bytes(corpus_manifest_path)
    try:
        corpus_payload = json.loads(
            corpus_manifest_bytes.decode("utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("corpus manifest is invalid JSON") from error
    corpus_manifest = _validate_task8_manifest(corpus_payload, "cath17287")
    if _canonical_task8_manifest(corpus_manifest) != corpus_manifest_bytes:
        raise ValueError("corpus manifest is not canonical")
    corpus_manifest_hash = hashlib.sha256(corpus_manifest_bytes).hexdigest()
    if corpus_manifest_hash != model_manifest["corpus_manifest_sha256"]:
        raise ValueError("source corpus manifest differs from frozen model graph")
    copied = model_manifest_path.parent / "cath17287_factorized_corpus_v1_manifest.json"
    if _stable_file_bytes(copied) != corpus_manifest_bytes:
        raise ValueError("source and copied corpus manifests differ")

    table_sets: dict[str, set[str]] = {}
    rejection_rows: list[dict[str, str]] = []
    for name, fields in _TABLE_FIELDS.items():
        path = corpus_root / f"{name}.csv"
        descriptor = corpus_manifest["tables"][name]
        if sha256_file(path) != descriptor["sha256"]:
            raise ValueError(f"corpus {name} hash mismatch")
        if descriptor["sha256"] != model_manifest["corpus_table_sha256s"][name]:
            raise ValueError(f"corpus {name} differs from model graph")
        ids, row_count, retained = _canonical_csv_ids(
            path,
            fields,
            unique=name == "chains",
        )
        if sha256_file(path) != descriptor["sha256"]:
            raise ValueError(f"corpus {name} changed while validating")
        if row_count != descriptor["rows"]:
            raise ValueError(f"corpus {name} row count mismatch")
        table_sets[name] = ids
        if retained is not None:
            rejection_rows = retained

    accepted = table_sets["chains"]
    rejected = table_sets["rejections"]
    if not accepted or len(accepted) != corpus_manifest["accepted_unique_chains"]:
        raise ValueError("accepted population count mismatch")
    if len(rejected) != corpus_manifest["rejected_unique_chains"]:
        raise ValueError("rejected population count mismatch")
    if accepted & rejected:
        raise ValueError("accepted and rejected populations overlap")
    if table_sets["counts"] != accepted or table_sets["candidates"] != accepted:
        raise ValueError("corpus tables disagree on accepted population")
    observed_candidate_codes = Counter(
        row["code"] for row in rejection_rows if row["scope"] == "candidate"
    )
    if dict(sorted(observed_candidate_codes.items())) != corpus_manifest[
        "candidate_rejection_counts_by_code"
    ]:
        raise ValueError("candidate rejection accounting mismatch")

    dataset_file = Path(dataset_path("cath17287"))
    if sha256_file(dataset_file) != corpus_manifest["dataset_sha256"]:
        raise ValueError("development dataset hash mismatch")
    metadata_ids = [_canonical_id(entry.entry_id) for entry in load_dataset("cath17287")]
    if sha256_file(dataset_file) != corpus_manifest["dataset_sha256"]:
        raise ValueError("development dataset changed while parsing")
    if len(metadata_ids) != len(set(metadata_ids)):
        raise ValueError("development metadata contains duplicate IDs")
    if set(metadata_ids) != accepted | rejected:
        raise ValueError("accepted/rejected accounting does not cover development metadata")

    fold_hash = sha256_file(fold_manifest_path)
    if fold_hash != model_manifest["fold_manifest_sha256"]:
        raise ValueError("fold manifest differs from model graph")
    fold_manifest = load_fold_manifest(
        fold_manifest_path,
        expected_sha256=fold_hash,
        expected_dataset_sha256=str(model_manifest["dataset_sha256"]),
        expected_corpus_sha256=corpus_manifest_hash,
        expected_chains_sha256=str(model_manifest["corpus_table_sha256s"]["chains"]),
    )
    if sha256_file(fold_manifest_path) != fold_hash:
        raise ValueError("fold manifest changed while validating")
    fold_ids = [_canonical_id(assignment.chain_id) for assignment in fold_manifest.assignments]
    if len(fold_ids) != len(set(fold_ids)) or set(fold_ids) != accepted:
        raise ValueError("fold IDs do not equal accepted population")
    if fold_manifest.n_folds != 5 or fold_manifest.seed != 37:
        raise ValueError("fold count or seed mismatch")
    if {assignment.fold for assignment in fold_manifest.assignments} != set(range(5)):
        raise ValueError("one or more frozen folds are empty")

    oof_hash = sha256_file(oof_predictions)
    if oof_hash != model_manifest["oof_predictions_sha256"]:
        raise ValueError("OOF predictions differ from model graph")
    oof_ids = _load_oof_ids(oof_predictions)
    if sha256_file(oof_predictions) != oof_hash:
        raise ValueError("OOF predictions changed while validating")
    if oof_ids != accepted:
        raise ValueError("OOF IDs do not equal accepted population")

    return {
        "accepted_count": len(accepted),
        "rejected_count": len(rejected),
        "accepted_id_set_sha256": canonical_id_set_hash(accepted),
        "fold_id_set_sha256": canonical_id_set_hash(fold_ids),
        "oof_id_set_sha256": canonical_id_set_hash(oof_ids),
        "rejected_id_set_sha256": canonical_id_set_hash(rejected),
        "dataset_sha256": _require_hash(model_manifest["dataset_sha256"], "dataset hash"),
        "corpus_manifest_sha256": corpus_manifest_hash,
        "corpus_table_sha256s": {
            name: descriptor["sha256"]
            for name, descriptor in sorted(corpus_manifest["tables"].items())
        },
        "fold_manifest_sha256": fold_hash,
        "oof_predictions_sha256": oof_hash,
        "n_folds": 5,
        "seed": 37,
    }


def _command_text(command: Sequence[str], cwd: Path) -> str:
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise ValueError(f"{' '.join(command)} failed")
    return completed.stdout.strip()


def _build_payload(
    *,
    model_manifest_path: Path,
    corpus_dir: Path,
    fold_manifest: Path,
    oof_predictions: Path,
    binary: Path,
    repo_root: Path,
    runtime_source_git_commit: str,
    schema_version: int = RUNTIME_FREEZE_SCHEMA_VERSION,
) -> dict[str, object]:
    repo_root = Path(repo_root).resolve(strict=True)
    model_info = Path(model_manifest_path).lstat()
    if stat.S_ISLNK(model_info.st_mode) or not stat.S_ISREG(model_info.st_mode):
        raise ValueError("model manifest is not a regular nonsymlink file")
    model_manifest_path = Path(model_manifest_path).resolve(strict=True)
    model_manifest_hash = sha256_file(model_manifest_path)
    model_manifest = verify_top_level_manifest(model_manifest_path)
    if sha256_file(model_manifest_path) != model_manifest_hash:
        raise ValueError("model manifest changed while validating")
    development = _verify_development_population(
        model_manifest=model_manifest,
        model_manifest_path=model_manifest_path,
        corpus_dir=corpus_dir,
        fold_manifest_path=fold_manifest,
        oof_predictions=oof_predictions,
    )

    runtime_paths = _runtime_paths(repo_root)
    runtime_hashes, runtime_tree_hash = _hash_git_inputs(
        repo_root,
        runtime_paths,
        runtime_roots=True,
    )
    evidence_tool_paths = _evidence_tool_paths(schema_version)
    evidence_hashes, evidence_tree_hash = _hash_git_inputs(
        repo_root,
        evidence_tool_paths,
        runtime_roots=False,
    )
    binary_info = Path(binary).lstat()
    if stat.S_ISLNK(binary_info.st_mode) or not stat.S_ISREG(binary_info.st_mode):
        raise ValueError("runtime binary is not a regular nonsymlink file")
    binary_path = Path(binary).resolve(strict=True)
    binary_hash, binary_size = _stable_file_hash_and_size(binary_path)
    if binary_path.name != "sword2" or not os.access(binary_path, os.X_OK):
        raise ValueError("runtime binary must be an executable named sword2")

    rustc_version = _command_text(["rustc", "-Vv"], repo_root)
    cargo_version = _command_text(["cargo", "-V"], repo_root)
    hosts = [line.split(":", 1)[1].strip() for line in rustc_version.splitlines() if line.startswith("host:")]
    if len(hosts) != 1 or not hosts[0]:
        raise ValueError("rustc -Vv did not expose one host target")

    artifact_hashes = {
        field: _require_hash(model_manifest[field], field)
        for field in MODEL_ARTIFACT_FIELDS
    }
    payload: dict[str, object] = {
        "schema_version": schema_version,
        "model_manifest_sha256": model_manifest_hash,
        "model_source_git_commit": model_manifest["source_git_commit"],
        "model_artifact_sha256s": artifact_hashes,
        "development_population": development,
        "runtime_source_git_commit": runtime_source_git_commit,
        "runtime_input_sha256s": runtime_hashes,
        "runtime_input_tree_sha256": runtime_tree_hash,
        "evidence_tool_sha256s": evidence_hashes,
        "evidence_tool_tree_sha256": evidence_tree_hash,
        "cargo_lock_sha256": runtime_hashes["Cargo.lock"],
        "rustc_version_verbose": rustc_version,
        "cargo_version": cargo_version,
        "target_triple": hosts[0],
        "build_profile": BUILD_PROFILE,
        "build_command": BUILD_COMMAND,
        "build_environment": {
            name: os.environ.get(name)
            for name in ("CARGO_INCREMENTAL", "RUSTFLAGS", "RUSTC_WRAPPER")
        },
        "binary_name": "sword2",
        "binary_size": binary_size,
        "binary_sha256": binary_hash,
        "selector_cache_contract": CACHE_CONTRACT,
    }
    current_runtime_hashes, current_runtime_tree = _hash_git_inputs(
        repo_root,
        runtime_paths,
        runtime_roots=True,
    )
    current_evidence_hashes, current_evidence_tree = _hash_git_inputs(
        repo_root,
        evidence_tool_paths,
        runtime_roots=False,
    )
    if (
        current_runtime_hashes != runtime_hashes
        or current_runtime_tree != runtime_tree_hash
        or current_evidence_hashes != evidence_hashes
        or current_evidence_tree != evidence_tree_hash
    ):
        raise ValueError("runtime or evidence-tool sources changed while freezing")
    if _stable_file_hash_and_size(binary_path) != (binary_hash, binary_size):
        raise ValueError("runtime binary changed while freezing")
    if sha256_file(model_manifest_path) != model_manifest_hash:
        raise ValueError("model manifest changed while freezing")
    _validate_runtime_payload(payload)
    return payload


def create_runtime_freeze(
    *,
    model_manifest: Path,
    corpus_dir: Path,
    fold_manifest: Path,
    oof_predictions: Path,
    binary: Path,
    repo_root: Path,
) -> dict[str, object]:
    root = Path(repo_root).resolve(strict=True)
    return _build_payload(
        model_manifest_path=model_manifest,
        corpus_dir=corpus_dir,
        fold_manifest=fold_manifest,
        oof_predictions=oof_predictions,
        binary=binary,
        repo_root=root,
        runtime_source_git_commit=_git_head(root),
    )


def _evidence_tool_paths(schema_version: int) -> tuple[str, ...]:
    if schema_version == 1:
        return EVIDENCE_TOOL_PATHS_V1
    if schema_version == 2:
        return EVIDENCE_TOOL_PATHS_V2
    raise ValueError("unsupported runtime freeze schema version")


def _validate_runtime_payload(payload: Mapping[str, object]) -> None:
    if set(payload) != RUNTIME_FREEZE_KEYS:
        raise ValueError("runtime freeze top-level schema mismatch")
    schema_version = _require_plain_int(payload["schema_version"], "schema version", 1)
    evidence_tool_paths = _evidence_tool_paths(schema_version)
    for field in (
        "model_manifest_sha256",
        "runtime_input_tree_sha256",
        "evidence_tool_tree_sha256",
        "cargo_lock_sha256",
        "binary_sha256",
    ):
        _require_hash(payload[field], field)
    for field in ("model_source_git_commit", "runtime_source_git_commit"):
        value = payload[field]
        if not isinstance(value, str) or len(value) != 40 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError(f"{field} must be a full lowercase Git SHA-1")
    for field in ("model_artifact_sha256s", "runtime_input_sha256s", "evidence_tool_sha256s"):
        value = payload[field]
        if not isinstance(value, dict) or not value:
            raise ValueError(f"{field} must be a nonempty mapping")
        if any(not isinstance(key, str) or not key or not isinstance(item, str) for key, item in value.items()):
            raise ValueError(f"{field} has invalid entries")
        for item in value.values():
            _require_hash(item, field)
    if payload["runtime_input_tree_sha256"] != _tree_hash(payload["runtime_input_sha256s"]):
        raise ValueError("runtime input tree hash disagrees with its mapping")
    if payload["evidence_tool_tree_sha256"] != _tree_hash(payload["evidence_tool_sha256s"]):
        raise ValueError("evidence tool tree hash disagrees with its mapping")
    if (
        payload["runtime_input_sha256s"].get("Cargo.lock")
        != payload["cargo_lock_sha256"]
    ):
        raise ValueError("Cargo.lock hash disagrees with the runtime input mapping")
    if set(payload["model_artifact_sha256s"]) != set(MODEL_ARTIFACT_FIELDS):
        raise ValueError("runtime model-artifact mapping mismatch")
    for field in ("runtime_input_sha256s", "evidence_tool_sha256s"):
        for name in payload[field]:
            if _canonical_relative_path(name) != name or Path(name).as_posix() != name:
                raise ValueError(f"{field} contains a noncanonical path")
    if set(payload["evidence_tool_sha256s"]) != set(evidence_tool_paths):
        raise ValueError("runtime evidence-tool path set mismatch")
    development = payload["development_population"]
    expected_development = {
        "accepted_count",
        "rejected_count",
        "accepted_id_set_sha256",
        "fold_id_set_sha256",
        "oof_id_set_sha256",
        "rejected_id_set_sha256",
        "dataset_sha256",
        "corpus_manifest_sha256",
        "corpus_table_sha256s",
        "fold_manifest_sha256",
        "oof_predictions_sha256",
        "n_folds",
        "seed",
    }
    if not isinstance(development, dict) or set(development) != expected_development:
        raise ValueError("runtime development-population schema mismatch")
    _require_plain_int(development["accepted_count"], "accepted count", 1)
    _require_plain_int(development["rejected_count"], "rejected count", 0)
    if development["n_folds"] != 5 or development["seed"] != 37:
        raise ValueError("runtime fold contract mismatch")
    for field in expected_development - {"accepted_count", "rejected_count", "n_folds", "seed", "corpus_table_sha256s"}:
        _require_hash(development[field], field)
    tables = development["corpus_table_sha256s"]
    if not isinstance(tables, dict) or set(tables) != set(_TABLE_FIELDS):
        raise ValueError("runtime corpus-table mapping mismatch")
    for value in tables.values():
        _require_hash(value, "corpus table hash")
    if not (
        development["accepted_id_set_sha256"]
        == development["fold_id_set_sha256"]
        == development["oof_id_set_sha256"]
    ):
        raise ValueError("runtime accepted/fold/OOF population hashes differ")
    if payload["build_profile"] != BUILD_PROFILE or payload["build_command"] != BUILD_COMMAND:
        raise ValueError("runtime build contract mismatch")
    if payload["binary_name"] != "sword2":
        raise ValueError("runtime binary name mismatch")
    _require_plain_int(payload["binary_size"], "binary size", 1)
    environment = payload["build_environment"]
    if not isinstance(environment, dict) or set(environment) != {
        "CARGO_INCREMENTAL",
        "RUSTFLAGS",
        "RUSTC_WRAPPER",
    } or any(value is not None and not isinstance(value, str) for value in environment.values()):
        raise ValueError("runtime build environment mismatch")
    for field in ("rustc_version_verbose", "cargo_version", "target_triple"):
        if not isinstance(payload[field], str) or not payload[field]:
            raise ValueError(f"runtime {field} is invalid")
    cache_contract = payload["selector_cache_contract"]
    if (
        not isinstance(cache_contract, dict)
        or set(cache_contract) != set(CACHE_CONTRACT)
        or cache_contract["default_promotion_compatible"] is not False
        or cache_contract["factorized_cache_hit_behavior"]
        != "whole_chain_legacy_fallback"
        or cache_contract["typed_context_reloaded"] is not False
    ):
        raise ValueError("runtime selector cache contract is not schema-v1 compatible")


def write_runtime_freeze(path: Path, payload: Mapping[str, object]) -> str:
    _validate_runtime_payload(payload)
    target = Path(path)
    if target.exists() or target.is_symlink():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    data = canonical_json_bytes(dict(payload))
    temporary = target.with_name(f".{target.name}.tmp-{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(temporary)
    try:
        with temporary.open("xb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, target)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    temporary.unlink()
    if target.read_bytes() != data:
        raise OSError("installed runtime freeze bytes differ")
    return hashlib.sha256(data).hexdigest()


def _require_ancestor(repo_root: Path, commit: str) -> None:
    exists = subprocess.run(
        ["git", "-C", os.fspath(repo_root), "cat-file", "-e", f"{commit}^{{commit}}"],
        check=False,
        capture_output=True,
    )
    if exists.returncode != 0:
        raise ValueError("runtime source commit does not exist")
    ancestor = subprocess.run(
        ["git", "-C", os.fspath(repo_root), "merge-base", "--is-ancestor", commit, "HEAD"],
        check=False,
        capture_output=True,
    )
    if ancestor.returncode != 0:
        raise ValueError("runtime source commit is not an ancestor of HEAD")


def verify_runtime_freeze(
    path: Path,
    *,
    binary: Path,
    repo_root: Path,
) -> dict[str, object]:
    manifest_info = Path(path).lstat()
    if stat.S_ISLNK(manifest_info.st_mode) or not stat.S_ISREG(manifest_info.st_mode):
        raise ValueError("runtime manifest is not a regular nonsymlink file")
    manifest_path = Path(path).resolve(strict=True)
    payload = load_canonical_json(manifest_path)
    _validate_runtime_payload(payload)
    root = Path(repo_root).resolve(strict=True)
    runtime_commit = str(payload["runtime_source_git_commit"])
    _require_ancestor(root, runtime_commit)
    expected = _build_payload(
        model_manifest_path=root / "benchmark/models/factorized_ranker_v1_manifest.json",
        corpus_dir=root / "benchmark/data/cath17287_factorized_corpus_v1",
        fold_manifest=root / "benchmark/models/cath17287_factorized_folds_v1.json",
        oof_predictions=root / "benchmark/models/factorized_ranker_v1_oof.csv",
        binary=binary,
        repo_root=root,
        runtime_source_git_commit=runtime_commit,
        schema_version=int(payload["schema_version"]),
    )
    if payload != expected:
        raise ValueError("runtime freeze no longer matches model/source/tool/binary inputs")
    return payload


__all__ = [
    "BUILD_COMMAND",
    "CACHE_CONTRACT",
    "EVIDENCE_TOOL_PATHS",
    "EVIDENCE_TOOL_PATHS_V1",
    "EVIDENCE_TOOL_PATHS_V2",
    "RUNTIME_FREEZE_SCHEMA_VERSION",
    "canonical_id_set_hash",
    "canonical_json_bytes",
    "create_runtime_freeze",
    "hash_directory_tree",
    "hash_file_or_tree",
    "load_canonical_json",
    "sha256_file",
    "verify_runtime_freeze",
    "write_runtime_freeze",
]
