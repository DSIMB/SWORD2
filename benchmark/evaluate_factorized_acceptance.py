"""Validate locked evidence and evaluate the frozen factorized-ranker gates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import stat
import subprocess
import sys
import unicodedata
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark.factorized_ranker.runtime_freeze import (
    CACHE_CONTRACT,
    canonical_id_set_hash,
    canonical_json_bytes,
    load_canonical_json,
    sha256_file,
    verify_runtime_freeze,
)
from benchmark.factorized_ranker.eligibility import (
    POLICY as FACTORIZED_ELIGIBILITY_POLICY,
    load_eligibility_manifest,
)
from benchmark.factorized_ranker.model_artifact import verify_top_level_manifest
from benchmark.datasets import CathEntry, parse_cath_domain_string, read_merizo_csv
from benchmark.score import FailureRow, LockedRunRow
from benchmark.stats import paired_chain_bootstrap


ACCEPTANCE_SCHEMA_VERSION = 2
COVERAGE_SCHEMA_VERSION = 2
LOCKED_POPULATION_SIZE = 663
BOOTSTRAP = {
    "method": "paired_chain_mean_delta",
    "quantile_method": "linear",
    "replicates": 10_000,
    "seed": 37,
}
GATE_ORDER = (
    "overall_ndo",
    "merizo_ndo_ci_low",
    "chainsaw_ndo_ci_low",
    "domain_count_accuracy",
    "boundary_f1_10",
    "contiguous_standalone_ndo_delta",
    "discontinuous_standalone_ndo_delta",
    "runtime_median_ratio",
    "rss_maxima_ratio",
)


class InvalidEvidence(ValueError):
    """A stable exit-2 validation error for locked benchmark evidence."""


LOCKED_MANIFEST_KEYS = {
    "schema_version",
    "status",
    "created_at",
    "completed_at",
    "hostname",
    "platform",
    "locked_role",
    "dataset",
    "dataset_sha256",
    "dataset_id_count",
    "dataset_ids",
    "dataset_id_set_sha256",
    "structure_sha256s",
    "structure_tree_sha256",
    "raw_evidence_sha256s",
    "raw_evidence_tree_sha256",
    "tools",
    "sword2_extra_args",
    "sword2_threads",
    "environment",
    "normalized_argv",
    "gnu_time_version",
    "external_artifacts",
    "expected_chainsaw",
    "order_design",
    "model_manifest_sha256",
    "runtime_manifest_sha256",
    "eligibility_manifest_sha256",
    "eligibility_policy",
    "factorized_eligible_count",
    "factorized_eligible_id_set_sha256",
    "structural_abstention_count",
    "structural_abstention_id_set_sha256",
    "locked_jobs",
    "binary_sha256",
    "runtime_source_git_commit",
    "runtime_input_tree_sha256",
    "evidence_tool_tree_sha256",
    "cargo_lock_sha256",
    "row_counts",
    "success_id_counts",
    "success_id_set_sha256s",
    "failure_id_counts",
    "failure_id_set_sha256s",
    "selector_counts",
    "fallback_count",
    "excluded_candidate_count",
    "scores_sha256",
    "runs_sha256",
    "failures_sha256",
}

SCORE_IDENTITY_COLUMNS = (
    "dataset",
    "entry_id",
    "pdb_id",
    "chain_id",
    "tool",
    "variant",
    "partition",
)

COVERAGE_KEYS = {
    "schema_version",
    "valid",
    "roles",
    "model_manifest_sha256",
    "runtime_manifest_sha256",
    "eligibility_manifest_sha256",
    "eligibility_policy",
    "binary_sha256",
    "dataset_sha256",
    "dataset_id_count",
    "dataset_id_set_sha256",
    "factorized_eligible_count",
    "factorized_eligible_id_set_sha256",
    "structural_abstention_count",
    "structural_abstention_id_set_sha256",
    "factorized_structural_coverage",
    "chainsaw_success_count",
    "chainsaw_success_id_set_sha256",
    "chainsaw_failure_count",
    "chainsaw_failure_id_set_sha256",
    "resource_pair_count",
    "resource_order_sha256",
    "selector_counts",
    "fallback_count",
    "excluded_candidate_count",
    "evidence_sha256s",
}

ACCEPTANCE_KEYS = {
    "schema_version",
    "gate_order",
    "bootstrap",
    "model_manifest_sha256",
    "runtime_manifest_sha256",
    "coverage_attestation_sha256",
    "evidence_sha256s",
    "denominators",
    "gates",
    "diagnostics",
    "all_gates_measured",
    "all_gates_pass",
}

LOCKED_RUN_HEADER = tuple(LockedRunRow.__dataclass_fields__)
FAILURE_HEADER = tuple(FailureRow.__dataclass_fields__)
EVIDENCE_ARGUMENT_NAMES = (
    "legacy_manifest",
    "legacy_scores",
    "legacy_runs",
    "legacy_failures",
    "factorized_manifest",
    "factorized_scores",
    "factorized_runs",
    "factorized_failures",
    "resource_manifest",
    "resource_runs",
    "resource_failures",
    "dataset_metadata",
    "chainsaw_expected_success_ids",
    "model_manifest",
    "runtime_manifest",
    "eligibility_manifest",
    "standalone_baseline",
)


def _hash_value(value: object, description: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise InvalidEvidence(f"{description} is not a lowercase SHA-256")
    return value


def _ordinary_int(value: object, description: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise InvalidEvidence(f"{description} is not an integer >= {minimum}")
    return value


def _length_framed_mapping_hash(mapping: Mapping[str, str]) -> str:
    digest = hashlib.sha256()
    for name, value in sorted(mapping.items()):
        name_bytes = name.encode("utf-8")
        value_bytes = value.encode("ascii")
        digest.update(len(name_bytes).to_bytes(8, "big"))
        digest.update(name_bytes)
        digest.update(len(value_bytes).to_bytes(8, "big"))
        digest.update(value_bytes)
    return digest.hexdigest()


def _require_regular_file(path: Path, description: str) -> Path:
    path = Path(path)
    try:
        info = path.lstat()
    except OSError as error:
        raise InvalidEvidence(f"{description} is unavailable") from error
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
        raise InvalidEvidence(f"{description} is not a regular nonsymlink file")
    return path.resolve(strict=True)


def _validate_input_paths(paths: Mapping[str, Path]) -> dict[str, Path]:
    resolved: dict[str, Path] = {}
    inodes: dict[tuple[int, int], str] = {}
    for role, path in paths.items():
        canonical = _require_regular_file(path, role)
        info = canonical.stat()
        identity = (info.st_dev, info.st_ino)
        if canonical in resolved.values():
            raise InvalidEvidence(f"evidence input paths alias at {role}")
        if identity in inodes:
            raise InvalidEvidence(
                f"evidence inputs {inodes[identity]} and {role} are hard-link aliases"
            )
        resolved[role] = canonical
        inodes[identity] = role
    return resolved


def _validate_absent_outputs(outputs: Mapping[str, Path], inputs: Mapping[str, Path]) -> None:
    resolved_outputs: dict[str, Path] = {}
    for role, raw in outputs.items():
        path = Path(raw).absolute()
        if path.exists() or path.is_symlink():
            raise InvalidEvidence(f"{role} must be absent")
        try:
            parent = path.parent.resolve(strict=True)
        except OSError as error:
            raise InvalidEvidence(f"{role} parent is unavailable") from error
        canonical = parent / path.name
        if canonical in resolved_outputs.values():
            raise InvalidEvidence("acceptance outputs alias")
        if canonical in inputs.values():
            raise InvalidEvidence(f"{role} aliases an evidence input")
        resolved_outputs[role] = canonical


def _parse_bool_cell(value: str, description: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise InvalidEvidence(f"{description} is not an exact CSV boolean")


def _canonical_ids(values: object, description: str) -> tuple[str, ...]:
    if not isinstance(values, list):
        raise InvalidEvidence(f"{description} is not an array")
    result: list[str] = []
    for value in values:
        if (
            not isinstance(value, str)
            or not value
            or value.strip() != value
            or unicodedata.normalize("NFC", value) != value
            or "\0" in value
            or any(ord(character) < 32 or ord(character) == 127 for character in value)
        ):
            raise InvalidEvidence(f"{description} contains a noncanonical ID")
        result.append(value)
    if result != sorted(result) or len(result) != len(set(result)):
        raise InvalidEvidence(f"{description} is not sorted and unique")
    return tuple(result)


def _load_expected_ids(path: Path) -> tuple[str, ...]:
    if Path(path).is_symlink() or not Path(path).is_file():
        raise InvalidEvidence("expected-success ID input is not a regular file")
    data = Path(path).read_bytes()
    if not data or not data.endswith(b"\n") or b"\r" in data:
        raise InvalidEvidence("expected-success ID input has noncanonical newlines")
    try:
        values = data.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise InvalidEvidence("expected-success ID input is not UTF-8") from error
    return _canonical_ids(values, "expected-success IDs")


def _read_locked_metadata(path: Path) -> list[CathEntry]:
    try:
        entries = read_merizo_csv(path, dataset="cath663")
    except (OSError, ValueError) as error:
        raise InvalidEvidence("locked dataset metadata is invalid") from error
    raw_ids: list[str] = []
    try:
        with Path(path).open(newline="", encoding="utf-8") as handle:
            for row in csv.reader(handle):
                if not row or row[0].startswith("#"):
                    continue
                if len(row) < 7:
                    raise InvalidEvidence("locked dataset metadata row is truncated")
                raw_ids.append(row[1])
    except (OSError, UnicodeError, csv.Error) as error:
        raise InvalidEvidence("locked dataset metadata bytes are invalid") from error
    if len(raw_ids) != len(set(raw_ids)):
        raise InvalidEvidence("dataset metadata IDs contain duplicates")
    for value in raw_ids:
        if (
            not value
            or value.strip() != value
            or unicodedata.normalize("NFC", value) != value
            or "\0" in value
            or any(ord(character) < 32 or ord(character) == 127 for character in value)
        ):
            raise InvalidEvidence("dataset metadata contains a noncanonical ID")
    if [entry.entry_id for entry in entries] != raw_ids:
        raise InvalidEvidence("dataset reader changed locked entry identities")
    return entries


def _read_dict_csv(
    path: Path,
    expected_header: tuple[str, ...] | None = None,
) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    if Path(path).is_symlink() or not Path(path).is_file():
        raise InvalidEvidence(f"evidence CSV is not a regular file: {path}")
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        header = tuple(reader.fieldnames or ())
        if not header or len(header) != len(set(header)):
            raise InvalidEvidence(f"evidence CSV header is empty or duplicated: {path}")
        if expected_header is not None and header != expected_header:
            raise InvalidEvidence(f"evidence CSV header mismatch: {path}")
        rows = list(reader)
    if any(set(row) != set(header) or any(value is None for value in row.values()) for row in rows):
        raise InvalidEvidence(f"evidence CSV row schema mismatch: {path}")
    return header, rows


def _load_locked_manifest(path: Path, role: str) -> dict[str, object]:
    try:
        manifest = load_canonical_json(path)
    except (OSError, ValueError) as error:
        raise InvalidEvidence(f"locked {role} manifest is invalid") from error
    if set(manifest) != LOCKED_MANIFEST_KEYS:
        raise InvalidEvidence(f"locked {role} manifest schema mismatch")
    if (
        type(manifest["schema_version"]) is not int
        or manifest["schema_version"] != 2
        or manifest["status"] != "complete"
    ):
        raise InvalidEvidence(f"locked {role} manifest is incomplete")
    if manifest["locked_role"] != role:
        raise InvalidEvidence(f"locked manifest role mismatch for {role}")
    ids = _canonical_ids(manifest["dataset_ids"], f"{role} dataset IDs")
    if manifest["dataset_id_count"] != len(ids):
        raise InvalidEvidence(f"locked {role} dataset count mismatch")
    if manifest["dataset_id_set_sha256"] != canonical_id_set_hash(ids):
        raise InvalidEvidence(f"locked {role} dataset ID hash mismatch")
    if manifest["sword2_threads"] != 1:
        raise InvalidEvidence(f"locked {role} thread count mismatch")
    expected_jobs = 1 if role == "paired-sword-resources" else 32
    if manifest["locked_jobs"] != expected_jobs:
        raise InvalidEvidence(f"locked {role} worker count mismatch")
    if manifest["eligibility_policy"] != FACTORIZED_ELIGIBILITY_POLICY:
        raise InvalidEvidence(f"locked {role} eligibility policy mismatch")
    if (
        _ordinary_int(
            manifest["factorized_eligible_count"],
            f"locked {role} eligible count",
        )
        + _ordinary_int(
            manifest["structural_abstention_count"],
            f"locked {role} structural-abstention count",
        )
        != len(ids)
    ):
        raise InvalidEvidence(f"locked {role} eligibility counts mismatch")
    if manifest["dataset"] != "cath663":
        raise InvalidEvidence(f"locked {role} dataset name mismatch")
    expected_tools = {
        "legacy-accuracy": ["sword2-rust", "merizo", "chainsaw"],
        "factorized-accuracy": ["sword2-rust"],
        "paired-sword-resources": ["sword2-rust"],
    }[role]
    if manifest["tools"] != expected_tools:
        raise InvalidEvidence(f"locked {role} tool contract mismatch")
    expected_extra = ["--use-factorized-ranker"] if role == "factorized-accuracy" else []
    if manifest["sword2_extra_args"] != expected_extra:
        raise InvalidEvidence(f"locked {role} SWORD arguments mismatch")
    environment = manifest["environment"]
    expected_environment = {
        "PYTHONHASHSEED": "0",
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "VECLIB_MAXIMUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
        "RAYON_NUM_THREADS": "1",
    }
    if environment != expected_environment:
        raise InvalidEvidence(f"locked {role} environment mismatch")
    if not isinstance(manifest["gnu_time_version"], str) or "gnu time" not in manifest[
        "gnu_time_version"
    ].casefold():
        raise InvalidEvidence(f"locked {role} GNU-time identity mismatch")
    for field in (
        "dataset_sha256",
        "dataset_id_set_sha256",
        "structure_tree_sha256",
        "raw_evidence_tree_sha256",
        "model_manifest_sha256",
        "runtime_manifest_sha256",
        "eligibility_manifest_sha256",
        "factorized_eligible_id_set_sha256",
        "structural_abstention_id_set_sha256",
        "binary_sha256",
        "runtime_input_tree_sha256",
        "evidence_tool_tree_sha256",
        "cargo_lock_sha256",
        "runs_sha256",
        "failures_sha256",
    ):
        _hash_value(manifest[field], f"{role} {field}")
    if manifest["scores_sha256"] is not None:
        _hash_value(manifest["scores_sha256"], f"{role} scores hash")
    for field in ("created_at", "completed_at", "hostname", "platform"):
        if not isinstance(manifest[field], str) or not manifest[field]:
            raise InvalidEvidence(f"locked {role} {field} is invalid")
    if not isinstance(manifest["normalized_argv"], list) or not all(
        isinstance(value, str) for value in manifest["normalized_argv"]
    ):
        raise InvalidEvidence(f"locked {role} normalized argv is invalid")
    normalized_argv = manifest["normalized_argv"]
    if any(
        not token
        or "\0" in token
        or token.startswith("/")
        or "=/" in token
        for token in normalized_argv
    ):
        raise InvalidEvidence(f"locked {role} normalized argv leaks a physical path")

    def option_value(name: str) -> str | None:
        exact = [
            normalized_argv[index + 1]
            for index, token in enumerate(normalized_argv[:-1])
            if token == name
        ]
        assigned = [
            token.split("=", 1)[1]
            for token in normalized_argv
            if token.startswith(name + "=")
        ]
        values = exact + assigned
        return values[0] if len(values) == 1 else None

    if (
        option_value("--locked-role") != role
        or option_value("--dataset") != "cath663"
        or option_value("--cache-dir") != "<cache>"
        or option_value("--results-dir") != "<results>"
        or option_value("--locked-factorized-manifest") != "<runtime_manifest>"
        or option_value("--locked-eligibility-manifest")
        != "<eligibility_manifest>"
        or option_value("--locked-jobs") != str(expected_jobs)
        or option_value("--sword2-threads") != "1"
        or option_value("--tools") != ",".join(expected_tools)
        or normalized_argv.count("--strict") != 1
        or normalized_argv.count("--no-download") != 1
    ):
        raise InvalidEvidence(f"locked {role} normalized invocation mismatch")
    forbidden = {
        "--allow-dl-cpu",
        "--limit",
        "--skip-existing",
        "--reuse-tool-results-dir",
        "--score-csv",
        "--sword2-experiments",
    }
    if any(
        token in forbidden
        or any(token.startswith(name + "=") for name in forbidden)
        for token in normalized_argv
    ):
        raise InvalidEvidence(f"locked {role} invocation contains a forbidden option")
    merizo_device = option_value("--merizo-device")
    if merizo_device not in {None, "cuda"}:
        raise InvalidEvidence(f"locked {role} invocation changes the Merizo device")
    extra_value = option_value("--sword2-extra-args")
    expected_extra_value = "--use-factorized-ranker" if role == "factorized-accuracy" else None
    if extra_value != expected_extra_value:
        raise InvalidEvidence(f"locked {role} normalized SWORD-extra option mismatch")
    expected_success_value = option_value("--chainsaw-expected-success-ids")
    if role == "legacy-accuracy":
        if expected_success_value != "<chainsaw_expected_success_ids>":
            raise InvalidEvidence("locked legacy invocation lacks frozen Chainsaw IDs")
    elif expected_success_value is not None:
        raise InvalidEvidence(f"locked {role} invocation unexpectedly binds Chainsaw IDs")

    artifact_tokens = [
        normalized_argv[index + 1]
        for index, token in enumerate(normalized_argv[:-1])
        if token == "--locked-artifact"
    ]
    expected_artifact_tokens = sorted(
        f"{artifact_role}=<locked_artifact:{artifact_role}>"
        for artifact_role in manifest["external_artifacts"]
    )
    if sorted(artifact_tokens) != expected_artifact_tokens:
        raise InvalidEvidence(f"locked {role} normalized artifact arguments mismatch")
    for field in (
        "row_counts",
        "success_id_counts",
        "success_id_set_sha256s",
        "failure_id_counts",
        "failure_id_set_sha256s",
        "selector_counts",
    ):
        if not isinstance(manifest[field], dict):
            raise InvalidEvidence(f"locked {role} {field} is invalid")
    if set(manifest["row_counts"]) != {"scores", "runs", "failures"}:
        raise InvalidEvidence(f"locked {role} row-count schema mismatch")
    for name, value in manifest["row_counts"].items():
        _ordinary_int(value, f"locked {role} {name} row count")
    if set(manifest["selector_counts"]) != {"legacy", "factorized"}:
        raise InvalidEvidence(f"locked {role} selector-count schema mismatch")
    for name, value in manifest["selector_counts"].items():
        _ordinary_int(value, f"locked {role} {name} selector count")
    _ordinary_int(manifest["fallback_count"], f"locked {role} fallback count")
    _ordinary_int(
        manifest["excluded_candidate_count"],
        f"locked {role} excluded-candidate count",
    )
    structure_hashes = manifest["structure_sha256s"]
    raw_hashes = manifest["raw_evidence_sha256s"]
    if not isinstance(structure_hashes, dict) or set(structure_hashes) != set(ids):
        raise InvalidEvidence(f"locked {role} structure mapping mismatch")
    if not isinstance(raw_hashes, dict) or not raw_hashes:
        raise InvalidEvidence(f"locked {role} raw evidence mapping is empty")
    for mapping, description in ((structure_hashes, "structure"), (raw_hashes, "raw evidence")):
        for key, value in mapping.items():
            if not isinstance(key, str) or not key:
                raise InvalidEvidence(f"locked {role} {description} role is invalid")
            _hash_value(value, f"locked {role} {description} hash")
    if manifest["structure_tree_sha256"] != _length_framed_mapping_hash(structure_hashes):
        raise InvalidEvidence(f"locked {role} structure-tree hash mismatch")
    if manifest["raw_evidence_tree_sha256"] != _length_framed_mapping_hash(raw_hashes):
        raise InvalidEvidence(f"locked {role} raw-evidence tree hash mismatch")
    return manifest


def _validate_manifest_files(
    manifest_path: Path,
    manifest: Mapping[str, object],
    *,
    scores_path: Path | None,
    runs_path: Path,
    failures_path: Path,
) -> None:
    root = Path(manifest_path).parent.resolve(strict=True)
    expected_paths = [runs_path, failures_path]
    if scores_path is not None:
        expected_paths.append(scores_path)
    if any(Path(path).parent.resolve(strict=True) != root for path in expected_paths):
        raise InvalidEvidence("locked evidence files do not share their manifest root")
    if sha256_file(runs_path) != manifest["runs_sha256"]:
        raise InvalidEvidence("locked runs bytes differ from manifest")
    if sha256_file(failures_path) != manifest["failures_sha256"]:
        raise InvalidEvidence("locked failures bytes differ from manifest")
    if scores_path is None:
        if manifest["scores_sha256"] is not None:
            raise InvalidEvidence("resource manifest unexpectedly binds scores")
    elif sha256_file(scores_path) != manifest["scores_sha256"]:
        raise InvalidEvidence("locked scores bytes differ from manifest")
    raw_hashes = manifest["raw_evidence_sha256s"]
    observed: dict[str, str] = {}
    observed_inodes: dict[tuple[int, int], str] = {}
    authority_inodes = {
        (Path(path).stat().st_dev, Path(path).stat().st_ino)
        for path in expected_paths
    }
    manifest_info = Path(manifest_path).stat()
    authority_inodes.add((manifest_info.st_dev, manifest_info.st_ino))
    raw_root = root / "raw"
    if not raw_root.is_dir() or raw_root.is_symlink():
        raise InvalidEvidence("locked raw evidence root is unavailable")
    for path in sorted(raw_root.rglob("*")):
        if path.is_symlink():
            raise InvalidEvidence("locked raw evidence contains a symlink")
        if path.is_dir():
            continue
        if not path.is_file():
            raise InvalidEvidence("locked raw evidence contains a nonregular file")
        info = path.stat()
        identity = (info.st_dev, info.st_ino)
        relative = path.relative_to(root).as_posix()
        if identity in authority_inodes or identity in observed_inodes:
            raise InvalidEvidence("locked raw evidence contains an aliased file")
        observed_inodes[identity] = relative
        observed[relative] = sha256_file(path)
    if observed != raw_hashes:
        raise InvalidEvidence("locked raw evidence bytes differ from manifest")


def _score_identity_rows(path: Path) -> list[dict[str, str]]:
    header, rows = _read_dict_csv(path)
    if not set(SCORE_IDENTITY_COLUMNS).issubset(header):
        raise InvalidEvidence("locked score identity columns are missing")
    return [
        {column: row[column] for column in SCORE_IDENTITY_COLUMNS}
        for row in rows
    ]


def _unique_exact_rows(
    rows: list[dict[str, str]],
    predicate,
    expected_ids: set[str],
    description: str,
) -> None:
    selected = [row for row in rows if predicate(row)]
    ids = [row["entry_id"] for row in selected]
    if len(ids) != len(set(ids)) or set(ids) != expected_ids:
        raise InvalidEvidence(f"{description} coverage is duplicated or incomplete")


def _parse_json_cell(value: str, description: str) -> list[str]:
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise InvalidEvidence(f"{description} is invalid JSON") from error
    if not isinstance(parsed, list) or not all(isinstance(token, str) for token in parsed):
        raise InvalidEvidence(f"{description} is not a string array")
    if json.dumps(parsed, separators=(",", ":"), ensure_ascii=False) != value:
        raise InvalidEvidence(f"{description} is not canonical JSON")
    return parsed


def _validate_sword_run(
    row: dict[str, str],
    selector: str,
    manifest: Mapping[str, object],
    result_root: Path,
    *,
    resource: bool,
    structural_abstention: bool = False,
) -> None:
    if row["tool"] != "sword2-rust" or row["selector_variant"] != selector:
        raise InvalidEvidence("locked SWORD selector row mismatch")
    if row["locked_role"] != manifest["locked_role"]:
        raise InvalidEvidence("locked SWORD row role mismatch")
    if _parse_bool_cell(row["reused_output"], "locked reused_output") or row["returncode"] != "0":
        raise InvalidEvidence("locked SWORD row is reused or unsuccessful")
    try:
        runtime = float(row["runtime_s"])
        rss = int(row["peak_rss_kb"])
        exclusions = int(row["excluded_candidate_count"])
    except ValueError as error:
        raise InvalidEvidence("locked SWORD runtime/RSS/status values are malformed") from error
    if not math.isfinite(runtime) or runtime <= 0 or rss <= 0 or exclusions < 0:
        raise InvalidEvidence("locked SWORD runtime/RSS/status values are invalid")
    fallback = _parse_bool_cell(row["fallback"], "locked fallback")
    if structural_abstention:
        if not (
            selector == "factorized"
            and row["requested_selector"] == "factorized"
            and row["selector_used"] == "legacy"
            and fallback is True
            and row["error_code"] == "structural_quality_abstention"
            and row["selector_warning_code"] == "factorized_fallback"
            and exclusions == 0
        ):
            raise InvalidEvidence(
                "locked structural abstention status is inconsistent"
            )
    elif (
        row["requested_selector"] != selector
        or row["selector_used"] != selector
        or fallback
        or row["error_code"]
        or row["selector_warning_code"]
    ):
        raise InvalidEvidence("locked SWORD status does not prove requested inference")
    if selector == "legacy" and exclusions != 0:
        raise InvalidEvidence("locked legacy row has candidate exclusions")
    if (
        row["binary_sha256"] != manifest["binary_sha256"]
        or row["model_manifest_sha256"] != manifest["model_manifest_sha256"]
        or row["runtime_manifest_sha256"] != manifest["runtime_manifest_sha256"]
    ):
        raise InvalidEvidence("locked SWORD row identity hashes disagree with manifest")
    entry_id = row["entry_id"]
    if row["input_structure_sha256"] != manifest["structure_sha256s"].get(entry_id):
        raise InvalidEvidence("locked SWORD input structure hash mismatch")
    normalized = _parse_json_cell(row["normalized_command_json"], "normalized SWORD command")
    command = _parse_json_cell(row["command_json"], "exact SWORD command")
    expected = ["<binary>", "-i", "<input_structure>", "-o", "<process_output>", "-j", "1"]
    if selector == "factorized":
        expected.append("--use-factorized-ranker")
    if normalized != expected:
        raise InvalidEvidence("locked SWORD normalized command mismatch")
    if len(command) != len(normalized) or not command or Path(command[0]).name != "sword2":
        raise InvalidEvidence("locked SWORD exact command mismatch")
    if row["cwd_role"] != "repository_root":
        raise InvalidEvidence("locked SWORD cwd role mismatch")
    summary_role = row["raw_summary_role"]
    if not summary_role.startswith("raw/") or Path(summary_role).is_absolute() or ".." in Path(summary_role).parts:
        raise InvalidEvidence("locked SWORD summary role is invalid")
    summary = result_root / summary_role
    if (
        manifest["raw_evidence_sha256s"].get(summary_role) != row["raw_summary_sha256"]
        or not summary.is_file()
        or sha256_file(summary) != row["raw_summary_sha256"]
    ):
        raise InvalidEvidence("locked SWORD summary hash mismatch")
    process_candidates = [
        parent
        for parent in summary.parents
        if parent != result_root
        and result_root in parent.parents
        and (parent / "selector_status.json").is_file()
        and (parent / "stdout.log").is_file()
        and (parent / "stderr.log").is_file()
    ]
    if len(process_candidates) != 1:
        raise InvalidEvidence("locked SWORD process evidence root is ambiguous")
    process_dir = process_candidates[0]
    status_path = process_dir / "selector_status.json"
    stdout_path = process_dir / "stdout.log"
    stderr_path = process_dir / "stderr.log"
    for path, field in (
        (status_path, "selector_status_sha256"),
        (stdout_path, "stdout_sha256"),
        (stderr_path, "stderr_sha256"),
    ):
        if not path.is_file() or sha256_file(path) != row[field]:
            raise InvalidEvidence("locked SWORD raw status/log hash mismatch")
        relative = path.relative_to(result_root).as_posix()
        if manifest["raw_evidence_sha256s"].get(relative) != row[field]:
            raise InvalidEvidence("locked SWORD row hash is not bound by raw manifest")
    from benchmark.run_benchmark import parse_selector_status

    try:
        parsed_status = parse_selector_status(
            status_path.read_bytes(),
            selector,
            require_success=not structural_abstention,
        )
    except ValueError as error:
        raise InvalidEvidence("locked SWORD status bytes are invalid") from error
    parsed_error = (
        "" if parsed_status["error_code"] is None else str(parsed_status["error_code"])
    )
    if (
        int(parsed_status["excluded_candidate_count"]) != exclusions
        or str(parsed_status["requested_selector"])
        != row["requested_selector"]
        or str(parsed_status["selector_used"]) != row["selector_used"]
        or bool(parsed_status["fallback"]) is not fallback
        or parsed_error != row["error_code"]
    ):
        raise InvalidEvidence("locked SWORD exclusion count disagrees with status bytes")
    if resource:
        if row["pair_order"] not in {"legacy_first", "factorized_first"} or row["pair_position"] not in {"1", "2"}:
            raise InvalidEvidence("locked SWORD resource pair metadata is invalid")
    elif row["pair_order"] or row["pair_position"]:
        raise InvalidEvidence("accuracy SWORD row contains resource pair metadata")


def _validate_competitor_run(
    row: dict[str, str],
    tool: str,
    manifest: Mapping[str, object],
    result_root: Path,
) -> None:
    if (
        row["tool"] != tool
        or row["locked_role"] != "legacy-accuracy"
        or _parse_bool_cell(row["reused_output"], f"locked {tool} reused_output")
        or row["returncode"] != "0"
    ):
        raise InvalidEvidence(f"locked {tool} run row is invalid")
    for field in (
        "selector_variant",
        "pair_order",
        "pair_position",
        "selector_status_sha256",
        "requested_selector",
        "selector_used",
        "fallback",
        "error_code",
        "selector_warning_code",
        "excluded_candidate_count",
        "binary_sha256",
        "model_manifest_sha256",
        "runtime_manifest_sha256",
    ):
        if row[field]:
            raise InvalidEvidence(f"locked {tool} row contains SWORD-only evidence")
    if row["input_structure_sha256"] != manifest["structure_sha256s"].get(row["entry_id"]):
        raise InvalidEvidence(f"locked {tool} input structure hash mismatch")
    command = _parse_json_cell(row["command_json"], f"locked {tool} exact command")
    normalized = _parse_json_cell(row["normalized_command_json"], f"locked {tool} normalized command")
    if (
        not normalized
        or len(command) != len(normalized)
        or any(token.startswith("/") or "=/" in token for token in normalized)
    ):
        raise InvalidEvidence(f"locked {tool} normalized command leaks a physical path")
    if tool == "merizo":
        if (
            normalized[:3]
            != ["<locked_artifact:merizo_python>", "predict.py", "-i"]
            or "-d" not in normalized
            or normalized[normalized.index("-d") + 1 :]
            != [
                "cuda",
                "--return_indices",
                "--output_headers",
                "--pdb_chain",
                row["chain_id"],
            ]
            or f"<input_structure:{row['entry_id']}>" not in normalized
        ):
            raise InvalidEvidence("locked Merizo normalized command mismatch")
    elif normalized != [
        "<locked_artifact:chainsaw_python>",
        "get_predictions.py",
        "--structure_directory",
        "<results>/raw/chainsaw/_batch_stage",
        "--output",
        "<results>/raw/chainsaw/_batch.tsv",
    ]:
        raise InvalidEvidence("locked Chainsaw normalized command mismatch")
    if row["cwd_role"] != f"locked_artifact:{tool}_source_tree":
        raise InvalidEvidence(f"locked {tool} cwd role mismatch")
    if not row["raw_summary_role"].startswith("raw/"):
        raise InvalidEvidence(f"locked {tool} raw role is invalid")
    raw = result_root / row["raw_summary_role"]
    if (
        manifest["raw_evidence_sha256s"].get(row["raw_summary_role"])
        != row["raw_summary_sha256"]
        or not raw.is_file()
        or sha256_file(raw) != row["raw_summary_sha256"]
    ):
        raise InvalidEvidence(f"locked {tool} raw prediction hash mismatch")
    raw_values = set(manifest["raw_evidence_sha256s"].values())
    if row["stdout_sha256"] not in raw_values or row["stderr_sha256"] not in raw_values:
        raise InvalidEvidence(f"locked {tool} log hashes are not bound by raw manifest")


def _validate_role_runs(
    rows: list[dict[str, str]],
    manifest: Mapping[str, object],
    expected_ids: set[str],
    factorized_eligible_ids: set[str],
    structural_abstention_ids: set[str],
    chainsaw_ids: set[str],
    result_root: Path,
) -> tuple[int, int]:
    role = str(manifest["locked_role"])
    if len(rows) != manifest["row_counts"]["runs"]:
        raise InvalidEvidence(f"locked {role} run-row count mismatch")
    for row in rows:
        if row["dataset"] != manifest["dataset"] or row["entry_id"] not in expected_ids:
            raise InvalidEvidence(f"locked {role} run identity mismatch")
    run_keys = [
        (row["entry_id"], row["tool"], row["selector_variant"], row["pair_position"])
        for row in rows
    ]
    if len(run_keys) != len(set(run_keys)):
        raise InvalidEvidence(f"locked {role} contains duplicate run rows")
    fallback_count = 0
    exclusions = 0
    if role == "paired-sword-resources":
        resource_ids = factorized_eligible_ids
        if len(rows) != 2 * len(resource_ids):
            raise InvalidEvidence("paired resource rows are incomplete")
        by_id: dict[str, list[dict[str, str]]] = {}
        for row in rows:
            by_id.setdefault(row["entry_id"], []).append(row)
        if set(by_id) != resource_ids:
            raise InvalidEvidence("paired resource ID set mismatch")
        assignments = manifest["order_design"]
        if not isinstance(assignments, dict):
            raise InvalidEvidence("paired resource order design is absent")
        assignment_rows = assignments.get("assignments", [])
        if not isinstance(assignment_rows, list) or not all(
            isinstance(item, list)
            and len(item) == 2
            and all(isinstance(value, str) for value in item)
            for item in assignment_rows
        ):
            raise InvalidEvidence("paired resource assignments are malformed")
        frozen = {entry_id: first for entry_id, first in assignment_rows}
        if set(frozen) != resource_ids:
            raise InvalidEvidence("paired resource order assignment mismatch")
        ordered = sorted(
            resource_ids,
            key=lambda entry_id: (
                hashlib.sha256(b"37\0" + entry_id.encode("utf-8")).digest(),
                entry_id,
            ),
        )
        expected_assignment_rows = [
            [entry_id, "legacy" if index % 2 == 0 else "factorized"]
            for index, entry_id in enumerate(ordered)
        ]
        if assignment_rows != expected_assignment_rows:
            raise InvalidEvidence("paired resource order is not the frozen seed-37 design")
        if (
            assignments.get("seed") != 37
            or assignments.get("assignments_sha256")
            != hashlib.sha256(canonical_json_bytes(assignment_rows)).hexdigest()
            or assignments.get("legacy_first_count")
            != sum(first == "legacy" for first in frozen.values())
            or assignments.get("factorized_first_count")
            != sum(first == "factorized" for first in frozen.values())
        ):
            raise InvalidEvidence("paired resource order hash/count mismatch")
        for entry_id, pair in by_id.items():
            if len(pair) != 2 or {row["selector_variant"] for row in pair} != {"legacy", "factorized"}:
                raise InvalidEvidence("paired resource selector pair is incomplete")
            first = frozen[entry_id]
            expected_order = f"{first}_first"
            for row in pair:
                _validate_sword_run(
                    row,
                    row["selector_variant"],
                    manifest,
                    result_root,
                    resource=True,
                )
                if row["pair_order"] != expected_order:
                    raise InvalidEvidence("paired resource row disagrees with frozen order")
                expected_position = "1" if row["selector_variant"] == first else "2"
                if row["pair_position"] != expected_position:
                    raise InvalidEvidence("paired resource position disagrees with frozen order")
                exclusions += int(row["excluded_candidate_count"])
    else:
        if manifest["order_design"] is not None:
            raise InvalidEvidence(f"locked {role} unexpectedly has a resource order")
        selector = "legacy" if role == "legacy-accuracy" else "factorized"
        sword_rows = [row for row in rows if row["tool"] == "sword2-rust"]
        if len(sword_rows) != len(expected_ids) or {row["entry_id"] for row in sword_rows} != expected_ids:
            raise InvalidEvidence(f"locked {role} SWORD run coverage mismatch")
        for row in sword_rows:
            structural_abstention = (
                role == "factorized-accuracy"
                and row["entry_id"] in structural_abstention_ids
            )
            _validate_sword_run(
                row,
                selector,
                manifest,
                result_root,
                resource=False,
                structural_abstention=structural_abstention,
            )
            fallback_count += int(
                _parse_bool_cell(row["fallback"], "locked fallback")
            )
            exclusions += int(row["excluded_candidate_count"])
        if role == "factorized-accuracy":
            observed_abstentions = {
                row["entry_id"]
                for row in sword_rows
                if _parse_bool_cell(row["fallback"], "locked fallback")
            }
            if observed_abstentions != structural_abstention_ids:
                raise InvalidEvidence(
                    "factorized abstention IDs disagree with eligibility"
                )
        other = [row for row in rows if row["tool"] != "sword2-rust"]
        if role == "factorized-accuracy":
            if other:
                raise InvalidEvidence("factorized accuracy contains competitor run rows")
        else:
            merizo = [row for row in other if row["tool"] == "merizo"]
            chainsaw = [row for row in other if row["tool"] == "chainsaw"]
            if len(merizo) != len(expected_ids) or {row["entry_id"] for row in merizo} != expected_ids:
                raise InvalidEvidence("locked Merizo run coverage mismatch")
            if len(chainsaw) != len(chainsaw_ids) or {row["entry_id"] for row in chainsaw} != chainsaw_ids:
                raise InvalidEvidence("locked Chainsaw run coverage mismatch")
            for row in merizo:
                _validate_competitor_run(row, "merizo", manifest, result_root)
            for row in chainsaw:
                _validate_competitor_run(row, "chainsaw", manifest, result_root)
            if len(other) != len(merizo) + len(chainsaw):
                raise InvalidEvidence("legacy accuracy contains an unknown tool")
    if fallback_count != manifest["fallback_count"] or exclusions != manifest["excluded_candidate_count"]:
        raise InvalidEvidence(f"locked {role} selector accounting mismatch")
    observed_selector_counts = {
        "legacy": sum(row["selector_variant"] == "legacy" for row in rows),
        "factorized": sum(row["selector_variant"] == "factorized" for row in rows),
    }
    if observed_selector_counts != manifest["selector_counts"]:
        raise InvalidEvidence(f"locked {role} selector-count mismatch")
    return fallback_count, exclusions


def _validate_external_artifacts(manifest: Mapping[str, object]) -> None:
    role = str(manifest["locked_role"])
    expected = (
        {
            "merizo_python",
            "merizo_source_tree",
            "merizo_model_tree",
            "chainsaw_python",
            "chainsaw_source_tree",
            "chainsaw_model_tree",
        }
        if role == "legacy-accuracy"
        else set()
    )
    artifacts = manifest["external_artifacts"]
    if not isinstance(artifacts, dict) or set(artifacts) != expected:
        raise InvalidEvidence(f"locked {role} external-artifact roles mismatch")
    for artifact_role, descriptor in artifacts.items():
        if not isinstance(descriptor, dict) or set(descriptor) != {
            "role",
            "kind",
            "byte_count",
            "file_count",
            "sha256",
        }:
            raise InvalidEvidence(f"locked artifact {artifact_role} schema mismatch")
        if descriptor["role"] != artifact_role or descriptor["kind"] not in {"file", "tree"}:
            raise InvalidEvidence(f"locked artifact {artifact_role} identity mismatch")
        _ordinary_int(descriptor["byte_count"], f"{artifact_role} byte count", 1)
        _ordinary_int(descriptor["file_count"], f"{artifact_role} file count", 1)
        _hash_value(descriptor["sha256"], f"{artifact_role} hash")


def _validate_manifest_accounting(
    manifest: Mapping[str, object],
    expected_ids: set[str],
    factorized_eligible_ids: set[str],
    chainsaw_ids: set[str],
    score_rows: list[dict[str, str]] | None,
    failure_rows: list[dict[str, str]],
) -> None:
    role = str(manifest["locked_role"])
    expected_success: dict[str, set[str]]
    if role == "paired-sword-resources":
        expected_success = {
            "sword2-rust:factorized": factorized_eligible_ids,
            "sword2-rust:legacy": factorized_eligible_ids,
        }
    elif role == "factorized-accuracy":
        expected_success = {"sword2-rust": factorized_eligible_ids}
    else:
        expected_success = {
            "sword2-rust": expected_ids,
            "merizo": expected_ids,
            "chainsaw": chainsaw_ids,
        }
    expected_failures = (
        {"chainsaw": expected_ids - chainsaw_ids}
        if role == "legacy-accuracy" and expected_ids - chainsaw_ids
        else {}
    )
    if manifest["success_id_counts"] != {
        name: len(values) for name, values in sorted(expected_success.items())
    }:
        raise InvalidEvidence(f"locked {role} success counts mismatch")
    if manifest["success_id_set_sha256s"] != {
        name: canonical_id_set_hash(values)
        for name, values in sorted(expected_success.items())
    }:
        raise InvalidEvidence(f"locked {role} success hashes mismatch")
    if manifest["failure_id_counts"] != {
        name: len(values) for name, values in sorted(expected_failures.items())
    }:
        raise InvalidEvidence(f"locked {role} failure counts mismatch")
    if manifest["failure_id_set_sha256s"] != {
        name: canonical_id_set_hash(values)
        for name, values in sorted(expected_failures.items())
    }:
        raise InvalidEvidence(f"locked {role} failure hashes mismatch")

    if len(failure_rows) != manifest["row_counts"]["failures"]:
        raise InvalidEvidence(f"locked {role} failure-row count mismatch")
    if role != "legacy-accuracy" and failure_rows:
        raise InvalidEvidence(f"locked {role} unexpectedly records failures")
    observed_failure_ids: list[str] = []
    for row in failure_rows:
        if (
            row["dataset"] != "cath663"
            or row["tool"] != "chainsaw"
            or row["stage"] != "expected_failure"
            or row["message"] != "chainsaw_expected_failure"
            or row["entry_id"] not in expected_ids - chainsaw_ids
        ):
            raise InvalidEvidence("locked Chainsaw failure row is invalid")
        observed_failure_ids.append(row["entry_id"])
    expected_failure_ids = (
        expected_ids - chainsaw_ids if role == "legacy-accuracy" else set()
    )
    if observed_failure_ids != sorted(expected_failure_ids) or len(
        observed_failure_ids
    ) != len(set(observed_failure_ids)):
        raise InvalidEvidence("locked Chainsaw failure coverage mismatch")

    if score_rows is None:
        if manifest["row_counts"]["scores"] != 0:
            raise InvalidEvidence("resource manifest has nonzero score rows")
        return
    if len(score_rows) != manifest["row_counts"]["scores"]:
        raise InvalidEvidence(f"locked {role} score-row count mismatch")


def _validate_score_coverage(
    rows: list[dict[str, str]],
    *,
    role: str,
    metadata: Mapping[str, object],
    expected_ids: set[str],
    chainsaw_ids: set[str],
) -> None:
    allowed_tools = {"sword2-rust"} if role == "factorized-accuracy" else {
        "sword2-rust",
        "merizo",
        "chainsaw",
    }
    for row in rows:
        entry_id = row["entry_id"]
        if entry_id not in expected_ids or row["dataset"] != "cath663" or row["tool"] not in allowed_tools:
            raise InvalidEvidence(f"locked {role} score identity is invalid")
        entry = metadata[entry_id]
        if row["pdb_id"] != entry.pdb_id or row["chain_id"] != entry.chain_id:
            raise InvalidEvidence(f"locked {role} score metadata identity mismatch")
    _unique_exact_rows(
        rows,
        lambda row: row["tool"] == "sword2-rust"
        and row["variant"] == "optimal"
        and row["partition"] == "Optimal partition",
        expected_ids,
        f"locked {role} rank-1 SWORD",
    )
    if role == "legacy-accuracy":
        _unique_exact_rows(
            rows,
            lambda row: row["tool"] == "merizo",
            expected_ids,
            "locked Merizo",
        )
        _unique_exact_rows(
            rows,
            lambda row: row["tool"] == "chainsaw",
            chainsaw_ids,
            "locked Chainsaw",
        )


def _validate_expected_chainsaw(
    manifest: Mapping[str, object],
    expected_path: Path,
    expected_ids: set[str],
    chainsaw_ids: set[str],
) -> None:
    role = str(manifest["locked_role"])
    record = manifest["expected_chainsaw"]
    if role != "legacy-accuracy":
        if record is not None:
            raise InvalidEvidence(f"locked {role} unexpectedly binds Chainsaw IDs")
        return
    if not isinstance(record, dict) or set(record) != {
        "file_sha256",
        "success_count",
        "success_id_set_sha256",
        "failure_count",
        "failure_id_set_sha256",
    }:
        raise InvalidEvidence("locked Chainsaw expectation schema mismatch")
    failures = expected_ids - chainsaw_ids
    if record != {
        "file_sha256": sha256_file(expected_path),
        "success_count": len(chainsaw_ids),
        "success_id_set_sha256": canonical_id_set_hash(chainsaw_ids),
        "failure_count": len(failures),
        "failure_id_set_sha256": canonical_id_set_hash(failures),
    }:
        raise InvalidEvidence("locked Chainsaw expectation differs from frozen input")


def _manifest_common_identity(manifest: Mapping[str, object]) -> dict[str, object]:
    return {
        field: manifest[field]
        for field in (
            "dataset",
            "dataset_sha256",
            "dataset_id_count",
            "dataset_ids",
            "dataset_id_set_sha256",
            "structure_sha256s",
            "structure_tree_sha256",
            "model_manifest_sha256",
            "runtime_manifest_sha256",
            "eligibility_manifest_sha256",
            "eligibility_policy",
            "factorized_eligible_count",
            "factorized_eligible_id_set_sha256",
            "structural_abstention_count",
            "structural_abstention_id_set_sha256",
            "binary_sha256",
            "runtime_source_git_commit",
            "runtime_input_tree_sha256",
            "evidence_tool_tree_sha256",
            "cargo_lock_sha256",
        )
    }


def _evidence_path_arguments(args: argparse.Namespace) -> dict[str, Path]:
    return {name: Path(getattr(args, name)) for name in EVIDENCE_ARGUMENT_NAMES}


def _build_coverage_attestation(
    args: argparse.Namespace,
    *,
    expected_population: int | None = None,
) -> dict[str, object]:
    if expected_population is None:
        expected_population = LOCKED_POPULATION_SIZE
    if type(expected_population) is not int or expected_population <= 0:
        raise InvalidEvidence("locked population size is invalid")
    raw_paths = _evidence_path_arguments(args)
    paths = _validate_input_paths(raw_paths)
    initial_evidence_hashes = {
        role: sha256_file(path) for role, path in sorted(paths.items())
    }
    result_roots = {
        paths["legacy_manifest"].parent,
        paths["factorized_manifest"].parent,
        paths["resource_manifest"].parent,
    }
    if len(result_roots) != 3:
        raise InvalidEvidence("locked result roots alias")

    try:
        model_manifest = verify_top_level_manifest(paths["model_manifest"])
    except (OSError, ValueError) as error:
        raise InvalidEvidence("model manifest graph is invalid") from error
    repo_root = Path(__file__).resolve().parents[1]
    binary = repo_root / "target/release/sword2"
    try:
        runtime_manifest = verify_runtime_freeze(
            paths["runtime_manifest"],
            binary=binary,
            repo_root=repo_root,
        )
    except (OSError, ValueError) as error:
        raise InvalidEvidence("runtime manifest graph is invalid") from error
    model_hash = sha256_file(paths["model_manifest"])
    runtime_hash = sha256_file(paths["runtime_manifest"])
    if runtime_manifest["model_manifest_sha256"] != model_hash:
        raise InvalidEvidence("runtime manifest binds another model manifest")
    if (
        runtime_manifest["model_artifact_sha256s"]["standalone_baseline_sha256"]
        != model_manifest["standalone_baseline_sha256"]
        or sha256_file(paths["standalone_baseline"])
        != model_manifest["standalone_baseline_sha256"]
    ):
        raise InvalidEvidence("standalone baseline differs from frozen model graph")

    metadata_entries = _read_locked_metadata(paths["dataset_metadata"])
    metadata_ids = [entry.entry_id for entry in metadata_entries]
    try:
        metadata_id_hash = canonical_id_set_hash(metadata_ids)
    except ValueError as error:
        raise InvalidEvidence("locked dataset identities are invalid") from error
    if len(metadata_ids) != expected_population:
        raise InvalidEvidence(
            f"locked population must contain exactly {expected_population} entries"
        )
    metadata = {entry.entry_id: entry for entry in metadata_entries}
    expected_ids = set(metadata)
    try:
        eligibility = load_eligibility_manifest(paths["eligibility_manifest"])
    except (OSError, ValueError) as error:
        raise InvalidEvidence("eligibility manifest is invalid") from error
    eligibility_hash = sha256_file(paths["eligibility_manifest"])
    if (
        eligibility["dataset"] != "cath663"
        or eligibility["dataset_sha256"]
        != sha256_file(paths["dataset_metadata"])
        or eligibility["dataset_ids"] != metadata_ids
        or eligibility["dataset_id_count"] != expected_population
        or eligibility["dataset_id_set_sha256"] != metadata_id_hash
        or eligibility["runtime_manifest_sha256"] != runtime_hash
        or eligibility["runtime_source_git_commit"]
        != runtime_manifest["runtime_source_git_commit"]
        or eligibility["binary_sha256"] != runtime_manifest["binary_sha256"]
        or eligibility["policy"] != FACTORIZED_ELIGIBILITY_POLICY
    ):
        raise InvalidEvidence("eligibility authority chain mismatch")
    factorized_eligible_ids = set(eligibility["eligible_ids"])
    structural_abstention_ids = set(eligibility["ineligible_ids"])
    if (
        factorized_eligible_ids & structural_abstention_ids
        or factorized_eligible_ids | structural_abstention_ids != expected_ids
    ):
        raise InvalidEvidence("eligibility IDs do not partition the dataset")
    chainsaw_tuple = _load_expected_ids(paths["chainsaw_expected_success_ids"])
    chainsaw_ids = set(chainsaw_tuple)
    if not chainsaw_ids or not chainsaw_ids.issubset(expected_ids):
        raise InvalidEvidence("frozen Chainsaw success IDs are not a nonempty dataset subset")

    manifests = {
        "legacy": _load_locked_manifest(paths["legacy_manifest"], "legacy-accuracy"),
        "factorized": _load_locked_manifest(
            paths["factorized_manifest"], "factorized-accuracy"
        ),
        "resource": _load_locked_manifest(
            paths["resource_manifest"], "paired-sword-resources"
        ),
    }
    common = _manifest_common_identity(manifests["legacy"])
    if any(_manifest_common_identity(manifest) != common for manifest in manifests.values()):
        raise InvalidEvidence("locked manifests do not bind one common input/runtime chain")
    if (
        common["dataset_sha256"] != sha256_file(paths["dataset_metadata"])
        or common["dataset_id_count"] != expected_population
        or common["dataset_ids"] != sorted(metadata_ids)
        or common["dataset_id_set_sha256"] != metadata_id_hash
        or common["model_manifest_sha256"] != model_hash
        or common["runtime_manifest_sha256"] != runtime_hash
        or common["eligibility_manifest_sha256"] != eligibility_hash
        or common["eligibility_policy"] != eligibility["policy"]
        or common["factorized_eligible_count"]
        != len(factorized_eligible_ids)
        or common["factorized_eligible_id_set_sha256"]
        != eligibility["eligible_id_set_sha256"]
        or common["structural_abstention_count"]
        != len(structural_abstention_ids)
        or common["structural_abstention_id_set_sha256"]
        != eligibility["ineligible_id_set_sha256"]
        or common["structure_sha256s"] != eligibility["structure_sha256s"]
        or common["structure_tree_sha256"]
        != eligibility["structure_tree_sha256"]
        or common["binary_sha256"] != runtime_manifest["binary_sha256"]
        or common["runtime_source_git_commit"]
        != runtime_manifest["runtime_source_git_commit"]
        or common["runtime_input_tree_sha256"]
        != runtime_manifest["runtime_input_tree_sha256"]
        or common["evidence_tool_tree_sha256"]
        != runtime_manifest["evidence_tool_tree_sha256"]
        or common["cargo_lock_sha256"] != runtime_manifest["cargo_lock_sha256"]
    ):
        raise InvalidEvidence("locked manifest authority chain mismatch")

    score_rows: dict[str, list[dict[str, str]]] = {}
    run_rows: dict[str, list[dict[str, str]]] = {}
    failure_rows: dict[str, list[dict[str, str]]] = {}
    for short_role, manifest in manifests.items():
        prefix = "resource" if short_role == "resource" else short_role
        scores_path = None if short_role == "resource" else paths[f"{prefix}_scores"]
        runs_path = paths[f"{prefix}_runs"]
        failures_path = paths[f"{prefix}_failures"]
        _validate_manifest_files(
            paths[f"{prefix}_manifest"],
            manifest,
            scores_path=scores_path,
            runs_path=runs_path,
            failures_path=failures_path,
        )
        if scores_path is not None:
            score_rows[short_role] = _score_identity_rows(scores_path)
        _, run_rows[short_role] = _read_dict_csv(runs_path, LOCKED_RUN_HEADER)
        _, failure_rows[short_role] = _read_dict_csv(failures_path, FAILURE_HEADER)
        for row in run_rows[short_role]:
            entry = metadata.get(row["entry_id"])
            if (
                entry is None
                or row["pdb_id"] != entry.pdb_id
                or row["chain_id"] != entry.chain_id
            ):
                raise InvalidEvidence(f"locked {short_role} run metadata identity mismatch")
        for row in failure_rows[short_role]:
            entry = metadata.get(row["entry_id"])
            if (
                entry is None
                or row["pdb_id"] != entry.pdb_id
                or row["chain_id"] != entry.chain_id
            ):
                raise InvalidEvidence(f"locked {short_role} failure metadata identity mismatch")
        _validate_external_artifacts(manifest)
        _validate_expected_chainsaw(
            manifest,
            paths["chainsaw_expected_success_ids"],
            expected_ids,
            chainsaw_ids,
        )
        _validate_manifest_accounting(
            manifest,
            expected_ids,
            factorized_eligible_ids,
            chainsaw_ids,
            score_rows.get(short_role),
            failure_rows[short_role],
        )
        _validate_role_runs(
            run_rows[short_role],
            manifest,
            expected_ids,
            factorized_eligible_ids,
            structural_abstention_ids,
            chainsaw_ids,
            paths[f"{prefix}_manifest"].parent,
        )

    _validate_score_coverage(
        score_rows["legacy"],
        role="legacy-accuracy",
        metadata=metadata,
        expected_ids=expected_ids,
        chainsaw_ids=chainsaw_ids,
    )
    _validate_score_coverage(
        score_rows["factorized"],
        role="factorized-accuracy",
        metadata=metadata,
        expected_ids=factorized_eligible_ids,
        chainsaw_ids=chainsaw_ids,
    )

    selector_counts = {
        name: manifest["selector_counts"] for name, manifest in manifests.items()
    }
    evidence_hashes = {
        role: sha256_file(path) for role, path in sorted(paths.items())
    }
    if evidence_hashes != initial_evidence_hashes:
        raise InvalidEvidence("coverage evidence changed while validating")
    resource_order = manifests["resource"]["order_design"]
    if not isinstance(resource_order, dict):
        raise InvalidEvidence("resource order design is invalid")
    return {
        "schema_version": COVERAGE_SCHEMA_VERSION,
        "valid": True,
        "roles": {
            "legacy": "legacy-accuracy",
            "factorized": "factorized-accuracy",
            "resource": "paired-sword-resources",
        },
        "model_manifest_sha256": model_hash,
        "runtime_manifest_sha256": runtime_hash,
        "eligibility_manifest_sha256": eligibility_hash,
        "eligibility_policy": str(eligibility["policy"]),
        "binary_sha256": runtime_manifest["binary_sha256"],
        "dataset_sha256": sha256_file(paths["dataset_metadata"]),
        "dataset_id_count": expected_population,
        "dataset_id_set_sha256": metadata_id_hash,
        "factorized_eligible_count": len(factorized_eligible_ids),
        "factorized_eligible_id_set_sha256": eligibility[
            "eligible_id_set_sha256"
        ],
        "structural_abstention_count": len(structural_abstention_ids),
        "structural_abstention_id_set_sha256": eligibility[
            "ineligible_id_set_sha256"
        ],
        "factorized_structural_coverage": (
            len(factorized_eligible_ids) / expected_population
        ),
        "chainsaw_success_count": len(chainsaw_ids),
        "chainsaw_success_id_set_sha256": canonical_id_set_hash(chainsaw_ids),
        "chainsaw_failure_count": len(expected_ids - chainsaw_ids),
        "chainsaw_failure_id_set_sha256": canonical_id_set_hash(
            expected_ids - chainsaw_ids
        ),
        "resource_pair_count": len(factorized_eligible_ids),
        "resource_order_sha256": resource_order["assignments_sha256"],
        "selector_counts": selector_counts,
        "fallback_count": sum(
            int(manifest["fallback_count"]) for manifest in manifests.values()
        ),
        "excluded_candidate_count": sum(
            int(manifest["excluded_candidate_count"])
            for manifest in manifests.values()
        ),
        "evidence_sha256s": evidence_hashes,
    }


def _validate_coverage_payload(payload: Mapping[str, object]) -> None:
    if set(payload) != COVERAGE_KEYS:
        raise InvalidEvidence("coverage attestation schema mismatch")
    if (
        type(payload["schema_version"]) is not int
        or payload["schema_version"] != COVERAGE_SCHEMA_VERSION
        or payload["valid"] is not True
    ):
        raise InvalidEvidence("coverage attestation is not valid schema v2")
    for field in (
        "model_manifest_sha256",
        "runtime_manifest_sha256",
        "eligibility_manifest_sha256",
        "binary_sha256",
        "dataset_sha256",
        "dataset_id_set_sha256",
        "factorized_eligible_id_set_sha256",
        "structural_abstention_id_set_sha256",
        "chainsaw_success_id_set_sha256",
        "chainsaw_failure_id_set_sha256",
        "resource_order_sha256",
    ):
        _hash_value(payload[field], f"coverage {field}")
    for field in (
        "dataset_id_count",
        "factorized_eligible_count",
        "structural_abstention_count",
        "chainsaw_success_count",
        "chainsaw_failure_count",
        "resource_pair_count",
        "fallback_count",
        "excluded_candidate_count",
    ):
        _ordinary_int(payload[field], f"coverage {field}")
    if payload["eligibility_policy"] != FACTORIZED_ELIGIBILITY_POLICY:
        raise InvalidEvidence("coverage eligibility policy mismatch")
    structural_coverage = _finite_float(
        payload["factorized_structural_coverage"],
        "coverage factorized structural coverage",
    )
    if payload["roles"] != {
        "legacy": "legacy-accuracy",
        "factorized": "factorized-accuracy",
        "resource": "paired-sword-resources",
    }:
        raise InvalidEvidence("coverage role mapping mismatch")
    evidence = payload["evidence_sha256s"]
    if not isinstance(evidence, dict) or set(evidence) != set(EVIDENCE_ARGUMENT_NAMES):
        raise InvalidEvidence("coverage evidence hash mapping mismatch")
    for value in evidence.values():
        _hash_value(value, "coverage evidence hash")
    if (
        payload["factorized_eligible_count"]
        + payload["structural_abstention_count"]
        != payload["dataset_id_count"]
        or payload["resource_pair_count"]
        != payload["factorized_eligible_count"]
        or payload["chainsaw_success_count"] + payload["chainsaw_failure_count"]
        != payload["dataset_id_count"]
        or payload["fallback_count"]
        != payload["structural_abstention_count"]
        or structural_coverage
        != payload["factorized_eligible_count"] / payload["dataset_id_count"]
    ):
        raise InvalidEvidence("coverage denominator/fallback accounting mismatch")
    selector_counts = payload["selector_counts"]
    if not isinstance(selector_counts, dict) or set(selector_counts) != {
        "legacy",
        "factorized",
        "resource",
    }:
        raise InvalidEvidence("coverage selector-count mapping mismatch")
    if any(
        not isinstance(record, dict) or set(record) != {"legacy", "factorized"}
        for record in selector_counts.values()
    ):
        raise InvalidEvidence("coverage selector-count records are invalid")
    population = payload["dataset_id_count"]
    eligible = payload["factorized_eligible_count"]
    if selector_counts != {
        "legacy": {"legacy": population, "factorized": 0},
        "factorized": {"legacy": 0, "factorized": population},
        "resource": {"legacy": eligible, "factorized": eligible},
    }:
        raise InvalidEvidence("coverage selector counts do not prove complete roles")


def _write_canonical_absent(path: Path, payload: Mapping[str, object]) -> None:
    target = Path(path).absolute()
    data = canonical_json_bytes(dict(payload))
    temporary = _write_owned_temporary(target, data)
    try:
        os.link(temporary, target)
    finally:
        temporary.unlink(missing_ok=True)
    if target.read_bytes() != data:
        raise OSError("installed canonical authority bytes differ")


def _load_revalidated_coverage(args: argparse.Namespace) -> dict[str, object]:
    try:
        stored = load_canonical_json(Path(args.coverage_attestation))
    except (OSError, ValueError) as error:
        raise InvalidEvidence("coverage attestation is not canonical") from error
    _validate_coverage_payload(stored)
    current = _build_coverage_attestation(args)
    _validate_coverage_payload(current)
    if stored != current:
        raise InvalidEvidence("coverage attestation no longer matches evidence bytes")
    return stored


def _finite_float(value: object, description: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{description} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{description} must be a finite number") from error
    if not math.isfinite(number):
        raise ValueError(f"{description} must be finite")
    return number


def _validated_index(frame: pd.DataFrame, description: str) -> list[str]:
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise ValueError(f"{description} must be a nonempty data frame")
    if frame.index.has_duplicates:
        raise ValueError(f"{description} contains duplicate entry IDs")
    ids: list[str] = []
    for value in frame.index.tolist():
        if (
            not isinstance(value, str)
            or not value
            or value.strip() != value
            or unicodedata.normalize("NFC", value) != value
            or "\0" in value
            or any(ord(character) < 32 or ord(character) == 127 for character in value)
        ):
            raise ValueError(f"{description} contains a noncanonical entry ID")
        ids.append(value)
    if ids != sorted(ids):
        frame.sort_index(inplace=True)
        ids = frame.index.tolist()
    return ids


def _numeric_column(frame: pd.DataFrame, column: str, description: str) -> pd.Series:
    if column not in frame.columns:
        raise ValueError(f"{description} lacks {column}")
    values = pd.to_numeric(frame[column], errors="raise").astype(float)
    if not np.isfinite(values.to_numpy(dtype=np.float64)).all():
        raise ValueError(f"{description} {column} is nonfinite")
    return values


def _source_hashes(*frames: pd.DataFrame) -> dict[str, str]:
    result: dict[str, str] = {}
    for index, frame in enumerate(frames):
        role = frame.attrs.get("source_role", f"input_{index}")
        digest = frame.attrs.get("source_sha256")
        if digest is None:
            continue
        if (
            not isinstance(role, str)
            or not role
            or not isinstance(digest, str)
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            raise ValueError("data-frame source hash metadata is invalid")
        if role in result and result[role] != digest:
            raise ValueError("data-frame source hash roles conflict")
        result[role] = digest
    return dict(sorted(result.items()))


def _gate(
    name: str,
    value: float,
    threshold: float,
    comparison: str,
    denominator: int,
    method: str,
    source_sha256s: Mapping[str, str],
) -> dict[str, object]:
    value = _finite_float(value, name)
    threshold = _finite_float(threshold, f"{name} threshold")
    if type(denominator) is not int or denominator <= 0:
        raise ValueError(f"{name} denominator must be positive")
    if comparison == ">":
        passed = value > threshold
    elif comparison == ">=":
        passed = value >= threshold
    elif comparison == "<=":
        passed = value <= threshold
    else:
        raise ValueError(f"{name} has an invalid comparison")
    return {
        "name": name,
        "measured": True,
        "value": value,
        "threshold": threshold,
        "comparison": comparison,
        "denominator": denominator,
        "method": method,
        "passed": bool(passed),
        "source_sha256s": dict(sorted(source_sha256s.items())),
    }


def evaluate_accuracy_gates(
    scores: pd.DataFrame,
    baseline_predictions: pd.DataFrame,
    competitors: Mapping[str, pd.DataFrame],
) -> dict[str, object]:
    """Evaluate the seven frozen accuracy gates on prevalidated exact rows."""
    factorized_ids = _validated_index(scores, "factorized scores")
    baseline_ids = _validated_index(baseline_predictions, "standalone baseline")
    if baseline_ids != factorized_ids:
        raise ValueError("standalone and factorized ID sets differ")
    if not isinstance(competitors, Mapping) or set(competitors) != {"merizo", "chainsaw"}:
        raise ValueError("competitors must contain exactly merizo and chainsaw")
    merizo = competitors["merizo"]
    chainsaw = competitors["chainsaw"]
    merizo_ids = _validated_index(merizo, "Merizo scores")
    chainsaw_ids = _validated_index(chainsaw, "Chainsaw scores")
    if merizo_ids != factorized_ids:
        raise ValueError("Merizo and factorized ID sets differ")
    if not set(chainsaw_ids).issubset(factorized_ids) or not chainsaw_ids:
        raise ValueError("Chainsaw IDs are not a nonempty frozen subset")

    ndo = _numeric_column(scores, "ndo", "factorized scores")
    standalone_ndo = _numeric_column(
        baseline_predictions,
        "ndo",
        "standalone baseline",
    )
    merizo_ndo = _numeric_column(merizo, "ndo", "Merizo scores")
    chainsaw_ndo = _numeric_column(chainsaw, "ndo", "Chainsaw scores")
    d_count = _numeric_column(scores, "d_count_acc", "factorized scores")
    predicted_count = _numeric_column(scores, "n_pred_domains", "factorized scores")
    true_count = _numeric_column(scores, "n_true_domains", "factorized scores")
    boundary = _numeric_column(scores, "boundary_f1_10", "factorized scores")
    if not d_count.isin([0.0, 1.0]).all():
        raise ValueError("d_count_acc must be exactly zero or one")
    expected_count = (predicted_count == true_count).astype(float)
    if not np.array_equal(d_count.to_numpy(), expected_count.to_numpy()):
        raise ValueError("d_count_acc disagrees with predicted/true count equality")
    if "continuity_cohort" not in scores.columns:
        raise ValueError("factorized scores lack truth continuity cohorts")
    cohorts = scores["continuity_cohort"]
    if not cohorts.isin(["contiguous", "discontinuous"]).all():
        raise ValueError("truth continuity cohorts are invalid")
    contiguous_ids = sorted(cohorts[cohorts == "contiguous"].index)
    discontinuous_ids = sorted(cohorts[cohorts == "discontinuous"].index)
    if (
        not contiguous_ids
        or not discontinuous_ids
        or set(contiguous_ids) & set(discontinuous_ids)
        or set(contiguous_ids) | set(discontinuous_ids) != set(factorized_ids)
    ):
        raise ValueError("truth continuity cohorts are not a nonempty partition")

    merizo_bootstrap = paired_chain_bootstrap(
        merizo_ndo.to_dict(),
        ndo.to_dict(),
        n_resamples=10_000,
        seed=37,
        higher_is_better=True,
    )
    chain_factorized = ndo.loc[chainsaw_ids]
    chainsaw_bootstrap = paired_chain_bootstrap(
        chainsaw_ndo.to_dict(),
        chain_factorized.to_dict(),
        n_resamples=10_000,
        seed=37,
        higher_is_better=True,
    )
    if merizo_bootstrap.n != len(factorized_ids) or chainsaw_bootstrap.n != len(chainsaw_ids):
        raise ValueError("paired bootstrap denominator mismatch")

    factorized_hashes = _source_hashes(scores)
    standalone_hashes = _source_hashes(scores, baseline_predictions)
    merizo_hashes = _source_hashes(scores, merizo)
    chainsaw_hashes = _source_hashes(scores, chainsaw)
    gates = {
        "overall_ndo": _gate(
            "overall_ndo",
            float(ndo.mean()),
            0.8389,
            ">",
            len(factorized_ids),
            "arithmetic_mean",
            factorized_hashes,
        ),
        "merizo_ndo_ci_low": _gate(
            "merizo_ndo_ci_low",
            merizo_bootstrap.ci_low,
            0.0,
            ">",
            merizo_bootstrap.n,
            "paired_chain_bootstrap_ci_low",
            merizo_hashes,
        ),
        "chainsaw_ndo_ci_low": _gate(
            "chainsaw_ndo_ci_low",
            chainsaw_bootstrap.ci_low,
            0.0,
            ">",
            chainsaw_bootstrap.n,
            "paired_chain_bootstrap_ci_low",
            chainsaw_hashes,
        ),
        "domain_count_accuracy": _gate(
            "domain_count_accuracy",
            float(d_count.mean()),
            0.745,
            ">=",
            len(factorized_ids),
            "arithmetic_mean",
            factorized_hashes,
        ),
        "boundary_f1_10": _gate(
            "boundary_f1_10",
            float(boundary.mean()),
            0.620,
            ">=",
            len(factorized_ids),
            "arithmetic_mean",
            factorized_hashes,
        ),
        "contiguous_standalone_ndo_delta": _gate(
            "contiguous_standalone_ndo_delta",
            float(ndo.loc[contiguous_ids].mean() - standalone_ndo.loc[contiguous_ids].mean()),
            -0.005,
            ">=",
            len(contiguous_ids),
            "paired_cohort_mean_delta",
            standalone_hashes,
        ),
        "discontinuous_standalone_ndo_delta": _gate(
            "discontinuous_standalone_ndo_delta",
            float(
                ndo.loc[discontinuous_ids].mean()
                - standalone_ndo.loc[discontinuous_ids].mean()
            ),
            -0.005,
            ">=",
            len(discontinuous_ids),
            "paired_cohort_mean_delta",
            standalone_hashes,
        ),
    }
    diagnostics = {
        "merizo_bootstrap": {
            "mean_delta": merizo_bootstrap.mean_delta,
            "ci_high": merizo_bootstrap.ci_high,
        },
        "chainsaw_bootstrap": {
            "mean_delta": chainsaw_bootstrap.mean_delta,
            "ci_high": chainsaw_bootstrap.ci_high,
        },
        "contiguous": {
            "denominator": len(contiguous_ids),
            "id_set_sha256": canonical_id_set_hash(contiguous_ids),
            "factorized_mean": float(ndo.loc[contiguous_ids].mean()),
            "standalone_mean": float(standalone_ndo.loc[contiguous_ids].mean()),
        },
        "discontinuous": {
            "denominator": len(discontinuous_ids),
            "id_set_sha256": canonical_id_set_hash(discontinuous_ids),
            "factorized_mean": float(ndo.loc[discontinuous_ids].mean()),
            "standalone_mean": float(standalone_ndo.loc[discontinuous_ids].mean()),
        },
    }
    return {
        "gates": gates,
        "diagnostics": diagnostics,
        "denominators": {
            "overall": len(factorized_ids),
            "chainsaw": len(chainsaw_ids),
            "contiguous": len(contiguous_ids),
            "discontinuous": len(discontinuous_ids),
        },
    }


def evaluate_resource_gates(
    baseline_runs: pd.DataFrame,
    factorized_runs: pd.DataFrame,
) -> dict[str, object]:
    """Evaluate the two frozen resource gates on exact paired SWORD rows."""
    baseline_ids = _validated_index(baseline_runs, "legacy resource rows")
    factorized_ids = _validated_index(factorized_runs, "factorized resource rows")
    if baseline_ids != factorized_ids:
        raise ValueError("resource selector ID sets differ")
    legacy_runtime = _numeric_column(baseline_runs, "runtime_s", "legacy resource rows")
    factorized_runtime = _numeric_column(
        factorized_runs,
        "runtime_s",
        "factorized resource rows",
    )
    legacy_rss = _numeric_column(baseline_runs, "peak_rss_kb", "legacy resource rows")
    factorized_rss = _numeric_column(
        factorized_runs,
        "peak_rss_kb",
        "factorized resource rows",
    )
    if any((values <= 0).any() for values in (legacy_runtime, factorized_runtime, legacy_rss, factorized_rss)):
        raise ValueError("resource inputs must be positive")
    runtime_ratios = factorized_runtime / legacy_runtime
    rss_ratios = factorized_rss / legacy_rss
    maxima_ratio = float(factorized_rss.max() / legacy_rss.max())
    hashes = _source_hashes(baseline_runs, factorized_runs)
    gates = {
        "runtime_median_ratio": _gate(
            "runtime_median_ratio",
            float(runtime_ratios.median()),
            1.15,
            "<=",
            len(baseline_ids),
            "median_of_paired_ratios",
            hashes,
        ),
        "rss_maxima_ratio": _gate(
            "rss_maxima_ratio",
            maxima_ratio,
            1.10,
            "<=",
            len(baseline_ids),
            "ratio_of_selector_maxima",
            hashes,
        ),
    }
    return {
        "gates": gates,
        "diagnostics": {
            "runtime_ratio_p95": float(
                np.quantile(runtime_ratios.to_numpy(), 0.95, method="linear")
            ),
            "rss_per_chain_ratio_median": float(rss_ratios.median()),
            "rss_per_chain_ratio_p95": float(
                np.quantile(rss_ratios.to_numpy(), 0.95, method="linear")
            ),
            "legacy_peak_rss_kb_max": float(legacy_rss.max()),
            "factorized_peak_rss_kb_max": float(factorized_rss.max()),
        },
        "denominators": {"resources": len(baseline_ids)},
    }


def _select_score_rows(
    path: Path,
    *,
    tool: str,
    expected_ids: set[str],
    rank_one: bool,
    description: str,
    source_role: str,
) -> pd.DataFrame:
    try:
        frame = pd.read_csv(
            path,
            dtype={column: "string" for column in SCORE_IDENTITY_COLUMNS},
        )
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise InvalidEvidence(f"{description} metrics CSV is invalid") from error
    required = set(SCORE_IDENTITY_COLUMNS) | {
        "ndo",
        "n_pred_domains",
        "n_true_domains",
        "d_count_acc",
        "boundary_f1_10",
    }
    if not required.issubset(frame.columns):
        raise InvalidEvidence(f"{description} metrics CSV lacks required columns")
    selected = frame[frame["tool"] == tool]
    if rank_one:
        selected = selected[
            (selected["variant"] == "optimal")
            & (selected["partition"] == "Optimal partition")
        ]
    ids = selected["entry_id"].tolist()
    if (
        len(ids) != len(set(ids))
        or set(ids) != expected_ids
        or any(not isinstance(entry_id, str) for entry_id in ids)
    ):
        raise InvalidEvidence(f"{description} metric row coverage mismatch")
    selected = selected.copy().set_index("entry_id").sort_index()
    selected.attrs["source_role"] = source_role
    selected.attrs["source_sha256"] = sha256_file(path)
    return selected


def _load_standalone_metrics(path: Path, expected_ids: set[str]) -> pd.DataFrame:
    try:
        frame = pd.read_csv(path, dtype={"chain_id": "string"})
    except (OSError, ValueError, pd.errors.ParserError) as error:
        raise InvalidEvidence("standalone baseline metrics are invalid") from error
    if "chain_id" not in frame.columns or "ndo" not in frame.columns:
        raise InvalidEvidence("standalone baseline lacks chain_id/ndo")
    ids = frame["chain_id"].tolist()
    if len(ids) != len(set(ids)) or set(ids) != expected_ids:
        raise InvalidEvidence("standalone baseline identity coverage mismatch")
    if len(frame) != len(expected_ids):
        raise InvalidEvidence("standalone baseline contains extra rows")
    result = frame.rename(columns={"chain_id": "entry_id"}).set_index("entry_id").sort_index()
    _numeric_column(result, "ndo", "standalone baseline")
    result.attrs["source_role"] = "standalone_baseline"
    result.attrs["source_sha256"] = sha256_file(path)
    return result


def _resource_metric_frames(
    path: Path,
    expected_ids: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    _, rows = _read_dict_csv(path, LOCKED_RUN_HEADER)
    frames: dict[str, pd.DataFrame] = {}
    digest = sha256_file(path)
    for selector in ("legacy", "factorized"):
        selected = [
            row
            for row in rows
            if row["tool"] == "sword2-rust" and row["selector_variant"] == selector
        ]
        ids = [row["entry_id"] for row in selected]
        if len(ids) != len(set(ids)) or set(ids) != expected_ids:
            raise InvalidEvidence(f"resource {selector} metric coverage mismatch")
        frame = pd.DataFrame(
            {
                "entry_id": ids,
                "runtime_s": [row["runtime_s"] for row in selected],
                "peak_rss_kb": [row["peak_rss_kb"] for row in selected],
            }
        ).set_index("entry_id").sort_index()
        frame.attrs["source_role"] = "resource_runs"
        frame.attrs["source_sha256"] = digest
        frames[selector] = frame
    return frames["legacy"], frames["factorized"]


def _truth_cohorts(path: Path, expected_ids: set[str]) -> dict[str, str]:
    entries = _read_locked_metadata(path)
    if {entry.entry_id for entry in entries} != expected_ids or len(entries) != len(expected_ids):
        raise InvalidEvidence("truth cohort metadata coverage mismatch")
    cohorts: dict[str, str] = {}
    for entry in entries:
        domains = parse_cath_domain_string(entry.chopping)
        if not domains or len(domains) != entry.n_domains:
            raise InvalidEvidence("truth metadata has an invalid domain chopping")
        cohorts[entry.entry_id] = (
            "discontinuous"
            if any(len(domain) > 1 for domain in domains)
            else "contiguous"
        )
    if set(cohorts.values()) != {"contiguous", "discontinuous"}:
        raise InvalidEvidence("truth continuity cohorts must both be nonempty")
    return cohorts


def _build_acceptance_payload(
    args: argparse.Namespace,
    coverage: Mapping[str, object],
) -> dict[str, object]:
    paths = _evidence_path_arguments(args)
    expected_evidence_hashes = coverage["evidence_sha256s"]
    if not isinstance(expected_evidence_hashes, dict):
        raise InvalidEvidence("coverage evidence hashes are invalid")
    if {
        role: sha256_file(path) for role, path in sorted(paths.items())
    } != expected_evidence_hashes:
        raise InvalidEvidence("metric evidence differs from validated coverage")
    expected_ids = set(_load_expected_ids(paths["chainsaw_expected_success_ids"]))
    metadata_entries = _read_locked_metadata(paths["dataset_metadata"])
    all_ids = {entry.entry_id for entry in metadata_entries}
    if len(all_ids) != coverage["dataset_id_count"]:
        raise InvalidEvidence("metric population differs from coverage")
    try:
        eligibility = load_eligibility_manifest(paths["eligibility_manifest"])
    except (OSError, ValueError) as error:
        raise InvalidEvidence("metric eligibility manifest is invalid") from error
    eligible_ids = set(eligibility["eligible_ids"])
    structural_abstention_ids = set(eligibility["ineligible_ids"])
    if (
        eligible_ids | structural_abstention_ids != all_ids
        or eligible_ids & structural_abstention_ids
        or sha256_file(paths["eligibility_manifest"])
        != coverage["eligibility_manifest_sha256"]
        or len(eligible_ids) != coverage["factorized_eligible_count"]
        or eligibility["eligible_id_set_sha256"]
        != coverage["factorized_eligible_id_set_sha256"]
        or len(structural_abstention_ids)
        != coverage["structural_abstention_count"]
        or eligibility["ineligible_id_set_sha256"]
        != coverage["structural_abstention_id_set_sha256"]
    ):
        raise InvalidEvidence("metric eligibility differs from validated coverage")
    chainsaw_ids = expected_ids & eligible_ids
    if not chainsaw_ids:
        raise InvalidEvidence(
            "no frozen Chainsaw success remains in the eligible population"
        )

    factorized = _select_score_rows(
        paths["factorized_scores"],
        tool="sword2-rust",
        expected_ids=eligible_ids,
        rank_one=True,
        description="factorized",
        source_role="factorized_scores",
    )
    legacy_full = _select_score_rows(
        paths["legacy_scores"],
        tool="sword2-rust",
        expected_ids=all_ids,
        rank_one=True,
        description="runtime legacy",
        source_role="legacy_scores",
    )
    merizo_full = _select_score_rows(
        paths["legacy_scores"],
        tool="merizo",
        expected_ids=all_ids,
        rank_one=False,
        description="merizo",
        source_role="legacy_scores",
    )
    chainsaw_full = _select_score_rows(
        paths["legacy_scores"],
        tool="chainsaw",
        expected_ids=expected_ids,
        rank_one=False,
        description="chainsaw",
        source_role="legacy_scores",
    )
    standalone_full = _load_standalone_metrics(
        paths["standalone_baseline"], all_ids
    )

    def restrict(frame: pd.DataFrame, ids: set[str]) -> pd.DataFrame:
        restricted = frame.loc[sorted(ids)].copy()
        restricted.attrs = dict(frame.attrs)
        return restricted

    legacy = restrict(legacy_full, eligible_ids)
    merizo = restrict(merizo_full, eligible_ids)
    chainsaw = restrict(chainsaw_full, chainsaw_ids)
    standalone = restrict(standalone_full, eligible_ids)
    factorized["continuity_cohort"] = pd.Series(
        _truth_cohorts(paths["dataset_metadata"], all_ids)
    ).reindex(factorized.index)
    legacy_resources, factorized_resources = _resource_metric_frames(
        paths["resource_runs"], eligible_ids
    )

    try:
        accuracy = evaluate_accuracy_gates(
            factorized,
            standalone,
            {"merizo": merizo, "chainsaw": chainsaw},
        )
        resources = evaluate_resource_gates(legacy_resources, factorized_resources)
    except (TypeError, ValueError) as error:
        raise InvalidEvidence("gate inputs or derived values are invalid") from error
    gates = {**accuracy["gates"], **resources["gates"]}
    if list(gates) != list(GATE_ORDER):
        raise InvalidEvidence("gate construction order mismatch")
    diagnostics = {
        "accuracy": accuracy["diagnostics"],
        "resources": resources["diagnostics"],
        "structural_coverage": {
            "full_denominator": len(all_ids),
            "factorized_eligible_count": len(eligible_ids),
            "structural_abstention_count": len(structural_abstention_ids),
            "factorized_structural_coverage": coverage[
                "factorized_structural_coverage"
            ],
            "ineligibility_reason_counts": eligibility[
                "ineligibility_reason_counts"
            ],
        },
        "conditional_mean_ndo": {
            "factorized": float(
                _numeric_column(factorized, "ndo", "factorized scores").mean()
            ),
            "runtime_legacy": float(
                _numeric_column(legacy, "ndo", "runtime legacy scores").mean()
            ),
            "merizo": float(
                _numeric_column(merizo, "ndo", "Merizo scores").mean()
            ),
            "chainsaw": float(
                _numeric_column(chainsaw, "ndo", "Chainsaw scores").mean()
            ),
            "chainsaw_denominator": len(chainsaw_ids),
        },
        "full_population_mean_ndo": {
            "runtime_legacy": float(
                _numeric_column(
                    legacy_full, "ndo", "full runtime legacy scores"
                ).mean()
            ),
            "merizo": float(
                _numeric_column(merizo_full, "ndo", "full Merizo scores").mean()
            ),
            "chainsaw": float(
                _numeric_column(
                    chainsaw_full, "ndo", "frozen-success Chainsaw scores"
                ).mean()
            ),
            "chainsaw_denominator": len(expected_ids),
        },
        "selector_counts": coverage["selector_counts"],
        "fallback_count": coverage["fallback_count"],
        "excluded_candidate_count": coverage["excluded_candidate_count"],
        "resource_order_sha256": coverage["resource_order_sha256"],
        "source_sha256s": dict(sorted(coverage["evidence_sha256s"].items())),
    }
    runtime_payload = load_canonical_json(paths["runtime_manifest"])
    diagnostics["default_promotion_compatible"] = bool(
        coverage["factorized_structural_coverage"] == 1.0
        and runtime_payload.get("selector_cache_contract") == CACHE_CONTRACT
        and CACHE_CONTRACT.get("default_promotion_compatible") is True
    )
    denominators = {
        **accuracy["denominators"],
        **resources["denominators"],
    }
    all_measured = all(record["measured"] is True for record in gates.values())
    all_pass = all(record["measured"] is True and record["passed"] is True for record in gates.values())
    if {
        role: sha256_file(path) for role, path in sorted(paths.items())
    } != expected_evidence_hashes:
        raise InvalidEvidence("metric evidence changed while evaluating gates")
    payload = {
        "schema_version": ACCEPTANCE_SCHEMA_VERSION,
        "gate_order": list(GATE_ORDER),
        "bootstrap": BOOTSTRAP,
        "model_manifest_sha256": coverage["model_manifest_sha256"],
        "runtime_manifest_sha256": coverage["runtime_manifest_sha256"],
        "coverage_attestation_sha256": sha256_file(Path(args.coverage_attestation)),
        "evidence_sha256s": dict(sorted(coverage["evidence_sha256s"].items())),
        "denominators": denominators,
        "gates": gates,
        "diagnostics": diagnostics,
        "all_gates_measured": all_measured,
        "all_gates_pass": all_pass,
    }
    _validate_acceptance_payload(payload)
    return payload


def _validate_acceptance_payload(payload: Mapping[str, object]) -> None:
    if set(payload) != ACCEPTANCE_KEYS:
        raise InvalidEvidence("acceptance top-level schema mismatch")
    if (
        type(payload["schema_version"]) is not int
        or payload["schema_version"] != ACCEPTANCE_SCHEMA_VERSION
    ):
        raise InvalidEvidence("acceptance schema version mismatch")
    if payload["gate_order"] != list(GATE_ORDER) or payload["bootstrap"] != BOOTSTRAP:
        raise InvalidEvidence("acceptance gate/bootstrap contract mismatch")
    for field in (
        "model_manifest_sha256",
        "runtime_manifest_sha256",
        "coverage_attestation_sha256",
    ):
        _hash_value(payload[field], f"acceptance {field}")
    evidence = payload["evidence_sha256s"]
    if not isinstance(evidence, dict) or set(evidence) != set(EVIDENCE_ARGUMENT_NAMES):
        raise InvalidEvidence("acceptance evidence mapping mismatch")
    for digest in evidence.values():
        _hash_value(digest, "acceptance evidence hash")
    denominators = payload["denominators"]
    if not isinstance(denominators, dict) or set(denominators) != {
        "overall",
        "chainsaw",
        "contiguous",
        "discontinuous",
        "resources",
    }:
        raise InvalidEvidence("acceptance denominator schema mismatch")
    for name, value in denominators.items():
        _ordinary_int(value, f"acceptance {name} denominator", 1)
    if (
        denominators["overall"] != denominators["resources"]
        or denominators["contiguous"] + denominators["discontinuous"]
        != denominators["overall"]
        or denominators["chainsaw"] > denominators["overall"]
    ):
        raise InvalidEvidence("acceptance denominators are inconsistent")
    gates = payload["gates"]
    if not isinstance(gates, dict) or set(gates) != set(GATE_ORDER):
        raise InvalidEvidence("acceptance gate mapping mismatch")
    expected_rules = {
        "overall_ndo": (0.8389, ">", "arithmetic_mean"),
        "merizo_ndo_ci_low": (0.0, ">", "paired_chain_bootstrap_ci_low"),
        "chainsaw_ndo_ci_low": (0.0, ">", "paired_chain_bootstrap_ci_low"),
        "domain_count_accuracy": (0.745, ">=", "arithmetic_mean"),
        "boundary_f1_10": (0.620, ">=", "arithmetic_mean"),
        "contiguous_standalone_ndo_delta": (-0.005, ">=", "paired_cohort_mean_delta"),
        "discontinuous_standalone_ndo_delta": (-0.005, ">=", "paired_cohort_mean_delta"),
        "runtime_median_ratio": (1.15, "<=", "median_of_paired_ratios"),
        "rss_maxima_ratio": (1.10, "<=", "ratio_of_selector_maxima"),
    }
    expected_denominator_roles = {
        "overall_ndo": "overall",
        "merizo_ndo_ci_low": "overall",
        "chainsaw_ndo_ci_low": "chainsaw",
        "domain_count_accuracy": "overall",
        "boundary_f1_10": "overall",
        "contiguous_standalone_ndo_delta": "contiguous",
        "discontinuous_standalone_ndo_delta": "discontinuous",
        "runtime_median_ratio": "resources",
        "rss_maxima_ratio": "resources",
    }
    for name in GATE_ORDER:
        record = gates[name]
        if not isinstance(record, dict) or set(record) != {
            "name",
            "measured",
            "value",
            "threshold",
            "comparison",
            "denominator",
            "method",
            "passed",
            "source_sha256s",
        }:
            raise InvalidEvidence(f"acceptance gate {name} schema mismatch")
        threshold, comparison, method = expected_rules[name]
        if (
            record["name"] != name
            or record["measured"] is not True
            or record["threshold"] != threshold
            or record["comparison"] != comparison
            or record["method"] != method
            or type(record["passed"]) is not bool
        ):
            raise InvalidEvidence(f"acceptance gate {name} contract mismatch")
        value = _finite_float(record["value"], f"acceptance {name}")
        denominator = _ordinary_int(
            record["denominator"], f"acceptance {name} denominator", 1
        )
        if denominator != denominators[expected_denominator_roles[name]]:
            raise InvalidEvidence(f"acceptance gate {name} denominator mismatch")
        expected_pass = (
            value > threshold
            if comparison == ">"
            else value >= threshold
            if comparison == ">="
            else value <= threshold
        )
        if record["passed"] is not expected_pass:
            raise InvalidEvidence(f"acceptance gate {name} pass boolean mismatch")
        hashes = record["source_sha256s"]
        if not isinstance(hashes, dict) or not hashes:
            raise InvalidEvidence(f"acceptance gate {name} lacks source hashes")
        for role, digest in hashes.items():
            _hash_value(digest, f"acceptance gate {name} source hash")
            if evidence.get(role) != digest:
                raise InvalidEvidence(
                    f"acceptance gate {name} source hash is not evidence-bound"
                )
    expected_measured = all(record["measured"] is True for record in gates.values())
    expected_pass = all(
        record["measured"] is True and record["passed"] is True
        for record in gates.values()
    )
    if (
        payload["all_gates_measured"] is not expected_measured
        or payload["all_gates_pass"] is not expected_pass
    ):
        raise InvalidEvidence("acceptance aggregate booleans are inconsistent")
    diagnostics = payload["diagnostics"]
    expected_diagnostic_keys = {
        "accuracy",
        "resources",
        "structural_coverage",
        "conditional_mean_ndo",
        "full_population_mean_ndo",
        "selector_counts",
        "fallback_count",
        "excluded_candidate_count",
        "resource_order_sha256",
        "source_sha256s",
        "default_promotion_compatible",
    }
    if not isinstance(diagnostics, dict) or set(diagnostics) != expected_diagnostic_keys:
        raise InvalidEvidence("acceptance denominator/diagnostic sections are invalid")
    if not isinstance(diagnostics["accuracy"], dict) or not isinstance(
        diagnostics["resources"], dict
    ):
        raise InvalidEvidence("acceptance metric diagnostics are invalid")
    structural = diagnostics["structural_coverage"]
    if not isinstance(structural, dict) or set(structural) != {
        "full_denominator",
        "factorized_eligible_count",
        "structural_abstention_count",
        "factorized_structural_coverage",
        "ineligibility_reason_counts",
    }:
        raise InvalidEvidence("acceptance structural-coverage diagnostics are invalid")
    full_count = _ordinary_int(
        structural["full_denominator"], "acceptance full denominator", 1
    )
    eligible_count = _ordinary_int(
        structural["factorized_eligible_count"],
        "acceptance eligible count",
        1,
    )
    abstention_count = _ordinary_int(
        structural["structural_abstention_count"],
        "acceptance abstention count",
    )
    structural_fraction = _finite_float(
        structural["factorized_structural_coverage"],
        "acceptance structural coverage",
    )
    reasons = structural["ineligibility_reason_counts"]
    if not isinstance(reasons, dict) or any(
        not isinstance(reason, str)
        or not reason
        or type(count) is not int
        or count <= 0
        for reason, count in reasons.items()
    ):
        raise InvalidEvidence("acceptance ineligibility reasons are invalid")
    if (
        eligible_count + abstention_count != full_count
        or eligible_count != denominators["overall"]
        or eligible_count != denominators["resources"]
        or structural_fraction != eligible_count / full_count
        or sum(reasons.values()) != abstention_count
    ):
        raise InvalidEvidence("acceptance structural-coverage accounting mismatch")
    conditional = diagnostics["conditional_mean_ndo"]
    if not isinstance(conditional, dict) or set(conditional) != {
        "factorized",
        "runtime_legacy",
        "merizo",
        "chainsaw",
        "chainsaw_denominator",
    }:
        raise InvalidEvidence("acceptance conditional means are invalid")
    full_means = diagnostics["full_population_mean_ndo"]
    if not isinstance(full_means, dict) or set(full_means) != {
        "runtime_legacy",
        "merizo",
        "chainsaw",
        "chainsaw_denominator",
    }:
        raise InvalidEvidence("acceptance full-population means are invalid")
    for description, record, names in (
        (
            "conditional",
            conditional,
            ("factorized", "runtime_legacy", "merizo", "chainsaw"),
        ),
        (
            "full-population",
            full_means,
            ("runtime_legacy", "merizo", "chainsaw"),
        ),
    ):
        for name in names:
            value = _finite_float(record[name], f"acceptance {description} {name}")
            if not 0.0 <= value <= 1.0:
                raise InvalidEvidence(f"acceptance {description} {name} is outside [0,1]")
    if (
        conditional["factorized"] != gates["overall_ndo"]["value"]
        or _ordinary_int(
            conditional["chainsaw_denominator"],
            "acceptance conditional Chainsaw denominator",
            1,
        )
        != denominators["chainsaw"]
        or _ordinary_int(
            full_means["chainsaw_denominator"],
            "acceptance full Chainsaw denominator",
            1,
        )
        < denominators["chainsaw"]
    ):
        raise InvalidEvidence("acceptance mean diagnostics disagree with gates")
    expected_selector_counts = {
        "legacy": {"legacy": full_count, "factorized": 0},
        "factorized": {"legacy": 0, "factorized": full_count},
        "resource": {
            "legacy": eligible_count,
            "factorized": eligible_count,
        },
    }
    if diagnostics["selector_counts"] != expected_selector_counts:
        raise InvalidEvidence("acceptance selector-count diagnostics are invalid")
    if (
        diagnostics["fallback_count"] != abstention_count
        or type(diagnostics["excluded_candidate_count"]) is not int
        or diagnostics["excluded_candidate_count"] < 0
    ):
        raise InvalidEvidence("acceptance selector diagnostics are inconsistent")
    _hash_value(
        diagnostics["resource_order_sha256"],
        "acceptance resource-order hash",
    )
    if diagnostics["source_sha256s"] != evidence:
        raise InvalidEvidence("acceptance diagnostic source hashes mismatch")
    promotion_compatible = diagnostics["default_promotion_compatible"]
    if type(promotion_compatible) is not bool or (
        structural_fraction < 1.0 and promotion_compatible
    ):
        raise InvalidEvidence("acceptance promotion compatibility is invalid")


def _render_acceptance_markdown(results: Mapping[str, object]) -> bytes:
    diagnostics = results.get("diagnostics", {})
    coverage = (
        diagnostics.get("structural_coverage", {})
        if isinstance(diagnostics, Mapping)
        else {}
    )
    if not isinstance(coverage, Mapping):
        coverage = {}
    eligible = coverage.get("factorized_eligible_count", "unknown")
    full = coverage.get("full_denominator", "unknown")
    abstentions = coverage.get("structural_abstention_count", "unknown")
    fraction = coverage.get("factorized_structural_coverage")
    fraction_text = (
        f"{_finite_float(fraction, 'structural coverage'):.6%}"
        if fraction is not None
        else "unknown"
    )
    lines = [
        "# Factorized structural ranker acceptance",
        "",
        f"Structural coverage: {eligible}/{full} ({fraction_text}); "
        f"abstentions: {abstentions}.",
        "",
        "All performance gates below are conditional on structurally eligible chains. "
        "Legacy fallbacks for abstained chains are not factorized scores.",
        "",
        "| Gate | Value | Rule | n | Result |",
        "|---|---:|:---:|---:|:---:|",
    ]
    gates = results.get("gates", {})
    if not isinstance(gates, Mapping):
        raise ValueError("acceptance gates are not a mapping")
    order = results.get("gate_order", [])
    if not isinstance(order, list) or not all(isinstance(name, str) for name in order):
        raise ValueError("acceptance gate order is invalid")
    for name in order:
        record = gates.get(name)
        if not isinstance(record, Mapping):
            continue
        value = _finite_float(record.get("value"), f"{name} value")
        threshold = _finite_float(record.get("threshold"), f"{name} threshold")
        comparison = record.get("comparison")
        denominator = record.get("denominator")
        passed = record.get("passed")
        lines.append(
            f"| {name} | {value:.12g} | {comparison} {threshold:.12g} | {denominator} | "
            f"{'PASS' if passed is True else 'FAIL'} |"
        )
    lines.extend(
        [
            "",
            f"All gates measured: {str(results.get('all_gates_measured')).lower()}",
            f"All gates pass: {str(results.get('all_gates_pass')).lower()}",
            "",
        ]
    )
    return "\n".join(lines).encode("utf-8")


def _write_owned_temporary(path: Path, data: bytes) -> Path:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(temporary)
    with temporary.open("xb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    return temporary


def write_acceptance_report(
    json_path: Path,
    markdown_path: Path,
    results: Mapping[str, object],
) -> None:
    json_path = Path(json_path).absolute()
    markdown_path = Path(markdown_path).absolute()
    if json_path == markdown_path:
        raise ValueError("acceptance JSON and Markdown outputs alias")
    if json_path.exists() or json_path.is_symlink():
        raise FileExistsError(json_path)
    if markdown_path.exists() or markdown_path.is_symlink():
        raise FileExistsError(markdown_path)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    payload = dict(results)
    _validate_acceptance_payload(payload)
    json_data = canonical_json_bytes(payload)
    if json.loads(json_data) != payload:
        raise ValueError("acceptance JSON does not round-trip")
    markdown_data = _render_acceptance_markdown(payload)
    json_temporary = _write_owned_temporary(json_path, json_data)
    markdown_temporary = _write_owned_temporary(markdown_path, markdown_data)
    installed_markdown = False
    try:
        os.link(markdown_temporary, markdown_path)
        installed_markdown = True
        os.link(json_temporary, json_path)
    except BaseException:
        if installed_markdown:
            markdown_path.unlink(missing_ok=True)
        raise
    finally:
        json_temporary.unlink(missing_ok=True)
        markdown_temporary.unlink(missing_ok=True)
    if json_path.read_bytes() != json_data or markdown_path.read_bytes() != markdown_data:
        raise OSError("installed acceptance report bytes differ")


def _add_evidence_arguments(parser: argparse.ArgumentParser) -> None:
    for prefix in ("legacy", "factorized"):
        parser.add_argument(f"--{prefix}-manifest", type=Path, required=True)
        parser.add_argument(f"--{prefix}-scores", type=Path, required=True)
        parser.add_argument(f"--{prefix}-runs", type=Path, required=True)
        parser.add_argument(f"--{prefix}-failures", type=Path, required=True)
    parser.add_argument("--resource-manifest", type=Path, required=True)
    parser.add_argument("--resource-runs", type=Path, required=True)
    parser.add_argument("--resource-failures", type=Path, required=True)
    parser.add_argument("--dataset-metadata", type=Path, required=True)
    parser.add_argument(
        "--chainsaw-expected-success-ids",
        type=Path,
        required=True,
    )
    parser.add_argument("--model-manifest", type=Path, required=True)
    parser.add_argument("--runtime-manifest", type=Path, required=True)
    parser.add_argument("--eligibility-manifest", type=Path, required=True)
    parser.add_argument("--standalone-baseline", type=Path, required=True)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    coverage = subparsers.add_parser("coverage")
    _add_evidence_arguments(coverage)
    coverage.add_argument("--coverage-out", type=Path, required=True)

    evaluate = subparsers.add_parser("evaluate")
    _add_evidence_arguments(evaluate)
    evaluate.add_argument("--coverage-attestation", type=Path, required=True)
    evaluate.add_argument("--json-out", type=Path, required=True)
    evaluate.add_argument("--markdown-out", type=Path, required=True)
    evaluate.add_argument("--bootstrap-replicates", type=int, required=True)
    evaluate.add_argument("--seed", type=int, required=True)

    promotion = subparsers.add_parser("validate-promotion")
    _add_evidence_arguments(promotion)
    promotion.add_argument("--coverage-attestation", type=Path, required=True)
    promotion.add_argument("--acceptance", type=Path, required=True)
    promotion.add_argument("--markdown", type=Path, required=True)
    promotion.add_argument("--task17-commit", required=True)
    return parser.parse_args(argv)


def _coverage_command(args: argparse.Namespace) -> int:
    inputs = _validate_input_paths(_evidence_path_arguments(args))
    _validate_absent_outputs({"coverage_out": args.coverage_out}, inputs)
    payload = _build_coverage_attestation(args)
    _validate_coverage_payload(payload)
    _write_canonical_absent(args.coverage_out, payload)
    return 0


def _evaluate_command(args: argparse.Namespace) -> int:
    if args.bootstrap_replicates != BOOTSTRAP["replicates"] or args.seed != BOOTSTRAP["seed"]:
        raise InvalidEvidence("evaluation bootstrap arguments differ from the frozen contract")
    inputs = _validate_input_paths(
        {
            **_evidence_path_arguments(args),
            "coverage_attestation": args.coverage_attestation,
        }
    )
    _validate_absent_outputs(
        {"json_out": args.json_out, "markdown_out": args.markdown_out},
        inputs,
    )
    coverage = _load_revalidated_coverage(args)
    payload = _build_acceptance_payload(args, coverage)
    write_acceptance_report(args.json_out, args.markdown_out, payload)
    return 0 if payload["all_gates_pass"] is True else 1


def _git_output(repo_root: Path, *arguments: str) -> bytes:
    completed = subprocess.run(
        ["git", "-C", os.fspath(repo_root), *arguments],
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise InvalidEvidence("Task 17 Git commit validation failed")
    return completed.stdout


def _verify_commit_binding(args: argparse.Namespace) -> None:
    commit = args.task17_commit
    if (
        not isinstance(commit, str)
        or len(commit) != 40
        or any(character not in "0123456789abcdef" for character in commit)
    ):
        raise InvalidEvidence("Task 17 commit must be a full lowercase Git SHA-1")
    repo_root = Path(
        _git_output(Path.cwd(), "rev-parse", "--show-toplevel").decode("utf-8").strip()
    ).resolve(strict=True)
    _git_output(repo_root, "cat-file", "-e", f"{commit}^{{commit}}")
    completed = subprocess.run(
        ["git", "-C", os.fspath(repo_root), "merge-base", "--is-ancestor", commit, "HEAD"],
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise InvalidEvidence("Task 17 commit is not an ancestor of HEAD")

    bound_paths = [
        *_evidence_path_arguments(args).values(),
        Path(args.coverage_attestation),
        Path(args.acceptance),
        Path(args.markdown),
        repo_root / "benchmark/models/factorized_ranker_v1_execution_receipt.json",
        repo_root / "benchmark/REPORT.md",
    ]
    # The model/runtime/baseline/dataset/expected-ID authority inputs are not
    # Task 17 evidence blobs. Remove them from the sixteen-file commit set.
    excluded_roles = {
        "dataset_metadata",
        "chainsaw_expected_success_ids",
        "model_manifest",
        "runtime_manifest",
        "eligibility_manifest",
        "standalone_baseline",
    }
    excluded = {
        Path(getattr(args, role)).resolve(strict=True) for role in excluded_roles
    }
    selected: list[Path] = []
    for raw_path in bound_paths:
        path = _require_regular_file(raw_path, "Task 17 committed evidence")
        if path in excluded:
            continue
        try:
            path.relative_to(repo_root)
        except ValueError as error:
            raise InvalidEvidence("Task 17 evidence lies outside the repository") from error
        if path not in selected:
            selected.append(path)
    if len(selected) != 16:
        raise InvalidEvidence("Task 17 commit binding does not contain exactly sixteen paths")
    for path in selected:
        relative = path.relative_to(repo_root).as_posix()
        committed = _git_output(repo_root, "show", f"{commit}:{relative}")
        if committed != path.read_bytes():
            raise InvalidEvidence(f"Task 17 committed blob differs for {relative}")


def _promotion_command(args: argparse.Namespace) -> int:
    _validate_input_paths(
        {
            **_evidence_path_arguments(args),
            "coverage_attestation": args.coverage_attestation,
            "acceptance": args.acceptance,
            "markdown": args.markdown,
        }
    )
    coverage = _load_revalidated_coverage(args)
    expected = _build_acceptance_payload(args, coverage)
    try:
        observed = load_canonical_json(Path(args.acceptance))
    except (OSError, ValueError) as error:
        raise InvalidEvidence("acceptance JSON is not canonical") from error
    _validate_acceptance_payload(observed)
    if observed != expected:
        raise InvalidEvidence("acceptance JSON does not match rederived gates")
    if Path(args.markdown).read_bytes() != _render_acceptance_markdown(observed):
        raise InvalidEvidence("acceptance Markdown does not match canonical JSON")
    _verify_commit_binding(args)
    runtime = load_canonical_json(Path(args.runtime_manifest))
    if runtime.get("selector_cache_contract") != CACHE_CONTRACT:
        raise InvalidEvidence("runtime cache contract is invalid")

    reasons: list[str] = []
    if observed["all_gates_pass"] is not True:
        reasons.append("performance_gates_failed")
    if coverage["factorized_structural_coverage"] < 1.0:
        reasons.append("structural_coverage_incomplete")
    if runtime["selector_cache_contract"]["default_promotion_compatible"] is not True:
        reasons.append("cache_context_not_promotable")
    if reasons:
        print("KEEP_OPT_IN " + " ".join(reasons))
        return 1
    print("PROMOTE")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.command == "coverage":
            return _coverage_command(args)
        if args.command == "evaluate":
            return _evaluate_command(args)
        if args.command == "validate-promotion":
            return _promotion_command(args)
        raise InvalidEvidence("unknown acceptance command")
    except InvalidEvidence:
        raise
    except (OSError, TypeError, ValueError) as error:
        raise InvalidEvidence(str(error)) from error


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except InvalidEvidence as error:
        print(f"invalid evidence: {error}", file=sys.stderr)
        raise SystemExit(2) from error
