"""Validate and export frozen factorized models, Rust arrays, and goldens."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import platform
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import scipy
import sklearn

import benchmark.factorized_ranker.pairs as pair_api
from benchmark.factorized_ranker.corpus import (
    CANDIDATE_FIELDS,
    CHAIN_FIELDS,
    COUNT_FIELDS,
    NORMALIZED_REJECTION_FIELDS,
    feature_schema_hash,
)
from benchmark.factorized_ranker.folds import (
    _read_canonical_table,
    _validate_task8_manifest,
    load_fold_manifest,
)
from benchmark.factorized_ranker.model_artifact import (
    ARTIFACT_SCHEMA_VERSION,
    MODEL_INPUT_DTYPE,
    THRESHOLD_POLICY,
    _canonical_json_bytes,
    _matrix_hash,
    _sha256,
    golden_vectors,
    load_artifact,
    reference_pair_from_batch,
    render_rust_models,
    validate_artifacts,
    validate_golden_vectors,
)
from benchmark.factorized_ranker.pairs import build_candidate_pairs, build_count_pairs
from benchmark.factorized_ranker.training import (
    FEATURE_FAMILY_ORDER,
    MODEL_GRID,
    OOF_FIELDS,
    OOFRow,
    CorpusTables,
    Hyperparameters,
    _id_hash,
    _params_payload,
    _render_oof,
    _spec_payload,
    head_feature_spec,
    validate_retained_families,
)


MANIFEST_SCHEMA_VERSION = 1
MANIFEST_KEYS = {
    "schema_version",
    "seed",
    "input_dtype",
    "threshold_policy",
    "retained_feature_families",
    "training_command",
    "export_command",
    "source_git_commit",
    "feature_dump_binary_sha256",
    "dataset_sha256",
    "corpus_manifest_sha256",
    "corpus_table_sha256s",
    "fold_manifest_sha256",
    "cv_report_sha256",
    "ablation_report_sha256",
    "oof_predictions_sha256",
    "count_model_sha256",
    "candidate_model_sha256",
    "golden_sha256",
    "generated_rust_sha256",
    "standalone_baseline_sha256",
    "feature_schema_sha256",
    "count_pair_feature_names_sha256",
    "candidate_pair_feature_names_sha256",
    "count_hyperparameters",
    "candidate_hyperparameters",
    "tree_counts",
    "node_counts",
    "combined_tree_count",
    "combined_node_count",
    "versions",
}

_PATH_ROLES = {
    "--count-model": "<COUNT_MODEL>",
    "--candidate-model": "<CANDIDATE_MODEL>",
    "--corpus-manifest": "<CORPUS_MANIFEST>",
    "--fold-manifest": "<FOLD_MANIFEST>",
    "--cv-report": "<CV_REPORT>",
    "--ablation-report": "<ABLATION_REPORT>",
    "--oof-predictions": "<OOF_PREDICTIONS>",
    "--standalone-baseline": "<STANDALONE_BASELINE>",
    "--rust-out": "<RUST_OUT>",
    "--golden-out": "<GOLDEN_OUT>",
    "--manifest-out": "<MANIFEST_OUT>",
}
_TRAINING_PATH_ROLES = {
    "--corpus-dir": "<CORPUS_DIR>",
    "--fold-manifest": "<FOLD_MANIFEST>",
    "--out-dir": "<OUT_DIR>",
    "--count-model-out": "<COUNT_MODEL_OUT>",
    "--candidate-model-out": "<CANDIDATE_MODEL_OUT>",
}
_CV_KEYS = {
    "schema_version",
    "seed",
    "evidence_role",
    "normalized_command",
    "dataset",
    "dataset_sha256",
    "corpus_manifest_sha256",
    "chains_sha256",
    "fold_manifest_sha256",
    "feature_schema_version",
    "feature_schema_sha256",
    "binary_sha256",
    "git_commit",
    "accepted_chain_count",
    "accepted_chain_id_sha256",
    "feature_family_order",
    "feature_family_mapping",
    "model_grid",
    "versions",
    "grid_history",
    "final",
}
_ABLATION_KEYS = (_CV_KEYS - {"grid_history"}) | {"stages"}
_FINAL_KEYS = {
    "retained_families",
    "selected_hyperparameters",
    "feature_specs",
    "oof_filename",
    "oof_sha256",
    "oof_chain_count",
    "oof_cohorts",
    "final_fit",
}
_OOF_INTEGER_FIELDS = {"fold", "selected_count", "n_true_domains", "count_correct"}
_OOF_FLOAT_FIELDS = {
    "count_borda_score",
    "candidate_borda_score",
    "ndo",
    "boundary_f1_10",
    "matched_dice",
    "total_regret",
    "count_regret",
    "within_count_regret",
}


@dataclass(frozen=True)
class ExportResult:
    manifest: dict[str, object]
    generated_rust_sha256: str
    golden_sha256: str
    manifest_sha256: str


def normalize_export_argv(argv: Sequence[str]) -> list[str]:
    normalized: list[str] = []
    index = 0
    while index < len(argv):
        token = str(argv[index])
        option, separator, _value = token.partition("=")
        if option in _PATH_ROLES and separator:
            normalized.append(f"{option}={_PATH_ROLES[option]}")
            index += 1
        elif token in _PATH_ROLES:
            if index + 1 >= len(argv):
                raise ValueError(f"path option {token} has no value")
            normalized.extend((token, _PATH_ROLES[token]))
            index += 2
        else:
            normalized.append(token)
            index += 1
    return normalized


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    output: dict[str, object] = {}
    for key, value in pairs:
        if key in output:
            raise ValueError(f"duplicate JSON key {key!r}")
        output[key] = value
    return output


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value}")


def _load_canonical_json(path: Path, description: str) -> tuple[dict[str, object], bytes]:
    data = Path(path).read_bytes()
    try:
        value = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=_unique_object,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"{description} is not valid canonical JSON") from error
    if not isinstance(value, dict):
        raise ValueError(f"{description} root must be an object")
    if _canonical_json_bytes(value) != data:
        raise ValueError(f"{description} is not canonically serialized")
    return value, data


def _integer(value: object, description: str, *, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{description} must be an integer >= {minimum}")
    return value


def _finite(value: object, description: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{description} must be a finite number")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{description} must be finite")
    return number


def _text(value: object, description: str) -> str:
    if not isinstance(value, str) or not value or "\0" in value:
        raise ValueError(f"{description} must be nonempty text without NUL")
    return value


def _hash(value: object, description: str) -> str:
    text = _text(value, description)
    if len(text) != 64 or any(character not in "0123456789abcdef" for character in text):
        raise ValueError(f"{description} must be a lowercase SHA-256")
    return text


def _command(value: object, description: str) -> tuple[str, ...]:
    if not isinstance(value, list) or not value:
        raise ValueError(f"{description} must be a nonempty token array")
    tokens = tuple(_text(token, description) for token in value)
    if any(token.startswith("/") or "=/" in token for token in tokens):
        raise ValueError(f"{description} contains a physical path")
    return tokens


def _validate_semantic_command(
    value: object,
    *,
    executable: str,
    path_roles: Mapping[str, str],
    allow_seed: bool,
    allow_jobs: bool = False,
    allowed_flags: frozenset[str] = frozenset(),
    description: str,
) -> tuple[str, ...]:
    tokens = _command(value, description)
    if tokens[0] != executable:
        raise ValueError(f"{description} executable mismatch")
    seen: set[str] = set()
    index = 1
    while index < len(tokens):
        token = tokens[index]
        option, separator, inline = token.partition("=")
        if option in path_roles:
            if option in seen:
                raise ValueError(f"{description} repeats {option}")
            if separator:
                observed = inline
                index += 1
            else:
                if index + 1 >= len(tokens):
                    raise ValueError(f"{description} path option {option} has no value")
                observed = tokens[index + 1]
                index += 2
            if observed != path_roles[option]:
                raise ValueError(f"{description} path role {option} mismatch")
            seen.add(option)
        elif allow_seed and option == "--seed":
            if option in seen:
                raise ValueError(f"{description} repeats --seed")
            if separator:
                observed = inline
                index += 1
            else:
                if index + 1 >= len(tokens):
                    raise ValueError(f"{description} --seed has no value")
                observed = tokens[index + 1]
                index += 2
            if observed != "37":
                raise ValueError(f"{description} seed is not 37")
            seen.add(option)
        elif allow_jobs and option == "--jobs":
            if option in seen:
                raise ValueError(f"{description} repeats --jobs")
            if separator:
                observed = inline
                index += 1
            else:
                if index + 1 >= len(tokens):
                    raise ValueError(f"{description} --jobs has no value")
                observed = tokens[index + 1]
                index += 2
            if observed not in ("1", "2", "3", "4", "5", "6", "7", "8"):
                raise ValueError(f"{description} jobs is not in 1..=8")
            seen.add(option)
        elif option in allowed_flags:
            if separator or option in seen:
                raise ValueError(f"{description} has invalid or repeated {option}")
            seen.add(option)
            index += 1
        else:
            raise ValueError(f"{description} contains unknown token {token!r}")
    if not set(path_roles).issubset(seen):
        raise ValueError(f"{description} omits a required semantic path role")
    return tokens


def _no_physical_or_self_data(value: object) -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            if str(key) in {"manifest_sha256", "self_sha256"}:
                raise ValueError("export manifest contains a self-hash field")
            _no_physical_or_self_data(child)
    elif isinstance(value, list):
        for child in value:
            _no_physical_or_self_data(child)
    elif isinstance(value, str) and (value.startswith("/") or "=/" in value):
        raise ValueError("export manifest contains a physical path")


def _manifest_hyperparameters(value: object, description: str) -> Hyperparameters:
    if not isinstance(value, Mapping) or set(value) != {
        "n_estimators",
        "learning_rate",
        "min_samples_leaf",
        "max_depth",
        "loss",
        "random_state",
    }:
        raise ValueError(f"{description} hyperparameter schema mismatch")
    params = Hyperparameters(
        _integer(value["n_estimators"], f"{description} n_estimators", minimum=1),
        _finite(value["learning_rate"], f"{description} learning_rate"),
        _integer(value["min_samples_leaf"], f"{description} min_samples_leaf", minimum=1),
        _integer(value["max_depth"], f"{description} max_depth", minimum=1),
    )
    if params not in MODEL_GRID or value["loss"] != "log_loss" or value["random_state"] != 37:
        raise ValueError(f"{description} hyperparameters are outside the frozen contract")
    return params


def validate_export_manifest(manifest: Mapping[str, object]) -> None:
    if not isinstance(manifest, Mapping) or set(manifest) != MANIFEST_KEYS:
        raise ValueError("export manifest keys do not match schema version 1")
    if _integer(manifest["schema_version"], "manifest schema version", minimum=1) != MANIFEST_SCHEMA_VERSION:
        raise ValueError("export manifest schema version is unsupported")
    if _integer(manifest["seed"], "manifest seed") != 37:
        raise ValueError("export manifest seed must be 37")
    if manifest["input_dtype"] != MODEL_INPUT_DTYPE or manifest["threshold_policy"] != THRESHOLD_POLICY:
        raise ValueError("export manifest numeric policy mismatch")
    validate_retained_families(manifest["retained_feature_families"])  # type: ignore[arg-type]
    _validate_semantic_command(
        manifest["training_command"],
        executable="benchmark.train_factorized_ranker",
        path_roles=_TRAINING_PATH_ROLES,
        allow_seed=True,
        allow_jobs=True,
        allowed_flags=frozenset({"--resume"}),
        description="manifest training command",
    )
    _validate_semantic_command(
        manifest["export_command"],
        executable="benchmark.export_factorized_ranker",
        path_roles=_PATH_ROLES,
        allow_seed=False,
        description="manifest export command",
    )
    source_commit = _text(manifest["source_git_commit"], "source Git commit")
    if len(source_commit) != 40 or any(character not in "0123456789abcdef" for character in source_commit):
        raise ValueError("source Git commit is not lowercase 40-hex")
    for field in (
        "feature_dump_binary_sha256",
        "dataset_sha256",
        "corpus_manifest_sha256",
        "fold_manifest_sha256",
        "cv_report_sha256",
        "ablation_report_sha256",
        "oof_predictions_sha256",
        "count_model_sha256",
        "candidate_model_sha256",
        "golden_sha256",
        "generated_rust_sha256",
        "standalone_baseline_sha256",
        "feature_schema_sha256",
        "count_pair_feature_names_sha256",
        "candidate_pair_feature_names_sha256",
    ):
        _hash(manifest[field], field)
    table_hashes = manifest["corpus_table_sha256s"]
    if not isinstance(table_hashes, Mapping) or set(table_hashes) != {
        "chains",
        "counts",
        "candidates",
        "rejections",
    }:
        raise ValueError("manifest corpus table hash schema mismatch")
    for name, value in table_hashes.items():
        _hash(value, f"{name} table hash")
    _manifest_hyperparameters(manifest["count_hyperparameters"], "count")
    _manifest_hyperparameters(manifest["candidate_hyperparameters"], "candidate")
    tree_counts = manifest["tree_counts"]
    node_counts = manifest["node_counts"]
    if not isinstance(tree_counts, Mapping) or set(tree_counts) != {"count", "candidate"}:
        raise ValueError("manifest tree counts schema mismatch")
    if not isinstance(node_counts, Mapping) or set(node_counts) != {"count", "candidate"}:
        raise ValueError("manifest node counts schema mismatch")
    count_trees = _integer(tree_counts["count"], "count tree count", minimum=1)
    candidate_trees = _integer(tree_counts["candidate"], "candidate tree count", minimum=1)
    count_nodes = _integer(node_counts["count"], "count node count", minimum=1)
    candidate_nodes = _integer(node_counts["candidate"], "candidate node count", minimum=1)
    if manifest["combined_tree_count"] != count_trees + candidate_trees:
        raise ValueError("manifest combined tree count mismatch")
    if manifest["combined_node_count"] != count_nodes + candidate_nodes:
        raise ValueError("manifest combined node count mismatch")
    versions = manifest["versions"]
    if not isinstance(versions, Mapping) or set(versions) != {
        "python",
        "numpy",
        "pandas",
        "scipy",
        "scikit_learn",
        "rustfmt",
    }:
        raise ValueError("manifest versions schema mismatch")
    for name, value in versions.items():
        _text(value, f"{name} version")
    _no_physical_or_self_data(manifest)


def _load_corpus(
    manifest_path: Path,
    expected_manifest_sha256: str,
    count_artifact: Mapping[str, object],
) -> tuple[dict[str, object], CorpusTables, dict[str, str], set[str]]:
    payload, data = _load_canonical_json(manifest_path, "corpus manifest")
    if _sha256(data) != expected_manifest_sha256:
        raise ValueError("corpus manifest hash disagrees with model artifacts")
    manifest = _validate_task8_manifest(payload, "cath17287")
    if manifest["binary_sha256"] != count_artifact["feature_dump_binary_sha256"]:
        raise ValueError("corpus feature-dump binary hash mismatch")
    if manifest["git_commit"] != count_artifact["source_git_commit"]:
        raise ValueError("corpus source Git commit mismatch")
    if manifest["feature_schema_hash"] != count_artifact["feature_schema_sha256"]:
        raise ValueError("corpus feature schema hash mismatch")
    root = Path(manifest_path).resolve(strict=True).parent
    fields = {
        "chains": CHAIN_FIELDS,
        "counts": COUNT_FIELDS,
        "candidates": CANDIDATE_FIELDS,
        "rejections": NORMALIZED_REJECTION_FIELDS,
    }
    frames: dict[str, pd.DataFrame] = {}
    hashes: dict[str, str] = {}
    for name, expected_fields in fields.items():
        path = (root / f"{name}.csv").resolve(strict=True)
        if path.parent != root:
            raise ValueError("corpus table escapes the manifest directory")
        table_data, rows = _read_canonical_table(path, expected_fields)
        descriptor = manifest["tables"][name]  # type: ignore[index]
        if _sha256(table_data) != descriptor["sha256"] or len(rows) != descriptor["rows"]:  # type: ignore[index]
            raise ValueError(f"corpus table {name} hash/row count mismatch")
        hashes[name] = _sha256(table_data)
        # Match Task 11's authoritative pandas CSV conversion exactly. The
        # canonical text rows above validate bytes/schema; rebuilding floats
        # through Python's scalar parser can differ by one ulp from pandas.
        frames[name] = pd.read_csv(path, keep_default_na=False)
        if tuple(frames[name].columns) != tuple(expected_fields):
            raise ValueError(f"corpus table {name} pandas header mismatch")
    corpus = CorpusTables(frames["chains"], frames["counts"], frames["candidates"])
    chain_rows = pair_api._validate_chains(corpus.chains)
    count_rows = pair_api._validate_counts(corpus.counts, chain_rows)
    candidate_rows = pair_api._validate_candidates(corpus.candidates, chain_rows)
    count_groups = {
        (chain_id, row.count)
        for chain_id, rows in count_rows.items()
        for row in rows
    }
    candidate_groups = {
        (chain_id, row.count)
        for chain_id, rows in candidate_rows.items()
        for row in rows
    }
    if not candidate_groups.issubset(count_groups):
        raise ValueError("candidate count groups are not covered by count rows")
    if manifest["accepted_unique_chains"] != len(chain_rows):
        raise ValueError("corpus accepted-chain count mismatch")
    return manifest, corpus, hashes, set(chain_rows)


def _validate_reports(
    cv_path: Path,
    ablation_path: Path,
    count_artifact: Mapping[str, object],
    candidate_artifact: Mapping[str, object],
    corpus_manifest: Mapping[str, object],
    corpus_manifest_sha256: str,
    fold_manifest_sha256: str,
    chains_sha256: str,
    accepted_ids: set[str],
    corpus: CorpusTables,
) -> tuple[dict[str, object], dict[str, object], PairBatchBundle]:
    cv, cv_bytes = _load_canonical_json(cv_path, "CV report")
    ablation, ablation_bytes = _load_canonical_json(ablation_path, "ablation report")
    if set(cv) != _CV_KEYS or set(ablation) != _ABLATION_KEYS:
        raise ValueError("training report schema mismatch")
    if _sha256(cv_bytes) != count_artifact["cv_report_sha256"]:
        raise ValueError("CV report hash mismatch")
    if _sha256(ablation_bytes) != count_artifact["ablation_report_sha256"]:
        raise ValueError("ablation report hash mismatch")
    common_fields = _CV_KEYS - {"grid_history", "final"}
    for field in common_fields:
        if cv[field] != ablation[field]:
            raise ValueError(f"CV/ablation common field {field} disagrees")
    if cv["final"] != ablation["final"]:
        raise ValueError("CV/ablation final decisions disagree")
    if cv["schema_version"] != 1 or cv["seed"] != 37:
        raise ValueError("training report schema/seed mismatch")
    expected_source = {
        "dataset": "cath17287",
        "dataset_sha256": corpus_manifest["dataset_sha256"],
        "corpus_manifest_sha256": corpus_manifest_sha256,
        "chains_sha256": chains_sha256,
        "fold_manifest_sha256": fold_manifest_sha256,
        "feature_schema_version": 1,
        "feature_schema_sha256": feature_schema_hash(),
        "binary_sha256": corpus_manifest["binary_sha256"],
        "git_commit": corpus_manifest["git_commit"],
        "accepted_chain_count": len(accepted_ids),
        "accepted_chain_id_sha256": _id_hash(sorted(accepted_ids)),
        "feature_family_order": list(FEATURE_FAMILY_ORDER),
        "model_grid": [_params_payload(params) for params in MODEL_GRID],
    }
    for field, expected in expected_source.items():
        if cv[field] != expected:
            raise ValueError(f"training report source field {field} mismatch")
    if cv["evidence_role"] != "development_model_selection_not_locked_acceptance":
        raise ValueError("training report evidence role mismatch")
    if cv["normalized_command"] != count_artifact["training_command"]:
        raise ValueError("training report command mismatch")
    versions = cv["versions"]
    expected_versions = {
        "python_implementation": platform.python_implementation(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
    }
    if versions != expected_versions:
        raise ValueError("training report versions mismatch current exporter")
    artifact_versions = {
        "python": count_artifact["python_version"],
        "numpy": count_artifact["numpy_version"],
        "pandas": count_artifact["pandas_version"],
        "scipy": count_artifact["scipy_version"],
        "scikit_learn": count_artifact["sklearn_version"],
    }
    if {name: versions[name] for name in artifact_versions} != artifact_versions:  # type: ignore[index]
        raise ValueError("model/report versions disagree")
    expected_family_mapping = {
        family: {
            head: _spec_payload(
                head_feature_spec(
                    head,
                    tuple(
                        candidate
                        for candidate in FEATURE_FAMILY_ORDER
                        if candidate == "base"
                        or FEATURE_FAMILY_ORDER.index(candidate)
                        <= FEATURE_FAMILY_ORDER.index(family)
                    ),
                )
            )
            for head in ("count", "candidate")
        }
        for family in FEATURE_FAMILY_ORDER
    }
    if cv["feature_family_mapping"] != expected_family_mapping:
        raise ValueError("training report feature-family mapping mismatch")
    grid_history = cv["grid_history"]
    if (
        not isinstance(grid_history, list)
        or len(grid_history) != len(FEATURE_FAMILY_ORDER)
        or [
            entry.get("family") if isinstance(entry, Mapping) else None
            for entry in grid_history
        ]
        != list(FEATURE_FAMILY_ORDER)
        or any(
            not isinstance(entry, Mapping)
            or set(entry) != {"family", "count", "candidate"}
            for entry in grid_history
        )
    ):
        raise ValueError("CV grid history is malformed or out of order")
    stages = ablation["stages"]
    if (
        not isinstance(stages, list)
        or len(stages) != len(FEATURE_FAMILY_ORDER)
        or [stage.get("family") if isinstance(stage, Mapping) else None for stage in stages]
        != list(FEATURE_FAMILY_ORDER)
        or any(
            not isinstance(stage, Mapping)
            or not {"family", "retained", "oof_sha256"}.issubset(stage)
            or type(stage["retained"]) is not bool
            for stage in stages
        )
    ):
        raise ValueError("ablation stages are missing")
    final = cv["final"]
    if not isinstance(final, Mapping) or set(final) != _FINAL_KEYS:
        raise ValueError("training final decision schema mismatch")
    retained = validate_retained_families(final["retained_families"])  # type: ignore[arg-type]
    if list(retained) != count_artifact["retained_feature_families"]:
        raise ValueError("training retained families mismatch artifacts")
    selected = final["selected_hyperparameters"]
    specs = final["feature_specs"]
    final_fit = final["final_fit"]
    if not isinstance(selected, Mapping) or set(selected) != {"count", "candidate"}:
        raise ValueError("training selected hyperparameters schema mismatch")
    if not isinstance(specs, Mapping) or set(specs) != {"count", "candidate"}:
        raise ValueError("training feature specs schema mismatch")
    if not isinstance(final_fit, Mapping) or set(final_fit) != {"count", "candidate"}:
        raise ValueError("training final-fit schema mismatch")
    batches: dict[str, Any] = {}
    for head, artifact in (("count", count_artifact), ("candidate", candidate_artifact)):
        params = Hyperparameters(
            int(artifact["n_estimators"]),
            float(artifact["learning_rate"]),
            int(artifact["min_samples_leaf"]),
            int(artifact["max_depth"]),
        )
        if selected[head] != _params_payload(params):  # type: ignore[index]
            raise ValueError(f"training {head} selected parameters mismatch artifact")
        spec = head_feature_spec(head, retained)  # type: ignore[arg-type]
        if specs[head] != _spec_payload(spec):  # type: ignore[index]
            raise ValueError(f"training {head} feature spec mismatch")
        if head == "count":
            batch = build_count_pairs(
                corpus.chains,
                corpus.counts,
                shared_features=spec.shared_features,
                item_features=spec.item_features,
                seed=37,
            )
        else:
            batch = build_candidate_pairs(
                corpus.chains,
                corpus.candidates,
                shared_features=spec.shared_features,
                item_features=spec.item_features,
                seed=37,
            )
        expected_fit = {
            "rows": len(batch.y),
            "features": batch.x.shape[1],
            "classes": [0, 1],
        }
        if final_fit[head] != expected_fit:  # type: ignore[index]
            raise ValueError(f"training {head} final-fit summary mismatch")
        audit = artifact["parity_audit"]
        if audit["row_count"] != len(batch.x) or audit["column_count"] != batch.x.shape[1]:  # type: ignore[index]
            raise ValueError(f"{head} artifact parity shape mismatch final pair matrix")
        if audit["matrix_sha256"] != _matrix_hash(batch.x, batch.feature_names):  # type: ignore[index]
            raise ValueError(f"{head} artifact parity matrix hash mismatch")
        if artifact["reference_pair"] != reference_pair_from_batch(
            corpus, batch, head, retained  # type: ignore[arg-type]
        ):
            raise ValueError(f"{head} artifact reference pair mismatch")
        batches[head] = batch
    if final["oof_filename"] != "oof_predictions.csv":
        raise ValueError("training OOF filename mismatch")
    if final["oof_sha256"] != count_artifact["oof_predictions_sha256"]:
        raise ValueError("training OOF hash mismatch artifacts")
    if final["oof_chain_count"] != len(accepted_ids):
        raise ValueError("training OOF chain count mismatch")
    retained_stages = [
        stage
        for stage in stages
        if isinstance(stage, Mapping) and stage.get("retained") is True
    ]
    if [stage["family"] for stage in retained_stages] != list(retained):
        raise ValueError("ablation retained-stage sequence disagrees with final families")
    if not retained_stages or retained_stages[-1].get("oof_sha256") != final["oof_sha256"]:
        raise ValueError("last retained ablation stage does not bind final OOF")
    return cv, ablation, PairBatchBundle(batches["count"], batches["candidate"])


@dataclass(frozen=True)
class PairBatchBundle:
    count: Any
    candidate: Any


def _load_oof(
    path: Path,
    expected_hash: str,
    accepted_ids: set[str],
    fold_assignments: Mapping[str, object],
    corpus: CorpusTables,
) -> tuple[tuple[OOFRow, ...], bytes]:
    data = Path(path).read_bytes()
    if _sha256(data) != expected_hash:
        raise ValueError("OOF predictions hash mismatch")
    if not data.endswith(b"\n") or b"\r" in data:
        raise ValueError("OOF predictions use noncanonical newlines")
    try:
        reader = csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""))
    except UnicodeDecodeError as error:
        raise ValueError("OOF predictions are not UTF-8") from error
    if tuple(reader.fieldnames or ()) != OOF_FIELDS:
        raise ValueError("OOF predictions header mismatch")
    rows: list[OOFRow] = []
    for raw in reader:
        if set(raw) != set(OOF_FIELDS) or any(value is None for value in raw.values()):
            raise ValueError("OOF row schema mismatch")
        values: dict[str, object] = {}
        for field in OOF_FIELDS:
            text = raw[field]
            if field in _OOF_INTEGER_FIELDS:
                try:
                    number = int(text)
                except ValueError as error:
                    raise ValueError(f"OOF integer field {field} is malformed") from error
                if str(number) != text:
                    raise ValueError(f"OOF integer field {field} is noncanonical")
                values[field] = number
            elif field in _OOF_FLOAT_FIELDS:
                try:
                    number = float(text)
                except ValueError as error:
                    raise ValueError(f"OOF float field {field} is malformed") from error
                if not math.isfinite(number) or format(number, ".17g") != text:
                    raise ValueError(f"OOF float field {field} is noncanonical")
                values[field] = number
            else:
                values[field] = _text(text, f"OOF {field}")
        rows.append(OOFRow(**values))  # type: ignore[arg-type]
    if _render_oof(rows) != data:
        raise ValueError("OOF predictions are not canonically ordered/serialized")
    if len(rows) != len(accepted_ids) or {row.chain_id for row in rows} != accepted_ids:
        raise ValueError("OOF predictions do not cover exact accepted chains")
    chains = pair_api._validate_chains(corpus.chains)
    candidates = pair_api._validate_candidates(corpus.candidates, chains)
    for row in rows:
        assignment = fold_assignments.get(row.chain_id)
        if assignment is None or row.fold != assignment.fold:  # type: ignore[union-attr]
            raise ValueError("OOF fold disagrees with fold manifest")
        if row.true_count_bin != assignment.true_count_bin or row.length_bin != assignment.length_bin:  # type: ignore[union-attr]
            raise ValueError("OOF strata disagree with fold manifest")
        if row.n_true_domains != chains[row.chain_id].n_true_domains:
            raise ValueError("OOF truth count disagrees with chains table")
        if row.count_correct != int(row.selected_count == row.n_true_domains):
            raise ValueError("OOF count correctness is inconsistent")
        matches = [
            item
            for item in candidates[row.chain_id]
            if item.item_id == row.selected_candidate_id
            and item.canonical_delineation == row.canonical_delineation
            and item.count == row.selected_count
        ]
        if len(matches) != 1:
            raise ValueError("OOF winner identity does not match a normalized candidate")
        candidate_frame = corpus.candidates[corpus.candidates["chain_id"] == row.chain_id]
        winner_frame = candidate_frame[
            candidate_frame["candidate_id"] == row.selected_candidate_id
        ]
        if len(winner_frame) != 1:
            raise ValueError("OOF winner row is missing or duplicated")
        winner = winner_frame.iloc[0]
        for field in ("ndo", "boundary_f1_10", "matched_dice"):
            if getattr(row, field) != float(winner[field]):
                raise ValueError(f"OOF winner metric {field} was not copied exactly")
        best_all = float(candidate_frame["ndo"].max())
        selected_rows = candidate_frame[candidate_frame["num_domains"] == row.selected_count]
        if selected_rows.empty:
            raise ValueError("OOF selected count has no candidate rows")
        best_selected = float(selected_rows["ndo"].max())
        expected = (
            best_all - row.ndo,
            best_all - best_selected,
            best_selected - row.ndo,
        )
        observed = (row.total_regret, row.count_regret, row.within_count_regret)
        if any(value < 0.0 for value in observed) or any(
            abs(left - right) > 1e-12 for left, right in zip(observed, expected, strict=True)
        ):
            raise ValueError("OOF regret decomposition mismatch")
        if row.continuity_cohort not in {"contiguous", "discontinuous"}:
            raise ValueError("OOF continuity cohort is invalid")
        if row.label_cohort not in {"seen", "unseen", "unknown"}:
            raise ValueError("OOF label cohort is invalid")
    return tuple(rows), data


def _rustfmt_version() -> str:
    result = subprocess.run(
        ["rustfmt", "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode != 0 or not result.stdout.strip():
        raise ValueError("rustfmt version could not be determined")
    return result.stdout.strip()


def _validate_generated_rust(source: bytes) -> None:
    try:
        text = source.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("generated Rust is not UTF-8") from error
    with tempfile.TemporaryDirectory(prefix="factorized-export-rust-") as raw_root:
        root = Path(raw_root)
        generated = root / "generated.rs"
        generated.write_bytes(source)
        check = subprocess.run(
            ["rustfmt", "--edition", "2021", "--check", str(generated)],
            capture_output=True,
            text=True,
            check=False,
        )
        if check.returncode != 0:
            raise ValueError(f"rustfmt rejected generated Rust: {check.stderr.strip()}")
        (root / "model.rs").write_text(
            """
pub(crate) struct StaticNode { pub feature: u16, pub threshold: f64, pub left: u16, pub right: u16, pub leaf_value: f64 }
pub(crate) struct StaticTree { pub root: u16 }
pub(crate) struct StaticBoostedModel { pub schema_version: u32, pub feature_names: &'static [&'static str], pub initial_log_odds: f64, pub learning_rate: f64, pub trees: &'static [StaticTree], pub nodes: &'static [StaticNode] }
""".lstrip(),
            encoding="utf-8",
        )
        (root / "lib.rs").write_text(
            '#[path = "model.rs"] mod model;\n#[path = "generated.rs"] mod generated;\n',
            encoding="utf-8",
        )
        compile_result = subprocess.run(
            ["rustc", "--edition", "2021", "--crate-type", "lib", "lib.rs"],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
        if compile_result.returncode != 0:
            raise ValueError(
                f"generated Rust stub compilation failed: {compile_result.stderr.strip()}"
            )
    if not text.endswith("\n"):
        raise ValueError("generated Rust lacks a trailing newline")


def _model_hyperparameters(artifact: Mapping[str, object]) -> dict[str, object]:
    return {
        "n_estimators": artifact["n_estimators"],
        "learning_rate": artifact["learning_rate"],
        "min_samples_leaf": artifact["min_samples_leaf"],
        "max_depth": artifact["max_depth"],
        "loss": artifact["loss"],
        "random_state": artifact["random_state"],
    }


def _validate_output_targets(inputs: Sequence[Path], outputs: Sequence[Path]) -> None:
    input_paths = {Path(path).resolve(strict=True) for path in inputs}
    output_paths = [Path(path).resolve(strict=False) for path in outputs]
    if len(set(output_paths)) != len(output_paths):
        raise ValueError("export output paths alias each other")
    if any(path in input_paths for path in output_paths):
        raise ValueError("export output path aliases an input")


def _stage(path: Path, data: bytes) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.tmp-{os.getpid()}")
    with temporary.open("xb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    if temporary.read_bytes() != data:
        temporary.unlink(missing_ok=True)
        raise ValueError("staged export bytes changed after write")
    return temporary


def export_factorized_ranker(
    *,
    count_model: Path,
    candidate_model: Path,
    corpus_manifest: Path,
    fold_manifest: Path,
    cv_report: Path,
    ablation_report: Path,
    oof_predictions: Path,
    standalone_baseline: Path,
    rust_out: Path,
    golden_out: Path,
    manifest_out: Path,
    normalized_command: Sequence[str],
) -> ExportResult:
    command = list(normalized_command)
    _validate_semantic_command(
        command,
        executable="benchmark.export_factorized_ranker",
        path_roles=_PATH_ROLES,
        allow_seed=False,
        description="normalized export command",
    )
    input_paths = [
        Path(count_model),
        Path(candidate_model),
        Path(corpus_manifest),
        Path(fold_manifest),
        Path(cv_report),
        Path(ablation_report),
        Path(oof_predictions),
        Path(standalone_baseline),
    ]
    output_paths = [Path(rust_out), Path(golden_out), Path(manifest_out)]
    _validate_output_targets(input_paths, output_paths)
    count_bytes = Path(count_model).read_bytes()
    candidate_bytes = Path(candidate_model).read_bytes()
    count_artifact = load_artifact(count_model, expected_head="count")
    candidate_artifact = load_artifact(candidate_model, expected_head="candidate")
    validate_artifacts(count_artifact, candidate_artifact)
    corpus_hash = str(count_artifact["corpus_manifest_sha256"])
    manifest_source, corpus, table_hashes, accepted_ids = _load_corpus(
        corpus_manifest, corpus_hash, count_artifact
    )
    fold_hash = _sha256(Path(fold_manifest).read_bytes())
    if fold_hash != count_artifact["fold_manifest_sha256"]:
        raise ValueError("fold manifest hash disagrees with model artifacts")
    folds = load_fold_manifest(
        fold_manifest,
        expected_sha256=fold_hash,
        expected_dataset_sha256=str(manifest_source["dataset_sha256"]),
        expected_corpus_sha256=corpus_hash,
        expected_chains_sha256=table_hashes["chains"],
    )
    if folds.seed != 37 or folds.n_folds != 5:
        raise ValueError("fold manifest seed/count mismatch")
    if {assignment.chain_id for assignment in folds.assignments} != accepted_ids:
        raise ValueError("fold/accepted chain populations disagree")
    cv, _ablation, _batches = _validate_reports(
        cv_report,
        ablation_report,
        count_artifact,
        candidate_artifact,
        manifest_source,
        corpus_hash,
        fold_hash,
        table_hashes["chains"],
        accepted_ids,
        corpus,
    )
    oof_rows, oof_bytes = _load_oof(
        oof_predictions,
        str(count_artifact["oof_predictions_sha256"]),
        accepted_ids,
        {assignment.chain_id: assignment for assignment in folds.assignments},
        corpus,
    )
    if len(oof_rows) != cv["final"]["oof_chain_count"]:  # type: ignore[index]
        raise ValueError("OOF rows disagree with report count")
    baseline_bytes = Path(standalone_baseline).read_bytes()
    baseline_hash = _sha256(baseline_bytes)
    rust_text = render_rust_models(count_artifact, candidate_artifact)
    rust_bytes = rust_text.encode("utf-8")
    golden = golden_vectors(count_artifact, candidate_artifact)
    validate_golden_vectors(golden, count_artifact, candidate_artifact)
    golden_bytes = _canonical_json_bytes(golden)
    rust_hash = _sha256(rust_bytes)
    golden_hash = _sha256(golden_bytes)
    rustfmt_version = _rustfmt_version()
    manifest: dict[str, object] = {
        "schema_version": MANIFEST_SCHEMA_VERSION,
        "seed": 37,
        "input_dtype": MODEL_INPUT_DTYPE,
        "threshold_policy": THRESHOLD_POLICY,
        "retained_feature_families": list(count_artifact["retained_feature_families"]),  # type: ignore[arg-type]
        "training_command": list(count_artifact["training_command"]),  # type: ignore[arg-type]
        "export_command": command,
        "source_git_commit": count_artifact["source_git_commit"],
        "feature_dump_binary_sha256": count_artifact["feature_dump_binary_sha256"],
        "dataset_sha256": manifest_source["dataset_sha256"],
        "corpus_manifest_sha256": corpus_hash,
        "corpus_table_sha256s": table_hashes,
        "fold_manifest_sha256": fold_hash,
        "cv_report_sha256": count_artifact["cv_report_sha256"],
        "ablation_report_sha256": count_artifact["ablation_report_sha256"],
        "oof_predictions_sha256": _sha256(oof_bytes),
        "count_model_sha256": _sha256(count_bytes),
        "candidate_model_sha256": _sha256(candidate_bytes),
        "golden_sha256": golden_hash,
        "generated_rust_sha256": rust_hash,
        "standalone_baseline_sha256": baseline_hash,
        "feature_schema_sha256": count_artifact["feature_schema_sha256"],
        "count_pair_feature_names_sha256": count_artifact["feature_names_sha256"],
        "candidate_pair_feature_names_sha256": candidate_artifact["feature_names_sha256"],
        "count_hyperparameters": _model_hyperparameters(count_artifact),
        "candidate_hyperparameters": _model_hyperparameters(candidate_artifact),
        "tree_counts": {
            "count": count_artifact["tree_count"],
            "candidate": candidate_artifact["tree_count"],
        },
        "node_counts": {
            "count": count_artifact["node_count"],
            "candidate": candidate_artifact["node_count"],
        },
        "combined_tree_count": int(count_artifact["tree_count"])
        + int(candidate_artifact["tree_count"]),
        "combined_node_count": int(count_artifact["node_count"])
        + int(candidate_artifact["node_count"]),
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "rustfmt": rustfmt_version,
        },
    }
    validate_export_manifest(manifest)
    manifest_bytes = _canonical_json_bytes(manifest)
    # Revalidate every derived byte string before any final target is replaced.
    _validate_generated_rust(rust_bytes)
    parsed_golden = json.loads(golden_bytes)
    validate_golden_vectors(parsed_golden, count_artifact, candidate_artifact)
    parsed_manifest = json.loads(manifest_bytes)
    validate_export_manifest(parsed_manifest)
    if parsed_manifest["generated_rust_sha256"] != _sha256(rust_bytes):
        raise ValueError("manifest generated Rust hash mismatch")
    if parsed_manifest["golden_sha256"] != _sha256(golden_bytes):
        raise ValueError("manifest golden hash mismatch")

    temporary_paths: list[Path] = []
    try:
        rust_temp = _stage(rust_out, rust_bytes)
        temporary_paths.append(rust_temp)
        golden_temp = _stage(golden_out, golden_bytes)
        temporary_paths.append(golden_temp)
        manifest_temp = _stage(manifest_out, manifest_bytes)
        temporary_paths.append(manifest_temp)
        os.replace(rust_temp, rust_out)
        temporary_paths.remove(rust_temp)
        if Path(rust_out).read_bytes() != rust_bytes:
            raise OSError(f"installed generated Rust differs at {rust_out}")
        os.replace(golden_temp, golden_out)
        temporary_paths.remove(golden_temp)
        if Path(golden_out).read_bytes() != golden_bytes:
            raise OSError(f"installed golden differs at {golden_out}")
        os.replace(manifest_temp, manifest_out)
        temporary_paths.remove(manifest_temp)
        if Path(manifest_out).read_bytes() != manifest_bytes:
            raise OSError(f"installed manifest differs at {manifest_out}")
    except BaseException as error:
        for temporary in temporary_paths:
            temporary.unlink(missing_ok=True)
        raise OSError(
            f"export installation failed for {rust_out}, {golden_out}, and {manifest_out}"
        ) from error
    return ExportResult(manifest, rust_hash, golden_hash, _sha256(manifest_bytes))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count-model", type=Path, required=True)
    parser.add_argument("--candidate-model", type=Path, required=True)
    parser.add_argument("--corpus-manifest", type=Path, required=True)
    parser.add_argument("--fold-manifest", type=Path, required=True)
    parser.add_argument("--cv-report", type=Path, required=True)
    parser.add_argument("--ablation-report", type=Path, required=True)
    parser.add_argument("--oof-predictions", type=Path, required=True)
    parser.add_argument("--standalone-baseline", type=Path, required=True)
    parser.add_argument("--rust-out", type=Path, required=True)
    parser.add_argument("--golden-out", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    args = parser.parse_args(argv)
    raw = [
        "benchmark.export_factorized_ranker",
        *(sys.argv[1:] if argv is None else argv),
    ]
    try:
        result = export_factorized_ranker(
            count_model=args.count_model,
            candidate_model=args.candidate_model,
            corpus_manifest=args.corpus_manifest,
            fold_manifest=args.fold_manifest,
            cv_report=args.cv_report,
            ablation_report=args.ablation_report,
            oof_predictions=args.oof_predictions,
            standalone_baseline=args.standalone_baseline,
            rust_out=args.rust_out,
            golden_out=args.golden_out,
            manifest_out=args.manifest_out,
            normalized_command=normalize_export_argv(raw),
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(f"generated_rust_sha256 {result.generated_rust_sha256}")
    print(f"golden_sha256 {result.golden_sha256}")
    print(f"manifest_sha256 {result.manifest_sha256}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
