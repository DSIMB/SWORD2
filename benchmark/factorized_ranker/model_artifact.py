"""Validated frozen tree artifacts and deterministic Rust/golden rendering."""

from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import struct
import subprocess
import unicodedata
from dataclasses import dataclass
from numbers import Integral, Real
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.ensemble import GradientBoostingClassifier

import benchmark.factorized_ranker.pairs as pair_api
from benchmark.factorized_ranker.pairs import PairBatch, pair_vector
from benchmark.factorized_ranker.ranking import select_candidate, select_count
from benchmark.factorized_ranker.training import (
    MODEL_GRID,
    CorpusTables,
    HeadFeatureSpec,
    Hyperparameters,
    head_feature_spec,
    validate_retained_families,
)


ARTIFACT_SCHEMA_VERSION = 1
MODEL_INPUT_DTYPE = "float32"
THRESHOLD_POLICY = "floor_to_f32"
MAX_TREES_PER_HEAD = 96
MAX_COMBINED_TREES = 192
MAX_COMBINED_NODES = 2_880
U16_SENTINEL = 2**16 - 1
GOLDEN_SCHEMA_VERSION = 1


@dataclass(frozen=True)
class TreeNode:
    feature_index: int
    threshold: float
    left_child: int
    right_child: int
    leaf_value: float


@dataclass(frozen=True)
class FrozenTree:
    nodes: tuple[TreeNode, ...]


@dataclass(frozen=True)
class ModelProvenance:
    corpus_manifest_sha256: str
    fold_manifest_sha256: str
    cv_report_sha256: str
    ablation_report_sha256: str
    oof_predictions_sha256: str
    feature_schema_sha256: str
    feature_dump_binary_sha256: str
    source_git_commit: str
    training_command: tuple[str, ...]
    python_version: str
    numpy_version: str
    pandas_version: str
    scipy_version: str
    sklearn_version: str


_ARTIFACT_KEYS = {
    "schema_version",
    "head",
    "input_dtype",
    "threshold_policy",
    "feature_names",
    "feature_count",
    "feature_names_sha256",
    "retained_feature_families",
    "seed",
    "training_command",
    "source_git_commit",
    "feature_dump_binary_sha256",
    "corpus_manifest_sha256",
    "fold_manifest_sha256",
    "cv_report_sha256",
    "ablation_report_sha256",
    "oof_predictions_sha256",
    "feature_schema_sha256",
    "python_version",
    "numpy_version",
    "pandas_version",
    "scipy_version",
    "sklearn_version",
    "n_estimators",
    "max_depth",
    "learning_rate",
    "min_samples_leaf",
    "loss",
    "random_state",
    "initial_log_odds",
    "tree_count",
    "node_count",
    "parity_audit",
    "reference_pair",
    "trees",
}
_TREE_KEYS = {"nodes"}
_NODE_KEYS = {
    "feature_index",
    "threshold",
    "left_child",
    "right_child",
    "leaf_value",
}
_PARITY_KEYS = {
    "row_count",
    "column_count",
    "matrix_sha256",
    "batch_size",
    "internal_threshold_count",
    "max_raw_score_delta",
    "max_probability_delta",
}
_REFERENCE_KEYS = {
    "chain_id",
    "left_id",
    "right_id",
    "shared_names",
    "item_names",
    "shared_values",
    "left_values",
    "right_values",
}
_SHARED_ARTIFACT_FIELDS = (
    "schema_version",
    "input_dtype",
    "threshold_policy",
    "retained_feature_families",
    "seed",
    "training_command",
    "source_git_commit",
    "feature_dump_binary_sha256",
    "corpus_manifest_sha256",
    "fold_manifest_sha256",
    "cv_report_sha256",
    "ablation_report_sha256",
    "oof_predictions_sha256",
    "feature_schema_sha256",
    "python_version",
    "numpy_version",
    "pandas_version",
    "scipy_version",
    "sklearn_version",
)


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _mapping(value: object, keys: set[str], description: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != keys:
        raise ValueError(f"{description} keys do not match the frozen schema")
    if any(not isinstance(key, str) for key in value):
        raise ValueError(f"{description} keys must be strings")
    return value


def _integer(
    value: object,
    description: str,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> int:
    if not isinstance(value, Integral) or isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{description} must be an integer")
    number = int(value)
    if minimum is not None and number < minimum:
        raise ValueError(f"{description} is below its minimum")
    if maximum is not None and number > maximum:
        raise ValueError(f"{description} exceeds its maximum")
    return number


def _finite(value: object, description: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{description} must be a real number")
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


def _git_hash(value: object) -> str:
    text = _text(value, "source Git commit")
    if len(text) != 40 or any(character not in "0123456789abcdef" for character in text):
        raise ValueError("source Git commit must be a lowercase 40-hex hash")
    return text


def _string_tuple(value: object, description: str) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{description} must be a string sequence")
    result = tuple(_text(item, description) for item in value)
    if len(set(result)) != len(result):
        raise ValueError(f"{description} contains duplicates")
    return result


def _validate_command(value: object, description: str) -> tuple[str, ...]:
    command = _string_tuple(value, description)
    if not command:
        raise ValueError(f"{description} is empty")
    if any(token.startswith("/") or "=/" in token for token in command):
        raise ValueError(f"{description} contains a physical path")
    return command


def _validate_versions(provenance: ModelProvenance) -> None:
    expected = {
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": scipy.__version__,
        "sklearn_version": sklearn.__version__,
    }
    for field, current in expected.items():
        recorded = _text(getattr(provenance, field), field)
        if recorded != current:
            raise ValueError(f"recorded {field} does not match the export environment")


def _validate_provenance(provenance: ModelProvenance) -> None:
    if type(provenance) is not ModelProvenance:
        raise ValueError("provenance must be ModelProvenance")
    for field in (
        "corpus_manifest_sha256",
        "fold_manifest_sha256",
        "cv_report_sha256",
        "ablation_report_sha256",
        "oof_predictions_sha256",
        "feature_schema_sha256",
        "feature_dump_binary_sha256",
    ):
        _hash(getattr(provenance, field), field)
    _git_hash(provenance.source_git_commit)
    _validate_command(provenance.training_command, "training command")
    _validate_versions(provenance)


def _feature_names_hash(names: Sequence[str]) -> str:
    encoded = json.dumps(
        list(names),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return _sha256(encoded)


def _matrix_hash(x: np.ndarray, feature_names: Sequence[str]) -> str:
    canonical = np.ascontiguousarray(x, dtype="<f8")
    digest = hashlib.sha256()
    digest.update(b"factorized-audit-matrix-v1\0")
    digest.update(canonical.shape[0].to_bytes(8, "big", signed=False))
    digest.update(canonical.shape[1].to_bytes(8, "big", signed=False))
    names = json.dumps(
        list(feature_names),
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    digest.update(len(names).to_bytes(8, "big", signed=False))
    digest.update(names)
    digest.update(b"<f8\0C\0")
    digest.update(canonical.tobytes(order="C"))
    return digest.hexdigest()


def _canonical_threshold(value: float) -> float:
    threshold = _finite(value, "sklearn threshold")
    with np.errstate(over="ignore", invalid="ignore"):
        threshold32 = np.float32(threshold)
    if not np.isfinite(threshold32):
        raise ValueError("sklearn threshold overflows float32")
    if float(threshold32) > threshold:
        threshold32 = np.nextafter(
            threshold32, np.float32(-np.inf), dtype=np.float32
        )
    result = float(threshold32)
    if not math.isfinite(result) or float(np.float32(result)) != result:
        raise ValueError("canonical threshold is not a finite float32 value")
    return result


def _canonical_inputs(x: np.ndarray | Sequence[Sequence[float]], feature_count: int) -> np.ndarray:
    try:
        values = np.asarray(x, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("model input cannot be converted to float64") from error
    if values.ndim != 2 or values.shape[1] != feature_count:
        raise ValueError("model input has the wrong two-dimensional shape")
    if not np.isfinite(values).all():
        raise ValueError("model input contains a non-finite f64 value")
    with np.errstate(over="ignore", invalid="ignore"):
        values32 = values.astype(np.float32)
    if not np.isfinite(values32).all():
        raise ValueError("model input overflows float32")
    return values32.astype(np.float64)


def _node_mapping(node: TreeNode) -> dict[str, object]:
    return {
        "feature_index": node.feature_index,
        "threshold": node.threshold,
        "left_child": node.left_child,
        "right_child": node.right_child,
        "leaf_value": node.leaf_value,
    }


def _tree_mapping(tree: FrozenTree) -> dict[str, object]:
    return {"nodes": [_node_mapping(node) for node in tree.nodes]}


def _validate_tree(
    value: object,
    feature_count: int,
    tree_index: int,
) -> tuple[int, int]:
    tree = _mapping(value, _TREE_KEYS, f"tree {tree_index}")
    nodes_value = tree["nodes"]
    if not isinstance(nodes_value, list) or not nodes_value:
        raise ValueError(f"tree {tree_index} nodes must be a nonempty list")
    nodes = [
        _mapping(node, _NODE_KEYS, f"tree {tree_index} node {index}")
        for index, node in enumerate(nodes_value)
    ]
    state = [0] * len(nodes)
    maximum_depth = 0

    def visit(index: int, depth: int) -> None:
        nonlocal maximum_depth
        if not 0 <= index < len(nodes):
            raise ValueError(f"tree {tree_index} child index is out of range")
        if state[index] == 1:
            raise ValueError(f"tree {tree_index} contains a cycle")
        if state[index] == 2:
            raise ValueError(f"tree {tree_index} contains a shared child")
        state[index] = 1
        maximum_depth = max(maximum_depth, depth)
        node = nodes[index]
        feature = _integer(node["feature_index"], "node feature index")
        left = _integer(node["left_child"], "node left child")
        right = _integer(node["right_child"], "node right child")
        threshold = _finite(node["threshold"], "node threshold")
        leaf = _finite(node["leaf_value"], "node leaf value")
        if feature == -1:
            if left != -1 or right != -1 or threshold != 0.0:
                raise ValueError("leaf sentinels or threshold are invalid")
        else:
            if not 0 <= feature < feature_count:
                raise ValueError("internal feature index is out of range")
            if left < 0 or right < 0 or left == right or left == index or right == index:
                raise ValueError("internal child indices are invalid")
            if leaf != 0.0:
                raise ValueError("internal node leaf value is not zero")
            if float(np.float32(threshold)) != threshold:
                raise ValueError("internal threshold is not canonical float32")
            visit(left, depth + 1)
            visit(right, depth + 1)
        state[index] = 2

    visit(0, 0)
    if any(value != 2 for value in state):
        raise ValueError(f"tree {tree_index} contains an unreachable node")
    if maximum_depth > 3:
        raise ValueError(f"tree {tree_index} exceeds maximum depth 3")
    internal_count = sum(
        1 for node in nodes if int(node["feature_index"]) != -1
    )
    return maximum_depth, internal_count


def _validate_reference(
    value: object,
    spec: HeadFeatureSpec,
) -> None:
    reference = _mapping(value, _REFERENCE_KEYS, "reference pair")
    chain_id = _text(reference["chain_id"], "reference chain ID")
    left_id = _text(reference["left_id"], "reference left ID")
    right_id = _text(reference["right_id"], "reference right ID")
    if left_id == right_id:
        raise ValueError("reference item IDs are duplicated")
    if not chain_id:
        raise ValueError("reference chain ID is empty")
    shared_names = _string_tuple(reference["shared_names"], "reference shared names")
    item_names = _string_tuple(reference["item_names"], "reference item names")
    if shared_names != spec.shared_features or item_names != spec.item_features:
        raise ValueError("reference raw feature names disagree with the head spec")

    def vector(name: str, length: int) -> list[float]:
        raw = reference[name]
        if not isinstance(raw, list) or len(raw) != length:
            raise ValueError(f"reference {name} has the wrong length")
        return [_finite(item, f"reference {name}") for item in raw]

    shared = vector("shared_values", len(shared_names))
    left = vector("left_values", len(item_names))
    right = vector("right_values", len(item_names))
    assembled = pair_vector(shared, left, right)
    if len(assembled) != len(spec.pair_features):
        raise ValueError("reference pair assembly has the wrong length")


def validate_artifact(
    artifact: Mapping[str, object],
    *,
    expected_head: Literal["count", "candidate"] | None = None,
) -> None:
    value = _mapping(artifact, _ARTIFACT_KEYS, "model artifact")
    if _integer(value["schema_version"], "artifact schema version") != ARTIFACT_SCHEMA_VERSION:
        raise ValueError("artifact schema version is unsupported")
    head = _text(value["head"], "artifact head")
    if head not in {"count", "candidate"}:
        raise ValueError("artifact head is invalid")
    if expected_head is not None and head != expected_head:
        raise ValueError("artifact head does not match the expected role")
    if value["input_dtype"] != MODEL_INPUT_DTYPE:
        raise ValueError("artifact input dtype is unsupported")
    if value["threshold_policy"] != THRESHOLD_POLICY:
        raise ValueError("artifact threshold policy is unsupported")
    families = validate_retained_families(
        _string_tuple(value["retained_feature_families"], "retained families")
    )
    spec = head_feature_spec(head, families)  # type: ignore[arg-type]
    names = _string_tuple(value["feature_names"], "artifact feature names")
    if names != spec.pair_features:
        raise ValueError("artifact feature names disagree with the exact head pair schema")
    feature_count = _integer(value["feature_count"], "artifact feature count", minimum=1)
    if feature_count != len(names):
        raise ValueError("artifact feature count disagrees with names")
    if _hash(value["feature_names_sha256"], "feature names hash") != _feature_names_hash(names):
        raise ValueError("artifact feature names hash mismatch")
    seed = _integer(value["seed"], "artifact seed")
    if seed != 37:
        raise ValueError("artifact seed must be 37")
    _validate_command(value["training_command"], "artifact training command")
    _git_hash(value["source_git_commit"])
    for field in (
        "feature_dump_binary_sha256",
        "corpus_manifest_sha256",
        "fold_manifest_sha256",
        "cv_report_sha256",
        "ablation_report_sha256",
        "oof_predictions_sha256",
        "feature_schema_sha256",
    ):
        _hash(value[field], field)
    for field in (
        "python_version",
        "numpy_version",
        "pandas_version",
        "scipy_version",
        "sklearn_version",
    ):
        _text(value[field], field)
    params = Hyperparameters(
        _integer(value["n_estimators"], "n_estimators"),
        _finite(value["learning_rate"], "learning_rate"),
        _integer(value["min_samples_leaf"], "min_samples_leaf"),
        _integer(value["max_depth"], "max_depth"),
    )
    if params not in MODEL_GRID:
        raise ValueError("artifact hyperparameters are outside the approved grid")
    if value["loss"] != "log_loss":
        raise ValueError("artifact loss is unsupported")
    if _integer(value["random_state"], "random state") != seed:
        raise ValueError("artifact random state disagrees with seed")
    _finite(value["initial_log_odds"], "initial log odds")
    tree_count = _integer(value["tree_count"], "tree count", minimum=1, maximum=MAX_TREES_PER_HEAD)
    if tree_count != params.n_estimators:
        raise ValueError("tree count disagrees with n_estimators")
    trees = value["trees"]
    if not isinstance(trees, list) or len(trees) != tree_count:
        raise ValueError("artifact trees have the wrong count")
    internal_count = 0
    node_count = 0
    for tree_index, tree in enumerate(trees):
        _depth, internals = _validate_tree(tree, feature_count, tree_index)
        internal_count += internals
        assert isinstance(tree, Mapping)
        node_count += len(tree["nodes"])  # type: ignore[arg-type]
    if node_count != _integer(value["node_count"], "node count", minimum=1):
        raise ValueError("artifact node count mismatch")
    if node_count >= U16_SENTINEL:
        raise ValueError("per-head flattened node count exceeds u16")
    parity = _mapping(value["parity_audit"], _PARITY_KEYS, "parity audit")
    rows = _integer(parity["row_count"], "parity row count", minimum=1)
    columns = _integer(parity["column_count"], "parity column count", minimum=1)
    if columns != feature_count:
        raise ValueError("parity audit feature count mismatch")
    _hash(parity["matrix_sha256"], "parity matrix hash")
    batch_size = _integer(parity["batch_size"], "parity batch size", minimum=1, maximum=4096)
    if batch_size > rows:
        raise ValueError("parity batch size exceeds audit rows")
    if _integer(parity["internal_threshold_count"], "internal threshold count", minimum=0) != internal_count:
        raise ValueError("parity internal threshold count mismatch")
    if _finite(parity["max_raw_score_delta"], "maximum raw-score delta") > 1e-12:
        raise ValueError("raw-score parity exceeds tolerance")
    if _finite(parity["max_probability_delta"], "maximum probability delta") > 1e-12:
        raise ValueError("probability parity exceeds tolerance")
    _validate_reference(value["reference_pair"], spec)


def validate_artifacts(
    count_artifact: Mapping[str, object],
    candidate_artifact: Mapping[str, object],
) -> None:
    validate_artifact(count_artifact, expected_head="count")
    validate_artifact(candidate_artifact, expected_head="candidate")
    for field in _SHARED_ARTIFACT_FIELDS:
        if count_artifact[field] != candidate_artifact[field]:
            raise ValueError(f"cross-head artifact field {field} disagrees")
    combined_trees = int(count_artifact["tree_count"]) + int(candidate_artifact["tree_count"])
    combined_nodes = int(count_artifact["node_count"]) + int(candidate_artifact["node_count"])
    if combined_trees > MAX_COMBINED_TREES:
        raise ValueError("combined tree count exceeds the frozen cap")
    if combined_nodes > MAX_COMBINED_NODES:
        raise ValueError("combined node count exceeds the frozen cap")


def raw_tree_value(tree: Mapping[str, object], canonical_x: np.ndarray) -> float:
    if not isinstance(canonical_x, np.ndarray) or canonical_x.dtype != np.float64 or canonical_x.ndim != 1:
        raise ValueError("tree input must be a canonical float64 vector")
    if not np.isfinite(canonical_x).all():
        raise ValueError("tree input contains a non-finite value")
    mapping = _mapping(tree, _TREE_KEYS, "tree")
    nodes = mapping["nodes"]
    if not isinstance(nodes, list) or not nodes:
        raise ValueError("tree nodes are invalid")
    index = 0
    visited: set[int] = set()
    while True:
        if index in visited or not 0 <= index < len(nodes):
            raise ValueError("tree traversal encountered invalid topology")
        visited.add(index)
        node = _mapping(nodes[index], _NODE_KEYS, "tree node")
        feature = _integer(node["feature_index"], "tree feature index")
        if feature == -1:
            return _finite(node["leaf_value"], "tree leaf value")
        if not 0 <= feature < len(canonical_x):
            raise ValueError("tree feature index is out of range")
        threshold = _finite(node["threshold"], "tree threshold")
        index = _integer(
            node["left_child"] if canonical_x[feature] <= threshold else node["right_child"],
            "tree child index",
        )


def _raw_score_validated(artifact: Mapping[str, object], canonical: np.ndarray) -> np.ndarray:
    trees = artifact["trees"]
    assert isinstance(trees, list)
    output = np.empty(canonical.shape[0], dtype=np.float64)
    initial = float(artifact["initial_log_odds"])
    rate = float(artifact["learning_rate"])
    for row_index, row in enumerate(canonical):
        ordered_sum = 0.0
        for tree in trees:
            assert isinstance(tree, Mapping)
            ordered_sum += raw_tree_value(tree, row)
        output[row_index] = initial + rate * ordered_sum
    if not np.isfinite(output).all():
        raise ValueError("artifact raw score is non-finite")
    return output


def raw_artifact_score(
    artifact: Mapping[str, object],
    x: np.ndarray,
) -> np.ndarray:
    validate_artifact(artifact)
    canonical = _canonical_inputs(x, int(artifact["feature_count"]))
    return _raw_score_validated(artifact, canonical)


def _sigmoid(raw: np.ndarray) -> np.ndarray:
    output = np.empty_like(raw, dtype=np.float64)
    positive = raw >= 0.0
    output[positive] = 1.0 / (1.0 + np.exp(-raw[positive]))
    exponential = np.exp(raw[~positive])
    output[~positive] = exponential / (1.0 + exponential)
    if not np.isfinite(output).all() or np.any(output < 0.0) or np.any(output > 1.0):
        raise ValueError("artifact probability is invalid")
    return output


def predict_artifact(
    artifact: Mapping[str, object],
    x: np.ndarray,
) -> np.ndarray:
    return _sigmoid(raw_artifact_score(artifact, x))


def _synthetic_reference(head: str, spec: HeadFeatureSpec) -> dict[str, object]:
    shared = [0.0] * len(spec.shared_features)
    left = [0.0] * len(spec.item_features)
    right = [0.0] * len(spec.item_features)
    if head == "count" and left:
        left[0] = 1.0
        right[0] = 2.0
    left_id, right_id = sorted(
        (
            _sha256(f"{head}:synthetic:left".encode("utf-8")),
            _sha256(f"{head}:synthetic:right".encode("utf-8")),
        )
    )
    return {
        "chain_id": f"synthetic-{head}-reference",
        "left_id": left_id,
        "right_id": right_id,
        "shared_names": list(spec.shared_features),
        "item_names": list(spec.item_features),
        "shared_values": shared,
        "left_values": left,
        "right_values": right,
    }


def reference_pair_from_batch(
    corpus: CorpusTables,
    batch: PairBatch,
    head_name: Literal["count", "candidate"],
    retained_feature_families: Sequence[str],
) -> dict[str, object]:
    """Reconstruct the lexically first canonical raw pair in a final Task 10 batch."""
    if not isinstance(corpus, CorpusTables) or not isinstance(batch, PairBatch):
        raise ValueError("reference construction requires CorpusTables and PairBatch")
    if head_name not in {"count", "candidate"}:
        raise ValueError("reference head must be count or candidate")
    families = validate_retained_families(retained_feature_families)
    spec = head_feature_spec(head_name, families)
    if batch.feature_names != spec.pair_features:
        raise ValueError("reference batch feature names disagree with the head spec")
    arrays = (batch.x, batch.chain_ids, batch.left_ids, batch.right_ids)
    if (
        batch.x.dtype != np.float64
        or batch.x.ndim != 2
        or batch.x.shape[1] != len(spec.pair_features)
        or any(array.ndim != 1 or len(array) != len(batch.x) for array in arrays[1:])
    ):
        raise ValueError("reference batch arrays are malformed")
    canonical_rows = sorted(
        (
            str(batch.chain_ids[index]),
            str(batch.left_ids[index]),
            str(batch.right_ids[index]),
            index,
        )
        for index in range(len(batch.x))
        if str(batch.left_ids[index]) < str(batch.right_ids[index])
    )
    if not canonical_rows:
        raise ValueError("final pair batch has no canonical orientation")
    chain_id, left_id, right_id, row_index = canonical_rows[0]
    chain_rows = pair_api._validate_chains(corpus.chains)
    if chain_id not in chain_rows:
        raise ValueError("reference batch chain does not join to corpus chains")
    chain = chain_rows[chain_id]
    if head_name == "count":
        grouped = pair_api._validate_counts(corpus.counts, chain_rows)
    else:
        grouped = pair_api._validate_candidates(corpus.candidates, chain_rows)
    items = {item.item_id: item for item in grouped[chain_id]}
    if set((left_id, right_id)) - set(items):
        raise ValueError("reference batch item IDs do not join to corpus rows")
    shared = [chain.values[name] for name in spec.shared_features]
    left = [items[left_id].values[name] for name in spec.item_features]
    right = [items[right_id].values[name] for name in spec.item_features]
    assembled = pair_vector(shared, left, right)
    if not np.array_equal(assembled, batch.x[row_index]):
        raise ValueError("reconstructed raw reference does not reproduce final batch row")
    matches = np.flatnonzero(
        (batch.chain_ids == chain_id)
        & (batch.left_ids == right_id)
        & (batch.right_ids == left_id)
    )
    if len(matches) != 1 or not np.array_equal(
        pair_vector(shared, right, left), batch.x[int(matches[0])]
    ):
        raise ValueError("reference reverse orientation is absent or inconsistent")
    reference = {
        "chain_id": chain_id,
        "left_id": left_id,
        "right_id": right_id,
        "shared_names": list(spec.shared_features),
        "item_names": list(spec.item_features),
        "shared_values": shared,
        "left_values": left,
        "right_values": right,
    }
    _validate_reference(reference, spec)
    return reference


def _extract_trees(
    model: GradientBoostingClassifier,
    feature_count: int,
) -> tuple[list[dict[str, object]], list[tuple[float, float]], int]:
    trees: list[dict[str, object]] = []
    threshold_pairs: list[tuple[float, float]] = []
    node_count = 0
    for stage in range(model.estimators_.shape[0]):
        native = model.estimators_[stage, 0].tree_
        if int(native.node_count) <= 0:
            raise ValueError("sklearn tree has no nodes")
        nodes: list[TreeNode] = []
        for index in range(int(native.node_count)):
            left = int(native.children_left[index])
            right = int(native.children_right[index])
            if left == -1 and right == -1:
                nodes.append(
                    TreeNode(
                        -1,
                        0.0,
                        -1,
                        -1,
                        _finite(native.value[index, 0, 0], "sklearn leaf value"),
                    )
                )
            else:
                feature = int(native.feature[index])
                if not 0 <= feature < feature_count:
                    raise ValueError("sklearn internal feature index is invalid")
                original = _finite(native.threshold[index], "sklearn threshold")
                canonical = _canonical_threshold(original)
                previous = np.nextafter(
                    np.float32(canonical), np.float32(-np.inf), dtype=np.float32
                )
                equal = np.float32(canonical)
                following = np.nextafter(
                    np.float32(canonical), np.float32(np.inf), dtype=np.float32
                )
                if not all(np.isfinite(value) for value in (previous, equal, following)):
                    raise ValueError("threshold predicate probes are non-finite")
                for probe in (previous, equal, following):
                    if (float(probe) <= original) != (float(probe) <= canonical):
                        raise ValueError("canonical threshold changes a float32 branch predicate")
                threshold_pairs.append((original, canonical))
                nodes.append(TreeNode(feature, canonical, left, right, 0.0))
        frozen = FrozenTree(tuple(nodes))
        mapping = _tree_mapping(frozen)
        depth, _internals = _validate_tree(mapping, feature_count, stage)
        if depth != int(native.max_depth):
            raise ValueError("frozen tree depth disagrees with sklearn")
        trees.append(mapping)
        node_count += len(nodes)
    return trees, threshold_pairs, node_count


def _validate_estimator(
    model: GradientBoostingClassifier,
    feature_count: int,
    params: Hyperparameters,
    seed: int,
) -> None:
    if type(model) is not GradientBoostingClassifier:
        raise ValueError("model must be exactly GradientBoostingClassifier")
    if type(params) is not Hyperparameters or params not in MODEL_GRID:
        raise ValueError("selected hyperparameters are outside the approved grid")
    if seed != 37 or isinstance(seed, bool):
        raise ValueError("artifact freeze seed must be 37")
    options = model.get_params(deep=False)
    expected = {
        "n_estimators": params.n_estimators,
        "learning_rate": params.learning_rate,
        "min_samples_leaf": params.min_samples_leaf,
        "max_depth": params.max_depth,
        "loss": "log_loss",
        "random_state": seed,
        "criterion": "friedman_mse",
        "min_samples_split": 2,
        "min_weight_fraction_leaf": 0.0,
        "subsample": 1.0,
        "max_features": None,
        "max_leaf_nodes": None,
        "min_impurity_decrease": 0.0,
        "init": None,
        "validation_fraction": 0.1,
        "n_iter_no_change": None,
        "tol": 0.0001,
        "ccp_alpha": 0.0,
        "warm_start": False,
        "verbose": 0,
    }
    for name, expected_value in expected.items():
        if name not in options or options[name] != expected_value:
            raise ValueError(f"sklearn option {name} is outside the frozen contract")
    classes = getattr(model, "classes_", None)
    if (
        not isinstance(classes, np.ndarray)
        or classes.ndim != 1
        or len(classes) != 2
        or classes.dtype.kind not in "iuf"
        or not np.array_equal(classes, np.asarray([0, 1]))
    ):
        raise ValueError("sklearn classes must be exact ordered numeric [0,1]")
    if int(getattr(model, "n_features_in_", -1)) != feature_count:
        raise ValueError("sklearn feature count disagrees with artifact names")
    estimators = getattr(model, "estimators_", None)
    if not isinstance(estimators, np.ndarray) or estimators.shape != (params.n_estimators, 1):
        raise ValueError("sklearn estimator array has the wrong shape")
    if not callable(getattr(model, "_raw_predict_init", None)):
        raise ValueError("sklearn private raw initializer is unavailable")


def _raw_initializer(model: GradientBoostingClassifier, feature_count: int) -> float:
    observed: list[float] = []
    for fill in (0.0, 1.0, -1.0):
        probe = np.full((1, feature_count), fill, dtype=np.float32)
        try:
            raw = model._raw_predict_init(probe)
        except Exception as error:  # sklearn private API: fail closed with context
            raise ValueError("sklearn private raw initializer failed") from error
        if not isinstance(raw, np.ndarray) or raw.shape != (1, 1):
            raise ValueError("sklearn private raw initializer returned the wrong shape")
        observed.append(_finite(raw[0, 0], "sklearn initial log odds"))
    if any(struct.pack(">d", value) != struct.pack(">d", observed[0]) for value in observed[1:]):
        raise ValueError("sklearn private raw initializer is input-dependent")
    return observed[0]


def _freeze_classifier(
    model: GradientBoostingClassifier,
    head_name: Literal["count", "candidate"],
    feature_names: Sequence[str],
    retained_feature_families: Sequence[str],
    selected_hyperparameters: Hyperparameters,
    provenance: ModelProvenance,
    audit_x: np.ndarray,
    *,
    seed: int = 37,
    reference_pair: Mapping[str, object] | None = None,
) -> dict[str, object]:
    if head_name not in {"count", "candidate"}:
        raise ValueError("head name must be count or candidate")
    families = validate_retained_families(retained_feature_families)
    spec = head_feature_spec(head_name, families)
    names = _string_tuple(feature_names, "feature names")
    if names != spec.pair_features:
        raise ValueError("freeze feature names disagree with the exact head pair schema")
    _validate_provenance(provenance)
    _validate_estimator(model, len(names), selected_hyperparameters, seed)
    if not isinstance(audit_x, np.ndarray) or audit_x.dtype != np.float64 or audit_x.ndim != 2:
        raise ValueError("audit matrix must be a two-dimensional float64 ndarray")
    if audit_x.shape[0] == 0 or audit_x.shape[1] != len(names):
        raise ValueError("audit matrix has the wrong shape")
    if not audit_x.flags.c_contiguous or not np.isfinite(audit_x).all():
        raise ValueError("audit matrix must be finite C-contiguous float64")
    canonical_audit = _canonical_inputs(audit_x, len(names))
    initial = _raw_initializer(model, len(names))
    trees, threshold_pairs, node_count = _extract_trees(model, len(names))
    reference = (
        dict(reference_pair)
        if reference_pair is not None
        else _synthetic_reference(head_name, spec)
    )
    artifact: dict[str, object] = {
        "schema_version": ARTIFACT_SCHEMA_VERSION,
        "head": head_name,
        "input_dtype": MODEL_INPUT_DTYPE,
        "threshold_policy": THRESHOLD_POLICY,
        "feature_names": list(names),
        "feature_count": len(names),
        "feature_names_sha256": _feature_names_hash(names),
        "retained_feature_families": list(families),
        "seed": seed,
        "training_command": list(provenance.training_command),
        "source_git_commit": provenance.source_git_commit,
        "feature_dump_binary_sha256": provenance.feature_dump_binary_sha256,
        "corpus_manifest_sha256": provenance.corpus_manifest_sha256,
        "fold_manifest_sha256": provenance.fold_manifest_sha256,
        "cv_report_sha256": provenance.cv_report_sha256,
        "ablation_report_sha256": provenance.ablation_report_sha256,
        "oof_predictions_sha256": provenance.oof_predictions_sha256,
        "feature_schema_sha256": provenance.feature_schema_sha256,
        "python_version": provenance.python_version,
        "numpy_version": provenance.numpy_version,
        "pandas_version": provenance.pandas_version,
        "scipy_version": provenance.scipy_version,
        "sklearn_version": provenance.sklearn_version,
        "n_estimators": selected_hyperparameters.n_estimators,
        "max_depth": selected_hyperparameters.max_depth,
        "learning_rate": selected_hyperparameters.learning_rate,
        "min_samples_leaf": selected_hyperparameters.min_samples_leaf,
        "loss": "log_loss",
        "random_state": seed,
        "initial_log_odds": initial,
        "tree_count": len(trees),
        "node_count": node_count,
        "parity_audit": {
            "row_count": len(audit_x),
            "column_count": audit_x.shape[1],
            "matrix_sha256": _matrix_hash(audit_x, names),
            "batch_size": min(4096, len(audit_x)),
            "internal_threshold_count": len(threshold_pairs),
            "max_raw_score_delta": 0.0,
            "max_probability_delta": 0.0,
        },
        "reference_pair": reference,
        "trees": trees,
    }
    validate_artifact(artifact, expected_head=head_name)
    maximum_raw_delta = 0.0
    maximum_probability_delta = 0.0
    batch_size = min(4096, len(audit_x))
    for start in range(0, len(audit_x), batch_size):
        stop = min(start + batch_size, len(audit_x))
        canonical64 = canonical_audit[start:stop]
        canonical32 = canonical64.astype(np.float32)
        expected_raw = np.asarray(model.decision_function(canonical32), dtype=np.float64)
        expected_probability = np.asarray(
            model.predict_proba(canonical32)[:, 1], dtype=np.float64
        )
        actual_raw = _raw_score_validated(artifact, canonical64)
        actual_probability = _sigmoid(actual_raw)
        if expected_raw.shape != actual_raw.shape or expected_probability.shape != actual_probability.shape:
            raise ValueError("sklearn parity output shape mismatch")
        maximum_raw_delta = max(
            maximum_raw_delta,
            float(np.max(np.abs(expected_raw - actual_raw))),
        )
        maximum_probability_delta = max(
            maximum_probability_delta,
            float(np.max(np.abs(expected_probability - actual_probability))),
        )
    if not math.isfinite(maximum_raw_delta) or maximum_raw_delta > 1e-12:
        raise ValueError("full audit raw-score parity exceeds 1e-12")
    if not math.isfinite(maximum_probability_delta) or maximum_probability_delta > 1e-12:
        raise ValueError("full audit probability parity exceeds 1e-12")
    artifact["parity_audit"]["max_raw_score_delta"] = maximum_raw_delta  # type: ignore[index]
    artifact["parity_audit"]["max_probability_delta"] = maximum_probability_delta  # type: ignore[index]
    validate_artifact(artifact, expected_head=head_name)
    return artifact


def freeze_classifier(
    model: GradientBoostingClassifier,
    head_name: Literal["count", "candidate"],
    feature_names: Sequence[str],
    retained_feature_families: Sequence[str],
    selected_hyperparameters: Hyperparameters,
    provenance: ModelProvenance,
    audit_x: np.ndarray,
    *,
    seed: int = 37,
) -> dict[str, object]:
    return _freeze_classifier(
        model,
        head_name,
        feature_names,
        retained_feature_families,
        selected_hyperparameters,
        provenance,
        audit_x,
        seed=seed,
    )


def _duplicate_rejecting_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value}")


def load_artifact(
    path: Path,
    *,
    expected_sha256: str | None = None,
    expected_head: Literal["count", "candidate"] | None = None,
) -> dict[str, object]:
    data = Path(path).read_bytes()
    if expected_sha256 is not None and _sha256(data) != _hash(expected_sha256, "expected artifact hash"):
        raise ValueError("artifact file hash mismatch")
    try:
        value = json.loads(
            data,
            object_pairs_hook=_duplicate_rejecting_object,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("artifact JSON is invalid") from error
    if not isinstance(value, dict):
        raise ValueError("artifact root must be an object")
    if _canonical_json_bytes(value) != data:
        raise ValueError("artifact bytes are not canonical JSON")
    validate_artifact(value, expected_head=expected_head)
    return value


def _write_temporary(path: Path, data: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("wb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    if temporary.read_bytes() != data:
        temporary.unlink(missing_ok=True)
        raise ValueError("temporary artifact bytes changed after writing")
    return temporary


def write_artifact(
    path: Path,
    artifact: Mapping[str, object],
) -> str:
    validate_artifact(artifact)
    data = _canonical_json_bytes(artifact)
    target = Path(path)
    temporary = _write_temporary(target, data)
    try:
        os.replace(temporary, target)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    installed = target.read_bytes()
    if installed != data:
        raise OSError(f"installed artifact bytes differ at {target}")
    return _sha256(installed)


def write_artifacts(
    count_path: Path,
    count_artifact: Mapping[str, object],
    candidate_path: Path,
    candidate_artifact: Mapping[str, object],
) -> tuple[str, str]:
    validate_artifacts(count_artifact, candidate_artifact)
    count_target = Path(count_path)
    candidate_target = Path(candidate_path)
    if count_target.resolve(strict=False) == candidate_target.resolve(strict=False):
        raise ValueError("count and candidate artifact paths alias")
    count_data = _canonical_json_bytes(count_artifact)
    candidate_data = _canonical_json_bytes(candidate_artifact)
    count_temp: Path | None = None
    candidate_temp: Path | None = None
    try:
        count_temp = _write_temporary(count_target, count_data)
        candidate_temp = _write_temporary(candidate_target, candidate_data)
        os.replace(count_temp, count_target)
        count_temp = None
        os.replace(candidate_temp, candidate_target)
        candidate_temp = None
    except BaseException as error:
        if count_temp is not None:
            count_temp.unlink(missing_ok=True)
        if candidate_temp is not None:
            candidate_temp.unlink(missing_ok=True)
        raise OSError(
            f"paired artifact installation failed for {count_target} and {candidate_target}"
        ) from error
    if count_target.read_bytes() != count_data or candidate_target.read_bytes() != candidate_data:
        raise OSError("installed paired artifact bytes differ")
    return _sha256(count_data), _sha256(candidate_data)


def _rust_f64_literal(value: object) -> str:
    number = _finite(value, "Rust f64 literal")
    token = repr(number)
    if "nan" in token.casefold() or "inf" in token.casefold():
        raise ValueError("Rust f64 literal is non-finite")
    parsed = float(token)
    if struct.pack(">d", parsed) != struct.pack(">d", number):
        raise ValueError("Rust f64 literal does not round-trip binary64")
    return f"{token}_f64"


def _rust_string_literal(value: object) -> str:
    if not isinstance(value, str):
        raise ValueError("Rust string value must be text")
    output = ['"']
    for character in value:
        if character == '"':
            output.append('\\"')
        elif character == "\\":
            output.append("\\\\")
        elif character == "\n":
            output.append("\\n")
        elif character == "\r":
            output.append("\\r")
        elif character == "\t":
            output.append("\\t")
        elif character == "\0":
            output.append("\\0")
        elif ord(character) == 127 or unicodedata.category(character) == "Cc":
            output.append(f"\\u{{{ord(character):x}}}")
        else:
            output.append(character)
    output.append('"')
    return "".join(output)


def _flatten_head(artifact: Mapping[str, object]) -> tuple[list[dict[str, object]], list[int]]:
    flat_nodes: list[dict[str, object]] = []
    roots: list[int] = []
    trees = artifact["trees"]
    assert isinstance(trees, list)
    for tree in trees:
        assert isinstance(tree, Mapping)
        nodes = tree["nodes"]
        assert isinstance(nodes, list)
        offset = len(flat_nodes)
        if offset >= U16_SENTINEL:
            raise ValueError("flattened tree root exceeds u16")
        roots.append(offset)
        for node in nodes:
            assert isinstance(node, Mapping)
            feature = int(node["feature_index"])
            if feature == -1:
                rendered = {
                    "feature": U16_SENTINEL,
                    "threshold": float(node["threshold"]),
                    "left": U16_SENTINEL,
                    "right": U16_SENTINEL,
                    "leaf_value": float(node["leaf_value"]),
                }
            else:
                left = offset + int(node["left_child"])
                right = offset + int(node["right_child"])
                if max(feature, left, right) >= U16_SENTINEL:
                    raise ValueError("flattened internal index exceeds u16")
                rendered = {
                    "feature": feature,
                    "threshold": float(node["threshold"]),
                    "left": left,
                    "right": right,
                    "leaf_value": float(node["leaf_value"]),
                }
            flat_nodes.append(rendered)
    return flat_nodes, roots


def _render_head(prefix: str, artifact: Mapping[str, object]) -> list[str]:
    nodes, roots = _flatten_head(artifact)
    names = artifact["feature_names"]
    assert isinstance(names, list)
    lines = [
        f"pub(crate) static {prefix}_FEATURE_NAMES: [&str; {len(names)}] = [",
        *(f"    {_rust_string_literal(name)}," for name in names),
        "];",
        "",
        f"pub(crate) static {prefix}_NODES: [StaticNode; {len(nodes)}] = [",
    ]
    for node in nodes:
        feature = "u16::MAX" if node["feature"] == U16_SENTINEL else f"{node['feature']}_u16"
        left = "u16::MAX" if node["left"] == U16_SENTINEL else f"{node['left']}_u16"
        right = "u16::MAX" if node["right"] == U16_SENTINEL else f"{node['right']}_u16"
        lines.extend(
            (
                "    StaticNode {",
                f"        feature: {feature},",
                f"        threshold: {_rust_f64_literal(node['threshold'])},",
                f"        left: {left},",
                f"        right: {right},",
                f"        leaf_value: {_rust_f64_literal(node['leaf_value'])},",
                "    },",
            )
        )
    lines.extend(("];", "", f"pub(crate) static {prefix}_TREES: [StaticTree; {len(roots)}] = ["))
    lines.extend(f"    StaticTree {{ root: {root}_u16 }}," for root in roots)
    lines.extend(
        (
            "];",
            "",
            f"pub(crate) static {prefix}_MODEL: StaticBoostedModel = StaticBoostedModel {{",
            f"    schema_version: {ARTIFACT_SCHEMA_VERSION}_u32,",
            f"    feature_names: &{prefix}_FEATURE_NAMES,",
            f"    initial_log_odds: {_rust_f64_literal(artifact['initial_log_odds'])},",
            f"    learning_rate: {_rust_f64_literal(artifact['learning_rate'])},",
            f"    trees: &{prefix}_TREES,",
            f"    nodes: &{prefix}_NODES,",
            "};",
            "",
        )
    )
    return lines


def _rustfmt_source(source: str) -> str:
    result = subprocess.run(
        ["rustfmt", "--edition", "2021"],
        input=source,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise ValueError(f"rustfmt rejected generated source: {result.stderr.strip()}")
    if not result.stdout.endswith("\n"):
        raise ValueError("rustfmt returned source without a trailing newline")
    return result.stdout


def render_rust_models(
    count_artifact: Mapping[str, object],
    candidate_artifact: Mapping[str, object],
) -> str:
    validate_artifacts(count_artifact, candidate_artifact)
    families = count_artifact["retained_feature_families"]
    assert isinstance(families, list)
    lines = [
        "// @generated by benchmark.export_factorized_ranker; do not edit.",
        "use super::model::{StaticBoostedModel, StaticNode, StaticTree};",
        "",
        f"pub(crate) const MODEL_INPUT_DTYPE: &str = {_rust_string_literal(MODEL_INPUT_DTYPE)};",
        f"pub(crate) const MODEL_THRESHOLD_POLICY: &str = {_rust_string_literal(THRESHOLD_POLICY)};",
        "pub(crate) static RETAINED_FEATURE_FAMILIES: &[&str] = &[",
        *(f"    {_rust_string_literal(family)}," for family in families),
        "];",
        "",
        *_render_head("COUNT", count_artifact),
        *_render_head("CANDIDATE", candidate_artifact),
    ]
    return _rustfmt_source("\n".join(lines))


def _canonical_vector(values: np.ndarray) -> list[float]:
    return _canonical_inputs(values.reshape(1, -1), len(values))[0].tolist()


def _representative_split(artifact: Mapping[str, object]) -> dict[str, object]:
    names = artifact["feature_names"]
    reference = artifact["reference_pair"]
    assert isinstance(names, list) and isinstance(reference, Mapping)
    shared_count = len(reference["shared_names"])  # type: ignore[arg-type]
    item_count = len(reference["item_names"])  # type: ignore[arg-type]
    trees = artifact["trees"]
    assert isinstance(trees, list)
    for tree_index, tree in enumerate(trees):
        assert isinstance(tree, Mapping)
        nodes = tree["nodes"]
        assert isinstance(nodes, list)
        root = nodes[0]
        assert isinstance(root, Mapping)
        feature = int(root["feature_index"])
        if feature < 0:
            continue
        threshold = float(root["threshold"])
        previous = float(
            np.nextafter(np.float32(threshold), np.float32(-np.inf), dtype=np.float32)
        )
        following = float(
            np.nextafter(np.float32(threshold), np.float32(np.inf), dtype=np.float32)
        )
        if not all(math.isfinite(value) for value in (previous, threshold, following)):
            continue
        if feature < shared_count:
            kind = "shared"
            raw_index = feature
        elif feature < shared_count + item_count:
            kind = "signed_difference"
            raw_index = feature - shared_count
        else:
            kind = "absolute_difference"
            raw_index = feature - shared_count - item_count
            if raw_index >= item_count or previous < 0.0:
                continue
        return {
            "tree_index": tree_index,
            "node_index": 0,
            "feature_index": feature,
            "feature_name": names[feature],
            "raw_kind": kind,
            "raw_index": raw_index,
            "threshold": threshold,
            "previous": previous,
            "following": following,
        }
    raise ValueError("head has no realizable non-stump root for golden probes")


def _golden_case(
    artifact: Mapping[str, object],
    representative: Mapping[str, object],
    kind: str,
    target: float | None,
) -> dict[str, object]:
    reference = artifact["reference_pair"]
    assert isinstance(reference, Mapping)
    shared = [float(value) for value in reference["shared_values"]]  # type: ignore[index]
    left = [float(value) for value in reference["left_values"]]  # type: ignore[index]
    right = [float(value) for value in reference["right_values"]]  # type: ignore[index]
    if target is not None:
        index = int(representative["raw_index"])
        raw_kind = representative["raw_kind"]
        if raw_kind == "shared":
            shared[index] = target
        else:
            left[index] = target
            right[index] = 0.0
    forward = pair_vector(shared, left, right)
    reverse = pair_vector(shared, right, left)
    canonical_forward = np.asarray(_canonical_vector(forward), dtype=np.float64)
    canonical_reverse = np.asarray(_canonical_vector(reverse), dtype=np.float64)
    trees = artifact["trees"]
    assert isinstance(trees, list)
    forward_leaves = [raw_tree_value(tree, canonical_forward) for tree in trees]  # type: ignore[arg-type]
    reverse_leaves = [raw_tree_value(tree, canonical_reverse) for tree in trees]  # type: ignore[arg-type]
    forward_raw = float(raw_artifact_score(artifact, forward.reshape(1, -1))[0])
    reverse_raw = float(raw_artifact_score(artifact, reverse.reshape(1, -1))[0])
    forward_probability = float(predict_artifact(artifact, forward.reshape(1, -1))[0])
    reverse_probability = float(predict_artifact(artifact, reverse.reshape(1, -1))[0])
    symmetrized = 0.5 * (forward_probability + 1.0 - reverse_probability)
    head = str(artifact["head"])
    left_id = str(reference["left_id"])
    right_id = str(reference["right_id"])
    if head == "count":
        left_key, right_key = 1, 2
        scores = {str(left_key): symmetrized, str(right_key): 1.0 - symmetrized}
        legacy_count = 1
        winner: int | str = select_count(
            {left_key: symmetrized, right_key: 1.0 - symmetrized}, legacy_count
        )
        tie_input: dict[str, object] = {
            "kind": "count",
            "left_count": left_key,
            "right_count": right_key,
            "legacy_count": legacy_count,
        }
    else:
        scores = {left_id: symmetrized, right_id: 1.0 - symmetrized}
        winner = select_candidate(scores)
        tie_input = {"kind": "candidate", "tie_break": "lexical"}
    return {
        "kind": kind,
        "chain_id": reference["chain_id"],
        "left_id": left_id,
        "right_id": right_id,
        "shared_names": list(reference["shared_names"]),  # type: ignore[arg-type]
        "item_names": list(reference["item_names"]),  # type: ignore[arg-type]
        "pair_feature_names": list(artifact["feature_names"]),  # type: ignore[arg-type]
        "shared_values": shared,
        "left_values": left,
        "right_values": right,
        "forward_pair_vector": forward.tolist(),
        "reverse_pair_vector": reverse.tolist(),
        "forward_canonical_vector": canonical_forward.tolist(),
        "reverse_canonical_vector": canonical_reverse.tolist(),
        "forward_tree_leaf_values": forward_leaves,
        "reverse_tree_leaf_values": reverse_leaves,
        "forward_raw_margin": forward_raw,
        "reverse_raw_margin": reverse_raw,
        "forward_probability": forward_probability,
        "reverse_probability": reverse_probability,
        "symmetrized_probability": symmetrized,
        "borda_scores": scores,
        "tie_input": tie_input,
        "winner": winner,
        "expected_equality_branch": "left" if kind == "threshold_equal" else "not_applicable",
    }


def golden_vectors(
    count_artifact: Mapping[str, object],
    candidate_artifact: Mapping[str, object],
) -> dict[str, object]:
    validate_artifacts(count_artifact, candidate_artifact)
    heads: dict[str, object] = {}
    for head, artifact in (
        ("count", count_artifact),
        ("candidate", candidate_artifact),
    ):
        representative = _representative_split(artifact)
        cases = [
            _golden_case(artifact, representative, "realistic", None),
            _golden_case(
                artifact,
                representative,
                "threshold_below",
                float(representative["previous"]),
            ),
            _golden_case(
                artifact,
                representative,
                "threshold_equal",
                float(representative["threshold"]),
            ),
            _golden_case(
                artifact,
                representative,
                "threshold_above",
                float(representative["following"]),
            ),
        ]
        shared_or_left = "shared" if artifact["reference_pair"]["shared_values"] else "left"  # type: ignore[index]
        errors = [
            {"kind": "nan_at", "vector": shared_or_left, "index": 0},
            {"kind": "positive_infinity_at", "vector": "left", "index": 0},
            {"kind": "negative_infinity_at", "vector": "right", "index": 0},
            {"kind": "length_mismatch", "vector": "right", "drop_last": 1},
            {"kind": "duplicate_item_id"},
            {"kind": "float32_overflow_at", "vector": "left", "index": 0},
        ]
        heads[head] = {
            "feature_names_sha256": artifact["feature_names_sha256"],
            "representative": representative,
            "cases": cases,
            "errors": errors,
        }
    result = {
        "schema_version": GOLDEN_SCHEMA_VERSION,
        "input_dtype": MODEL_INPUT_DTYPE,
        "threshold_policy": THRESHOLD_POLICY,
        "retained_feature_families": list(count_artifact["retained_feature_families"]),  # type: ignore[arg-type]
        "heads": heads,
    }
    validate_golden_vectors(result, count_artifact, candidate_artifact)
    return result


def validate_golden_vectors(
    value: Mapping[str, object],
    count_artifact: Mapping[str, object],
    candidate_artifact: Mapping[str, object],
) -> None:
    validate_artifacts(count_artifact, candidate_artifact)
    golden = _mapping(
        value,
        {"schema_version", "input_dtype", "threshold_policy", "retained_feature_families", "heads"},
        "golden vectors",
    )
    if _integer(golden["schema_version"], "golden schema version") != GOLDEN_SCHEMA_VERSION:
        raise ValueError("golden schema version is unsupported")
    if golden["input_dtype"] != MODEL_INPUT_DTYPE or golden["threshold_policy"] != THRESHOLD_POLICY:
        raise ValueError("golden numeric policy mismatch")
    if golden["retained_feature_families"] != count_artifact["retained_feature_families"]:
        raise ValueError("golden retained families mismatch")
    heads = _mapping(golden["heads"], {"count", "candidate"}, "golden heads")
    for head, artifact in (("count", count_artifact), ("candidate", candidate_artifact)):
        section = _mapping(
            heads[head],
            {"feature_names_sha256", "representative", "cases", "errors"},
            f"golden {head} section",
        )
        if section["feature_names_sha256"] != artifact["feature_names_sha256"]:
            raise ValueError("golden feature hash mismatch")
        cases = section["cases"]
        errors = section["errors"]
        if not isinstance(cases, list) or [case.get("kind") if isinstance(case, Mapping) else None for case in cases] != [
            "realistic",
            "threshold_below",
            "threshold_equal",
            "threshold_above",
        ]:
            raise ValueError("golden cases are not the exact frozen set")
        if not isinstance(errors, list) or [error.get("kind") if isinstance(error, Mapping) else None for error in errors] != [
            "nan_at",
            "positive_infinity_at",
            "negative_infinity_at",
            "length_mismatch",
            "duplicate_item_id",
            "float32_overflow_at",
        ]:
            raise ValueError("golden errors are not the exact tagged set")
        regenerated_representative = _representative_split(artifact)
        if section["representative"] != regenerated_representative:
            raise ValueError("golden representative split mismatch")
        expected_cases = [
            _golden_case(artifact, regenerated_representative, "realistic", None),
            _golden_case(artifact, regenerated_representative, "threshold_below", float(regenerated_representative["previous"])),
            _golden_case(artifact, regenerated_representative, "threshold_equal", float(regenerated_representative["threshold"])),
            _golden_case(artifact, regenerated_representative, "threshold_above", float(regenerated_representative["following"])),
        ]
        if cases != expected_cases:
            raise ValueError("golden case values do not reproduce the artifacts")
        reference = artifact["reference_pair"]
        assert isinstance(reference, Mapping)
        shared_or_left = "shared" if reference["shared_values"] else "left"
        expected_errors = [
            {"kind": "nan_at", "vector": shared_or_left, "index": 0},
            {"kind": "positive_infinity_at", "vector": "left", "index": 0},
            {"kind": "negative_infinity_at", "vector": "right", "index": 0},
            {"kind": "length_mismatch", "vector": "right", "drop_last": 1},
            {"kind": "duplicate_item_id"},
            {"kind": "float32_overflow_at", "vector": "left", "index": 0},
        ]
        if errors != expected_errors:
            raise ValueError("golden tagged error payloads are not exact")
        encoded = _canonical_json_bytes(errors)
        if b"NaN" in encoded or b"Infinity" in encoded:
            raise ValueError("golden errors contain noncanonical numeric values")


def write_golden_vectors(
    path: Path,
    count_artifact: Mapping[str, object],
    candidate_artifact: Mapping[str, object],
) -> str:
    value = golden_vectors(count_artifact, candidate_artifact)
    data = _canonical_json_bytes(value)
    target = Path(path)
    temporary = _write_temporary(target, data)
    try:
        os.replace(temporary, target)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    if target.read_bytes() != data:
        raise OSError("installed golden bytes differ")
    return _sha256(data)


def verify_top_level_manifest(path: Path) -> dict[str, object]:
    """Verify the complete conventionally installed Task 13 model-only graph."""
    manifest_path = Path(path).resolve(strict=True)
    data = manifest_path.read_bytes()
    try:
        manifest = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=_duplicate_rejecting_object,
            parse_constant=_reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("top-level manifest is not valid JSON") from error
    if not isinstance(manifest, dict) or _canonical_json_bytes(manifest) != data:
        raise ValueError("top-level manifest is not canonical JSON")
    # Import lazily to avoid a module cycle while the exporter imports this file.
    from benchmark.export_factorized_ranker import validate_export_manifest
    from benchmark.factorized_ranker.folds import (
        _validate_task8_manifest,
        load_fold_manifest,
    )

    validate_export_manifest(manifest)
    model_dir = manifest_path.parent
    benchmark_dir = model_dir.parent
    repo_root = benchmark_dir.parent
    paths = {
        "corpus_manifest_sha256": model_dir
        / "cath17287_factorized_corpus_v1_manifest.json",
        "fold_manifest_sha256": model_dir / "cath17287_factorized_folds_v1.json",
        "cv_report_sha256": model_dir / "factorized_ranker_v1_cv.json",
        "ablation_report_sha256": model_dir
        / "factorized_ranker_v1_ablations.json",
        "oof_predictions_sha256": model_dir / "factorized_ranker_v1_oof.csv",
        "count_model_sha256": model_dir / "factorized_count_v1.json",
        "candidate_model_sha256": model_dir / "factorized_candidate_v1.json",
        "golden_sha256": model_dir / "factorized_ranker_v1_golden.json",
        "standalone_baseline_sha256": model_dir
        / "cath663_standalone_structural_baseline.csv",
        "generated_rust_sha256": repo_root
        / "sword2-lib/src/sword/factorized_ranker/generated_model.rs",
    }
    for field, artifact_path in paths.items():
        try:
            artifact_bytes = artifact_path.read_bytes()
        except OSError as error:
            raise ValueError(f"top-level dependency is unavailable: {artifact_path}") from error
        if _sha256(artifact_bytes) != manifest[field]:
            raise ValueError(f"top-level dependency hash mismatch for {field}")

    count_artifact = load_artifact(
        paths["count_model_sha256"],
        expected_sha256=str(manifest["count_model_sha256"]),
        expected_head="count",
    )
    candidate_artifact = load_artifact(
        paths["candidate_model_sha256"],
        expected_sha256=str(manifest["candidate_model_sha256"]),
        expected_head="candidate",
    )
    validate_artifacts(count_artifact, candidate_artifact)
    for artifact in (count_artifact, candidate_artifact):
        for artifact_field, manifest_field in (
            ("corpus_manifest_sha256", "corpus_manifest_sha256"),
            ("fold_manifest_sha256", "fold_manifest_sha256"),
            ("cv_report_sha256", "cv_report_sha256"),
            ("ablation_report_sha256", "ablation_report_sha256"),
            ("oof_predictions_sha256", "oof_predictions_sha256"),
            ("feature_schema_sha256", "feature_schema_sha256"),
            ("feature_dump_binary_sha256", "feature_dump_binary_sha256"),
            ("source_git_commit", "source_git_commit"),
            ("training_command", "training_command"),
            ("retained_feature_families", "retained_feature_families"),
        ):
            if artifact[artifact_field] != manifest[manifest_field]:
                raise ValueError(
                    f"top-level manifest disagrees with model field {artifact_field}"
                )
    corpus_bytes = paths["corpus_manifest_sha256"].read_bytes()
    corpus_payload = json.loads(
        corpus_bytes,
        object_pairs_hook=_duplicate_rejecting_object,
        parse_constant=_reject_json_constant,
    )
    if _canonical_json_bytes(corpus_payload) != corpus_bytes:
        raise ValueError("copied corpus manifest is not canonical JSON")
    corpus_manifest = _validate_task8_manifest(corpus_payload, "cath17287")
    if (
        corpus_manifest["dataset_sha256"] != manifest["dataset_sha256"]
        or corpus_manifest["binary_sha256"]
        != manifest["feature_dump_binary_sha256"]
        or corpus_manifest["feature_schema_hash"]
        != manifest["feature_schema_sha256"]
        or corpus_manifest["git_commit"] != manifest["source_git_commit"]
    ):
        raise ValueError("copied corpus manifest disagrees with top-level provenance")
    if {
        name: descriptor["sha256"]
        for name, descriptor in corpus_manifest["tables"].items()
    } != manifest["corpus_table_sha256s"]:
        raise ValueError("copied corpus table hashes disagree with top-level manifest")
    load_fold_manifest(
        paths["fold_manifest_sha256"],
        expected_sha256=str(manifest["fold_manifest_sha256"]),
        expected_dataset_sha256=str(manifest["dataset_sha256"]),
        expected_corpus_sha256=str(manifest["corpus_manifest_sha256"]),
        expected_chains_sha256=str(manifest["corpus_table_sha256s"]["chains"]),  # type: ignore[index]
    )
    for report_field in ("cv_report_sha256", "ablation_report_sha256"):
        report_bytes = paths[report_field].read_bytes()
        try:
            report = json.loads(
                report_bytes.decode("utf-8"),
                object_pairs_hook=_duplicate_rejecting_object,
            )
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
            raise ValueError(f"installed {report_field} is invalid JSON") from error
        if _canonical_json_bytes(report) != report_bytes:
            raise ValueError(f"installed {report_field} is not canonical JSON")
    golden_bytes = paths["golden_sha256"].read_bytes()
    try:
        golden = json.loads(
            golden_bytes.decode("utf-8"),
            object_pairs_hook=_duplicate_rejecting_object,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("installed golden vectors are invalid JSON") from error
    if _canonical_json_bytes(golden) != golden_bytes:
        raise ValueError("installed golden vectors are not canonical JSON")
    validate_golden_vectors(golden, count_artifact, candidate_artifact)
    if not paths["generated_rust_sha256"].read_bytes().endswith(b"\n"):
        raise ValueError("installed generated Rust lacks a trailing newline")
    return manifest
