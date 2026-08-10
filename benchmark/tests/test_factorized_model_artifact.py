from __future__ import annotations

import copy
import json
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy
import sklearn
from sklearn.ensemble import GradientBoostingClassifier

from benchmark.factorized_ranker.model_artifact import (
    ARTIFACT_SCHEMA_VERSION,
    MODEL_INPUT_DTYPE,
    THRESHOLD_POLICY,
    ModelProvenance,
    freeze_classifier,
    load_artifact,
    predict_artifact,
    raw_artifact_score,
    raw_tree_value,
    validate_artifact,
    validate_artifacts,
    write_artifact,
    write_artifacts,
)
from benchmark.factorized_ranker.training import MODEL_GRID, head_feature_spec


ARTIFACT_KEYS = {
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


def _provenance() -> ModelProvenance:
    return ModelProvenance(
        corpus_manifest_sha256="1" * 64,
        fold_manifest_sha256="2" * 64,
        cv_report_sha256="3" * 64,
        ablation_report_sha256="4" * 64,
        oof_predictions_sha256="5" * 64,
        feature_schema_sha256="6" * 64,
        feature_dump_binary_sha256="7" * 64,
        source_git_commit="8" * 40,
        training_command=(
            "benchmark.train_factorized_ranker",
            "--corpus-dir",
            "<CORPUS_DIR>",
            "--count-model-out",
            "<COUNT_MODEL_OUT>",
        ),
        python_version=platform.python_version(),
        numpy_version=np.__version__,
        pandas_version=pd.__version__,
        scipy_version=scipy.__version__,
        sklearn_version=sklearn.__version__,
    )


def _fitted(head: str = "count") -> tuple[GradientBoostingClassifier, np.ndarray, tuple[str, ...]]:
    names = head_feature_spec(head, ("base",)).pair_features
    rng = np.random.default_rng(37 if head == "count" else 73)
    x = np.ascontiguousarray(rng.normal(size=(256, len(names))), dtype=np.float64)
    y = (x[:, 0] + 0.2 * x[:, 1] > 0.0).astype(np.int8)
    params = MODEL_GRID[0]
    model = GradientBoostingClassifier(
        n_estimators=params.n_estimators,
        learning_rate=params.learning_rate,
        min_samples_leaf=params.min_samples_leaf,
        max_depth=params.max_depth,
        criterion="friedman_mse",
        random_state=37,
        loss="log_loss",
    )
    model.fit(x, y, sample_weight=np.ones(len(y), dtype=np.float64))
    return model, x, names


@pytest.fixture(scope="module")
def count_artifact() -> dict[str, object]:
    model, x, names = _fitted("count")
    return freeze_classifier(
        model,
        head_name="count",
        feature_names=names,
        retained_feature_families=("base",),
        selected_hyperparameters=MODEL_GRID[0],
        provenance=_provenance(),
        audit_x=x,
    )


@pytest.fixture(scope="module")
def candidate_artifact() -> dict[str, object]:
    model, x, names = _fitted("candidate")
    return freeze_classifier(
        model,
        head_name="candidate",
        feature_names=names,
        retained_feature_families=("base",),
        selected_hyperparameters=MODEL_GRID[0],
        provenance=_provenance(),
        audit_x=x,
    )


def test_frozen_tree_probability_matches_sklearn() -> None:
    model, x, names = _fitted("count")
    artifact = freeze_classifier(
        model,
        head_name="count",
        feature_names=names,
        retained_feature_families=("base",),
        selected_hyperparameters=MODEL_GRID[0],
        provenance=_provenance(),
        audit_x=x,
    )
    canonical = x.astype(np.float32)
    np.testing.assert_allclose(
        raw_artifact_score(artifact, x),
        model.decision_function(canonical),
        rtol=0.0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        predict_artifact(artifact, x),
        model.predict_proba(canonical)[:, 1],
        rtol=0.0,
        atol=1e-12,
    )


def test_threshold_equality_takes_left_branch() -> None:
    tree = {
        "nodes": [
            {
                "feature_index": 0,
                "threshold": 1.25,
                "left_child": 1,
                "right_child": 2,
                "leaf_value": 0.0,
            },
            {
                "feature_index": -1,
                "threshold": 0.0,
                "left_child": -1,
                "right_child": -1,
                "leaf_value": -1.0,
            },
            {
                "feature_index": -1,
                "threshold": 0.0,
                "left_child": -1,
                "right_child": -1,
                "leaf_value": 1.0,
            },
        ]
    }
    assert raw_tree_value(tree, np.asarray([1.25], dtype=np.float64)) == -1.0


def test_exact_schema_and_float32_threshold_policy(count_artifact) -> None:
    assert set(count_artifact) == ARTIFACT_KEYS
    assert count_artifact["schema_version"] == ARTIFACT_SCHEMA_VERSION == 1
    assert count_artifact["input_dtype"] == MODEL_INPUT_DTYPE == "float32"
    assert count_artifact["threshold_policy"] == THRESHOLD_POLICY == "floor_to_f32"
    assert tuple(count_artifact["feature_names"]) == head_feature_spec(
        "count", ("base",)
    ).pair_features
    for tree in count_artifact["trees"]:
        for node in tree["nodes"]:
            if node["feature_index"] >= 0:
                assert float(np.float32(node["threshold"])) == node["threshold"]
    validate_artifact(count_artifact, expected_head="count")


def test_nonfinite_and_float32_overflow_inputs_fail(count_artifact) -> None:
    width = count_artifact["feature_count"]
    for value in (np.nan, np.inf, -np.inf, np.finfo(np.float64).max):
        x = np.zeros((1, width), dtype=np.float64)
        x[0, 0] = value
        with pytest.raises(ValueError):
            predict_artifact(count_artifact, x)


@pytest.mark.parametrize(
    "change",
    [
        "invalid_child",
        "nan_threshold",
        "nan_leaf",
        "wrong_feature_count",
        "unreachable_node",
        "shared_child",
        "cycle",
        "bad_leaf_sentinel",
        "bad_internal_leaf_value",
        "noncanonical_threshold",
        "too_many_trees",
    ],
)
def test_invalid_artifacts_are_rejected(count_artifact, change: str) -> None:
    artifact = copy.deepcopy(count_artifact)
    tree = artifact["trees"][0]
    root = tree["nodes"][0]
    if change == "invalid_child":
        root["left_child"] = len(tree["nodes"])
    elif change == "nan_threshold":
        root["threshold"] = float("nan")
    elif change == "nan_leaf":
        next(node for node in tree["nodes"] if node["feature_index"] == -1)[
            "leaf_value"
        ] = float("nan")
    elif change == "wrong_feature_count":
        artifact["feature_count"] += 1
    elif change == "unreachable_node":
        tree["nodes"].append(
            {
                "feature_index": -1,
                "threshold": 0.0,
                "left_child": -1,
                "right_child": -1,
                "leaf_value": 0.0,
            }
        )
        artifact["node_count"] += 1
    elif change == "shared_child":
        root["right_child"] = root["left_child"]
    elif change == "cycle":
        root["left_child"] = 0
    elif change == "bad_leaf_sentinel":
        next(node for node in tree["nodes"] if node["feature_index"] == -1)[
            "left_child"
        ] = 0
    elif change == "bad_internal_leaf_value":
        root["leaf_value"] = 0.25
    elif change == "noncanonical_threshold":
        root["threshold"] = float(root["threshold"]) + 2.0**-40
    else:
        artifact["trees"] = artifact["trees"] * 97
        artifact["tree_count"] = len(artifact["trees"])
        artifact["n_estimators"] = len(artifact["trees"])
        artifact["node_count"] = sum(
            len(value["nodes"]) for value in artifact["trees"]
        )
    with pytest.raises(ValueError):
        validate_artifact(artifact)


def test_cross_head_provenance_is_exact(count_artifact, candidate_artifact) -> None:
    validate_artifacts(count_artifact, candidate_artifact)
    changed = copy.deepcopy(candidate_artifact)
    changed["oof_predictions_sha256"] = "a" * 64
    with pytest.raises(ValueError):
        validate_artifacts(count_artifact, changed)


def test_artifact_writes_are_canonical_and_paired(tmp_path: Path, count_artifact, candidate_artifact) -> None:
    count_path = tmp_path / "a" / "count.json"
    candidate_path = tmp_path / "a" / "candidate.json"
    hashes = write_artifacts(
        count_path, count_artifact, candidate_path, candidate_artifact
    )
    assert hashes == (
        write_artifact(tmp_path / "b" / "count.json", dict(reversed(list(count_artifact.items())))),
        write_artifact(tmp_path / "b" / "candidate.json", dict(reversed(list(candidate_artifact.items())))),
    )
    assert count_path.read_bytes() == (tmp_path / "b" / "count.json").read_bytes()
    assert candidate_path.read_bytes() == (tmp_path / "b" / "candidate.json").read_bytes()
    assert count_path.read_bytes().endswith(b"\n")
    assert load_artifact(count_path, expected_sha256=hashes[0]) == count_artifact
    noncanonical = tmp_path / "noncanonical.json"
    noncanonical.write_text(json.dumps(count_artifact, indent=2), encoding="utf-8")
    with pytest.raises(ValueError):
        load_artifact(noncanonical)


@pytest.mark.parametrize("defect", ["type", "subsample", "init", "random_state", "features", "version"])
def test_freeze_rejects_nonapproved_estimator_or_provenance(defect: str) -> None:
    model, x, names = _fitted("count")
    provenance = _provenance()
    if defect == "type":
        model = object()  # type: ignore[assignment]
    elif defect == "subsample":
        model.subsample = 0.5
    elif defect == "init":
        model.init = "zero"
    elif defect == "random_state":
        model.random_state = 38
    elif defect == "features":
        names = names[:-1]
    else:
        provenance = ModelProvenance(
            **{**provenance.__dict__, "sklearn_version": "0.0"}
        )
    with pytest.raises(ValueError):
        freeze_classifier(
            model,
            head_name="count",
            feature_names=names,
            retained_feature_families=("base",),
            selected_hyperparameters=MODEL_GRID[0],
            provenance=provenance,
            audit_x=x,
        )


@pytest.mark.parametrize("defect", ["missing", "shape", "nonfinite", "input_dependent"])
def test_freeze_rejects_unusable_private_initializer(defect: str) -> None:
    model, x, names = _fitted("count")
    if defect == "missing":
        model._raw_predict_init = None  # type: ignore[method-assign]
    elif defect == "shape":
        model._raw_predict_init = lambda values: np.zeros(1, dtype=np.float64)  # type: ignore[method-assign]
    elif defect == "nonfinite":
        model._raw_predict_init = lambda values: np.asarray([[np.nan]])  # type: ignore[method-assign]
    else:
        model._raw_predict_init = lambda values: np.asarray(  # type: ignore[method-assign]
            [[float(values[0, 0])]], dtype=np.float64
        )
    with pytest.raises(ValueError):
        freeze_classifier(
            model,
            head_name="count",
            feature_names=names,
            retained_feature_families=("base",),
            selected_hyperparameters=MODEL_GRID[0],
            provenance=_provenance(),
            audit_x=x,
        )


@pytest.mark.parametrize("defect", ["classes", "estimators", "feature_count", "tree_value"])
def test_freeze_rejects_malformed_fitted_state(defect: str) -> None:
    model, x, names = _fitted("count")
    if defect == "classes":
        model.classes_ = np.asarray([1, 0])
    elif defect == "estimators":
        model.estimators_ = model.estimators_[:-1]
    elif defect == "feature_count":
        model.n_features_in_ += 1
    else:
        tree = model.estimators_[0, 0].tree_
        leaf = int(np.flatnonzero(tree.children_left == -1)[0])
        tree.value[leaf, 0, 0] = np.nan
    with pytest.raises(ValueError):
        freeze_classifier(
            model,
            head_name="count",
            feature_names=names,
            retained_feature_families=("base",),
            selected_hyperparameters=MODEL_GRID[0],
            provenance=_provenance(),
            audit_x=x,
        )


def test_paired_validation_failure_preserves_both_existing_files(
    tmp_path: Path, count_artifact, candidate_artifact
) -> None:
    count_path = tmp_path / "count.json"
    candidate_path = tmp_path / "candidate.json"
    count_path.write_bytes(b"old count\n")
    candidate_path.write_bytes(b"old candidate\n")
    invalid = copy.deepcopy(candidate_artifact)
    invalid["seed"] = 38
    with pytest.raises(ValueError):
        write_artifacts(count_path, count_artifact, candidate_path, invalid)
    assert count_path.read_bytes() == b"old count\n"
    assert candidate_path.read_bytes() == b"old candidate\n"


def test_recorded_parity_is_full_matrix_and_finite(count_artifact) -> None:
    audit = count_artifact["parity_audit"]
    assert audit["row_count"] == 256
    assert audit["column_count"] == count_artifact["feature_count"]
    assert audit["batch_size"] <= 4096
    assert audit["max_raw_score_delta"] <= 1e-12
    assert audit["max_probability_delta"] <= 1e-12
    assert len(audit["matrix_sha256"]) == 64
