"""Grouped development CV and ordered ablations for factorized ranker heads."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import multiprocessing
import os
import platform
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from io import StringIO
from numbers import Integral, Real
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.ensemble import GradientBoostingClassifier

import benchmark.factorized_ranker.pairs as pair_api
from benchmark.datasets import CathEntry, parse_cath_domain_string
from benchmark.factorized_ranker.corpus import (
    CANDIDATE_FIELDS,
    CHAIN_FIELDS,
    COUNT_FIELDS,
    CorpusPaths,
    feature_schema_hash,
)
from benchmark.factorized_ranker.folds import (
    FoldAssignment,
    _load_accepted_entries,
    chain_label_cohorts,
    load_fold_manifest,
    validate_folds,
)
from benchmark.factorized_ranker.pairs import (
    PairBatch,
    build_candidate_pairs,
    build_count_pairs,
)
from benchmark.factorized_ranker.ranking import (
    normalized_borda,
    select_candidate,
    select_count,
    symmetrized_probability,
)
from benchmark.factorized_ranker.schema import (
    BASE_CANDIDATE_FEATURES,
    BOUNDARY_LOCAL_FEATURES,
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    COUNT_SUMMARY_SOURCES,
    DISCONTINUITY_FEATURES,
    DOMAIN_CONDITIONED_FEATURES,
    GLOBAL_FEATURES,
    RELATIVE_HIERARCHY_FEATURES,
    SCHEMA_VERSION,
    pair_feature_names,
)
from benchmark.stats import paired_chain_bootstrap


FEATURE_FAMILY_ORDER = (
    "base",
    "global_count",
    "domain_conditioned",
    "boundary_local",
    "relative_hierarchy",
    "discontinuity",
)

BASE_COUNT_SUMMARY_SOURCES = (
    "min_size",
    "max_cr",
    "density_min",
    "mean_density",
    "contact_q_mean",
    "contact_q_max",
    "n_segments",
    "n_discontinuous",
    "boundary_coil_fraction",
)

BASE_COUNT_ITEM_FEATURES = (
    "count_num_domains",
    "count_n_candidates",
    "count_candidate_fraction",
    "count_modal_distance",
    *(
        f"count_{source}_{stat}"
        for source in BASE_COUNT_SUMMARY_SOURCES
        for stat in ("min", "mean", "max")
    ),
)

if BASE_COUNT_SUMMARY_SOURCES != tuple(
    source for source in COUNT_SUMMARY_SOURCES if source != "legacy_distance"
):
    raise RuntimeError("base count summary sources disagree with the frozen schema")
if any(source not in BASE_CANDIDATE_FEATURES for source in BASE_COUNT_SUMMARY_SOURCES):
    raise RuntimeError("base count summary source is not a base candidate feature")
_COUNT_POSITIONS = {name: index for index, name in enumerate(COUNT_ITEM_FEATURES)}
if any(name not in _COUNT_POSITIONS for name in BASE_COUNT_ITEM_FEATURES) or [
    _COUNT_POSITIONS[name] for name in BASE_COUNT_ITEM_FEATURES
] != sorted(_COUNT_POSITIONS[name] for name in BASE_COUNT_ITEM_FEATURES):
    raise RuntimeError("base count features are not an ordered count-schema subsequence")


@dataclass(frozen=True)
class HeadFeatureSpec:
    shared_features: tuple[str, ...]
    item_features: tuple[str, ...]
    pair_features: tuple[str, ...]


@dataclass(frozen=True, order=True)
class Hyperparameters:
    n_estimators: Literal[64, 96]
    learning_rate: Literal[0.03, 0.05]
    min_samples_leaf: Literal[32, 64]
    max_depth: Literal[3] = 3


MODEL_GRID = tuple(
    Hyperparameters(n_estimators, learning_rate, min_samples_leaf)
    for n_estimators in (64, 96)
    for learning_rate in (0.03, 0.05)
    for min_samples_leaf in (32, 64)
)


def validate_training_jobs(jobs: int) -> int:
    if (
        isinstance(jobs, bool)
        or not isinstance(jobs, Integral)
        or not 1 <= int(jobs) <= 8
    ):
        raise ValueError("training jobs must be a non-boolean integer in 1..=8")
    return int(jobs)


def _fork_context() -> multiprocessing.context.BaseContext:
    if "fork" not in multiprocessing.get_all_start_methods():
        raise ValueError("parallel factorized training requires the fork start method")
    return multiprocessing.get_context("fork")


@dataclass(frozen=True)
class CorpusTables:
    chains: pd.DataFrame
    counts: pd.DataFrame
    candidates: pd.DataFrame


@dataclass(frozen=True)
class FoldConfigurationResult:
    fold: int
    training_chain_ids: tuple[str, ...]
    validation_chain_ids: tuple[str, ...]
    count_accuracy: float
    count_abs_error: float
    candidate_ndo: float


@dataclass(frozen=True)
class ConfigurationResult:
    retained_families: tuple[str, ...]
    count_params: Hyperparameters
    candidate_params: Hyperparameters
    count_accuracy: float
    count_abs_error: float
    candidate_ndo: float
    folds: tuple[FoldConfigurationResult, ...]


@dataclass(frozen=True)
class AblationGateReport:
    family: str
    ndo_delta: float
    ci_low: float
    count_delta: float
    contiguous_delta: float


OOF_FIELDS = (
    "chain_id",
    "fold",
    "selected_count",
    "selected_candidate_id",
    "canonical_delineation",
    "count_borda_score",
    "candidate_borda_score",
    "n_true_domains",
    "count_correct",
    "ndo",
    "boundary_f1_10",
    "matched_dice",
    "total_regret",
    "count_regret",
    "within_count_regret",
    "true_count_bin",
    "length_bin",
    "continuity_cohort",
    "label_cohort",
    "count_tie_break",
    "candidate_tie_break",
)


@dataclass(frozen=True)
class VerifiedTrainingData:
    corpus: CorpusTables
    assignments: tuple[FoldAssignment, ...]
    metadata: tuple[CathEntry, ...]
    corpus_manifest: dict[str, Any]
    dataset_sha256: str
    corpus_manifest_sha256: str
    chains_sha256: str
    fold_manifest_sha256: str


@dataclass(frozen=True)
class CountDecision:
    selected_count: int
    borda_score: float
    tie_break: Literal["none", "legacy_count", "lower_count"]


@dataclass(frozen=True)
class OOFRow:
    chain_id: str
    fold: int
    selected_count: int
    selected_candidate_id: str
    canonical_delineation: str
    count_borda_score: float
    candidate_borda_score: float
    n_true_domains: int
    count_correct: int
    ndo: float
    boundary_f1_10: float
    matched_dice: float
    total_regret: float
    count_regret: float
    within_count_regret: float
    true_count_bin: str
    length_bin: str
    continuity_cohort: str
    label_cohort: str
    count_tie_break: str
    candidate_tie_break: str


@dataclass(frozen=True)
class OOFResult:
    rows: tuple[OOFRow, ...]
    count_decisions: dict[str, CountDecision]
    csv_bytes: bytes
    sha256: str


@dataclass(frozen=True)
class _OOFFoldTask:
    fold: int
    retained_families: tuple[str, ...]
    count_params: Hyperparameters
    candidate_params: Hyperparameters
    seed: int
    reuse_count_decisions: bool


@dataclass(frozen=True)
class _OOFFoldWorkResult:
    task: _OOFFoldTask
    rows: tuple[OOFRow, ...]
    count_decisions: tuple[tuple[str, CountDecision], ...]


@dataclass(frozen=True)
class _OOFWorkerContext:
    data: VerifiedTrainingData
    reused_count_decisions: Mapping[str, CountDecision] | None


_OOF_WORKER_CONTEXT: _OOFWorkerContext | None = None


@dataclass(frozen=True)
class TrainingRun:
    count_model: GradientBoostingClassifier
    candidate_model: GradientBoostingClassifier
    count_batch: PairBatch
    candidate_batch: PairBatch
    retained_families: tuple[str, ...]
    count_params: Hyperparameters
    candidate_params: Hyperparameters
    oof: OOFResult
    cv_report: dict[str, Any]
    ablation_report: dict[str, Any]
    artifact_hashes: dict[str, str]


def validate_retained_families(retained_families: Sequence[str]) -> tuple[str, ...]:
    try:
        families = tuple(retained_families)
    except TypeError as error:
        raise ValueError("retained families must be a sequence") from error
    if not families or families[0] != "base":
        raise ValueError("base must be the first retained family")
    if any(not isinstance(family, str) for family in families):
        raise ValueError("retained family names must be strings")
    if len(set(families)) != len(families):
        raise ValueError("retained families contain duplicates")
    positions = {family: index for index, family in enumerate(FEATURE_FAMILY_ORDER)}
    if any(family not in positions for family in families):
        raise ValueError("retained families contain an unknown name")
    observed = [positions[family] for family in families]
    if observed != sorted(observed):
        raise ValueError("retained families do not preserve approved order")
    return families


def head_feature_spec(
    head: Literal["count", "candidate"],
    retained_families: Sequence[str],
) -> HeadFeatureSpec:
    families = validate_retained_families(retained_families)
    if head not in {"count", "candidate"}:
        raise ValueError("head must be count or candidate")
    has_global = "global_count" in families
    shared = GLOBAL_FEATURES if has_global else ()
    if head == "count":
        item = COUNT_ITEM_FEATURES if has_global else BASE_COUNT_ITEM_FEATURES
    else:
        additions = {
            "domain_conditioned": DOMAIN_CONDITIONED_FEATURES,
            "boundary_local": BOUNDARY_LOCAL_FEATURES,
            "relative_hierarchy": RELATIVE_HIERARCHY_FEATURES,
            "discontinuity": DISCONTINUITY_FEATURES,
        }
        item = (
            *BASE_CANDIDATE_FEATURES,
            *(
                name
                for family in FEATURE_FAMILY_ORDER
                if family in families and family in additions
                for name in additions[family]
            ),
        )
    schema = COUNT_ITEM_FEATURES if head == "count" else CANDIDATE_FEATURES
    positions = {name: index for index, name in enumerate(schema)}
    if len(set(item)) != len(item) or any(name not in positions for name in item):
        raise ValueError("derived head item features are invalid")
    if [positions[name] for name in item] != sorted(positions[name] for name in item):
        raise ValueError("derived head item features are out of order")
    pair = tuple(pair_feature_names(shared, item))
    if len(set(pair)) != len(pair):
        raise ValueError("derived pair feature names are not unique")
    return HeadFeatureSpec(tuple(shared), tuple(item), pair)


def _validate_pair_batch(batch: PairBatch) -> None:
    if not isinstance(batch, PairBatch):
        raise ValueError("fit batch must be PairBatch")
    if (
        not isinstance(batch.feature_names, tuple)
        or not batch.feature_names
        or any(not isinstance(name, str) or not name for name in batch.feature_names)
        or len(set(batch.feature_names)) != len(batch.feature_names)
    ):
        raise ValueError("fit batch feature names are invalid")
    if batch.x.dtype != np.float64 or batch.x.ndim != 2:
        raise ValueError("fit batch x must be a float64 matrix")
    row_count, feature_count = batch.x.shape
    if row_count == 0 or feature_count != len(batch.feature_names):
        raise ValueError("fit batch matrix shape is invalid")
    arrays = (
        batch.y,
        batch.sample_weight,
        batch.chain_ids,
        batch.left_ids,
        batch.right_ids,
    )
    if any(array.ndim != 1 or len(array) != row_count for array in arrays):
        raise ValueError("fit batch metadata shape is invalid")
    if batch.y.dtype != np.int8 or batch.sample_weight.dtype != np.float64:
        raise ValueError("fit batch label/weight dtype is invalid")
    if batch.chain_ids.dtype.kind != "U" or batch.left_ids.dtype.kind != "U" or batch.right_ids.dtype.kind != "U":
        raise ValueError("fit batch ID arrays must use fixed Unicode dtypes")
    if not np.isfinite(batch.x).all() or not np.isfinite(batch.sample_weight).all():
        raise ValueError("fit batch contains non-finite values")
    if np.any(batch.sample_weight < 0.0) or not float(batch.sample_weight.sum()) > 0.0:
        raise ValueError("fit batch weights are invalid")
    if set(batch.y.tolist()) != {0, 1}:
        raise ValueError("fit batch must contain both binary classes")


def fit_head(
    batch: PairBatch,
    params: Hyperparameters,
    seed: int = 37,
) -> GradientBoostingClassifier:
    _validate_pair_batch(batch)
    if type(params) is not Hyperparameters or params not in MODEL_GRID:
        raise ValueError("hyperparameters are outside the approved grid")
    if not isinstance(seed, Integral) or isinstance(seed, (bool, np.bool_)):
        raise ValueError("seed must be a non-boolean integer")
    seed_value = int(seed)
    if not 0 <= seed_value <= 2**32 - 1:
        raise ValueError("seed is outside uint32 range")
    model = GradientBoostingClassifier(
        n_estimators=params.n_estimators,
        learning_rate=params.learning_rate,
        min_samples_leaf=params.min_samples_leaf,
        max_depth=3,
        criterion="friedman_mse",
        random_state=seed_value,
        loss="log_loss",
    )
    fitted = model.fit(batch.x, batch.y, sample_weight=batch.sample_weight)
    if fitted is not model:
        raise ValueError("sklearn fit returned a different estimator")
    if not np.array_equal(model.classes_, np.asarray([0, 1])):
        raise ValueError("fitted head classes are not exact [0,1]")
    if int(model.n_features_in_) != len(batch.feature_names):
        raise ValueError("fitted head feature count mismatch")
    if model.estimators_.shape != (params.n_estimators, 1):
        raise ValueError("fitted head estimator shape mismatch")
    for estimator in model.estimators_.ravel():
        if not np.isfinite(estimator.tree_.threshold).all() or not np.isfinite(
            estimator.tree_.value
        ).all():
            raise ValueError("fitted head contains non-finite tree state")
    probabilities = model.predict_proba(batch.x[:1])
    if probabilities.shape != (1, 2) or not np.isfinite(probabilities).all():
        raise ValueError("fitted head prediction state is invalid")
    return model


def _finite(value: object, description: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{description} must be a finite real")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{description} must be finite")
    return number


def _gate_value(report: object, name: str) -> object:
    if isinstance(report, Mapping):
        if name not in report:
            raise ValueError(f"ablation report is missing {name}")
        return report[name]
    if not hasattr(report, name):
        raise ValueError(f"ablation report is missing {name}")
    return getattr(report, name)


def retain_ablation(report: AblationGateReport | Mapping[str, object]) -> bool:
    family = _gate_value(report, "family")
    if family not in FEATURE_FAMILY_ORDER[1:]:
        raise ValueError("ablation gate family is invalid")
    ndo_delta = _finite(_gate_value(report, "ndo_delta"), "NDO delta")
    ci_low = _finite(_gate_value(report, "ci_low"), "NDO CI low")
    count_delta = _finite(_gate_value(report, "count_delta"), "count delta")
    contiguous_delta = _finite(
        _gate_value(report, "contiguous_delta"), "contiguous delta"
    )
    if not (ndo_delta > 0.0 and ci_low >= 0.0):
        return False
    if family == "global_count" and count_delta < 0.0:
        return False
    if family == "discontinuity" and contiguous_delta < -0.005:
        return False
    return True


def _model_probability(model: GradientBoostingClassifier, vector: np.ndarray) -> float:
    probability = model.predict_proba(vector.reshape(1, -1))
    if probability.shape != (1, 2) or not np.isfinite(probability).all():
        raise ValueError("head returned invalid probabilities")
    value = float(probability[0, 1])
    if not 0.0 <= value <= 1.0:
        raise ValueError("head probability is outside [0,1]")
    return value


def _row_float(row: Mapping[str, Any], name: str) -> float:
    value = float(row[name])
    if not math.isfinite(value):
        raise ValueError(f"validation field {name} is non-finite")
    return value


def _row_int(row: Mapping[str, Any], name: str) -> int:
    value = _row_float(row, name)
    if not value.is_integer():
        raise ValueError(f"validation field {name} is not integral")
    return int(value)


def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    columns = list(frame.columns)
    return [
        dict(zip(columns, row, strict=True))
        for row in frame.itertuples(index=False, name=None)
    ]


def _legacy_count(candidate_rows: Sequence[Mapping[str, Any]]) -> int:
    representatives: dict[int, Mapping[str, Any]] = {}
    for row in candidate_rows:
        count = _row_int(row, "num_domains")
        source_index = _row_int(row, "source_index")
        previous = representatives.get(count)
        if previous is None or source_index < _row_int(previous, "source_index"):
            representatives[count] = row
    if not representatives:
        raise ValueError("validation chain has no candidate representatives")
    winner = max(
        representatives.values(),
        key=lambda row: (
            _row_float(row, "legacy_distance"),
            -_row_int(row, "source_index"),
            -_row_int(row, "num_domains"),
        ),
    )
    return _row_int(winner, "num_domains")


def _count_decision(
    model: GradientBoostingClassifier,
    chain_row: Mapping[str, Any],
    count_rows: Sequence[Mapping[str, Any]],
    candidate_rows: Sequence[Mapping[str, Any]],
    spec: HeadFeatureSpec,
) -> tuple[int, float, str]:
    items = {_row_int(row, "count_num_domains"): row for row in count_rows}
    if not items:
        raise ValueError("validation chain has no count rows")
    shared = [_row_float(chain_row, name) for name in spec.shared_features]

    def probability(left: int | str, right: int | str) -> float:
        if not isinstance(left, int) or not isinstance(right, int):
            raise ValueError("count Borda IDs are not integers")
        return symmetrized_probability(
            lambda vector: _model_probability(model, vector),
            shared,
            [_row_float(items[left], name) for name in spec.item_features],
            [_row_float(items[right], name) for name in spec.item_features],
        )

    scores = normalized_borda(sorted(items), probability)
    legacy = _legacy_count(candidate_rows)
    selected = select_count(scores, legacy)
    maximum = max(scores.values())
    tied = [count for count, score in scores.items() if score == maximum]
    if len(tied) == 1:
        tie_break = "none"
    elif legacy in tied and selected == legacy:
        tie_break = "legacy_count"
    else:
        tie_break = "lower_count"
    return selected, float(scores[selected]), tie_break


def _candidate_decisions(
    model: GradientBoostingClassifier,
    chain_row: Mapping[str, Any],
    candidate_rows: Sequence[Mapping[str, Any]],
    spec: HeadFeatureSpec,
) -> dict[int, tuple[Mapping[str, Any], float, str]]:
    groups: dict[int, list[Mapping[str, Any]]] = {}
    for row in candidate_rows:
        groups.setdefault(_row_int(row, "num_domains"), []).append(row)
    shared = [_row_float(chain_row, name) for name in spec.shared_features]
    decisions: dict[int, tuple[Mapping[str, Any], float, str]] = {}
    for count in sorted(groups):
        rows = groups[count]
        by_canonical = {str(row["canonical_delineation"]): row for row in rows}
        if len(by_canonical) != len(rows):
            raise ValueError("validation candidates contain duplicate canonical IDs")

        def probability(left: int | str, right: int | str) -> float:
            if not isinstance(left, str) or not isinstance(right, str):
                raise ValueError("candidate Borda IDs are not strings")
            return symmetrized_probability(
                lambda vector: _model_probability(model, vector),
                shared,
                [
                    _row_float(by_canonical[left], name)
                    for name in spec.item_features
                ],
                [
                    _row_float(by_canonical[right], name)
                    for name in spec.item_features
                ],
            )

        scores = normalized_borda(sorted(by_canonical), probability)
        selected = select_candidate(scores)
        maximum = max(scores.values())
        tied = [candidate for candidate, score in scores.items() if score == maximum]
        tie_break = "lexical" if len(tied) > 1 else "none"
        decisions[count] = (by_canonical[selected], float(scores[selected]), tie_break)
    if not decisions:
        raise ValueError("validation chain has no candidate groups")
    return decisions


def _subset(frame: pd.DataFrame, chain_ids: set[str]) -> pd.DataFrame:
    return frame[frame["chain_id"].isin(chain_ids)].copy().reset_index(drop=True)


def _validate_cv_inputs(
    corpus: CorpusTables, assignments: Sequence[FoldAssignment]
) -> tuple[FoldAssignment, ...]:
    if not isinstance(corpus, CorpusTables):
        raise ValueError("corpus must be CorpusTables")
    if (
        tuple(corpus.chains.columns) != CHAIN_FIELDS
        or tuple(corpus.counts.columns) != COUNT_FIELDS
        or tuple(corpus.candidates.columns) != CANDIDATE_FIELDS
    ):
        raise ValueError("in-memory corpus headers do not match Task 8")
    # Validate the complete population without constructing target-labelled pairs;
    # actual pair batches are built only after each explicit fold split.
    validated_chains = pair_api._validate_chains(corpus.chains)
    pair_api._validate_counts(corpus.counts, validated_chains)
    pair_api._validate_candidates(corpus.candidates, validated_chains)
    canonical_assignments = tuple(sorted(assignments, key=lambda value: value.chain_id))
    validate_folds(canonical_assignments)
    chain_ids = set(str(value) for value in corpus.chains["chain_id"].tolist())
    if chain_ids != {assignment.chain_id for assignment in canonical_assignments}:
        raise ValueError("corpus/fold chain sets disagree")
    if {assignment.fold for assignment in canonical_assignments} != set(range(5)):
        raise ValueError("training requires exactly five populated folds")
    return canonical_assignments


def cross_validate_configuration(
    corpus: CorpusTables,
    assignments: Sequence[FoldAssignment],
    count_params: Hyperparameters,
    candidate_params: Hyperparameters,
    retained_families: Sequence[str],
    *,
    seed: int = 37,
) -> ConfigurationResult:
    families = validate_retained_families(retained_families)
    canonical_assignments = _validate_cv_inputs(corpus, assignments)
    count_spec = head_feature_spec("count", families)
    candidate_spec = head_feature_spec("candidate", families)
    fold_results: list[FoldConfigurationResult] = []
    all_count_correct: list[float] = []
    all_count_errors: list[float] = []
    all_candidate_ndo: list[float] = []
    all_ids = {assignment.chain_id for assignment in canonical_assignments}

    for fold in range(5):
        validation_ids = {
            assignment.chain_id
            for assignment in canonical_assignments
            if assignment.fold == fold
        }
        training_ids = all_ids - validation_ids
        training_chains = _subset(corpus.chains, training_ids)
        count_batch = build_count_pairs(
            training_chains,
            _subset(corpus.counts, training_ids),
            shared_features=count_spec.shared_features,
            item_features=count_spec.item_features,
            seed=seed,
        )
        candidate_batch = build_candidate_pairs(
            training_chains,
            _subset(corpus.candidates, training_ids),
            shared_features=candidate_spec.shared_features,
            item_features=candidate_spec.item_features,
            seed=seed,
        )
        count_model = fit_head(count_batch, count_params, seed)
        candidate_model = fit_head(candidate_batch, candidate_params, seed)

        fold_correct: list[float] = []
        fold_errors: list[float] = []
        fold_candidate_ndo: list[float] = []
        for chain_id in sorted(validation_ids):
            chain_rows = _records(_subset(corpus.chains, {chain_id}))
            count_rows = _records(_subset(corpus.counts, {chain_id}))
            candidate_rows = _records(_subset(corpus.candidates, {chain_id}))
            if len(chain_rows) != 1:
                raise ValueError("validation chain row is missing or duplicated")
            selected_count, _score, _tie = _count_decision(
                count_model,
                chain_rows[0],
                count_rows,
                candidate_rows,
                count_spec,
            )
            true_count = _row_int(chain_rows[0], "n_true_domains")
            fold_correct.append(float(selected_count == true_count))
            fold_errors.append(float(abs(selected_count - true_count)))
            decisions = _candidate_decisions(
                candidate_model, chain_rows[0], candidate_rows, candidate_spec
            )
            fold_candidate_ndo.append(
                sum(_row_float(row, "ndo") for row, _score, _tie in decisions.values())
                / len(decisions)
            )
        if not fold_correct or not fold_candidate_ndo:
            raise ValueError("validation fold is empty")
        all_count_correct.extend(fold_correct)
        all_count_errors.extend(fold_errors)
        all_candidate_ndo.extend(fold_candidate_ndo)
        fold_results.append(
            FoldConfigurationResult(
                fold,
                tuple(sorted(training_ids)),
                tuple(sorted(validation_ids)),
                float(np.mean(fold_correct)),
                float(np.mean(fold_errors)),
                float(np.mean(fold_candidate_ndo)),
            )
        )
    result = ConfigurationResult(
        families,
        count_params,
        candidate_params,
        _finite(float(np.mean(all_count_correct)), "count accuracy"),
        _finite(float(np.mean(all_count_errors)), "count absolute error"),
        _finite(float(np.mean(all_candidate_ndo)), "candidate NDO"),
        tuple(fold_results),
    )
    return result


@dataclass(frozen=True)
class HeadFoldResult:
    fold: int
    validation_chain_ids: tuple[str, ...]
    primary: float
    secondary: float | None


@dataclass(frozen=True)
class _HeadFoldTask:
    head: Literal["count", "candidate"]
    params: Hyperparameters
    retained_families: tuple[str, ...]
    fold: int
    seed: int


@dataclass(frozen=True)
class _HeadFoldWorkResult:
    task: _HeadFoldTask
    fold_result: HeadFoldResult
    primary_values: tuple[float, ...]
    secondary_values: tuple[float, ...]


@dataclass(frozen=True)
class _GridWorkerContext:
    corpus: CorpusTables
    assignments: tuple[FoldAssignment, ...]


_GRID_WORKER_CONTEXT: _GridWorkerContext | None = None


@dataclass(frozen=True)
class HeadGridEvaluation:
    head: Literal["count", "candidate"]
    params: Hyperparameters
    primary: float
    secondary: float | None
    folds: tuple[HeadFoldResult, ...]


@dataclass(frozen=True)
class GridSelection:
    head: Literal["count", "candidate"]
    selected: Hyperparameters
    evaluations: tuple[HeadGridEvaluation, ...]


@dataclass(frozen=True)
class TrainingStageState:
    completed_family_index: int
    retained_families: tuple[str, ...]
    count_params: Hyperparameters
    candidate_params: Hyperparameters
    retained_oof: OOFResult
    stage_records: tuple[dict[str, object], ...]
    grid_history: tuple[dict[str, object], ...]


class TrainingCheckpointStore:
    """Canonical exact-context checkpoints for individual grid/fold fits."""

    _MANIFEST_KEYS = {"schema_version", "kind", "context", "context_sha256"}
    _FOLD_KEYS = {
        "schema_version",
        "kind",
        "context",
        "context_sha256",
        "stage_family",
        "retained_feature_families",
        "head",
        "hyperparameters",
        "fold",
        "seed",
        "training_chain_ids",
        "training_chain_id_sha256",
        "validation_chain_ids",
        "validation_chain_id_sha256",
        "pair_feature_names",
        "pair_feature_names_sha256",
        "result",
    }
    _STATE_KEYS = {
        "schema_version",
        "kind",
        "context",
        "context_sha256",
        "completed_family_index",
        "retained_feature_families",
        "count_hyperparameters",
        "candidate_hyperparameters",
        "retained_oof",
        "stage_records",
        "grid_history",
    }
    _OOF_STATE_KEYS = {
        "sha256",
        "chain_count",
        "chain_id_sha256",
        "rows",
        "count_decisions",
    }

    def __init__(self, root: Path, context: Mapping[str, object]):
        if not isinstance(context, Mapping) or not context or any(
            not isinstance(key, str) or not key for key in context
        ):
            raise ValueError("checkpoint context must be a nonempty string-key mapping")
        self.root = Path(root)
        self.grid_dir = self.root / "grid"
        self.context = dict(context)
        context_bytes = _canonical_json_bytes(self.context)
        self.context_sha256 = _sha256_bytes(context_bytes)
        manifest = {
            "schema_version": 1,
            "kind": "factorized_training_checkpoints",
            "context": self.context,
            "context_sha256": self.context_sha256,
        }
        manifest_bytes = _canonical_json_bytes(manifest)
        manifest_path = self.root / "checkpoint_manifest.json"
        if manifest_path.exists():
            observed = self._read_json(manifest_path, self._MANIFEST_KEYS)
            if observed != manifest or manifest_path.read_bytes() != manifest_bytes:
                raise ValueError("checkpoint context disagrees with existing manifest")
        else:
            self.root.mkdir(parents=True, exist_ok=True)
            self._atomic_write(manifest_path, manifest_bytes)
        self.grid_dir.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _read_json(path: Path, keys: set[str]) -> dict[str, object]:
        data = Path(path).read_bytes()
        try:
            value = json.loads(data)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ValueError(f"checkpoint {path.name} is invalid JSON") from error
        if not isinstance(value, dict) or set(value) != keys:
            raise ValueError(f"checkpoint {path.name} schema mismatch")
        if _canonical_json_bytes(value) != data:
            raise ValueError(f"checkpoint {path.name} is not canonical JSON")
        return value

    @staticmethod
    def _atomic_write(path: Path, data: bytes) -> None:
        temporary = path.with_name(f".{path.name}.tmp")
        # A prior interrupted write is task-owned and explicitly incomplete.
        temporary.unlink(missing_ok=True)
        try:
            with temporary.open("xb") as handle:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
            if temporary.read_bytes() != data:
                raise ValueError("checkpoint temporary bytes changed after writing")
            os.replace(temporary, path)
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise

    def _fold_payload(
        self,
        *,
        stage_family: str,
        retained_families: tuple[str, ...],
        head: Literal["count", "candidate"],
        params: Hyperparameters,
        fold: int,
        seed: int,
        training_ids: set[str],
        validation_ids: set[str],
        pair_feature_names: tuple[str, ...],
    ) -> dict[str, object]:
        if stage_family not in FEATURE_FAMILY_ORDER:
            raise ValueError("checkpoint stage family is invalid")
        if stage_family not in retained_families:
            raise ValueError("checkpoint stage family is absent from proposed families")
        training = tuple(sorted(training_ids))
        validation = tuple(sorted(validation_ids))
        if not training or not validation or set(training) & set(validation):
            raise ValueError("checkpoint train/validation IDs are invalid")
        feature_bytes = _canonical_json_bytes(list(pair_feature_names))
        return {
            "schema_version": 1,
            "kind": "grid_fold",
            "context": self.context,
            "context_sha256": self.context_sha256,
            "stage_family": stage_family,
            "retained_feature_families": list(retained_families),
            "head": head,
            "hyperparameters": _params_payload(params),
            "fold": fold,
            "seed": seed,
            "training_chain_ids": list(training),
            "training_chain_id_sha256": _id_hash(training),
            "validation_chain_ids": list(validation),
            "validation_chain_id_sha256": _id_hash(validation),
            "pair_feature_names": list(pair_feature_names),
            "pair_feature_names_sha256": _sha256_bytes(feature_bytes),
        }

    def _fold_path(self, expected: Mapping[str, object]) -> Path:
        identifier = _sha256_bytes(_canonical_json_bytes(expected))
        return self.grid_dir / f"{identifier}.json"

    def load_fold(
        self,
        **kwargs: Any,
    ) -> tuple[HeadFoldResult, tuple[float, ...], tuple[float, ...]] | None:
        expected = self._fold_payload(**kwargs)
        path = self._fold_path(expected)
        if not path.exists():
            return None
        payload = self._read_json(path, self._FOLD_KEYS)
        for name, value in expected.items():
            if payload[name] != value:
                raise ValueError(f"grid checkpoint field {name} mismatch")
        result = payload["result"]
        if not isinstance(result, dict) or set(result) != {
            "primary",
            "secondary",
            "primary_values",
            "secondary_values",
        }:
            raise ValueError("grid checkpoint result schema mismatch")
        primary = _finite(result["primary"], "checkpoint fold primary")
        secondary_raw = result["secondary"]
        secondary = (
            None
            if secondary_raw is None
            else _finite(secondary_raw, "checkpoint fold secondary")
        )
        primary_raw = result["primary_values"]
        secondary_values_raw = result["secondary_values"]
        if not isinstance(primary_raw, list) or not isinstance(
            secondary_values_raw, list
        ):
            raise ValueError("grid checkpoint per-chain values are malformed")
        primary_values = tuple(
            _finite(value, "checkpoint primary value") for value in primary_raw
        )
        secondary_values = tuple(
            _finite(value, "checkpoint secondary value")
            for value in secondary_values_raw
        )
        validation_ids = tuple(expected["validation_chain_ids"])
        if len(primary_values) != len(validation_ids):
            raise ValueError("grid checkpoint primary population mismatch")
        if expected["head"] == "count":
            if secondary is None or len(secondary_values) != len(validation_ids):
                raise ValueError("count checkpoint secondary population mismatch")
        elif secondary is not None or secondary_values:
            raise ValueError("candidate checkpoint has unexpected secondary values")
        if primary != float(np.mean(primary_values)):
            raise ValueError("grid checkpoint primary aggregate mismatch")
        if secondary is not None and secondary != float(np.mean(secondary_values)):
            raise ValueError("grid checkpoint secondary aggregate mismatch")
        return (
            HeadFoldResult(
                int(expected["fold"]),
                validation_ids,
                primary,
                secondary,
            ),
            primary_values,
            secondary_values,
        )

    def write_fold(
        self,
        result: HeadFoldResult,
        primary_values: Sequence[float],
        secondary_values: Sequence[float],
        **kwargs: Any,
    ) -> None:
        expected = self._fold_payload(**kwargs)
        validated_primary = tuple(
            _finite(value, "checkpoint primary value") for value in primary_values
        )
        validated_secondary = tuple(
            _finite(value, "checkpoint secondary value")
            for value in secondary_values
        )
        if (
            result.fold != expected["fold"]
            or list(result.validation_chain_ids)
            != expected["validation_chain_ids"]
            or len(validated_primary) != len(result.validation_chain_ids)
            or result.primary != float(np.mean(validated_primary))
        ):
            raise ValueError("grid checkpoint result identity mismatch")
        if expected["head"] == "count":
            if (
                result.secondary is None
                or len(validated_secondary) != len(result.validation_chain_ids)
                or result.secondary != float(np.mean(validated_secondary))
            ):
                raise ValueError("count grid checkpoint secondary result mismatch")
        elif result.secondary is not None or validated_secondary:
            raise ValueError("candidate grid checkpoint has unexpected secondary values")
        payload = {
            **expected,
            "result": {
                "primary": _finite(result.primary, "checkpoint primary"),
                "secondary": (
                    None
                    if result.secondary is None
                    else _finite(result.secondary, "checkpoint secondary")
                ),
                "primary_values": list(validated_primary),
                "secondary_values": list(validated_secondary),
            },
        }
        path = self._fold_path(expected)
        data = _canonical_json_bytes(payload)
        if path.exists():
            if path.read_bytes() != data:
                raise ValueError("existing grid checkpoint bytes disagree")
            return
        self._atomic_write(path, data)

    @staticmethod
    def _checkpoint_params(value: object, description: str) -> Hyperparameters:
        if not isinstance(value, Mapping) or set(value) != {
            "n_estimators",
            "learning_rate",
            "min_samples_leaf",
            "max_depth",
        }:
            raise ValueError(f"{description} checkpoint hyperparameters are malformed")
        if any(
            type(value[name]) is not int
            for name in ("n_estimators", "min_samples_leaf", "max_depth")
        ):
            raise ValueError(
                f"{description} checkpoint integer hyperparameters are malformed"
            )
        if type(value["learning_rate"]) not in {int, float}:
            raise ValueError(f"{description} checkpoint learning rate is malformed")
        params = Hyperparameters(
            value["n_estimators"],  # type: ignore[arg-type]
            float(value["learning_rate"]),  # type: ignore[arg-type]
            value["min_samples_leaf"],  # type: ignore[arg-type]
            value["max_depth"],  # type: ignore[arg-type]
        )
        if params not in MODEL_GRID:
            raise ValueError(
                f"{description} checkpoint hyperparameters are outside the grid"
            )
        return params

    @staticmethod
    def _oof_state(oof: OOFResult, accepted_ids: Sequence[str]) -> dict[str, object]:
        accepted = tuple(sorted(accepted_ids))
        if not accepted or len(set(accepted)) != len(accepted):
            raise ValueError("stage checkpoint accepted IDs are invalid")
        ordered_rows = tuple(sorted(oof.rows, key=lambda row: row.chain_id))
        csv_bytes = _render_oof(ordered_rows)
        if (
            csv_bytes != oof.csv_bytes
            or _sha256_bytes(csv_bytes) != oof.sha256
            or tuple(row.chain_id for row in ordered_rows) != accepted
            or set(oof.count_decisions) != set(accepted)
        ):
            raise ValueError("stage checkpoint OOF population or hash mismatch")
        decisions: list[dict[str, object]] = []
        for row in ordered_rows:
            decision = oof.count_decisions[row.chain_id]
            if (
                decision.selected_count != row.selected_count
                or decision.borda_score != row.count_borda_score
                or decision.tie_break != row.count_tie_break
            ):
                raise ValueError("stage checkpoint count decision disagrees with OOF")
            decisions.append(
                {
                    "chain_id": row.chain_id,
                    "selected_count": decision.selected_count,
                    "borda_score": _finite(
                        decision.borda_score, "checkpoint count Borda score"
                    ),
                    "tie_break": decision.tie_break,
                }
            )
        return {
            "sha256": oof.sha256,
            "chain_count": len(accepted),
            "chain_id_sha256": _id_hash(accepted),
            "rows": [_oof_mapping(row) for row in ordered_rows],
            "count_decisions": decisions,
        }

    @classmethod
    def _restore_oof_state(
        cls, value: object, accepted_ids: Sequence[str]
    ) -> OOFResult:
        if not isinstance(value, Mapping) or set(value) != cls._OOF_STATE_KEYS:
            raise ValueError("stage checkpoint OOF schema mismatch")
        accepted = tuple(sorted(accepted_ids))
        if (
            type(value["chain_count"]) is not int
            or value["chain_count"] != len(accepted)
            or value["chain_id_sha256"] != _id_hash(accepted)
        ):
            raise ValueError("stage checkpoint OOF population mismatch")
        raw_rows = value["rows"]
        raw_decisions = value["count_decisions"]
        if not isinstance(raw_rows, list) or not isinstance(raw_decisions, list):
            raise ValueError("stage checkpoint OOF rows are malformed")
        rows: list[OOFRow] = []
        for raw in raw_rows:
            if not isinstance(raw, Mapping) or set(raw) != set(OOF_FIELDS):
                raise ValueError("stage checkpoint OOF row schema mismatch")
            rows.append(OOFRow(**raw))  # type: ignore[arg-type]
        ordered = tuple(sorted(rows, key=lambda row: row.chain_id))
        csv_bytes = _render_oof(ordered)
        sha256 = _sha256_bytes(csv_bytes)
        if (
            tuple(row.chain_id for row in ordered) != accepted
            or value["sha256"] != sha256
        ):
            raise ValueError("stage checkpoint OOF hash or chain IDs mismatch")
        decisions: dict[str, CountDecision] = {}
        for raw in raw_decisions:
            if not isinstance(raw, Mapping) or set(raw) != {
                "chain_id",
                "selected_count",
                "borda_score",
                "tie_break",
            }:
                raise ValueError("stage checkpoint count decision schema mismatch")
            chain_id = raw["chain_id"]
            selected_count = raw["selected_count"]
            tie_break = raw["tie_break"]
            if (
                not isinstance(chain_id, str)
                or not chain_id
                or chain_id in decisions
                or type(selected_count) is not int
                or selected_count < 1
                or tie_break not in {"none", "legacy_count", "lower_count"}
            ):
                raise ValueError("stage checkpoint count decision is malformed")
            decisions[chain_id] = CountDecision(
                selected_count,
                _finite(raw["borda_score"], "checkpoint count Borda score"),
                tie_break,  # type: ignore[arg-type]
            )
        if set(decisions) != set(accepted):
            raise ValueError("stage checkpoint count decisions have wrong population")
        for row in ordered:
            decision = decisions[row.chain_id]
            if (
                decision.selected_count != row.selected_count
                or decision.borda_score != row.count_borda_score
                or decision.tie_break != row.count_tie_break
            ):
                raise ValueError("stage checkpoint count decision disagrees with OOF")
        return OOFResult(ordered, decisions, csv_bytes, sha256)

    @staticmethod
    def _validate_stage_lists(
        completed_family_index: int,
        retained_families: tuple[str, ...],
        stage_records: Sequence[Mapping[str, object]],
        grid_history: Sequence[Mapping[str, object]],
    ) -> None:
        expected_families = FEATURE_FAMILY_ORDER[: completed_family_index + 1]
        if (
            len(stage_records) != len(expected_families)
            or len(grid_history) != len(expected_families)
            or tuple(record.get("family") for record in stage_records)
            != expected_families
            or tuple(record.get("family") for record in grid_history)
            != expected_families
            or any(type(record.get("retained")) is not bool for record in stage_records)
        ):
            raise ValueError("stage checkpoint history is incomplete or out of order")
        observed_retained = tuple(
            str(record["family"])
            for record in stage_records
            if record["retained"] is True
        )
        if observed_retained != retained_families:
            raise ValueError("stage checkpoint retained-family history mismatch")

    def write_stage(
        self,
        *,
        completed_family_index: int,
        retained_families: Sequence[str],
        count_params: Hyperparameters,
        candidate_params: Hyperparameters,
        retained_oof: OOFResult,
        stage_records: Sequence[Mapping[str, object]],
        grid_history: Sequence[Mapping[str, object]],
        accepted_ids: Sequence[str],
    ) -> None:
        if (
            type(completed_family_index) is not int
            or not 0 <= completed_family_index < len(FEATURE_FAMILY_ORDER)
        ):
            raise ValueError("completed checkpoint family index is invalid")
        families = validate_retained_families(retained_families)
        if count_params not in MODEL_GRID or candidate_params not in MODEL_GRID:
            raise ValueError("stage checkpoint hyperparameters are outside the grid")
        records = [dict(record) for record in stage_records]
        grids = [dict(record) for record in grid_history]
        self._validate_stage_lists(
            completed_family_index, families, records, grids
        )
        payload = {
            "schema_version": 1,
            "kind": "stage_state",
            "context": self.context,
            "context_sha256": self.context_sha256,
            "completed_family_index": completed_family_index,
            "retained_feature_families": list(families),
            "count_hyperparameters": _params_payload(count_params),
            "candidate_hyperparameters": _params_payload(candidate_params),
            "retained_oof": self._oof_state(retained_oof, accepted_ids),
            "stage_records": records,
            "grid_history": grids,
        }
        path = self.root / "stage_state.json"
        data = _canonical_json_bytes(payload)
        if path.exists():
            current = self._read_json(path, self._STATE_KEYS)
            current_index = current["completed_family_index"]
            if type(current_index) is not int:
                raise ValueError("existing stage checkpoint index is malformed")
            if current_index == completed_family_index:
                if path.read_bytes() != data:
                    raise ValueError("existing completed stage checkpoint disagrees")
                return
            if current_index > completed_family_index:
                raise ValueError("refusing to replace a later completed stage checkpoint")
        self._atomic_write(path, data)

    def load_stage(self, accepted_ids: Sequence[str]) -> TrainingStageState | None:
        path = self.root / "stage_state.json"
        if not path.exists():
            return None
        payload = self._read_json(path, self._STATE_KEYS)
        if (
            payload["schema_version"] != 1
            or payload["kind"] != "stage_state"
            or payload["context"] != self.context
            or payload["context_sha256"] != self.context_sha256
        ):
            raise ValueError("stage checkpoint context mismatch")
        completed = payload["completed_family_index"]
        if type(completed) is not int or not 0 <= completed < len(FEATURE_FAMILY_ORDER):
            raise ValueError("stage checkpoint completed-family index is invalid")
        raw_families = payload["retained_feature_families"]
        if not isinstance(raw_families, list):
            raise ValueError("stage checkpoint retained families are malformed")
        families = validate_retained_families(raw_families)
        raw_records = payload["stage_records"]
        raw_grids = payload["grid_history"]
        if (
            not isinstance(raw_records, list)
            or not isinstance(raw_grids, list)
            or any(not isinstance(value, Mapping) for value in (*raw_records, *raw_grids))
        ):
            raise ValueError("stage checkpoint histories are malformed")
        records = tuple(dict(value) for value in raw_records)
        grids = tuple(dict(value) for value in raw_grids)
        self._validate_stage_lists(completed, families, records, grids)
        return TrainingStageState(
            completed,
            families,
            self._checkpoint_params(payload["count_hyperparameters"], "count"),
            self._checkpoint_params(
                payload["candidate_hyperparameters"], "candidate"
            ),
            self._restore_oof_state(payload["retained_oof"], accepted_ids),
            records,
            grids,
        )


def _head_fold_populations(
    assignments: tuple[FoldAssignment, ...], fold: int
) -> tuple[set[str], tuple[str, ...]]:
    if type(fold) is not int or not 0 <= fold < 5:
        raise ValueError("head fold is invalid")
    all_ids = {assignment.chain_id for assignment in assignments}
    validation = tuple(
        sorted(
            assignment.chain_id
            for assignment in assignments
            if assignment.fold == fold
        )
    )
    if not validation:
        raise ValueError("head validation fold is empty")
    training = all_ids - set(validation)
    if not training or training & set(validation):
        raise ValueError("head train/validation populations are invalid")
    return training, validation


def _evaluate_head_fold(
    corpus: CorpusTables,
    assignments: tuple[FoldAssignment, ...],
    task: _HeadFoldTask,
) -> _HeadFoldWorkResult:
    if task.head not in {"count", "candidate"}:
        raise ValueError("head must be count or candidate")
    if task.params not in MODEL_GRID:
        raise ValueError("head fold hyperparameters are outside the frozen grid")
    if task.seed != 37 or isinstance(task.seed, bool):
        raise ValueError("head fold seed is frozen at 37")
    families = validate_retained_families(task.retained_families)
    training_ids, validation_ids = _head_fold_populations(assignments, task.fold)
    spec = head_feature_spec(task.head, families)
    training_chains = _subset(corpus.chains, training_ids)
    if task.head == "count":
        batch = build_count_pairs(
            training_chains,
            _subset(corpus.counts, training_ids),
            shared_features=spec.shared_features,
            item_features=spec.item_features,
            seed=task.seed,
        )
    else:
        batch = build_candidate_pairs(
            training_chains,
            _subset(corpus.candidates, training_ids),
            shared_features=spec.shared_features,
            item_features=spec.item_features,
            seed=task.seed,
        )
    model = fit_head(batch, task.params, task.seed)
    fold_primary: list[float] = []
    fold_secondary: list[float] = []
    for chain_id in validation_ids:
        chain_rows = _records(_subset(corpus.chains, {chain_id}))
        candidate_rows = _records(_subset(corpus.candidates, {chain_id}))
        if len(chain_rows) != 1:
            raise ValueError("validation chain row is missing or duplicated")
        if task.head == "count":
            count_rows = _records(_subset(corpus.counts, {chain_id}))
            selected, _score, _tie = _count_decision(
                model, chain_rows[0], count_rows, candidate_rows, spec
            )
            truth = _row_int(chain_rows[0], "n_true_domains")
            fold_primary.append(float(selected == truth))
            fold_secondary.append(float(abs(selected - truth)))
        else:
            decisions = _candidate_decisions(
                model, chain_rows[0], candidate_rows, spec
            )
            fold_primary.append(
                sum(
                    _row_float(row, "ndo")
                    for row, _score, _tie in decisions.values()
                )
                / len(decisions)
            )
    primary_values = tuple(fold_primary)
    secondary_values = tuple(fold_secondary)
    fold_result = HeadFoldResult(
        task.fold,
        validation_ids,
        _finite(float(np.mean(primary_values)), "fold primary objective"),
        (
            _finite(float(np.mean(secondary_values)), "fold secondary objective")
            if task.head == "count"
            else None
        ),
    )
    return _HeadFoldWorkResult(
        task,
        fold_result,
        primary_values,
        secondary_values,
    )


def _run_head_fold_task(task: _HeadFoldTask) -> _HeadFoldWorkResult:
    context = _GRID_WORKER_CONTEXT
    if context is None:
        raise RuntimeError("grid worker context is unavailable")
    return _evaluate_head_fold(context.corpus, context.assignments, task)


def _validate_head_fold_work_result(
    expected_task: _HeadFoldTask,
    expected_validation_ids: tuple[str, ...],
    result: _HeadFoldWorkResult,
) -> _HeadFoldWorkResult:
    if not isinstance(result, _HeadFoldWorkResult) or result.task != expected_task:
        raise ValueError("grid worker task identity mismatch")
    fold_result = result.fold_result
    if (
        fold_result.fold != expected_task.fold
        or fold_result.validation_chain_ids != expected_validation_ids
        or len(result.primary_values) != len(expected_validation_ids)
    ):
        raise ValueError("grid worker fold population mismatch")
    primary_values = tuple(
        _finite(value, "grid worker primary value") for value in result.primary_values
    )
    if fold_result.primary != float(np.mean(primary_values)):
        raise ValueError("grid worker primary aggregate mismatch")
    if expected_task.head == "count":
        secondary_values = tuple(
            _finite(value, "grid worker secondary value")
            for value in result.secondary_values
        )
        if (
            fold_result.secondary is None
            or len(secondary_values) != len(expected_validation_ids)
            or fold_result.secondary != float(np.mean(secondary_values))
        ):
            raise ValueError("count grid worker secondary aggregate mismatch")
    elif fold_result.secondary is not None or result.secondary_values:
        raise ValueError("candidate grid worker has unexpected secondary values")
    return result


def _assemble_head_evaluation(
    head: Literal["count", "candidate"],
    params: Hyperparameters,
    fold_results: Sequence[_HeadFoldWorkResult],
) -> HeadGridEvaluation:
    ordered = tuple(sorted(fold_results, key=lambda value: value.task.fold))
    if (
        len(ordered) != 5
        or {result.task.fold for result in ordered} != set(range(5))
        or any(result.task.head != head or result.task.params != params for result in ordered)
    ):
        raise ValueError("head grid fold results are incomplete")
    primary_values = tuple(
        value for result in ordered for value in result.primary_values
    )
    secondary_values = tuple(
        value for result in ordered for value in result.secondary_values
    )
    primary = _finite(float(np.mean(primary_values)), "head primary objective")
    secondary = (
        _finite(float(np.mean(secondary_values)), "head secondary objective")
        if head == "count"
        else None
    )
    return HeadGridEvaluation(
        head,
        params,
        primary,
        secondary,
        tuple(result.fold_result for result in ordered),
    )


def _assemble_grid_evaluations(
    head: Literal["count", "candidate"],
    results: Mapping[_HeadFoldTask, _HeadFoldWorkResult],
) -> tuple[HeadGridEvaluation, ...]:
    if not results:
        raise ValueError("head grid results are empty")
    first_task = next(iter(results))
    expected_tasks = tuple(
        _HeadFoldTask(
            head,
            params,
            first_task.retained_families,
            fold,
            first_task.seed,
        )
        for params in MODEL_GRID
        for fold in range(5)
    )
    if set(results) != set(expected_tasks):
        raise ValueError("head grid task results are incomplete or unexpected")
    return tuple(
        _assemble_head_evaluation(
            head,
            params,
            tuple(
                results[
                    _HeadFoldTask(
                        head,
                        params,
                        first_task.retained_families,
                        fold,
                        first_task.seed,
                    )
                ]
                for fold in range(5)
            ),
        )
        for params in MODEL_GRID
    )


def _evaluate_head_configuration(
    corpus: CorpusTables,
    assignments: Sequence[FoldAssignment],
    head: Literal["count", "candidate"],
    params: Hyperparameters,
    retained_families: Sequence[str],
    *,
    seed: int = 37,
    checkpoint_store: TrainingCheckpointStore | None = None,
    stage_family: str | None = None,
) -> HeadGridEvaluation:
    families = validate_retained_families(retained_families)
    canonical_assignments = _validate_cv_inputs(corpus, assignments)
    if (checkpoint_store is None) != (stage_family is None):
        raise ValueError("checkpoint store and stage family must be supplied together")
    spec = head_feature_spec(head, families)
    work_results: list[_HeadFoldWorkResult] = []
    for fold in range(5):
        task = _HeadFoldTask(head, params, families, fold, seed)
        training_ids, validation_ids = _head_fold_populations(
            canonical_assignments, fold
        )
        checkpoint_kwargs = {
            "stage_family": stage_family,
            "retained_families": families,
            "head": head,
            "params": params,
            "fold": fold,
            "seed": seed,
            "training_ids": training_ids,
            "validation_ids": set(validation_ids),
            "pair_feature_names": spec.pair_features,
        }
        cached = (
            checkpoint_store.load_fold(**checkpoint_kwargs)
            if checkpoint_store is not None
            else None
        )
        if cached is None:
            work = _evaluate_head_fold(corpus, canonical_assignments, task)
        else:
            fold_result, primary_values, secondary_values = cached
            work = _HeadFoldWorkResult(
                task, fold_result, primary_values, secondary_values
            )
        work = _validate_head_fold_work_result(task, validation_ids, work)
        if cached is None and checkpoint_store is not None:
            checkpoint_store.write_fold(
                work.fold_result,
                work.primary_values,
                work.secondary_values,
                **checkpoint_kwargs,
            )
        work_results.append(work)
    return _assemble_head_evaluation(head, params, work_results)


def _evaluate_head_grid_parallel(
    corpus: CorpusTables,
    assignments: Sequence[FoldAssignment],
    head: Literal["count", "candidate"],
    retained_families: Sequence[str],
    *,
    seed: int,
    jobs: int,
    mp_context: multiprocessing.context.BaseContext,
    checkpoint_store: TrainingCheckpointStore | None,
    stage_family: str | None,
) -> tuple[HeadGridEvaluation, ...]:
    global _GRID_WORKER_CONTEXT

    families = validate_retained_families(retained_families)
    canonical_assignments = _validate_cv_inputs(corpus, assignments)
    if (checkpoint_store is None) != (stage_family is None):
        raise ValueError("checkpoint store and stage family must be supplied together")
    spec = head_feature_spec(head, families)
    tasks = tuple(
        _HeadFoldTask(head, params, families, fold, seed)
        for params in MODEL_GRID
        for fold in range(5)
    )
    results: dict[_HeadFoldTask, _HeadFoldWorkResult] = {}
    checkpoint_kwargs: dict[_HeadFoldTask, dict[str, object]] = {}
    cached_count = 0
    for task in tasks:
        training_ids, validation_ids = _head_fold_populations(
            canonical_assignments, task.fold
        )
        kwargs: dict[str, object] = {
            "stage_family": stage_family,
            "retained_families": families,
            "head": head,
            "params": task.params,
            "fold": task.fold,
            "seed": seed,
            "training_ids": training_ids,
            "validation_ids": set(validation_ids),
            "pair_feature_names": spec.pair_features,
        }
        checkpoint_kwargs[task] = kwargs
        cached = (
            checkpoint_store.load_fold(**kwargs)
            if checkpoint_store is not None
            else None
        )
        if cached is not None:
            fold_result, primary_values, secondary_values = cached
            work = _HeadFoldWorkResult(
                task, fold_result, primary_values, secondary_values
            )
            results[task] = _validate_head_fold_work_result(
                task, validation_ids, work
            )
            cached_count += 1

    missing = tuple(task for task in tasks if task not in results)
    if missing:
        if _GRID_WORKER_CONTEXT is not None:
            raise RuntimeError("grid worker context is already active")
        _GRID_WORKER_CONTEXT = _GridWorkerContext(corpus, canonical_assignments)
        executor = ProcessPoolExecutor(max_workers=jobs, mp_context=mp_context)
        futures = {}
        fitted_count = 0
        try:
            futures = {
                executor.submit(_run_head_fold_task, task): task for task in missing
            }
            try:
                for future in as_completed(futures):
                    expected_task = futures[future]
                    work = future.result()
                    _training_ids, validation_ids = _head_fold_populations(
                        canonical_assignments, expected_task.fold
                    )
                    work = _validate_head_fold_work_result(
                        expected_task, validation_ids, work
                    )
                    if work.task in results:
                        raise ValueError("grid worker returned a duplicate task")
                    if checkpoint_store is not None:
                        checkpoint_store.write_fold(
                            work.fold_result,
                            work.primary_values,
                            work.secondary_values,
                            **checkpoint_kwargs[work.task],
                        )
                    results[work.task] = work
                    fitted_count += 1
                    print(
                        "factorized-training grid "
                        f"head={head} family={stage_family or families[-1]} "
                        f"completed={len(results)}/40 cached={cached_count} "
                        f"fitted={fitted_count}",
                        file=sys.stderr,
                        flush=True,
                    )
            except BaseException:
                for future in futures:
                    future.cancel()
                executor.shutdown(wait=True, cancel_futures=True)
                raise
            else:
                executor.shutdown(wait=True)
        finally:
            _GRID_WORKER_CONTEXT = None
    return _assemble_grid_evaluations(head, results)


def select_head_hyperparameters(
    corpus: CorpusTables,
    assignments: Sequence[FoldAssignment],
    head: Literal["count", "candidate"],
    retained_families: Sequence[str],
    *,
    seed: int = 37,
    jobs: int = 1,
    checkpoint_store: TrainingCheckpointStore | None = None,
    stage_family: str | None = None,
) -> GridSelection:
    jobs = validate_training_jobs(jobs)
    if head not in {"count", "candidate"}:
        raise ValueError("head must be count or candidate")
    resume_kwargs: dict[str, object] = {}
    if checkpoint_store is not None or stage_family is not None:
        resume_kwargs = {
            "checkpoint_store": checkpoint_store,
            "stage_family": stage_family,
        }
    if jobs == 1:
        evaluations = tuple(
            _evaluate_head_configuration(
                corpus,
                assignments,
                head,
                params,
                retained_families,
                seed=seed,
                **resume_kwargs,  # type: ignore[arg-type]
            )
            for params in MODEL_GRID
        )
    else:
        evaluations = _evaluate_head_grid_parallel(
            corpus,
            assignments,
            head,
            retained_families,
            seed=seed,
            jobs=jobs,
            mp_context=_fork_context(),
            checkpoint_store=checkpoint_store,
            stage_family=stage_family,
        )
    if len(evaluations) != len(MODEL_GRID):
        raise ValueError("head grid evaluation count mismatch")
    if head == "count":
        selected_evaluation = min(
            evaluations,
            key=lambda evaluation: (
                -evaluation.primary,
                evaluation.secondary,
                evaluation.params.n_estimators,
                -evaluation.params.min_samples_leaf,
                evaluation.params.learning_rate,
            ),
        )
    else:
        selected_evaluation = min(
            evaluations,
            key=lambda evaluation: (
                -evaluation.primary,
                evaluation.params.n_estimators,
                -evaluation.params.min_samples_leaf,
                evaluation.params.learning_rate,
            ),
        )
    return GridSelection(head, selected_evaluation.params, evaluations)


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _read_normalized_csv(path: Path, expected_fields: Sequence[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, keep_default_na=False)
    if tuple(frame.columns) != tuple(expected_fields):
        raise ValueError(f"normalized table {path.name} header mismatch")
    return frame


def load_verified_training_data(
    corpus_dir: Path,
    fold_manifest_path: Path,
) -> VerifiedTrainingData:
    corpus_dir = Path(corpus_dir)
    fold_manifest_path = Path(fold_manifest_path)
    metadata, dataset_hash, corpus_hash, chains_hash = _load_accepted_entries(
        corpus_dir, "cath17287"
    )
    corpus_paths = CorpusPaths.at(corpus_dir)
    manifest_bytes = corpus_paths.manifest.read_bytes()
    if _sha256_bytes(manifest_bytes) != corpus_hash:
        raise ValueError("corpus manifest hash changed after verification")
    corpus_manifest = json.loads(manifest_bytes)
    fold_bytes = fold_manifest_path.read_bytes()
    fold_hash = _sha256_bytes(fold_bytes)
    fold_manifest = load_fold_manifest(
        fold_manifest_path,
        expected_sha256=fold_hash,
        expected_dataset_sha256=dataset_hash,
        expected_corpus_sha256=corpus_hash,
        expected_chains_sha256=chains_hash,
    )
    chains = _read_normalized_csv(corpus_paths.chains, CHAIN_FIELDS)
    counts = _read_normalized_csv(corpus_paths.counts, COUNT_FIELDS)
    candidates = _read_normalized_csv(corpus_paths.candidates, CANDIDATE_FIELDS)
    corpus = CorpusTables(chains, counts, candidates)
    assignments = tuple(fold_manifest.assignments)
    _validate_cv_inputs(corpus, assignments)
    chain_ids = set(str(value) for value in chains["chain_id"].tolist())
    count_chain_ids = set(str(value) for value in counts["chain_id"].tolist())
    candidate_chain_ids = set(str(value) for value in candidates["chain_id"].tolist())
    if chain_ids != count_chain_ids or chain_ids != candidate_chain_ids:
        raise ValueError("normalized table chain populations disagree")
    count_groups = {
        (str(row.chain_id), int(row.count_num_domains))
        for row in counts[["chain_id", "count_num_domains"]].itertuples(index=False)
    }
    candidate_groups = {
        (str(row.chain_id), int(row.num_domains))
        for row in candidates[["chain_id", "num_domains"]].itertuples(index=False)
    }
    if not candidate_groups.issubset(count_groups):
        raise ValueError("candidate count groups are not covered by complete count rows")
    # _validate_cv_inputs has applied Task 10's complete row/join validators but
    # deliberately has not constructed any target-labelled pre-split pairs.
    if {entry.entry_id for entry in metadata} != chain_ids:
        raise ValueError("verified metadata/accepted chain populations disagree")
    return VerifiedTrainingData(
        corpus,
        tuple(sorted(assignments, key=lambda value: value.chain_id)),
        tuple(sorted(metadata, key=lambda value: value.entry_id)),
        corpus_manifest,
        dataset_hash,
        corpus_hash,
        chains_hash,
        fold_hash,
    )


_OOF_INTEGER_FIELDS = {
    "fold",
    "selected_count",
    "n_true_domains",
    "count_correct",
}
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


def _oof_mapping(row: OOFRow) -> dict[str, object]:
    return {field: getattr(row, field) for field in OOF_FIELDS}


def _render_oof(rows: Sequence[OOFRow]) -> bytes:
    ordered = sorted(rows, key=lambda row: row.chain_id)
    if len({row.chain_id for row in ordered}) != len(ordered):
        raise ValueError("OOF rows contain duplicate chain IDs")
    output = StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=OOF_FIELDS, lineterminator="\n")
    writer.writeheader()
    for row in ordered:
        values = _oof_mapping(row)
        formatted: dict[str, str] = {}
        for field in OOF_FIELDS:
            value = values[field]
            if field in _OOF_INTEGER_FIELDS:
                if not isinstance(value, Integral) or isinstance(value, bool):
                    raise ValueError(f"OOF integer field {field} is invalid")
                formatted[field] = str(int(value))
            elif field in _OOF_FLOAT_FIELDS:
                formatted[field] = format(_finite(value, f"OOF {field}"), ".17g")
            else:
                if not isinstance(value, str) or not value or "\0" in value:
                    raise ValueError(f"OOF string field {field} is invalid")
                formatted[field] = value
        writer.writerow(formatted)
    return output.getvalue().encode("utf-8")


def _truth_continuity(entry: CathEntry) -> str:
    domains = parse_cath_domain_string(entry.chopping)
    if len(domains) != entry.n_domains or not domains:
        raise ValueError("truth topology disagrees with verified domain count")
    if any(not segments for segments in domains):
        raise ValueError("truth topology contains an empty domain")
    return "contiguous" if all(len(segments) == 1 for segments in domains) else "discontinuous"


def _generate_oof_fold(
    data: VerifiedTrainingData,
    task: _OOFFoldTask,
    reused_count_decisions: Mapping[str, CountDecision] | None,
) -> _OOFFoldWorkResult:
    if task.count_params not in MODEL_GRID or task.candidate_params not in MODEL_GRID:
        raise ValueError("OOF fold hyperparameters are outside the frozen grid")
    if task.seed != 37 or isinstance(task.seed, bool):
        raise ValueError("OOF fold seed is frozen at 37")
    if type(task.reuse_count_decisions) is not bool or task.reuse_count_decisions != (
        reused_count_decisions is not None
    ):
        raise ValueError("OOF fold reused-count identity mismatch")
    families = validate_retained_families(task.retained_families)
    count_spec = head_feature_spec("count", families)
    candidate_spec = head_feature_spec("candidate", families)
    assignments = data.assignments
    all_ids = {assignment.chain_id for assignment in assignments}
    metadata = {entry.entry_id: entry for entry in data.metadata}
    assignment_by_id = {
        assignment.chain_id: assignment for assignment in assignments
    }
    if set(metadata) != all_ids or set(assignment_by_id) != all_ids:
        raise ValueError("OOF metadata/fold populations disagree")
    if reused_count_decisions is not None and set(reused_count_decisions) != all_ids:
        raise ValueError("reused count decisions do not cover exact accepted chains")
    training_ids, validation_ids = _head_fold_populations(assignments, task.fold)
    training_chains = _subset(data.corpus.chains, training_ids)
    if reused_count_decisions is None:
        count_batch = build_count_pairs(
            training_chains,
            _subset(data.corpus.counts, training_ids),
            shared_features=count_spec.shared_features,
            item_features=count_spec.item_features,
            seed=task.seed,
        )
        count_model: GradientBoostingClassifier | None = fit_head(
            count_batch, task.count_params, task.seed
        )
    else:
        count_model = None
    candidate_batch = build_candidate_pairs(
        training_chains,
        _subset(data.corpus.candidates, training_ids),
        shared_features=candidate_spec.shared_features,
        item_features=candidate_spec.item_features,
        seed=task.seed,
    )
    candidate_model = fit_head(candidate_batch, task.candidate_params, task.seed)
    label_cohorts = chain_label_cohorts(assignments, task.fold)
    rows: list[OOFRow] = []
    count_decisions: dict[str, CountDecision] = {}
    for chain_id in validation_ids:
        chain_values = _records(_subset(data.corpus.chains, {chain_id}))
        count_values = _records(_subset(data.corpus.counts, {chain_id}))
        candidate_values = _records(_subset(data.corpus.candidates, {chain_id}))
        if len(chain_values) != 1:
            raise ValueError("OOF chain row is missing or duplicated")
        if reused_count_decisions is None:
            assert count_model is not None
            selected_count, count_score, count_tie = _count_decision(
                count_model,
                chain_values[0],
                count_values,
                candidate_values,
                count_spec,
            )
            decision = CountDecision(selected_count, count_score, count_tie)
        else:
            decision = reused_count_decisions[chain_id]
        count_decisions[chain_id] = decision
        candidate_decisions = _candidate_decisions(
            candidate_model,
            chain_values[0],
            candidate_values,
            candidate_spec,
        )
        if decision.selected_count not in candidate_decisions:
            raise ValueError("selected count has no candidate group")
        winner, candidate_score, candidate_tie = candidate_decisions[
            decision.selected_count
        ]
        winner_id = str(winner["candidate_id"])
        canonical = str(winner["canonical_delineation"])
        matching = [
            row
            for row in candidate_values
            if str(row["candidate_id"]) == winner_id
            and str(row["canonical_delineation"]) == canonical
        ]
        if len(matching) != 1 or matching[0] is not winner:
            raise ValueError("selected candidate identity does not match normalized row")
        chosen_ndo = _row_float(winner, "ndo")
        best_all = max(_row_float(row, "ndo") for row in candidate_values)
        selected_count_rows = [
            row
            for row in candidate_values
            if _row_int(row, "num_domains") == decision.selected_count
        ]
        if not selected_count_rows:
            raise ValueError("selected count has no normalized candidates")
        best_selected = max(_row_float(row, "ndo") for row in selected_count_rows)
        total_regret = best_all - chosen_ndo
        count_regret = best_all - best_selected
        within_regret = best_selected - chosen_ndo
        if (
            any(
                not math.isfinite(value) or value < 0.0
                for value in (total_regret, count_regret, within_regret)
            )
            or abs(total_regret - (count_regret + within_regret)) > 1e-12
        ):
            raise ValueError("OOF regret decomposition is invalid")
        assignment = assignment_by_id[chain_id]
        true_count = _row_int(chain_values[0], "n_true_domains")
        if true_count != metadata[chain_id].n_domains:
            raise ValueError("OOF truth count disagrees with verified metadata")
        if chain_id not in label_cohorts:
            raise ValueError("OOF chain label cohort is missing")
        rows.append(
            OOFRow(
                chain_id=chain_id,
                fold=task.fold,
                selected_count=decision.selected_count,
                selected_candidate_id=winner_id,
                canonical_delineation=canonical,
                count_borda_score=decision.borda_score,
                candidate_borda_score=candidate_score,
                n_true_domains=true_count,
                count_correct=int(decision.selected_count == true_count),
                ndo=chosen_ndo,
                boundary_f1_10=_row_float(winner, "boundary_f1_10"),
                matched_dice=_row_float(winner, "matched_dice"),
                total_regret=total_regret,
                count_regret=count_regret,
                within_count_regret=within_regret,
                true_count_bin=assignment.true_count_bin,
                length_bin=assignment.length_bin,
                continuity_cohort=_truth_continuity(metadata[chain_id]),
                label_cohort=label_cohorts[chain_id],
                count_tie_break=decision.tie_break,
                candidate_tie_break=candidate_tie,
            )
        )
    return _OOFFoldWorkResult(
        task,
        tuple(sorted(rows, key=lambda row: row.chain_id)),
        tuple(sorted(count_decisions.items())),
    )


def _run_oof_fold_task(task: _OOFFoldTask) -> _OOFFoldWorkResult:
    context = _OOF_WORKER_CONTEXT
    if context is None:
        raise RuntimeError("OOF worker context is unavailable")
    return _generate_oof_fold(
        data=context.data,
        task=task,
        reused_count_decisions=context.reused_count_decisions,
    )


def _validate_oof_fold_work_result(
    data: VerifiedTrainingData,
    expected_task: _OOFFoldTask,
    result: _OOFFoldWorkResult,
    *,
    reused_count_decisions: Mapping[str, CountDecision] | None = None,
) -> _OOFFoldWorkResult:
    if not isinstance(result, _OOFFoldWorkResult) or result.task != expected_task:
        raise ValueError("OOF worker task identity mismatch")
    _training_ids, validation_ids = _head_fold_populations(
        data.assignments, expected_task.fold
    )
    row_ids = tuple(row.chain_id for row in result.rows)
    decision_ids = tuple(chain_id for chain_id, _decision in result.count_decisions)
    if (
        row_ids != validation_ids
        or decision_ids != validation_ids
        or any(row.fold != expected_task.fold for row in result.rows)
    ):
        raise ValueError("OOF worker fold population mismatch")
    decisions = dict(result.count_decisions)
    if len(decisions) != len(result.count_decisions):
        raise ValueError("OOF worker returned duplicate count decisions")
    if expected_task.reuse_count_decisions != (reused_count_decisions is not None):
        raise ValueError("OOF worker reused count decision mode mismatch")
    for row in result.rows:
        decision = decisions[row.chain_id]
        if (
            not isinstance(decision, CountDecision)
            or row.selected_count != decision.selected_count
            or row.count_borda_score != decision.borda_score
            or row.count_tie_break != decision.tie_break
        ):
            raise ValueError("OOF worker row/count decision mismatch")
        if (
            reused_count_decisions is not None
            and decision != reused_count_decisions[row.chain_id]
        ):
            raise ValueError("OOF worker changed a reused count decision")
    _render_oof(result.rows)
    return result


def _combine_oof_fold_results(
    data: VerifiedTrainingData,
    tasks: Sequence[_OOFFoldTask],
    results: Mapping[_OOFFoldTask, _OOFFoldWorkResult],
    *,
    reused_count_decisions: Mapping[str, CountDecision] | None = None,
) -> OOFResult:
    ordered_tasks = tuple(sorted(tasks, key=lambda task: task.fold))
    if (
        len(ordered_tasks) != 5
        or {task.fold for task in ordered_tasks} != set(range(5))
        or set(results) != set(ordered_tasks)
    ):
        raise ValueError("OOF fold task results are incomplete or unexpected")
    rows: list[OOFRow] = []
    count_decisions: dict[str, CountDecision] = {}
    for task in ordered_tasks:
        result = _validate_oof_fold_work_result(
            data,
            task,
            results[task],
            reused_count_decisions=reused_count_decisions,
        )
        rows.extend(result.rows)
        for chain_id, decision in result.count_decisions:
            if chain_id in count_decisions:
                raise ValueError("OOF count decisions contain duplicate chain IDs")
            count_decisions[chain_id] = decision
    all_ids = {assignment.chain_id for assignment in data.assignments}
    ordered = tuple(sorted(rows, key=lambda row: row.chain_id))
    if {row.chain_id for row in ordered} != all_ids or len(ordered) != len(all_ids):
        raise ValueError("OOF rows do not cover exact accepted chains")
    if set(count_decisions) != all_ids:
        raise ValueError("OOF count decisions do not cover exact accepted chains")
    csv_bytes = _render_oof(ordered)
    return OOFResult(ordered, count_decisions, csv_bytes, _sha256_bytes(csv_bytes))


def generate_oof(
    data: VerifiedTrainingData,
    retained_families: Sequence[str],
    count_params: Hyperparameters,
    candidate_params: Hyperparameters,
    *,
    seed: int = 37,
    jobs: int = 1,
    reused_count_decisions: Mapping[str, CountDecision] | None = None,
) -> OOFResult:
    global _OOF_WORKER_CONTEXT

    jobs = validate_training_jobs(jobs)
    families = validate_retained_families(retained_families)
    assignments = data.assignments
    all_ids = {assignment.chain_id for assignment in assignments}
    metadata = {entry.entry_id: entry for entry in data.metadata}
    if set(metadata) != all_ids or len(assignments) != len(all_ids):
        raise ValueError("OOF metadata/fold populations disagree")
    if reused_count_decisions is not None and set(reused_count_decisions) != all_ids:
        raise ValueError("reused count decisions do not cover exact accepted chains")
    tasks = tuple(
        _OOFFoldTask(
            fold,
            families,
            count_params,
            candidate_params,
            seed,
            reused_count_decisions is not None,
        )
        for fold in range(5)
    )
    results: dict[_OOFFoldTask, _OOFFoldWorkResult] = {}
    if jobs == 1:
        for task in tasks:
            result = _generate_oof_fold(data, task, reused_count_decisions)
            results[task] = _validate_oof_fold_work_result(
                data,
                task,
                result,
                reused_count_decisions=reused_count_decisions,
            )
    else:
        mp_context = _fork_context()
        if _OOF_WORKER_CONTEXT is not None:
            raise RuntimeError("OOF worker context is already active")
        _OOF_WORKER_CONTEXT = _OOFWorkerContext(data, reused_count_decisions)
        executor = ProcessPoolExecutor(
            max_workers=min(jobs, len(tasks)), mp_context=mp_context
        )
        futures = {}
        try:
            futures = {
                executor.submit(_run_oof_fold_task, task): task for task in tasks
            }
            try:
                for future in as_completed(futures):
                    expected_task = futures[future]
                    result = _validate_oof_fold_work_result(
                        data,
                        expected_task,
                        future.result(),
                        reused_count_decisions=reused_count_decisions,
                    )
                    if result.task in results:
                        raise ValueError("OOF worker returned a duplicate task")
                    results[result.task] = result
                    print(
                        "factorized-training oof "
                        f"family={families[-1]} completed={len(results)}/5",
                        file=sys.stderr,
                        flush=True,
                    )
            except BaseException:
                for future in futures:
                    future.cancel()
                executor.shutdown(wait=True, cancel_futures=True)
                raise
            else:
                executor.shutdown(wait=True)
        finally:
            _OOF_WORKER_CONTEXT = None
    return _combine_oof_fold_results(
        data,
        tasks,
        results,
        reused_count_decisions=reused_count_decisions,
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


def _id_hash(chain_ids: Sequence[str]) -> str:
    ordered = sorted(chain_ids)
    if len(set(ordered)) != len(ordered) or any(not value or "\0" in value for value in ordered):
        raise ValueError("cohort chain IDs are invalid")
    return hashlib.sha256("\0".join(ordered).encode("utf-8")).hexdigest()


def _params_payload(params: Hyperparameters) -> dict[str, object]:
    return {
        "n_estimators": params.n_estimators,
        "learning_rate": params.learning_rate,
        "min_samples_leaf": params.min_samples_leaf,
        "max_depth": params.max_depth,
    }


def _spec_payload(spec: HeadFeatureSpec) -> dict[str, object]:
    payload = {
        "shared_features": list(spec.shared_features),
        "item_features": list(spec.item_features),
        "pair_features": list(spec.pair_features),
    }
    payload["sha256"] = _sha256_bytes(_canonical_json_bytes(payload).rstrip(b"\n"))
    return payload


def _grid_payload(selection: GridSelection) -> dict[str, object]:
    return {
        "head": selection.head,
        "selected": _params_payload(selection.selected),
        "evaluations": [
            {
                "params": _params_payload(evaluation.params),
                "primary": evaluation.primary,
                "secondary": evaluation.secondary,
                "folds": [
                    {
                        "fold": fold.fold,
                        "validation_chain_ids": list(fold.validation_chain_ids),
                        "validation_chain_id_sha256": _id_hash(
                            fold.validation_chain_ids
                        ),
                        "primary": fold.primary,
                        "secondary": fold.secondary,
                    }
                    for fold in evaluation.folds
                ],
            }
            for evaluation in selection.evaluations
        ],
    }


def _cohort_summary(rows: Sequence[OOFRow], chain_ids: Sequence[str]) -> dict[str, object]:
    ids = sorted(chain_ids)
    if not ids:
        return {"status": "unavailable", "n": 0}
    by_id = {row.chain_id: row for row in rows}
    if any(chain_id not in by_id for chain_id in ids):
        raise ValueError("cohort contains a chain absent from OOF")
    selected = [by_id[chain_id] for chain_id in ids]
    return {
        "status": "available",
        "n": len(ids),
        "chain_id_sha256": _id_hash(ids),
        "mean_ndo": _finite(float(np.mean([row.ndo for row in selected])), "cohort NDO"),
        "count_accuracy": _finite(
            float(np.mean([row.count_correct for row in selected])),
            "cohort count accuracy",
        ),
        "mean_boundary_f1_10": _finite(
            float(np.mean([row.boundary_f1_10 for row in selected])),
            "cohort boundary F1",
        ),
        "mean_matched_dice": _finite(
            float(np.mean([row.matched_dice for row in selected])),
            "cohort matched Dice",
        ),
        "mean_total_regret": _finite(
            float(np.mean([row.total_regret for row in selected])),
            "cohort total regret",
        ),
        "mean_count_regret": _finite(
            float(np.mean([row.count_regret for row in selected])),
            "cohort count regret",
        ),
        "mean_within_count_regret": _finite(
            float(np.mean([row.within_count_regret for row in selected])),
            "cohort within-count regret",
        ),
    }


def _all_cohorts(oof: OOFResult) -> dict[str, object]:
    rows = oof.rows
    reports: dict[str, object] = {
        "overall": _cohort_summary(rows, [row.chain_id for row in rows]),
        "fold": {},
        "true_count_bin": {},
        "length_bin": {},
        "continuity": {},
        "label": {},
    }
    group_specs = (
        ("fold", tuple(str(value) for value in range(5)), lambda row: str(row.fold)),
        ("true_count_bin", ("1", "2", "3", "4", "5+"), lambda row: row.true_count_bin),
        (
            "length_bin",
            ("<250", "250-349", "350-449", "450+"),
            lambda row: row.length_bin,
        ),
        (
            "continuity",
            ("contiguous", "discontinuous"),
            lambda row: row.continuity_cohort,
        ),
        ("label", ("seen", "unseen", "unknown"), lambda row: row.label_cohort),
    )
    for report_name, values, key in group_specs:
        report = reports[report_name]
        assert isinstance(report, dict)
        for value in values:
            report[value] = _cohort_summary(
                rows, [row.chain_id for row in rows if key(row) == value]
            )
    return reports


def _cohort_delta(
    previous: OOFResult,
    proposed: OOFResult,
    cohort: str,
) -> dict[str, object]:
    previous_by_id = {row.chain_id: row for row in previous.rows}
    proposed_by_id = {row.chain_id: row for row in proposed.rows}
    ids = sorted(
        chain_id
        for chain_id, row in previous_by_id.items()
        if row.continuity_cohort == cohort
    )
    if not ids:
        return {"status": "unavailable", "n": 0}
    if any(proposed_by_id[chain_id].continuity_cohort != cohort for chain_id in ids):
        raise ValueError("truth continuity cohort changed between stages")
    delta = float(
        np.mean(
            [
                proposed_by_id[chain_id].ndo - previous_by_id[chain_id].ndo
                for chain_id in ids
            ]
        )
    )
    return {
        "status": "available",
        "n": len(ids),
        "chain_id_sha256": _id_hash(ids),
        "mean_ndo_delta": _finite(delta, f"{cohort} NDO delta"),
        "diagnostic_threshold": -0.005,
        "passes_diagnostic": delta >= -0.005,
    }


def _reason(report: AblationGateReport, retained: bool) -> str:
    if retained:
        return "retained_all_declared_gates"
    if report.ndo_delta <= 0.0:
        return "rejected_nonpositive_ndo_delta"
    if report.ci_low < 0.0:
        return "rejected_negative_ndo_ci_low"
    if report.family == "global_count" and report.count_delta < 0.0:
        return "rejected_negative_count_accuracy_delta"
    if report.family == "discontinuity" and report.contiguous_delta < -0.005:
        return "rejected_contiguous_guard"
    raise ValueError("ablation rejection has no stable reason")


def _common_report(
    data: VerifiedTrainingData,
    normalized_command: Sequence[str],
    accepted_ids: Sequence[str],
) -> dict[str, object]:
    command = list(normalized_command)
    if not command or any(not isinstance(token, str) for token in command):
        raise ValueError("normalized training command is invalid")
    if any(token.startswith("/") or "=/" in token for token in command):
        raise ValueError("normalized training command contains a physical path")
    family_mapping = {
        family: {
            head: _spec_payload(head_feature_spec(head, tuple(
                candidate
                for candidate in FEATURE_FAMILY_ORDER
                if candidate == "base" or (
                    candidate in FEATURE_FAMILY_ORDER[1 : FEATURE_FAMILY_ORDER.index(family) + 1]
                )
            )))
            for head in ("count", "candidate")
        }
        for family in FEATURE_FAMILY_ORDER
    }
    return {
        "schema_version": 1,
        "seed": 37,
        "evidence_role": "development_model_selection_not_locked_acceptance",
        "normalized_command": command,
        "dataset": "cath17287",
        "dataset_sha256": data.dataset_sha256,
        "corpus_manifest_sha256": data.corpus_manifest_sha256,
        "chains_sha256": data.chains_sha256,
        "fold_manifest_sha256": data.fold_manifest_sha256,
        "feature_schema_version": SCHEMA_VERSION,
        "feature_schema_sha256": feature_schema_hash(),
        "binary_sha256": data.corpus_manifest["binary_sha256"],
        "git_commit": data.corpus_manifest["git_commit"],
        "accepted_chain_count": len(accepted_ids),
        "accepted_chain_id_sha256": _id_hash(accepted_ids),
        "feature_family_order": list(FEATURE_FAMILY_ORDER),
        "feature_family_mapping": family_mapping,
        "model_grid": [_params_payload(params) for params in MODEL_GRID],
        "versions": {
            "python_implementation": platform.python_implementation(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
    }


_TRAINING_SOURCE_FILES = (
    "benchmark/train_factorized_ranker.py",
    "benchmark/export_factorized_ranker.py",
    "benchmark/factorized_ranker/training.py",
    "benchmark/factorized_ranker/model_artifact.py",
    "benchmark/factorized_ranker/pairs.py",
    "benchmark/factorized_ranker/ranking.py",
    "benchmark/factorized_ranker/schema.py",
    "benchmark/factorized_ranker/corpus.py",
    "benchmark/factorized_ranker/folds.py",
    "benchmark/datasets.py",
    "benchmark/stats.py",
)


def build_training_checkpoint_context(
    data: VerifiedTrainingData,
    normalized_command: Sequence[str],
) -> dict[str, object]:
    """Bind resumable work to the exact code, data, folds, and invocation."""

    command = list(normalized_command)
    if not command or any(not isinstance(token, str) or not token for token in command):
        raise ValueError("checkpoint training command is invalid")
    if any(token.startswith("/") or "=/" in token for token in command):
        raise ValueError("checkpoint training command contains a physical path")
    repo_root = Path(__file__).resolve().parents[2]
    source_files: dict[str, str] = {}
    for relative in _TRAINING_SOURCE_FILES:
        path = repo_root / relative
        try:
            source_files[relative] = _sha256_bytes(path.read_bytes())
        except OSError as error:
            raise ValueError(f"training source file is unavailable: {relative}") from error
    git_result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        capture_output=True,
        text=True,
        check=False,
    )
    source_commit = git_result.stdout.strip()
    if (
        git_result.returncode != 0
        or len(source_commit) != 40
        or any(character not in "0123456789abcdef" for character in source_commit)
    ):
        raise ValueError("training source Git commit is unavailable")
    accepted_ids = tuple(
        sorted(assignment.chain_id for assignment in data.assignments)
    )
    source_bytes = _canonical_json_bytes(source_files)
    return {
        "schema_version": 1,
        "kind": "factorized_training_resume_context",
        "source_git_commit": source_commit,
        "source_files": source_files,
        "training_source_sha256": _sha256_bytes(source_bytes),
        "dataset_sha256": data.dataset_sha256,
        "corpus_manifest_sha256": data.corpus_manifest_sha256,
        "chains_sha256": data.chains_sha256,
        "fold_manifest_sha256": data.fold_manifest_sha256,
        "feature_schema_sha256": feature_schema_hash(),
        "feature_family_order": list(FEATURE_FAMILY_ORDER),
        "model_grid": [_params_payload(params) for params in MODEL_GRID],
        "accepted_chain_count": len(accepted_ids),
        "accepted_chain_id_sha256": _id_hash(accepted_ids),
        "seed": 37,
        "normalized_command": command,
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
    }


def _atomic_write(path: Path, data: bytes) -> None:
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


def run_grouped_training(
    data: VerifiedTrainingData,
    out_dir: Path,
    normalized_command: Sequence[str],
    *,
    seed: int = 37,
    jobs: int = 1,
    checkpoint_store: TrainingCheckpointStore | None = None,
) -> TrainingRun:
    if seed != 37 or isinstance(seed, bool):
        raise ValueError("grouped training seed is frozen at 37")
    jobs = validate_training_jobs(jobs)
    accepted_ids = tuple(sorted(assignment.chain_id for assignment in data.assignments))
    common = _common_report(data, normalized_command, accepted_ids)
    restored = (
        checkpoint_store.load_stage(accepted_ids)
        if checkpoint_store is not None
        else None
    )
    if restored is None:
        retained = ("base",)
        count_grid = select_head_hyperparameters(
            data.corpus,
            data.assignments,
            "count",
            retained,
            seed=seed,
            jobs=jobs,
            checkpoint_store=checkpoint_store,
            stage_family="base" if checkpoint_store is not None else None,
        )
        candidate_grid = select_head_hyperparameters(
            data.corpus,
            data.assignments,
            "candidate",
            retained,
            seed=seed,
            jobs=jobs,
            checkpoint_store=checkpoint_store,
            stage_family="base" if checkpoint_store is not None else None,
        )
        retained_count_params = count_grid.selected
        retained_candidate_params = candidate_grid.selected
        retained_oof = generate_oof(
            data,
            retained,
            retained_count_params,
            retained_candidate_params,
            seed=seed,
            jobs=jobs,
        )
        stage_records: list[dict[str, object]] = [
            {
                "family": "base",
                "previous_retained": [],
                "proposed": ["base"],
                "next_retained": ["base"],
                "affected_heads": ["count", "candidate"],
                "reused_heads": [],
                "count_selection": _grid_payload(count_grid),
                "candidate_selection": _grid_payload(candidate_grid),
                "count_feature_spec": _spec_payload(
                    head_feature_spec("count", retained)
                ),
                "candidate_feature_spec": _spec_payload(
                    head_feature_spec("candidate", retained)
                ),
                "oof_sha256": retained_oof.sha256,
                "oof_chain_count": len(retained_oof.rows),
                "oof_chain_id_sha256": _id_hash(
                    [row.chain_id for row in retained_oof.rows]
                ),
                "gates": {"status": "mandatory_base_no_bootstrap"},
                "diagnostics": _all_cohorts(retained_oof),
                "retained": True,
                "reason": "mandatory_base",
            }
        ]
        grid_history: list[dict[str, object]] = [
            {
                "family": "base",
                "count": _grid_payload(count_grid),
                "candidate": _grid_payload(candidate_grid),
            }
        ]
        completed_family_index = 0
        if checkpoint_store is not None:
            checkpoint_store.write_stage(
                completed_family_index=completed_family_index,
                retained_families=retained,
                count_params=retained_count_params,
                candidate_params=retained_candidate_params,
                retained_oof=retained_oof,
                stage_records=stage_records,
                grid_history=grid_history,
                accepted_ids=accepted_ids,
            )
    else:
        retained = restored.retained_families
        retained_count_params = restored.count_params
        retained_candidate_params = restored.candidate_params
        retained_oof = restored.retained_oof
        stage_records = [dict(record) for record in restored.stage_records]
        grid_history = [dict(record) for record in restored.grid_history]
        completed_family_index = restored.completed_family_index

    for family_index in range(
        completed_family_index + 1, len(FEATURE_FAMILY_ORDER)
    ):
        family = FEATURE_FAMILY_ORDER[family_index]
        previous = retained
        proposal = (*retained, family)
        if family == "global_count":
            proposed_count_grid = select_head_hyperparameters(
                data.corpus,
                data.assignments,
                "count",
                proposal,
                seed=seed,
                jobs=jobs,
                checkpoint_store=checkpoint_store,
                stage_family=family if checkpoint_store is not None else None,
            )
            proposed_count_params = proposed_count_grid.selected
            reused_count = None
            affected_heads = ["count", "candidate"]
            reused_heads: list[str] = []
        else:
            proposed_count_grid = None
            proposed_count_params = retained_count_params
            reused_count = retained_oof.count_decisions
            affected_heads = ["candidate"]
            reused_heads = ["count"]
        proposed_candidate_grid = select_head_hyperparameters(
            data.corpus,
            data.assignments,
            "candidate",
            proposal,
            seed=seed,
            jobs=jobs,
            checkpoint_store=checkpoint_store,
            stage_family=family if checkpoint_store is not None else None,
        )
        proposed_candidate_params = proposed_candidate_grid.selected
        proposed_oof = generate_oof(
            data,
            proposal,
            proposed_count_params,
            proposed_candidate_params,
            seed=seed,
            jobs=jobs,
            reused_count_decisions=reused_count,
        )
        previous_by_id = {row.chain_id: row for row in retained_oof.rows}
        proposed_by_id = {row.chain_id: row for row in proposed_oof.rows}
        if set(previous_by_id) != set(proposed_by_id) or set(previous_by_id) != set(
            accepted_ids
        ):
            raise ValueError("ablation OOF common chain IDs are incomplete")
        bootstrap = paired_chain_bootstrap(
            {chain_id: previous_by_id[chain_id].ndo for chain_id in accepted_ids},
            {chain_id: proposed_by_id[chain_id].ndo for chain_id in accepted_ids},
            n_resamples=10_000,
            seed=37,
        )
        count_delta = _finite(
            float(
                np.mean(
                    [
                        proposed_by_id[chain_id].count_correct
                        - previous_by_id[chain_id].count_correct
                        for chain_id in accepted_ids
                    ]
                )
            ),
            "ablation count accuracy delta",
        )
        contiguous = _cohort_delta(retained_oof, proposed_oof, "contiguous")
        discontinuous = _cohort_delta(retained_oof, proposed_oof, "discontinuous")
        if family == "discontinuity" and contiguous["status"] != "available":
            raise ValueError("discontinuity gate requires a nonempty contiguous cohort")
        contiguous_value = (
            float(contiguous["mean_ndo_delta"])
            if contiguous["status"] == "available"
            else 0.0
        )
        gate_report = AblationGateReport(
            family,
            bootstrap.mean_delta,
            bootstrap.ci_low,
            count_delta,
            contiguous_value,
        )
        retained_stage = retain_ablation(gate_report)
        next_retained = proposal if retained_stage else previous
        stage_record = {
            "family": family,
            "previous_retained": list(previous),
            "proposed": list(proposal),
            "next_retained": list(next_retained),
            "affected_heads": affected_heads,
            "reused_heads": reused_heads,
            "count_selection": (
                _grid_payload(proposed_count_grid)
                if proposed_count_grid is not None
                else {
                    "status": "reused_last_retained",
                    "selected": _params_payload(proposed_count_params),
                }
            ),
            "candidate_selection": _grid_payload(proposed_candidate_grid),
            "count_feature_spec": _spec_payload(
                head_feature_spec("count", proposal)
            ),
            "candidate_feature_spec": _spec_payload(
                head_feature_spec("candidate", proposal)
            ),
            "oof_sha256": proposed_oof.sha256,
            "oof_chain_count": len(proposed_oof.rows),
            "oof_chain_id_sha256": _id_hash(accepted_ids),
            "gates": {
                "overall_ndo": {
                    "method": "paired_chain_bootstrap_mean_delta",
                    "n": bootstrap.n,
                    "mean_delta": bootstrap.mean_delta,
                    "ci_low": bootstrap.ci_low,
                    "ci_high": bootstrap.ci_high,
                    "mean_threshold": 0.0,
                    "mean_strictly_greater": True,
                    "ci_low_threshold": 0.0,
                    "passes": bootstrap.mean_delta > 0.0
                    and bootstrap.ci_low >= 0.0,
                },
                "count_accuracy": {
                    "declared_gate": family == "global_count",
                    "delta": count_delta,
                    "threshold": 0.0,
                    "passes": count_delta >= 0.0,
                },
                "contiguous_guard": {
                    "declared_gate": family == "discontinuity",
                    "threshold": -0.005,
                    "diagnostic": contiguous,
                    "passes": contiguous_value >= -0.005,
                },
            },
            "diagnostics": {
                "contiguous": contiguous,
                "discontinuous": discontinuous,
                "proposed_cohorts": _all_cohorts(proposed_oof),
            },
            "retained": retained_stage,
            "reason": _reason(gate_report, retained_stage),
        }
        stage_records.append(stage_record)
        grid_history.append(
            {
                "family": family,
                "count": stage_record["count_selection"],
                "candidate": stage_record["candidate_selection"],
            }
        )
        if retained_stage:
            retained = proposal
            retained_count_params = proposed_count_params
            retained_candidate_params = proposed_candidate_params
            retained_oof = proposed_oof
        if checkpoint_store is not None:
            checkpoint_store.write_stage(
                completed_family_index=family_index,
                retained_families=retained,
                count_params=retained_count_params,
                candidate_params=retained_candidate_params,
                retained_oof=retained_oof,
                stage_records=stage_records,
                grid_history=grid_history,
                accepted_ids=accepted_ids,
            )

    final_oof = generate_oof(
        data,
        retained,
        retained_count_params,
        retained_candidate_params,
        seed=seed,
        jobs=jobs,
    )
    if final_oof.csv_bytes != retained_oof.csv_bytes:
        raise ValueError("fresh final OOF does not match last retained stage")
    count_spec = head_feature_spec("count", retained)
    candidate_spec = head_feature_spec("candidate", retained)
    count_batch = build_count_pairs(
        data.corpus.chains,
        data.corpus.counts,
        shared_features=count_spec.shared_features,
        item_features=count_spec.item_features,
        seed=seed,
    )
    candidate_batch = build_candidate_pairs(
        data.corpus.chains,
        data.corpus.candidates,
        shared_features=candidate_spec.shared_features,
        item_features=candidate_spec.item_features,
        seed=seed,
    )
    count_model = fit_head(count_batch, retained_count_params, seed)
    candidate_model = fit_head(candidate_batch, retained_candidate_params, seed)

    final_payload = {
        "retained_families": list(retained),
        "selected_hyperparameters": {
            "count": _params_payload(retained_count_params),
            "candidate": _params_payload(retained_candidate_params),
        },
        "feature_specs": {
            "count": _spec_payload(count_spec),
            "candidate": _spec_payload(candidate_spec),
        },
        "oof_filename": "oof_predictions.csv",
        "oof_sha256": final_oof.sha256,
        "oof_chain_count": len(final_oof.rows),
        "oof_cohorts": _all_cohorts(final_oof),
        "final_fit": {
            "count": {
                "rows": len(count_batch.y),
                "features": count_batch.x.shape[1],
                "classes": count_model.classes_.tolist(),
            },
            "candidate": {
                "rows": len(candidate_batch.y),
                "features": candidate_batch.x.shape[1],
                "classes": candidate_model.classes_.tolist(),
            },
        },
    }
    cv_report = {
        **common,
        "grid_history": grid_history,
        "final": final_payload,
    }
    ablation_report = {
        **common,
        "stages": stage_records,
        "final": final_payload,
    }
    cv_bytes = _canonical_json_bytes(cv_report)
    ablation_bytes = _canonical_json_bytes(ablation_report)
    # All payloads are fully rendered and validated before the first install.
    json.loads(cv_bytes)
    json.loads(ablation_bytes)
    out_dir = Path(out_dir)
    outputs = {
        "oof_predictions.csv": final_oof.csv_bytes,
        "cv_report.json": cv_bytes,
        "ablation_report.json": ablation_bytes,
    }
    for name, payload in outputs.items():
        _atomic_write(out_dir / name, payload)
    hashes = {name: _sha256_bytes(payload) for name, payload in outputs.items()}
    return TrainingRun(
        count_model,
        candidate_model,
        count_batch,
        candidate_batch,
        retained,
        retained_count_params,
        retained_candidate_params,
        final_oof,
        cv_report,
        ablation_report,
        hashes,
    )
