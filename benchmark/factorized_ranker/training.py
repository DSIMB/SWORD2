"""Grouped development CV and ordered ablations for factorized ranker heads."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import sys
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


def _evaluate_head_configuration(
    corpus: CorpusTables,
    assignments: Sequence[FoldAssignment],
    head: Literal["count", "candidate"],
    params: Hyperparameters,
    retained_families: Sequence[str],
    *,
    seed: int = 37,
) -> HeadGridEvaluation:
    families = validate_retained_families(retained_families)
    canonical_assignments = _validate_cv_inputs(corpus, assignments)
    spec = head_feature_spec(head, families)
    all_ids = {assignment.chain_id for assignment in canonical_assignments}
    primary_values: list[float] = []
    secondary_values: list[float] = []
    fold_results: list[HeadFoldResult] = []
    for fold in range(5):
        validation_ids = {
            assignment.chain_id
            for assignment in canonical_assignments
            if assignment.fold == fold
        }
        training_ids = all_ids - validation_ids
        training_chains = _subset(corpus.chains, training_ids)
        if head == "count":
            batch = build_count_pairs(
                training_chains,
                _subset(corpus.counts, training_ids),
                shared_features=spec.shared_features,
                item_features=spec.item_features,
                seed=seed,
            )
        else:
            batch = build_candidate_pairs(
                training_chains,
                _subset(corpus.candidates, training_ids),
                shared_features=spec.shared_features,
                item_features=spec.item_features,
                seed=seed,
            )
        model = fit_head(batch, params, seed)
        fold_primary: list[float] = []
        fold_secondary: list[float] = []
        for chain_id in sorted(validation_ids):
            chain_rows = _records(_subset(corpus.chains, {chain_id}))
            candidate_rows = _records(_subset(corpus.candidates, {chain_id}))
            if len(chain_rows) != 1:
                raise ValueError("validation chain row is missing or duplicated")
            if head == "count":
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
        if not fold_primary:
            raise ValueError("head validation fold is empty")
        fold_primary_mean = _finite(float(np.mean(fold_primary)), "fold primary objective")
        primary_values.extend(fold_primary)
        if head == "count":
            fold_secondary_mean: float | None = _finite(
                float(np.mean(fold_secondary)), "fold secondary objective"
            )
            secondary_values.extend(fold_secondary)
        else:
            fold_secondary_mean = None
        fold_results.append(
            HeadFoldResult(
                fold,
                tuple(sorted(validation_ids)),
                fold_primary_mean,
                fold_secondary_mean,
            )
        )
    primary = _finite(float(np.mean(primary_values)), "head primary objective")
    secondary = (
        _finite(float(np.mean(secondary_values)), "head secondary objective")
        if head == "count"
        else None
    )
    return HeadGridEvaluation(head, params, primary, secondary, tuple(fold_results))


def select_head_hyperparameters(
    corpus: CorpusTables,
    assignments: Sequence[FoldAssignment],
    head: Literal["count", "candidate"],
    retained_families: Sequence[str],
    *,
    seed: int = 37,
) -> GridSelection:
    if head not in {"count", "candidate"}:
        raise ValueError("head must be count or candidate")
    evaluations = tuple(
        _evaluate_head_configuration(
            corpus,
            assignments,
            head,
            params,
            retained_families,
            seed=seed,
        )
        for params in MODEL_GRID
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


def generate_oof(
    data: VerifiedTrainingData,
    retained_families: Sequence[str],
    count_params: Hyperparameters,
    candidate_params: Hyperparameters,
    *,
    seed: int = 37,
    reused_count_decisions: Mapping[str, CountDecision] | None = None,
) -> OOFResult:
    families = validate_retained_families(retained_families)
    count_spec = head_feature_spec("count", families)
    candidate_spec = head_feature_spec("candidate", families)
    assignments = data.assignments
    all_ids = {assignment.chain_id for assignment in assignments}
    metadata = {entry.entry_id: entry for entry in data.metadata}
    assignment_by_id = {assignment.chain_id: assignment for assignment in assignments}
    if set(metadata) != all_ids or set(assignment_by_id) != all_ids:
        raise ValueError("OOF metadata/fold populations disagree")
    if reused_count_decisions is not None and set(reused_count_decisions) != all_ids:
        raise ValueError("reused count decisions do not cover exact accepted chains")

    rows: list[OOFRow] = []
    count_decisions: dict[str, CountDecision] = {}
    for fold in range(5):
        validation_ids = {
            assignment.chain_id for assignment in assignments if assignment.fold == fold
        }
        training_ids = all_ids - validation_ids
        training_chains = _subset(data.corpus.chains, training_ids)
        if reused_count_decisions is None:
            count_batch = build_count_pairs(
                training_chains,
                _subset(data.corpus.counts, training_ids),
                shared_features=count_spec.shared_features,
                item_features=count_spec.item_features,
                seed=seed,
            )
            count_model: GradientBoostingClassifier | None = fit_head(
                count_batch, count_params, seed
            )
        else:
            count_model = None
        candidate_batch = build_candidate_pairs(
            training_chains,
            _subset(data.corpus.candidates, training_ids),
            shared_features=candidate_spec.shared_features,
            item_features=candidate_spec.item_features,
            seed=seed,
        )
        candidate_model = fit_head(candidate_batch, candidate_params, seed)
        label_cohorts = chain_label_cohorts(assignments, fold)

        for chain_id in sorted(validation_ids):
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
                    fold=fold,
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
    ordered = tuple(sorted(rows, key=lambda row: row.chain_id))
    if {row.chain_id for row in ordered} != all_ids or len(ordered) != len(all_ids):
        raise ValueError("OOF rows do not cover exact accepted chains")
    if set(count_decisions) != all_ids:
        raise ValueError("OOF count decisions do not cover exact accepted chains")
    csv_bytes = _render_oof(ordered)
    return OOFResult(ordered, count_decisions, csv_bytes, _sha256_bytes(csv_bytes))


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
) -> TrainingRun:
    if seed != 37 or isinstance(seed, bool):
        raise ValueError("grouped training seed is frozen at 37")
    accepted_ids = tuple(sorted(assignment.chain_id for assignment in data.assignments))
    common = _common_report(data, normalized_command, accepted_ids)
    retained = ("base",)
    count_grid = select_head_hyperparameters(
        data.corpus, data.assignments, "count", retained, seed=seed
    )
    candidate_grid = select_head_hyperparameters(
        data.corpus, data.assignments, "candidate", retained, seed=seed
    )
    retained_count_params = count_grid.selected
    retained_candidate_params = candidate_grid.selected
    retained_oof = generate_oof(
        data,
        retained,
        retained_count_params,
        retained_candidate_params,
        seed=seed,
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
            "count_feature_spec": _spec_payload(head_feature_spec("count", retained)),
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

    for family in FEATURE_FAMILY_ORDER[1:]:
        previous = retained
        proposal = (*retained, family)
        if family == "global_count":
            proposed_count_grid = select_head_hyperparameters(
                data.corpus, data.assignments, "count", proposal, seed=seed
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
            data.corpus, data.assignments, "candidate", proposal, seed=seed
        )
        proposed_candidate_params = proposed_candidate_grid.selected
        proposed_oof = generate_oof(
            data,
            proposal,
            proposed_count_params,
            proposed_candidate_params,
            seed=seed,
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

    final_oof = generate_oof(
        data,
        retained,
        retained_count_params,
        retained_candidate_params,
        seed=seed,
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
