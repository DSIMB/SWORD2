from __future__ import annotations

import hashlib
import json
import multiprocessing
import shutil
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

import benchmark.factorized_ranker.training as training
from benchmark.factorized_ranker.corpus import CANDIDATE_FIELDS, CHAIN_FIELDS, COUNT_FIELDS
from benchmark.factorized_ranker.folds import FoldAssignment
from benchmark.factorized_ranker.pairs import PairBatch
from benchmark.factorized_ranker.schema import (
    BASE_CANDIDATE_FEATURES,
    BOUNDARY_LOCAL_FEATURES,
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    DISCONTINUITY_FEATURES,
    DOMAIN_CONDITIONED_FEATURES,
    GLOBAL_FEATURES,
    RELATIVE_HIERARCHY_FEATURES,
    pair_feature_names,
)
from benchmark.factorized_ranker.training import (
    BASE_COUNT_ITEM_FEATURES,
    BASE_COUNT_SUMMARY_SOURCES,
    FEATURE_FAMILY_ORDER,
    MODEL_GRID,
    AblationGateReport,
    CorpusTables,
    Hyperparameters,
    OOFResult,
    OOFRow,
    TrainingCheckpointStore,
    build_training_checkpoint_context,
    cross_validate_configuration,
    fit_head,
    generate_oof,
    head_feature_spec,
    load_verified_training_data,
    retain_ablation,
    select_head_hyperparameters,
    validate_retained_families,
)


FIXTURES = Path(__file__).parent / "fixtures"


def _checkpoint_bytes(root: Path) -> dict[str, bytes]:
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in sorted(root.rglob("*.json"))
    }


def _always_fail_fit(*_args: object, **_kwargs: object) -> object:
    raise RuntimeError("injected worker failure")


@pytest.mark.parametrize("value", [True, False, 0, -1, 9, 1.0, "8", None])
def test_parallel_jobs_reject_invalid_library_values(value: object) -> None:
    with pytest.raises(ValueError, match="jobs"):
        training.validate_training_jobs(value)  # type: ignore[arg-type]


def test_parallel_jobs_accept_exact_supported_range() -> None:
    assert [training.validate_training_jobs(value) for value in range(1, 9)] == list(
        range(1, 9)
    )


def test_model_grid_is_exact_and_bounded() -> None:
    assert len(MODEL_GRID) == 8
    assert {
        (p.n_estimators, p.learning_rate, p.min_samples_leaf, p.max_depth)
        for p in MODEL_GRID
    } == {
        (n, lr, leaf, 3)
        for n in (64, 96)
        for lr in (0.03, 0.05)
        for leaf in (32, 64)
    }


def test_feature_family_mapping_is_exact_and_cumulative() -> None:
    assert FEATURE_FAMILY_ORDER == (
        "base",
        "global_count",
        "domain_conditioned",
        "boundary_local",
        "relative_hierarchy",
        "discontinuity",
    )
    assert BASE_COUNT_SUMMARY_SOURCES == (
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
    assert BASE_COUNT_ITEM_FEATURES[:4] == (
        "count_num_domains",
        "count_n_candidates",
        "count_candidate_fraction",
        "count_modal_distance",
    )
    assert not any("legacy_distance" in name for name in BASE_COUNT_ITEM_FEATURES)

    base_count = head_feature_spec("count", ("base",))
    base_candidate = head_feature_spec("candidate", ("base",))
    assert base_count.shared_features == ()
    assert base_count.item_features == BASE_COUNT_ITEM_FEATURES
    assert base_candidate.shared_features == ()
    assert base_candidate.item_features == BASE_CANDIDATE_FEATURES
    assert base_candidate.pair_features == pair_feature_names((), BASE_CANDIDATE_FEATURES)

    global_count = head_feature_spec("count", ("base", "global_count"))
    global_candidate = head_feature_spec("candidate", ("base", "global_count"))
    assert global_count.shared_features == GLOBAL_FEATURES
    assert global_count.item_features == COUNT_ITEM_FEATURES
    assert global_candidate.shared_features == GLOBAL_FEATURES
    assert global_candidate.item_features == BASE_CANDIDATE_FEATURES

    all_candidate = head_feature_spec("candidate", FEATURE_FAMILY_ORDER)
    assert all_candidate.item_features == (
        *BASE_CANDIDATE_FEATURES,
        *DOMAIN_CONDITIONED_FEATURES,
        *BOUNDARY_LOCAL_FEATURES,
        *RELATIVE_HIERARCHY_FEATURES,
        *DISCONTINUITY_FEATURES,
    ) == CANDIDATE_FEATURES


@pytest.mark.parametrize(
    "families",
    [
        (),
        ("global_count",),
        ("base", "base"),
        ("base", "unknown"),
        ("base", "boundary_local", "domain_conditioned"),
    ],
)
def test_invalid_retained_family_sequences_fail(families) -> None:
    with pytest.raises(ValueError):
        validate_retained_families(families)


def test_rejected_families_may_be_skipped_but_order_is_preserved() -> None:
    assert validate_retained_families(
        ("base", "boundary_local", "discontinuity")
    ) == ("base", "boundary_local", "discontinuity")


def _batch() -> PairBatch:
    arrays = PairBatch(
        feature_names=("f0", "f1"),
        x=np.asarray([[0.0, 1.0], [0.0, -1.0], [1.0, 2.0], [1.0, -2.0]], dtype=np.float64),
        y=np.asarray([1, 0, 1, 0], dtype=np.int8),
        sample_weight=np.asarray([0.25] * 4, dtype=np.float64),
        chain_ids=np.asarray(["a", "a", "b", "b"], dtype="<U1"),
        left_ids=np.asarray(["a" * 64, "b" * 64, "c" * 64, "d" * 64], dtype="<U64"),
        right_ids=np.asarray(["b" * 64, "a" * 64, "d" * 64, "c" * 64], dtype="<U64"),
    )
    return arrays


def test_fit_head_uses_exact_approved_model_contract() -> None:
    params = Hyperparameters(64, 0.03, 32)
    model = fit_head(_batch(), params, seed=37)
    observed = model.get_params()
    assert observed["n_estimators"] == 64
    assert observed["learning_rate"] == 0.03
    assert observed["min_samples_leaf"] == 32
    assert observed["max_depth"] == 3
    assert observed["loss"] == "log_loss"
    assert observed["criterion"] == "friedman_mse"
    assert observed["random_state"] == 37
    assert model.classes_.tolist() == [0, 1]
    assert model.n_features_in_ == 2


@pytest.mark.parametrize("defect", ["empty", "one_class", "shape", "dtype", "nan", "weight", "feature"])
def test_fit_head_rejects_malformed_batches(defect: str) -> None:
    batch = _batch()
    if defect == "empty":
        batch = replace(
            batch,
            x=np.empty((0, 2), dtype=np.float64),
            y=np.empty(0, dtype=np.int8),
            sample_weight=np.empty(0, dtype=np.float64),
            chain_ids=np.empty(0, dtype="<U1"),
            left_ids=np.empty(0, dtype="<U64"),
            right_ids=np.empty(0, dtype="<U64"),
        )
    elif defect == "one_class":
        batch = replace(batch, y=np.ones(4, dtype=np.int8))
    elif defect == "shape":
        batch = replace(batch, y=np.asarray([0, 1], dtype=np.int8))
    elif defect == "dtype":
        batch = replace(batch, x=batch.x.astype(np.float32))
    elif defect == "nan":
        x = batch.x.copy()
        x[0, 0] = np.nan
        batch = replace(batch, x=x)
    elif defect == "weight":
        batch = replace(batch, sample_weight=np.zeros(4, dtype=np.float64))
    else:
        batch = replace(batch, feature_names=("only_one",))
    with pytest.raises(ValueError):
        fit_head(batch, MODEL_GRID[0])


@pytest.mark.parametrize(
    "params,seed",
    [
        (Hyperparameters(1, 0.03, 32), 37),
        (MODEL_GRID[0], True),
        (MODEL_GRID[0], -1),
        (MODEL_GRID[0], 2**32),
    ],
)
def test_fit_head_rejects_unapproved_params_or_seed(params, seed) -> None:
    with pytest.raises(ValueError):
        fit_head(_batch(), params, seed=seed)


def test_ablation_guards_are_exact() -> None:
    assert retain_ablation(AblationGateReport("boundary_local", 0.01, 0.0, 0.0, 0.0))
    assert not retain_ablation(
        AblationGateReport("boundary_local", 0.01, -1e-6, 0.0, 0.0)
    )
    assert not retain_ablation(
        AblationGateReport("global_count", 0.01, 0.0, -1e-6, 0.0)
    )
    assert not retain_ablation(
        AblationGateReport("discontinuity", 0.01, 0.0, 0.0, -0.0051)
    )
    assert retain_ablation(
        AblationGateReport("discontinuity", 0.01, 0.0, 0.0, -0.005)
    )


def _synthetic_corpus() -> tuple[CorpusTables, tuple[FoldAssignment, ...]]:
    chain_rows = []
    count_rows = []
    candidate_rows = []
    assignments = []
    for index in range(10):
        chain_id = f"chain{index}"
        true_count = index % 5 + 1
        chain_rows.append(
            {
                "chain_id": chain_id,
                "n_true_domains": true_count,
                **{name: index + offset / 100.0 for offset, name in enumerate(GLOBAL_FEATURES)},
            }
        )
        counts = sorted({true_count, max(1, true_count - 1), true_count + 1})
        source_index = 0
        for count in counts:
            count_row = {
                "chain_id": chain_id,
                **{
                    name: count + offset / 100.0
                    for offset, name in enumerate(COUNT_ITEM_FEATURES)
                },
            }
            count_row["count_num_domains"] = count
            count_rows.append(count_row)
            for candidate_index, ndo in enumerate((0.2, 0.8)):
                canonical = f"{count}:{candidate_index}:{chain_id}"
                candidate = {
                    "candidate_id": hashlib.sha256(
                        (chain_id + "\0" + canonical).encode()
                    ).hexdigest(),
                    "chain_id": chain_id,
                    "canonical_delineation": canonical,
                    "source_index": source_index,
                    "legacy_distance": float(count),
                    **{
                        name: count + candidate_index + offset / 1000.0
                        for offset, name in enumerate(CANDIDATE_FEATURES)
                    },
                    "n_true_domains": true_count,
                    "n_pred_domains": count,
                    "ndo": ndo,
                    "iou": ndo,
                    "boundary_f1_10": ndo,
                    "matched_dice": ndo,
                    "d_count_acc": float(count == true_count),
                    "S": ndo,
                    "is_oracle_s": int(candidate_index == 1),
                }
                candidate["num_domains"] = count
                candidate_rows.append(candidate)
                source_index += 1
        component_id = hashlib.sha256(chain_id.encode()).hexdigest()
        assignments.append(
            FoldAssignment(
                chain_id=chain_id,
                pdb_id=f"{index:04x}",
                family_combination=(f"family{index}",),
                family_labels=(f"family{index}",),
                true_count_bin="5+" if true_count >= 5 else str(true_count),
                length_bin="<250",
                component_id=component_id,
                fold=index % 5,
            )
        )
    return (
        CorpusTables(
            pd.DataFrame(chain_rows, columns=CHAIN_FIELDS),
            pd.DataFrame(count_rows, columns=COUNT_FIELDS),
            pd.DataFrame(candidate_rows, columns=CANDIDATE_FIELDS),
        ),
        tuple(assignments),
    )


def test_training_never_sees_validation_chain(monkeypatch: pytest.MonkeyPatch) -> None:
    corpus, folds = _synthetic_corpus()
    original_fit = training.fit_head
    original_count_pairs = training.build_count_pairs
    original_candidate_pairs = training.build_candidate_pairs
    fitted: list[set[str]] = []
    built: list[tuple[str, set[str]]] = []

    def capture(batch, params, seed=37):
        fitted.append(set(batch.chain_ids.tolist()))
        return original_fit(batch, params, seed)

    def capture_count_pairs(chains, counts, **kwargs):
        built.append(("count", set(chains["chain_id"].tolist())))
        return original_count_pairs(chains, counts, **kwargs)

    def capture_candidate_pairs(chains, candidates, **kwargs):
        built.append(("candidate", set(chains["chain_id"].tolist())))
        return original_candidate_pairs(chains, candidates, **kwargs)

    monkeypatch.setattr(training, "fit_head", capture)
    monkeypatch.setattr(training, "build_count_pairs", capture_count_pairs)
    monkeypatch.setattr(training, "build_candidate_pairs", capture_candidate_pairs)
    result = cross_validate_configuration(
        corpus,
        folds,
        MODEL_GRID[0],
        MODEL_GRID[0],
        ("base",),
    )
    assert len(built) == len(fitted) == 10
    all_ids = {assignment.chain_id for assignment in folds}
    for fold in range(5):
        validation = {assignment.chain_id for assignment in folds if assignment.fold == fold}
        expected_training = all_ids - validation
        assert built[2 * fold] == ("count", expected_training)
        assert built[2 * fold + 1] == ("candidate", expected_training)
        assert fitted[2 * fold] == expected_training
        assert fitted[2 * fold + 1] == expected_training
    assert result.count_accuracy >= 0.0
    assert result.count_abs_error >= 0.0
    assert 0.2 <= result.candidate_ndo <= 0.8


def test_cross_validation_is_source_order_independent() -> None:
    corpus, folds = _synthetic_corpus()
    first = cross_validate_configuration(
        corpus, folds, MODEL_GRID[0], MODEL_GRID[0], ("base",)
    )
    reversed_corpus = CorpusTables(
        corpus.chains.iloc[::-1].reset_index(drop=True),
        corpus.counts.iloc[::-1].reset_index(drop=True),
        corpus.candidates.iloc[::-1].reset_index(drop=True),
    )
    second = cross_validate_configuration(
        reversed_corpus,
        tuple(reversed(folds)),
        MODEL_GRID[0],
        MODEL_GRID[0],
        ("base",),
    )
    assert first == second


def test_grid_ties_use_exact_approved_parameter_order(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[tuple[str, Hyperparameters]] = []

    def evaluate(_corpus, _assignments, head, params, _families, seed=37):
        calls.append((head, params))
        return training.HeadGridEvaluation(
            head=head,
            params=params,
            primary=0.5,
            secondary=1.0 if head == "count" else None,
            folds=(),
        )

    monkeypatch.setattr(training, "_evaluate_head_configuration", evaluate)
    corpus, assignments = _synthetic_corpus()
    count = training.select_head_hyperparameters(
        corpus, assignments, "count", ("base",)
    )
    candidate = training.select_head_hyperparameters(
        corpus, assignments, "candidate", ("base",)
    )
    expected = Hyperparameters(64, 0.03, 64)
    assert count.selected == candidate.selected == expected
    assert [head for head, _params in calls].count("count") == 8
    assert [head for head, _params in calls].count("candidate") == 8


def test_grid_fold_checkpoints_are_exact_and_reused(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    corpus, assignments = _synthetic_corpus()
    context = {
        "source_git_commit": "a" * 40,
        "training_source_sha256": "b" * 64,
        "corpus_manifest_sha256": "c" * 64,
        "fold_manifest_sha256": "d" * 64,
        "feature_schema_sha256": "e" * 64,
        "seed": 37,
    }
    store = TrainingCheckpointStore(tmp_path, context)
    original = training.fit_head
    fits = 0

    def capture(batch, params, seed=37):
        nonlocal fits
        fits += 1
        return original(batch, params, seed)

    monkeypatch.setattr(training, "fit_head", capture)
    first = select_head_hyperparameters(
        corpus,
        assignments,
        "count",
        ("base",),
        checkpoint_store=store,
        stage_family="base",
    )
    assert fits == 40
    payloads = sorted((tmp_path / "grid").glob("*.json"))
    assert len(payloads) == 40
    payload = json.loads(payloads[0].read_bytes())
    assert set(payload) == {
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
    assert payload["context"] == context

    def forbidden_fit(*_args, **_kwargs):
        raise AssertionError("valid fold checkpoint was not reused")

    monkeypatch.setattr(training, "fit_head", forbidden_fit)
    second = select_head_hyperparameters(
        corpus,
        tuple(reversed(assignments)),
        "count",
        ("base",),
        checkpoint_store=store,
        stage_family="base",
    )
    assert second == first
    with pytest.raises(ValueError, match="context"):
        TrainingCheckpointStore(tmp_path, {**context, "seed": 38})


@pytest.mark.skipif(
    "fork" not in multiprocessing.get_all_start_methods(), reason="fork required"
)
def test_parallel_grid_matches_sequential_selection_and_checkpoint_bytes(
    tmp_path: Path,
) -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    context = {"source_git_commit": "a" * 40, "seed": 37}
    sequential_root = tmp_path / "sequential"
    parallel_root = tmp_path / "parallel"
    sequential = select_head_hyperparameters(
        data.corpus,
        data.assignments,
        "count",
        ("base",),
        jobs=1,
        checkpoint_store=TrainingCheckpointStore(sequential_root, context),
        stage_family="base",
    )
    parallel = select_head_hyperparameters(
        data.corpus,
        data.assignments,
        "count",
        ("base",),
        jobs=8,
        checkpoint_store=TrainingCheckpointStore(parallel_root, context),
        stage_family="base",
    )
    assert parallel == sequential
    assert _checkpoint_bytes(parallel_root) == _checkpoint_bytes(sequential_root)


def test_grid_aggregation_ignores_worker_completion_order() -> None:
    canonical: dict[training._HeadFoldTask, training._HeadFoldWorkResult] = {}
    for params_index, params in enumerate(MODEL_GRID):
        for fold in range(5):
            task = training._HeadFoldTask("count", params, ("base",), fold, 37)
            primary = (float((params_index + fold) % 2),)
            secondary = (float(fold),)
            canonical[task] = training._HeadFoldWorkResult(
                task,
                training.HeadFoldResult(
                    fold,
                    (f"chain{fold}",),
                    float(np.mean(primary)),
                    float(np.mean(secondary)),
                ),
                primary,
                secondary,
            )
    reversed_results = dict(reversed(tuple(canonical.items())))
    assert training._assemble_grid_evaluations(
        "count", reversed_results
    ) == training._assemble_grid_evaluations("count", canonical)


@pytest.mark.skipif(
    "fork" not in multiprocessing.get_all_start_methods(), reason="fork required"
)
def test_parallel_grid_reuses_partial_checkpoints_without_rewriting_them(
    tmp_path: Path,
) -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    context = {"source_git_commit": "a" * 40, "seed": 37}
    complete_root = tmp_path / "complete"
    complete = select_head_hyperparameters(
        data.corpus,
        data.assignments,
        "count",
        ("base",),
        jobs=1,
        checkpoint_store=TrainingCheckpointStore(complete_root, context),
        stage_family="base",
    )
    partial_root = tmp_path / "partial"
    partial_grid = partial_root / "grid"
    partial_grid.mkdir(parents=True)
    shutil.copy2(
        complete_root / "checkpoint_manifest.json",
        partial_root / "checkpoint_manifest.json",
    )
    for source in sorted((complete_root / "grid").glob("*.json"))[:13]:
        shutil.copy2(source, partial_grid / source.name)
    cached_paths = tuple(sorted(partial_grid.glob("*.json")))
    cached_bytes = {path.name: path.read_bytes() for path in cached_paths}
    cached_mtimes = {path.name: path.stat().st_mtime_ns for path in cached_paths}

    resumed = select_head_hyperparameters(
        data.corpus,
        tuple(reversed(data.assignments)),
        "count",
        ("base",),
        jobs=8,
        checkpoint_store=TrainingCheckpointStore(partial_root, context),
        stage_family="base",
    )

    assert resumed == complete
    assert len(list(partial_grid.glob("*.json"))) == 40
    assert {path.name: path.read_bytes() for path in cached_paths} == cached_bytes
    assert {path.name: path.stat().st_mtime_ns for path in cached_paths} == cached_mtimes


@pytest.mark.skipif(
    "fork" not in multiprocessing.get_all_start_methods(), reason="fork required"
)
def test_parallel_grid_worker_failure_never_writes_stage_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    store = TrainingCheckpointStore(
        tmp_path, {"source_git_commit": "a" * 40, "seed": 37}
    )
    monkeypatch.setattr(training, "fit_head", _always_fail_fit)
    with pytest.raises(RuntimeError, match="injected worker failure"):
        select_head_hyperparameters(
            data.corpus,
            data.assignments,
            "count",
            ("base",),
            jobs=8,
            checkpoint_store=store,
            stage_family="base",
        )
    assert not (tmp_path / "stage_state.json").exists()
    assert not tuple((tmp_path / "grid").glob(".*.tmp"))


def test_completed_stage_checkpoint_round_trips_exact_oof_state(tmp_path: Path) -> None:
    context = {"source_git_commit": "a" * 40, "seed": 37}
    store = TrainingCheckpointStore(tmp_path, context)
    row = OOFRow(
        chain_id="chain0",
        fold=0,
        selected_count=2,
        selected_candidate_id="candidate0",
        canonical_delineation="1-10 11-20",
        count_borda_score=0.75,
        candidate_borda_score=0.625,
        n_true_domains=2,
        count_correct=1,
        ndo=0.8,
        boundary_f1_10=0.9,
        matched_dice=0.85,
        total_regret=0.1,
        count_regret=0.05,
        within_count_regret=0.05,
        true_count_bin="2",
        length_bin="<250",
        continuity_cohort="contiguous",
        label_cohort="unseen",
        count_tie_break="none",
        candidate_tie_break="none",
    )
    csv_bytes = training._render_oof((row,))
    oof = OOFResult(
        (row,),
        {"chain0": training.CountDecision(2, 0.75, "none")},
        csv_bytes,
        hashlib.sha256(csv_bytes).hexdigest(),
    )
    stage_records = [{"family": "base", "retained": True}]
    grid_history = [{"family": "base", "count": {}, "candidate": {}}]
    store.write_stage(
        completed_family_index=0,
        retained_families=("base",),
        count_params=MODEL_GRID[0],
        candidate_params=MODEL_GRID[1],
        retained_oof=oof,
        stage_records=stage_records,
        grid_history=grid_history,
        accepted_ids=("chain0",),
    )
    restored = store.load_stage(("chain0",))
    assert restored is not None
    assert restored.completed_family_index == 0
    assert restored.retained_families == ("base",)
    assert restored.count_params == MODEL_GRID[0]
    assert restored.candidate_params == MODEL_GRID[1]
    assert restored.retained_oof == oof
    assert list(restored.stage_records) == stage_records
    assert list(restored.grid_history) == grid_history
    state_path = tmp_path / "stage_state.json"
    payload = json.loads(state_path.read_bytes())
    assert set(payload) == {
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

    payload["retained_oof"]["sha256"] = "0" * 64
    state_path.write_bytes(training._canonical_json_bytes(payload))
    with pytest.raises(ValueError, match="OOF hash"):
        store.load_stage(("chain0",))


def test_committed_fixture_verifies_exact_hash_chain_and_source_cohorts() -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus",
        FIXTURES / "factorized_folds.json",
    )
    assert len(data.assignments) == 15
    assert {assignment.fold for assignment in data.assignments} == set(range(5))
    assert {assignment.true_count_bin for assignment in data.assignments} == {
        "2",
        "3",
        "4",
        "5+",
    }
    assert {assignment.length_bin for assignment in data.assignments} == {
        "<250",
        "250-349",
        "350-449",
        "450+",
    }
    assert not any(
        token in column.casefold()
        for column in data.corpus.candidates.columns
        for token in ("merizo", "chainsaw")
    )


def test_fixture_oof_is_complete_and_regrets_decompose() -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus",
        FIXTURES / "factorized_folds.json",
    )
    oof = generate_oof(data, ("base",), MODEL_GRID[0], MODEL_GRID[0])
    assert len(oof.rows) == len(data.assignments) == 15
    assert {row.chain_id for row in oof.rows} == {
        assignment.chain_id for assignment in data.assignments
    }
    assert {row.continuity_cohort for row in oof.rows} == {
        "contiguous",
        "discontinuous",
    }
    assert {row.label_cohort for row in oof.rows} == {"seen", "unseen"}
    unavailable_truth = next(row for row in oof.rows if row.chain_id == "1aroP")
    assert unavailable_truth.n_true_domains == 5
    assert unavailable_truth.count_correct == 0
    for row in oof.rows:
        assert row.total_regret >= 0.0
        assert row.count_regret >= 0.0
        assert row.within_count_regret >= 0.0
        assert row.total_regret == pytest.approx(
            row.count_regret + row.within_count_regret, abs=1e-12
        )
    assert hashlib.sha256(oof.csv_bytes).hexdigest() == oof.sha256


@pytest.mark.skipif(
    "fork" not in multiprocessing.get_all_start_methods(), reason="fork required"
)
def test_parallel_oof_matches_sequential_bytes_and_decisions() -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    sequential = generate_oof(
        data, ("base",), MODEL_GRID[0], MODEL_GRID[0], jobs=1
    )
    parallel = generate_oof(
        data, ("base",), MODEL_GRID[0], MODEL_GRID[0], jobs=8
    )
    assert parallel == sequential
    assert parallel.csv_bytes == sequential.csv_bytes
    assert parallel.sha256 == sequential.sha256
    assert parallel.count_decisions == sequential.count_decisions


@pytest.mark.skipif(
    "fork" not in multiprocessing.get_all_start_methods(), reason="fork required"
)
def test_parallel_oof_is_repeatable_with_reused_count_decisions() -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    reference = generate_oof(
        data, ("base",), MODEL_GRID[0], MODEL_GRID[0], jobs=1
    )
    sequential = generate_oof(
        data,
        ("base",),
        MODEL_GRID[0],
        MODEL_GRID[0],
        jobs=1,
        reused_count_decisions=reference.count_decisions,
    )
    first = generate_oof(
        data,
        ("base",),
        MODEL_GRID[0],
        MODEL_GRID[0],
        jobs=8,
        reused_count_decisions=reference.count_decisions,
    )
    second = generate_oof(
        data,
        ("base",),
        MODEL_GRID[0],
        MODEL_GRID[0],
        jobs=8,
        reused_count_decisions=reference.count_decisions,
    )
    assert first == second == sequential
    assert first.csv_bytes == second.csv_bytes == sequential.csv_bytes


def test_oof_combination_ignores_worker_completion_order() -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    tasks = tuple(
        training._OOFFoldTask(
            fold, ("base",), MODEL_GRID[0], MODEL_GRID[0], 37, False
        )
        for fold in range(5)
    )
    canonical = {
        task: training._generate_oof_fold(data, task, None) for task in tasks
    }
    reversed_results = dict(reversed(tuple(canonical.items())))
    assert training._combine_oof_fold_results(
        data, tasks, reversed_results
    ) == training._combine_oof_fold_results(data, tasks, canonical)


def test_oof_parent_rejects_changed_reused_count_decision() -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    reference = generate_oof(
        data, ("base",), MODEL_GRID[0], MODEL_GRID[0], jobs=1
    )
    task = training._OOFFoldTask(
        0, ("base",), MODEL_GRID[0], MODEL_GRID[0], 37, True
    )
    work = training._generate_oof_fold(data, task, reference.count_decisions)
    chain_id, decision = work.count_decisions[0]
    changed_decision = replace(
        decision,
        selected_count=decision.selected_count + 1,
        borda_score=decision.borda_score + 0.125,
    )
    changed_rows = tuple(
        replace(
            row,
            selected_count=changed_decision.selected_count,
            count_borda_score=changed_decision.borda_score,
        )
        if row.chain_id == chain_id
        else row
        for row in work.rows
    )
    changed_decisions = tuple(
        (observed_id, changed_decision if observed_id == chain_id else observed)
        for observed_id, observed in work.count_decisions
    )
    changed = replace(
        work, rows=changed_rows, count_decisions=changed_decisions
    )
    with pytest.raises(ValueError, match="reused count decision"):
        training._validate_oof_fold_work_result(
            data,
            task,
            changed,
            reused_count_decisions=reference.count_decisions,
        )


@pytest.mark.skipif(
    "fork" not in multiprocessing.get_all_start_methods(), reason="fork required"
)
def test_parallel_oof_worker_failure_returns_no_partial_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    monkeypatch.setattr(training, "fit_head", _always_fail_fit)
    with pytest.raises(RuntimeError, match="injected worker failure"):
        generate_oof(
            data, ("base",), MODEL_GRID[0], MODEL_GRID[0], jobs=8
        )


def test_grouped_training_forwards_jobs_to_final_oof(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus", FIXTURES / "factorized_folds.json"
    )
    retained_oof = generate_oof(
        data, ("base",), MODEL_GRID[0], MODEL_GRID[0], jobs=1
    )
    restored = training.TrainingStageState(
        completed_family_index=len(FEATURE_FAMILY_ORDER) - 1,
        retained_families=("base",),
        count_params=MODEL_GRID[0],
        candidate_params=MODEL_GRID[0],
        retained_oof=retained_oof,
        stage_records=(),
        grid_history=(),
    )

    class RestoredStore:
        def load_stage(self, _accepted_ids: tuple[str, ...]) -> object:
            return restored

    observed: list[int] = []

    def capture_oof(
        _data: object,
        _families: tuple[str, ...],
        _count_params: Hyperparameters,
        _candidate_params: Hyperparameters,
        *,
        seed: int,
        jobs: int,
        reused_count_decisions: object = None,
    ) -> OOFResult:
        assert seed == 37
        assert reused_count_decisions is None
        observed.append(jobs)
        return retained_oof

    monkeypatch.setattr(training, "generate_oof", capture_oof)
    training.run_grouped_training(
        data,
        tmp_path / "out",
        ["benchmark.train_factorized_ranker", "--jobs", "8"],
        jobs=8,
        checkpoint_store=RestoredStore(),  # type: ignore[arg-type]
    )
    assert observed == [8]


def test_training_argv_roles_and_model_outputs_are_all_or_none(
    capsys: pytest.CaptureFixture[str],
) -> None:
    from benchmark.train_factorized_ranker import main, normalize_training_argv

    assert normalize_training_argv(
        [
            "benchmark.train_factorized_ranker",
            "--corpus-dir=/tmp/a",
            "--fold-manifest",
            "/tmp/b",
            "--out-dir",
            "/tmp/c",
            "--seed",
            "37",
        ]
    ) == [
        "benchmark.train_factorized_ranker",
        "--corpus-dir=<CORPUS_DIR>",
        "--fold-manifest",
        "<FOLD_MANIFEST>",
        "--out-dir",
        "<OUT_DIR>",
        "--seed",
        "37",
    ]
    with pytest.raises(SystemExit) as error:
        main(
            [
                "--corpus-dir",
                "unused",
                "--fold-manifest",
                "unused",
                "--out-dir",
                "unused",
                "--count-model-out",
                "count.json",
            ]
        )
    assert error.value.code == 2
    assert "must be supplied together" in capsys.readouterr().err


def test_training_cli_exposes_resume(capsys: pytest.CaptureFixture[str]) -> None:
    from benchmark.train_factorized_ranker import main

    with pytest.raises(SystemExit) as error:
        main(["--help"])
    assert error.value.code == 0
    assert "--resume" in capsys.readouterr().out


def test_training_cli_exposes_bounded_jobs(
    capsys: pytest.CaptureFixture[str],
) -> None:
    from benchmark.train_factorized_ranker import main

    with pytest.raises(SystemExit) as help_exit:
        main(["--help"])
    assert help_exit.value.code == 0
    assert "--jobs" in capsys.readouterr().out
    with pytest.raises(SystemExit) as invalid_exit:
        main(
            [
                "--corpus-dir",
                "unused",
                "--fold-manifest",
                "unused",
                "--out-dir",
                "unused",
                "--jobs",
                "9",
            ]
        )
    assert invalid_exit.value.code == 2
    assert "1..=8" in capsys.readouterr().err


def test_training_cli_forwards_jobs_and_binds_normalized_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import benchmark.train_factorized_ranker as cli

    sentinel = object()
    observed: dict[str, object] = {}

    def fake_run(
        data: object,
        out_dir: Path,
        normalized_command: list[str],
        *,
        seed: int,
        jobs: int,
        checkpoint_store: object,
    ) -> SimpleNamespace:
        observed.update(
            data=data,
            out_dir=out_dir,
            normalized_command=normalized_command,
            seed=seed,
            jobs=jobs,
            checkpoint_store=checkpoint_store,
        )
        return SimpleNamespace(
            artifact_hashes={
                "oof_predictions.csv": "a" * 64,
                "cv_report.json": "b" * 64,
                "ablation_report.json": "c" * 64,
            }
        )

    monkeypatch.setattr(cli, "load_verified_training_data", lambda *_args: sentinel)
    monkeypatch.setattr(cli, "run_grouped_training", fake_run)
    assert (
        cli.main(
            [
                "--corpus-dir",
                str(tmp_path / "corpus"),
                "--fold-manifest",
                str(tmp_path / "folds.json"),
                "--out-dir",
                str(tmp_path / "out"),
                "--jobs",
                "8",
            ]
        )
        == 0
    )
    assert observed["data"] is sentinel
    assert observed["jobs"] == 8
    assert observed["seed"] == 37
    assert observed["checkpoint_store"] is None
    assert observed["normalized_command"] == [
        "benchmark.train_factorized_ranker",
        "--corpus-dir",
        "<CORPUS_DIR>",
        "--fold-manifest",
        "<FOLD_MANIFEST>",
        "--out-dir",
        "<OUT_DIR>",
        "--jobs",
        "8",
    ]


def test_training_checkpoint_context_binds_source_data_and_command() -> None:
    data = load_verified_training_data(
        FIXTURES / "factorized_corpus",
        FIXTURES / "factorized_folds.json",
    )
    command = [
        "benchmark.train_factorized_ranker",
        "--corpus-dir",
        "<CORPUS_DIR>",
        "--fold-manifest",
        "<FOLD_MANIFEST>",
        "--out-dir",
        "<OUT_DIR>",
        "--count-model-out",
        "<COUNT_MODEL_OUT>",
        "--candidate-model-out",
        "<CANDIDATE_MODEL_OUT>",
        "--resume",
    ]
    context = build_training_checkpoint_context(data, command)
    assert set(context) == {
        "schema_version",
        "kind",
        "source_git_commit",
        "source_files",
        "training_source_sha256",
        "dataset_sha256",
        "corpus_manifest_sha256",
        "chains_sha256",
        "fold_manifest_sha256",
        "feature_schema_sha256",
        "feature_family_order",
        "model_grid",
        "accepted_chain_count",
        "accepted_chain_id_sha256",
        "seed",
        "normalized_command",
        "versions",
    }
    assert len(context["source_git_commit"]) == 40
    assert context["normalized_command"] == command
    assert context["accepted_chain_count"] == 15
    assert context["feature_schema_sha256"] == training.feature_schema_hash()
    assert all(not path.startswith("/") for path in context["source_files"])
    assert build_training_checkpoint_context(data, command) == context
    assert build_training_checkpoint_context(data, [*command, "--seed", "37"]) != context
