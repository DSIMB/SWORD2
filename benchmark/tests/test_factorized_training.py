from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path

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
    cross_validate_configuration,
    fit_head,
    generate_oof,
    head_feature_spec,
    load_verified_training_data,
    retain_ablation,
    validate_retained_families,
)


FIXTURES = Path(__file__).parent / "fixtures"


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


def test_training_argv_roles_and_task12_guard_are_exact(capsys: pytest.CaptureFixture[str]) -> None:
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
                "--candidate-model-out",
                "candidate.json",
            ]
        )
    assert error.value.code == 2
    assert "model artifact export requires Task 12" in capsys.readouterr().err
