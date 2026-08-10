from __future__ import annotations

import json
import hashlib
import math
import platform
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy
import sklearn
from sklearn.ensemble import GradientBoostingClassifier

from benchmark.factorized_ranker.model_artifact import (
    ModelProvenance,
    _freeze_classifier,
    _rust_f64_literal,
    _rust_string_literal,
    freeze_classifier,
    reference_pair_from_batch,
    render_rust_models,
    validate_golden_vectors,
    verify_top_level_manifest,
    write_golden_vectors,
    write_artifacts,
)
from benchmark.factorized_ranker.corpus import feature_schema_hash
from benchmark.factorized_ranker.folds import load_fold_manifest
from benchmark.factorized_ranker.pairs import build_candidate_pairs, build_count_pairs
from benchmark.factorized_ranker.training import (
    FEATURE_FAMILY_ORDER,
    MODEL_GRID,
    OOFRow,
    _id_hash,
    _params_payload,
    _render_oof,
    _spec_payload,
    fit_head,
    head_feature_spec,
)
from benchmark.export_factorized_ranker import (
    MANIFEST_KEYS,
    export_factorized_ranker,
    normalize_export_argv,
    validate_export_manifest,
)


def _artifact(head: str) -> dict[str, object]:
    names = head_feature_spec(head, ("base",)).pair_features
    rng = np.random.default_rng(137 if head == "count" else 173)
    x = np.ascontiguousarray(rng.normal(size=(256, len(names))), dtype=np.float64)
    y = (x[:, 0] - 0.1 * x[:, 1] > 0.0).astype(np.int8)
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
    provenance = ModelProvenance(
        corpus_manifest_sha256="1" * 64,
        fold_manifest_sha256="2" * 64,
        cv_report_sha256="3" * 64,
        ablation_report_sha256="4" * 64,
        oof_predictions_sha256="5" * 64,
        feature_schema_sha256="6" * 64,
        feature_dump_binary_sha256="7" * 64,
        source_git_commit="8" * 40,
        training_command=("train", "--corpus-dir", "<CORPUS_DIR>"),
        python_version=platform.python_version(),
        numpy_version=np.__version__,
        pandas_version=pd.__version__,
        scipy_version=scipy.__version__,
        sklearn_version=sklearn.__version__,
    )
    return freeze_classifier(
        model,
        head_name=head,
        feature_names=names,
        retained_feature_families=("base",),
        selected_hyperparameters=params,
        provenance=provenance,
        audit_x=x,
    )


@pytest.fixture(scope="module")
def count_artifact() -> dict[str, object]:
    return _artifact("count")


@pytest.fixture(scope="module")
def candidate_artifact() -> dict[str, object]:
    return _artifact("candidate")


def _canonical_json(value: object) -> bytes:
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


def _write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


@pytest.fixture(scope="module")
def export_inputs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("factorized-export-inputs")
    fixture_root = Path(__file__).parent / "fixtures"
    corpus_dir = fixture_root / "factorized_corpus"
    corpus_manifest_path = corpus_dir / "corpus_manifest.json"
    fold_path = fixture_root / "factorized_folds.json"
    corpus_manifest_bytes = corpus_manifest_path.read_bytes()
    corpus_manifest = json.loads(corpus_manifest_bytes)
    fold_bytes = fold_path.read_bytes()
    folds = load_fold_manifest(fold_path)
    chains = pd.read_csv(corpus_dir / "chains.csv", keep_default_na=False)
    counts = pd.read_csv(corpus_dir / "counts.csv", keep_default_na=False)
    candidates = pd.read_csv(corpus_dir / "candidates.csv", keep_default_na=False)
    retained = ("base",)
    count_spec = head_feature_spec("count", retained)
    candidate_spec = head_feature_spec("candidate", retained)
    count_batch = build_count_pairs(
        chains,
        counts,
        shared_features=count_spec.shared_features,
        item_features=count_spec.item_features,
        seed=37,
    )
    candidate_batch = build_candidate_pairs(
        chains,
        candidates,
        shared_features=candidate_spec.shared_features,
        item_features=candidate_spec.item_features,
        seed=37,
    )
    params = MODEL_GRID[0]
    count_model = fit_head(count_batch, params, 37)
    candidate_model = fit_head(candidate_batch, params, 37)

    assignment_by_id = {value.chain_id: value for value in folds.assignments}
    chain_truth = {
        str(row.chain_id): int(row.n_true_domains)
        for row in chains[["chain_id", "n_true_domains"]].itertuples(index=False)
    }
    oof_rows: list[OOFRow] = []
    for chain_id in sorted(chain_truth):
        group = candidates[candidates["chain_id"] == chain_id].sort_values(
            ["num_domains", "canonical_delineation"]
        )
        winner = group.iloc[0]
        selected_count = int(winner["num_domains"])
        chosen_ndo = float(winner["ndo"])
        best_all = float(group["ndo"].max())
        best_selected = float(
            group[group["num_domains"] == selected_count]["ndo"].max()
        )
        assignment = assignment_by_id[chain_id]
        oof_rows.append(
            OOFRow(
                chain_id=chain_id,
                fold=assignment.fold,
                selected_count=selected_count,
                selected_candidate_id=str(winner["candidate_id"]),
                canonical_delineation=str(winner["canonical_delineation"]),
                count_borda_score=1.0,
                candidate_borda_score=1.0,
                n_true_domains=chain_truth[chain_id],
                count_correct=int(selected_count == chain_truth[chain_id]),
                ndo=chosen_ndo,
                boundary_f1_10=float(winner["boundary_f1_10"]),
                matched_dice=float(winner["matched_dice"]),
                total_regret=best_all - chosen_ndo,
                count_regret=best_all - best_selected,
                within_count_regret=best_selected - chosen_ndo,
                true_count_bin=assignment.true_count_bin,
                length_bin=assignment.length_bin,
                continuity_cohort="contiguous",
                label_cohort="unknown",
                count_tie_break="none",
                candidate_tie_break="none",
            )
        )
    oof_bytes = _render_oof(oof_rows)
    oof_path = root / "oof_predictions.csv"
    _write(oof_path, oof_bytes)
    oof_hash = hashlib.sha256(oof_bytes).hexdigest()
    training_command = (
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
    )
    versions = {
        "python_implementation": platform.python_implementation(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
    }
    family_mapping = {
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
    final = {
        "retained_families": list(retained),
        "selected_hyperparameters": {
            "count": _params_payload(params),
            "candidate": _params_payload(params),
        },
        "feature_specs": {
            "count": _spec_payload(count_spec),
            "candidate": _spec_payload(candidate_spec),
        },
        "oof_filename": "oof_predictions.csv",
        "oof_sha256": oof_hash,
        "oof_chain_count": len(oof_rows),
        "oof_cohorts": {},
        "final_fit": {
            "count": {
                "rows": len(count_batch.y),
                "features": count_batch.x.shape[1],
                "classes": [0, 1],
            },
            "candidate": {
                "rows": len(candidate_batch.y),
                "features": candidate_batch.x.shape[1],
                "classes": [0, 1],
            },
        },
    }
    common = {
        "schema_version": 1,
        "seed": 37,
        "evidence_role": "development_model_selection_not_locked_acceptance",
        "normalized_command": list(training_command),
        "dataset": "cath17287",
        "dataset_sha256": corpus_manifest["dataset_sha256"],
        "corpus_manifest_sha256": hashlib.sha256(corpus_manifest_bytes).hexdigest(),
        "chains_sha256": corpus_manifest["tables"]["chains"]["sha256"],
        "fold_manifest_sha256": hashlib.sha256(fold_bytes).hexdigest(),
        "feature_schema_version": 1,
        "feature_schema_sha256": feature_schema_hash(),
        "binary_sha256": corpus_manifest["binary_sha256"],
        "git_commit": corpus_manifest["git_commit"],
        "accepted_chain_count": len(oof_rows),
        "accepted_chain_id_sha256": _id_hash([row.chain_id for row in oof_rows]),
        "feature_family_order": list(FEATURE_FAMILY_ORDER),
        "feature_family_mapping": family_mapping,
        "model_grid": [_params_payload(value) for value in MODEL_GRID],
        "versions": versions,
        "final": final,
    }
    cv_path = root / "cv_report.json"
    ablation_path = root / "ablation_report.json"
    cv_bytes = _canonical_json(
        {
            **common,
            "grid_history": [
                {"family": family, "count": {}, "candidate": {}}
                for family in FEATURE_FAMILY_ORDER
            ],
        }
    )
    ablation_bytes = _canonical_json(
        {
            **common,
            "stages": [
                {
                    "family": family,
                    "retained": family == "base",
                    "oof_sha256": oof_hash,
                }
                for family in FEATURE_FAMILY_ORDER
            ],
        }
    )
    _write(cv_path, cv_bytes)
    _write(ablation_path, ablation_bytes)
    provenance = ModelProvenance(
        corpus_manifest_sha256=hashlib.sha256(corpus_manifest_bytes).hexdigest(),
        fold_manifest_sha256=hashlib.sha256(fold_bytes).hexdigest(),
        cv_report_sha256=hashlib.sha256(cv_bytes).hexdigest(),
        ablation_report_sha256=hashlib.sha256(ablation_bytes).hexdigest(),
        oof_predictions_sha256=oof_hash,
        feature_schema_sha256=feature_schema_hash(),
        feature_dump_binary_sha256=corpus_manifest["binary_sha256"],
        source_git_commit=corpus_manifest["git_commit"],
        training_command=training_command,
        python_version=platform.python_version(),
        numpy_version=np.__version__,
        pandas_version=pd.__version__,
        scipy_version=scipy.__version__,
        sklearn_version=sklearn.__version__,
    )
    from benchmark.factorized_ranker.training import CorpusTables

    corpus = CorpusTables(chains, counts, candidates)
    count_artifact = _freeze_classifier(
        count_model,
        "count",
        count_batch.feature_names,
        retained,
        params,
        provenance,
        count_batch.x,
        reference_pair=reference_pair_from_batch(
            corpus, count_batch, "count", retained
        ),
    )
    candidate_artifact = _freeze_classifier(
        candidate_model,
        "candidate",
        candidate_batch.feature_names,
        retained,
        params,
        provenance,
        candidate_batch.x,
        reference_pair=reference_pair_from_batch(
            corpus, candidate_batch, "candidate", retained
        ),
    )
    count_model_path = root / "count.json"
    candidate_model_path = root / "candidate.json"
    write_artifacts(
        count_model_path, count_artifact, candidate_model_path, candidate_artifact
    )
    baseline = root / "standalone.bin"
    baseline.write_bytes(b"opaque standalone baseline\n")
    return {
        "count_model": count_model_path,
        "candidate_model": candidate_model_path,
        "corpus_manifest": corpus_manifest_path,
        "fold_manifest": fold_path,
        "cv_report": cv_path,
        "ablation_report": ablation_path,
        "oof_predictions": oof_path,
        "standalone_baseline": baseline,
    }


def test_rust_literals_preserve_bits_and_escape_source() -> None:
    assert _rust_f64_literal(-0.0).startswith("-0.0")
    assert math.copysign(1.0, float(_rust_f64_literal(-0.0).removesuffix("_f64"))) < 0
    escaped = _rust_string_literal('"\\\n\r\t\0\x01é')
    assert escaped == '"\\"\\\\\\n\\r\\t\\0\\u{1}é"'


def test_generated_rust_is_ordered_formatted_and_compiles(
    tmp_path: Path, count_artifact, candidate_artifact
) -> None:
    source = render_rust_models(count_artifact, candidate_artifact)
    declarations = (
        "MODEL_INPUT_DTYPE",
        "MODEL_THRESHOLD_POLICY",
        "RETAINED_FEATURE_FAMILIES",
        "COUNT_FEATURE_NAMES",
        "COUNT_NODES",
        "COUNT_TREES",
        "COUNT_MODEL",
        "CANDIDATE_FEATURE_NAMES",
        "CANDIDATE_NODES",
        "CANDIDATE_TREES",
        "CANDIDATE_MODEL",
    )
    positions = [source.index(name) for name in declarations]
    assert positions == sorted(positions)
    assert "u16::MAX" in source
    assert "timestamp" not in source.casefold()
    generated = tmp_path / "generated.rs"
    generated.write_text(source, encoding="utf-8", newline="")
    subprocess.run(
        ["rustfmt", "--edition", "2021", "--check", str(generated)],
        check=True,
        capture_output=True,
        text=True,
    )
    model = tmp_path / "model.rs"
    model.write_text(
        """
pub(crate) struct StaticNode { pub feature: u16, pub threshold: f64, pub left: u16, pub right: u16, pub leaf_value: f64 }
pub(crate) struct StaticTree { pub root: u16 }
pub(crate) struct StaticBoostedModel { pub schema_version: u32, pub feature_names: &'static [&'static str], pub initial_log_odds: f64, pub learning_rate: f64, pub trees: &'static [StaticTree], pub nodes: &'static [StaticNode] }
""".lstrip(),
        encoding="utf-8",
    )
    harness = tmp_path / "lib.rs"
    harness.write_text(
        '#[path = "model.rs"] mod model;\n#[path = "generated.rs"] mod generated;\n',
        encoding="utf-8",
    )
    subprocess.run(
        ["rustc", "--edition", "2021", "--crate-type", "lib", str(harness)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )


def test_golden_vectors_have_exact_cases_and_tagged_errors(
    tmp_path: Path, count_artifact, candidate_artifact
) -> None:
    path = tmp_path / "golden.json"
    digest = write_golden_vectors(path, count_artifact, candidate_artifact)
    assert len(digest) == 64
    payload = json.loads(path.read_bytes())
    validate_golden_vectors(payload, count_artifact, candidate_artifact)
    for head in ("count", "candidate"):
        section = payload["heads"][head]
        assert [case["kind"] for case in section["cases"]] == [
            "realistic",
            "threshold_below",
            "threshold_equal",
            "threshold_above",
        ]
        assert [case["expected_equality_branch"] for case in section["cases"]][2] == "left"
        assert [error["kind"] for error in section["errors"]] == [
            "nan_at",
            "positive_infinity_at",
            "negative_infinity_at",
            "length_mismatch",
            "duplicate_item_id",
            "float32_overflow_at",
        ]
        for case in section["cases"]:
            assert len(case["forward_pair_vector"]) == len(
                case["forward_canonical_vector"]
            )
            assert 0.0 <= case["symmetrized_probability"] <= 1.0
            assert sum(case["borda_scores"].values()) == pytest.approx(1.0)
    assert b"NaN" not in path.read_bytes()
    assert b"Infinity" not in path.read_bytes()


def _export(
    inputs: dict[str, Path],
    output_root: Path,
    *,
    cv_report: Path | None = None,
):
    argv = [
        "benchmark.export_factorized_ranker",
        "--count-model",
        str(inputs["count_model"]),
        "--candidate-model",
        str(inputs["candidate_model"]),
        "--corpus-manifest",
        str(inputs["corpus_manifest"]),
        "--fold-manifest",
        str(inputs["fold_manifest"]),
        "--cv-report",
        str(cv_report or inputs["cv_report"]),
        "--ablation-report",
        str(inputs["ablation_report"]),
        "--oof-predictions",
        str(inputs["oof_predictions"]),
        "--standalone-baseline",
        str(inputs["standalone_baseline"]),
        "--rust-out",
        str(output_root / "generated.rs"),
        "--golden-out",
        str(output_root / "golden.json"),
        "--manifest-out",
        str(output_root / "manifest.json"),
    ]
    return export_factorized_ranker(
        count_model=inputs["count_model"],
        candidate_model=inputs["candidate_model"],
        corpus_manifest=inputs["corpus_manifest"],
        fold_manifest=inputs["fold_manifest"],
        cv_report=cv_report or inputs["cv_report"],
        ablation_report=inputs["ablation_report"],
        oof_predictions=inputs["oof_predictions"],
        standalone_baseline=inputs["standalone_baseline"],
        rust_out=output_root / "generated.rs",
        golden_out=output_root / "golden.json",
        manifest_out=output_root / "manifest.json",
        normalized_command=normalize_export_argv(argv),
    )


def test_export_path_roles_are_exact() -> None:
    first = normalize_export_argv(
        [
            "export",
            "--count-model=/one/count.json",
            "--rust-out",
            "/one/generated.rs",
            "--seed",
            "37",
        ]
    )
    second = normalize_export_argv(
        [
            "export",
            "--count-model=/two/count.json",
            "--rust-out",
            "/two/generated.rs",
            "--seed",
            "37",
        ]
    )
    assert first == second == [
        "export",
        "--count-model=<COUNT_MODEL>",
        "--rust-out",
        "<RUST_OUT>",
        "--seed",
        "37",
    ]


def test_complete_export_is_byte_identical_and_manifest_bound(
    tmp_path: Path, export_inputs: dict[str, Path]
) -> None:
    first = _export(export_inputs, tmp_path / "first")
    second = _export(export_inputs, tmp_path / "second")
    for name in ("generated.rs", "golden.json", "manifest.json"):
        assert (tmp_path / "first" / name).read_bytes() == (
            tmp_path / "second" / name
        ).read_bytes()
    manifest = json.loads((tmp_path / "first" / "manifest.json").read_bytes())
    assert set(manifest) == MANIFEST_KEYS
    validate_export_manifest(manifest)
    assert first.manifest == second.manifest == manifest
    assert first.generated_rust_sha256 == hashlib.sha256(
        (tmp_path / "first" / "generated.rs").read_bytes()
    ).hexdigest()
    assert first.golden_sha256 == hashlib.sha256(
        (tmp_path / "first" / "golden.json").read_bytes()
    ).hexdigest()
    assert "manifest_sha256" not in manifest
    encoded = (tmp_path / "first" / "manifest.json").read_text(encoding="utf-8")
    assert "/home/" not in encoded and "/tmp/" not in encoded

    install_root = tmp_path / "installed"
    model_dir = install_root / "benchmark" / "models"
    installed = {
        model_dir / "cath17287_factorized_corpus_v1_manifest.json": export_inputs[
            "corpus_manifest"
        ],
        model_dir / "cath17287_factorized_folds_v1.json": export_inputs[
            "fold_manifest"
        ],
        model_dir / "factorized_ranker_v1_cv.json": export_inputs["cv_report"],
        model_dir / "factorized_ranker_v1_ablations.json": export_inputs[
            "ablation_report"
        ],
        model_dir / "factorized_ranker_v1_oof.csv": export_inputs[
            "oof_predictions"
        ],
        model_dir / "factorized_count_v1.json": export_inputs["count_model"],
        model_dir / "factorized_candidate_v1.json": export_inputs[
            "candidate_model"
        ],
        model_dir / "factorized_ranker_v1_golden.json": tmp_path
        / "first"
        / "golden.json",
        model_dir / "cath663_standalone_structural_baseline.csv": export_inputs[
            "standalone_baseline"
        ],
        install_root
        / "sword2-lib/src/sword/factorized_ranker/generated_model.rs": tmp_path
        / "first"
        / "generated.rs",
        model_dir / "factorized_ranker_v1_manifest.json": tmp_path
        / "first"
        / "manifest.json",
    }
    for destination, source in installed.items():
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(source.read_bytes())
    assert verify_top_level_manifest(
        model_dir / "factorized_ranker_v1_manifest.json"
    ) == manifest


def test_input_failure_leaves_all_existing_outputs_untouched(
    tmp_path: Path, export_inputs: dict[str, Path]
) -> None:
    output = tmp_path / "output"
    output.mkdir()
    sentinels = {
        output / "generated.rs": b"old rust\n",
        output / "golden.json": b"old golden\n",
        output / "manifest.json": b"old manifest\n",
    }
    for path, data in sentinels.items():
        path.write_bytes(data)
    bad_cv = tmp_path / "bad-cv.json"
    payload = json.loads(export_inputs["cv_report"].read_bytes())
    payload["seed"] = 38
    bad_cv.write_bytes(_canonical_json(payload))
    with pytest.raises(ValueError):
        _export(export_inputs, output, cv_report=bad_cv)
    assert {path: path.read_bytes() for path in sentinels} == sentinels
