"""Train factorized structural heads with grouped development-only ablations."""

from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path
from typing import Sequence

from benchmark.factorized_ranker.training import (
    FEATURE_FAMILY_ORDER,
    TrainingCheckpointStore,
    VerifiedTrainingData,
    build_training_checkpoint_context,
    load_verified_training_data,
    run_grouped_training,
)


_PATH_ROLES = {
    "--corpus-dir": "<CORPUS_DIR>",
    "--fold-manifest": "<FOLD_MANIFEST>",
    "--out-dir": "<OUT_DIR>",
    "--count-model-out": "<COUNT_MODEL_OUT>",
    "--candidate-model-out": "<CANDIDATE_MODEL_OUT>",
}


def _jobs_argument(value: str) -> int:
    try:
        jobs = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("jobs must be an integer in 1..=8") from error
    if not 1 <= jobs <= 8:
        raise argparse.ArgumentTypeError("jobs must be an integer in 1..=8")
    return jobs


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _reuse_completed_outputs(
    data: VerifiedTrainingData,
    out_dir: Path,
    count_model_path: Path,
    candidate_model_path: Path,
    normalized_command: Sequence[str],
    checkpoint_store: TrainingCheckpointStore,
) -> tuple[dict[str, str], tuple[str, str]] | None:
    accepted_ids = {assignment.chain_id for assignment in data.assignments}
    state = checkpoint_store.load_stage(tuple(sorted(accepted_ids)))
    if (
        state is None
        or state.completed_family_index != len(FEATURE_FAMILY_ORDER) - 1
    ):
        return None
    report_paths = {
        "oof_predictions.csv": Path(out_dir) / "oof_predictions.csv",
        "cv_report.json": Path(out_dir) / "cv_report.json",
        "ablation_report.json": Path(out_dir) / "ablation_report.json",
    }
    required = [
        *report_paths.values(),
        Path(count_model_path),
        Path(candidate_model_path),
    ]
    if not all(path.exists() for path in required):
        return None

    from benchmark.export_factorized_ranker import _load_oof, _validate_reports
    from benchmark.factorized_ranker.model_artifact import (
        load_artifact,
        validate_artifacts,
    )

    count_artifact = load_artifact(count_model_path, expected_head="count")
    candidate_artifact = load_artifact(
        candidate_model_path, expected_head="candidate"
    )
    validate_artifacts(count_artifact, candidate_artifact)
    if (
        list(state.retained_families)
        != count_artifact["retained_feature_families"]
        or state.retained_oof.sha256
        != count_artifact["oof_predictions_sha256"]
    ):
        raise ValueError("completed checkpoint state disagrees with frozen models")
    for artifact, params, head in (
        (count_artifact, state.count_params, "count"),
        (candidate_artifact, state.candidate_params, "candidate"),
    ):
        if (
            artifact["n_estimators"] != params.n_estimators
            or artifact["learning_rate"] != params.learning_rate
            or artifact["min_samples_leaf"] != params.min_samples_leaf
            or artifact["max_depth"] != params.max_depth
        ):
            raise ValueError(
                f"completed checkpoint {head} hyperparameters disagree with model"
            )
    expected_provenance = {
        "training_command": list(normalized_command),
        "corpus_manifest_sha256": data.corpus_manifest_sha256,
        "fold_manifest_sha256": data.fold_manifest_sha256,
        "feature_dump_binary_sha256": data.corpus_manifest["binary_sha256"],
        "source_git_commit": data.corpus_manifest["git_commit"],
    }
    for field, expected in expected_provenance.items():
        if count_artifact[field] != expected:
            raise ValueError(f"completed model provenance field {field} mismatch")
    _validate_reports(
        report_paths["cv_report.json"],
        report_paths["ablation_report.json"],
        count_artifact,
        candidate_artifact,
        data.corpus_manifest,
        data.corpus_manifest_sha256,
        data.fold_manifest_sha256,
        data.chains_sha256,
        accepted_ids,
        data.corpus,
    )
    _load_oof(
        report_paths["oof_predictions.csv"],
        str(count_artifact["oof_predictions_sha256"]),
        accepted_ids,
        {assignment.chain_id: assignment for assignment in data.assignments},
        data.corpus,
    )
    report_hashes = {
        name: _sha256(path) for name, path in report_paths.items()
    }
    model_hashes = (_sha256(count_model_path), _sha256(candidate_model_path))
    return report_hashes, model_hashes


def normalize_training_argv(argv: Sequence[str]) -> list[str]:
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


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus-dir", type=Path, required=True)
    parser.add_argument("--fold-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=37)
    parser.add_argument(
        "--jobs",
        type=_jobs_argument,
        default=1,
        help="independent fold workers (1..=8)",
    )
    parser.add_argument("--count-model-out", type=Path, default=None)
    parser.add_argument("--candidate-model-out", type=Path, default=None)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="reuse exact-context fold and completed-stage checkpoints",
    )
    args = parser.parse_args(argv)
    if args.seed != 37:
        parser.error("factorized grouped training seed is frozen at 37")
    if (args.count_model_out is None) != (args.candidate_model_out is None):
        parser.error("count and candidate model outputs must be supplied together")
    raw_argv = [
        "benchmark.train_factorized_ranker",
        *(sys.argv[1:] if argv is None else argv),
    ]
    try:
        normalized = normalize_training_argv(raw_argv)
        data = load_verified_training_data(args.corpus_dir, args.fold_manifest)
        checkpoint_store = (
            TrainingCheckpointStore(
                args.out_dir / ".factorized_training_checkpoints",
                build_training_checkpoint_context(data, normalized),
            )
            if args.resume
            else None
        )
        completed = (
            _reuse_completed_outputs(
                data,
                args.out_dir,
                args.count_model_out,
                args.candidate_model_out,
                normalized,
                checkpoint_store,
            )
            if checkpoint_store is not None
            and args.count_model_out is not None
            and args.candidate_model_out is not None
            else None
        )
        if completed is not None:
            completed_reports, completed_models = completed
            for name in (
                "oof_predictions.csv",
                "cv_report.json",
                "ablation_report.json",
            ):
                print(f"{name} {completed_reports[name]}")
            print(f"count_model.json {completed_models[0]}")
            print(f"candidate_model.json {completed_models[1]}")
            return 0
        result = run_grouped_training(
            data,
            args.out_dir,
            normalized,
            seed=args.seed,
            jobs=args.jobs,
            checkpoint_store=checkpoint_store,
        )
        model_hashes: tuple[str, str] | None = None
        if args.count_model_out is not None and args.candidate_model_out is not None:
            from benchmark.factorized_ranker.model_artifact import (
                ModelProvenance,
                _freeze_classifier,
                reference_pair_from_batch,
                write_artifacts,
            )

            versions = result.cv_report["versions"]
            if not isinstance(versions, dict):
                raise ValueError("training report versions are malformed")
            provenance = ModelProvenance(
                corpus_manifest_sha256=data.corpus_manifest_sha256,
                fold_manifest_sha256=data.fold_manifest_sha256,
                cv_report_sha256=result.artifact_hashes["cv_report.json"],
                ablation_report_sha256=result.artifact_hashes[
                    "ablation_report.json"
                ],
                oof_predictions_sha256=result.artifact_hashes[
                    "oof_predictions.csv"
                ],
                feature_schema_sha256=str(
                    result.cv_report["feature_schema_sha256"]
                ),
                feature_dump_binary_sha256=str(data.corpus_manifest["binary_sha256"]),
                source_git_commit=str(data.corpus_manifest["git_commit"]),
                training_command=tuple(normalized),
                python_version=str(versions["python"]),
                numpy_version=str(versions["numpy"]),
                pandas_version=str(versions["pandas"]),
                scipy_version=str(versions["scipy"]),
                sklearn_version=str(versions["scikit_learn"]),
            )
            count_reference = reference_pair_from_batch(
                data.corpus,
                result.count_batch,
                "count",
                result.retained_families,
            )
            candidate_reference = reference_pair_from_batch(
                data.corpus,
                result.candidate_batch,
                "candidate",
                result.retained_families,
            )
            count_artifact = _freeze_classifier(
                result.count_model,
                "count",
                result.count_batch.feature_names,
                result.retained_families,
                result.count_params,
                provenance,
                result.count_batch.x,
                seed=args.seed,
                reference_pair=count_reference,
            )
            candidate_artifact = _freeze_classifier(
                result.candidate_model,
                "candidate",
                result.candidate_batch.feature_names,
                result.retained_families,
                result.candidate_params,
                provenance,
                result.candidate_batch.x,
                seed=args.seed,
                reference_pair=candidate_reference,
            )
            model_hashes = write_artifacts(
                args.count_model_out,
                count_artifact,
                args.candidate_model_out,
                candidate_artifact,
            )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    for name in ("oof_predictions.csv", "cv_report.json", "ablation_report.json"):
        print(f"{name} {result.artifact_hashes[name]}")
    if model_hashes is not None:
        print(f"count_model.json {model_hashes[0]}")
        print(f"candidate_model.json {model_hashes[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
