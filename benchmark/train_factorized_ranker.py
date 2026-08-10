"""Train factorized structural heads with grouped development-only ablations."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from benchmark.factorized_ranker.training import (
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
    parser.add_argument("--count-model-out", type=Path, default=None)
    parser.add_argument("--candidate-model-out", type=Path, default=None)
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
        result = run_grouped_training(
            data,
            args.out_dir,
            normalized,
            seed=args.seed,
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
