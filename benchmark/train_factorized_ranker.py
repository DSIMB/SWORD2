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
    if args.count_model_out is not None:
        parser.error("model artifact export requires Task 12")
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
    except (OSError, ValueError) as error:
        parser.error(str(error))
    for name in ("oof_predictions.csv", "cv_report.json", "ablation_report.json"):
        print(f"{name} {result.artifact_hashes[name]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
