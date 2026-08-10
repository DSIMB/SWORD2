"""Build a verified leakage-resistant fold manifest for the Task 8 corpus."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from benchmark.factorized_ranker.folds import (
    _load_accepted_entries,
    assign_folds,
    load_fold_manifest,
    write_fold_manifest,
)


_CORPUS_INPUT_NAMES = (
    "corpus_manifest.json",
    "chains.csv",
    "counts.csv",
    "candidates.csv",
    "rejections.csv",
)


def output_aliases_corpus_input(output: Path, corpus_dir: Path) -> bool:
    output_path = Path(output).resolve(strict=False)
    for name in _CORPUS_INPUT_NAMES:
        input_path = (Path(corpus_dir) / name).resolve(strict=False)
        if output_path == input_path:
            return True
        if output_path.exists() and input_path.exists():
            try:
                if output_path.samefile(input_path):
                    return True
            except OSError:
                pass
    return False


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--corpus-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=37)
    args = parser.parse_args(argv)

    if args.dataset != "cath17287":
        parser.error("factorized folds accept only dataset cath17287")
    if args.seed != 37:
        parser.error("the real factorized fold seed is frozen at 37")
    if output_aliases_corpus_input(args.out, args.corpus_dir):
        parser.error("fold output must not alias a corpus input")

    try:
        entries, dataset_hash, corpus_hash, chains_hash = _load_accepted_entries(
            args.corpus_dir, args.dataset
        )
        assignments = assign_folds(entries, n_folds=args.n_folds, seed=args.seed)
        manifest_hash = write_fold_manifest(
            args.out,
            assignments,
            dataset_hash,
            corpus_hash,
            chains_sha256=chains_hash,
            seed=args.seed,
        )
        load_fold_manifest(
            args.out,
            expected_sha256=manifest_hash,
            expected_dataset_sha256=dataset_hash,
            expected_corpus_sha256=corpus_hash,
            expected_chains_sha256=chains_hash,
        )
    except (OSError, ValueError) as error:
        parser.error(str(error))
    print(manifest_hash)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
