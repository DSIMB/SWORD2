"""Create or verify the frozen factorized-selector runtime identity."""

from __future__ import annotations

import argparse
import stat
import sys
from pathlib import Path

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark.factorized_ranker.runtime_freeze import (
    create_runtime_freeze,
    verify_runtime_freeze,
    write_runtime_freeze,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    create = subparsers.add_parser("create")
    create.add_argument("--model-manifest", type=Path, required=True)
    create.add_argument("--corpus-dir", type=Path, required=True)
    create.add_argument("--fold-manifest", type=Path, required=True)
    create.add_argument("--oof-predictions", type=Path, required=True)
    create.add_argument("--binary", type=Path, required=True)
    create.add_argument("--repo-root", type=Path, required=True)
    create.add_argument("--out", type=Path, required=True)

    verify = subparsers.add_parser("verify")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--binary", type=Path, required=True)
    verify.add_argument("--repo-root", type=Path, required=True)
    return parser.parse_args(argv)


def _validate_create_paths(args: argparse.Namespace) -> None:
    file_inputs = {
        "model manifest": args.model_manifest,
        "fold manifest": args.fold_manifest,
        "OOF predictions": args.oof_predictions,
        "binary": args.binary,
    }
    resolved_files: dict[str, Path] = {}
    inodes: set[tuple[int, int]] = set()
    for role, raw_path in file_inputs.items():
        info = Path(raw_path).lstat()
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
            raise ValueError(f"{role} is not a regular nonsymlink file")
        resolved = Path(raw_path).resolve(strict=True)
        identity = (info.st_dev, info.st_ino)
        if identity in inodes:
            raise ValueError("runtime-freeze inputs alias")
        inodes.add(identity)
        resolved_files[role] = resolved

    corpus_info = Path(args.corpus_dir).lstat()
    if stat.S_ISLNK(corpus_info.st_mode) or not stat.S_ISDIR(corpus_info.st_mode):
        raise ValueError("corpus input is not a regular nonsymlink directory")
    corpus = Path(args.corpus_dir).resolve(strict=True)
    repo = Path(args.repo_root).resolve(strict=True)
    output = Path(args.out).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    parent = output.parent.resolve(strict=True)
    output = parent / output.name
    if output in resolved_files.values():
        raise ValueError("runtime-freeze output aliases an input")
    try:
        output.relative_to(corpus)
    except ValueError:
        pass
    else:
        raise ValueError("runtime-freeze output is inside the read-only corpus")
    try:
        corpus.relative_to(output)
    except ValueError:
        pass
    else:
        raise ValueError("runtime-freeze output contains the read-only corpus")
    if not (repo / "benchmark").is_dir():
        raise ValueError("repository root is invalid")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.command == "create":
            _validate_create_paths(args)
            payload = create_runtime_freeze(
                model_manifest=args.model_manifest,
                corpus_dir=args.corpus_dir,
                fold_manifest=args.fold_manifest,
                oof_predictions=args.oof_predictions,
                binary=args.binary,
                repo_root=args.repo_root,
            )
            digest = write_runtime_freeze(args.out, payload)
        else:
            verify_runtime_freeze(
                args.manifest,
                binary=args.binary,
                repo_root=args.repo_root,
            )
            from benchmark.factorized_ranker.runtime_freeze import sha256_file

            digest = sha256_file(args.manifest)
    except (OSError, ValueError) as error:
        print(f"invalid runtime freeze: {error}", file=sys.stderr)
        return 2
    print(digest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
