"""Create or verify the frozen input-only factorized eligibility manifest."""

from __future__ import annotations

import argparse
import hashlib
import stat
import sys
from pathlib import Path

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark.factorized_ranker.eligibility import (
    create_eligibility_manifest,
    stable_file_bytes,
    verify_eligibility_manifest,
    write_eligibility_manifest,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    create = subparsers.add_parser("create")
    create.add_argument("--dataset-metadata", type=Path, required=True)
    create.add_argument("--cache-dir", type=Path, required=True)
    create.add_argument("--runtime-manifest", type=Path, required=True)
    create.add_argument("--binary", type=Path, required=True)
    create.add_argument("--repo-root", type=Path, required=True)
    create.add_argument("--jobs", type=int, required=True)
    create.add_argument("--out", type=Path, required=True)

    verify = subparsers.add_parser("verify")
    verify.add_argument("--manifest", type=Path, required=True)
    verify.add_argument("--dataset-metadata", type=Path, required=True)
    verify.add_argument("--cache-dir", type=Path, required=True)
    verify.add_argument("--runtime-manifest", type=Path, required=True)
    verify.add_argument("--binary", type=Path, required=True)
    verify.add_argument("--repo-root", type=Path, required=True)
    return parser.parse_args(argv)


def _regular_file(path: Path, description: str) -> tuple[Path, tuple[int, int]]:
    info = Path(path).lstat()
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
        raise ValueError(f"{description} is not a regular nonsymlink file")
    return Path(path).resolve(strict=True), (info.st_dev, info.st_ino)


def _regular_directory(path: Path, description: str) -> Path:
    info = Path(path).lstat()
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise ValueError(f"{description} is not a regular nonsymlink directory")
    return Path(path).resolve(strict=True)


def _validate_authority_paths(args: argparse.Namespace) -> None:
    files = {
        "dataset metadata": args.dataset_metadata,
        "runtime manifest": args.runtime_manifest,
        "binary": args.binary,
    }
    identities: set[tuple[int, int]] = set()
    for description, path in files.items():
        _resolved, identity = _regular_file(path, description)
        if identity in identities:
            raise ValueError("eligibility authority inputs alias")
        identities.add(identity)
    cache = _regular_directory(args.cache_dir, "cache directory")
    _regular_directory(cache / "chains", "chain cache")
    repo = _regular_directory(args.repo_root, "repository root")
    if not (repo / "benchmark").is_dir():
        raise ValueError("repository root is invalid")


def _validate_create_paths(args: argparse.Namespace) -> None:
    _validate_authority_paths(args)
    if type(args.jobs) is not int or args.jobs < 1:
        raise ValueError("eligibility jobs must be a positive integer")
    output = Path(args.out).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    parent = _regular_directory(output.parent, "output parent")
    output = parent / output.name
    for source in (args.dataset_metadata, args.runtime_manifest, args.binary):
        if output == Path(source).resolve(strict=True):
            raise ValueError("eligibility output aliases an authority input")


def _validate_verify_paths(args: argparse.Namespace) -> None:
    _validate_authority_paths(args)
    _regular_file(args.manifest, "eligibility manifest")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.command == "create":
            _validate_create_paths(args)
            payload = create_eligibility_manifest(
                dataset_metadata=args.dataset_metadata,
                cache_dir=args.cache_dir,
                runtime_manifest=args.runtime_manifest,
                binary=args.binary,
                repo_root=args.repo_root,
                jobs=args.jobs,
            )
            digest = write_eligibility_manifest(args.out, payload)
        else:
            _validate_verify_paths(args)
            verify_eligibility_manifest(
                args.manifest,
                dataset_metadata=args.dataset_metadata,
                cache_dir=args.cache_dir,
                runtime_manifest=args.runtime_manifest,
                binary=args.binary,
                repo_root=args.repo_root,
            )
            digest = hashlib.sha256(stable_file_bytes(args.manifest)).hexdigest()
    except (OSError, ValueError) as error:
        print(f"invalid eligibility freeze: {error}", file=sys.stderr)
        return 2
    print(digest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
