from __future__ import annotations

import copy
import hashlib
import os
from pathlib import Path

import pytest

import benchmark.factorized_ranker.eligibility as eligibility_module
import benchmark.freeze_factorized_eligibility as eligibility_cli
from benchmark.factorized_ranker.eligibility import (
    DatasetIdentity,
    build_eligibility_manifest,
    create_eligibility_manifest,
    inspect_quality_record,
    load_quality_record_bytes,
    read_dataset_identity_prefixes,
    validate_eligibility_manifest,
    verify_eligibility_manifest,
    write_eligibility_manifest,
)
from benchmark.factorized_ranker.runtime_freeze import canonical_json_bytes


def _report(*, chain_id: str, eligible: bool) -> bytes:
    missing = (
        []
        if eligible
        else [
            {
                "author_residue_number": 7,
                "chain_id": chain_id,
                "missing_atoms": ["O"],
            }
        ]
    )
    return canonical_json_bytes(
        {
            "candidate_residue_count": 10,
            "chain_id": chain_id,
            "complete_backbone_residue_count": 10 if eligible else 9,
            "eligible": eligible,
            "incomplete_residues": missing,
            "policy": "strict_complete_backbone_v1",
            "reason_code": None,
            "schema_version": 1,
            "structural_coverage": 1.0 if eligible else 0.9,
        }
    )


def _synthetic_runtime() -> dict[str, object]:
    return {
        "runtime_source_git_commit": "a" * 40,
        "binary_sha256": "b" * 64,
    }


def _manifest_fixture(tmp_path: Path) -> dict[str, object]:
    chain_root = tmp_path / "chains"
    chain_root.mkdir()
    (chain_root / "oneA.pdb").write_bytes(b"first structure\n")
    (chain_root / "twoB.pdb").write_bytes(b"second structure\n")
    identities = [
        DatasetIdentity("1aaa", "oneA", "A"),
        DatasetIdentity("2bbb", "twoB", "B"),
    ]
    reports = {
        "oneA": _report(chain_id="A", eligible=True),
        "twoB": _report(chain_id="B", eligible=False),
    }
    return build_eligibility_manifest(
        identities=identities,
        chain_root=chain_root,
        dataset="synthetic",
        dataset_sha256="d" * 64,
        runtime=_synthetic_runtime(),
        runtime_manifest_sha256="c" * 64,
        binary_sha256="b" * 64,
        inspector=lambda identity, _path: reports[identity.entry_id],
        jobs=2,
    )


def test_build_manifest_is_order_stable_and_partitions_all_ids(tmp_path: Path) -> None:
    manifest = _manifest_fixture(tmp_path)

    assert manifest["dataset_ids"] == ["oneA", "twoB"]
    assert manifest["eligible_ids"] == ["oneA"]
    assert manifest["ineligible_ids"] == ["twoB"]
    assert manifest["eligible_count"] + manifest["ineligible_count"] == 2
    assert manifest["ineligibility_reason_counts"] == {"incomplete_backbone": 1}
    validate_eligibility_manifest(manifest)


def test_structure_evidence_is_streamed_without_retaining_whole_pdb_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    chain_root = tmp_path / "chains"
    chain_root.mkdir()
    (chain_root / "oneA.pdb").write_bytes(b"structure\n")
    original = eligibility_module.stable_file_bytes

    def reject_whole_pdb_read(path: Path) -> bytes:
        if Path(path).suffix == ".pdb":
            raise AssertionError("PDB evidence must be hashed as a stream")
        return original(path)

    monkeypatch.setattr(eligibility_module, "stable_file_bytes", reject_whole_pdb_read)
    manifest = build_eligibility_manifest(
        identities=[DatasetIdentity("1aaa", "oneA", "A")],
        chain_root=chain_root,
        dataset="synthetic",
        dataset_sha256="d" * 64,
        runtime=_synthetic_runtime(),
        runtime_manifest_sha256="c" * 64,
        binary_sha256="b" * 64,
        inspector=lambda _identity, _path: _report(chain_id="A", eligible=True),
        jobs=1,
    )
    assert manifest["structure_sha256s"] == {
        "oneA": hashlib.sha256(b"structure\n").hexdigest()
    }


def test_identity_prefix_parser_never_retains_label_remainder(tmp_path: Path) -> None:
    sentinel = "DO_NOT_PARSE_THIS_CATH_LABEL"
    data = f"1aaa,oneA,A,2,10,99,1-99:{sentinel}\n".encode("ascii")
    path = tmp_path / "dataset.csv"
    path.write_bytes(data)

    identities, digest = read_dataset_identity_prefixes(path)

    assert identities == [DatasetIdentity("1aaa", "oneA", "A")]
    assert sentinel not in repr(identities)
    assert digest == hashlib.sha256(data).hexdigest()


@pytest.mark.parametrize(
    "data",
    [
        b"1aaa,oneA,A\n",
        b"1aaa,oneA,AA,ignored\n",
        b"1aaa, oneA,A,ignored\n",
        b"1aaa,oneA,A,ignored\n2bbb,oneA,B,ignored\n",
        b"\xff,oneA,A,ignored\n",
    ],
)
def test_identity_prefix_parser_rejects_truncation_noncanonical_and_duplicates(
    tmp_path: Path, data: bytes
) -> None:
    path = tmp_path / "dataset.csv"
    path.write_bytes(data)
    with pytest.raises(ValueError):
        read_dataset_identity_prefixes(path)


def test_quality_record_validation_accepts_exact_rust_schema() -> None:
    record = load_quality_record_bytes(
        _report(chain_id="A", eligible=False), expected_chain_id="A"
    )
    assert record["eligible"] is False
    assert record["incomplete_residues"][0]["missing_atoms"] == ["O"]


def test_quality_record_rejects_duplicate_keys_noncanonical_bytes_and_constants() -> None:
    valid = _report(chain_id="A", eligible=True)
    with pytest.raises(ValueError):
        load_quality_record_bytes(valid.replace(b'{"candidate', b'{"eligible":true,"candidate'))
    with pytest.raises(ValueError):
        load_quality_record_bytes(valid[:-1] + b" \n")
    with pytest.raises(ValueError):
        load_quality_record_bytes(valid.replace(b"1.0", b"NaN"))


def test_quality_record_rejects_each_semantic_inconsistency() -> None:
    base = {
        "candidate_residue_count": 10,
        "chain_id": "A",
        "complete_backbone_residue_count": 9,
        "eligible": False,
        "incomplete_residues": [
            {
                "author_residue_number": 7,
                "chain_id": "A",
                "missing_atoms": ["CA", "O"],
            }
        ],
        "policy": "strict_complete_backbone_v1",
        "reason_code": None,
        "schema_version": 1,
        "structural_coverage": 0.9,
    }
    mutations = [
        lambda row: row.update(policy="other"),
        lambda row: row.update(schema_version=2),
        lambda row: row.update(candidate_residue_count=True),
        lambda row: row.update(complete_backbone_residue_count=10),
        lambda row: row.update(structural_coverage=0.8),
        lambda row: row.update(eligible=True),
        lambda row: row.update(reason_code="empty_candidate_population"),
        lambda row: row["incomplete_residues"][0].update(chain_id="B"),
        lambda row: row["incomplete_residues"][0].update(missing_atoms=["O", "CA"]),
        lambda row: row["incomplete_residues"][0].update(missing_atoms=["CA", "CA"]),
    ]
    for mutate in mutations:
        payload = copy.deepcopy(base)
        mutate(payload)
        with pytest.raises(ValueError):
            load_quality_record_bytes(canonical_json_bytes(payload))

    with pytest.raises(ValueError):
        load_quality_record_bytes(
            _report(chain_id="A", eligible=True), expected_chain_id="B"
        )


def test_empty_candidate_record_has_one_stable_reason() -> None:
    payload = {
        "candidate_residue_count": 0,
        "chain_id": "A",
        "complete_backbone_residue_count": 0,
        "eligible": False,
        "incomplete_residues": [],
        "policy": "strict_complete_backbone_v1",
        "reason_code": "empty_candidate_population",
        "schema_version": 1,
        "structural_coverage": 0.0,
    }
    assert load_quality_record_bytes(canonical_json_bytes(payload)) == payload
    payload["reason_code"] = None
    with pytest.raises(ValueError):
        load_quality_record_bytes(canonical_json_bytes(payload))


def test_manifest_validation_rejects_one_field_tampering(tmp_path: Path) -> None:
    valid = _manifest_fixture(tmp_path)
    expected_dataset_ids = tuple(valid["dataset_ids"])

    def duplicate_id(row: dict[str, object]) -> None:
        row["dataset_ids"].append("oneA")

    def structure_hash(row: dict[str, object]) -> None:
        row["structure_sha256s"]["oneA"] = "e" * 64

    def quality_record(row: dict[str, object]) -> None:
        row["quality_records"]["oneA"]["eligible"] = False

    def quality_hash(row: dict[str, object]) -> None:
        row["quality_record_sha256s"]["oneA"] = "e" * 64

    def bad_policy(row: dict[str, object]) -> None:
        row["policy"] = "other"

    def bad_runtime_hash(row: dict[str, object]) -> None:
        row["runtime_manifest_sha256"] = "bad"

    def bad_binary_hash(row: dict[str, object]) -> None:
        row["binary_sha256"] = "bad"

    def overlap(row: dict[str, object]) -> None:
        row["ineligible_ids"] = ["oneA", "twoB"]

    def omission(row: dict[str, object]) -> None:
        row["ineligible_ids"] = []

    def reorder(row: dict[str, object]) -> None:
        row["dataset_ids"] = ["twoB", "oneA"]

    for mutate in (
        duplicate_id,
        structure_hash,
        quality_record,
        quality_hash,
        bad_policy,
        bad_runtime_hash,
        bad_binary_hash,
        overlap,
        omission,
        reorder,
    ):
        payload = copy.deepcopy(valid)
        mutate(payload)
        with pytest.raises(ValueError):
            validate_eligibility_manifest(
                payload, expected_dataset_ids=expected_dataset_ids
            )


def test_build_rejects_report_chain_mismatch_and_changed_input(tmp_path: Path) -> None:
    chain_root = tmp_path / "chains"
    chain_root.mkdir()
    structure = chain_root / "oneA.pdb"
    structure.write_bytes(b"structure\n")
    identity = DatasetIdentity("1aaa", "oneA", "A")
    common = {
        "identities": [identity],
        "chain_root": chain_root,
        "dataset": "synthetic",
        "dataset_sha256": "d" * 64,
        "runtime": _synthetic_runtime(),
        "runtime_manifest_sha256": "c" * 64,
        "binary_sha256": "b" * 64,
        "jobs": 1,
    }
    with pytest.raises(ValueError):
        build_eligibility_manifest(
            **common,
            inspector=lambda _identity, _path: _report(chain_id="B", eligible=True),
        )

    def mutate_input(_identity: DatasetIdentity, path: Path) -> bytes:
        path.write_bytes(b"changed\n")
        return _report(chain_id="A", eligible=True)

    with pytest.raises(ValueError):
        build_eligibility_manifest(**common, inspector=mutate_input)


@pytest.mark.parametrize(("stderr", "returncode"), [(b"warning\n", 0), (b"", 7)])
def test_inspector_rejects_stderr_or_nonzero_exit(
    tmp_path: Path, stderr: bytes, returncode: int
) -> None:
    binary = tmp_path / "sword2"
    stderr_command = "printf 'warning\\n' >&2" if stderr else ":"
    binary.write_text(
        "#!/bin/sh\n"
        "printf '%s' '{\"not\":\"a quality record\"}\\n'\n"
        f"{stderr_command}\n"
        f"exit {returncode}\n",
        encoding="utf-8",
    )
    binary.chmod(0o755)
    structure = tmp_path / "oneA.pdb"
    structure.write_bytes(b"structure\n")

    with pytest.raises(ValueError):
        inspect_quality_record(
            binary=binary,
            identity=DatasetIdentity("1aaa", "oneA", "A"),
            structure_path=structure,
        )


def test_manifest_write_is_canonical_absent_only_and_refuses_symlink(tmp_path: Path) -> None:
    manifest = _manifest_fixture(tmp_path)
    target = tmp_path / "eligibility.json"
    digest = write_eligibility_manifest(target, manifest)
    data = target.read_bytes()
    assert data == canonical_json_bytes(manifest)
    assert digest == hashlib.sha256(data).hexdigest()
    with pytest.raises(FileExistsError):
        write_eligibility_manifest(target, manifest)

    symlink = tmp_path / "symlink.json"
    symlink.symlink_to(target)
    with pytest.raises(FileExistsError):
        write_eligibility_manifest(symlink, manifest)


def _authority_fixture(tmp_path: Path) -> dict[str, Path]:
    repo = tmp_path / "repo"
    (repo / "benchmark").mkdir(parents=True)
    cache = tmp_path / "cache"
    chains = cache / "chains"
    chains.mkdir(parents=True)
    (chains / "oneA.pdb").write_bytes(b"first structure\n")
    (chains / "twoB.pdb").write_bytes(b"second structure\n")
    dataset = tmp_path / "CATH-663.csv"
    dataset.write_bytes(
        b"1aaa,oneA,A,DO_NOT_PARSE_LABELS\n2bbb,twoB,B,DO_NOT_PARSE_LABELS\n"
    )
    runtime = tmp_path / "runtime.json"
    runtime.write_bytes(canonical_json_bytes(_synthetic_runtime()))
    binary = tmp_path / "sword2"
    binary.write_bytes(b"synthetic binary\n")
    binary.chmod(0o755)
    return {
        "repo": repo,
        "cache": cache,
        "dataset": dataset,
        "runtime": runtime,
        "binary": binary,
    }


def test_create_and_verify_rederive_every_authority_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    paths = _authority_fixture(tmp_path)
    reports = {
        "oneA": _report(chain_id="A", eligible=True),
        "twoB": _report(chain_id="B", eligible=False),
    }
    binary_hash = hashlib.sha256(paths["binary"].read_bytes()).hexdigest()
    runtime_payload = {
        "runtime_source_git_commit": "a" * 40,
        "binary_sha256": binary_hash,
    }
    monkeypatch.setattr(
        eligibility_module,
        "verify_runtime_freeze",
        lambda *_args, **_kwargs: runtime_payload,
    )
    monkeypatch.setattr(
        eligibility_module,
        "inspect_quality_record",
        lambda *, identity, **_kwargs: reports[identity.entry_id],
    )

    manifest = create_eligibility_manifest(
        dataset_metadata=paths["dataset"],
        cache_dir=paths["cache"],
        runtime_manifest=paths["runtime"],
        binary=paths["binary"],
        repo_root=paths["repo"],
        jobs=2,
    )
    assert manifest["dataset"] == "cath663"
    target = tmp_path / "eligibility.json"
    write_eligibility_manifest(target, manifest)
    verified = verify_eligibility_manifest(
        target,
        dataset_metadata=paths["dataset"],
        cache_dir=paths["cache"],
        runtime_manifest=paths["runtime"],
        binary=paths["binary"],
        repo_root=paths["repo"],
    )
    assert verified == manifest

    paths["cache"].joinpath("chains/oneA.pdb").write_bytes(b"tampered\n")
    with pytest.raises(ValueError):
        verify_eligibility_manifest(
            target,
            dataset_metadata=paths["dataset"],
            cache_dir=paths["cache"],
            runtime_manifest=paths["runtime"],
            binary=paths["binary"],
            repo_root=paths["repo"],
        )


def test_cli_parses_exact_create_and_verify_interfaces(tmp_path: Path) -> None:
    create = eligibility_cli.parse_args(
        [
            "create",
            "--dataset-metadata",
            "dataset.csv",
            "--cache-dir",
            "cache",
            "--runtime-manifest",
            "runtime.json",
            "--binary",
            "sword2",
            "--repo-root",
            ".",
            "--jobs",
            "32",
            "--out",
            "eligibility.json",
        ]
    )
    assert create.command == "create"
    assert create.jobs == 32

    verify = eligibility_cli.parse_args(
        [
            "verify",
            "--manifest",
            "eligibility.json",
            "--dataset-metadata",
            "dataset.csv",
            "--cache-dir",
            "cache",
            "--runtime-manifest",
            "runtime.json",
            "--binary",
            "sword2",
            "--repo-root",
            ".",
        ]
    )
    assert verify.command == "verify"
    assert not hasattr(verify, "jobs")


def test_cli_create_rejects_existing_output_and_input_alias(
    tmp_path: Path,
) -> None:
    paths = _authority_fixture(tmp_path)
    args = eligibility_cli.parse_args(
        [
            "create",
            "--dataset-metadata",
            os.fspath(paths["dataset"]),
            "--cache-dir",
            os.fspath(paths["cache"]),
            "--runtime-manifest",
            os.fspath(paths["runtime"]),
            "--binary",
            os.fspath(paths["binary"]),
            "--repo-root",
            os.fspath(paths["repo"]),
            "--jobs",
            "2",
            "--out",
            os.fspath(paths["runtime"]),
        ]
    )
    with pytest.raises((FileExistsError, ValueError)):
        eligibility_cli._validate_create_paths(args)


def test_cli_failure_is_concise_stderr_exit_two_and_no_stdout(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(
        eligibility_cli,
        "_validate_create_paths",
        lambda _args: (_ for _ in ()).throw(ValueError("fixture")),
    )
    exit_code = eligibility_cli.main(
        [
            "create",
            "--dataset-metadata",
            "dataset.csv",
            "--cache-dir",
            "cache",
            "--runtime-manifest",
            "runtime.json",
            "--binary",
            "sword2",
            "--repo-root",
            ".",
            "--jobs",
            "1",
            "--out",
            "out.json",
        ]
    )
    captured = capsys.readouterr()
    assert exit_code == 2
    assert captured.out == ""
    assert captured.err == "invalid eligibility freeze: fixture\n"
