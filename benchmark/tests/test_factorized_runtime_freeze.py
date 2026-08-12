from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from pathlib import Path

import pytest

import benchmark.factorized_ranker.runtime_freeze as runtime_freeze_module
from benchmark.factorized_ranker.runtime_freeze import (
    BUILD_COMMAND,
    CACHE_CONTRACT,
    EVIDENCE_TOOL_PATHS_V1,
    EVIDENCE_TOOL_PATHS_V2,
    MODEL_ARTIFACT_FIELDS,
    RUNTIME_FREEZE_SCHEMA_VERSION,
    canonical_id_set_hash,
    canonical_json_bytes,
    hash_directory_tree,
    hash_file_or_tree,
    verify_runtime_freeze,
    write_runtime_freeze,
)


def _tree_hash(mapping: dict[str, str]) -> str:
    value = hashlib.sha256()
    for name, item in sorted(mapping.items()):
        name_bytes = name.encode()
        item_bytes = item.encode()
        value.update(len(name_bytes).to_bytes(8, "big"))
        value.update(name_bytes)
        value.update(len(item_bytes).to_bytes(8, "big"))
        value.update(item_bytes)
    return value.hexdigest()


def _runtime_payload(schema_version: int = 1) -> dict[str, object]:
    digest = "a" * 64
    runtime_inputs = {"Cargo.lock": digest, "Cargo.toml": "b" * 64}
    evidence_paths = (
        EVIDENCE_TOOL_PATHS_V1 if schema_version == 1 else EVIDENCE_TOOL_PATHS_V2
    )
    evidence_inputs = {
        name: hashlib.sha256(name.encode()).hexdigest()
        for name in evidence_paths
    }

    population_hash = "d" * 64
    return {
        "schema_version": schema_version,
        "model_manifest_sha256": "e" * 64,
        "model_source_git_commit": "1" * 40,
        "model_artifact_sha256s": {
            name: hashlib.sha256(name.encode()).hexdigest()
            for name in MODEL_ARTIFACT_FIELDS
        },
        "development_population": {
            "accepted_count": 2,
            "rejected_count": 1,
            "accepted_id_set_sha256": population_hash,
            "fold_id_set_sha256": population_hash,
            "oof_id_set_sha256": population_hash,
            "rejected_id_set_sha256": "f" * 64,
            "dataset_sha256": "0" * 64,
            "corpus_manifest_sha256": "1" * 64,
            "corpus_table_sha256s": {
                "chains": "2" * 64,
                "counts": "3" * 64,
                "candidates": "4" * 64,
                "rejections": "5" * 64,
            },
            "fold_manifest_sha256": "6" * 64,
            "oof_predictions_sha256": "7" * 64,
            "n_folds": 5,
            "seed": 37,
        },
        "runtime_source_git_commit": "2" * 40,
        "runtime_input_sha256s": runtime_inputs,
        "runtime_input_tree_sha256": _tree_hash(runtime_inputs),
        "evidence_tool_sha256s": evidence_inputs,
        "evidence_tool_tree_sha256": _tree_hash(evidence_inputs),
        "cargo_lock_sha256": digest,
        "rustc_version_verbose": "rustc test\nhost: x86_64-unknown-linux-gnu",
        "cargo_version": "cargo test",
        "target_triple": "x86_64-unknown-linux-gnu",
        "build_profile": "release",
        "build_command": BUILD_COMMAND,
        "build_environment": {
            "CARGO_INCREMENTAL": None,
            "RUSTFLAGS": None,
            "RUSTC_WRAPPER": None,
        },
        "binary_name": "sword2",
        "binary_size": 1,
        "binary_sha256": "8" * 64,
        "selector_cache_contract": CACHE_CONTRACT,
    }


def test_runtime_v2_binds_eligibility_tools_and_v1_keeps_its_original_path_set(
    tmp_path: Path,
):
    assert RUNTIME_FREEZE_SCHEMA_VERSION == 2
    assert set(EVIDENCE_TOOL_PATHS_V2) == {
        *EVIDENCE_TOOL_PATHS_V1,
        "benchmark/factorized_ranker/eligibility.py",
        "benchmark/freeze_factorized_eligibility.py",
    }

    created = _runtime_payload(schema_version=2)
    assert created["schema_version"] == 2
    assert set(created["evidence_tool_sha256s"]) == set(EVIDENCE_TOOL_PATHS_V2)
    write_runtime_freeze(tmp_path / "v2.json", created)

    legacy = _runtime_payload(schema_version=1)
    assert set(legacy["evidence_tool_sha256s"]) == set(EVIDENCE_TOOL_PATHS_V1)
    write_runtime_freeze(tmp_path / "v1.json", legacy)


def test_runtime_schema_and_new_tool_hash_tampering_are_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    legacy_relabelled = _runtime_payload(schema_version=1)
    legacy_relabelled["schema_version"] = 2
    with pytest.raises(ValueError, match="path set"):
        write_runtime_freeze(tmp_path / "legacy-relabelled.json", legacy_relabelled)

    expected = _runtime_payload(schema_version=2)
    repo = tmp_path / "repo"
    repo.mkdir()
    binary = tmp_path / "sword2"
    binary.write_bytes(b"binary\n")
    monkeypatch.setattr(runtime_freeze_module, "_require_ancestor", lambda *_args: None)
    monkeypatch.setattr(
        runtime_freeze_module, "_build_payload", lambda **_kwargs: expected
    )
    for tool in (
        "benchmark/factorized_ranker/eligibility.py",
        "benchmark/freeze_factorized_eligibility.py",
    ):
        tampered = _runtime_payload(schema_version=2)
        tampered["evidence_tool_sha256s"][tool] = "f" * 64
        tampered["evidence_tool_tree_sha256"] = _tree_hash(
            tampered["evidence_tool_sha256s"]
        )
        path = tmp_path / f"tampered-{Path(tool).name}.json"
        write_runtime_freeze(path, tampered)
        with pytest.raises(ValueError, match="no longer matches"):
            verify_runtime_freeze(path, binary=binary, repo_root=repo)


def test_canonical_runtime_helpers_are_root_independent_and_strict(tmp_path: Path):
    assert canonical_json_bytes({"z": 1, "a": ["é"]}) == b'{"a":["\xc3\xa9"],"z":1}\n'
    expected = hashlib.sha256(b'["a","b"]\n').hexdigest()
    assert canonical_id_set_hash(["b", "a"]) == expected

    for invalid in (["a", "a"], [""], [" a"], ["a\n"], [1]):
        with pytest.raises(ValueError):
            canonical_id_set_hash(invalid)  # type: ignore[arg-type]

    left = tmp_path / "left"
    right = tmp_path / "right"
    for root in (left, right):
        (root / "nested").mkdir(parents=True)
        (root / "a.txt").write_bytes(b"a\n")
        (root / "nested" / "b.bin").write_bytes(b"b")
    assert hash_directory_tree(left) == hash_directory_tree(right)

    (left / "link").symlink_to(left / "a.txt")
    with pytest.raises(ValueError, match="symlink"):
        hash_directory_tree(left)


def test_opaque_artifact_tree_hashes_symlink_target_without_following(tmp_path: Path):
    root = tmp_path / "artifact"
    root.mkdir()
    (root / "payload.txt").write_bytes(b"payload\n")
    (root / "relative-link").symlink_to("payload.txt")

    assert hash_file_or_tree(root) == {
        "kind": "tree",
        "byte_count": 19,
        "file_count": 2,
        "sha256": "f9e30f8564dee7116cb1e8e4c0b439ce80738670caa7aca2065dac0936d14971",
    }
    assert hash_file_or_tree(root / "relative-link") == {
        "kind": "file",
        "byte_count": 11,
        "file_count": 1,
        "sha256": "af3ff413ed80f409410ac36d2a0e22a3af2f0cde2db0f44594a292df9cf3fac6",
    }


def test_runtime_manifest_write_is_canonical_absent_only_and_cache_is_frozen(tmp_path: Path):
    assert CACHE_CONTRACT == {
        "default_promotion_compatible": False,
        "factorized_cache_hit_behavior": "whole_chain_legacy_fallback",
        "typed_context_reloaded": False,
    }
    path = tmp_path / "runtime.json"
    payload = _runtime_payload()
    digest = write_runtime_freeze(path, payload)
    assert path.read_bytes() == canonical_json_bytes(payload)
    assert digest == hashlib.sha256(path.read_bytes()).hexdigest()

    with pytest.raises(FileExistsError):
        write_runtime_freeze(path, payload)
    assert json.loads(path.read_bytes())["selector_cache_contract"] == CACHE_CONTRACT

    incompatible = deepcopy(payload)
    incompatible["selector_cache_contract"]["default_promotion_compatible"] = True
    invalid_path = tmp_path / "invalid.json"
    with pytest.raises(ValueError, match="cache contract"):
        write_runtime_freeze(invalid_path, incompatible)
    assert not invalid_path.exists()
