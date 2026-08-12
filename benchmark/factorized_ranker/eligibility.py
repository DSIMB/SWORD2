"""Canonical input-only structural eligibility for the frozen factorized ranker."""

from __future__ import annotations

import hashlib
import json
import math
import os
import stat
import subprocess
import unicodedata
from collections import Counter
from collections.abc import Callable, Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from benchmark.factorized_ranker.runtime_freeze import (
    canonical_id_set_hash,
    canonical_json_bytes,
    verify_runtime_freeze,
)


ELIGIBILITY_SCHEMA_VERSION = 1
QUALITY_REPORT_SCHEMA_VERSION = 1
POLICY = "strict_complete_backbone_v1"
REQUIRED_ATOM_ORDER = ("N", "CA", "C", "O")

QUALITY_REPORT_KEYS = {
    "candidate_residue_count",
    "chain_id",
    "complete_backbone_residue_count",
    "eligible",
    "incomplete_residues",
    "policy",
    "reason_code",
    "schema_version",
    "structural_coverage",
}
INCOMPLETE_RESIDUE_KEYS = {
    "author_residue_number",
    "chain_id",
    "missing_atoms",
}
ELIGIBILITY_MANIFEST_KEYS = {
    "binary_sha256",
    "dataset",
    "dataset_id_count",
    "dataset_id_set_sha256",
    "dataset_ids",
    "dataset_sha256",
    "eligible_count",
    "eligible_id_set_sha256",
    "eligible_ids",
    "ineligibility_reason_counts",
    "ineligible_count",
    "ineligible_id_set_sha256",
    "ineligible_ids",
    "policy",
    "quality_record_sha256s",
    "quality_record_tree_sha256",
    "quality_records",
    "runtime_manifest_sha256",
    "runtime_source_git_commit",
    "schema_version",
    "structure_sha256s",
    "structure_tree_sha256",
}


@dataclass(frozen=True, slots=True)
class DatasetIdentity:
    pdb_id: str
    entry_id: str
    chain_id: str


def _reject_duplicate_keys(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value}")


def _canonical_text(value: object, description: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value.strip() != value
        or unicodedata.normalize("NFC", value) != value
        or "\0" in value
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise ValueError(f"{description} is not canonical text")
    try:
        value.encode("ascii")
    except UnicodeEncodeError as error:
        raise ValueError(f"{description} must be ASCII") from error
    return value


def _canonical_chain_id(value: object) -> str:
    chain_id = _canonical_text(value, "chain ID")
    if len(chain_id) != 1:
        raise ValueError("chain ID must contain exactly one ASCII character")
    return chain_id


def _validate_identity(identity: DatasetIdentity) -> DatasetIdentity:
    if not isinstance(identity, DatasetIdentity):
        raise ValueError("dataset identity has the wrong type")
    pdb_id = _canonical_text(identity.pdb_id, "PDB ID")
    entry_id = _canonical_text(identity.entry_id, "entry ID")
    if Path(entry_id).name != entry_id or entry_id in {".", ".."} or "\\" in entry_id:
        raise ValueError("entry ID is not a safe file identity")
    chain_id = _canonical_chain_id(identity.chain_id)
    return DatasetIdentity(pdb_id, entry_id, chain_id)


def _require_plain_int(
    value: object, description: str, minimum: int, maximum: int
) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{description} is outside its integer contract")
    return value


def _require_hash(value: object, description: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{description} must be a lowercase SHA-256")
    return value


def _require_git_commit(value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 40
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError("runtime source commit must be a full lowercase Git SHA-1")
    return value


def stable_file_bytes(path: Path) -> bytes:
    path = Path(path)
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError(f"input is not a regular nonsymlink file: {path}")
    data = path.read_bytes()
    after = path.lstat()
    before_identity = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
    )
    if before_identity != after_identity or len(data) != before.st_size:
        raise ValueError(f"input changed while reading: {path}")
    return data


def _stable_file_hash_and_identity(
    path: Path,
) -> tuple[str, tuple[int, int, int, int, int]]:
    path = Path(path)
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError(f"input is not a regular nonsymlink file: {path}")
    digest = hashlib.sha256()
    byte_count = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
            byte_count += len(chunk)
    after = path.lstat()
    before_identity = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if before_identity != after_identity or byte_count != before.st_size:
        raise ValueError(f"input changed while hashing: {path}")
    return digest.hexdigest(), before_identity


def _mapping_tree_hash(mapping: Mapping[str, str]) -> str:
    digest = hashlib.sha256()
    for name, value in sorted(mapping.items()):
        name_bytes = name.encode("utf-8")
        value_bytes = value.encode("ascii")
        digest.update(len(name_bytes).to_bytes(8, "big"))
        digest.update(name_bytes)
        digest.update(len(value_bytes).to_bytes(8, "big"))
        digest.update(value_bytes)
    return digest.hexdigest()


def read_dataset_identity_prefixes(
    path: Path,
) -> tuple[list[DatasetIdentity], str]:
    data = stable_file_bytes(path)
    identities: list[DatasetIdentity] = []
    seen_entry_ids: set[str] = set()
    for line in data.splitlines():
        if not line or line.startswith(b"#"):
            continue
        prefix = line.split(b",", 3)
        if len(prefix) != 4:
            raise ValueError("dataset identity row is truncated")
        try:
            pdb_id, entry_id, chain_id = (
                field.decode("ascii") for field in prefix[:3]
            )
        except UnicodeDecodeError as error:
            raise ValueError("dataset identity prefix must be ASCII") from error
        identity = _validate_identity(DatasetIdentity(pdb_id, entry_id, chain_id))
        if identity.entry_id in seen_entry_ids:
            raise ValueError("dataset identity prefixes contain duplicate entry IDs")
        seen_entry_ids.add(identity.entry_id)
        identities.append(identity)
    if not identities:
        raise ValueError("dataset contains no identities")
    return identities, hashlib.sha256(data).hexdigest()


def _validate_quality_payload(
    payload: object, *, expected_chain_id: str | None = None
) -> dict[str, object]:
    if not isinstance(payload, dict) or set(payload) != QUALITY_REPORT_KEYS:
        raise ValueError("quality record top-level schema mismatch")
    if (
        _require_plain_int(payload["schema_version"], "quality schema version", 1, 1)
        != QUALITY_REPORT_SCHEMA_VERSION
    ):
        raise ValueError("unsupported quality record schema")
    if payload["policy"] != POLICY:
        raise ValueError("quality record policy mismatch")
    chain_id = _canonical_chain_id(payload["chain_id"])
    if expected_chain_id is not None and chain_id != _canonical_chain_id(expected_chain_id):
        raise ValueError("quality record chain does not match dataset identity")

    candidate_count = _require_plain_int(
        payload["candidate_residue_count"], "candidate residue count", 0, 2**64 - 1
    )
    complete_count = _require_plain_int(
        payload["complete_backbone_residue_count"],
        "complete-backbone residue count",
        0,
        2**64 - 1,
    )
    eligible = payload["eligible"]
    if type(eligible) is not bool:
        raise ValueError("quality eligibility must be Boolean")
    coverage = payload["structural_coverage"]
    if type(coverage) is not float or not math.isfinite(coverage):
        raise ValueError("structural coverage must be a finite JSON float")

    raw_incomplete = payload["incomplete_residues"]
    if not isinstance(raw_incomplete, list):
        raise ValueError("incomplete residues must be a list")
    seen_residues: set[tuple[str, int]] = set()
    for raw_residue in raw_incomplete:
        if not isinstance(raw_residue, dict) or set(raw_residue) != INCOMPLETE_RESIDUE_KEYS:
            raise ValueError("incomplete residue schema mismatch")
        residue_chain = _canonical_chain_id(raw_residue["chain_id"])
        if residue_chain != chain_id:
            raise ValueError("incomplete residue chain does not match its quality record")
        residue_number = _require_plain_int(
            raw_residue["author_residue_number"],
            "author residue number",
            -(2**31),
            2**31 - 1,
        )
        residue_identity = (residue_chain, residue_number)
        if residue_identity in seen_residues:
            raise ValueError("quality record contains duplicate incomplete residues")
        seen_residues.add(residue_identity)
        missing_atoms = raw_residue["missing_atoms"]
        if not isinstance(missing_atoms, list) or not missing_atoms:
            raise ValueError("missing atoms must be a nonempty list")
        if any(type(atom) is not str or atom not in REQUIRED_ATOM_ORDER for atom in missing_atoms):
            raise ValueError("missing atoms contain an unknown atom name")
        positions = [REQUIRED_ATOM_ORDER.index(atom) for atom in missing_atoms]
        if positions != sorted(set(positions)):
            raise ValueError("missing atoms are not a canonical backbone subsequence")

    if complete_count > candidate_count:
        raise ValueError("complete-backbone count exceeds candidate count")
    if candidate_count - complete_count != len(raw_incomplete):
        raise ValueError("quality residue counts disagree with incomplete rows")
    expected_coverage = 0.0 if candidate_count == 0 else complete_count / candidate_count
    if coverage != expected_coverage:
        raise ValueError("structural coverage disagrees with residue counts")
    expected_eligible = candidate_count > 0 and not raw_incomplete
    if eligible is not expected_eligible:
        raise ValueError("quality eligibility disagrees with residue evidence")
    expected_reason = "empty_candidate_population" if candidate_count == 0 else None
    if payload["reason_code"] != expected_reason:
        raise ValueError("quality reason code disagrees with residue evidence")
    return payload


def load_quality_record_bytes(
    data: bytes, *, expected_chain_id: str | None = None
) -> dict[str, object]:
    if type(data) is not bytes:
        raise ValueError("quality record must be bytes")
    try:
        payload = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("quality record is not valid canonical JSON") from error
    if canonical_json_bytes(payload) != data:
        raise ValueError("quality record is not canonically serialized")
    return _validate_quality_payload(payload, expected_chain_id=expected_chain_id)


def inspect_quality_record(
    *, binary: Path, identity: DatasetIdentity, structure_path: Path
) -> bytes:
    identity = _validate_identity(identity)
    binary_info = Path(binary).lstat()
    if stat.S_ISLNK(binary_info.st_mode) or not stat.S_ISREG(binary_info.st_mode):
        raise ValueError("inspector binary is not a regular nonsymlink file")
    binary_path = Path(binary).resolve(strict=True)
    if not os.access(binary_path, os.X_OK):
        raise ValueError("inspector binary is not executable")
    structure = Path(structure_path).resolve(strict=True)
    completed = subprocess.run(
        [
            os.fspath(binary_path),
            "inspect-factorized-eligibility",
            "--input",
            os.fspath(structure),
            "--chain",
            identity.chain_id,
            "--nmr-model",
            "1",
        ],
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise ValueError("factorized eligibility inspector exited nonzero")
    if completed.stderr:
        raise ValueError("factorized eligibility inspector wrote to stderr")
    load_quality_record_bytes(completed.stdout, expected_chain_id=identity.chain_id)
    return completed.stdout


def _ineligibility_reason(record: Mapping[str, object]) -> str:
    reason = record["reason_code"]
    return str(reason) if reason is not None else "incomplete_backbone"


def build_eligibility_manifest(
    *,
    identities: Sequence[DatasetIdentity],
    chain_root: Path,
    dataset: str,
    dataset_sha256: str,
    runtime: Mapping[str, object],
    runtime_manifest_sha256: str,
    binary_sha256: str,
    inspector: Callable[[DatasetIdentity, Path], bytes],
    jobs: int,
) -> dict[str, object]:
    dataset = _canonical_text(dataset, "dataset name")
    dataset_sha256 = _require_hash(dataset_sha256, "dataset hash")
    runtime_manifest_sha256 = _require_hash(
        runtime_manifest_sha256, "runtime manifest hash"
    )
    binary_sha256 = _require_hash(binary_sha256, "binary hash")
    runtime_binary_hash = _require_hash(runtime.get("binary_sha256"), "runtime binary hash")
    if binary_sha256 != runtime_binary_hash:
        raise ValueError("binary hash disagrees with runtime manifest")
    runtime_source_git_commit = _require_git_commit(
        runtime.get("runtime_source_git_commit")
    )
    _require_plain_int(jobs, "eligibility worker count", 1, 2**31 - 1)
    if not callable(inspector):
        raise ValueError("eligibility inspector must be callable")

    normalized_identities = [_validate_identity(identity) for identity in identities]
    dataset_ids = [identity.entry_id for identity in normalized_identities]
    if not dataset_ids:
        raise ValueError("eligibility population is empty")
    if len(dataset_ids) != len(set(dataset_ids)):
        raise ValueError("eligibility population contains duplicate entry IDs")

    root_info = Path(chain_root).lstat()
    if stat.S_ISLNK(root_info.st_mode) or not stat.S_ISDIR(root_info.st_mode):
        raise ValueError("chain root must be a regular nonsymlink directory")
    root = Path(chain_root).resolve(strict=True)
    structure_paths = {
        identity.entry_id: root / f"{identity.entry_id}.pdb"
        for identity in normalized_identities
    }
    initial_structure_evidence = {
        entry_id: _stable_file_hash_and_identity(path)
        for entry_id, path in structure_paths.items()
    }
    structure_sha256s = {
        entry_id: evidence[0]
        for entry_id, evidence in initial_structure_evidence.items()
    }

    def inspect_one(identity: DatasetIdentity) -> tuple[bytes, dict[str, object]]:
        raw = inspector(identity, structure_paths[identity.entry_id])
        record = load_quality_record_bytes(raw, expected_chain_id=identity.chain_id)
        return raw, record

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [executor.submit(inspect_one, identity) for identity in normalized_identities]
        inspected = [future.result() for future in futures]

    for entry_id, path in structure_paths.items():
        if _stable_file_hash_and_identity(path) != initial_structure_evidence[entry_id]:
            raise ValueError("structure input changed during eligibility inspection")

    quality_records = {
        identity.entry_id: record
        for identity, (_raw, record) in zip(normalized_identities, inspected, strict=True)
    }
    quality_record_sha256s = {
        identity.entry_id: hashlib.sha256(raw).hexdigest()
        for identity, (raw, _record) in zip(normalized_identities, inspected, strict=True)
    }
    eligible_ids = [
        entry_id for entry_id in dataset_ids if quality_records[entry_id]["eligible"] is True
    ]
    ineligible_ids = [
        entry_id for entry_id in dataset_ids if quality_records[entry_id]["eligible"] is False
    ]
    reason_counts = Counter(
        _ineligibility_reason(quality_records[entry_id]) for entry_id in ineligible_ids
    )

    manifest: dict[str, object] = {
        "binary_sha256": binary_sha256,
        "dataset": dataset,
        "dataset_id_count": len(dataset_ids),
        "dataset_id_set_sha256": canonical_id_set_hash(dataset_ids),
        "dataset_ids": dataset_ids,
        "dataset_sha256": dataset_sha256,
        "eligible_count": len(eligible_ids),
        "eligible_id_set_sha256": canonical_id_set_hash(eligible_ids),
        "eligible_ids": eligible_ids,
        "ineligibility_reason_counts": dict(sorted(reason_counts.items())),
        "ineligible_count": len(ineligible_ids),
        "ineligible_id_set_sha256": canonical_id_set_hash(ineligible_ids),
        "ineligible_ids": ineligible_ids,
        "policy": POLICY,
        "quality_record_sha256s": quality_record_sha256s,
        "quality_record_tree_sha256": _mapping_tree_hash(quality_record_sha256s),
        "quality_records": quality_records,
        "runtime_manifest_sha256": runtime_manifest_sha256,
        "runtime_source_git_commit": runtime_source_git_commit,
        "schema_version": ELIGIBILITY_SCHEMA_VERSION,
        "structure_sha256s": structure_sha256s,
        "structure_tree_sha256": _mapping_tree_hash(structure_sha256s),
    }
    validate_eligibility_manifest(manifest, expected_dataset_ids=dataset_ids)
    return manifest


def _require_id_list(value: object, description: str) -> list[str]:
    if not isinstance(value, list):
        raise ValueError(f"{description} must be a list")
    ids = [_canonical_text(item, description) for item in value]
    if len(ids) != len(set(ids)):
        raise ValueError(f"{description} contains duplicates")
    return ids


def _require_hash_mapping(
    value: object, expected_ids: Iterable[str], description: str
) -> dict[str, str]:
    if not isinstance(value, dict) or set(value) != set(expected_ids):
        raise ValueError(f"{description} identity mapping mismatch")
    return {
        _canonical_text(key, description): _require_hash(item, description)
        for key, item in value.items()
    }


def validate_eligibility_manifest(
    payload: Mapping[str, object], *, expected_dataset_ids: Sequence[str] | None = None
) -> None:
    if not isinstance(payload, dict) or set(payload) != ELIGIBILITY_MANIFEST_KEYS:
        raise ValueError("eligibility manifest top-level schema mismatch")
    if (
        _require_plain_int(payload["schema_version"], "eligibility schema version", 1, 1)
        != ELIGIBILITY_SCHEMA_VERSION
    ):
        raise ValueError("unsupported eligibility manifest schema")
    if payload["policy"] != POLICY:
        raise ValueError("eligibility manifest policy mismatch")
    _canonical_text(payload["dataset"], "dataset name")
    for field in (
        "binary_sha256",
        "dataset_sha256",
        "dataset_id_set_sha256",
        "eligible_id_set_sha256",
        "ineligible_id_set_sha256",
        "quality_record_tree_sha256",
        "runtime_manifest_sha256",
        "structure_tree_sha256",
    ):
        _require_hash(payload[field], field)
    _require_git_commit(payload["runtime_source_git_commit"])

    dataset_ids = _require_id_list(payload["dataset_ids"], "dataset IDs")
    if expected_dataset_ids is not None and dataset_ids != list(expected_dataset_ids):
        raise ValueError("eligibility dataset order disagrees with its authority")
    dataset_id_count = _require_plain_int(
        payload["dataset_id_count"], "dataset ID count", 1, 2**64 - 1
    )
    if dataset_id_count != len(dataset_ids):
        raise ValueError("eligibility dataset count mismatch")
    if payload["dataset_id_set_sha256"] != canonical_id_set_hash(dataset_ids):
        raise ValueError("eligibility dataset ID-set hash mismatch")

    eligible_ids = _require_id_list(payload["eligible_ids"], "eligible IDs")
    ineligible_ids = _require_id_list(payload["ineligible_ids"], "ineligible IDs")
    eligible_count = _require_plain_int(
        payload["eligible_count"], "eligible count", 0, dataset_id_count
    )
    ineligible_count = _require_plain_int(
        payload["ineligible_count"], "ineligible count", 0, dataset_id_count
    )
    if eligible_count != len(eligible_ids) or ineligible_count != len(ineligible_ids):
        raise ValueError("eligibility subset counts mismatch")
    if set(eligible_ids) & set(ineligible_ids):
        raise ValueError("eligible and ineligible identity sets overlap")
    if set(eligible_ids) | set(ineligible_ids) != set(dataset_ids):
        raise ValueError("eligible and ineligible identities do not partition the dataset")
    if payload["eligible_id_set_sha256"] != canonical_id_set_hash(eligible_ids):
        raise ValueError("eligible ID-set hash mismatch")
    if payload["ineligible_id_set_sha256"] != canonical_id_set_hash(ineligible_ids):
        raise ValueError("ineligible ID-set hash mismatch")

    structure_hashes = _require_hash_mapping(
        payload["structure_sha256s"], dataset_ids, "structure hashes"
    )
    quality_hashes = _require_hash_mapping(
        payload["quality_record_sha256s"], dataset_ids, "quality record hashes"
    )
    if payload["structure_tree_sha256"] != _mapping_tree_hash(structure_hashes):
        raise ValueError("structure tree hash mismatch")
    if payload["quality_record_tree_sha256"] != _mapping_tree_hash(quality_hashes):
        raise ValueError("quality-record tree hash mismatch")

    raw_records = payload["quality_records"]
    if not isinstance(raw_records, dict) or set(raw_records) != set(dataset_ids):
        raise ValueError("quality-record identity mapping mismatch")
    records: dict[str, dict[str, object]] = {}
    for entry_id in dataset_ids:
        record_bytes = canonical_json_bytes(raw_records[entry_id])
        record = load_quality_record_bytes(record_bytes)
        if hashlib.sha256(record_bytes).hexdigest() != quality_hashes[entry_id]:
            raise ValueError("embedded quality record hash mismatch")
        records[entry_id] = record
    expected_eligible = [entry_id for entry_id in dataset_ids if records[entry_id]["eligible"]]
    expected_ineligible = [entry_id for entry_id in dataset_ids if not records[entry_id]["eligible"]]
    if eligible_ids != expected_eligible or ineligible_ids != expected_ineligible:
        raise ValueError("eligibility subsets disagree with quality records")

    raw_reason_counts = payload["ineligibility_reason_counts"]
    if not isinstance(raw_reason_counts, dict):
        raise ValueError("ineligibility reason counts must be a mapping")
    reason_counts: dict[str, int] = {}
    for reason, count in raw_reason_counts.items():
        reason = _canonical_text(reason, "ineligibility reason")
        reason_counts[reason] = _require_plain_int(
            count, "ineligibility reason count", 1, dataset_id_count
        )
    expected_reason_counts = dict(
        sorted(
            Counter(
                _ineligibility_reason(records[entry_id])
                for entry_id in ineligible_ids
            ).items()
        )
    )
    if reason_counts != expected_reason_counts:
        raise ValueError("ineligibility reason counts disagree with quality records")


def write_eligibility_manifest(path: Path, payload: Mapping[str, object]) -> str:
    validate_eligibility_manifest(payload)
    target = Path(path)
    if target.exists() or target.is_symlink():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    data = canonical_json_bytes(dict(payload))
    temporary = target.with_name(f".{target.name}.tmp-{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FileExistsError(temporary)
    try:
        with temporary.open("xb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, target)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    temporary.unlink()
    if stable_file_bytes(target) != data:
        raise OSError("installed eligibility manifest bytes differ")
    return hashlib.sha256(data).hexdigest()


def _dataset_name(path: Path) -> str:
    stem = _canonical_text(Path(path).stem, "dataset file stem").lower()
    normalized = "".join(character for character in stem if character.isalnum())
    if not normalized or any(ord(character) > 127 for character in normalized):
        raise ValueError("dataset file name does not define a canonical dataset identity")
    return normalized


def load_eligibility_manifest(path: Path) -> dict[str, object]:
    data = stable_file_bytes(path)
    try:
        payload = json.loads(
            data.decode("utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("eligibility manifest is not valid canonical JSON") from error
    if not isinstance(payload, dict) or canonical_json_bytes(payload) != data:
        raise ValueError("eligibility manifest is not canonically serialized")
    validate_eligibility_manifest(payload)
    return payload


def create_eligibility_manifest(
    *,
    dataset_metadata: Path,
    cache_dir: Path,
    runtime_manifest: Path,
    binary: Path,
    repo_root: Path,
    jobs: int,
) -> dict[str, object]:
    dataset_path = Path(dataset_metadata)
    runtime_path = Path(runtime_manifest)
    binary_path = Path(binary)
    dataset_bytes = stable_file_bytes(dataset_path)
    runtime_bytes = stable_file_bytes(runtime_path)
    binary_bytes = stable_file_bytes(binary_path)
    runtime = verify_runtime_freeze(
        runtime_path,
        binary=binary_path,
        repo_root=repo_root,
    )
    binary_sha256 = hashlib.sha256(binary_bytes).hexdigest()
    if runtime.get("binary_sha256") != binary_sha256:
        raise ValueError("runtime manifest does not bind the inspector binary")
    identities, dataset_sha256 = read_dataset_identity_prefixes(dataset_path)
    if hashlib.sha256(dataset_bytes).hexdigest() != dataset_sha256:
        raise ValueError("dataset metadata changed while parsing identities")

    cache_info = Path(cache_dir).lstat()
    if stat.S_ISLNK(cache_info.st_mode) or not stat.S_ISDIR(cache_info.st_mode):
        raise ValueError("cache directory must be a regular nonsymlink directory")
    chain_root = Path(cache_dir).resolve(strict=True) / "chains"

    manifest = build_eligibility_manifest(
        identities=identities,
        chain_root=chain_root,
        dataset=_dataset_name(dataset_path),
        dataset_sha256=dataset_sha256,
        runtime=runtime,
        runtime_manifest_sha256=hashlib.sha256(runtime_bytes).hexdigest(),
        binary_sha256=binary_sha256,
        inspector=lambda identity, structure_path: inspect_quality_record(
            binary=binary_path,
            identity=identity,
            structure_path=structure_path,
        ),
        jobs=jobs,
    )

    if stable_file_bytes(dataset_path) != dataset_bytes:
        raise ValueError("dataset metadata changed during eligibility inspection")
    if stable_file_bytes(runtime_path) != runtime_bytes:
        raise ValueError("runtime manifest changed during eligibility inspection")
    if stable_file_bytes(binary_path) != binary_bytes:
        raise ValueError("inspector binary changed during eligibility inspection")
    if (
        verify_runtime_freeze(
            runtime_path,
            binary=binary_path,
            repo_root=repo_root,
        )
        != runtime
    ):
        raise ValueError("runtime authority changed during eligibility inspection")
    return manifest


def verify_eligibility_manifest(
    path: Path,
    *,
    dataset_metadata: Path,
    cache_dir: Path,
    runtime_manifest: Path,
    binary: Path,
    repo_root: Path,
) -> dict[str, object]:
    manifest_bytes = stable_file_bytes(path)
    payload = load_eligibility_manifest(path)
    identities, _dataset_sha256 = read_dataset_identity_prefixes(dataset_metadata)
    expected_ids = [identity.entry_id for identity in identities]
    validate_eligibility_manifest(payload, expected_dataset_ids=expected_ids)
    expected = create_eligibility_manifest(
        dataset_metadata=dataset_metadata,
        cache_dir=cache_dir,
        runtime_manifest=runtime_manifest,
        binary=binary,
        repo_root=repo_root,
        jobs=min(32, len(expected_ids)),
    )
    if payload != expected:
        raise ValueError("eligibility manifest no longer matches its authority inputs")
    if stable_file_bytes(path) != manifest_bytes:
        raise ValueError("eligibility manifest changed while verifying")
    return payload


__all__ = [
    "DatasetIdentity",
    "ELIGIBILITY_MANIFEST_KEYS",
    "ELIGIBILITY_SCHEMA_VERSION",
    "POLICY",
    "QUALITY_REPORT_KEYS",
    "REQUIRED_ATOM_ORDER",
    "build_eligibility_manifest",
    "create_eligibility_manifest",
    "inspect_quality_record",
    "load_eligibility_manifest",
    "load_quality_record_bytes",
    "read_dataset_identity_prefixes",
    "stable_file_bytes",
    "validate_eligibility_manifest",
    "verify_eligibility_manifest",
    "write_eligibility_manifest",
]
