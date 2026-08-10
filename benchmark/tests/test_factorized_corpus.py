from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

import benchmark.build_factorized_corpus as build_factorized_corpus
from benchmark.build_training_table import ScoredChain
from benchmark.datasets import CathEntry
from benchmark.dump_candidate_corpus import (
    required_dump_fields,
    valid_resume_part,
    validate_dump_header,
    validate_dump_rows,
)
from benchmark.build_factorized_corpus import validate_raw_population
from benchmark.factorized_ranker.integrity import RejectionCode, RejectionRecord
from benchmark.factorized_ranker.corpus import (
    ACQUISITION_REJECTION_FIELDS,
    NORMALIZED_REJECTION_FIELDS,
    CorpusPaths,
    feature_schema_hash,
    manifest_hash,
    normalize_argv,
    write_corpus,
)
from benchmark.factorized_ranker.schema import (
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    GLOBAL_FEATURES,
)


def _raw_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for chain_id, count, delineation, source in (
        ("chainB", 2, "0-1 2-3", 9),
        ("chainA", 1, "0-3", 4),
        ("chainA", 2, "0-1 2-3", 7),
    ):
        row: dict[str, object] = {
            "chain_id": chain_id,
            "canonical_delineation": delineation,
            "source_index": source,
            "legacy_distance": 0.125 + source / 1000,
            "n_true_domains": 2,
            "n_pred_domains": count,
            "ndo": 0.1,
            "iou": 0.2,
            "boundary_f1_10": 0.3,
            "matched_dice": 0.4,
            "d_count_acc": 1 if count == 2 else 0,
            "S": 0.5,
            "is_oracle_s": 1 if count == 2 else 0,
        }
        row.update({name: index + 0.25 for index, name in enumerate(GLOBAL_FEATURES)})
        row.update({name: index + count / 10 for index, name in enumerate(COUNT_ITEM_FEATURES)})
        row.update({name: index + source / 100 for index, name in enumerate(CANDIDATE_FEATURES)})
        row["num_domains"] = count
        row["count_num_domains"] = count
        rows.append(row)
    # Repeated chain/count vectors are authoritative and byte-identical.
    for row in rows:
        if row["chain_id"] == "chainA":
            for index, name in enumerate(GLOBAL_FEATURES):
                row[name] = index + 0.25
    return rows


def _exact_raw_row(row: dict[str, object], chain_id: str | None = None) -> dict[str, str]:
    result = {field: str(row[field]) for field in required_dump_fields()}
    if chain_id is not None:
        result["chain_id"] = chain_id
    return result


def _reference(*aliases: str) -> dict[str, CathEntry]:
    entry = CathEntry(
        pdb_id="1abc",
        chain_id="A",
        entry_id="canonicalA",
        n_domains=2,
        n_residues=4,
        chopping="1-2:A|3-4:A",
        dataset="cath17287",
    )
    return {alias: entry for alias in (*aliases, entry.entry_id)}


def _scored_row(canonical: str, n_pred_domains: int) -> dict[str, object]:
    return {
        "delineation": canonical,
        "n_true_domains": 2,
        "n_pred_domains": n_pred_domains,
        "ndo": 0.1,
        "iou": 0.2,
        "boundary_f1_10": 0.3,
        "matched_dice": 0.4,
        "d_count_acc": 1.0,
        "S": 0.5,
        "is_oracle_s": 1,
    }


def _build_fixture_corpus(directory: Path, rows: list[dict[str, object]]) -> CorpusPaths:
    paths = CorpusPaths.at(directory)
    chain_rows: list[dict[str, object]] = []
    count_rows: list[dict[str, object]] = []
    for row in rows:
        chain_id = str(row["chain_id"])
        chain_rows.append(
            {"chain_id": chain_id, "n_true_domains": row["n_true_domains"], **{name: row[name] for name in GLOBAL_FEATURES}}
        )
        count_rows.append(
            {"chain_id": chain_id, **{name: row[name] for name in COUNT_ITEM_FEATURES}}
        )
    write_corpus(
        paths,
        chain_rows,
        count_rows,
        rows,
        [],
        {
            "dataset": "cath17287",
            "dataset_sha256": "a" * 64,
            "binary_sha256": "b" * 64,
            "git_commit": "c" * 40,
            "dump_argv_normalized": ["dump", "--out", "$OUT"],
            "build_argv_normalized": ["build", "--out-dir", "$OUT_DIR"],
            "seed": 37,
            "jobs": 1,
            "threads": 1,
            "feature_schema_hash": feature_schema_hash(),
        },
    )
    return paths


def test_shuffled_input_produces_identical_corpus_bytes(tmp_path: Path) -> None:
    rows = _raw_rows()
    first = _build_fixture_corpus(tmp_path / "first", rows)
    second = _build_fixture_corpus(tmp_path / "second", list(reversed(rows)))
    for name in ("chains", "counts", "candidates", "rejections", "manifest"):
        assert getattr(first, name).read_bytes() == getattr(second, name).read_bytes()


def test_manifest_changes_when_one_feature_changes(tmp_path: Path) -> None:
    rows = _raw_rows()
    changed = [dict(row) for row in rows]
    changed[0][CANDIDATE_FEATURES[-1]] = 999.125
    a = _build_fixture_corpus(tmp_path / "a", rows)
    b = _build_fixture_corpus(tmp_path / "b", changed)
    assert manifest_hash(a.manifest) != manifest_hash(b.manifest)


def test_exact_raw_header_has_one_frozen_num_domains() -> None:
    fields = required_dump_fields()
    assert len(fields) == 293
    assert fields == [
        "chain_id",
        "canonical_delineation",
        "source_index",
        "legacy_distance",
        *GLOBAL_FEATURES,
        *COUNT_ITEM_FEATURES,
        *CANDIDATE_FEATURES,
    ]
    assert fields.count("num_domains") == 1
    validate_dump_header(fields)


def test_merged_raw_identity_uniqueness_is_scoped_by_chain() -> None:
    first = _exact_raw_row(_raw_rows()[1], "chainA")
    second = dict(first)
    second["chain_id"] = "chainB"
    assert [row["chain_id"] for row in validate_dump_rows([second, first])] == [
        "chainA",
        "chainB",
    ]
    with pytest.raises(ValueError, match="duplicate candidate identity"):
        validate_dump_rows([first, dict(first)])


@pytest.mark.parametrize(
    "fields,match",
    [
        (lambda fields: [*fields, fields[-1]], "duplicate"),
        (lambda fields: fields[:-1], "missing|exact"),
        (lambda fields: [fields[1], fields[0], *fields[2:]], "order|exact"),
        (lambda fields: [*fields, "surprise"], "unknown|exact"),
        (lambda fields: [*fields, "merizo_score"], "external predictor"),
        (lambda fields: [*fields, "chainsaw_confidence"], "external predictor"),
    ],
)
def test_dump_header_fails_closed(fields, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        validate_dump_header(fields(required_dump_fields()))


def test_rejection_roles_have_same_exact_task1_schema() -> None:
    assert ACQUISITION_REJECTION_FIELDS == NORMALIZED_REJECTION_FIELDS == (
        "chain_id",
        "scope",
        "code",
        "detail",
        "delineation",
    )


def test_normalize_argv_removes_absolute_paths_and_preserves_order() -> None:
    first = normalize_argv(
        ["tool", "--dataset", "cath17287", "--dump", "/a/raw.csv", "--out-dir", "/a/out"]
    )
    second = normalize_argv(
        ["tool", "--dataset", "cath17287", "--dump", "/b/raw.csv", "--out-dir", "/b/out"]
    )
    assert first == second == [
        "tool",
        "--dataset",
        "cath17287",
        "--dump",
        "$DUMP",
        "--out-dir",
        "$OUT_DIR",
    ]


def test_corpus_newlines_float_format_ids_and_manifest_contract(tmp_path: Path) -> None:
    paths = _build_fixture_corpus(tmp_path / "corpus", _raw_rows())
    for path in (paths.chains, paths.counts, paths.candidates, paths.rejections, paths.manifest):
        data = path.read_bytes()
        assert b"\r\n" not in data
        assert data.endswith(b"\n")
    with paths.candidates.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert rows == sorted(
        rows,
        key=lambda row: (
            row["chain_id"],
            int(row["num_domains"]),
            row["canonical_delineation"],
            int(row["source_index"]),
        ),
    )
    for row in rows:
        expected = hashlib.sha256(
            (row["chain_id"] + "\0" + row["canonical_delineation"]).encode()
        ).hexdigest()
        assert row["candidate_id"] == expected
    manifest = json.loads(paths.manifest.read_text())
    encoded = paths.manifest.read_text()
    assert "timestamp" not in encoded
    assert "manifest_sha256" not in encoded
    assert str(tmp_path) not in encoded
    assert manifest["feature_schema_hash"] == feature_schema_hash()


def test_conflicting_repeated_vectors_fail_closed(tmp_path: Path) -> None:
    rows = _raw_rows()
    conflict = [dict(row) for row in rows]
    conflict[2][GLOBAL_FEATURES[0]] = -1
    with pytest.raises(ValueError, match="global"):
        _build_fixture_corpus(tmp_path / "global", conflict)

    conflict = [dict(row) for row in rows]
    extra = dict(conflict[1])
    extra["source_index"] = 99
    extra["canonical_delineation"] = "0;2 1;3"
    extra[COUNT_ITEM_FEATURES[1]] = -1
    conflict.append(extra)
    with pytest.raises(ValueError, match="count"):
        _build_fixture_corpus(tmp_path / "count", conflict)


def test_raw_population_conflicts_are_order_independent() -> None:
    rows = [{field: str(value) for field, value in row.items()} for row in _raw_rows()[1:]]
    rows[1][GLOBAL_FEATURES[0]] = "-1"
    for ordered in (rows, list(reversed(rows))):
        with pytest.raises(ValueError, match="conflicting repeated global"):
            validate_raw_population("chainA", ordered)

    rows = [{field: str(value) for field, value in row.items()} for row in _raw_rows()[1:]]
    duplicate = dict(rows[0])
    duplicate["canonical_delineation"] = "0;2 1;3"
    duplicate["source_index"] = "99"
    duplicate[COUNT_ITEM_FEATURES[1]] = "-1"
    rows.append(duplicate)
    for ordered in (rows, list(reversed(rows))):
        with pytest.raises(ValueError, match="conflicting repeated count"):
            validate_raw_population("chainA", ordered)


def test_manifest_rejects_nested_absolute_path(tmp_path: Path) -> None:
    rows = _raw_rows()
    paths = CorpusPaths.at(tmp_path / "bad-manifest")
    chain_rows = [
        {"chain_id": row["chain_id"], "n_true_domains": row["n_true_domains"], **{name: row[name] for name in GLOBAL_FEATURES}}
        for row in rows
    ]
    count_rows = [
        {"chain_id": row["chain_id"], **{name: row[name] for name in COUNT_ITEM_FEATURES}}
        for row in rows
    ]
    with pytest.raises(ValueError, match="absolute path"):
        write_corpus(
            paths,
            chain_rows,
            count_rows,
            rows,
            [],
            {
                "dataset": "cath17287",
                "dataset_sha256": "a" * 64,
                "binary_sha256": "b" * 64,
                "git_commit": "c" * 40,
                "dump_argv_normalized": ["dump", {"nested": "/private/raw.csv"}],
                "build_argv_normalized": ["build"],
                "seed": 37,
                "jobs": 1,
                "threads": 1,
            },
        )


def test_nonfinite_and_missing_values_are_rejected(tmp_path: Path) -> None:
    rows = _raw_rows()
    rows[0][CANDIDATE_FEATURES[-1]] = float("nan")
    with pytest.raises(ValueError, match="finite"):
        _build_fixture_corpus(tmp_path / "bad", rows)


def test_resume_part_requires_exact_nonempty_schema_and_identity(tmp_path: Path) -> None:
    fields = required_dump_fields()
    raw = {field: str(value) for field, value in _raw_rows()[1].items() if field in fields}
    raw["chain_id"] = "chainA"

    def write(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    valid = tmp_path / "chainA.csv"
    write(valid, fields, [raw])
    assert valid_resume_part(valid, "chainA")

    header_only = tmp_path / "header.csv"
    write(header_only, fields, [])
    assert not valid_resume_part(header_only, "header")

    stale = tmp_path / "stale.csv"
    write(stale, fields[:-1], [{key: value for key, value in raw.items() if key in fields[:-1]}])
    assert not valid_resume_part(stale, "chainA")

    duplicate = tmp_path / "duplicate.csv"
    write(duplicate, [*fields, fields[-1]], [raw])
    assert not valid_resume_part(duplicate, "chainA")

    wrong = tmp_path / "wrong.csv"
    write(wrong, fields, [raw])
    assert not valid_resume_part(wrong, "other")

    malformed = tmp_path / "malformed.csv"
    bad = dict(raw)
    bad[CANDIDATE_FEATURES[-1]] = "nan"
    write(malformed, fields, [bad])
    assert not valid_resume_part(malformed, "chainA")

    surplus = tmp_path / "surplus.csv"
    write(surplus, fields, [raw])
    lines = surplus.read_text().splitlines()
    surplus.write_text(lines[0] + "\n" + lines[1] + ",extra\n")
    assert not valid_resume_part(surplus, "chainA")


def test_acquisition_only_failure_writes_no_accepted_rows(tmp_path: Path) -> None:
    paths = CorpusPaths.at(tmp_path / "only-rejection")
    write_corpus(
        paths,
        [],
        [],
        [],
        [{"chain_id": "bad", "scope": "chain", "code": "schema_mismatch", "detail": "header-only", "delineation": ""}],
        {
            "dataset": "cath17287",
            "dataset_sha256": "a" * 64,
            "binary_sha256": "b" * 64,
            "git_commit": "c" * 40,
            "dump_argv_normalized": ["dump"],
            "build_argv_normalized": ["build"],
            "seed": 37,
            "jobs": 1,
            "threads": 1,
        },
    )
    assert paths.chains.read_text().count("\n") == 1
    assert paths.counts.read_text().count("\n") == 1
    assert paths.candidates.read_text().count("\n") == 1
    assert paths.rejections.read_text().count("\n") == 2


def test_candidate_rejection_keeps_complete_chain_count_population(tmp_path: Path) -> None:
    raw_rows = [row for row in _raw_rows() if row["chain_id"] == "chainA"]
    paths = CorpusPaths.at(tmp_path / "partial-candidates")
    write_corpus(
        paths,
        [
            {"chain_id": row["chain_id"], "n_true_domains": row["n_true_domains"], **{name: row[name] for name in GLOBAL_FEATURES}}
            for row in raw_rows
        ],
        [
            {"chain_id": row["chain_id"], **{name: row[name] for name in COUNT_ITEM_FEATURES}}
            for row in raw_rows
        ],
        raw_rows[:1],
        [{"chain_id": "chainA", "scope": "candidate", "code": "candidate_parse_failed", "detail": "invalid partition", "delineation": raw_rows[1]["canonical_delineation"]}],
        {
            "dataset": "cath17287",
            "dataset_sha256": "a" * 64,
            "binary_sha256": "b" * 64,
            "git_commit": "c" * 40,
            "dump_argv_normalized": ["dump"],
            "build_argv_normalized": ["build"],
            "seed": 37,
            "jobs": 1,
            "threads": 1,
        },
    )
    with paths.counts.open(newline="") as handle:
        assert len(list(csv.DictReader(handle))) == 2
    with paths.candidates.open(newline="") as handle:
        assert len(list(csv.DictReader(handle))) == 1


def test_normalize_raw_rows_preserves_raw_precision_and_complete_counts(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    rows = [_exact_raw_row(row, "rawA") for row in _raw_rows()[1:]]
    precise_field = CANDIDATE_FEATURES[-1]
    rows[0][precise_field] = "0.12345678901234567"
    rejected = RejectionRecord(
        chain_id="canonicalA",
        scope="candidate",
        code=RejectionCode.CANDIDATE_PARSE_FAILED,
        detail="invalid candidate",
        delineation=rows[1]["canonical_delineation"],
    )

    def score(*_args, **_kwargs) -> ScoredChain:
        return ScoredChain(
            rows=[_scored_row(rows[0]["canonical_delineation"], 1)],
            rejections=[rejected],
        )

    monkeypatch.setattr(build_factorized_corpus, "_score_candidates", score)
    chains, counts, candidates, rejections = build_factorized_corpus.normalize_raw_rows(
        rows, _reference("rawA"), tmp_path
    )
    assert len(chains) == 1
    assert [int(float(row["count_num_domains"])) for row in counts] == [1, 2]
    assert len(candidates) == 1
    assert candidates[0][precise_field] == "0.12345678901234567"
    assert rejections == [rejected]


@pytest.mark.parametrize("defect", ["missing", "extra", "duplicate"])
def test_normalize_raw_rows_rejects_non_bijective_scoring_identity(
    defect: str, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    rows = [_exact_raw_row(row, "rawA") for row in _raw_rows()[1:]]
    first = _scored_row(rows[0]["canonical_delineation"], 1)
    second = _scored_row(rows[1]["canonical_delineation"], 2)
    if defect == "missing":
        scored_rows = [first]
    elif defect == "extra":
        scored_rows = [first, second, _scored_row("0;2 1;3", 2)]
    else:
        scored_rows = [first, second, dict(first)]

    monkeypatch.setattr(
        build_factorized_corpus,
        "_score_candidates",
        lambda *_args, **_kwargs: ScoredChain(rows=scored_rows, rejections=[]),
    )
    chains, counts, candidates, rejections = build_factorized_corpus.normalize_raw_rows(
        rows, _reference("rawA"), tmp_path
    )
    assert chains == counts == candidates == []
    assert len(rejections) == 1
    assert rejections[0]["code"] == "schema_mismatch"
    assert rejections[0]["detail"] == "Task 1 scoring identity join is not one-to-one"


def test_normalize_raw_rows_rejects_canonical_alias_mixture_before_scoring(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    rows = [
        _exact_raw_row(_raw_rows()[1], "rawA"),
        _exact_raw_row(_raw_rows()[2], "rawAlias"),
    ]
    calls = 0

    def score(*_args, **_kwargs) -> ScoredChain:
        nonlocal calls
        calls += 1
        return ScoredChain(rows=[], rejections=[])

    monkeypatch.setattr(build_factorized_corpus, "_score_candidates", score)
    outcomes = []
    for ordered in (rows, list(reversed(rows))):
        result = build_factorized_corpus.normalize_raw_rows(
            ordered, _reference("rawA", "rawAlias"), tmp_path
        )
        outcomes.append(result)
    assert calls == 0
    assert outcomes[0] == outcomes[1]
    chains, counts, candidates, rejections = outcomes[0]
    assert chains == counts == candidates == []
    assert rejections == [
        {
            "chain_id": "canonicalA",
            "scope": "chain",
            "code": "schema_mismatch",
            "detail": "canonical chain is sourced from multiple raw chain identities",
            "delineation": "",
        }
    ]


def _provenance_fixture(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    dataset = tmp_path / "dataset.csv"
    dataset.write_bytes(b"dataset\n")
    binary = tmp_path / "sword2"
    binary.write_bytes(b"binary\n")
    expected_commit = "a" * 40
    monkeypatch.setattr(build_factorized_corpus, "dataset_path", lambda _name: dataset)
    monkeypatch.setattr(
        build_factorized_corpus,
        "_expected_git_commit",
        lambda: expected_commit,
        raising=False,
    )
    provenance = {
        "dataset": "cath17287",
        "dataset_sha256": hashlib.sha256(dataset.read_bytes()).hexdigest(),
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "git_commit": expected_commit,
        "seed": 37,
        "jobs": 3,
        "threads": 4,
        "timeout": 123,
        "limit": 2,
        "resume": False,
        "feature_schema_hash": feature_schema_hash(),
        "dump_argv_normalized": [
            "benchmark.dump_candidate_corpus",
            "--dataset=cath17287",
            "--limit=2",
            "--jobs=3",
            "--threads=4",
            "--timeout=123",
            "--seed=37",
            "--binary=$BINARY",
            "--chain-cache-dir=$CHAIN_CACHE_DIR",
            "--out=$OUT",
        ],
    }
    path = tmp_path / "dump.provenance.json"
    path.write_text(json.dumps(provenance, sort_keys=True, separators=(",", ":")) + "\n")
    return path, binary, provenance


def test_provenance_accepts_strict_canonical_equal_form(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, binary, provenance = _provenance_fixture(tmp_path, monkeypatch)
    assert build_factorized_corpus._load_and_verify_provenance(
        path, "cath17287", binary
    ) == provenance


def test_dump_argv_parser_accepts_separate_form_and_applies_defaults() -> None:
    parsed = build_factorized_corpus._parse_dump_argv_normalized(
        [
            "benchmark.dump_candidate_corpus",
            "--dataset",
            "cath17287",
            "--jobs",
            "8",
            "--out",
            "$OUT",
        ]
    )
    assert parsed == {
        "dataset": "cath17287",
        "out": "$OUT",
        "parts_dir": None,
        "chain_cache_dir": "$CHAIN_CACHE_DIR",
        "binary": "$BINARY",
        "jobs": 8,
        "threads": 1,
        "timeout": 300,
        "limit": None,
        "seed": 37,
        "resume": False,
    }


@pytest.mark.parametrize(
    "tamper",
    [
        lambda value: value.update(git_commit="b" * 40),
        lambda value: value.update(seed=38),
        lambda value: value.update(jobs=5),
        lambda value: value.update(threads=5),
        lambda value: value.update(timeout=124),
        lambda value: value.update(limit=3),
        lambda value: value.update(resume=0),
        lambda value: value["dump_argv_normalized"].append("--unknown=value"),
        lambda value: value["dump_argv_normalized"].extend(["--jobs", "3"]),
        lambda value: value["dump_argv_normalized"].__setitem__(1, "--dataset=other"),
        lambda value: value["dump_argv_normalized"].__setitem__(-1, "--out=/absolute/raw.csv"),
        lambda value: value["dump_argv_normalized"].__setitem__(-3, "--binary=/absolute/sword2"),
    ],
)
def test_provenance_rejects_each_sidecar_or_argv_tamper(
    tamper, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path, binary, provenance = _provenance_fixture(tmp_path, monkeypatch)
    tamper(provenance)
    path.write_text(json.dumps(provenance, sort_keys=True, separators=(",", ":")) + "\n")
    with pytest.raises(ValueError, match="provenance|argv|option|role|commit|mismatch"):
        build_factorized_corpus._load_and_verify_provenance(path, "cath17287", binary)
