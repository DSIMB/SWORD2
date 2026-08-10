from __future__ import annotations

import csv
import hashlib
import io
import json
from dataclasses import replace
from pathlib import Path

import pytest

import benchmark.factorized_ranker.folds as folds
from benchmark.datasets import CathEntry, cath_family_combination, cath_family_labels
from benchmark.factorized_ranker.corpus import CorpusPaths, write_corpus
from benchmark.factorized_ranker.folds import (
    FoldAssignment,
    assign_folds,
    chain_label_cohorts,
    connected_components,
    individual_label_seen,
    load_fold_manifest,
    validate_folds,
    write_fold_manifest,
)
from benchmark.factorized_ranker.schema import (
    CANDIDATE_FEATURES,
    COUNT_ITEM_FEATURES,
    GLOBAL_FEATURES,
)


def _entry(
    entry_id: str,
    pdb_id: str,
    labels: tuple[str, ...] = (),
    *,
    n_domains: int = 2,
    n_residues: int = 300,
    dataset: str = "cath17287",
) -> CathEntry:
    if labels:
        chopping = "|".join(
            f"{index * 10 + 1}-{index * 10 + 10}:{label}"
            for index, label in enumerate(labels)
        )
    else:
        chopping = "1-10:999_999"
    return CathEntry(
        pdb_id=pdb_id,
        chain_id=entry_id[-1:],
        entry_id=entry_id,
        n_domains=n_domains,
        n_residues=n_residues,
        chopping=chopping,
        dataset=dataset,
    )


def _independent_entries() -> list[CathEntry]:
    return [
        _entry("entryA", "1aaa", ("fam.a",), n_domains=1, n_residues=200),
        _entry("entryB", "2bbb", ("fam.b",), n_domains=2, n_residues=250),
        _entry("entryC", "3ccc", ("fam.c",), n_domains=3, n_residues=350),
        _entry("entryD", "4ddd", ("fam.d",), n_domains=4, n_residues=450),
        _entry("entryE", "5eee", ("fam.e",), n_domains=5, n_residues=449),
        _entry("entryF", "6fff", (), n_domains=6, n_residues=349),
    ]


def test_family_combination_is_sorted_multiset_without_sentinel() -> None:
    chopping = "1-10:2.40.50.140|11-20:1.10.8.10|21-30:2.40.50.140|31-40:999_999"
    assert cath_family_labels(chopping) == (
        "2.40.50.140",
        "1.10.8.10",
        "2.40.50.140",
    )
    assert cath_family_combination(chopping) == (
        "1.10.8.10",
        "2.40.50.140",
        "2.40.50.140",
    )
    assert cath_family_combination("1-2|3-4: |5-6:999_999") == ()


def test_transitive_pdb_and_combination_links_share_component() -> None:
    entries = [
        _entry("A", "1abc", ("x",)),
        _entry("B", " 1ABC ", ("y",)),
        _entry("C", "2def", ("y",)),
    ]
    components = connected_components(entries)
    assert [[entry.entry_id for entry in component] for component in components] == [["A", "B", "C"]]


def test_empty_combinations_do_not_union_but_pdb_links_still_do() -> None:
    entries = [
        _entry("A", "1abc"),
        _entry("B", "2def"),
        _entry("C", "1ABC", ()),
    ]
    components = connected_components(entries)
    assert sorted(tuple(entry.entry_id for entry in component) for component in components) == [
        ("A", "C"),
        ("B",),
    ]


def test_no_pdb_or_exact_combination_crosses_folds() -> None:
    entries = _independent_entries() + [
        _entry("entryG", "1AAA", ("fam.g",), n_domains=2, n_residues=240),
        _entry("entryH", "8hhh", ("fam.b",), n_domains=3, n_residues=360),
    ]
    assignments = assign_folds(entries, n_folds=5, seed=37)
    validate_folds(assignments)
    assert len({assignment.fold for assignment in assignments}) == 5
    by_pdb: dict[str, set[int]] = {}
    by_combo: dict[tuple[str, ...], set[int]] = {}
    for assignment in assignments:
        by_pdb.setdefault(assignment.pdb_id, set()).add(assignment.fold)
        if assignment.family_combination:
            by_combo.setdefault(assignment.family_combination, set()).add(assignment.fold)
    assert all(len(folds_) == 1 for folds_ in (*by_pdb.values(), *by_combo.values()))


def test_assignments_are_order_independent_and_cover_exact_bins() -> None:
    entries = _independent_entries()
    forward = assign_folds(entries, n_folds=5, seed=37)
    reverse = assign_folds(list(reversed(entries)), n_folds=5, seed=37)
    assert forward == reverse
    assert {assignment.true_count_bin for assignment in forward} == {"1", "2", "3", "4", "5+"}
    assert {assignment.length_bin for assignment in forward} == {
        "<250",
        "250-349",
        "350-449",
        "450+",
    }
    ordered_components = sorted(
        (assignment.component_id for assignment in forward),
        key=lambda component_id: hashlib.sha256(f"37:{component_id}".encode()).hexdigest(),
    )
    by_component = {assignment.component_id: assignment.fold for assignment in forward}
    assert by_component[ordered_components[0]] == 0


def test_greedy_assignment_reserves_enough_components_to_populate_every_fold() -> None:
    strata = [
        (5, 106),
        (7, 491),
        (3, 93),
        (2, 21),
        (7, 563),
        (5, 784),
        (1, 228),
        (6, 284),
        (3, 109),
        (5, 220),
        (1, 657),
        (5, 279),
    ]
    entries = [
        _entry(
            f"skew{index}",
            f"{index:04x}",
            (f"unique.{index}",),
            n_domains=count,
            n_residues=length,
        )
        for index, (count, length) in enumerate(strata)
    ]
    assignments = assign_folds(entries, n_folds=5, seed=37)
    assert {assignment.fold for assignment in assignments} == set(range(5))


def test_greedy_assignment_balances_many_independent_equal_components() -> None:
    entries = [
        _entry(
            f"equal{index}",
            f"{index:04x}",
            (f"unique.{index}",),
            n_domains=2,
            n_residues=300,
        )
        for index in range(50)
    ]

    assignments = assign_folds(entries, n_folds=5, seed=37)
    fold_sizes = [
        sum(assignment.fold == fold for assignment in assignments)
        for fold in range(5)
    ]

    assert fold_sizes == [10, 10, 10, 10, 10]


@pytest.mark.parametrize(
    "entry,match",
    [
        (_entry("", "1abc"), "entry"),
        (_entry("bad\x00id", "1abc"), "entry"),
        (_entry("bad\nid", "1abc"), "entry"),
        (_entry("A", "abc"), "PDB"),
        (_entry("A", "1-bc"), "PDB"),
        (_entry("A", "åbcd"), "PDB"),
        (_entry("A", "1abc", n_domains=0), "count"),
        (_entry("A", "1abc", n_residues=0), "length"),
    ],
)
def test_invalid_entry_metadata_fails_closed(entry: CathEntry, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        connected_components([entry])


def test_invalid_fold_parameters_and_crossed_assignments_fail_closed() -> None:
    entries = _independent_entries()
    for n_folds in (True, 1, 7):
        with pytest.raises(ValueError):
            assign_folds(entries, n_folds=n_folds, seed=37)
    for seed in (True, -1):
        with pytest.raises(ValueError):
            assign_folds(entries, n_folds=5, seed=seed)

    assignments = list(assign_folds(entries, n_folds=5, seed=37))
    with pytest.raises(ValueError, match="duplicate"):
        validate_folds([*assignments, assignments[0]])
    crossed = [
        replace(assignments[0], pdb_id=assignments[1].pdb_id),
        *assignments[1:],
    ]
    if crossed[0].fold == crossed[1].fold:
        crossed[0] = replace(crossed[0], fold=(crossed[1].fold + 1) % 5)
    with pytest.raises(ValueError):
        validate_folds(crossed)
    with pytest.raises(ValueError):
        validate_folds([replace(assignments[0], fold=9), *assignments[1:]])


def test_duplicate_and_nonstring_entry_ids_fail_as_validation_errors() -> None:
    duplicate = [_entry("same", "1abc"), _entry("same", "2def")]
    with pytest.raises(ValueError, match="duplicate"):
        connected_components(duplicate)
    malformed = replace(_entry("valid", "1abc"), entry_id=None)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="entry"):
        connected_components([malformed])


def test_seen_flags_and_chain_cohorts_use_only_other_folds() -> None:
    assignments = (
        FoldAssignment("train", "1aaa", ("known",), ("known",), "1", "<250", "", 0),
        FoldAssignment("seen", "2bbb", ("known", "known"), ("known", "known"), "2", "250-349", "", 1),
        FoldAssignment("unseen", "3ccc", ("known", "new"), ("known", "new"), "3", "350-449", "", 1),
        FoldAssignment("unknown", "4ddd", (), (), "4", "450+", "", 1),
    )
    # Component IDs depend only on each singleton chain ID.
    assignments = tuple(
        replace(
            assignment,
            component_id=hashlib.sha256(assignment.chain_id.encode()).hexdigest(),
        )
        for assignment in assignments
    )
    validate_folds(assignments)
    assert individual_label_seen(assignments, 1) == {"known": True, "new": False}
    assert chain_label_cohorts(assignments, 1) == {
        "seen": "seen",
        "unknown": "unknown",
        "unseen": "unseen",
    }


def _valid_assignments() -> tuple[FoldAssignment, ...]:
    return assign_folds(_independent_entries(), n_folds=5, seed=37)


def test_fold_manifest_is_canonical_external_hash_and_order_independent(tmp_path: Path) -> None:
    assignments = _valid_assignments()
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    expected = {
        "dataset_sha256": "a" * 64,
        "corpus_sha256": "b" * 64,
        "chains_sha256": "c" * 64,
    }
    first_hash = write_fold_manifest(first, assignments, seed=37, **expected)
    second_hash = write_fold_manifest(second, list(reversed(assignments)), seed=37, **expected)
    assert first.read_bytes() == second.read_bytes()
    assert first_hash == second_hash == hashlib.sha256(first.read_bytes()).hexdigest()
    assert first.read_bytes().endswith(b"\n") and b"\n " not in first.read_bytes()
    payload = json.loads(first.read_text())
    assert "manifest_sha256" not in payload
    loaded = load_fold_manifest(
        first,
        expected_sha256=first_hash,
        expected_dataset_sha256="a" * 64,
        expected_corpus_sha256="b" * 64,
        expected_chains_sha256="c" * 64,
    )
    assert loaded.assignments == assignments


@pytest.mark.parametrize("defect", ["whitespace", "unknown", "record_order", "hash"])
def test_fold_manifest_rejects_noncanonical_or_tampered_bytes(
    defect: str, tmp_path: Path
) -> None:
    path = tmp_path / "folds.json"
    manifest_hash = write_fold_manifest(
        path,
        _valid_assignments(),
        "a" * 64,
        "b" * 64,
        chains_sha256="c" * 64,
    )
    payload = json.loads(path.read_text())
    if defect == "whitespace":
        path.write_text(json.dumps(payload, indent=2) + "\n")
    elif defect == "unknown":
        payload["unknown"] = 1
        path.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")
    elif defect == "record_order":
        payload["assignments"].reverse()
        path.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")
    else:
        payload["dataset_sha256"] = "z" * 64
        path.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n")
    with pytest.raises(ValueError):
        load_fold_manifest(path, expected_sha256=manifest_hash if defect != "hash" else None)


def _corpus_rows(entries: list[CathEntry]):
    chains = []
    counts = []
    candidates = []
    for index, entry in enumerate(entries):
        chains.append(
            {
                "chain_id": entry.entry_id,
                "n_true_domains": entry.n_domains,
                **{field: 0.0 for field in GLOBAL_FEATURES},
            }
        )
        count_row = {field: 0.0 for field in COUNT_ITEM_FEATURES}
        count_row["chain_id"] = entry.entry_id
        count_row["count_num_domains"] = entry.n_domains
        counts.append(count_row)
        candidate_row = {field: 0.0 for field in CANDIDATE_FEATURES}
        candidate_row.update(
            {
                "chain_id": entry.entry_id,
                "canonical_delineation": f"0-{entry.n_residues - 1}",
                "source_index": index,
                "legacy_distance": 0.0,
                "num_domains": entry.n_domains,
                "n_true_domains": entry.n_domains,
                "n_pred_domains": entry.n_domains,
                "ndo": 1.0,
                "iou": 1.0,
                "boundary_f1_10": 1.0,
                "matched_dice": 1.0,
                "d_count_acc": 1.0,
                "S": 1.0,
                "is_oracle_s": 1,
            }
        )
        candidates.append(candidate_row)
    return chains, counts, candidates


def _build_corpus(
    tmp_path: Path,
    entries: list[CathEntry],
    dataset_file: Path,
    *,
    reverse_input_rows: bool = False,
) -> Path:
    corpus_dir = tmp_path / "corpus"
    provenance = {
        "dataset": "cath17287",
        "dataset_sha256": hashlib.sha256(dataset_file.read_bytes()).hexdigest(),
        "binary_sha256": "a" * 64,
        "git_commit": "b" * 40,
        "dump_argv_normalized": ["benchmark.dump_candidate_corpus"],
        "build_argv_normalized": ["benchmark.build_factorized_corpus"],
        "seed": 37,
        "jobs": 1,
        "threads": 1,
        "timeout": 300,
        "limit": None,
        "resume": False,
    }
    chains, counts, candidates = _corpus_rows(entries)
    if reverse_input_rows:
        chains.reverse()
        counts.reverse()
        candidates.reverse()
    write_corpus(
        CorpusPaths.at(corpus_dir),
        chains,
        counts,
        candidates,
        [],
        provenance,
    )
    return corpus_dir


def _rewrite_table_manifest(
    corpus_dir: Path,
    table_name: str,
    fields: tuple[str, ...],
    rows: list[dict[str, str]],
    *,
    accepted_unique_chains: int | None = None,
) -> None:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = output.getvalue().encode()
    (corpus_dir / f"{table_name}.csv").write_bytes(data)
    manifest_path = corpus_dir / "corpus_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["tables"][table_name] = {
        "sha256": hashlib.sha256(data).hexdigest(),
        "rows": len(rows),
    }
    if accepted_unique_chains is not None:
        manifest["accepted_unique_chains"] = accepted_unique_chains
    manifest_path.write_text(
        json.dumps(manifest, sort_keys=True, separators=(",", ":")) + "\n"
    )


def test_load_accepted_entries_verifies_hashes_and_excludes_source_only_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset_file = tmp_path / "CATH-17287.csv"
    dataset_file.write_text("source bytes\n")
    accepted = _independent_entries()[:5]
    source = [*accepted, _entry("sourceOnly", "9zzz", ("unused",))]
    corpus_dir = _build_corpus(tmp_path, accepted, dataset_file)
    monkeypatch.setattr(folds, "dataset_path", lambda _dataset: dataset_file)
    monkeypatch.setattr(folds, "load_dataset", lambda _dataset: list(reversed(source)))

    loaded, dataset_hash, corpus_hash, chains_hash = folds._load_accepted_entries(
        corpus_dir, "cath17287"
    )
    assert [entry.entry_id for entry in loaded] == sorted(entry.entry_id for entry in accepted)
    assert "sourceOnly" not in {entry.entry_id for entry in loaded}
    assert dataset_hash == hashlib.sha256(dataset_file.read_bytes()).hexdigest()
    assert corpus_hash == hashlib.sha256((corpus_dir / "corpus_manifest.json").read_bytes()).hexdigest()
    assert chains_hash == hashlib.sha256((corpus_dir / "chains.csv").read_bytes()).hexdigest()

    with (corpus_dir / "counts.csv").open("a") as handle:
        handle.write("tamper\n")
    with pytest.raises(ValueError, match="hash"):
        folds._load_accepted_entries(corpus_dir, "cath17287")


def test_reversed_chains_and_metadata_produce_identical_fold_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset_file = tmp_path / "CATH-17287.csv"
    dataset_file.write_text("source bytes\n")
    accepted = _independent_entries()
    first_corpus = _build_corpus(tmp_path / "first", accepted, dataset_file)
    second_corpus = _build_corpus(
        tmp_path / "second", accepted, dataset_file, reverse_input_rows=True
    )
    monkeypatch.setattr(folds, "dataset_path", lambda _dataset: dataset_file)
    monkeypatch.setattr(folds, "load_dataset", lambda _dataset: accepted)
    loaded, dataset_hash, corpus_hash, chains_hash = folds._load_accepted_entries(
        first_corpus, "cath17287"
    )
    first = tmp_path / "first-folds.json"
    write_fold_manifest(
        first,
        assign_folds(loaded),
        dataset_hash,
        corpus_hash,
        chains_sha256=chains_hash,
    )
    monkeypatch.setattr(folds, "load_dataset", lambda _dataset: list(reversed(accepted)))
    loaded, dataset_hash, corpus_hash, chains_hash = folds._load_accepted_entries(
        second_corpus, "cath17287"
    )
    second = tmp_path / "second-folds.json"
    write_fold_manifest(
        second,
        assign_folds(loaded),
        dataset_hash,
        corpus_hash,
        chains_sha256=chains_hash,
    )
    assert first.read_bytes() == second.read_bytes()


@pytest.mark.parametrize("defect", ["duplicate", "empty", "dataset_hash"])
def test_load_accepted_entries_rejects_accepted_population_integrity_defects(
    defect: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset_file = tmp_path / "CATH-17287.csv"
    dataset_file.write_text("source bytes\n")
    accepted = _independent_entries()[:5]
    corpus_dir = _build_corpus(tmp_path, accepted, dataset_file)
    monkeypatch.setattr(folds, "dataset_path", lambda _dataset: dataset_file)
    monkeypatch.setattr(folds, "load_dataset", lambda _dataset: accepted)
    if defect == "dataset_hash":
        dataset_file.write_text("changed source bytes\n")
    else:
        with (corpus_dir / "chains.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        if defect == "duplicate":
            rows.append(dict(rows[0]))
            accepted_count = len(accepted)
        else:
            rows = []
            accepted_count = 0
        _rewrite_table_manifest(
            corpus_dir,
            "chains",
            tuple(("chain_id", "n_true_domains", *GLOBAL_FEATURES)),
            rows,
            accepted_unique_chains=accepted_count,
        )
    with pytest.raises(ValueError):
        folds._load_accepted_entries(corpus_dir, "cath17287")


@pytest.mark.parametrize("defect", ["missing", "duplicate", "count", "dataset"])
def test_load_accepted_entries_rejects_metadata_join_defects(
    defect: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset_file = tmp_path / "CATH-17287.csv"
    dataset_file.write_text("source bytes\n")
    accepted = _independent_entries()[:5]
    corpus_dir = _build_corpus(tmp_path, accepted, dataset_file)
    metadata = list(accepted)
    if defect == "missing":
        metadata.pop()
    elif defect == "duplicate":
        metadata.append(metadata[0])
    elif defect == "count":
        metadata[0] = replace(metadata[0], n_domains=metadata[0].n_domains + 1)
    else:
        metadata[0] = replace(metadata[0], dataset="other")
    monkeypatch.setattr(folds, "dataset_path", lambda _dataset: dataset_file)
    monkeypatch.setattr(folds, "load_dataset", lambda _dataset: metadata)
    with pytest.raises(ValueError):
        folds._load_accepted_entries(corpus_dir, "cath17287")


def test_fold_cli_rejects_output_alias(tmp_path: Path) -> None:
    from benchmark.build_factorized_folds import output_aliases_corpus_input

    corpus_dir = tmp_path / "corpus"
    corpus_dir.mkdir()
    manifest = corpus_dir / "corpus_manifest.json"
    manifest.write_text("{}\n")
    assert output_aliases_corpus_input(manifest, corpus_dir)
    assert not output_aliases_corpus_input(tmp_path / "folds.json", corpus_dir)
