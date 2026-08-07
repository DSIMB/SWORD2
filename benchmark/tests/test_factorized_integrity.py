from pathlib import Path

import pytest

from benchmark.datasets import CathEntry
from benchmark.factorized_ranker.integrity import (
    RejectionCode,
    map_cath_reference,
    validate_partition,
)


def _write_pdb(path: Path, author_numbers: list[int], chain_id: str = "A") -> None:
    lines = [
        (
            f"ATOM  {serial:5d}  CA  ALA {chain_id}{author_number:4d}    "
            f"{3.8 * (serial - 1):8.3f}{0.0:8.3f}{0.0:8.3f}"
            "  1.00 20.00           C\n"
        )
        for serial, author_number in enumerate(author_numbers, start=1)
    ]
    path.write_text("".join(lines) + "END\n")


def test_partition_requires_exact_coverage():
    with pytest.raises(
        ValueError,
        match=RejectionCode.CANDIDATE_INCOMPLETE_COVERAGE.value,
    ):
        validate_partition("0-2 4-5", n_residues=6, declared_domains=2)


@pytest.mark.parametrize(
    ("delineation", "code"),
    [
        ("0-3 3-5", RejectionCode.CANDIDATE_OVERLAP),
        ("0-2 3-6", RejectionCode.CANDIDATE_OUT_OF_RANGE),
        ("2-0 1-5", RejectionCode.CANDIDATE_PARSE_FAILED),
        ("0-2 3-5", RejectionCode.CANDIDATE_DOMAIN_COUNT_MISMATCH),
    ],
)
def test_partition_rejects_invalid_candidates(delineation, code):
    declared = 3 if code is RejectionCode.CANDIDATE_DOMAIN_COUNT_MISMATCH else 2
    with pytest.raises(ValueError, match=code.value):
        validate_partition(delineation, n_residues=6, declared_domains=declared)


@pytest.mark.parametrize("delineation", ["", "0-2;", "0-1 nope"])
def test_partition_rejects_empty_or_malformed_tokens(delineation):
    with pytest.raises(ValueError, match=RejectionCode.CANDIDATE_PARSE_FAILED.value):
        validate_partition(delineation, n_residues=3, declared_domains=1)


def test_discontinuous_partition_has_stable_canonical_key():
    parsed = validate_partition("4-5 0-1;3-3 2-2", 6, 3)

    assert parsed.canonical_delineation == "0-1;3 2 4-5"
    assert parsed.domains == (((0, 1), (3, 3)), ((2, 2),), ((4, 5),))
    assert parsed.residue_to_domain == (0, 0, 1, 0, 2, 2)


def test_adjacent_segments_have_one_topology_canonical_key():
    split = validate_partition("0-1;2-3 4-5", 6, 2)
    merged = validate_partition("0-3 4-5", 6, 2)

    assert split.canonical_delineation == merged.canonical_delineation == "0-3 4-5"
    assert split.domains == merged.domains
    assert split.residue_to_domain == merged.residue_to_domain


@pytest.mark.parametrize(
    "delineation",
    ["+0-+1 2-3", "０-１ ２-３", "0_0-1 2-3"],
)
def test_partition_accepts_only_ascii_unsigned_integer_tokens(delineation):
    with pytest.raises(ValueError, match=RejectionCode.CANDIDATE_PARSE_FAILED.value):
        validate_partition(delineation, n_residues=4, declared_domains=2)


def test_cath_reference_is_mapped_to_zero_based_chain_coordinates(tmp_path):
    pdb_path = tmp_path / "example.pdb"
    _write_pdb(pdb_path, [10, 11, 20])
    entry = CathEntry(
        pdb_id="test",
        chain_id="A",
        entry_id="testA",
        n_domains=2,
        n_residues=3,
        chopping="10-11:1_1|20:2_2",
        dataset="test",
    )

    reference = map_cath_reference(entry, pdb_path)

    assert reference.chopping == "0-1,2"
    assert reference.n_residues == 3
    assert reference.n_domains == 2


@pytest.mark.parametrize("chopping", ["10||11", "10|11|", "10__11"])
def test_cath_reference_rejects_empty_domains_or_segments(tmp_path, chopping):
    pdb_path = tmp_path / "example.pdb"
    _write_pdb(pdb_path, [10, 11])
    entry = CathEntry(
        pdb_id="test",
        chain_id="A",
        entry_id="testA",
        n_domains=2 if "|" in chopping else 1,
        n_residues=2,
        chopping=chopping,
        dataset="test",
    )

    with pytest.raises(ValueError, match=RejectionCode.AUTHOR_MAPPING_FAILED.value):
        map_cath_reference(entry, pdb_path)


def test_cath_reference_rejects_overlapping_mapped_domains(tmp_path):
    pdb_path = tmp_path / "example.pdb"
    _write_pdb(pdb_path, [10, 11])
    entry = CathEntry(
        pdb_id="test",
        chain_id="A",
        entry_id="testA",
        n_domains=2,
        n_residues=2,
        chopping="10-11:1_1|10-11:2_2",
        dataset="test",
    )

    with pytest.raises(ValueError, match=RejectionCode.AUTHOR_MAPPING_FAILED.value):
        map_cath_reference(entry, pdb_path)
