import csv
from pathlib import Path

import pytest

from benchmark.build_training_table import (
    BASE_CANDIDATE_FIELDS,
    FIELDNAMES,
    _load_reference,
    _score_chain_file,
    _score_candidates,
    build_from_dump_csv,
    write_rejections,
)
from benchmark.candidate_geometry import CANDIDATE_GEOMETRY_FIELDS
from benchmark.datasets import CathEntry
from benchmark.factorized_ranker.integrity import RejectionCode, RejectionRecord


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


def _entry(
    *,
    n_residues: int,
    n_domains: int,
    chopping: str,
    chain_id: str = "A",
) -> CathEntry:
    return CathEntry(
        pdb_id="test",
        chain_id=chain_id,
        entry_id="testchain",
        n_domains=n_domains,
        n_residues=n_residues,
        chopping=chopping,
        dataset="test",
    )


def _candidate(
    *,
    delineation: str,
    num_domains: int,
    output_dir: Path,
) -> dict[str, str]:
    row = {
        "output_dir": str(output_dir),
        "num_domains": str(num_domains),
        "min_size": "1",
        "max_cr": "0.1",
        "density_min": "1.0",
        "mean_density": "2.0",
        "delineation": delineation,
        "boundary_coil_fraction": "0.75",
        "energy_z": "",
        "modal_count_distance": "0.0",
    }
    row.update({field: "0.0" for field in CANDIDATE_GEOMETRY_FIELDS})
    return row


def test_fieldnames_include_new_reranker_features():
    for col in ("boundary_coil_fraction", "energy_z", "modal_count_distance"):
        assert col in FIELDNAMES, f"missing column: {col}"


def test_fieldnames_include_structural_candidate_features():
    for col in CANDIDATE_GEOMETRY_FIELDS:
        assert col in FIELDNAMES, f"missing column: {col}"


def test_score_candidates_passes_through_new_features(tmp_path):
    _write_pdb(tmp_path / "testchain.pdb", list(range(1, 21)))
    reference = {
        "testchain": _entry(n_residues=20, n_domains=1, chopping="1-20:1_1")
    }
    candidate = _candidate(
        delineation="0-19", num_domains=1, output_dir=tmp_path / "missing-output"
    )
    candidate["min_size"] = "20"
    candidate["energy_z"] = "-2.5"

    scored = _score_candidates(
        "testchain", [candidate], reference, chain_cache_dir=tmp_path
    )

    assert scored.rejections == []
    assert len(scored.rows) == 1
    assert scored.rows[0]["boundary_coil_fraction"] == "0.75"
    assert scored.rows[0]["energy_z"] == "-2.5"
    assert scored.rows[0]["modal_count_distance"] == "0.0"


def test_score_candidates_maps_author_reference_from_chain_cache(tmp_path):
    # The candidate uses its native 0-based residue indices, while CATH uses
    # author numbers.  The temporary SWORD output is intentionally absent: the
    # durable cache must be used to map the CATH reference before scoring.
    _write_pdb(tmp_path / "testchain.pdb", [10, 11])
    reference = {
        "testchain": _entry(n_residues=2, n_domains=2, chopping="10:1_1|11:2_2")
    }
    candidates = [
        _candidate(
            delineation="0 1",
            num_domains=2,
            output_dir=tmp_path / "missing-output",
        )
    ]

    scored = _score_candidates(
        "testchain", candidates, reference, chain_cache_dir=tmp_path
    )

    assert scored.rejections == []
    assert len(scored.rows) == 1
    assert scored.rows[0]["n_true_domains"] == 2
    assert scored.rows[0]["n_pred_domains"] == 2
    assert scored.rows[0]["ndo"] == 1.0


def test_reference_alias_uses_canonical_entry_id_for_chain_cache(tmp_path):
    _write_pdb(tmp_path / "testchain.pdb", [10, 11])
    entry = _entry(n_residues=2, n_domains=1, chopping="10-11:1_1")
    candidate = _candidate(
        delineation="0-1", num_domains=1, output_dir=tmp_path / "missing-output"
    )

    scored = _score_candidates(
        "TEST_A", [candidate], {"TEST_A": entry}, chain_cache_dir=tmp_path
    )

    assert scored.rejections == []
    assert len(scored.rows) == 1
    assert scored.rows[0]["chain_id"] == "testchain"


def test_reference_loader_rejects_alias_collisions(tmp_path):
    reference_path = tmp_path / "reference.csv"
    reference_path.write_text(
        "test,shared,A,1,unused,1,1\n"
        "other,shared,B,1,unused,1,1\n"
    )

    with pytest.raises(ValueError, match="reference alias collision"):
        _load_reference(str(reference_path))


def test_missing_output_dir_does_not_search_the_process_working_directory(
    tmp_path, monkeypatch
):
    accidental = tmp_path / "intermediate"
    accidental.mkdir()
    _write_pdb(accidental / "testchain.pdb", [10, 11])
    monkeypatch.chdir(tmp_path)
    entry = _entry(n_residues=2, n_domains=1, chopping="10-11:1_1")
    candidate = _candidate(delineation="0-1", num_domains=1, output_dir=tmp_path)
    candidate.pop("output_dir")

    scored = _score_candidates("testchain", [candidate], {"testchain": entry})

    assert scored.rows == []
    assert [record.code for record in scored.rejections] == [
        RejectionCode.MISSING_CHAIN_PDB
    ]


def test_score_candidates_rejects_missing_reference(tmp_path):
    scored = _score_candidates(
        "unknown",
        [_candidate(delineation="0", num_domains=1, output_dir=tmp_path)],
        {},
        chain_cache_dir=tmp_path,
    )

    assert scored.rows == []
    assert [record.code for record in scored.rejections] == [
        RejectionCode.MISSING_REFERENCE
    ]


def test_score_candidates_rejects_missing_chain_pdb(tmp_path):
    reference = {
        "testchain": _entry(n_residues=2, n_domains=1, chopping="10-11:1_1")
    }

    scored = _score_candidates(
        "testchain",
        [_candidate(delineation="0-1", num_domains=1, output_dir=tmp_path)],
        reference,
        chain_cache_dir=tmp_path,
    )

    assert scored.rows == []
    assert [record.code for record in scored.rejections] == [
        RejectionCode.MISSING_CHAIN_PDB
    ]


def test_score_candidates_rejects_reference_residue_count_mismatch(tmp_path):
    _write_pdb(tmp_path / "testchain.pdb", [10, 11])
    reference = {
        "testchain": _entry(n_residues=3, n_domains=1, chopping="10-11:1_1")
    }

    scored = _score_candidates(
        "testchain",
        [_candidate(delineation="0-1", num_domains=1, output_dir=tmp_path)],
        reference,
        chain_cache_dir=tmp_path,
    )

    assert scored.rows == []
    assert [record.code for record in scored.rejections] == [
        RejectionCode.RESIDUE_COUNT_MISMATCH
    ]


def test_score_candidates_rejects_reference_domain_count_mismatch(tmp_path):
    _write_pdb(tmp_path / "testchain.pdb", [10, 11])
    reference = {
        "testchain": _entry(n_residues=2, n_domains=3, chopping="10:1_1|11:2_2")
    }

    scored = _score_candidates(
        "testchain",
        [_candidate(delineation="0 1", num_domains=2, output_dir=tmp_path)],
        reference,
        chain_cache_dir=tmp_path,
    )

    assert scored.rows == []
    assert [record.code for record in scored.rejections] == [
        RejectionCode.TRUE_DOMAIN_COUNT_MISMATCH
    ]


def test_score_candidates_rejects_only_the_invalid_candidate(tmp_path):
    _write_pdb(tmp_path / "testchain.pdb", [10, 11, 12, 13])
    reference = {
        "testchain": _entry(n_residues=4, n_domains=2, chopping="10-11:1_1|12-13:2_2")
    }
    candidates = [
        _candidate(delineation="0-2 2-3", num_domains=2, output_dir=tmp_path),
        _candidate(delineation="0-1 2-3", num_domains=2, output_dir=tmp_path),
    ]

    scored = _score_candidates(
        "testchain", candidates, reference, chain_cache_dir=tmp_path
    )

    assert len(scored.rows) == 1
    assert scored.rows[0]["delineation"] == "0-1 2-3"
    assert [record.code for record in scored.rejections] == [
        RejectionCode.CANDIDATE_OVERLAP
    ]


def test_nonfinite_required_feature_rejects_candidate(tmp_path):
    _write_pdb(tmp_path / "testchain.pdb", [10, 11])
    reference = {
        "testchain": _entry(n_residues=2, n_domains=1, chopping="10-11:1_1")
    }
    candidate = _candidate(delineation="0-1", num_domains=1, output_dir=tmp_path)
    candidate["density_min"] = "nan"

    scored = _score_candidates(
        "testchain", [candidate], reference, chain_cache_dir=tmp_path
    )

    assert scored.rows == []
    assert [record.code for record in scored.rejections] == [
        RejectionCode.NONFINITE_FEATURE
    ]


def test_nonfinite_optional_energy_is_blank_but_does_not_reject_candidate(tmp_path):
    _write_pdb(tmp_path / "testchain.pdb", [10, 11])
    reference = {
        "testchain": _entry(n_residues=2, n_domains=1, chopping="10-11:1_1")
    }
    candidate = _candidate(delineation="0-1", num_domains=1, output_dir=tmp_path)
    candidate["energy_z"] = "nan"

    scored = _score_candidates(
        "testchain", [candidate], reference, chain_cache_dir=tmp_path
    )

    assert scored.rejections == []
    assert len(scored.rows) == 1
    assert scored.rows[0]["energy_z"] == ""


def test_write_rejections_is_deterministically_sorted(tmp_path):
    path = tmp_path / "rejections.csv"
    records = [
        RejectionRecord(
            chain_id="b",
            scope="candidate",
            code=RejectionCode.CANDIDATE_OVERLAP,
            detail="later",
            delineation="0-2 2-3",
        ),
        RejectionRecord(
            chain_id="a",
            scope="chain",
            code=RejectionCode.MISSING_CHAIN_PDB,
            detail="missing",
        ),
        RejectionRecord(
            chain_id="b",
            scope="candidate",
            code=RejectionCode.CANDIDATE_OUT_OF_RANGE,
            detail="first",
            delineation="0-5",
        ),
    ]

    write_rejections(path, records)

    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [(row["chain_id"], row["scope"], row["code"]) for row in rows] == [
        ("a", "chain", "missing_chain_pdb"),
        ("b", "candidate", "candidate_out_of_range"),
        ("b", "candidate", "candidate_overlap"),
    ]
    assert tuple(rows[0]) == ("chain_id", "scope", "code", "detail", "delineation")


def test_part_file_rejects_row_identity_that_disagrees_with_filename(tmp_path):
    part = tmp_path / "wrong-chain.csv"
    fieldnames = ["chain_id", "entry_id", *BASE_CANDIDATE_FIELDS, *CANDIDATE_GEOMETRY_FIELDS]
    with part.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        candidate = _candidate(delineation="0", num_domains=1, output_dir=tmp_path)
        candidate.pop("output_dir")
        writer.writerow(
            {
                "chain_id": "testchain",
                "entry_id": "testchain",
                **candidate,
            }
        )

    scored = _score_chain_file(part)

    assert scored.rows == []
    assert [record.code for record in scored.rejections] == [
        RejectionCode.SCHEMA_MISMATCH
    ]


def test_merged_dump_without_identity_column_writes_schema_rejection(tmp_path):
    dump_path = tmp_path / "dump.csv"
    dump_path.write_text("num_domains,delineation\n1,0\n")
    output_path = tmp_path / "training.csv"
    rejections_path = tmp_path / "rejections.csv"

    build_from_dump_csv(
        dump_path,
        reference={},
        output_path=output_path,
        rejections_path=rejections_path,
    )

    with rejections_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["code"] for row in rows] == ["schema_mismatch"]
    assert "chain_id" in rows[0]["detail"]
