from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

from benchmark.numbering import numbering_from_pdb
from benchmark.structures import download_pdb, extract_single_chain


def test_extract_single_chain_preserves_header_for_legacy_tools(tmp_path):
    raw = tmp_path / "raw.pdb"
    raw.write_text(
        "HEADER    TEST STRUCTURE                           01-JAN-00   1ABC\n"
        "TITLE     EXAMPLE\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N  \n"
        "ATOM      2  N   GLY B   1       1.000   1.000   1.000  1.00 10.00           N  \n"
        "END\n"
    )
    output = tmp_path / "chain.pdb"

    extract_single_chain(raw, "A", output)

    lines = output.read_text().splitlines()
    assert lines[0].startswith("HEADER")
    assert lines[1].startswith("TITLE")
    assert sum(line.startswith("ATOM") for line in lines) == 1
    assert " A   1 " in lines[2]


def test_extract_single_chain_renames_single_available_chain_to_requested_chain(tmp_path):
    raw = tmp_path / "raw.pdb"
    raw.write_text(
        "HEADER    TEST STRUCTURE                           01-JAN-00   1ABC\n"
        "TITLE     EXAMPLE\n"
        "ATOM      1  N   ALA B   1       0.000   0.000   0.000  1.00 10.00           N  \n"
        "ATOM      2  CA  ALA B   1       1.000   1.000   1.000  1.00 10.00           C  \n"
        "END\n"
    )
    output = tmp_path / "chain.pdb"

    extract_single_chain(raw, "A", output)

    atom_lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
    assert len(atom_lines) == 2
    assert {line[21] for line in atom_lines} == {"A"}


def test_extract_single_chain_resolves_alias_from_atom_chains_not_ligand_chains(tmp_path):
    raw = tmp_path / "raw.pdb"
    raw.write_text(
        "HEADER    TEST STRUCTURE                           01-JAN-00   1ABC\n"
        "TITLE     EXAMPLE\n"
        "HETATM    1  C1  LIG A 900       9.000   9.000   9.000  1.00 10.00           C  \n"
        "ATOM      2  N   ALA B   1       0.000   0.000   0.000  1.00 10.00           N  \n"
        "END\n"
    )
    output = tmp_path / "chain.pdb"

    extract_single_chain(raw, "A", output)

    lines = output.read_text().splitlines()
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    assert len(atom_lines) == 1
    assert atom_lines[0][21] == "A"
    assert not any(line.startswith("HETATM") for line in lines)


def test_numbering_uses_only_standard_protein_atom_residues(tmp_path):
    raw = tmp_path / "raw.pdb"
    raw.write_text(
        "HEADER    TEST STRUCTURE                           01-JAN-00   1ABC\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N  \n"
        "ATOM      2  CA  MSE A   2       1.000   1.000   1.000  1.00 10.00           C  \n"
        "HETATM    3  C1  LIG A 900       9.000   9.000   9.000  1.00 10.00           C  \n"
        "ATOM      4  N   GLY A   3       2.000   2.000   2.000  1.00 10.00           N  \n"
        "END\n"
    )

    numbering = numbering_from_pdb(raw, chain_id="A")

    assert [residue.auth_seq for residue in numbering.residues] == ["1", "3"]


def test_extract_single_chain_writes_only_standard_protein_atoms(tmp_path):
    raw = tmp_path / "raw.pdb"
    raw.write_text(
        "HEADER    TEST STRUCTURE                           01-JAN-00   1ABC\n"
        "TITLE     EXAMPLE\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N  \n"
        "ATOM      2  CA  MSE A   2       1.000   1.000   1.000  1.00 10.00           C  \n"
        "HETATM    3  C1  LIG A 900       9.000   9.000   9.000  1.00 10.00           C  \n"
        "ATOM      4  N   GLY A   3       2.000   2.000   2.000  1.00 10.00           N  \n"
        "END\n"
    )
    output = tmp_path / "chain.pdb"

    extract_single_chain(raw, "A", output)

    atom_lines = [line for line in output.read_text().splitlines() if line.startswith("ATOM")]
    assert len(atom_lines) == 2
    assert {line[17:20].strip() for line in atom_lines} == {"ALA", "GLY"}
    assert not any(line.startswith("HETATM") for line in output.read_text().splitlines())


@contextmanager
def _fake_response(body: bytes):
    class _Response:
        def read(self_inner) -> bytes:
            return body

        def __enter__(self_inner):
            return self_inner

        def __exit__(self_inner, *exc_info):
            return False

    yield _Response()


def test_download_pdb_skips_network_when_already_cached(tmp_path):
    cache_dir = tmp_path / "raw"
    cache_dir.mkdir()
    existing = cache_dir / "1abc.pdb"
    existing.write_text("cached content")

    with patch("urllib.request.urlopen") as mock_urlopen:
        result = download_pdb("1ABC", cache_dir)

    mock_urlopen.assert_not_called()
    assert result == existing
    assert result.read_text() == "cached content"


def test_download_pdb_writes_atomically_leaving_no_temp_file(tmp_path):
    cache_dir = tmp_path / "raw"

    with patch("urllib.request.urlopen", return_value=_fake_response(b"ATOM data")):
        result = download_pdb("1ABC", cache_dir)

    assert result == cache_dir / "1abc.pdb"
    assert result.read_bytes() == b"ATOM data"
    assert list(cache_dir.glob("*.tmp")) == []
