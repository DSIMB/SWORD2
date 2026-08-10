from __future__ import annotations

import os
import threading
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from benchmark.datasets import CathEntry
from benchmark.numbering import StructureNumbering, is_standard_protein_atom_line, numbering_from_pdb


@dataclass(frozen=True)
class CanonicalStructure:
    entry: CathEntry
    pdb_path: Path
    numbering: StructureNumbering


def download_pdb(pdb_id: str, cache_dir: Path) -> Path:
    """Download a PDB file, writing it atomically so concurrent callers for the
    same pdb_id (shared across multiple chain entries) never observe a partially
    written file."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{pdb_id.lower()}.pdb"
    if path.exists():
        return path
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    with urllib.request.urlopen(url, timeout=60) as response:
        data = response.read()
    tmp_path = cache_dir / f".{pdb_id.lower()}.{os.getpid()}.{threading.get_ident()}.pdb.tmp"
    tmp_path.write_bytes(data)
    os.replace(tmp_path, path)
    return path


def _atom_chains(input_pdb: Path) -> set[str]:
    atom_chains: set[str] = set()
    with input_pdb.open() as handle:
        for line in handle:
            if is_standard_protein_atom_line(line):
                atom_chains.add(line[21].strip())
    return atom_chains


def _replace_pdb_chain_id(line: str, chain_id: str) -> str:
    if len(chain_id) != 1:
        return line
    return f"{line[:21]}{chain_id}{line[22:]}"


def extract_single_chain(input_pdb: Path, chain_id: str, output_pdb: Path) -> Path:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    available_chains = _atom_chains(input_pdb)
    source_chain_id = chain_id
    if chain_id not in available_chains:
        if len(available_chains) == 1:
            source_chain_id = next(iter(available_chains))
        else:
            available = ", ".join(repr(chain) for chain in sorted(available_chains)) or "none"
            raise ValueError(
                f"Chain {chain_id!r} was not found in {input_pdb}; available atom chains: {available}"
            )

    header_records = {
        "HEADER",
        "TITLE ",
        "COMPND",
        "SOURCE",
        "KEYWDS",
        "EXPDTA",
        "AUTHOR",
        "REVDAT",
        "JRNL  ",
        "REMARK",
        "SEQRES",
        "DBREF ",
    }
    wrote_atom = False
    with input_pdb.open() as src, output_pdb.open("w") as dst:
        for line in src:
            record = line[:6]
            if not wrote_atom and record in header_records:
                dst.write(line)
                continue
            if record == "ATOM  ":
                if not is_standard_protein_atom_line(line):
                    continue
                if line[21].strip() != source_chain_id:
                    continue
                if source_chain_id != chain_id:
                    line = _replace_pdb_chain_id(line, chain_id)
                dst.write(line)
                wrote_atom = True
            elif record == "TER   " and wrote_atom:
                if source_chain_id != chain_id:
                    line = _replace_pdb_chain_id(line, chain_id)
                dst.write(line)
                break
        dst.write("END\n")
    return output_pdb


def canonicalize_entry(entry: CathEntry, cache_dir: Path, download: bool = True) -> CanonicalStructure:
    raw_dir = cache_dir / "raw"
    chain_dir = cache_dir / "chains"
    raw_pdb = raw_dir / f"{entry.pdb_id.lower()}.pdb"
    if download:
        raw_pdb = download_pdb(entry.pdb_id, raw_dir)
    elif not raw_pdb.exists():
        raise FileNotFoundError(raw_pdb)

    chain_pdb = chain_dir / f"{entry.entry_id}.pdb"
    extract_single_chain(raw_pdb, entry.chain_id, chain_pdb)
    return CanonicalStructure(
        entry=entry,
        pdb_path=chain_pdb,
        numbering=numbering_from_pdb(chain_pdb, chain_id=entry.chain_id),
    )
