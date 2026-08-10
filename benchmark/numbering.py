from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


_RESIDUE_TOKEN = re.compile(r"^\s*(-?\d+)([A-Za-z]?)\s*$")
_SEGMENT_TOKEN = re.compile(r"^\s*(-?\d+[A-Za-z]?)(?:-(-?\d+[A-Za-z]?))?\s*$")

STANDARD_PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


def is_standard_protein_atom_line(line: str) -> bool:
    return line.startswith("ATOM  ") and line[17:20].strip() in STANDARD_PROTEIN_RESIDUES


@dataclass(frozen=True, order=True)
class ResidueKey:
    chain_id: str
    auth_seq: str
    insertion_code: str = ""

    @classmethod
    def from_pdb_line(cls, line: str) -> "ResidueKey":
        return cls(
            chain_id=line[21].strip(),
            auth_seq=line[22:26].strip(),
            insertion_code=line[26].strip(),
        )

    @classmethod
    def from_token(cls, token: str, chain_id: str = "") -> "ResidueKey":
        match = _RESIDUE_TOKEN.match(token)
        if not match:
            raise ValueError(f"Unsupported residue token {token!r}")
        return cls(chain_id=chain_id, auth_seq=match.group(1), insertion_code=match.group(2))


@dataclass(frozen=True)
class StructureNumbering:
    residues: list[ResidueKey]

    def __post_init__(self) -> None:
        if not self.residues:
            raise ValueError("Structure numbering needs at least one residue")

    @property
    def n_residues(self) -> int:
        return len(self.residues)

    def index_of(self, key: ResidueKey) -> int:
        exact = {residue: idx for idx, residue in enumerate(self.residues)}
        if key in exact:
            return exact[key]

        chainless = [
            idx
            for idx, residue in enumerate(self.residues)
            if residue.auth_seq == key.auth_seq and residue.insertion_code == key.insertion_code
        ]
        if len(chainless) == 1:
            return chainless[0]
        if not chainless:
            raise KeyError(f"Residue {key.auth_seq}{key.insertion_code} is not present")
        raise KeyError(f"Residue {key.auth_seq}{key.insertion_code} is ambiguous across chains")

    def indices_in_author_range(self, start: ResidueKey, end: ResidueKey) -> tuple[int, int]:
        """Return first/last present residue indices covered by an author-numbered range."""
        if start.insertion_code or end.insertion_code:
            raise KeyError(
                f"Residue range {start.auth_seq}{start.insertion_code}-"
                f"{end.auth_seq}{end.insertion_code} is not present"
            )
        lo, hi = sorted((int(start.auth_seq), int(end.auth_seq)))
        chain = start.chain_id or end.chain_id
        matches = [
            idx
            for idx, residue in enumerate(self.residues)
            if (not chain or residue.chain_id == chain)
            and not residue.insertion_code
            and lo <= int(residue.auth_seq) <= hi
        ]
        if not matches:
            raise KeyError(f"No residues from {start.auth_seq}-{end.auth_seq} are present")
        return min(matches), max(matches)


def residues_from_pdb(path: Path, chain_id: str | None = None) -> list[ResidueKey]:
    residues: list[ResidueKey] = []
    seen: set[ResidueKey] = set()
    with path.open() as handle:
        for line in handle:
            if not is_standard_protein_atom_line(line):
                continue
            residue = ResidueKey.from_pdb_line(line)
            if chain_id is not None and residue.chain_id != chain_id:
                continue
            if residue not in seen:
                seen.add(residue)
                residues.append(residue)
    return residues


def numbering_from_pdb(path: Path, chain_id: str | None = None) -> StructureNumbering:
    return StructureNumbering(residues_from_pdb(path, chain_id=chain_id))


def split_domains(chopping: str) -> list[list[str]]:
    delimiter = "|" if "|" in chopping else ","
    domains: list[list[str]] = []
    for domain in chopping.strip().split(delimiter):
        domain = domain.strip()
        if not domain:
            continue
        domains.append([segment.strip() for segment in domain.split("_") if segment.strip()])
    return domains


def format_domains(domains: list[list[str]]) -> str:
    return ",".join("_".join(domain) for domain in domains)


def _split_segment(segment: str) -> tuple[str, str]:
    match = _SEGMENT_TOKEN.match(segment)
    if not match:
        raise ValueError(f"Unsupported residue segment {segment!r}")
    start = match.group(1)
    end = match.group(2) or start
    return start.strip(), end.strip()


def _format_zero_based_segment(start_idx: int, end_idx: int) -> str:
    lo, hi = sorted((start_idx, end_idx))
    return f"{lo}-{hi}" if lo != hi else str(lo)


def infer_n_res(choppings: list[str]) -> int:
    max_end = -1
    for chopping in choppings:
        for domain in split_domains(chopping):
            for segment in domain:
                _, end = _split_segment(segment)
                max_end = max(max_end, int(end))
    return max_end + 1


def map_one_based_chopping(chopping: str) -> str:
    mapped: list[list[str]] = []
    for domain in split_domains(chopping):
        mapped_domain: list[str] = []
        for segment in domain:
            start, end = _split_segment(segment)
            mapped_domain.append(_format_zero_based_segment(int(start) - 1, int(end) - 1))
        mapped.append(mapped_domain)
    return format_domains(mapped)


def map_zero_based_chopping(chopping: str) -> str:
    return format_domains(split_domains(chopping))


def map_author_chopping(
    chopping: str,
    numbering: StructureNumbering,
    chain_id: str | None = None,
) -> str:
    mapped: list[list[str]] = []
    chain = chain_id or ""
    for domain in split_domains(chopping):
        mapped_domain: list[str] = []
        for segment in domain:
            start, end = _split_segment(segment)
            start_key = ResidueKey.from_token(start, chain_id=chain)
            end_key = ResidueKey.from_token(end, chain_id=chain)
            try:
                start_idx = numbering.index_of(start_key)
                end_idx = numbering.index_of(end_key)
            except KeyError:
                start_idx, end_idx = numbering.indices_in_author_range(start_key, end_key)
            mapped_domain.append(_format_zero_based_segment(start_idx, end_idx))
        mapped.append(mapped_domain)
    return format_domains(mapped)


def read_sword2_residue_mapping(path: Path, chain_id: str | None = None) -> dict[int, ResidueKey]:
    """Read SWORD2 residue_mapping.txt into renumbered 1-based index -> author residue."""
    mapping: dict[int, ResidueKey] = {}
    chain = chain_id or ""
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.upper().startswith("ORIGINAL"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            author, renum = parts[0], int(parts[1])
            mapping[renum] = ResidueKey.from_token(author, chain_id=chain)
    return mapping


def map_sword2_chopping(
    chopping: str,
    numbering: StructureNumbering | None,
    residue_mapping: dict[int, ResidueKey] | None = None,
) -> str:
    if numbering is None:
        return map_one_based_chopping(chopping)

    mapped: list[list[str]] = []
    for domain in split_domains(chopping):
        mapped_domain: list[str] = []
        for segment in domain:
            start, end = _split_segment(segment)
            start_renum, end_renum = int(start), int(end)
            if residue_mapping:
                start_idx = numbering.index_of(residue_mapping[start_renum])
                end_idx = numbering.index_of(residue_mapping[end_renum])
            else:
                start_idx = start_renum - 1
                end_idx = end_renum - 1
            mapped_domain.append(_format_zero_based_segment(start_idx, end_idx))
        mapped.append(mapped_domain)
    return format_domains(mapped)
