"""Fail-closed reference mapping and candidate-partition validation."""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Literal

from benchmark.datasets import CathEntry, strip_cath_labels
from benchmark.numbering import map_author_chopping, numbering_from_pdb, split_domains


_CANDIDATE_SEGMENT = re.compile(r"[0-9]+(?:-[0-9]+)?", flags=re.ASCII)


class RejectionCode(str, Enum):
    MISSING_REFERENCE = "missing_reference"
    MISSING_CHAIN_PDB = "missing_chain_pdb"
    AUTHOR_MAPPING_FAILED = "author_mapping_failed"
    RESIDUE_COUNT_MISMATCH = "residue_count_mismatch"
    TRUE_DOMAIN_COUNT_MISMATCH = "true_domain_count_mismatch"
    CANDIDATE_PARSE_FAILED = "candidate_parse_failed"
    CANDIDATE_OUT_OF_RANGE = "candidate_out_of_range"
    CANDIDATE_OVERLAP = "candidate_overlap"
    CANDIDATE_INCOMPLETE_COVERAGE = "candidate_incomplete_coverage"
    CANDIDATE_DOMAIN_COUNT_MISMATCH = "candidate_domain_count_mismatch"
    SCORING_FAILED = "scoring_failed"
    NONFINITE_FEATURE = "nonfinite_feature"
    SCHEMA_MISMATCH = "schema_mismatch"


class IntegrityError(ValueError):
    """A rejected reference or candidate with a stable machine-readable code."""

    def __init__(self, code: RejectionCode, detail: str):
        self.code = code
        self.detail = detail
        super().__init__(f"{code.value}: {detail}")


@dataclass(frozen=True)
class RejectionRecord:
    chain_id: str
    scope: Literal["chain", "candidate"]
    code: RejectionCode
    detail: str
    delineation: str | None = None


@dataclass(frozen=True)
class CanonicalReference:
    chopping: str
    n_residues: int
    n_domains: int


@dataclass(frozen=True)
class ValidatedPartition:
    canonical_delineation: str
    domains: tuple[tuple[tuple[int, int], ...], ...]
    residue_to_domain: tuple[int, ...]


def map_cath_reference(entry: CathEntry, pdb_path: Path) -> CanonicalReference:
    """Map a CATH author-numbered reference onto the cached clean chain."""
    raw_domains = entry.chopping.strip().split("|")
    malformed_source = not entry.chopping.strip() or any(
        not raw_domain.strip()
        or not raw_domain.split(":", 1)[0].strip()
        or any(
            not segment.strip()
            for segment in raw_domain.split(":", 1)[0].split("_")
        )
        for raw_domain in raw_domains
    )
    if malformed_source:
        raise IntegrityError(
            RejectionCode.AUTHOR_MAPPING_FAILED,
            "CATH chopping contains an empty domain or segment",
        )

    try:
        numbering = numbering_from_pdb(pdb_path, chain_id=entry.chain_id)
    except Exception as exc:
        raise IntegrityError(
            RejectionCode.AUTHOR_MAPPING_FAILED,
            f"could not read chain numbering ({type(exc).__name__})",
        ) from exc

    if numbering.n_residues != entry.n_residues:
        raise IntegrityError(
            RejectionCode.RESIDUE_COUNT_MISMATCH,
            f"reference declares {entry.n_residues} residues but PDB has "
            f"{numbering.n_residues}",
        )

    try:
        chopping = map_author_chopping(
            strip_cath_labels(entry.chopping),
            numbering,
            chain_id=entry.chain_id,
        )
    except Exception as exc:
        raise IntegrityError(
            RejectionCode.AUTHOR_MAPPING_FAILED,
            f"could not map author residue ranges ({type(exc).__name__})",
        ) from exc

    mapped_domains = split_domains(chopping)
    if len(mapped_domains) != entry.n_domains:
        raise IntegrityError(
            RejectionCode.TRUE_DOMAIN_COUNT_MISMATCH,
            f"reference declares {entry.n_domains} domains but mapped chopping has "
            f"{len(mapped_domains)}",
        )

    occupied: set[int] = set()
    for domain in mapped_domains:
        for segment in domain:
            pieces = segment.split("-")
            try:
                start = int(pieces[0])
                end = int(pieces[1]) if len(pieces) == 2 else start
            except (IndexError, ValueError) as exc:
                raise IntegrityError(
                    RejectionCode.AUTHOR_MAPPING_FAILED,
                    "mapped reference contains a malformed segment",
                ) from exc
            residues = set(range(start, end + 1))
            if occupied.intersection(residues):
                raise IntegrityError(
                    RejectionCode.AUTHOR_MAPPING_FAILED,
                    "mapped reference domains overlap",
                )
            occupied.update(residues)

    return CanonicalReference(
        chopping=chopping,
        n_residues=numbering.n_residues,
        n_domains=len(mapped_domains),
    )


def _reject(code: RejectionCode, detail: str) -> None:
    raise IntegrityError(code, detail)


def validate_partition(
    delineation: str,
    n_residues: int,
    declared_domains: int,
) -> ValidatedPartition:
    """Parse and canonicalize a complete zero-based inclusive partition."""
    if n_residues <= 0 or not delineation.strip():
        _reject(RejectionCode.CANDIDATE_PARSE_FAILED, "empty candidate partition")

    parsed: list[list[tuple[int, int]]] = []
    for raw_domain in delineation.strip().split():
        raw_segments = raw_domain.split(";")
        if not raw_segments or any(not segment for segment in raw_segments):
            _reject(
                RejectionCode.CANDIDATE_PARSE_FAILED,
                "candidate contains an empty segment",
            )

        segments: list[tuple[int, int]] = []
        for raw_segment in raw_segments:
            if _CANDIDATE_SEGMENT.fullmatch(raw_segment) is None:
                _reject(
                    RejectionCode.CANDIDATE_PARSE_FAILED,
                    f"malformed segment {raw_segment!r}",
                )
            pieces = raw_segment.split("-")
            try:
                start = int(pieces[0])
                end = int(pieces[1]) if len(pieces) == 2 else start
            except ValueError:
                _reject(
                    RejectionCode.CANDIDATE_PARSE_FAILED,
                    f"malformed segment {raw_segment!r}",
                )
            if start < 0 or end < start:
                _reject(
                    RejectionCode.CANDIDATE_PARSE_FAILED,
                    f"invalid segment bounds {raw_segment!r}",
                )
            if end >= n_residues:
                _reject(
                    RejectionCode.CANDIDATE_OUT_OF_RANGE,
                    f"segment {start}-{end} exceeds chain length {n_residues}",
                )
            segments.append((start, end))
        parsed.append(sorted(segments))

    if len(parsed) != declared_domains:
        _reject(
            RejectionCode.CANDIDATE_DOMAIN_COUNT_MISMATCH,
            f"row declares {declared_domains} domains but delineation has {len(parsed)}",
        )

    parsed.sort(key=lambda domain: min(start for start, _ in domain))
    owner = [-1] * n_residues
    for domain_index, segments in enumerate(parsed):
        for start, end in segments:
            for residue in range(start, end + 1):
                if owner[residue] != -1:
                    _reject(
                        RejectionCode.CANDIDATE_OVERLAP,
                        f"residue {residue} occurs more than once",
                    )
                owner[residue] = domain_index

    if any(value == -1 for value in owner):
        _reject(
            RejectionCode.CANDIDATE_INCOMPLETE_COVERAGE,
            "candidate does not cover every chain residue",
        )

    coalesced: list[list[tuple[int, int]]] = []
    for domain_index in range(len(parsed)):
        residues = [
            residue for residue, owner_index in enumerate(owner) if owner_index == domain_index
        ]
        runs: list[tuple[int, int]] = []
        start = end = residues[0]
        for residue in residues[1:]:
            if residue == end + 1:
                end = residue
            else:
                runs.append((start, end))
                start = end = residue
        runs.append((start, end))
        coalesced.append(runs)

    def render(segment: tuple[int, int]) -> str:
        start, end = segment
        return str(start) if start == end else f"{start}-{end}"

    frozen = tuple(tuple(domain) for domain in coalesced)
    canonical = " ".join(";".join(map(render, domain)) for domain in frozen)
    return ValidatedPartition(canonical, frozen, tuple(owner))
