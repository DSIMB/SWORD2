from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


DEFAULT_MERIZO_DATASET_DIR = Path("/home/chili/cretin/PROJECTS/Merizo/datasets/merizo_domains")


@dataclass(frozen=True)
class CathEntry:
    pdb_id: str
    chain_id: str
    entry_id: str
    n_domains: int
    n_residues: int
    chopping: str
    dataset: str


DATASETS = {
    "cath663": DEFAULT_MERIZO_DATASET_DIR / "CATH-663.csv",
    "cath17287": DEFAULT_MERIZO_DATASET_DIR / "CATH-17287.csv",
    "afdb1195": DEFAULT_MERIZO_DATASET_DIR / "AFDB-1195.csv",
    "afdb7502": DEFAULT_MERIZO_DATASET_DIR / "AFDB-7502.csv",
}


def dataset_path(name: str) -> Path:
    try:
        return DATASETS[name.lower()]
    except KeyError as exc:
        known = ", ".join(sorted(DATASETS))
        raise ValueError(f"Unknown dataset {name!r}; known datasets: {known}") from exc


def parse_cath_domain_string(chopping: str) -> list[list[str]]:
    """Return domains as author-numbered segments, dropping CATH label suffixes."""
    domains: list[list[str]] = []
    for raw_domain in chopping.strip().split("|"):
        raw_domain = raw_domain.strip()
        if not raw_domain:
            continue
        domain_part = raw_domain.split(":", 1)[0]
        segments = [segment.strip() for segment in domain_part.split("_") if segment.strip()]
        if segments:
            domains.append(segments)
    return domains


def strip_cath_labels(chopping: str) -> str:
    return "|".join("_".join(domain) for domain in parse_cath_domain_string(chopping))


def cath_family_labels(chopping: str) -> tuple[str, ...]:
    """Return known CATH family suffixes in domain order, retaining multiplicity."""
    labels: list[str] = []
    for raw_domain in chopping.split("|"):
        if ":" not in raw_domain:
            continue
        label = raw_domain.rsplit(":", 1)[1].strip()
        if label and label != "999_999":
            labels.append(label)
    return tuple(labels)


def cath_family_combination(chopping: str) -> tuple[str, ...]:
    """Return the exact sorted multiset of known CATH family suffixes."""
    return tuple(sorted(cath_family_labels(chopping)))


def read_merizo_csv(path: Path, dataset: str, limit: int | None = None) -> list[CathEntry]:
    """Read Merizo benchmark CSVs with their headerless seven-column layout."""
    entries: list[CathEntry] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 7:
                raise ValueError(f"Expected at least 7 columns in {path}, got {len(row)}: {row!r}")
            entries.append(
                CathEntry(
                    pdb_id=row[0].strip().lower(),
                    entry_id=row[1].strip(),
                    chain_id=row[2].strip(),
                    n_domains=int(row[3]),
                    n_residues=int(row[5]),
                    chopping=row[6].strip(),
                    dataset=dataset,
                )
            )
            if limit is not None and len(entries) >= limit:
                break
    return entries


def load_dataset(name: str, limit: int | None = None) -> list[CathEntry]:
    return read_merizo_csv(dataset_path(name), dataset=name.lower(), limit=limit)
