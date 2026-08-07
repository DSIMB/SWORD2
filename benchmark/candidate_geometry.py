"""Fast structural features for a proposed SWORD domain partition.

The functions mirror the analytical domain statistics used during the original
selection investigation, but operate on one SWORD candidate delineation.  They
are deliberately independent of CATH labels and can therefore be computed at
training or inference time from only C-alpha coordinates and candidate bounds.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np

from benchmark.factorized_ranker.integrity import validate_partition
from benchmark.numbering import is_standard_protein_atom_line


# Average residue volume in A^3, matching compile_domain_statistics.py.
V0 = 141.0
RHO_IDEAL = 1.0 / V0
C_IDEAL = (3.0 * V0 / (4.0 * math.pi)) ** (1.0 / 3.0)
C_GYR_IDEAL = C_IDEAL / math.sqrt(5.0)
I3 = 780.69552
I4 = 6137.0618
K_MIN = (2.0 / 3.0) * (math.pi**2) * (RHO_IDEAL**2) * I4
K_MAX = (math.pi**2) * (RHO_IDEAL**2) * I3
_SEGMENT = re.compile(r"[0-9]+(?:-[0-9]+)?", flags=re.ASCII)


CANDIDATE_GEOMETRY_FIELDS = [
    "domain_q1_mean",
    "domain_q2_mean",
    "domain_q3_mean",
    "domain_q3_min",
    "domain_vol_ratio_mean",
    "domain_density_mean",
    "domain_density_min",
    "contact_q_mean",
    "contact_q_max",
    "n_segments",
    "n_discontinuous",
    "size_balance",
    "largest_domain_fraction",
    "min_segment_size",
    "mean_segment_size",
]


def load_ca_coordinates(pdb_path: Path) -> np.ndarray:
    """Load standard-residue C-alpha coordinates in sequential chain order."""
    coordinates: list[tuple[float, float, float]] = []
    with pdb_path.open() as handle:
        for line in handle:
            if not is_standard_protein_atom_line(line) or line[12:16].strip() != "CA":
                continue
            coordinates.append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            )
    return np.asarray(coordinates, dtype=float)


def parse_delineation(delineation: str) -> list[list[tuple[int, int]]]:
    """Parse SWORD's 0-based ``space``/``;`` candidate delineation syntax."""
    if not delineation.strip():
        return []
    domains: list[list[tuple[int, int]]] = []
    for domain in delineation.strip().split():
        segments: list[tuple[int, int]] = []
        for raw_segment in domain.split(";"):
            if not raw_segment or _SEGMENT.fullmatch(raw_segment) is None:
                return []
            try:
                pieces = raw_segment.split("-")
                if len(pieces) not in (1, 2):
                    return []
                start = int(pieces[0])
                end = int(pieces[1]) if len(pieces) == 2 else start
            except (ValueError, IndexError):
                return []
            if start < 0 or end < start:
                return []
            segments.append((start, end))
        if segments:
            domains.append(segments)
    return domains


def contact_probability_matrix(coordinates: np.ndarray) -> np.ndarray:
    """Return the SWORD logistic C-alpha contact-probability matrix."""
    coordinates = np.asarray(coordinates, dtype=float)
    if (
        coordinates.ndim != 2
        or coordinates.shape[1:] != (3,)
        or len(coordinates) == 0
        or not np.isfinite(coordinates).all()
    ):
        raise ValueError("invalid coordinate array")
    deltas = coordinates[:, None, :] - coordinates[None, :, :]
    distance = np.sqrt(np.einsum("ijk,ijk->ij", deltas, deltas))
    return 1.0 / (1.0 + np.exp(np.clip((distance - 6.0) / 1.5, -50.0, 50.0)))


def _segment_indices(
    segments: list[tuple[int, int]], n_residues: int
) -> np.ndarray | None:
    indices: list[int] = []
    for start, end in segments:
        if start < 0 or end >= n_residues:
            return None
        indices.extend(range(start, end + 1))
    return np.asarray(indices, dtype=np.intp)


def _domain_shape(points: np.ndarray) -> tuple[float, float, float, float, float]:
    """Return q1/q2/q3, ellipsoid volume ratio, and relative C-alpha density."""
    n_points = len(points)
    if n_points < 4:
        return (0.0, 0.0, 0.0, 0.0, 0.0)
    centered = points - points.mean(axis=0)
    tensor = centered.T @ centered / n_points
    eigenvalues = np.maximum(np.linalg.eigvalsh(tensor)[::-1], 0.0)
    radii = np.sqrt(eigenvalues)
    ideal_radius = C_GYR_IDEAL * n_points ** (1.0 / 3.0)
    q1, q2, q3 = radii / max(ideal_radius, 1e-12)
    volume_ratio = float(q1 * q2 * q3)
    clamped_radii = np.maximum(radii, 1e-3)
    ellipsoid_volume = (20.0 * math.sqrt(5.0) / 3.0) * math.pi * float(
        np.prod(clamped_radii)
    )
    relative_density = (n_points / ellipsoid_volume) / RHO_IDEAL
    return (float(q1), float(q2), float(q3), volume_ratio, relative_density)


def _contact_q(
    probability: float, size_a: int, size_b: int
) -> float:
    """Normalize an inter-domain contact mass by the sphere-model bounds."""
    if size_a < 4 or size_b < 4:
        return 0.0
    radius_a = C_IDEAL * size_a ** (1.0 / 3.0)
    radius_b = C_IDEAL * size_b ** (1.0 / 3.0)
    min_contacts = K_MIN * (radius_a * radius_b) / (radius_a + radius_b)
    radius_ab = C_IDEAL * (size_a + size_b) ** (1.0 / 3.0)
    fraction = min(size_a, size_b) / (size_a + size_b)
    angle = math.acos(max(-1.0, min(1.0, 2.0 * fraction - 1.0))) / 3.0
    x_value = 2.0 * math.cos(angle + 4.0 * math.pi / 3.0)
    max_contacts = K_MAX * radius_ab**2 * (1.0 - x_value**2)
    denominator = max_contacts - min_contacts
    if abs(denominator) <= 1e-5:
        return 0.0
    return float(np.clip((probability - min_contacts) / denominator, 0.0, 1.0))


def candidate_geometry_features(
    coordinates: np.ndarray,
    delineation: str,
    contact_matrix: np.ndarray | None = None,
) -> dict[str, float]:
    """Compute shape, contact, and fragmentation descriptors for one candidate."""
    coordinates = np.asarray(coordinates, dtype=float)
    if (
        coordinates.ndim != 2
        or coordinates.shape[1:] != (3,)
        or len(coordinates) == 0
        or not np.isfinite(coordinates).all()
    ):
        raise ValueError("invalid candidate partition")

    domains = parse_delineation(delineation)
    try:
        validated = validate_partition(
            delineation,
            n_residues=len(coordinates),
            declared_domains=len(domains),
        )
    except ValueError as exc:
        raise ValueError("invalid candidate partition") from exc
    domains = [list(domain) for domain in validated.domains]

    indices = [_segment_indices(segments, len(coordinates)) for segments in domains]
    if any(value is None or len(value) == 0 for value in indices):
        raise ValueError("invalid candidate partition")
    domain_indices = [value for value in indices if value is not None]
    sizes = [len(value) for value in domain_indices]
    segment_sizes = [end - start + 1 for segments in domains for start, end in segments]
    shapes = [_domain_shape(coordinates[value]) for value in domain_indices]
    q1, q2, q3, volume_ratio, density = (np.asarray(shapes, dtype=float).T)

    if contact_matrix is None:
        contact_matrix = contact_probability_matrix(coordinates)
    else:
        contact_matrix = np.asarray(contact_matrix, dtype=float)
        if (
            contact_matrix.shape != (len(coordinates), len(coordinates))
            or not np.isfinite(contact_matrix).all()
            or np.any(contact_matrix < 0.0)
            or np.any(contact_matrix > 1.0)
            or not np.allclose(
                contact_matrix,
                contact_matrix.T,
                rtol=0.0,
                atol=1e-12,
            )
        ):
            raise ValueError("invalid candidate partition")
    contact_q_values: list[float] = []
    for left in range(len(domain_indices)):
        for right in range(left + 1, len(domain_indices)):
            probability = float(
                contact_matrix[np.ix_(domain_indices[left], domain_indices[right])].sum()
            )
            contact_q_values.append(
                _contact_q(probability, sizes[left], sizes[right])
            )

    total_size = sum(sizes)
    features = {
        "domain_q1_mean": float(q1.mean()),
        "domain_q2_mean": float(q2.mean()),
        "domain_q3_mean": float(q3.mean()),
        "domain_q3_min": float(q3.min()),
        "domain_vol_ratio_mean": float(volume_ratio.mean()),
        "domain_density_mean": float(density.mean()),
        "domain_density_min": float(density.min()),
        "contact_q_mean": float(np.mean(contact_q_values)) if contact_q_values else 0.0,
        "contact_q_max": float(np.max(contact_q_values)) if contact_q_values else 0.0,
        "n_segments": float(len(segment_sizes)),
        "n_discontinuous": float(len(segment_sizes) - len(domains)),
        "size_balance": float(min(sizes) / max(sizes)),
        "largest_domain_fraction": float(max(sizes) / total_size),
        "min_segment_size": float(min(segment_sizes)),
        "mean_segment_size": float(np.mean(segment_sizes)),
    }
    if not all(math.isfinite(value) for value in features.values()):
        raise ValueError("invalid candidate partition")
    return features
