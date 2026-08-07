"""Independent NumPy reference formulas for structural ranker features."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal, Mapping, Sequence

import numpy as np

from benchmark.candidate_geometry import _domain_shape
from benchmark.factorized_ranker.integrity import ValidatedPartition
from benchmark.factorized_ranker.schema import (
    BOUNDARY_LOCAL_FEATURES,
    BOUNDARY_MEASURES,
    DISCONTINUITY_FEATURES,
    DISCONTINUITY_MEASURES,
    DOMAIN_CONDITIONED_FEATURES,
    DOMAIN_SIDE_MEASURES,
    GLOBAL_FEATURES,
    MAX_COUNT_HISTOGRAM_BIN,
)


@dataclass(frozen=True)
class DsspHydrogenBond:
    partner: int
    energy_kcal: float


@dataclass(frozen=True)
class DsspResidue:
    index: int
    ss: Literal["helix", "strand", "coil"]
    bridge_partners: tuple[int, ...]
    sheet_label: str
    hydrogen_bonds: tuple[DsspHydrogenBond, ...]
    kappa: float
    alpha: float


@dataclass(frozen=True)
class PeelingSummary:
    n_levels: int
    finest_pu_count: int


def _coordinates(value: np.ndarray) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if (
        result.ndim != 2
        or result.shape[1:] != (3,)
        or len(result) == 0
        or not np.isfinite(result).all()
    ):
        raise ValueError("invalid coordinates")
    return result


def _contacts(value: np.ndarray, n_residues: int) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if (
        result.shape != (n_residues, n_residues)
        or not np.isfinite(result).all()
        or np.any(result < 0.0)
        or np.any(result > 1.0)
        or not np.allclose(result, result.T, rtol=0.0, atol=1e-12)
    ):
        raise ValueError("invalid contact matrix")
    return result


def _dssp_records(
    value: Sequence[DsspResidue], n_residues: int
) -> tuple[DsspResidue, ...]:
    if len(value) != n_residues:
        raise ValueError("invalid DSSP evidence")
    by_index = {residue.index: residue for residue in value}
    if set(by_index) != set(range(n_residues)):
        raise ValueError("invalid DSSP evidence")
    result = tuple(by_index[index] for index in range(n_residues))
    for residue in result:
        scalars = (residue.kappa, residue.alpha)
        partners = (*residue.bridge_partners, *(bond.partner for bond in residue.hydrogen_bonds))
        if (
            not all(math.isfinite(float(item)) for item in scalars)
            or not all(0 <= partner < n_residues for partner in partners)
            or not all(math.isfinite(bond.energy_kcal) for bond in residue.hydrogen_bonds)
        ):
            raise ValueError("invalid DSSP evidence")
    return result


def _partition_length(partition: ValidatedPartition, n_residues: int) -> None:
    if len(partition.residue_to_domain) != n_residues:
        raise ValueError("partition length does not match structure")


def _ss_class(code: str) -> str:
    normalized = code.strip()
    if normalized in {"helix", "H", "G", "I"}:
        return "helix"
    if normalized in {"strand", "E", "B"}:
        return "strand"
    return "coil"


def _ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator != 0.0 else 0.0


def _summary(values: Sequence[float]) -> tuple[float, float, float]:
    if not values:
        return (0.0, 0.0, 0.0)
    array = np.asarray(values, dtype=float)
    return (float(array.min()), float(array.mean()), float(array.max()))


def _ordered_finite(names: Sequence[str], values: Mapping[str, float]) -> dict[str, float]:
    if set(values) != set(names):
        raise ValueError("feature schema mismatch")
    result = {name: float(values[name]) for name in names}
    if not all(math.isfinite(value) for value in result.values()):
        raise ValueError("non-finite structural feature")
    return result


def compute_global_features(
    coordinates: np.ndarray,
    contacts: np.ndarray,
    dssp: Sequence[DsspResidue],
    peeling: PeelingSummary,
    candidate_counts: Sequence[int],
) -> dict[str, float]:
    coordinates = _coordinates(coordinates)
    n_residues = len(coordinates)
    contacts = _contacts(contacts, n_residues)
    dssp = _dssp_records(dssp, n_residues)
    counts = tuple(int(value) for value in candidate_counts)
    if (
        any(value <= 0 for value in counts)
        or any(float(raw) != value for raw, value in zip(candidate_counts, counts))
        or peeling.n_levels < 0
        or peeling.finest_pu_count < 0
    ):
        raise ValueError("invalid global feature evidence")

    centered = coordinates - coordinates.mean(axis=0)
    radius_of_gyration = math.sqrt(float(np.mean(np.sum(centered * centered, axis=1))))
    tensor = centered.T @ centered / n_residues
    eigenvalues = np.maximum(np.linalg.eigvalsh(tensor)[::-1], 0.0)
    first = float(eigenvalues[0])

    nonlocal_values = [
        float(contacts[left, right])
        for left in range(n_residues)
        for right in range(left + 8, n_residues)
    ]
    pair_weights = np.triu(contacts, k=1)
    weight_sum = float(pair_weights.sum())
    weighted_separation = sum(
        float(contacts[left, right]) * (right - left)
        for left in range(n_residues)
        for right in range(left + 1, n_residues)
    )
    classes = tuple(_ss_class(residue.ss) for residue in dssp)

    def block_count(target: str) -> int:
        return sum(
            current == target and (index == 0 or classes[index - 1] != target)
            for index, current in enumerate(classes)
        )

    histogram = {value: 0 for value in range(1, MAX_COUNT_HISTOGRAM_BIN + 1)}
    overflow = 0
    for count in counts:
        if count <= MAX_COUNT_HISTOGRAM_BIN:
            histogram[count] += 1
        else:
            overflow += 1
    frequencies = {count: counts.count(count) for count in set(counts)}
    modal_count = min(frequencies, key=lambda count: (-frequencies[count], count)) if counts else 0

    values: dict[str, float] = {
        "chain_n_residues": n_residues,
        "chain_rg_normalized": radius_of_gyration / n_residues ** (1.0 / 3.0),
        "chain_inertia_ratio_21": math.sqrt(float(eigenvalues[1]) / first) if first else 0.0,
        "chain_inertia_ratio_31": math.sqrt(float(eigenvalues[2]) / first) if first else 0.0,
        "chain_nonlocal_contact_density": float(np.mean(nonlocal_values)) if nonlocal_values else 0.0,
        "chain_contact_order": weighted_separation / weight_sum / n_residues if weight_sum else 0.0,
        "chain_helix_fraction": classes.count("helix") / n_residues,
        "chain_strand_fraction": classes.count("strand") / n_residues,
        "chain_coil_fraction": classes.count("coil") / n_residues,
        "chain_helix_blocks": block_count("helix"),
        "chain_strand_blocks": block_count("strand"),
        "chain_peeling_levels": peeling.n_levels,
        "chain_finest_pus": peeling.finest_pu_count,
        "chain_candidate_total": len(counts),
        "chain_available_count_total": len(set(counts)),
        **{
            f"chain_count_hist_{value}": histogram[value] / len(counts) if counts else 0.0
            for value in histogram
        },
        "chain_count_hist_21_plus": overflow / len(counts) if counts else 0.0,
        "chain_modal_count": modal_count,
    }
    return _ordered_finite(GLOBAL_FEATURES, values)


def _domain_indices(partition: ValidatedPartition) -> list[np.ndarray]:
    return [
        np.asarray(
            [index for start, end in segments for index in range(start, end + 1)],
            dtype=np.intp,
        )
        for segments in partition.domains
    ]


def compute_domain_conditioned_features(
    partition: ValidatedPartition,
    coordinates: np.ndarray,
    contacts: np.ndarray,
) -> dict[str, float]:
    coordinates = _coordinates(coordinates)
    n_residues = len(coordinates)
    contacts = _contacts(contacts, n_residues)
    _partition_length(partition, n_residues)
    domains = _domain_indices(partition)
    if not domains or any(len(domain) == 0 for domain in domains):
        raise ValueError("invalid partition")

    measures: list[dict[str, float]] = []
    conductances: list[float] = []
    all_indices = np.arange(n_residues, dtype=np.intp)
    for indices in domains:
        q1, q2, q3, _volume_ratio, relative_density = _domain_shape(coordinates[indices])
        internal_pairs = [
            (int(left), int(right))
            for offset, left in enumerate(indices)
            for right in indices[offset + 1 :]
        ]
        nonlocal_contacts = [
            float(contacts[left, right])
            for left, right in internal_pairs
            if right - left >= 8
        ]
        internal_mass = sum(float(contacts[left, right]) for left, right in internal_pairs)
        external = all_indices[~np.isin(all_indices, indices)]
        external_mass = float(contacts[np.ix_(indices, external)].sum())
        incident_mass = internal_mass + external_mass
        weighted_separation = sum(
            float(contacts[left, right]) * (right - left)
            for left, right in internal_pairs
        )
        internal_weight = sum(float(contacts[left, right]) for left, right in internal_pairs)
        internal_fraction = internal_mass / incident_mass if incident_mass else 0.0
        denominator = 2.0 * internal_mass + external_mass
        conductance = external_mass / denominator if denominator else 0.0
        measures.append(
            {
                "size_fraction": len(indices) / n_residues,
                "q1": q1,
                "q2": q2,
                "q3": q3,
                "relative_density": relative_density,
                "internal_contact_density": float(np.mean(nonlocal_contacts)) if nonlocal_contacts else 0.0,
                "contact_order": weighted_separation / internal_weight / n_residues if internal_weight else 0.0,
                "internal_contact_fraction": internal_fraction,
            }
        )
        conductances.append(conductance)

    sizes = [len(domain) for domain in domains]
    smallest = measures[int(np.argmin(sizes))]
    largest = measures[int(np.argmax(sizes))]
    values = {
        **{f"smallest_{name}": smallest[name] for name in DOMAIN_SIDE_MEASURES},
        **{f"largest_{name}": largest[name] for name in DOMAIN_SIDE_MEASURES},
        "smallest_to_largest_density_ratio": _ratio(
            smallest["relative_density"], largest["relative_density"]
        ),
        "smallest_to_largest_internal_contact_density_ratio": _ratio(
            smallest["internal_contact_density"], largest["internal_contact_density"]
        ),
        "smallest_to_largest_q1_ratio": _ratio(smallest["q1"], largest["q1"]),
        "smallest_to_largest_q2_ratio": _ratio(smallest["q2"], largest["q2"]),
        "smallest_to_largest_q3_ratio": _ratio(smallest["q3"], largest["q3"]),
    }
    for stat, value in zip(
        ("min", "mean", "max"),
        _summary([measure["internal_contact_fraction"] for measure in measures]),
    ):
        values[f"domain_internal_contact_fraction_{stat}"] = value
    for stat, value in zip(("min", "mean", "max"), _summary(conductances)):
        values[f"domain_conductance_{stat}"] = value
    return _ordered_finite(DOMAIN_CONDITIONED_FEATURES, values)


def aggregate_measurements(
    measurements: Sequence[Mapping[str, float]],
) -> dict[str, float]:
    """Aggregate supplied boundary measurements in frozen schema order."""
    present = {
        name for name in BOUNDARY_MEASURES
        if not measurements or any(name in measurement for measurement in measurements)
    }
    result: dict[str, float] = {}
    for name in BOUNDARY_MEASURES:
        if name not in present:
            continue
        if any(name not in measurement for measurement in measurements):
            raise ValueError(f"missing boundary measurement {name}")
        summary = _summary([float(measurement[name]) for measurement in measurements])
        for stat, value in zip(("min", "mean", "max"), summary):
            result[f"boundary_{name}_{stat}"] = value
    if not all(math.isfinite(value) for value in result.values()):
        raise ValueError("non-finite boundary measurement")
    return result


def _crosses(left: int, right: int, boundary: int) -> bool:
    return (left <= boundary < right) or (right <= boundary < left)


def _boundary_measurement(
    boundary: int,
    contacts: np.ndarray,
    dssp: tuple[DsspResidue, ...],
) -> dict[str, float]:
    n_residues = len(dssp)
    classes = tuple(_ss_class(residue.ss) for residue in dssp)
    same_class = classes[boundary] == classes[boundary + 1]
    inside = classes[boundary] if same_class else "coil"
    terminus_distance = 0.0
    if inside in {"helix", "strand"} and same_class:
        start = boundary
        while start > 0 and classes[start - 1] == inside:
            start -= 1
        end = boundary + 1
        while end + 1 < n_residues and classes[end + 1] == inside:
            end += 1
        terminus_distance = float(min(boundary - start + 1, end - boundary))

    crossing_bonds: dict[tuple[int, int], float] = {}
    for residue in dssp:
        for bond in residue.hydrogen_bonds:
            if _crosses(residue.index, bond.partner, boundary):
                crossing_bonds.setdefault((residue.index, bond.partner), bond.energy_kcal)

    crossing_bridges = {
        tuple(sorted((residue.index, partner)))
        for residue in dssp
        for partner in residue.bridge_partners
        if residue.index != partner and _crosses(residue.index, partner, boundary)
    }
    sheet_pairs = {
        (left, right)
        for left in range(boundary + 1)
        for right in range(boundary + 1, n_residues)
        if dssp[left].sheet_label.strip()
        and dssp[left].sheet_label == dssp[right].sheet_label
    }

    def insulation(window: int) -> float:
        left = np.arange(max(0, boundary - window + 1), boundary + 1)
        right = np.arange(boundary + 1, min(n_residues, boundary + window + 1))
        return 1.0 - float(contacts[np.ix_(left, right)].mean())

    long_range = [
        float(contacts[left, right])
        for left in range(boundary + 1)
        for right in range(boundary + 1, n_residues)
        if right - left >= 8
    ]
    alpha_left = dssp[boundary].alpha
    alpha_right = dssp[boundary + 1].alpha
    if alpha_left == 360.0 or alpha_right == 360.0:
        dihedral_change = 0.0
    else:
        dihedral_change = abs((alpha_left - alpha_right + 180.0) % 360.0 - 180.0) / 180.0

    return {
        "inside_helix": float(inside == "helix" and same_class),
        "inside_strand": float(inside == "strand" and same_class),
        "inside_coil": float(inside == "coil"),
        "sse_terminus_distance": terminus_distance,
        "hbond_count": float(len(crossing_bonds)),
        "hbond_energy_kcal": float(sum(crossing_bonds.values())),
        "bridge_count": float(len(crossing_bridges)),
        "sheet_link_count": float(len(sheet_pairs)),
        "insulation_w8": insulation(8),
        "insulation_w16": insulation(16),
        "insulation_w32": insulation(32),
        "long_range_contact_density": float(np.mean(long_range)) if long_range else 0.0,
        "bend_change": abs(dssp[boundary].kappa - dssp[boundary + 1].kappa) / 180.0,
        "virtual_dihedral_change": dihedral_change,
    }


def compute_boundary_local_features(
    partition: ValidatedPartition,
    coordinates: np.ndarray,
    contacts: np.ndarray,
    dssp: Sequence[DsspResidue],
) -> dict[str, float]:
    coordinates = _coordinates(coordinates)
    n_residues = len(coordinates)
    contacts = _contacts(contacts, n_residues)
    dssp = _dssp_records(dssp, n_residues)
    _partition_length(partition, n_residues)
    owners = partition.residue_to_domain
    boundaries = [
        index for index in range(n_residues - 1) if owners[index] != owners[index + 1]
    ]
    measurements = [
        _boundary_measurement(boundary, contacts, dssp) for boundary in boundaries
    ]
    return _ordered_finite(BOUNDARY_LOCAL_FEATURES, aggregate_measurements(measurements))


def _separation_bin(separation: int) -> int:
    if separation <= 15:
        return 0
    if separation <= 31:
        return 1
    if separation <= 63:
        return 2
    return 3


def compute_discontinuity_features(
    partition: ValidatedPartition,
    contacts: np.ndarray,
    dssp: Sequence[DsspResidue],
) -> dict[str, float]:
    n_residues = len(partition.residue_to_domain)
    contacts = _contacts(contacts, n_residues)
    dssp = _dssp_records(dssp, n_residues)
    discontinuous_domains = {
        domain_index
        for domain_index, segments in enumerate(partition.domains)
        if len(segments) > 1
    }
    if not discontinuous_domains:
        return {name: 0.0 for name in DISCONTINUITY_FEATURES}

    domains = _domain_indices(partition)
    segment_of = [-1] * n_residues
    segments: list[tuple[int, int, np.ndarray]] = []
    for domain_index, domain_segments in enumerate(partition.domains):
        for start, end in domain_segments:
            segment_index = len(segments)
            indices = np.arange(start, end + 1, dtype=np.intp)
            segments.append((domain_index, segment_index, indices))
            for index in indices:
                segment_of[int(index)] = segment_index

    sheet_edges = {
        tuple(sorted((residue.index, partner)))
        for residue in dssp
        for partner in residue.bridge_partners
        if residue.index != partner
        and partition.residue_to_domain[residue.index]
        == partition.residue_to_domain[partner]
        and segment_of[residue.index] != segment_of[partner]
    }

    measurements: list[dict[str, float]] = []
    for domain_index, segment_index, segment in segments:
        if domain_index not in discontinuous_domains:
            continue
        segment_set = set(int(index) for index in segment)
        same_domain_other = np.asarray(
            [index for index in domains[domain_index] if index not in segment_set],
            dtype=np.intp,
        )
        affinity = float(contacts[np.ix_(segment, same_domain_other)].mean())
        competitor_affinities = [
            float(contacts[np.ix_(segment, domain)].mean())
            for other_index, domain in enumerate(domains)
            if other_index != domain_index
        ]
        margin = affinity - max(competitor_affinities) if competitor_affinities else 0.0

        internal_mass = 0.0
        all_mass = 0.0
        bin_masses = np.zeros(4, dtype=float)
        same_domain_set = set(int(index) for index in same_domain_other)
        for left in segment:
            for right in range(n_residues):
                if right in segment_set:
                    continue
                separation = abs(int(left) - right)
                if separation < 8:
                    continue
                mass = float(contacts[int(left), right])
                all_mass += mass
                if right in same_domain_set:
                    internal_mass += mass
                    bin_masses[_separation_bin(separation)] += mass
        capture = internal_mass / all_mass if all_mass else 0.0
        if internal_mass:
            probabilities = bin_masses[bin_masses > 0.0] / internal_mass
            entropy = -float(np.sum(probabilities * np.log(probabilities))) / math.log(4.0)
        else:
            entropy = 0.0
        sheet_links = sum(
            1
            for left, right in sheet_edges
            if segment_of[left] == segment_index or segment_of[right] == segment_index
        )
        measurements.append(
            {
                "same_domain_affinity": affinity,
                "affinity_margin": margin,
                "long_range_internal_capture": capture,
                "interface_span_entropy": entropy,
                "same_domain_sheet_links": float(sheet_links),
            }
        )

    values: dict[str, float] = {"has_discontinuity": 1.0}
    for measure in DISCONTINUITY_MEASURES:
        for stat, value in zip(
            ("min", "mean", "max"),
            _summary([measurement[measure] for measurement in measurements]),
        ):
            values[f"segment_{measure}_{stat}"] = value
    return _ordered_finite(DISCONTINUITY_FEATURES, values)
