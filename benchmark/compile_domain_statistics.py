#!/usr/bin/env python3
"""Compile domain shape, compactness, and inter-domain contact statistics over CATH.

Fits Gaussians to:
1. Domain principal radii of gyration ratios to ideal sphere: q1, q2, q3
2. Domain volume ratio to ideal sphere: vol_ratio
3. Domain C-alpha residue density: density
4. Inter-domain contact ratio: contact_q

Saves the fitted mean and std to benchmark/data/domain_statistics_gaussians.json.
Generates histograms of the distributions to verify if they are Gaussian.
"""
from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

from benchmark.datasets import load_dataset, strip_cath_labels
from benchmark.numbering import numbering_from_pdb, map_author_chopping

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
CHAINS = REPO / "benchmark/cache/chains"

# Constants for ideal protein domain sphere model
V0 = 141.0  # Average residue volume in Å^3
RHO_IDEAL = 1.0 / V0  # Ideal density in residues/Å^3
C_IDEAL = (3.0 * V0 / (4.0 * np.pi))**(1.0 / 3.0)  # Ideal sphere radius scaling factor
C_GYR_IDEAL = C_IDEAL / np.sqrt(5.0)  # Ideal principal radius of gyration scaling factor

# Integration constants for contact probability function
# p(d) = 1 / (1 + exp((d - D0)/delta)), D0=6.0, delta=1.5
# I3 = integral_0^inf w^3 p(w) dw = 780.69552
# I4 = integral_0^inf w^4 p(w) dw = 6137.0618
I3 = 780.69552
I4 = 6137.0618

K_MIN = (2.0 / 3.0) * (np.pi**2) * (RHO_IDEAL**2) * I4  # ~2.03097 Å^-1
K_MAX = (np.pi**2) * (RHO_IDEAL**2) * I3  # ~0.38755 Å^-2


def load_ca_coords(pdb_path: Path) -> list[np.ndarray]:
    """Parse CA atom coordinates from PDB file."""
    coords = []
    with pdb_path.open() as f:
        for line in f:
            if line.startswith("ATOM  ") and line[17:20].strip() in {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
                "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
                "TYR", "VAL"
            }:
                if line[12:16].strip() == "CA":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append(np.array([x, y, z]))
    return coords


def solve_eigenvalues(S: np.ndarray) -> list[float]:
    """Analytical eigenvalue solver for a symmetric 3x3 matrix S."""
    m = np.trace(S)
    n = S[0,0]*S[1,1] + S[1,1]*S[2,2] + S[0,0]*S[2,2] - S[0,1]**2 - S[1,2]**2 - S[0,2]**2
    d = np.linalg.det(S)

    p = m**2 - 3.0 * n
    q = 2.0 * m**3 - 9.0 * m * n + 27.0 * d

    if p <= 1e-9:
        val = max(0.0, m / 3.0)
        return [val, val, val]

    arg = q / (2.0 * p**1.5)
    arg = np.clip(arg, -1.0, 1.0)
    theta = np.arccos(arg) / 3.0

    l1 = (m + 2.0 * np.sqrt(p) * np.cos(theta)) / 3.0
    l2 = (m + 2.0 * np.sqrt(p) * np.cos(theta + 2.0 * np.pi / 3.0)) / 3.0
    l3 = (m + 2.0 * np.sqrt(p) * np.cos(theta + 4.0 * np.pi / 3.0)) / 3.0

    # Return sorted descending and non-negative
    vals = sorted([max(0.0, l1), max(0.0, l2), max(0.0, l3)], reverse=True)
    return vals


def parse_domains(chopping: str) -> list[list[tuple[int, int]]]:
    """Parse sequential 0-based chopping into domain segment lists."""
    domains = []
    # delimiter is comma
    for raw_domain in chopping.strip().split(","):
        if not raw_domain.strip():
            continue
        segments = []
        # segments within a domain are separated by underscore
        for raw_seg in raw_domain.split("_"):
            if not raw_seg.strip():
                continue
            if "-" in raw_seg:
                start, end = raw_seg.split("-")
                segments.append((int(start), int(end)))
            else:
                idx = int(raw_seg)
                segments.append((idx, idx))
        domains.append(segments)
    return domains


def get_domain_indices(segments: list[tuple[int, int]]) -> list[int]:
    """Get flat list of residue indices for a domain's segments."""
    indices = []
    for start, end in segments:
        indices.extend(range(start, end + 1))
    return indices


def solve_x(f: float) -> float:
    """Solve the cubic equation for domain volume fraction cut parameter x."""
    val = np.clip(2.0 * f - 1.0, -1.0, 1.0)
    theta = np.arccos(val) / 3.0
    x = 2.0 * np.cos(theta + 4.0 * np.pi / 3.0)
    return x


def compute_contacts(ca_coords: list[np.ndarray], idx_a: list[int], idx_b: list[int]) -> float:
    """Calculate actual contact count sum between two sets of C-alpha coordinates."""
    contacts = 0.0
    for i in idx_a:
        for j in idx_b:
            if i >= len(ca_coords) or j >= len(ca_coords):
                continue
            d = np.linalg.norm(ca_coords[i] - ca_coords[j])
            p = 1.0 / (1.0 + np.exp((d - 6.0) / 1.5))
            contacts += p
    return contacts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", default="cath663")
    ap.add_argument("--out-json", type=Path, default=REPO / "benchmark/data/domain_statistics_gaussians.json")
    ap.add_argument("--out-plot", type=Path, default=REPO / "benchmark/data/domain_statistics_histograms.png")
    args = ap.parse_args()

    entries = load_dataset(args.dataset)
    log.info("Processing %d entries in dataset %s", len(entries), args.dataset)

    q1_vals = []
    q2_vals = []
    q3_vals = []
    vol_ratio_vals = []
    density_vals = []
    contact_q_vals = []

    success_count = 0
    for entry in entries:
        pdb_path = CHAINS / f"{entry.entry_id}.pdb"
        if not pdb_path.exists():
            continue

        try:
            numbering = numbering_from_pdb(pdb_path, chain_id=entry.chain_id)
            true_chop = map_author_chopping(strip_cath_labels(entry.chopping), numbering, chain_id=entry.chain_id)
            ca_coords = load_ca_coords(pdb_path)
            if not ca_coords:
                continue
        except Exception:
            continue

        domains = parse_domains(true_chop)
        if len(domains) == 0:
            continue

        success_count += 1

        # Domain level features
        domain_res_indices = []
        for dom_segments in domains:
            idx = get_domain_indices(dom_segments)
            domain_res_indices.append(idx)

            # Select coordinates
            dom_coords = [ca_coords[i] for i in idx if i < len(ca_coords)]
            N = len(dom_coords)
            if N < 4:
                continue

            # Centroid
            centroid = np.mean(dom_coords, axis=0)

            # Covariance matrix S
            S = np.zeros((3, 3))
            for v in dom_coords:
                diff = v - centroid
                S += np.outer(diff, diff)
            S /= N

            # Eigenvalues
            eigenvals = solve_eigenvalues(S)
            r1, r2, r3 = np.sqrt(eigenvals[0]), np.sqrt(eigenvals[1]), np.sqrt(eigenvals[2])

            # Ideal principal radius of gyration
            r_ideal = C_GYR_IDEAL * (N**(1.0 / 3.0))

            # Clamp to prevent divide-by-zero
            r1_clamped = max(r1, 1e-3)
            r2_clamped = max(r2, 1e-3)
            r3_clamped = max(r3, 1e-3)

            q1 = r1_clamped / r_ideal
            q2 = r2_clamped / r_ideal
            q3 = r3_clamped / r_ideal
            vol_ratio = (r1_clamped * r2_clamped * r3_clamped) / (r_ideal**3)

            # Equivalent Volume V_real = 20 * sqrt(5) * pi * r1 * r2 * r3 / 3
            V_real = (20.0 * np.sqrt(5.0) / 3.0) * np.pi * r1_clamped * r2_clamped * r3_clamped
            density = N / V_real
            rel_density = density / RHO_IDEAL

            q1_vals.append(q1)
            q2_vals.append(q2)
            q3_vals.append(q3)
            vol_ratio_vals.append(vol_ratio)
            density_vals.append(rel_density)

        # Pair level features
        for a in range(len(domains)):
            for b in range(a + 1, len(domains)):
                idx_a = domain_res_indices[a]
                idx_b = domain_res_indices[b]
                N_A = len([i for i in idx_a if i < len(ca_coords)])
                N_B = len([i for i in idx_b if i < len(ca_coords)])
                if N_A < 4 or N_B < 4:
                    continue

                C_real = compute_contacts(ca_coords, idx_a, idx_b)
                if C_real < 1.0:
                    continue  # Not in contact

                # R_A and R_B ideal
                R_A = C_IDEAL * (N_A**(1.0 / 3.0))
                R_B = C_IDEAL * (N_B**(1.0 / 3.0))
                C_min = K_MIN * (R_A * R_B) / (R_A + R_B)

                R_AB = C_IDEAL * ((N_A + N_B)**(1.0 / 3.0))
                f = min(N_A, N_B) / (N_A + N_B)
                x = solve_x(f)
                C_max = K_MAX * (R_AB**2) * (1.0 - x**2)

                denom = C_max - C_min
                if abs(denom) > 1e-5:
                    contact_q = (C_real - C_min) / denom
                    # Clamp to [0, 1]
                    contact_q = np.clip(contact_q, 0.0, 1.0)
                    contact_q_vals.append(contact_q)

    log.info("Successfully analyzed %d/%d entries", success_count, len(entries))
    log.info("Collected %d domain-level samples and %d domain-pair samples", len(q1_vals), len(contact_q_vals))

    # Compute parameters and print
    features = {
        "q1": q1_vals,
        "q2": q2_vals,
        "q3": q3_vals,
        "vol_ratio": vol_ratio_vals,
        "rel_density": density_vals,
        "contact_q": contact_q_vals
    }

    fitted_params = {}
    print("\n--- Feature Distributions and Gaussian Fits ---")
    for name, vals in features.items():
        mean = float(np.mean(vals))
        std = float(np.std(vals))
        skew = float(stats.skew(vals))
        kurt = float(stats.kurtosis(vals))
        fitted_params[name] = {"mean": mean, "std": std}
        print(f"Feature: {name:12} | Mean: {mean:8.4f} | Std: {std:8.4f} | Skewness: {skew:8.4f} | Kurtosis: {kurt:8.4f}")

    # Save to JSON
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    with args.out_json.open("w") as f:
        json.dump(fitted_params, f, indent=2)
    log.info("Saved fitted parameters to %s", args.out_json)

    # Plot histograms
    fig, axes = plt.subplots(3, 2, figsize=(12, 15))
    axes = axes.flatten()
    for i, (name, vals) in enumerate(features.items()):
        ax = axes[i]
        sns.histplot(vals, kde=True, ax=ax, color="teal", stat="density")
        # Overlay fitted normal distribution
        mean = fitted_params[name]["mean"]
        std = fitted_params[name]["std"]
        x_plot = np.linspace(min(vals), max(vals), 100)
        y_plot = stats.norm.pdf(x_plot, mean, std)
        ax.plot(x_plot, y_plot, color="red", linestyle="--", linewidth=2, label="Normal Fit")
        ax.set_title(f"{name} (mean={mean:.3f}, std={std:.3f})")
        ax.legend()

    plt.tight_layout()
    args.out_plot.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out_plot, dpi=300)
    log.info("Saved feature histograms to %s", args.out_plot)


if __name__ == "__main__":
    main()
