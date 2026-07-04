"""Fit reference Gaussians for the analytical "ideal sphere" geometry criteria.

Mirrors benchmark/fit_count_calibration.py: fits on CATH-17287, evaluates
(without leaking into the fit) on held-out CATH-663, and writes the fitted
constants that Rust embeds via sword2-lib/src/sword/geometry_metrics.rs. See
that module's doc comment for the full mathematical derivation (kappa^2 shape
anisotropy, ideal-sphere density scaling law, and the tangent/half-sphere
interface-fraction bound). The geometry functions below are line-for-line
mirrors of the Rust ones so the fitted constants mean the same thing on both
sides.
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
from scipy import stats

from benchmark.datasets import CathEntry, load_dataset, strip_cath_labels
from benchmark.numbering import (
    is_standard_protein_atom_line,
    map_author_chopping,
    numbering_from_pdb,
    split_domains,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger(__name__)

REPO = Path(__file__).resolve().parents[1]
CHAINS = REPO / "benchmark/cache/chains"

# Must match CONTACT_CUTOFF_ANGSTROM in sword2-lib/src/sword/geometry_metrics.rs.
CONTACT_CUTOFF_ANGSTROM = 8.0
MIN_DOMAIN_SIZE = 4
SHAPIRO_MAX_SAMPLE = 5000


def load_ca_coords(pdb_path: Path) -> np.ndarray:
    """Parse Cα coordinates from a PDB file, in file order (matches the Rust
    pipeline's ca_coords, both built from the same cleaned-chain ATOM records)."""
    coords: list[list[float]] = []
    with pdb_path.open() as f:
        for line in f:
            if is_standard_protein_atom_line(line) and line[12:16].strip() == "CA":
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(coords)


def gyration_tensor(coords: np.ndarray, idx: list[int]) -> np.ndarray:
    pts = coords[idx]
    centered = pts - pts.mean(axis=0)
    return (centered.T @ centered) / len(idx)


def sym3x3_eigenvalues(s: np.ndarray) -> np.ndarray:
    """Mirror of geometry_metrics.rs::sym3x3_eigenvalues (Cardano's trig solution)."""
    p1 = s[0, 1] ** 2 + s[0, 2] ** 2 + s[1, 2] ** 2
    if p1 <= 1e-12:
        return np.sort(np.array([s[0, 0], s[1, 1], s[2, 2]]))[::-1]

    q = np.trace(s) / 3.0
    p2 = (s[0, 0] - q) ** 2 + (s[1, 1] - q) ** 2 + (s[2, 2] - q) ** 2 + 2.0 * p1
    p = np.sqrt(p2 / 6.0)
    b = (s - q * np.eye(3)) / p
    det_b = np.linalg.det(b)
    r = np.clip(det_b / 2.0, -1.0, 1.0)
    phi = np.arccos(r) / 3.0

    eig1 = q + 2.0 * p * np.cos(phi)
    eig3 = q + 2.0 * p * np.cos(phi + 2.0 * np.pi / 3.0)
    eig2 = 3.0 * q - eig1 - eig3
    return np.sort(np.array([eig1, eig2, eig3]))[::-1]


def relative_shape_anisotropy(eig: np.ndarray) -> float:
    """kappa^2: mirror of geometry_metrics.rs::relative_shape_anisotropy."""
    total = eig.sum()
    if total <= 0:
        return 0.0
    pair_sum = eig[0] * eig[1] + eig[1] * eig[2] + eig[2] * eig[0]
    return float(np.clip(1.0 - 3.0 * pair_sum / (total * total), 0.0, 1.0))


def effective_sphere_radius(rg: float) -> float:
    return rg * np.sqrt(5.0 / 3.0)


def ca_number_density(n: int, r_eff: float) -> float:
    if r_eff <= 0:
        return 0.0
    v = (4.0 / 3.0) * np.pi * r_eff**3
    return n / v


def interdomain_contacts(coords: np.ndarray, idx_a: list[int], idx_b: list[int]) -> int:
    a = coords[idx_a]
    b = coords[idx_b]
    d2 = ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)
    return int((d2 <= CONTACT_CUTOFF_ANGSTROM**2).sum())


def parse_zero_based_domains(chopping: str) -> list[list[int]]:
    """0-based chopping (comma-separated domains, underscore-separated
    segments — the format `map_author_chopping` returns) into flat index lists."""
    domains: list[list[int]] = []
    for domain in split_domains(chopping):
        idx: list[int] = []
        for segment in domain:
            if "-" in segment:
                start, end = segment.split("-")
                idx.extend(range(int(start), int(end) + 1))
            else:
                idx.append(int(segment))
        domains.append(idx)
    return domains


class Samples:
    """Accumulated raw observations across the training dataset."""

    def __init__(self) -> None:
        self.kappa2: list[float] = []
        self.density_n: list[int] = []
        self.density_rho: list[float] = []
        self.r_eff_n: list[int] = []
        self.r_eff: list[float] = []
        self.interface_contacts: list[int] = []
        self.interface_r_a: list[float] = []
        self.interface_r_b: list[float] = []
        self.n_chains_used = 0


def collect_samples(entries: list[CathEntry]) -> Samples:
    samples = Samples()

    for entry in entries:
        pdb_path = CHAINS / f"{entry.entry_id}.pdb"
        if not pdb_path.exists():
            continue
        try:
            numbering = numbering_from_pdb(pdb_path, chain_id=entry.chain_id)
            chopping0 = map_author_chopping(
                strip_cath_labels(entry.chopping), numbering, chain_id=entry.chain_id
            )
            coords = load_ca_coords(pdb_path)
        except Exception:
            continue
        if coords.size == 0:
            continue

        domains = parse_zero_based_domains(chopping0)
        domains = [[i for i in dom if i < len(coords)] for dom in domains]
        if not domains:
            continue

        r_effs: list[float | None] = []
        for idx in domains:
            if len(idx) < MIN_DOMAIN_SIZE:
                r_effs.append(None)
                continue
            eig = sym3x3_eigenvalues(gyration_tensor(coords, idx))
            kappa2 = relative_shape_anisotropy(eig)
            rg = np.sqrt(max(float(eig.sum()), 0.0))
            r_eff = effective_sphere_radius(rg)
            density = ca_number_density(len(idx), r_eff)

            samples.kappa2.append(kappa2)
            samples.r_eff_n.append(len(idx))
            samples.r_eff.append(r_eff)
            samples.density_n.append(len(idx))
            samples.density_rho.append(density)
            r_effs.append(r_eff)

        for i in range(len(domains) - 1):
            r_a, r_b = r_effs[i], r_effs[i + 1]
            if r_a is None or r_b is None:
                continue
            contacts = interdomain_contacts(coords, domains[i], domains[i + 1])
            samples.interface_contacts.append(contacts)
            samples.interface_r_a.append(r_a)
            samples.interface_r_b.append(r_b)

        samples.n_chains_used += 1

    return samples


def fit_power_law(ns: np.ndarray, rs: np.ndarray) -> tuple[float, float]:
    """Least-squares log-log fit of R_ideal(N) = a * N^b."""
    mask = (ns > 0) & (rs > 0)
    log_n = np.log(ns[mask].astype(float))
    log_r = np.log(rs[mask].astype(float))
    design = np.vstack([np.ones_like(log_n), log_n]).T
    coef, *_ = np.linalg.lstsq(design, log_r, rcond=None)
    log_a, b = coef
    return float(np.exp(log_a)), float(b)


def fit_gamma(contacts: np.ndarray, r_a: np.ndarray, r_b: np.ndarray) -> float:
    """Surface contact density gamma: ratio-of-sums estimator for
    contacts ~= gamma * pi * min(r_a, r_b)^2 (regression through the origin)."""
    denom = np.pi * np.minimum(r_a, r_b) ** 2
    total_denom = denom.sum()
    if total_denom <= 0:
        return 0.0
    return float(contacts.sum() / total_denom)


def shapiro_p(values: np.ndarray, seed: int = 0) -> float:
    if len(values) < 3:
        return 0.0
    sample = values
    if len(values) > SHAPIRO_MAX_SAMPLE:
        rng = np.random.default_rng(seed)
        sample = rng.choice(values, SHAPIRO_MAX_SAMPLE, replace=False)
    try:
        return float(stats.shapiro(sample).pvalue)
    except Exception:
        return 0.0


def choose_transform(values: np.ndarray, allow_log: bool = True) -> tuple[str, np.ndarray, float]:
    """Pick identity vs log by whichever sample is closer to normal (higher
    Shapiro-Wilk p-value). log is only tried when every value is positive."""
    candidates = {"identity": values}
    if allow_log and np.all(values > 0):
        candidates["log"] = np.log(values)

    best_name, best_vals, best_p = "identity", values, -1.0
    for name, vals in candidates.items():
        p = shapiro_p(vals)
        log.info(
            "  transform=%-8s n=%d skew=%+.4f kurtosis=%+.4f shapiro_p=%.4g",
            name, len(vals), stats.skew(vals), stats.kurtosis(vals), p,
        )
        if p > best_p:
            best_name, best_vals, best_p = name, vals, p
    return best_name, best_vals, best_p


def fit_gaussian_metric(name: str, values: np.ndarray, allow_log: bool = True) -> dict:
    log.info("Fitting %s (n=%d):", name, len(values))
    transform, transformed, _ = choose_transform(values, allow_log=allow_log)
    mu = float(np.mean(transformed))
    sigma = float(np.std(transformed))
    log.info("  chosen transform=%s mu=%.6f sigma=%.6f", transform, mu, sigma)
    return {"transform": transform, "mu": mu, "sigma": sigma}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--train", default="cath17287")
    ap.add_argument("--eval", default="cath663")
    ap.add_argument("--out", type=Path, default=REPO / "benchmark/data/geometry_reference.json")
    ap.add_argument("--plot", type=Path, default=None, help="optional diagnostic histogram PNG")
    args = ap.parse_args()

    train = load_dataset(args.train)
    eval_entries = load_dataset(args.eval)
    eval_ids = {e.entry_id for e in eval_entries}
    # keep the fit clean of the held-out benchmark (same discipline as
    # fit_count_calibration.py)
    train = [e for e in train if e.entry_id not in eval_ids]

    log.info(
        "Fitting on %s (n=%d, %s excluded); cache has %d chains available",
        args.train, len(train), args.eval, sum(1 for _ in CHAINS.glob("*.pdb")),
    )
    samples = collect_samples(train)
    log.info(
        "Used %d/%d chains: %d domain samples, %d adjacent-pair samples",
        samples.n_chains_used, len(train), len(samples.kappa2), len(samples.interface_contacts),
    )
    if not samples.kappa2 or not samples.interface_contacts:
        raise SystemExit(
            "No domain/pair samples collected — is benchmark/cache/chains populated "
            "for this dataset?"
        )

    kappa2 = np.array(samples.kappa2)
    sphericity = fit_gaussian_metric("sphericity (kappa^2)", kappa2, allow_log=True)

    a, b = fit_power_law(np.array(samples.r_eff_n), np.array(samples.r_eff))
    log.info("Density power law: R_ideal(N) = %.6f * N^%.6f", a, b)
    ideal_rho = np.array([ca_number_density(n, a * n**b) for n in samples.density_n])
    log_residual = np.log(np.maximum(samples.density_rho, 1e-12)) - np.log(np.maximum(ideal_rho, 1e-12))
    # already log-domain and can be negative -> identity only
    log_density_residual = fit_gaussian_metric("log_density_residual", log_residual, allow_log=False)

    contacts = np.array(samples.interface_contacts, dtype=float)
    r_a = np.array(samples.interface_r_a)
    r_b = np.array(samples.interface_r_b)
    gamma = fit_gamma(contacts, r_a, r_b)
    log.info("Interface surface contact density gamma = %.6f", gamma)
    denom = np.pi * np.minimum(r_a, r_b) ** 2
    fraction = np.divide(contacts, denom, out=np.zeros_like(contacts), where=denom > 0)
    interface_fraction = fit_gaussian_metric("interface_fraction", fraction, allow_log=True)

    reference = {
        "_comment": (
            f"Fitted on {args.train} (n={samples.n_chains_used} chains, "
            f"{args.eval} held out) by benchmark/fit_geometry_reference.py."
        ),
        "sphericity": sphericity,
        "log_density_residual": log_density_residual,
        "interface_fraction": interface_fraction,
        "density_power_law": {"a": a, "b": b},
        "gamma": gamma,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(reference, indent=2) + "\n")
    log.info("Wrote %s", args.out)

    if args.plot:
        _plot_diagnostics(
            args.plot,
            {
                "sphericity (kappa^2)": kappa2,
                "log_density_residual": log_residual,
                "interface_fraction": fraction,
            },
            reference,
        )
        log.info("Wrote diagnostic histograms to %s", args.plot)


def _plot_diagnostics(path: Path, raw_values: dict[str, np.ndarray], reference: dict) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(raw_values), figsize=(5 * len(raw_values), 4))
    for ax, (name, values) in zip(np.atleast_1d(axes), raw_values.items()):
        ax.hist(values, bins=40, density=True, color="teal", alpha=0.7)
        ax.set_title(name)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)


if __name__ == "__main__":
    main()
