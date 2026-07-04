//! Analytical "ideal sphere" geometric criteria for domain-partition selection.
//!
//! An ideal protein domain is modeled as a compact sphere. This module computes,
//! from raw Cα coordinates alone (no contact-probability matrix, no learned
//! weights), three purely geometric descriptors of how close a candidate domain
//! partition is to that ideal:
//!
//! 1. **Sphericity** — relative shape anisotropy κ² of the gyration tensor
//!    (0 = perfect sphere, 1 = a line).
//! 2. **Compactness** — Cα number density compared to the density expected for
//!    a domain of that size (fitted offline as `R_ideal(N) = a * N^b`).
//! 3. **Inter-domain contacts** — the fraction of the maximum possible contact
//!    interface (as if the two domains were touching half-spheres) that a pair
//!    of adjacent candidate domains actually realizes.
//!
//! Each descriptor is compared against a Gaussian fitted on true CATH domains
//! (see `benchmark/fit_geometry_reference.py`), producing a z-score and
//! p-value. This module deliberately does not implement any predictive model:
//! it is a self-contained companion to (not a replacement for) the
//! distance_model / candidate_features / reranker selection machinery, used
//! only to report and optionally bias candidate selection toward
//! geometrically plausible domains.

use std::f64::consts::PI;
use std::sync::LazyLock;

/// Hard Cα-Cα contact cutoff (Å) used for inter-domain contact counting.
/// Distinct from — and independent of — the sigmoid contact-probability
/// matrix used by peeling; this is a simple geometric distance threshold.
pub const CONTACT_CUTOFF_ANGSTROM: f64 = 8.0;

/// Centroid of the Cα coordinates at `indices`.
pub fn centroid(coords: &[[f64; 3]], indices: &[usize]) -> [f64; 3] {
    let n = indices.len() as f64;
    let mut c = [0.0; 3];
    for &i in indices {
        c[0] += coords[i][0];
        c[1] += coords[i][1];
        c[2] += coords[i][2];
    }
    if n > 0.0 {
        for v in &mut c {
            *v /= n;
        }
    }
    c
}

/// Gyration tensor `S = (1/N) * sum (x_i - c)(x_i - c)^T` for the residues at
/// `indices`, a symmetric 3x3 matrix returned as `[row][col]`.
pub fn gyration_tensor(coords: &[[f64; 3]], indices: &[usize]) -> [[f64; 3]; 3] {
    let c = centroid(coords, indices);
    let n = indices.len() as f64;
    let mut s = [[0.0; 3]; 3];
    for &i in indices {
        let d = [
            coords[i][0] - c[0],
            coords[i][1] - c[1],
            coords[i][2] - c[2],
        ];
        for a in 0..3 {
            for b in 0..3 {
                s[a][b] += d[a] * d[b];
            }
        }
    }
    if n > 0.0 {
        for row in &mut s {
            for v in row {
                *v /= n;
            }
        }
    }
    s
}

/// Analytical eigenvalues of a symmetric 3x3 matrix (Cardano's trigonometric
/// solution), sorted descending. Standard closed form for the real symmetric
/// case — see e.g. Smith, "Eigenvalues of a symmetric 3x3 matrix" (1961).
pub fn sym3x3_eigenvalues(m: [[f64; 3]; 3]) -> [f64; 3] {
    let p1 = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];
    if p1 <= 1e-12 {
        let mut eig = [m[0][0], m[1][1], m[2][2]];
        eig.sort_by(|a, b| b.partial_cmp(a).unwrap());
        return eig;
    }

    let q = (m[0][0] + m[1][1] + m[2][2]) / 3.0;
    let p2 = (m[0][0] - q).powi(2) + (m[1][1] - q).powi(2) + (m[2][2] - q).powi(2) + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();

    let mut b = [[0.0; 3]; 3];
    for a in 0..3 {
        for c in 0..3 {
            let diag = if a == c { q } else { 0.0 };
            b[a][c] = (m[a][c] - diag) / p;
        }
    }
    let det_b = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
        - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
        + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);

    let r = (det_b / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;

    let eig1 = q + 2.0 * p * phi.cos();
    let eig3 = q + 2.0 * p * (phi + 2.0 * PI / 3.0).cos();
    let eig2 = 3.0 * q - eig1 - eig3;

    let mut eig = [eig1, eig2, eig3];
    eig.sort_by(|a, b| b.partial_cmp(a).unwrap());
    eig
}

/// Principal radii `r_k = sqrt(lambda_k)`, from descending-sorted eigenvalues.
pub fn principal_radii(eigenvalues: [f64; 3]) -> [f64; 3] {
    [
        eigenvalues[0].max(0.0).sqrt(),
        eigenvalues[1].max(0.0).sqrt(),
        eigenvalues[2].max(0.0).sqrt(),
    ]
}

/// Radius of gyration `Rg = sqrt(lambda1 + lambda2 + lambda3)`.
pub fn radius_of_gyration(eigenvalues: [f64; 3]) -> f64 {
    (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]).max(0.0).sqrt()
}

/// Relative shape anisotropy kappa^2 in `[0, 1]`: 0 for a perfect sphere, 1
/// for a line. `kappa^2 = 1 - 3*(l1*l2 + l2*l3 + l3*l1) / (l1+l2+l3)^2`.
pub fn relative_shape_anisotropy(eigenvalues: [f64; 3]) -> f64 {
    let sum = eigenvalues[0] + eigenvalues[1] + eigenvalues[2];
    if sum <= 0.0 {
        return 0.0;
    }
    let pair_sum = eigenvalues[0] * eigenvalues[1]
        + eigenvalues[1] * eigenvalues[2]
        + eigenvalues[2] * eigenvalues[0];
    (1.0 - 3.0 * pair_sum / (sum * sum)).clamp(0.0, 1.0)
}

/// Radius of the sphere with the same radius of gyration as the domain
/// (uniform-density sphere: `Rg^2 = (3/5) R^2`).
pub fn effective_sphere_radius(rg: f64) -> f64 {
    rg * (5.0f64 / 3.0).sqrt()
}

/// Cα number density `N / V` of the effective sphere of radius `r_eff`.
pub fn ca_number_density(n: usize, r_eff: f64) -> f64 {
    if r_eff <= 0.0 {
        return 0.0;
    }
    let v = (4.0 / 3.0) * PI * r_eff.powi(3);
    n as f64 / v
}

/// Number of Cα-Cα pairs between `idx_a` and `idx_b` within `cutoff` Å.
pub fn interdomain_contacts(
    coords: &[[f64; 3]],
    idx_a: &[usize],
    idx_b: &[usize],
    cutoff: f64,
) -> usize {
    let cutoff_sq = cutoff * cutoff;
    let mut count = 0usize;
    for &i in idx_a {
        for &j in idx_b {
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            if dx * dx + dy * dy + dz * dz <= cutoff_sq {
                count += 1;
            }
        }
    }
    count
}

/// Maximum contact count if the two domains met across a shared great-circle
/// disk of radius `min(R_a, R_b)` (the "two half-spheres" bound), given a
/// fitted surface contact density `gamma` (residues per Å^2).
pub fn max_half_sphere_contacts(r_a: f64, r_b: f64, gamma: f64) -> f64 {
    let r_min = r_a.min(r_b);
    gamma * PI * r_min * r_min
}

/// Observed contacts as a fraction of the half-sphere maximum: ~0 means the
/// domains barely touch (tangent spheres — a clean cut), ~1 means the
/// interface is as wide as a full bisection (a likely over-split).
pub fn interface_fraction(contacts: usize, r_a: f64, r_b: f64, gamma: f64) -> f64 {
    let c_half = max_half_sphere_contacts(r_a, r_b, gamma);
    if c_half <= 0.0 {
        return 0.0;
    }
    contacts as f64 / c_half
}

/// Parse a raw (0-based) delineation string into per-domain residue index
/// lists. Domains are whitespace-separated; discontinuous segments within a
/// domain are `;`-separated; each segment is `start-end` (inclusive). This is
/// a local, self-contained parser — deliberately not shared with
/// `candidate_features.rs`.
pub fn parse_domain_indices(raw_delineation: &str) -> Vec<Vec<usize>> {
    raw_delineation
        .split_whitespace()
        .map(|domain_tok| {
            domain_tok
                .split(';')
                .filter(|seg| !seg.is_empty())
                .flat_map(|seg| {
                    let mut parts = seg.splitn(2, '-');
                    let start = parts.next().and_then(|s| s.parse::<usize>().ok());
                    let end = parts.next().and_then(|s| s.parse::<usize>().ok());
                    match (start, end) {
                        (Some(s), Some(e)) if s <= e => (s..=e).collect::<Vec<_>>(),
                        (Some(s), None) => vec![s],
                        _ => Vec::new(),
                    }
                })
                .collect()
        })
        .collect()
}

/// A transform applied to a raw metric value before comparing it against a
/// fitted Gaussian, chosen offline (per metric) by whichever makes the
/// reference sample closer to normal — see `benchmark/fit_geometry_reference.py`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Transform {
    Identity,
    /// `ln(max(x, epsilon))` — only sensible for non-negative raw metrics.
    Log,
}

impl Transform {
    fn apply(self, x: f64) -> f64 {
        match self {
            Transform::Identity => x,
            Transform::Log => x.max(1e-12).ln(),
        }
    }
}

/// Fitted Gaussian `N(mu, sigma)` for one (possibly transformed) geometric
/// residual.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct GaussianParams {
    pub transform: Transform,
    pub mu: f64,
    pub sigma: f64,
}

impl GaussianParams {
    /// z-score of `raw_value` after applying this metric's fitted transform.
    pub fn z(&self, raw_value: f64) -> f64 {
        (self.transform.apply(raw_value) - self.mu) / self.sigma
    }
}

/// Fitted `R_ideal(N) = a * N^b` domain-size-to-radius scaling law.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct DensityPowerLaw {
    pub a: f64,
    pub b: f64,
}

/// Reference Gaussians and density scaling law fitted offline on true CATH
/// domains by `benchmark/fit_geometry_reference.py`, embedded at compile time
/// from `benchmark/data/geometry_reference.json`.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct ReferenceDistributions {
    pub sphericity: GaussianParams,
    pub log_density_residual: GaussianParams,
    pub interface_fraction: GaussianParams,
    pub density_power_law: DensityPowerLaw,
    pub gamma: f64,
}

const REFERENCE_JSON: &str = include_str!("../../../benchmark/data/geometry_reference.json");

static DEFAULT_REFERENCE: LazyLock<ReferenceDistributions> = LazyLock::new(|| {
    serde_json::from_str(REFERENCE_JSON)
        .expect("embedded benchmark/data/geometry_reference.json must be valid")
});

impl ReferenceDistributions {
    /// The compiled-in reference fitted on CATH (see module docs).
    pub fn embedded() -> &'static ReferenceDistributions {
        &DEFAULT_REFERENCE
    }

    /// z-score of a domain's kappa^2 against the reference sphericity Gaussian.
    pub fn sphericity_z(&self, kappa2: f64) -> f64 {
        self.sphericity.z(kappa2)
    }

    /// Expected Ca density for a domain of `n` residues under the fitted
    /// ideal-sphere scaling law `R_ideal(N) = a * N^b`.
    pub fn ideal_density(&self, n: usize) -> f64 {
        let r_ideal = self.density_power_law.a * (n as f64).powf(self.density_power_law.b);
        ca_number_density(n, r_ideal)
    }

    /// z-score of a domain's log-density residual (`ln density - ln ideal`)
    /// against the reference. The residual is already log-domain by
    /// construction and can be negative, so its own `transform` is expected
    /// to always be `identity`.
    pub fn density_z(&self, n: usize, density: f64) -> f64 {
        let ideal = self.ideal_density(n).max(1e-12);
        let delta = density.max(1e-12).ln() - ideal.ln();
        self.log_density_residual.z(delta)
    }

    /// z-score of an inter-domain interface fraction against the reference.
    pub fn interface_z(&self, fraction: f64) -> f64 {
        self.interface_fraction.z(fraction)
    }
}

/// Standard normal CDF via the Abramowitz & Stegun 7.1.26 erf approximation
/// (max absolute error 1.5e-7) — sufficient precision for reporting p-values.
fn std_normal_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf(z / std::f64::consts::SQRT_2))
}

fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;
    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
    sign * y
}

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

/// Per-partition aggregate of the three geometric criteria: worst-case (most
/// anomalous) and mean z-score for sphericity and density across domains, and
/// for the interface fraction across sequentially adjacent domain pairs
/// (`None` for single-domain partitions). `geometry_score` is the fixed-weight
/// analytical penalty `G = kappa_z+ + density_deficit_z+ + interface_z+`
/// (positive parts only — only "worse than ideal" deviations are penalized;
/// weights are all 1, documented here rather than learned).
#[derive(Debug, Clone, Copy)]
pub struct GeometryReport {
    pub worst_sphericity_z: f64,
    pub mean_sphericity_z: f64,
    pub worst_density_z: f64,
    pub mean_density_z: f64,
    pub worst_interface_z: Option<f64>,
    pub mean_interface_z: Option<f64>,
    pub geometry_score: f64,
}

impl GeometryReport {
    /// Upper-tail p-value: probability of a domain this elongated or more.
    pub fn sphericity_p(&self) -> f64 {
        1.0 - std_normal_cdf(self.worst_sphericity_z)
    }

    /// Lower-tail p-value: probability of a domain this loose (or looser).
    pub fn density_p(&self) -> f64 {
        std_normal_cdf(self.worst_density_z)
    }

    /// Upper-tail p-value: probability of an interface this wide or wider.
    pub fn interface_p(&self) -> Option<f64> {
        self.worst_interface_z.map(|z| 1.0 - std_normal_cdf(z))
    }
}

/// Compute the full geometry report for one candidate partition. `domains` is
/// the per-domain list of (0-based) Cα indices, in the same order as they
/// appear in the raw delineation (see `parse_domain_indices`); consecutive
/// entries are treated as sequentially adjacent for the interface criterion.
pub fn geometry_report(
    coords: &[[f64; 3]],
    domains: &[Vec<usize>],
    reference: &ReferenceDistributions,
) -> GeometryReport {
    let mut sphericity_zs = Vec::with_capacity(domains.len());
    let mut density_zs = Vec::with_capacity(domains.len());
    let mut r_effs: Vec<Option<f64>> = Vec::with_capacity(domains.len());

    for dom in domains {
        if dom.len() < 2 {
            r_effs.push(None);
            continue;
        }
        let eig = sym3x3_eigenvalues(gyration_tensor(coords, dom));
        let kappa2 = relative_shape_anisotropy(eig);
        let r_eff = effective_sphere_radius(radius_of_gyration(eig));
        let density = ca_number_density(dom.len(), r_eff);
        sphericity_zs.push(reference.sphericity_z(kappa2));
        density_zs.push(reference.density_z(dom.len(), density));
        r_effs.push(Some(r_eff));
    }

    let mut interface_zs = Vec::new();
    for i in 0..domains.len().saturating_sub(1) {
        let (Some(r_a), Some(r_b)) = (r_effs[i], r_effs[i + 1]) else {
            continue;
        };
        let contacts =
            interdomain_contacts(coords, &domains[i], &domains[i + 1], CONTACT_CUTOFF_ANGSTROM);
        let fraction = interface_fraction(contacts, r_a, r_b, reference.gamma);
        interface_zs.push(reference.interface_z(fraction));
    }

    let (worst_sphericity_z, mean_sphericity_z) = if sphericity_zs.is_empty() {
        (0.0, 0.0)
    } else {
        (
            sphericity_zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            mean(&sphericity_zs),
        )
    };
    let (worst_density_z, mean_density_z) = if density_zs.is_empty() {
        (0.0, 0.0)
    } else {
        (
            density_zs.iter().cloned().fold(f64::INFINITY, f64::min),
            mean(&density_zs),
        )
    };
    let (worst_interface_z, mean_interface_z) = if interface_zs.is_empty() {
        (None, None)
    } else {
        (
            Some(interface_zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max)),
            Some(mean(&interface_zs)),
        )
    };

    let geometry_score = worst_sphericity_z.max(0.0)
        + (-worst_density_z).max(0.0)
        + worst_interface_z.map(|z| z.max(0.0)).unwrap_or(0.0);

    GeometryReport {
        worst_sphericity_z,
        mean_sphericity_z,
        worst_density_z,
        mean_density_z,
        worst_interface_z,
        mean_interface_z,
        geometry_score,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn test_reference() -> ReferenceDistributions {
        ReferenceDistributions {
            sphericity: GaussianParams { transform: Transform::Identity, mu: 0.3, sigma: 0.1 },
            log_density_residual: GaussianParams {
                transform: Transform::Identity,
                mu: 0.0,
                sigma: 0.2,
            },
            interface_fraction: GaussianParams {
                transform: Transform::Identity,
                mu: 0.2,
                sigma: 0.1,
            },
            density_power_law: DensityPowerLaw { a: 3.0, b: 0.4 },
            gamma: 0.05,
        }
    }

    #[test]
    fn log_transform_applies_ln_with_epsilon_floor() {
        let g = GaussianParams { transform: Transform::Log, mu: 0.0, sigma: 1.0 };
        // ln(1.0) == 0.0 == mu -> z == 0
        assert!(g.z(1.0).abs() < EPS);
        // ln(e) == 1.0 -> z == 1
        assert!((g.z(std::f64::consts::E) - 1.0).abs() < EPS);
        // non-positive input is floored, not -inf/NaN
        assert!(g.z(0.0).is_finite());
        assert!(g.z(-5.0).is_finite());
    }

    #[test]
    fn octahedron_is_a_perfect_sphere() {
        let coords = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let idx: Vec<usize> = (0..6).collect();
        let s = gyration_tensor(&coords, &idx);
        let eig = sym3x3_eigenvalues(s);
        assert!((eig[0] - eig[1]).abs() < EPS);
        assert!((eig[1] - eig[2]).abs() < EPS);
        let kappa2 = relative_shape_anisotropy(eig);
        assert!(kappa2 < EPS, "expected kappa^2 ~ 0, got {kappa2}");
    }

    #[test]
    fn collinear_points_have_unit_anisotropy() {
        let coords: Vec<[f64; 3]> = (-2..=2).map(|x| [x as f64, 0.0, 0.0]).collect();
        let idx: Vec<usize> = (0..coords.len()).collect();
        let s = gyration_tensor(&coords, &idx);
        let eig = sym3x3_eigenvalues(s);
        let kappa2 = relative_shape_anisotropy(eig);
        assert!((kappa2 - 1.0).abs() < EPS, "expected kappa^2 ~ 1, got {kappa2}");
    }

    #[test]
    fn eigenvalues_match_known_symmetric_matrix() {
        // [[2,1,0],[1,2,0],[0,0,3]] has eigenvalues {3, 3, 1}.
        let m = [[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let eig = sym3x3_eigenvalues(m);
        assert!((eig[0] - 3.0).abs() < 1e-6);
        assert!((eig[1] - 3.0).abs() < 1e-6);
        assert!((eig[2] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn diagonal_matrix_shortcut_sorts_descending() {
        let m = [[1.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 2.0]];
        let eig = sym3x3_eigenvalues(m);
        assert_eq!(eig, [3.0, 2.0, 1.0]);
    }

    #[test]
    fn principal_radii_and_rg_are_consistent() {
        let eig = [4.0, 1.0, 0.0];
        let radii = principal_radii(eig);
        assert!((radii[0] - 2.0).abs() < EPS);
        assert!((radii[1] - 1.0).abs() < EPS);
        assert!((radii[2] - 0.0).abs() < EPS);
        let rg = radius_of_gyration(eig);
        assert!((rg - 5.0f64.sqrt()).abs() < EPS);
    }

    #[test]
    fn density_matches_sphere_volume_formula() {
        let r_eff = effective_sphere_radius(1.0);
        assert!((r_eff - (5.0 / 3.0f64).sqrt()).abs() < EPS);
        let n = 100;
        let density = ca_number_density(n, r_eff);
        let expected = n as f64 / ((4.0 / 3.0) * PI * r_eff.powi(3));
        assert!((density - expected).abs() < EPS);
    }

    #[test]
    fn zero_radius_gives_zero_density() {
        assert_eq!(ca_number_density(10, 0.0), 0.0);
    }

    #[test]
    fn interdomain_contacts_counts_pairs_within_cutoff() {
        let coords = vec![
            [0.0, 0.0, 0.0],  // 0: domain a
            [10.0, 0.0, 0.0], // 1: domain a
            [3.0, 0.0, 0.0],  // 2: domain b, 3A from 0 -> contact
            [3.0, 0.0, 0.0],  // 3: domain b, 3A from 0 -> contact, 7A from 1 -> contact
            [50.0, 0.0, 0.0], // 4: domain b, far from both -> no contact
        ];
        let idx_a = [0usize, 1];
        let idx_b = [2usize, 3, 4];
        let count = interdomain_contacts(&coords, &idx_a, &idx_b, CONTACT_CUTOFF_ANGSTROM);
        // 0-2 (3A, yes), 0-3 (3A, yes), 1-2 (7A, yes), 1-3 (7A, yes), rest no.
        assert_eq!(count, 4);
    }

    #[test]
    fn interface_fraction_zero_when_no_contacts() {
        assert_eq!(interface_fraction(0, 5.0, 5.0, 0.1), 0.0);
    }

    #[test]
    fn interface_fraction_scales_with_contacts() {
        let r_a = 5.0;
        let r_b = 8.0;
        let gamma = 0.1;
        let c_half = max_half_sphere_contacts(r_a, r_b, gamma);
        let half = interface_fraction(c_half.round() as usize, r_a, r_b, gamma);
        // c_half itself is not an integer, so rounding to the nearest contact
        // count introduces a small (~1/c_half) relative error.
        assert!((half - 1.0).abs() < 0.05, "expected ~1.0, got {half}");
        let none = interface_fraction(0, r_a, r_b, gamma);
        assert_eq!(none, 0.0);
    }

    #[test]
    fn parse_domain_indices_simple() {
        let domains = parse_domain_indices("0-2 3-5");
        assert_eq!(domains, vec![vec![0, 1, 2], vec![3, 4, 5]]);
    }

    #[test]
    fn parse_domain_indices_discontinuous_segments() {
        let domains = parse_domain_indices("0-2;7-8 3-6");
        assert_eq!(domains, vec![vec![0, 1, 2, 7, 8], vec![3, 4, 5, 6]]);
    }

    #[test]
    fn parse_domain_indices_empty_string() {
        let domains = parse_domain_indices("");
        assert!(domains.is_empty());
    }

    #[test]
    fn sphericity_z_zero_at_reference_mean() {
        let r = test_reference();
        assert!(r.sphericity_z(0.3).abs() < EPS);
    }

    #[test]
    fn interface_z_zero_at_reference_mean() {
        let r = test_reference();
        assert!(r.interface_z(0.2).abs() < EPS);
    }

    #[test]
    fn density_z_zero_when_matching_ideal_law() {
        let r = test_reference();
        let ideal = r.ideal_density(120);
        assert!(r.density_z(120, ideal).abs() < 1e-9);
    }

    #[test]
    fn std_normal_cdf_matches_known_values() {
        assert!((std_normal_cdf(0.0) - 0.5).abs() < 1e-6);
        assert!((std_normal_cdf(1.0) - 0.8413).abs() < 1e-3);
        assert!((std_normal_cdf(-1.0) - 0.1587).abs() < 1e-3);
    }

    #[test]
    fn geometry_report_two_domain_partition_computes_all_metrics() {
        // Two small tetrahedra, far apart -> low/near-zero interface fraction.
        let offsets = [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, -1.0],
            [-1.0, 1.0, -1.0],
            [-1.0, -1.0, 1.0],
        ];
        let mut coords: Vec<[f64; 3]> = offsets.to_vec();
        coords.extend(offsets.iter().map(|o| [o[0] + 100.0, o[1], o[2]]));
        let domains = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]];
        let r = test_reference();
        let report = geometry_report(&coords, &domains, &r);

        assert!(report.worst_interface_z.is_some());
        assert!(report.mean_interface_z.is_some());
        assert!(report.geometry_score >= 0.0 && report.geometry_score.is_finite());
        assert!((0.0..=1.0).contains(&report.sphericity_p()));
        assert!((0.0..=1.0).contains(&report.density_p()));
        assert!((0.0..=1.0).contains(&report.interface_p().unwrap()));
    }

    #[test]
    fn geometry_report_single_domain_has_no_interface() {
        let coords = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let domains = vec![vec![0, 1, 2, 3]];
        let r = test_reference();
        let report = geometry_report(&coords, &domains, &r);
        assert!(report.worst_interface_z.is_none());
        assert!(report.mean_interface_z.is_none());
    }

    #[test]
    fn embedded_reference_json_parses_and_computes() {
        let reference = ReferenceDistributions::embedded();
        let coords = vec![
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 3.0],
            [10.0, 10.0, 10.0],
        ];
        let domains = vec![vec![0, 1, 2, 3]];
        let report = geometry_report(&coords, &domains, reference);
        assert!(report.geometry_score.is_finite());
    }
}
