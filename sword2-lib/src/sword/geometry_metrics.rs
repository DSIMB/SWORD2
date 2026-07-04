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

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

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
}
