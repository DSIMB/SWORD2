//! Backbone angle calculations: phi, psi, kappa, alpha, tco.

use crate::pdb::types::Point3D;

use super::types::{DsspChain, RADIAN};

/// Calculate all backbone angles (kappa, alpha, chirality) for the chain.
pub fn calculate_angles(chain: &mut DsspChain) {
    let len = chain.len;

    // Alpha (virtual dihedral CA[i-1]..CA[i+2]) and chirality
    for i in 2..=len.saturating_sub(2) {
        if chain.no_chain_break(i - 1, i + 2) {
            let alpha = dihedral_angle(
                &chain.get(i - 1).ca,
                &chain.get(i).ca,
                &chain.get(i + 1).ca,
                &chain.get(i + 2).ca,
            );
            chain.get_mut(i).alpha = alpha;
            chain.get_mut(i).ss[5] = if alpha < 0.0 { '-' } else { '+' };
        }
    }

    // Kappa (virtual bend angle CA[i-2]..CA[i]..CA[i+2])
    for i in 3..=len.saturating_sub(2) {
        if chain.no_chain_break(i - 2, i + 2) {
            let ckap = cos_angle(
                &chain.get(i).ca,
                &chain.get(i - 2).ca,
                &chain.get(i + 2).ca,
                &chain.get(i).ca,
            );
            let skap = (1.0 - ckap * ckap).max(0.0).sqrt();
            chain.get_mut(i).kappa = RADIAN * atan2(skap, ckap);
        }
    }
}

/// Dihedral angle (torsion) of four points v1-v2-v3-v4, in degrees.
pub fn dihedral_angle(v1: &Point3D, v2: &Point3D, v3: &Point3D, v4: &Point3D) -> f64 {
    let v12 = *v1 - *v2;
    let v43 = *v4 - *v3;
    let z = *v2 - *v3;

    let p = cross(&z, &v12);
    let x = cross(&z, &v43);
    let y = cross(&z, &x);

    let u_sq = dot(&x, &x);
    let v_sq = dot(&y, &y);

    if u_sq <= 0.0 || v_sq <= 0.0 {
        return 360.0;
    }

    let u = dot(&p, &x) / u_sq.sqrt();
    let v = dot(&p, &y) / v_sq.sqrt();

    if u != 0.0 || v != 0.0 {
        atan2(v, u) * RADIAN
    } else {
        360.0
    }
}

/// Cosine of angle between vectors (v1-v2) and (v3-v4).
pub fn cos_angle(v1: &Point3D, v2: &Point3D, v3: &Point3D, v4: &Point3D) -> f64 {
    let u = *v1 - *v2;
    let v = *v3 - *v4;
    let x = dot(&u, &u) * dot(&v, &v);
    if x > 0.0 {
        dot(&u, &v) / x.sqrt()
    } else {
        0.0
    }
}

/// Cosine of C=O[i] vs C=O[i-1] (TCO value).
pub fn tco(chain: &DsspChain, i: usize) -> f64 {
    if chain.no_chain_break(i - 1, i) {
        cos_angle(
            &chain.get(i).c,
            &chain.get(i).o,
            &chain.get(i - 1).c,
            &chain.get(i - 1).o,
        )
    } else {
        0.0
    }
}

/// Phi angle for residue i.
pub fn phi(chain: &DsspChain, i: usize) -> f64 {
    if chain.no_chain_break(i - 1, i) {
        dihedral_angle(
            &chain.get(i - 1).c,
            &chain.get(i).n,
            &chain.get(i).ca,
            &chain.get(i).c,
        )
    } else {
        360.0
    }
}

/// Psi angle for residue i.
pub fn psi(chain: &DsspChain, i: usize) -> f64 {
    if i < chain.len && chain.no_chain_break(i, i + 1) {
        dihedral_angle(
            &chain.get(i).n,
            &chain.get(i).ca,
            &chain.get(i).c,
            &chain.get(i + 1).n,
        )
    } else {
        360.0
    }
}

// Vector helpers operating on Point3D

fn cross(a: &Point3D, b: &Point3D) -> Point3D {
    Point3D::new(
        a.y * b.z - b.y * a.z,
        a.z * b.x - b.z * a.x,
        a.x * b.y - b.x * a.y,
    )
}

fn dot(a: &Point3D, b: &Point3D) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

fn atan2(y: f64, x: f64) -> f64 {
    // Match the original DSSP Atan2 behavior exactly
    let z;
    if x != 0.0 {
        z = (y / x).atan();
    } else if y > 0.0 {
        z = std::f64::consts::FRAC_PI_2;
    } else if y < 0.0 {
        z = -std::f64::consts::FRAC_PI_2;
    } else {
        return std::f64::consts::TAU; // 2*PI
    }

    if x >= 0.0 {
        z
    } else if y > 0.0 {
        z + std::f64::consts::PI
    } else {
        z - std::f64::consts::PI
    }
}
