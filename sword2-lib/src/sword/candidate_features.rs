//! Per-candidate features for domain-partition reranking.
//!
//! These operate on the same 0-based delineation strings produced by
//! `compute_measure::MeasureLine` and consumed by `parse_measure` — i.e.
//! *before* `remap_residue_numbers` converts indices to original PDB
//! numbering.

use crate::peeling::algorithm::SsType;

/// Fraction of residues within `window` positions of a domain-boundary
/// junction that DSSP classifies as coil (loop/turn). Domain linkers are
/// structurally coil far more often than not, so a low fraction is a signal
/// of a boundary cutting through a helix or strand.
pub fn boundary_coil_fraction(raw_delineation: &str, ss_types: &[SsType]) -> f64 {
    const WINDOW: i64 = 2;
    let mut total = 0usize;
    let mut coil = 0usize;

    for domain in raw_delineation.trim().split_whitespace() {
        for segment in domain.split(';') {
            let parts: Vec<&str> = segment.split('-').collect();
            if parts.len() != 2 {
                continue;
            }
            let (Ok(start), Ok(end)) = (parts[0].parse::<i64>(), parts[1].parse::<i64>()) else {
                continue;
            };
            for junction in [start, end] {
                for offset in -WINDOW..=WINDOW {
                    let idx = junction + offset;
                    if idx >= 0 && (idx as usize) < ss_types.len() {
                        total += 1;
                        if ss_types[idx as usize] == SsType::Coil {
                            coil += 1;
                        }
                    }
                }
            }
        }
    }

    if total == 0 {
        0.0
    } else {
        coil as f64 / total as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boundary_coil_fraction_all_coil() {
        let ss_types = vec![SsType::Coil; 20];
        let frac = boundary_coil_fraction("0-9 10-19", &ss_types);
        assert_eq!(frac, 1.0);
    }

    #[test]
    fn test_boundary_coil_fraction_boundary_in_helix() {
        // Junction at index 9/10 is deep inside a helix run (indices 5..15).
        let mut ss_types = vec![SsType::Coil; 20];
        for s in ss_types.iter_mut().take(15).skip(5) {
            *s = SsType::Helix;
        }
        let frac = boundary_coil_fraction("0-9 10-19", &ss_types);
        assert!(frac < 0.5, "expected low coil fraction, got {frac}");
    }

    #[test]
    fn test_boundary_coil_fraction_discontinuous_segment() {
        let ss_types = vec![SsType::Coil; 20];
        // Domain 1 is discontinuous (two segments); still four junctions total.
        let frac = boundary_coil_fraction("0-4;15-19 5-14", &ss_types);
        assert_eq!(frac, 1.0);
    }

    #[test]
    fn test_boundary_coil_fraction_empty_delineation() {
        let ss_types = vec![SsType::Coil; 20];
        assert_eq!(boundary_coil_fraction("", &ss_types), 0.0);
    }
}
