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

/// Most frequent `num_domains` value across one chain's own candidate set.
/// Ties broken toward the smaller count for determinism.
pub fn modal_num_domains(counts: &[usize]) -> usize {
    let mut freq: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new();
    for &nd in counts {
        *freq.entry(nd).or_insert(0) += 1;
    }
    freq.into_iter()
        .max_by_key(|&(nd, count)| (count, std::cmp::Reverse(nd)))
        .map(|(nd, _)| nd)
        .unwrap_or(0)
}

/// Absolute distance of a candidate's domain count from its chain's modal count.
pub fn modal_count_distance(num_domains: usize, modal: usize) -> f64 {
    (num_domains as f64 - modal as f64).abs()
}

/// Pseudo-energy Z-score for one candidate, computed per-domain and averaged.
///
/// Every candidate for a given protein partitions the *same* full residue
/// set, just grouped into different domains — scoring all of a candidate's
/// residues as a single concatenated blob would therefore score identically
/// for every candidate regardless of how the boundaries are drawn. Each
/// domain must be scored on its own (does *this* proposed domain look like a
/// real fold?), matching how `energy::calculate_all_energies` scores the
/// final selected partition. `remapped_delineation` must use original PDB
/// residue numbers (i.e. already passed through `remap_residue_numbers`),
/// unlike `boundary_coil_fraction`'s `raw_delineation`.
pub fn candidate_energy_z_score(
    energy_config: &crate::energy::EnergyConfig,
    pdb_path: &str,
    chain: &str,
    remapped_delineation: &str,
) -> Option<f64> {
    let mut z_scores: Vec<f64> = Vec::new();
    for domain in remapped_delineation.trim().split_whitespace() {
        let mut residues = String::new();
        for segment in domain.split(';') {
            let parts: Vec<&str> = segment.split('-').collect();
            if parts.len() != 2 {
                continue;
            }
            let (Ok(start), Ok(end)) = (parts[0].parse::<i32>(), parts[1].parse::<i32>()) else {
                continue;
            };
            let list = crate::energy::build_residue_list((start, end), chain);
            if !residues.is_empty() {
                residues.push(',');
            }
            residues.push_str(&list);
        }
        if residues.is_empty() {
            continue;
        }
        if let Ok(result) =
            crate::energy::get_energy_and_z_score(energy_config, pdb_path, Some(&residues))
        {
            if let Some(z) = result.z_score {
                z_scores.push(z);
            }
        }
    }
    if z_scores.is_empty() {
        None
    } else {
        Some(z_scores.iter().sum::<f64>() / z_scores.len() as f64)
    }
}

/// One row of dump/training-table features for a single candidate, computed
/// at the same shortlist stage the reranker scores at inference time.
pub struct DumpRow {
    pub num_domains: usize,
    pub min_size: usize,
    pub max_cr: f64,
    pub density_min: f64,
    pub mean_density: f64,
    pub delineation: String,
    pub boundary_coil_fraction: f64,
    pub energy_z: Option<f64>,
    pub modal_count_distance: f64,
}

/// Build a `DumpRow` from one `relevant_measure2` pipe-delimited line.
///
/// `raw_delineation` uses 0-based indices (matches `ss_types`); the caller
/// supplies `remapped_delineation` (original PDB numbering) separately since
/// only `candidate_energy_z_score` needs it. `modal` is
/// `modal_num_domains(...)` computed once across the whole shortlist by the
/// caller (not per-row) and passed in.
#[allow(clippy::too_many_arguments)]
pub fn build_dump_row(
    num_domains: usize,
    min_size: usize,
    max_cr: f64,
    density_min: f64,
    mean_density: f64,
    raw_delineation: &str,
    remapped_delineation: &str,
    ss_types: &[SsType],
    energy_config: Option<&crate::energy::EnergyConfig>,
    pdb_path: &str,
    chain: &str,
    modal: usize,
) -> DumpRow {
    let coil_fraction = boundary_coil_fraction(raw_delineation, ss_types);
    let energy_z = energy_config
        .and_then(|ec| candidate_energy_z_score(ec, pdb_path, chain, remapped_delineation));
    DumpRow {
        num_domains,
        min_size,
        max_cr,
        density_min,
        mean_density,
        delineation: raw_delineation.to_string(),
        boundary_coil_fraction: coil_fraction,
        energy_z,
        modal_count_distance: modal_count_distance(num_domains, modal),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_dump_row_without_energy_config() {
        let ss_types = vec![SsType::Coil; 20];
        let row = build_dump_row(
            2, 10, 0.5, 1.0, 2.0,
            "0-9 10-19", "1-10 11-20",
            &ss_types, None, "/nonexistent.pdb", "A", 2,
        );
        assert_eq!(row.num_domains, 2);
        assert_eq!(row.boundary_coil_fraction, 1.0);
        assert_eq!(row.energy_z, None);
        assert_eq!(row.modal_count_distance, 0.0);
    }

    #[test]
    fn test_candidate_energy_z_score_returns_finite_value() {
        use std::io::Write as _;
        let dir = tempfile::tempdir().unwrap();
        let pdb_path = dir.path().join("frag.pdb");
        let mut f = std::fs::File::create(&pdb_path).unwrap();
        // Residues 1-2 of 1JX4_A backbone atoms — see sword2-lib/tests/dssp_golden.rs
        // for the full fixture and provenance.
        writeln!(f, "ATOM      1  N   ILE A   1      46.170  17.543  13.913  1.00 33.83           N").unwrap();
        writeln!(f, "ATOM      2  CA  ILE A   1      45.665  16.827  12.751  1.00 32.01           C").unwrap();
        writeln!(f, "ATOM      3  C   ILE A   1      44.396  16.160  13.297  1.00 29.09           C").unwrap();
        writeln!(f, "ATOM      4  O   ILE A   1      44.466  15.419  14.274  1.00 27.71           O").unwrap();
        writeln!(f, "ATOM      9  N   VAL A   2      43.247  16.494  12.713  1.00 28.04           N").unwrap();
        writeln!(f, "ATOM     10  CA  VAL A   2      41.973  15.941  13.134  1.00 26.25           C").unwrap();
        writeln!(f, "ATOM     11  C   VAL A   2      41.475  15.044  12.028  1.00 27.46           C").unwrap();
        writeln!(f, "ATOM     12  O   VAL A   2      41.543  15.426  10.857  1.00 26.95           O").unwrap();
        drop(f);

        let bin_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../bin");
        let mut ec = crate::energy::EnergyConfig::from_bin_dir(bin_dir);
        ec.num_shuffles = 20; // keep the test fast; production reranking uses 200
        ec.preload().unwrap();

        let z = candidate_energy_z_score(&ec, pdb_path.to_str().unwrap(), "A", "1-2");
        assert!(z.is_some(), "expected a Z-score for a 2-residue fragment");
        assert!(z.unwrap().is_finite());
    }

    #[test]
    fn test_candidate_energy_z_score_empty_delineation_is_none() {
        let bin_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../bin");
        let ec = crate::energy::EnergyConfig::from_bin_dir(bin_dir);
        assert_eq!(candidate_energy_z_score(&ec, "/nonexistent.pdb", "A", ""), None);
    }

    /// Regression test: different domain splits of the *same* residue range
    /// must score differently. A prior version concatenated every domain's
    /// residues into one combined list before scoring, which made every
    /// candidate covering the full protein score identically regardless of
    /// where the boundaries were drawn (caught via manual smoke-testing the
    /// dump path on 1JX4_A — every one of 23 candidates got the same
    /// energy_z). Scoring must happen per-domain and be aggregated.
    #[test]
    fn test_candidate_energy_z_score_differs_by_domain_split() {
        // Residues 1-20 of 1JX4_A backbone atoms, matching
        // sword2-lib/tests/dssp_golden.rs::FIXTURE_PDB.
        const FIXTURE_PDB_20: &str = "\
ATOM      1  N   ILE A   1      46.170  17.543  13.913  1.00 33.83           N
ATOM      2  CA  ILE A   1      45.665  16.827  12.751  1.00 32.01           C
ATOM      3  C   ILE A   1      44.396  16.160  13.297  1.00 29.09           C
ATOM      4  O   ILE A   1      44.466  15.419  14.274  1.00 27.71           O
ATOM      9  N   VAL A   2      43.247  16.494  12.713  1.00 28.04           N
ATOM     10  CA  VAL A   2      41.973  15.941  13.134  1.00 26.25           C
ATOM     11  C   VAL A   2      41.475  15.044  12.028  1.00 27.46           C
ATOM     12  O   VAL A   2      41.543  15.426  10.857  1.00 26.95           O
ATOM     16  N   LEU A   3      40.986  13.867  12.417  1.00 26.01           N
ATOM     17  CA  LEU A   3      40.415  12.890  11.503  1.00 25.53           C
ATOM     18  C   LEU A   3      38.998  12.690  12.000  1.00 24.10           C
ATOM     19  O   LEU A   3      38.814  12.233  13.138  1.00 24.68           O
ATOM     24  N   PHE A   4      38.029  13.007  11.125  1.00 23.52           N
ATOM     25  CA  PHE A   4      36.591  12.947  11.382  1.00 23.38           C
ATOM     26  C   PHE A   4      36.034  11.771  10.603  1.00 23.86           C
ATOM     27  O   PHE A   4      36.392  11.560   9.449  1.00 25.30           O
ATOM     35  N   VAL A   5      35.163  10.999  11.241  1.00 24.37           N
ATOM     36  CA  VAL A   5      34.553   9.829  10.616  1.00 23.47           C
ATOM     37  C   VAL A   5      33.065  10.084  10.590  1.00 22.58           C
ATOM     38  O   VAL A   5      32.478  10.418  11.618  1.00 22.07           O
ATOM     42  N   ASP A   6      32.446   9.899   9.424  1.00 21.80           N
ATOM     43  CA  ASP A   6      31.022  10.191   9.273  1.00 21.96           C
ATOM     44  C   ASP A   6      30.440   8.962   8.554  1.00 21.90           C
ATOM     45  O   ASP A   6      30.778   8.694   7.402  1.00 20.77           O
ATOM     50  N   PHE A   7      29.621   8.183   9.249  1.00 20.49           N
ATOM     51  CA  PHE A   7      29.087   6.961   8.633  1.00 23.13           C
ATOM     52  C   PHE A   7      28.154   7.314   7.469  1.00 24.00           C
ATOM     53  O   PHE A   7      27.353   8.254   7.550  1.00 22.75           O
ATOM     61  N   ASP A   8      28.212   6.533   6.393  1.00 23.11           N
ATOM     62  CA  ASP A   8      27.370   6.885   5.240  1.00 21.98           C
ATOM     63  C   ASP A   8      25.887   6.510   5.392  1.00 21.02           C
ATOM     64  O   ASP A   8      25.564   5.394   5.796  1.00 21.79           O
ATOM     69  N   TYR A   9      25.007   7.441   5.023  1.00 20.58           N
ATOM     70  CA  TYR A   9      23.558   7.272   5.091  1.00 20.92           C
ATOM     71  C   TYR A   9      23.224   6.235   6.153  1.00 22.62           C
ATOM     72  O   TYR A   9      22.422   5.300   5.958  1.00 21.24           O
ATOM     81  N   PHE A  10      23.748   6.501   7.340  1.00 20.52           N
ATOM     82  CA  PHE A  10      23.718   5.478   8.361  1.00 19.69           C
ATOM     83  C   PHE A  10      22.537   4.569   8.612  1.00 20.06           C
ATOM     84  O   PHE A  10      22.721   3.356   8.506  1.00 21.07           O
ATOM     92  N   TYR A  11      21.339   5.083   8.936  1.00 18.38           N
ATOM     93  CA  TYR A  11      20.265   4.179   9.278  1.00 20.38           C
ATOM     94  C   TYR A  11      19.883   3.292   8.080  1.00 21.64           C
ATOM     95  O   TYR A  11      19.577   2.084   8.233  1.00 20.87           O
ATOM    104  N   ALA A  12      19.887   3.898   6.896  1.00 17.70           N
ATOM    105  CA  ALA A  12      19.539   3.124   5.725  1.00 19.35           C
ATOM    106  C   ALA A  12      20.628   2.115   5.446  1.00 20.93           C
ATOM    107  O   ALA A  12      20.310   1.001   5.042  1.00 19.80           O
ATOM    109  N   GLN A  13      21.916   2.467   5.650  1.00 21.40           N
ATOM    110  CA  GLN A  13      22.973   1.495   5.407  1.00 21.95           C
ATOM    111  C   GLN A  13      22.889   0.317   6.387  1.00 21.72           C
ATOM    112  O   GLN A  13      23.142  -0.839   5.990  1.00 23.52           O
ATOM    118  N   VAL A  14      22.554   0.581   7.649  1.00 20.36           N
ATOM    119  CA  VAL A  14      22.387  -0.510   8.587  1.00 20.59           C
ATOM    120  C   VAL A  14      21.278  -1.456   8.093  1.00 23.68           C
ATOM    121  O   VAL A  14      21.447  -2.697   8.139  1.00 23.53           O
ATOM    125  N   GLU A  15      20.159  -0.892   7.629  1.00 21.78           N
ATOM    126  CA  GLU A  15      19.099  -1.754   7.132  1.00 23.52           C
ATOM    127  C   GLU A  15      19.604  -2.634   5.957  1.00 22.81           C
ATOM    128  O   GLU A  15      19.224  -3.810   5.896  1.00 24.56           O
ATOM    134  N   GLU A  16      20.481  -2.092   5.093  1.00 23.01           N
ATOM    135  CA  GLU A  16      21.056  -2.832   3.949  1.00 25.82           C
ATOM    136  C   GLU A  16      21.997  -3.923   4.460  1.00 28.76           C
ATOM    137  O   GLU A  16      22.128  -4.972   3.816  1.00 29.72           O
ATOM    143  N   VAL A  17      22.691  -3.675   5.571  1.00 27.94           N
ATOM    144  CA  VAL A  17      23.572  -4.729   6.072  1.00 28.79           C
ATOM    145  C   VAL A  17      22.731  -5.870   6.647  1.00 27.56           C
ATOM    146  O   VAL A  17      23.083  -7.034   6.448  1.00 31.25           O
ATOM    150  N   LEU A  18      21.658  -5.551   7.356  1.00 27.86           N
ATOM    151  CA  LEU A  18      20.755  -6.546   7.935  1.00 27.86           C
ATOM    152  C   LEU A  18      19.855  -7.246   6.854  1.00 27.81           C
ATOM    153  O   LEU A  18      19.358  -8.353   7.092  1.00 27.83           O
ATOM    158  N   ASN A  19      19.684  -6.612   5.693  1.00 26.87           N
ATOM    159  CA  ASN A  19      18.928  -7.212   4.583  1.00 27.50           C
ATOM    160  C   ASN A  19      19.573  -6.791   3.260  1.00 27.77           C
ATOM    161  O   ASN A  19      19.141  -5.853   2.585  1.00 27.17           O
ATOM    166  N   PRO A  20      20.663  -7.485   2.868  1.00 28.99           N
ATOM    167  CA  PRO A  20      21.435  -7.239   1.642  1.00 29.83           C
ATOM    168  C   PRO A  20      20.602  -7.208   0.342  1.00 30.16           C
ATOM    169  O   PRO A  20      21.043  -6.664  -0.654  1.00 32.03           O
";
        let dir = tempfile::tempdir().unwrap();
        let pdb_path = dir.path().join("frag20.pdb");
        std::fs::write(&pdb_path, FIXTURE_PDB_20).unwrap();

        let bin_dir = concat!(env!("CARGO_MANIFEST_DIR"), "/../bin");
        let mut ec = crate::energy::EnergyConfig::from_bin_dir(bin_dir);
        ec.num_shuffles = 20;
        ec.preload().unwrap();

        let path = pdb_path.to_str().unwrap();
        let one_domain = candidate_energy_z_score(&ec, path, "A", "1-20").unwrap();
        let split_middle = candidate_energy_z_score(&ec, path, "A", "1-10 11-20").unwrap();
        let split_early = candidate_energy_z_score(&ec, path, "A", "1-5 6-20").unwrap();

        assert_ne!(one_domain, split_middle, "single-domain and two-domain scoring must differ");
        assert_ne!(split_middle, split_early, "different domain splits must score differently");
    }

    #[test]
    fn test_modal_num_domains_picks_most_frequent() {
        assert_eq!(modal_num_domains(&[2, 3, 3, 3, 4, 5]), 3);
    }

    #[test]
    fn test_modal_num_domains_ties_break_smaller() {
        assert_eq!(modal_num_domains(&[2, 2, 5, 5]), 2);
    }

    #[test]
    fn test_modal_num_domains_empty_is_zero() {
        assert_eq!(modal_num_domains(&[]), 0);
    }

    #[test]
    fn test_modal_count_distance() {
        assert_eq!(modal_count_distance(5, 3), 2.0);
        assert_eq!(modal_count_distance(3, 3), 0.0);
    }

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
