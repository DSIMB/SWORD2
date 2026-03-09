//! Pseudo-energy calculations for protein domains.
//!
//! This module runs the external `scoring_omp` binary (from mypmfs) to compute
//! pseudo-energy and Z-score for protein domains and protein units.

use std::process::Command;

use anyhow::{Context, Result};
use rayon::prelude::*;
use regex::Regex;

/// Result of an energy calculation for a domain or PU.
#[derive(Debug, Clone)]
pub struct EnergyResult {
    /// The pseudo-energy score.
    pub energy: Option<f64>,
    /// The Z-score (number of standard deviations from random).
    pub z_score: Option<f64>,
}

impl EnergyResult {
    /// Calculate the AUL (Autonomy Unit Level) percentage from the Z-score.
    ///
    /// AUL = (1 - 1/Z²) × 100 when |Z| >= 1, else 0.
    pub fn aul_percent(&self) -> i32 {
        match self.z_score {
            Some(z) if z.abs() >= 1.0 => ((1.0 - 1.0 / (z * z)) * 100.0) as i32,
            _ => 0,
        }
    }
}

/// Path configuration for energy calculations.
#[derive(Debug, Clone)]
pub struct EnergyConfig {
    /// Path to the scoring_omp binary.
    pub scoring_bin: String,
    /// Path to the potential directory (025_30_100_potential).
    pub potential_dir: String,
    /// Number of random shuffles for Z-score.
    pub num_shuffles: usize,
}

impl EnergyConfig {
    /// Create config from the bin directory path.
    pub fn from_bin_dir(bin_dir: &str) -> Self {
        Self {
            scoring_bin: format!("{}/mypmfs-master/scoring_omp", bin_dir),
            potential_dir: format!("{}/mypmfs-master/025_30_100_potential", bin_dir),
            num_shuffles: 2000,
        }
    }
}

/// Calculate pseudo-energy and Z-score for a set of residues.
///
/// Runs the external `scoring_omp` binary:
/// ```text
/// scoring_omp -i <pdb> -d <potential_dir> [-q <residue_list>] -z -s <num_shuffles>
/// ```
///
/// Ported from Python `get_energy_and_z_score()`.
pub fn get_energy_and_z_score(
    config: &EnergyConfig,
    pdb_path: &str,
    residue_list: Option<&str>,
) -> Result<EnergyResult> {
    let mut cmd = Command::new(&config.scoring_bin);
    cmd.arg("-i").arg(pdb_path);
    cmd.arg("-d").arg(&config.potential_dir);
    if let Some(res_list) = residue_list {
        cmd.arg("-q").arg(res_list);
    }
    cmd.arg("-z");
    cmd.arg("-s").arg(config.num_shuffles.to_string());

    let output = cmd
        .output()
        .with_context(|| format!("Failed to execute scoring binary: {}", config.scoring_bin))?;

    if !output.status.success() {
        tracing::warn!(
            "Scoring binary failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return Ok(EnergyResult {
            energy: None,
            z_score: None,
        });
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let energy_re = Regex::new(r"^Pseudo-energy = (.+)$").unwrap();
    let zscore_re = Regex::new(r"^Z-score = (.+)$").unwrap();

    let mut energy = None;
    let mut z_score = None;

    for line in stdout.lines() {
        if let Some(caps) = energy_re.captures(line) {
            energy = caps[1].trim().parse::<f64>().ok();
        }
        if let Some(caps) = zscore_re.captures(line) {
            z_score = caps[1].trim().parse::<f64>().ok();
        }
    }

    Ok(EnergyResult { energy, z_score })
}

/// Build the residue list string for a set of residue numbers and a chain.
///
/// Produces comma-separated entries like "1A,2A,3A,...".
pub fn build_residue_list(residue_range: (i32, i32), chain: &str) -> String {
    (residue_range.0..=residue_range.1)
        .map(|r| format!("{}{}", r, chain))
        .collect::<Vec<_>>()
        .join(",")
}

/// Calculate energies for all domains and PUs in a partitioning, in parallel.
///
/// Returns a map from keys to EnergyResult:
/// - Domain key: `(partition_idx, domain_idx)`
/// - PU key: `(partition_idx, domain_idx, start, end)`
///
/// Ported from Python `multiprocess_get_energy()`.
pub fn calculate_all_energies(
    config: &EnergyConfig,
    pdb_path: &str,
    chain: &str,
    partitions: &[crate::sword::SwordPartition],
) -> std::collections::HashMap<EnergyKey, EnergyResult> {
    use std::collections::HashMap;
    use std::sync::Mutex;

    let results = Mutex::new(HashMap::new());

    // Collect all work items
    let mut work_items: Vec<(usize, usize, Vec<(i32, i32)>)> = Vec::new();
    for (i, part) in partitions.iter().enumerate() {
        for (j, domain) in part.boundaries.iter().enumerate() {
            work_items.push((i, j, domain.clone()));
        }
    }

    // Process in parallel using rayon
    work_items.par_iter().for_each(|(i, j, domain)| {
        let mut dom_residues = String::new();

        // Calculate energy for each PU in this domain
        for &(start, end) in domain {
            let pu_res_list = build_residue_list((start, end), chain);
            if !dom_residues.is_empty() {
                dom_residues.push(',');
            }
            dom_residues.push_str(&pu_res_list);

            if let Ok(pu_result) =
                get_energy_and_z_score(config, pdb_path, Some(&pu_res_list))
            {
                results
                    .lock()
                    .unwrap()
                    .insert(EnergyKey::Pu(*i, *j, start, end), pu_result);
            }
        }

        // Calculate energy for the entire domain
        if let Ok(dom_result) =
            get_energy_and_z_score(config, pdb_path, Some(&dom_residues))
        {
            results
                .lock()
                .unwrap()
                .insert(EnergyKey::Domain(*i, *j), dom_result);
        }
    });

    results.into_inner().unwrap()
}

/// Key for energy results lookup.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum EnergyKey {
    /// Energy for a domain: (partition_index, domain_index).
    Domain(usize, usize),
    /// Energy for a PU: (partition_index, domain_index, start_residue, end_residue).
    Pu(usize, usize, i32, i32),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_residue_list() {
        let result = build_residue_list((1, 3), "A");
        assert_eq!(result, "1A,2A,3A");
    }

    #[test]
    fn test_aul_percent() {
        let r = EnergyResult {
            energy: Some(-10.0),
            z_score: Some(2.0),
        };
        assert_eq!(r.aul_percent(), 75);

        let r2 = EnergyResult {
            energy: Some(-1.0),
            z_score: Some(0.5),
        };
        assert_eq!(r2.aul_percent(), 0);

        let r3 = EnergyResult {
            energy: None,
            z_score: None,
        };
        assert_eq!(r3.aul_percent(), 0);
    }
}
