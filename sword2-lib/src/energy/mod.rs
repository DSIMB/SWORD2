//! Pseudo-energy calculations for protein domains.
//!
//! Pure-Rust pseudo-energy and Z-score computation (see [`score`]) for protein
//! domains and protein units, using the precomputed mypmfs potentials.

use std::sync::{Arc, OnceLock};

use anyhow::Result;
use rayon::prelude::*;

mod score;
pub use score::{score as score_structure, Potentials};

type WorkItem = (usize, usize, Vec<(i32, i32)>);

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
    /// AUL = (1 - 1/Z²) × 100 when Z <= -1, else 0.
    pub fn aul_percent(&self) -> i32 {
        match self.z_score {
            Some(z) => aul_percent_for_z(z),
            None => 0,
        }
    }
}

pub(crate) fn aul_percent_for_z(z: f64) -> i32 {
    if z <= -1.0 {
        ((1.0 - 1.0 / (z * z)) * 100.0) as i32
    } else {
        0
    }
}

/// Path configuration for energy calculations.
#[derive(Debug, Clone)]
pub struct EnergyConfig {
    /// Path to the potential directory (025_30_100_potential).
    pub potential_dir: String,
    /// Number of random shuffles for Z-score.
    pub num_shuffles: usize,
    /// Lazily-loaded, shared potentials (loaded once, reused across all calls).
    potentials: Arc<OnceLock<Potentials>>,
}

impl EnergyConfig {
    /// Create config from the bin directory path.
    pub fn from_bin_dir(bin_dir: &str) -> Self {
        Self {
            potential_dir: format!("{}/mypmfs-master/025_30_100_potential", bin_dir),
            num_shuffles: 2000,
            potentials: Arc::new(OnceLock::new()),
        }
    }

    /// Get the loaded potentials, loading them on first access.
    ///
    /// Safe to call from multiple threads: a losing racer simply discards its
    /// redundant load. Call [`EnergyConfig::preload`] once up front to avoid that
    /// race entirely before a parallel batch.
    fn potentials(&self) -> Result<&Potentials> {
        if let Some(p) = self.potentials.get() {
            return Ok(p);
        }
        let loaded = Potentials::load(&self.potential_dir)?;
        let _ = self.potentials.set(loaded);
        Ok(self.potentials.get().expect("potentials set"))
    }

    /// Eagerly load the potentials (call once before parallel scoring).
    pub fn preload(&self) -> Result<()> {
        self.potentials().map(|_| ())
    }
}

/// Calculate pseudo-energy and Z-score for a set of residues.
///
/// Pure-Rust scoring (see [`score`]), CA representation with linear interpolation.
/// `residue_list` (comma-separated `numchain` tokens) restricts the calculation to
/// a residue subset.
pub fn get_energy_and_z_score(
    config: &EnergyConfig,
    pdb_path: &str,
    residue_list: Option<&str>,
) -> Result<EnergyResult> {
    let pot = config.potentials()?;
    let result = score::score(pot, pdb_path, residue_list, config.num_shuffles, true)?;
    tracing::trace!(
        "Energy result: energy={:?}, z_score={:?}",
        result.energy,
        result.z_score
    );
    Ok(result)
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
    let mut work_items: Vec<WorkItem> = Vec::new();
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

            if let Ok(pu_result) = get_energy_and_z_score(config, pdb_path, Some(&pu_res_list)) {
                results
                    .lock()
                    .unwrap()
                    .insert(EnergyKey::Pu(*i, *j, start, end), pu_result);
            }
        }

        // Calculate energy for the entire domain
        if let Ok(dom_result) = get_energy_and_z_score(config, pdb_path, Some(&dom_residues)) {
            results
                .lock()
                .unwrap()
                .insert(EnergyKey::Domain(*i, *j), dom_result);
        }
    });

    results.into_inner().unwrap()
}

/// Calculate energies for all domains, reusing pre-computed PU energy cache.
///
/// PU-level energies are looked up from `pu_cache` instead of re-invoking
/// the external binary. Only domain-level (multi-PU) energies are computed fresh.
pub fn calculate_all_energies_with_cache(
    config: &EnergyConfig,
    pdb_path: &str,
    chain: &str,
    partitions: &[crate::sword::SwordPartition],
    pu_cache: &std::collections::HashMap<(i32, i32), EnergyResult>,
) -> std::collections::HashMap<EnergyKey, EnergyResult> {
    use std::collections::HashMap;
    use std::sync::Mutex;

    let results = Mutex::new(HashMap::new());

    // Collect all work items
    let mut work_items: Vec<WorkItem> = Vec::new();
    for (i, part) in partitions.iter().enumerate() {
        for (j, domain) in part.boundaries.iter().enumerate() {
            work_items.push((i, j, domain.clone()));
        }
    }

    // Process in parallel using rayon
    work_items.par_iter().for_each(|(i, j, domain)| {
        let mut dom_residues = String::new();

        // Use cached PU energies
        for &(start, end) in domain {
            let pu_res_list = build_residue_list((start, end), chain);
            if !dom_residues.is_empty() {
                dom_residues.push(',');
            }
            dom_residues.push_str(&pu_res_list);

            if let Some(pu_result) = pu_cache.get(&(start, end)) {
                results
                    .lock()
                    .unwrap()
                    .insert(EnergyKey::Pu(*i, *j, start, end), pu_result.clone());
            }
        }

        // Calculate energy for the entire domain (not cached — unique per partition)
        if let Ok(dom_result) = get_energy_and_z_score(config, pdb_path, Some(&dom_residues)) {
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

/// Pre-compute energies for a set of unique PU ranges in parallel.
///
/// Returns a map from (start, end) to EnergyResult. This avoids redundant
/// external binary invocations when the same PU ranges appear across
/// SWORD partitions and peeling levels.
pub fn compute_pu_energies_batch(
    config: &EnergyConfig,
    pdb_path: &str,
    chain: &str,
    pu_ranges: &[(i32, i32)],
) -> std::collections::HashMap<(i32, i32), EnergyResult> {
    tracing::debug!(
        "Computing energies for {} unique PU ranges",
        pu_ranges.len()
    );
    pu_ranges
        .par_iter()
        .filter_map(|&(start, end)| {
            let res_list = build_residue_list((start, end), chain);
            get_energy_and_z_score(config, pdb_path, Some(&res_list))
                .ok()
                .map(|r| ((start, end), r))
        })
        .collect()
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
        assert_eq!(r.aul_percent(), 0);

        let favorable = EnergyResult {
            energy: Some(-10.0),
            z_score: Some(-2.0),
        };
        assert_eq!(favorable.aul_percent(), 75);

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

    #[test]
    fn aul_formula_only_rewards_favorable_negative_z_scores() {
        assert_eq!(aul_percent_for_z(2.0), 0);
        assert_eq!(aul_percent_for_z(-0.5), 0);
        assert_eq!(aul_percent_for_z(-2.0), 75);
    }
}
