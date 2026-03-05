//! SWORD algorithm orchestration.
//!
//! This module coordinates the overall SWORD2 pipeline:
//! 1. Structure input and validation
//! 2. Running the SWORD binary
//! 3. Parsing results
//! 4. Energy calculations
//! 5. Output generation

use std::path::Path;
use std::process::Command;

use anyhow::{Context, Result};

use crate::output::Partitioning;

/// Configuration for a SWORD2 run.
#[derive(Debug, Clone)]
pub struct SwordConfig {
    /// Path to the SWORD binary.
    pub sword_bin: String,
    /// Path to the energy binary (mypmfs).
    pub energy_bin: Option<String>,
    /// Whether to compute energies.
    pub compute_energies: bool,
    /// Number of threads for parallel computation.
    pub num_threads: usize,
    /// Output directory.
    pub output_dir: String,
}

impl Default for SwordConfig {
    fn default() -> Self {
        Self {
            sword_bin: "SWORD".to_string(),
            energy_bin: None,
            compute_energies: true,
            num_threads: num_cpus::get(),
            output_dir: ".".to_string(),
        }
    }
}

/// Run the SWORD binary on a PDB file and return the raw output.
pub fn run_sword_binary(pdb_path: &Path, config: &SwordConfig) -> Result<String> {
    let output = Command::new(&config.sword_bin)
        .arg(pdb_path)
        .output()
        .with_context(|| format!("Failed to execute SWORD binary: {}", config.sword_bin))?;

    if !output.status.success() {
        anyhow::bail!(
            "SWORD binary exited with status {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        );
    }

    Ok(String::from_utf8_lossy(&output.stdout).to_string())
}

/// Parse the output from the SWORD binary into partitionings.
///
/// Ported from Python parse_sword() function.
pub fn parse_sword_output(_raw_output: &str) -> Result<Vec<Partitioning>> {
    // TODO: Port from Python parse_sword()
    todo!("SWORD output parsing not yet implemented")
}
