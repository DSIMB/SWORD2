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
use regex::Regex;

use crate::output::Partitioning;

/// Configuration for a SWORD2 run.
#[derive(Debug, Clone)]
pub struct SwordConfig {
    /// Path to the SWORD binary.
    pub sword_bin: String,
    /// Path to the display script (display_SWORD2_output.pl).
    pub display_script: Option<String>,
    /// Path to the energy binary (mypmfs).
    pub energy_bin: Option<String>,
    /// Whether to compute energies.
    pub compute_energies: bool,
    /// Whether to generate plots.
    pub generate_plots: bool,
    /// Number of threads for parallel computation.
    pub num_threads: usize,
    /// Output directory.
    pub output_dir: String,
}

impl Default for SwordConfig {
    fn default() -> Self {
        Self {
            sword_bin: "SWORD".to_string(),
            display_script: None,
            energy_bin: None,
            compute_energies: true,
            generate_plots: true,
            num_threads: num_cpus::get(),
            output_dir: ".".to_string(),
        }
    }
}

/// Raw parsed results from the SWORD binary output.
#[derive(Debug, Clone)]
pub struct SwordResults {
    /// Ambiguity index (e.g., "+++", "n/a").
    pub ambiguity: String,
    /// Domain partitionings keyed by alternative index (0 = optimal).
    pub domains: Vec<SwordPartition>,
}

/// A single SWORD partition (optimal or alternative).
#[derive(Debug, Clone)]
pub struct SwordPartition {
    /// Number of domains declared by SWORD.
    pub nb_domains: usize,
    /// Minimum domain size.
    pub min_size: usize,
    /// Domain boundaries: each domain is a list of (start, end) PU segments.
    pub boundaries: Vec<Vec<(i32, i32)>>,
    /// Average K value.
    pub average_k: f64,
    /// Quality indicator (e.g., "*****").
    pub quality: String,
}

/// Run the SWORD binary on a PDB file and return the raw output lines.
pub fn run_sword_binary(pdb_path: &Path, config: &SwordConfig) -> Result<Vec<String>> {
    let output = if let Some(display_script) = &config.display_script {
        // Use the display script wrapper: display_SWORD2_output.pl 'SWORD -i ...'
        let sword_cmd = format!(
            "{} -i {} --dir {} -max 9 -nbcpu {}",
            config.sword_bin,
            pdb_path.display(),
            config.output_dir,
            config.num_threads
        );
        Command::new(display_script)
            .arg(&sword_cmd)
            .output()
            .with_context(|| {
                format!(
                    "Failed to execute display script: {}",
                    display_script
                )
            })?
    } else {
        Command::new(&config.sword_bin)
            .arg(pdb_path)
            .output()
            .with_context(|| {
                format!("Failed to execute SWORD binary: {}", config.sword_bin)
            })?
    };

    if !output.status.success() {
        anyhow::bail!(
            "SWORD binary exited with status {}: {}",
            output.status,
            String::from_utf8_lossy(&output.stderr)
        );
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    Ok(stdout.lines().map(|s| s.to_string()).collect())
}

/// Parse the output from the SWORD binary into structured results.
///
/// Ported from Python `parse_sword()` function.
///
/// The SWORD output format has lines like:
/// ```text
/// A-index = +++
///   2 | 30 | 1-100 101-200 | 3.5 | *****
/// ```
/// Where fields separated by `|` are:
///   nb_domains | min_size | boundaries | average_k | quality
pub fn parse_sword_output(output: &[String]) -> Result<SwordResults> {
    let amb_re = Regex::new(r"^A-index = (\++)$")?;
    let ass_re = Regex::new(r"\d{1,}\s+\|")?;

    let mut ambiguity = "n/a".to_string();
    let mut domains = Vec::new();

    for line in output {
        // Check for ambiguity index
        if let Some(caps) = amb_re.captures(line) {
            ambiguity = caps[1].to_string();
            continue;
        }

        // Check for a domain assignment line
        if ass_re.is_match(line) {
            let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
            if parts.len() < 5 {
                continue;
            }

            let nb_domains: usize = parts[0].trim().parse().unwrap_or(0);
            let min_size: usize = parts[1].trim().parse().unwrap_or(0);

            // Parse boundaries: "1-100 101-200;201-250" etc.
            let boundary_strs: Vec<&str> = parts[2].split_whitespace().collect();
            let mut boundaries = Vec::new();
            for boundary in &boundary_strs {
                // Each domain can have multiple PU segments separated by ';'
                let segments: Vec<(i32, i32)> = boundary
                    .split(';')
                    .filter_map(|seg| {
                        let range: Vec<&str> = seg.split('-').collect();
                        if range.len() == 2 {
                            let start = range[0].parse::<i32>().ok()?;
                            let end = range[1].parse::<i32>().ok()?;
                            Some((start, end))
                        } else {
                            None
                        }
                    })
                    .collect();
                if !segments.is_empty() {
                    boundaries.push(segments);
                }
            }

            let average_k: f64 = parts[3].trim().parse().unwrap_or(0.0);
            let quality = parts[4].trim().to_string();

            domains.push(SwordPartition {
                nb_domains,
                min_size,
                boundaries,
                average_k,
                quality,
            });
        }
    }

    Ok(SwordResults { ambiguity, domains })
}

/// Convert SwordResults into the output Partitioning format.
pub fn results_to_partitionings(results: &SwordResults) -> Vec<Partitioning> {
    results
        .domains
        .iter()
        .map(|part| {
            let mut domains = Vec::new();
            for (i, boundary) in part.boundaries.iter().enumerate() {
                for &(start, end) in boundary {
                    domains.push(crate::output::Domain {
                        id: format!("d{}", i + 1),
                        chain: ' ',
                        start,
                        end,
                    });
                }
            }
            Partitioning {
                num_domains: part.boundaries.len(),
                domains,
                energy: None,
                z_score: None,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sword_output_basic() {
        let output = vec![
            "Some header line".to_string(),
            "A-index = +++".to_string(),
            "  2 | 30 | 1-100 101-200 | 3.5 | *****".to_string(),
            "  3 | 20 | 1-50 51-150 151-200 | 2.1 | ***".to_string(),
        ];
        let results = parse_sword_output(&output).unwrap();
        assert_eq!(results.ambiguity, "+++");
        assert_eq!(results.domains.len(), 2);

        assert_eq!(results.domains[0].nb_domains, 2);
        assert_eq!(results.domains[0].boundaries.len(), 2);
        assert_eq!(results.domains[0].boundaries[0], vec![(1, 100)]);
        assert_eq!(results.domains[0].boundaries[1], vec![(101, 200)]);
        assert_eq!(results.domains[0].quality, "*****");

        assert_eq!(results.domains[1].nb_domains, 3);
        assert_eq!(results.domains[1].boundaries.len(), 3);
    }

    #[test]
    fn test_parse_sword_output_multiseg() {
        let output = vec![
            "  2 | 30 | 1-50;151-200 51-150 | 3.5 | ****".to_string(),
        ];
        let results = parse_sword_output(&output).unwrap();
        assert_eq!(results.domains.len(), 1);
        // First domain has two segments
        assert_eq!(results.domains[0].boundaries[0], vec![(1, 50), (151, 200)]);
        // Second domain has one segment
        assert_eq!(results.domains[0].boundaries[1], vec![(51, 150)]);
    }

    #[test]
    fn test_parse_sword_no_ambiguity() {
        let output = vec![
            "  2 | 30 | 1-100 101-200 | 3.5 | *****".to_string(),
        ];
        let results = parse_sword_output(&output).unwrap();
        assert_eq!(results.ambiguity, "n/a");
    }
}
