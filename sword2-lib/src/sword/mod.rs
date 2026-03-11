//! SWORD algorithm orchestration.
//!
//! This module coordinates the full SWORD2 pipeline entirely in Rust:
//! 1. DSSP secondary structure assignment (via dsspcmbi C binary)
//! 2. Protein Peeling (via Peeling_omp C binary)
//! 3. ComputeMeasure (PU merging — pure Rust)
//! 4. ParseMeasure + prediction model (domain selection — pure Rust)
//! 5. Quality scoring and display (pure Rust)

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::LazyLock;

use anyhow::{Context, Result};
use regex::Regex;

static DIGITS_RE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d+").unwrap());
static AMB_RE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^A-index = (\++)$").unwrap());
static ASS_RE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d{1,}\s+\|").unwrap());

use crate::output::Partitioning;

pub mod compute_jones;
pub mod compute_measure;
pub mod distance_model;
pub mod junctions;
pub mod parse_measure;

/// Configuration for a SWORD2 run.
#[derive(Debug, Clone)]
pub struct SwordConfig {
    /// Path to the Peeling_omp binary.
    pub peeling_bin: String,
    /// Whether to compute energies.
    pub compute_energies: bool,
    /// Whether to generate plots.
    pub generate_plots: bool,
    /// Number of threads for parallel computation.
    pub num_threads: usize,
    /// Output directory.
    pub output_dir: String,
    /// Max alternative assignments (3, 9, or 15).
    pub max_alternatives: usize,
}

impl Default for SwordConfig {
    fn default() -> Self {
        Self {
            peeling_bin: "Peeling_omp".to_string(),
            compute_energies: true,
            generate_plots: true,
            num_threads: num_cpus::get(),
            output_dir: ".".to_string(),
            max_alternatives: 9,
        }
    }
}

/// Raw parsed results from the SWORD pipeline.
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

/// Run the complete SWORD pipeline (pure Rust, no Perl).
///
/// `pdb_path` - path to the cleaned PDB file (without .pdb extension).
/// `config` - pipeline configuration.
///
/// Returns raw output lines (same format as old SWORD binary) and the
/// structured SwordResults.
pub fn run_pipeline(
    pdb_path: &Path,
    config: &SwordConfig,
) -> Result<(Vec<String>, SwordResults)> {
    let pdb_name = pdb_path
        .file_name()
        .and_then(|f| f.to_str())
        .unwrap_or("unknown");

    let results_dir = PathBuf::from(&config.output_dir);
    let clean_dir = results_dir.join("PDBs_Clean").join(pdb_name);
    std::fs::create_dir_all(&clean_dir)?;

    // Copy PDB file to PDBs_Clean structure (SWORD convention)
    let pdb_file_src = pdb_path.to_path_buf();
    let pdb_file_dst = clean_dir.join(format!("{}.pdb", pdb_name));
    if pdb_file_src.exists() && !pdb_file_dst.exists() {
        std::fs::copy(&pdb_file_src, &pdb_file_dst)?;
    }

    // Step 1: Run DSSP (pure Rust)
    let dssp_file = clean_dir.join(format!("{}.dssp", pdb_name));
    if !dssp_file.exists() {
        tracing::debug!("Running DSSP on {}", pdb_file_dst.display());
        let s2d_file = clean_dir.join(format!("{}.s2d", pdb_name));
        crate::dssp::run_dssp(&pdb_file_dst, &dssp_file, &s2d_file, pdb_name)?;
    }

    // Step 2: Run Peeling
    let pu_delineation = clean_dir.join("file_pu_delineation.mtx");
    if !pu_delineation.exists() {
        tracing::debug!("Running Peeling on {}", pdb_file_dst.display());
        let peeling_dir = clean_dir.join("Peeling");
        std::fs::create_dir_all(&peeling_dir)?;

        let peeling_opts = format!(
            "-r 98 -s 8 -l 30 -m 0 -0 6.0 -t 1.5 -o 0 -g 0 -c 0 -n 30 -O {} -C {}",
            clean_dir.display(),
            config.num_threads,
        );

        let peeling_output = Command::new(&config.peeling_bin)
            .arg("-p")
            .arg(&pdb_file_dst)
            .arg("-d")
            .arg(&dssp_file)
            .args(peeling_opts.split_whitespace())
            .output()
            .with_context(|| format!("Failed to run Peeling: {}", config.peeling_bin))?;

        // Save peeling log
        let peeling_log = String::from_utf8_lossy(&peeling_output.stdout);
        std::fs::write(peeling_dir.join("Peeling.log"), peeling_log.as_ref())?;
    }

    // Step 3: Prepare .num file (residue number mapping)
    // The .num file should already exist from main.rs clean_chain_for_sword
    // but if not, create a default one from the PDB
    let num_file = clean_dir.join(format!("{}.num", pdb_name));
    if !num_file.exists() {
        // Create a simple 1-based num file
        let pdb_content = std::fs::read_to_string(&pdb_file_dst).unwrap_or_default();
        let resnums: Vec<String> = pdb_content
            .lines()
            .filter(|l| l.starts_with("ATOM"))
            .filter(|l| l.len() >= 26)
            .filter_map(|l| {
                let atom = l[12..16].trim();
                if atom == "CA" {
                    Some(l[22..26].trim().to_string())
                } else {
                    None
                }
            })
            .collect();
        std::fs::write(&num_file, resnums.join(" "))?;
    }

    // Read num file for original residue numbering
    let num_content = std::fs::read_to_string(&num_file).unwrap_or_default();
    let tab_num: Vec<i32> = num_content
        .split_whitespace()
        .filter_map(|s| s.parse::<i32>().ok())
        .collect();

    // Step 4: Compute measures and reconstruct domains
    let contact_matrix = clean_dir.join("file_matrix_pu_contact.mtx");
    if !pu_delineation.exists() {
        // No peeling result → single domain
        tracing::debug!("No peeling for chain, treating as single domain");
        let first = tab_num.first().copied().unwrap_or(1);
        let last = tab_num.last().copied().unwrap_or(1);
        let line = format!(
            "1 | {} | {}-{} | 0.000000 | n/a |",
            tab_num.len(),
            first,
            last
        );
        let results = SwordResults {
            ambiguity: "n/a".to_string(),
            domains: vec![SwordPartition {
                nb_domains: 1,
                min_size: tab_num.len(),
                boundaries: vec![vec![(first, last)]],
                average_k: 0.0,
                quality: "n/a".to_string(),
            }],
        };
        return Ok((vec![line], results));
    }

    // Determine alt_l and alt_b from max_alternatives
    let (alt_l, alt_b): (usize, usize) = match config.max_alternatives {
        9 => (3, 3),
        15 => (5, 3),
        _ => (3, 1),
    };

    // Run ComputeMeasure
    tracing::debug!("Computing criteria for PUs merging");
    let measure_lines = compute_measure::compute_measure(
        &contact_matrix,
        &pu_delineation,
        0.0001,
    );

    let measure_strings: Vec<String> = measure_lines.iter().map(|ml| ml.to_line()).collect();

    // First ParseMeasure pass — get all relevant measures
    let relevant_measure = parse_measure::parse_measure(
        &measure_strings,
        &results_dir.join("PDBs_Clean").to_string_lossy(),
        pdb_name,
        false,
        0,
        alt_b,
        alt_l,
        true,
    );

    // Prediction model
    tracing::debug!("Predicting structural domains from {} measures", relevant_measure.len());
    let predictions = prediction_model(&relevant_measure);
    let predictions_rev: Vec<i32> = predictions.iter().rev().cloned().collect();
    let mut n_dom = predictions_rev.len() + 1;
    for (i, &pred) in predictions_rev.iter().enumerate() {
        if pred == 0 {
            n_dom = i + 1;
            break;
        }
    }

    // Second ParseMeasure pass — select assignments around predicted N_dom
    tracing::debug!("Selecting domain assignments around N_dom={}", n_dom);
    let relevant_measure2 = parse_measure::parse_measure(
        &measure_strings,
        &results_dir.join("PDBs_Clean").to_string_lossy(),
        pdb_name,
        true,
        n_dom,
        alt_b,
        alt_l,
        true,
    );

    // Find the measure line for the predicted N_dom
    let mut to_print = String::new();
    for rm in &relevant_measure2 {
        let fields: Vec<&str> = rm.split('|').collect();
        if !fields.is_empty() {
            let n: usize = fields[0].trim().parse().unwrap_or(0);
            if n == n_dom {
                to_print = rm.clone();
                break;
            }
        }
    }

    // Quality and display
    let output_lines = quality_and_display(
        pdb_name,
        &tab_num,
        &to_print,
        n_dom,
        &relevant_measure2,
        alt_b,
        alt_l,
    );

    // Parse into structured results
    let results = parse_sword_output(&output_lines)?;

    Ok((output_lines, results))
}

/// Linear prediction model for optimal number of domains.
///
/// Port of `prediction_model()` from SWORD Perl script.
fn prediction_model(relevant_measure: &[String]) -> Vec<i32> {
    let mut predictions = Vec::new();
    let mut max_dom: i32 = -1;

    for idline in 0..relevant_measure.len().saturating_sub(1) {
        let fields: Vec<&str> = relevant_measure[idline].split('|').collect();
        if fields.is_empty() {
            continue;
        }

        let nd: i32 = fields[0].trim().parse().unwrap_or(0);
        if max_dom == -1 {
            max_dom = nd;
        }

        if nd == max_dom {
            max_dom -= 1;

            let cr: f64 = if fields.len() > 3 {
                fields[3].trim().parse().unwrap_or(0.0)
            } else {
                0.0
            };
            let obs_cpd: f64 = if fields.len() > 5 {
                fields[5].trim().parse().unwrap_or(0.0)
            } else {
                0.0
            };

            // Model parameters
            let diag_intercept: f64 = 2.818831;
            let diag_slope: f64 = 3.582524;
            let diag_inter_v: f64 = 0.09434462;
            let horizontal_lim: f64 = 3.166823;
            let vertical_lim: f64 = 0.231845;

            let theo_cpd = if cr <= diag_inter_v {
                horizontal_lim
            } else if cr >= vertical_lim {
                10000.0
            } else {
                cr * diag_slope + diag_intercept
            };

            if theo_cpd - obs_cpd > 0.0 {
                predictions.push(0);
            } else {
                predictions.push(1);
            }
        }
    }

    predictions
}

/// Generate quality and display output lines.
///
/// Port of `quality_and_display()` from SWORD Perl script.
fn quality_and_display(
    _pdb_name: &str,
    tab_num: &[i32],
    to_print: &str,
    n_dom: usize,
    relevant_measure: &[String],
    alt_b: usize,
    alt_l: usize,
) -> Vec<String> {
    let mut output = Vec::new();
    let mut globqual: Vec<usize> = Vec::new();

    // Parse the optimal assignment
    let to_print_fields: Vec<&str> = to_print.split('|').collect();

    if to_print_fields.len() < 7 {
        return output;
    }

    // Get quality for optimal
    let calc_print = if to_print_fields[0].trim().parse::<usize>().unwrap_or(0) != 1 {
        let cr: f64 = to_print_fields[3].trim().parse().unwrap_or(0.0);
        let cpd: f64 = to_print_fields[5].trim().parse().unwrap_or(0.0);
        let dist = distance_model::distance_model(cr, cpd, 0);
        let stars = distance_model::step_function(dist);
        let star_str = &"*****"[..stars];
        star_str.to_string()
    } else {
        "n/a".to_string()
    };

    // Remap residue numbers
    let delineation = remap_residue_numbers(to_print_fields[2].trim(), tab_num);
    let avg_k: f64 = to_print_fields[6].trim().parse().unwrap_or(0.0);

    // Output: ambiguity will be computed at the end
    // Header
    output.push(format!(
        "{:<2}|{:<3}|{:>60}|{:>12.6}|{:>10}|",
        to_print_fields[0].trim(),
        to_print_fields[1].trim(),
        delineation,
        avg_k,
        calc_print,
    ));

    if calc_print != "n/a" {
        globqual.push(calc_print.len());
    }

    // Alternatives
    if alt_b > 1 || alt_l >= 1 {
        let optimal_k = avg_k;

        for rm in relevant_measure {
            let fields: Vec<&str> = rm.split('|').collect();
            if fields.len() < 7 {
                continue;
            }

            let nd: usize = fields[0].trim().parse().unwrap_or(0);
            let this_k: f64 = fields[6].trim().parse().unwrap_or(0.0);

            if nd <= n_dom + alt_l && nd >= n_dom.saturating_sub(alt_l) {
                let alt_quality = if nd != 1 {
                    let cr: f64 = fields[3].trim().parse().unwrap_or(0.0);
                    let cpd: f64 = fields[5].trim().parse().unwrap_or(0.0);
                    let dist = distance_model::distance_model(cr, cpd, 0);
                    let stars = distance_model::step_function(dist);
                    let star_str = &"*****"[..stars];
                    star_str.to_string()
                } else {
                    "n/a".to_string()
                };

                if alt_quality != "n/a" && (this_k - optimal_k).abs() > 1e-10 {
                    let alt_del = remap_residue_numbers(fields[2].trim(), tab_num);
                    output.push(format!(
                        "{:<2}|{:<3}|{:>60}|{:>12.6}|{:>10}|",
                        fields[0].trim(),
                        fields[1].trim(),
                        alt_del,
                        this_k,
                        alt_quality,
                    ));

                    if alt_quality != "n/a" {
                        globqual.push(alt_quality.len());
                    }
                }
            }
        }
    }

    // Compute A-index
    let cindex = compute_cindex(&globqual);

    // Prepend A-index line
    let mut final_output = vec![format!("A-index = {}", cindex)];
    final_output.extend(output);

    final_output
}

/// Remap renumbered residue positions in a delineation string to original numbering.
fn remap_residue_numbers(delineation: &str, tab_num: &[i32]) -> String {
    DIGITS_RE.replace_all(delineation, |caps: &regex::Captures| {
        let idx: usize = caps[0].parse().unwrap_or(0);
        if idx > 0 && idx <= tab_num.len() {
            tab_num[idx - 1].to_string()
        } else if idx < tab_num.len() {
            tab_num[idx].to_string()
        } else {
            caps[0].to_string()
        }
    })
    .to_string()
}

/// Compute the complexity (ambiguity) index from quality scores.
///
/// Port of `cindex()` from SWORD Perl script.
fn compute_cindex(globqual: &[usize]) -> String {
    let c5 = globqual.iter().filter(|&&q| q >= 5).count();
    let c4 = globqual.iter().filter(|&&q| q >= 4).count();
    let c3 = globqual.iter().filter(|&&q| q >= 3).count();
    let c2 = globqual.iter().filter(|&&q| q >= 2).count();
    let c1 = globqual.iter().filter(|&&q| q >= 1).count();

    if c5 >= 5 {
        "+++++".to_string()
    } else if c4 >= 4 {
        "++++".to_string()
    } else if c3 >= 3 {
        "+++".to_string()
    } else if c2 >= 2 {
        "++".to_string()
    } else if c1 >= 1 {
        "+".to_string()
    } else {
        "+".to_string()
    }
}


/// Parse the output from the SWORD pipeline into structured results.
///
/// The output format has lines like:
/// ```text
/// A-index = +++
///   2 | 30 | 1-100 101-200 | 3.5 | *****
/// ```
pub fn parse_sword_output(output: &[String]) -> Result<SwordResults> {
    let mut ambiguity = "n/a".to_string();
    let mut domains = Vec::new();

    for line in output {
        // Check for ambiguity index
        if let Some(caps) = AMB_RE.captures(line) {
            ambiguity = caps[1].to_string();
            continue;
        }

        // Check for a domain assignment line
        if ASS_RE.is_match(line) {
            let parts: Vec<&str> = line.split('|').map(|s| s.trim()).collect();
            if parts.len() < 5 {
                continue;
            }

            let nb_domains: usize = parts[0].trim().parse().unwrap_or(0);
            let min_size: usize = parts[1].trim().parse().unwrap_or(0);

            // Parse boundaries: "1-100 101-200" etc.
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
        assert_eq!(results.domains[0].boundaries[0], vec![(1, 50), (151, 200)]);
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

    #[test]
    fn test_cindex() {
        assert_eq!(compute_cindex(&[5, 5, 5, 5, 5]), "+++++");
        assert_eq!(compute_cindex(&[4, 4, 4, 4]), "++++");
        assert_eq!(compute_cindex(&[3, 3, 3]), "+++");
        assert_eq!(compute_cindex(&[2, 2]), "++");
        assert_eq!(compute_cindex(&[1]), "+");
    }

    #[test]
    fn test_prediction_model() {
        // With CR below threshold and obs_CPD below theo → prediction = 0
        let lines = vec![
            "5|30|...|0.05|...|2.5|0.5".to_string(),
            "4|30|...|0.05|...|2.5|0.5".to_string(),
            "".to_string(),
        ];
        let preds = prediction_model(&lines);
        assert!(!preds.is_empty());
    }
}
