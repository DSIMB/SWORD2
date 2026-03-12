//! Protein peeling algorithm for hierarchical domain decomposition.
//!
//! This module implements and orchestrates the Protein Peeling algorithm
//! (Gelly et al., 2006) for iterative domain decomposition. The native
//! Rust implementation replaces the external `Peeling_omp` C binary.

pub mod algorithm;
pub mod contact_matrix;

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

pub use algorithm::{PeelingConfig, PeelingOutput, run_peeling};
pub use contact_matrix::ContactMatrix;

/// The backend that produced peeling results.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum PeelingBackend {
    /// Legacy external Peeling_omp binary.
    LegacyBinary,
    /// Future native Rust implementation.
    NativeRust,
}

/// A residue interval in original residue numbering.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct ResidueRange {
    /// First residue in the interval.
    pub start: i32,
    /// Last residue in the interval.
    pub end: i32,
}

impl ResidueRange {
    /// Create a new residue range.
    pub fn new(start: i32, end: i32) -> Self {
        Self { start, end }
    }
}

/// A Protein Unit (PU) identified by peeling.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinUnit {
    /// PU identifier.
    pub id: String,
    /// Chain ID.
    pub chain: char,
    /// Residue ranges.
    pub segments: Vec<ResidueRange>,
    /// Peeling level at which this PU was identified.
    pub level: usize,
}

/// Results from a single peeling level.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeelingLevel {
    /// Peeling level number (1-based).
    pub level: usize,
    /// Maximum contact ratio (Max_CR in legacy output).
    pub max_cr: f64,
    /// Minimum density (Min_Density in legacy output).
    pub min_density: f64,
    /// Compaction Index (CI).
    pub ci: f64,
    /// R metric from the legacy output.
    pub r: f64,
    /// Number of PUs at this level.
    pub num_pus: usize,
    /// Protein Unit boundaries in original residue numbering.
    pub pus: Vec<ResidueRange>,
}

/// Structured peeling results independent of the producing backend.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PeelingResults {
    /// Backend used to generate the results.
    pub backend: PeelingBackend,
    /// Original residue numbers for the cleaned chain.
    pub original_resnums: Vec<i32>,
    /// All peeling levels in ascending order.
    pub levels: Vec<PeelingLevel>,
}

/// Parse the .num file that maps renumbered residues to original numbers.
///
/// The .num file contains a single line with space-separated original residue numbers.
pub fn parse_num_file(num_path: &Path) -> Result<Vec<i32>> {
    let content = fs::read_to_string(num_path)
        .with_context(|| format!("Cannot read num file: {}", num_path.display()))?;

    let resnums: Vec<i32> = content
        .split_whitespace()
        .filter_map(|s| s.parse::<i32>().ok())
        .collect();

    Ok(resnums)
}

/// Parse the Peeling.log file to extract peeling levels and PU boundaries.
///
/// The log file format (after skipping the header line):
/// ```text
/// # header line
/// Max_CR  Min_Density  CI  R  N  start1 end1 start2 end2 ...
/// ```
///
/// The residue numbers in the log are 1-based indices into the renumbered sequence.
/// They are converted back to original numbering using `ori_resnums`.
pub fn parse_peeling_log(
    peeling_log: &Path,
    ori_resnums: &[i32],
) -> Result<Vec<PeelingLevel>> {
    let content = fs::read_to_string(peeling_log)
        .with_context(|| format!("Cannot read peeling log: {}", peeling_log.display()))?;

    let mut levels = Vec::new();
    let mut nb_lvl: usize = 1;

    for (i, line) in content.lines().enumerate() {
        // Skip the first line (header)
        if i == 0 {
            continue;
        }
        // Skip comment lines and empty lines
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 6 {
            continue;
        }

        let max_cr: f64 = fields[0].parse().unwrap_or(0.0);
        let min_density: f64 = fields[1].parse().unwrap_or(0.0);
        let ci: f64 = fields[2].parse().unwrap_or(0.0);
        let r: f64 = fields[3].parse().unwrap_or(0.0);
        let num_pus: usize = fields[4].parse().unwrap_or(0);

        // Parse PU boundaries from remaining fields
        // Fields 5+ are pairs: start1 end1 start2 end2 ...
        let boundary_fields = &fields[5..];
        let mut pus = Vec::new();
        for chunk in boundary_fields.chunks(2) {
            if chunk.len() == 2 {
                let start_idx: usize = chunk[0].parse().unwrap_or(0);
                let end_idx: usize = chunk[1].parse().unwrap_or(0);
                // Convert 1-based index to original residue numbering
                if start_idx > 0
                    && end_idx > 0
                    && start_idx <= ori_resnums.len()
                    && end_idx <= ori_resnums.len()
                {
                    let start = ori_resnums[start_idx - 1];
                    let end = ori_resnums[end_idx - 1];
                    pus.push(ResidueRange::new(start, end));
                }
            }
        }
        pus.sort_by_key(|range| range.start);

        levels.push(PeelingLevel {
            level: nb_lvl,
            max_cr,
            min_density,
            ci,
            r,
            num_pus,
            pus,
        });
        nb_lvl += 1;
    }

    Ok(levels)
}

/// Convert peeling levels to a vector of ProteinUnit vectors (one per level).
pub fn levels_to_protein_units(
    levels: &[PeelingLevel],
    chain_id: char,
) -> Vec<Vec<ProteinUnit>> {
    levels
        .iter()
        .map(|level| {
            level
                .pus
                .iter()
                .enumerate()
                .map(|(i, range)| ProteinUnit {
                    id: format!("PU{}", i + 1),
                    chain: chain_id,
                    segments: vec![*range],
                    level: level.level,
                })
                .collect()
        })
        .collect()
}

/// Load legacy peeling outputs into the structured result model.
pub fn load_legacy_results(
    peeling_log: &Path,
    num_path: &Path,
) -> Result<PeelingResults> {
    let original_resnums = parse_num_file(num_path)?;
    let levels = parse_peeling_log(peeling_log, &original_resnums)?;

    Ok(PeelingResults {
        backend: PeelingBackend::LegacyBinary,
        original_resnums,
        levels,
    })
}

fn format_metric(value: f64) -> String {
    let rounded = (value * 100.0).round() / 100.0;
    let trimmed = format!("{rounded:.2}")
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string();

    if trimmed.contains('.') {
        trimmed
    } else {
        format!("{trimmed}.0")
    }
}

#[derive(Serialize)]
struct PeelingSummaryEntry {
    boundary: ResidueRange,
    aul_percent: Option<i32>,
    z_score: Option<f64>,
}

#[derive(Serialize)]
struct PeelingLevelSummary<'a> {
    level: usize,
    max_cr: f64,
    min_density: f64,
    ci: f64,
    r: f64,
    num_pus: usize,
    pus: Vec<PeelingSummaryEntry>,
    #[serde(skip_serializing_if = "Option::is_none")]
    backend: Option<&'a PeelingBackend>,
}

#[derive(Serialize)]
struct PeelingResultsSummary<'a> {
    backend: PeelingBackend,
    original_resnums: &'a [i32],
    levels: Vec<PeelingLevelSummary<'a>>,
}

fn aul_percent(z_score: f64) -> i32 {
    if z_score.abs() >= 1.0 {
        ((1.0 - (1.0 / (z_score * z_score))) * 100.0) as i32
    } else {
        0
    }
}

fn build_summary<'a>(
    results: &'a PeelingResults,
    energies: Option<&std::collections::HashMap<(i32, i32), (Option<f64>, Option<f64>)>>,
) -> PeelingResultsSummary<'a> {
    let levels = results
        .levels
        .iter()
        .enumerate()
        .map(|(index, level)| {
            let pus = level
                .pus
                .iter()
                .map(|range| {
                    let (_, z_score) = energies
                        .and_then(|values| values.get(&(range.start, range.end)).copied())
                        .unwrap_or((None, None));

                    PeelingSummaryEntry {
                        boundary: *range,
                        aul_percent: z_score.map(aul_percent),
                        z_score,
                    }
                })
                .collect();

            PeelingLevelSummary {
                level: level.level,
                max_cr: level.max_cr,
                min_density: level.min_density,
                ci: level.ci,
                r: level.r,
                num_pus: level.num_pus,
                pus,
                backend: (index == 0).then_some(&results.backend),
            }
        })
        .collect();

    PeelingResultsSummary {
        backend: results.backend,
        original_resnums: &results.original_resnums,
        levels,
    }
}

/// Write peeling results summary to a text file.
pub fn write_peeling_summary(
    results: &PeelingResults,
    output_path: &Path,
    energies: Option<&std::collections::HashMap<(i32, i32), (Option<f64>, Option<f64>)>>,
) -> Result<()> {
    use std::io::Write;
    let mut f = fs::File::create(output_path)
        .with_context(|| format!("Cannot create {}", output_path.display()))?;

    writeln!(f, "Peeling backend: {:?}", results.backend)?;

    for level in &results.levels {
        writeln!(
            f,
            "\nPeeling level {}\n  Protein Units : {}\n  Max CR        : {}\n  Min Density   : {}\n  CI            : {}\n  R             : {}\n  Boundaries",
            level.level,
            level.num_pus,
            format_metric(level.max_cr),
            format_metric(level.min_density),
            format_metric(level.ci),
            format_metric(level.r),
        )?;

        for (index, range) in level.pus.iter().enumerate() {
            if let Some(energies) = energies {
                if let Some(&(_, Some(z_score))) = energies.get(&(range.start, range.end)) {
                    writeln!(
                        f,
                        "    {:>2}. {:>7}  AUL={:3}%  Z-score={:.1}",
                        index + 1,
                        format!("{}-{}", range.start, range.end),
                        aul_percent(z_score),
                        z_score
                    )?;
                } else {
                    writeln!(
                        f,
                        "    {:>2}. {:>7}",
                        index + 1,
                        format!("{}-{}", range.start, range.end)
                    )?;
                }
            } else {
                writeln!(
                    f,
                    "    {:>2}. {:>7}",
                    index + 1,
                    format!("{}-{}", range.start, range.end)
                )?;
            }
        }
    }

    Ok(())
}

/// Write peeling results summary to JSON.
pub fn write_peeling_summary_json(
    results: &PeelingResults,
    output_path: &Path,
    energies: Option<&std::collections::HashMap<(i32, i32), (Option<f64>, Option<f64>)>>,
) -> Result<()> {
    let summary = build_summary(results, energies);
    let json = serde_json::to_string_pretty(&summary)?;
    fs::write(output_path, json)
        .with_context(|| format!("Cannot create {}", output_path.display()))?;
    Ok(())
}

/// Convert native peeling output into the structured `PeelingResults` model.
///
/// The `PeelingOutput` from the native algorithm uses 0-based indexing internally.
/// The `true_nums` field maps those indices to DSSP residue numbers (which are
/// 1-based sequential for the cleaned chain). We use those as our "original" numbers.
pub fn native_peeling_to_results(output: &PeelingOutput) -> PeelingResults {
    let levels = output
        .iterations
        .iter()
        .enumerate()
        .map(|(i, iter)| {
            let mut pus: Vec<ResidueRange> = iter
                .pu_boundaries
                .iter()
                .map(|pu| {
                    let start = output.true_nums[pu[0]];
                    let end = output.true_nums[pu[1]];
                    ResidueRange::new(start, end)
                })
                .collect();
            pus.sort_by_key(|r| r.start);

            PeelingLevel {
                level: i + 1,
                max_cr: iter.max_cr,
                min_density: iter.min_density,
                ci: iter.ci,
                r: iter.r,
                num_pus: iter.num_pus,
                pus,
            }
        })
        .collect();

    PeelingResults {
        backend: PeelingBackend::NativeRust,
        original_resnums: output.true_nums.clone(),
        levels,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_peeling_log() {
        let dir = tempfile::tempdir().unwrap();
        let log_path = dir.path().join("Peeling.log");
        fs::write(
            &log_path,
            "# header\n0.5 0.3 0.95 1.2 2 1 50 51 100\n0.4 0.2 0.90 1.5 3 1 30 31 70 71 100\n",
        )
        .unwrap();

        // Original residue numbering: 1..100 mapped to 10..109
        let ori_resnums: Vec<i32> = (10..110).collect();
        let levels = parse_peeling_log(&log_path, &ori_resnums).unwrap();

        assert_eq!(levels.len(), 2);
        assert_eq!(levels[0].level, 1);
        assert_eq!(levels[0].max_cr, 0.5);
        assert_eq!(levels[0].min_density, 0.3);
        assert_eq!(levels[0].ci, 0.95);
        assert_eq!(levels[0].r, 1.2);
        assert_eq!(levels[0].num_pus, 2);
        assert_eq!(
            levels[0].pus,
            vec![ResidueRange::new(10, 59), ResidueRange::new(60, 109)]
        );
        assert_eq!(levels[1].level, 2);
        assert_eq!(levels[1].num_pus, 3);
    }

    #[test]
    fn test_load_legacy_results() {
        let dir = tempfile::tempdir().unwrap();
        let log_path = dir.path().join("Peeling.log");
        let num_path = dir.path().join("test.num");

        fs::write(&log_path, "Max_CR Min_Density CI R Num_PUs PU_Delineations\n0.5 0.3 0.95 1.2 2 1 2 3 4\n").unwrap();
        fs::write(&num_path, "10 11 12 13").unwrap();

        let results = load_legacy_results(&log_path, &num_path).unwrap();
        assert_eq!(results.backend, PeelingBackend::LegacyBinary);
        assert_eq!(results.original_resnums, vec![10, 11, 12, 13]);
        assert_eq!(results.levels[0].pus, vec![ResidueRange::new(10, 11), ResidueRange::new(12, 13)]);
    }
}
