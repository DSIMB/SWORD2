//! Protein peeling algorithm for hierarchical domain decomposition.
//!
//! This module parses the Peeling.log output from the SWORD binary
//! and the .num file mapping renumbered residues to original numbers.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};

/// A Protein Unit (PU) identified by peeling.
#[derive(Debug, Clone)]
pub struct ProteinUnit {
    /// PU identifier.
    pub id: String,
    /// Chain ID.
    pub chain: char,
    /// Residue ranges as (start, end) pairs.
    pub segments: Vec<(i32, i32)>,
    /// Peeling level at which this PU was identified.
    pub level: usize,
}

/// Results from a single peeling level.
#[derive(Debug, Clone)]
pub struct PeelingLevel {
    /// Peeling level number (1-based).
    pub level: usize,
    /// Internal/external contact ratio.
    pub ie_ratio: f64,
    /// Internal/(internal+external) ratio.
    pub ii_plus_e_ratio: f64,
    /// R-squared value.
    pub r2: f64,
    /// Compaction Index (CI).
    pub ci: f64,
    /// Number of PUs at this level.
    pub num_pus: usize,
    /// Protein Unit boundaries as (start, end) pairs in original residue numbering.
    pub pus: Vec<(i32, i32)>,
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
/// ie_ratio  ii+e_ratio  R2  CI  N  start1 end1 start2 end2 ...
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

        let ie_ratio: f64 = fields[0].parse().unwrap_or(0.0);
        let ii_plus_e_ratio: f64 = fields[1].parse().unwrap_or(0.0);
        let r2: f64 = fields[2].parse().unwrap_or(0.0);
        let ci: f64 = fields[3].parse().unwrap_or(0.0);
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
                    pus.push((start, end));
                }
            }
        }
        pus.sort_by_key(|&(s, _)| s);

        levels.push(PeelingLevel {
            level: nb_lvl,
            ie_ratio,
            ii_plus_e_ratio,
            r2,
            ci,
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
                .map(|(i, &(start, end))| ProteinUnit {
                    id: format!("PU{}", i + 1),
                    chain: chain_id,
                    segments: vec![(start, end)],
                    level: level.level,
                })
                .collect()
        })
        .collect()
}

/// Write peeling results summary to a text file.
pub fn write_peeling_summary(
    levels: &[PeelingLevel],
    output_path: &Path,
    energies: Option<&std::collections::HashMap<(i32, i32), (Option<f64>, Option<f64>)>>,
) -> Result<()> {
    use std::io::Write;
    let mut f = fs::File::create(output_path)
        .with_context(|| format!("Cannot create {}", output_path.display()))?;

    for level in levels {
        writeln!(
            f,
            "Peeling level {}\n    Number of Protein Units: {}\n    Compaction Index: {:.2}",
            level.level, level.num_pus, level.ci
        )?;

        for &(start, end) in &level.pus {
            if let Some(energies) = energies {
                if let Some(&(_, Some(z_score))) = energies.get(&(start, end)) {
                    let aul = if z_score.abs() >= 1.0 {
                        ((1.0 - (1.0 / (z_score * z_score))) * 100.0) as i32
                    } else {
                        0
                    };
                    writeln!(
                        f,
                        "    {:>7}: AUL={:3}% Z-score={:.1}",
                        format!("{}-{}", start, end),
                        aul,
                        z_score
                    )?;
                } else {
                    writeln!(f, "    {:>7}", format!("{}-{}", start, end))?;
                }
            } else {
                writeln!(f, "    {:>7}", format!("{}-{}", start, end))?;
            }
        }
    }

    Ok(())
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
        assert_eq!(levels[0].num_pus, 2);
        assert_eq!(levels[0].pus, vec![(10, 59), (60, 109)]);
        assert_eq!(levels[1].level, 2);
        assert_eq!(levels[1].num_pus, 3);
    }
}
