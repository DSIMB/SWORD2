//! Output formatting for SWORD2 results.
//!
//! This module handles writing domain partitioning results in text and JSON formats.
//! Ported from Python functions write_partitionings() and write_partitionings_json().

use std::io::Write;
use std::path::Path;

use anyhow::Result;
use serde::Serialize;

/// A single domain in a partitioning result.
#[derive(Debug, Clone, Serialize)]
pub struct Domain {
    /// Domain identifier (e.g., "d1", "d2").
    pub id: String,
    /// Chain ID.
    pub chain: char,
    /// Start residue number.
    pub start: i32,
    /// End residue number.
    pub end: i32,
}

/// A complete partitioning result.
#[derive(Debug, Clone, Serialize)]
pub struct Partitioning {
    /// Number of domains in this partitioning.
    pub num_domains: usize,
    /// The domains.
    pub domains: Vec<Domain>,
    /// Energy score (if calculated).
    pub energy: Option<f64>,
    /// Z-score (if calculated).
    pub z_score: Option<f64>,
}

/// Write partitioning results as JSON to a file.
pub fn write_json(partitionings: &[Partitioning], path: &Path) -> Result<()> {
    let json = serde_json::to_string_pretty(partitionings)?;
    std::fs::write(path, json)?;
    Ok(())
}

/// Write partitioning results as plain text to a writer.
pub fn write_text<W: Write>(partitionings: &[Partitioning], writer: &mut W) -> Result<()> {
    for (i, part) in partitionings.iter().enumerate() {
        writeln!(writer, "Partitioning {} ({} domains):", i + 1, part.num_domains)?;
        for domain in &part.domains {
            writeln!(
                writer,
                "  {} chain {} : {}-{}",
                domain.id, domain.chain, domain.start, domain.end
            )?;
        }
        if let Some(energy) = part.energy {
            writeln!(writer, "  Energy: {:.4}", energy)?;
        }
        if let Some(z_score) = part.z_score {
            writeln!(writer, "  Z-score: {:.4}", z_score)?;
        }
        writeln!(writer)?;
    }
    Ok(())
}
