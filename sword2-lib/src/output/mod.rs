//! Output formatting for SWORD2 results.
//!
//! This module handles writing domain partitioning results in text and JSON formats.
//! Ported from Python functions write_partitionings() and write_partitionings_json().

use std::io::Write;
use std::path::Path;

use anyhow::Result;
use serde::Serialize;

/// Format for machine-readable stdout output produced by `write_stdout_summary`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OutputFormat {
    #[default]
    Text,
    Tsv,
    Json,
}

impl std::str::FromStr for OutputFormat {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self> {
        match s.to_ascii_lowercase().as_str() {
            "text" => Ok(Self::Text),
            "tsv" => Ok(Self::Tsv),
            "json" => Ok(Self::Json),
            other => anyhow::bail!("Unknown output format '{}'. Valid: text, tsv, json", other),
        }
    }
}

use crate::energy::{EnergyKey, EnergyResult};
use crate::sword::SwordResults;

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
        writeln!(
            writer,
            "Partitioning {} ({} domains):",
            i + 1,
            part.num_domains
        )?;
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

/// Get quality indicator as number of bars.
fn quality_as_nb_bars(quality: &str) -> usize {
    if quality == "n/a" {
        0
    } else {
        quality.len()
    }
}

/// Write the full SWORD2 summary text file.
///
/// Ported from Python `write_partitionings()`.
pub fn write_sword_summary(
    results: &SwordResults,
    energies: &std::collections::HashMap<EnergyKey, EnergyResult>,
    disable_energies: bool,
    output_path: &Path,
) -> Result<()> {
    let mut f = std::fs::File::create(output_path)?;

    // Ambiguity index
    let nb_bars = quality_as_nb_bars(&results.ambiguity);
    writeln!(f, "Ambiguity index: {}", "*".repeat(nb_bars))?;

    for (i, part) in results.domains.iter().enumerate() {
        writeln!(f, "-----------------------")?;
        if i == 0 {
            writeln!(f, "Optimal partition")?;
        } else {
            writeln!(f, "Alternative partition {}", i)?;
        }
        let nb_bars = quality_as_nb_bars(&part.quality);
        writeln!(f, "Quality: {}", "*".repeat(nb_bars))?;
        writeln!(f, "Nb. domains: {}", part.boundaries.len())?;

        for (j, domain) in part.boundaries.iter().enumerate() {
            if !disable_energies {
                let dom_energy = energies.get(&EnergyKey::Domain(i, j));
                let (dom_aul, dom_z_str) = format_energy(dom_energy);
                writeln!(
                    f,
                    "Domain:{}       AUL={:3}% Z-score={}",
                    j + 1,
                    dom_aul,
                    dom_z_str
                )?;
            } else {
                writeln!(f, "Domain:{}", j + 1)?;
            }

            for &(start, end) in domain {
                if !disable_energies {
                    let pu_energy = energies.get(&EnergyKey::Pu(i, j, start, end));
                    let (pu_aul, pu_z_str) = format_energy(pu_energy);
                    writeln!(
                        f,
                        "    PU:{:>7} AUL={:3}% Z-score={}",
                        format!("{}-{}", start, end),
                        pu_aul,
                        pu_z_str
                    )?;
                } else {
                    writeln!(f, "    PU:{:>7}", format!("{}-{}", start, end))?;
                }
            }
        }
    }

    Ok(())
}

/// Write the full SWORD2 summary JSON file.
///
/// Ported from Python `write_partitionings_json()`.
pub fn write_sword_summary_json(
    results: &SwordResults,
    energies: &std::collections::HashMap<EnergyKey, EnergyResult>,
    disable_energies: bool,
    output_path: &Path,
) -> Result<()> {
    let mut json_results = serde_json::Map::new();

    // Ambiguity index
    let nb_bars = quality_as_nb_bars(&results.ambiguity);
    json_results.insert(
        "Ambiguity index".to_string(),
        serde_json::Value::String("*".repeat(nb_bars)),
    );

    for (i, part) in results.domains.iter().enumerate() {
        let mut alt_part_json = serde_json::Map::new();

        let partition_name = if i == 0 {
            "Optimal partition".to_string()
        } else {
            format!("Alternative partition {}", i)
        };
        alt_part_json.insert(
            "Partition".to_string(),
            serde_json::Value::String(partition_name.clone()),
        );

        let nb_bars = quality_as_nb_bars(&part.quality);
        alt_part_json.insert(
            "Quality".to_string(),
            serde_json::Value::String("*".repeat(nb_bars)),
        );
        alt_part_json.insert(
            "Nb. domains".to_string(),
            serde_json::Value::Number(serde_json::Number::from(part.boundaries.len())),
        );

        let mut domains_json = serde_json::Map::new();
        for (j, domain) in part.boundaries.iter().enumerate() {
            let mut domain_json = serde_json::Map::new();

            if !disable_energies {
                let dom_energy = energies.get(&EnergyKey::Domain(i, j));
                let (dom_aul, dom_z_str) = format_energy(dom_energy);
                domain_json.insert(
                    "AUL".to_string(),
                    serde_json::Value::Number(serde_json::Number::from(dom_aul)),
                );
                domain_json.insert("Z-score".to_string(), serde_json::Value::String(dom_z_str));
            }

            let mut pus_json = serde_json::Map::new();
            for &(start, end) in domain {
                let pu_key = format!("{}-{}", start, end);
                if !disable_energies {
                    let pu_energy = energies.get(&EnergyKey::Pu(i, j, start, end));
                    let (pu_aul, pu_z_str) = format_energy(pu_energy);
                    let mut pu_json = serde_json::Map::new();
                    pu_json.insert(
                        "AUL".to_string(),
                        serde_json::Value::Number(serde_json::Number::from(pu_aul)),
                    );
                    pu_json.insert("Z-score".to_string(), serde_json::Value::String(pu_z_str));
                    pus_json.insert(pu_key, serde_json::Value::Object(pu_json));
                } else {
                    pus_json.insert(pu_key, serde_json::Value::Object(serde_json::Map::new()));
                }
            }
            domain_json.insert("PUs".to_string(), serde_json::Value::Object(pus_json));
            domains_json.insert(
                format!("Domain {}", j + 1),
                serde_json::Value::Object(domain_json),
            );
        }
        alt_part_json.insert(
            "Domains".to_string(),
            serde_json::Value::Object(domains_json),
        );
        json_results.insert(partition_name, serde_json::Value::Object(alt_part_json));
    }

    let json_str = serde_json::to_string_pretty(&json_results)?;
    std::fs::write(output_path, json_str)?;

    Ok(())
}

/// Write a single-run result summary line to stdout in the requested format.
///
/// `Text` is a no-op — the spinner/finish line already handles human output.
/// `Tsv` emits one tab-separated data line (print header separately via `write_tsv_header`).
/// `Json` emits one compact JSON object.
pub fn write_stdout_summary(
    id_chain: &str,
    chain: char,
    results: &crate::sword::SwordResults,
    energies: &std::collections::HashMap<EnergyKey, EnergyResult>,
    format: OutputFormat,
) {
    if format == OutputFormat::Text {
        return;
    }

    let best = match results.domains.first() {
        Some(p) => p,
        None => return,
    };

    let n_domains = best.boundaries.len();

    // "1-34,227-256;35-91;92-144" — comma separates PU segments, semicolon separates domains
    let domains_str: String = best
        .boundaries
        .iter()
        .map(|segs| {
            segs.iter()
                .map(|(s, e)| format!("{}-{}", s, e))
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect::<Vec<_>>()
        .join(";");

    let quality = &best.quality;
    let ambiguity = &results.ambiguity;

    let dom_energy = energies.get(&EnergyKey::Domain(0, 0));
    let energy_val = dom_energy.and_then(|e| e.energy);
    let zscore_val = dom_energy.and_then(|e| e.z_score);

    match format {
        OutputFormat::Tsv => {
            let energy_str = energy_val.map(|v| format!("{:.4}", v)).unwrap_or_else(|| "NA".to_string());
            let zscore_str = zscore_val.map(|v| format!("{:.2}", v)).unwrap_or_else(|| "NA".to_string());
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                id_chain, chain, n_domains, domains_str, quality, ambiguity, energy_str, zscore_str
            );
        }
        OutputFormat::Json => {
            let energy_json = energy_val
                .map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "null".to_string());
            let zscore_json = zscore_val
                .map(|v| format!("{:.2}", v))
                .unwrap_or_else(|| "null".to_string());
            println!(
                "{{\"id\":\"{}\",\"chain\":\"{}\",\"n_domains\":{},\"domains\":\"{}\",\"quality\":\"{}\",\"ambiguity\":\"{}\",\"energy\":{},\"z_score\":{}}}",
                id_chain, chain, n_domains, domains_str, quality, ambiguity, energy_json, zscore_json,
            );
        }
        OutputFormat::Text => {}
    }
}

/// Print the TSV header line. Call once before processing any entries.
pub fn write_tsv_header() {
    println!("id\tchain\tn_domains\tdomains\tquality\tambiguity\tenergy\tz_score");
}

/// Format energy result into (AUL%, Z-score string).
fn format_energy(energy: Option<&EnergyResult>) -> (i32, String) {
    match energy {
        Some(e) => {
            let aul = e.aul_percent();
            let z_str = match e.z_score {
                Some(z) => format!("{:.1}", z),
                None => "n/a".to_string(),
            };
            (aul, z_str)
        }
        None => (0, "n/a".to_string()),
    }
}
