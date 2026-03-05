//! SWORD2 CLI: SWift and Optimized Recognition of protein Domains.
//!
//! Command-line interface for running the SWORD2 protein domain recognition pipeline.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;

use sword2_lib::pdb;

/// SWORD2: SWift and Optimized Recognition of protein Domains
#[derive(Parser, Debug)]
#[command(name = "sword2", version, about)]
struct Cli {
    /// PDB code to fetch and analyze (e.g., "1TIM")
    #[arg(short = 'p', long)]
    pdb_id: Option<String>,

    /// PDB chain to analyze (e.g., "A")
    #[arg(short = 'c', long)]
    chain: Option<String>,

    /// AlphaFold UniProt Accession ID
    #[arg(short = 'u', long)]
    uniprot_id: Option<String>,

    /// ESM Metagenomic Atlas MGnify ID
    #[arg(short = 'm', long)]
    mgnify_id: Option<String>,

    /// Path to a local PDB/mmCIF file
    #[arg(short = 'i', long)]
    input_file: Option<PathBuf>,

    /// Output directory
    #[arg(short = 'o', long, default_value = ".")]
    output_dir: PathBuf,

    /// Structure model number to parse (for NMR structures)
    #[arg(long, default_value = "1")]
    model: i32,

    /// Disable energy calculations
    #[arg(long)]
    disable_energies: bool,

    /// Number of threads for parallel computation
    #[arg(short = 't', long)]
    threads: Option<usize>,
}

fn main() -> Result<()> {
    // Initialize logging
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::from_default_env()
                .add_directive("sword2=info".parse().unwrap()),
        )
        .init();

    let cli = Cli::parse();

    // Validate that at least one input source is provided
    if cli.pdb_id.is_none()
        && cli.uniprot_id.is_none()
        && cli.mgnify_id.is_none()
        && cli.input_file.is_none()
    {
        anyhow::bail!(
            "Please provide an input source: --pdb-id, --uniprot-id, --mgnify-id, or --input-file"
        );
    }

    // Parse the structure
    if let Some(input_file) = &cli.input_file {
        tracing::info!("Parsing structure from {}", input_file.display());
        let structure = pdb::parse_pdb(input_file)
            .with_context(|| format!("Failed to parse {}", input_file.display()))?;

        tracing::info!("{}", structure);

        if let Some(model) = structure.first_model() {
            for chain in &model.chains {
                tracing::info!(
                    "  Chain {}: {} residues, sequence: {}",
                    chain.id,
                    chain.len(),
                    chain.get_sequence()
                );
            }
        }
    } else {
        // TODO: Implement fetching from PDB/AlphaFold/ESM
        tracing::warn!("Remote structure fetching not yet implemented");
    }

    Ok(())
}
