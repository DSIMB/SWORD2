//! SWORD2 CLI: SWift and Optimized Recognition of protein Domains.
//!
//! Command-line interface for running the SWORD2 protein domain recognition pipeline.

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;

use sword2_lib::{energy, fetch, output, pdb, peeling, plot, sword};

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
    #[arg(short = 'e', long)]
    disable_energies: bool,

    /// Disable generation of contact probability matrix plots
    #[arg(short = 'l', long)]
    disable_plots: bool,

    /// Number of threads for parallel computation (0 = all CPUs)
    #[arg(short = 'x', long, default_value = "0")]
    cpu: usize,
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

    // Determine number of threads
    let num_threads = if cli.cpu == 0 {
        num_cpus::get()
    } else {
        cli.cpu
    };

    // Ensure output directory exists
    std::fs::create_dir_all(&cli.output_dir)?;

    // Step 1: Obtain the structure file
    let (input_path, pdb_id_base) = resolve_input(&cli)?;

    // Step 2: Parse the structure
    tracing::info!("Parsing structure from {}", input_path.display());
    let structure = pdb::parse_pdb(&input_path)
        .with_context(|| format!("Failed to parse {}", input_path.display()))?;
    tracing::info!("{}", structure);

    // Determine chain
    let chain_id = if let Some(ref chain_str) = cli.chain {
        chain_str.chars().next().unwrap_or('A')
    } else if let Some(model) = structure.first_model() {
        if let Some(chain) = model.chains.first() {
            tracing::info!("No chain specified. Using first chain '{}'", chain.id);
            chain.id
        } else {
            anyhow::bail!("No chains found in the structure");
        }
    } else {
        anyhow::bail!("No models found in the structure");
    };

    // Verify chain exists
    if let Some(model) = structure.first_model() {
        if model.get_chain(chain_id).is_none() {
            let available: Vec<String> =
                model.chains.iter().map(|c| c.id.to_string()).collect();
            anyhow::bail!(
                "Chain '{}' not found. Available chains: {}",
                chain_id,
                available.join(", ")
            );
        }
        let chain = model.get_chain(chain_id).unwrap();
        tracing::info!(
            "Chain {}: {} residues, sequence: {}",
            chain.id,
            chain.len(),
            chain.get_sequence()
        );
    }

    let pdb_id_chain = format!("{}_{}", pdb_id_base, chain_id);
    let results_dir = cli.output_dir.join(&pdb_id_chain);
    std::fs::create_dir_all(&results_dir)?;

    tracing::info!("Results directory: {}", results_dir.display());
    tracing::info!("Using {} CPUs", num_threads);

    // Step 3: Run SWORD binary
    let exe_dir = std::env::current_exe()
        .ok()
        .and_then(|p| p.parent().map(|p| p.to_path_buf()))
        .unwrap_or_else(|| PathBuf::from("."));

    let bin_dir = exe_dir.join("bin");
    let sword_bin = bin_dir.join("SWORD/bin/SWORD/SWORD");
    let display_script = bin_dir.join("display_SWORD2_output.pl");

    let config = sword::SwordConfig {
        sword_bin: sword_bin.to_string_lossy().to_string(),
        display_script: if display_script.exists() {
            Some(display_script.to_string_lossy().to_string())
        } else {
            None
        },
        energy_bin: None,
        compute_energies: !cli.disable_energies,
        generate_plots: !cli.disable_plots,
        num_threads,
        output_dir: results_dir.to_string_lossy().to_string(),
    };

    tracing::info!("Launching SWORD");
    let sword_output = sword::run_sword_binary(&input_path, &config)
        .context("Failed to run SWORD binary")?;

    // Save raw SWORD output
    std::fs::write(
        results_dir.join("sword.txt"),
        sword_output.join("\n"),
    )?;

    // Step 4: Parse SWORD output
    let sword_results = sword::parse_sword_output(&sword_output)
        .context("Failed to parse SWORD output")?;

    tracing::info!(
        "Found {} partitioning(s), ambiguity: {}",
        sword_results.domains.len(),
        sword_results.ambiguity
    );

    // Step 5: Calculate energies
    let energies: HashMap<energy::EnergyKey, energy::EnergyResult> =
        if !cli.disable_energies {
            tracing::info!("Calculating pseudo-energies of Domains");
            let energy_config =
                energy::EnergyConfig::from_bin_dir(&bin_dir.to_string_lossy());
            energy::calculate_all_energies(
                &energy_config,
                &results_dir.join(&pdb_id_chain).to_string_lossy(),
                &chain_id.to_string(),
                &sword_results.domains,
            )
        } else {
            HashMap::new()
        };

    // Step 6: Write SWORD partitionings
    output::write_sword_summary(
        &sword_results,
        &energies,
        cli.disable_energies,
        &results_dir.join("SWORD2_summary.txt"),
    )?;
    output::write_sword_summary_json(
        &sword_results,
        &energies,
        cli.disable_energies,
        &results_dir.join("SWORD2_summary.json"),
    )?;

    // Step 7: Write Peeling results
    let peeling_num = results_dir
        .join("PDBs_Clean")
        .join(&pdb_id_chain)
        .join(format!("{}.num", pdb_id_chain));
    let peeling_log = results_dir
        .join("PDBs_Clean")
        .join(&pdb_id_chain)
        .join("Peeling")
        .join("Peeling.log");

    if peeling_log.exists() {
        let ori_resnums = if peeling_num.exists() {
            peeling::parse_num_file(&peeling_num)?
        } else {
            // Default: 1..N identity mapping
            Vec::new()
        };

        if !ori_resnums.is_empty() {
            let peeling_levels = peeling::parse_peeling_log(&peeling_log, &ori_resnums)?;
            peeling::write_peeling_summary(
                &peeling_levels,
                &results_dir.join("PEELING_summary.txt"),
                None,
            )?;
            tracing::info!("Wrote peeling summary");
        }
    }

    // Step 8: Generate plots
    if !cli.disable_plots {
        let contact_matrix_dir = results_dir.join("Contact_Probability_Matrix");
        std::fs::create_dir_all(&contact_matrix_dir)?;

        let proba_mat_file = results_dir
            .join("PDBs_Clean")
            .join(&pdb_id_chain)
            .join("file_proba_contact.mat");

        if proba_mat_file.exists() {
            tracing::info!("Generating contact probability matrices");
            let mat = plot::load_contact_matrix(&proba_mat_file)?;

            for (i, part) in sword_results.domains.iter().enumerate() {
                // All PUs for this partition
                let all_pus: Vec<(i32, i32)> = part
                    .boundaries
                    .iter()
                    .flatten()
                    .copied()
                    .collect();

                let title = if i == 0 {
                    "Contact Probability Map of the\noptimal partition (all Protein Units)"
                        .to_string()
                } else {
                    format!(
                        "Contact Probability Map of the alternative\npartition n.{} (all Protein Units)",
                        i
                    )
                };

                plot::write_contact_matrix(
                    &mat,
                    contact_matrix_dir
                        .join(format!("contact_probability_matrix_alternative_{}.svg", i))
                        .to_str()
                        .unwrap(),
                    &title,
                    &all_pus,
                )?;

                // Per-domain plots
                for (j, domain) in part.boundaries.iter().enumerate() {
                    let dom_title = if i == 0 {
                        format!(
                            "Contact Probability Map of the domain {}\nof the optimal partition",
                            j + 1
                        )
                    } else {
                        format!(
                            "Contact Probability Map of the domain {}\nof the alternative partition n.{}",
                            j + 1, i
                        )
                    };

                    plot::write_contact_matrix(
                        &mat,
                        contact_matrix_dir
                            .join(format!(
                                "contact_probability_matrix_alternative_{}_domain_{}.svg",
                                i, j
                            ))
                            .to_str()
                            .unwrap(),
                        &dom_title,
                        domain,
                    )?;

                    // Per-PU plots
                    for &(start, end) in domain {
                        let pu_title = if i == 0 {
                            format!(
                                "Contact Probability Map of PU {}-{} of the domain {}\nof the optimal partition",
                                start, end, j + 1
                            )
                        } else {
                            format!(
                                "Contact Probability Map of PU {}-{} of the domain {}\nof the alternative partition n.{}",
                                start, end, j + 1, i
                            )
                        };

                        plot::write_contact_matrix(
                            &mat,
                            contact_matrix_dir
                                .join(format!(
                                    "contact_probability_matrix_alternative_{}_domain_{}_pu_{}_{}.svg",
                                    i, j, start, end
                                ))
                                .to_str()
                                .unwrap(),
                            &pu_title,
                            &[(start, end)],
                        )?;
                    }
                }
            }
        }

        // Domain histogram
        let domain_counts = plot::count_domains(&sword_results.domains);
        plot::write_domain_histogram(
            &domain_counts,
            results_dir
                .join("domains_histogram.svg")
                .to_str()
                .unwrap(),
        )?;
        tracing::info!("Generated domain histogram");
    }

    tracing::info!("Results can be found here: {}", results_dir.display());
    Ok(())
}

/// Resolve the input source to a file path and base name.
fn resolve_input(cli: &Cli) -> Result<(PathBuf, String)> {
    if let Some(ref input_file) = cli.input_file {
        let base = input_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        Ok((input_file.clone(), base))
    } else if let Some(ref uniprot_id) = cli.uniprot_id {
        let path = fetch::fetch_alphafold(uniprot_id, &cli.output_dir)?;
        Ok((path, uniprot_id.clone()))
    } else if let Some(ref mgnify_id) = cli.mgnify_id {
        let path = fetch::fetch_esm(mgnify_id, &cli.output_dir)?;
        Ok((path, mgnify_id.clone()))
    } else if let Some(ref pdb_id) = cli.pdb_id {
        let path = fetch::fetch_pdb(pdb_id, &cli.output_dir)?;
        Ok((path, pdb_id.to_uppercase()))
    } else {
        anyhow::bail!("No input source specified");
    }
}
