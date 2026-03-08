//! SWORD2 CLI: SWift and Optimized Recognition of protein Domains.
//!
//! Command-line interface for running the SWORD2 protein domain recognition pipeline.
//! Replicates the Python SWORD2.py pipeline faithfully.

use std::collections::HashMap;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;

use sword2_lib::{energy, fetch, output, pdb, peeling, sword};

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

    /// Path to SWORD2 base directory (defaults to parent of binary location)
    #[arg(long)]
    base_dir: Option<PathBuf>,
}

fn main() -> Result<()> {
    let start = Instant::now();

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

    // Resolve base directory (where bin/ lives)
    let base_dir = if let Some(ref bd) = cli.base_dir {
        bd.clone()
    } else {
        // Try: parent of the binary, then current directory
        std::env::current_exe()
            .ok()
            .and_then(|p| p.parent().map(|pp| pp.to_path_buf()))
            .and_then(|p| {
                if p.join("bin").exists() {
                    Some(p)
                } else {
                    p.parent().and_then(|pp| {
                        if pp.join("bin").exists() {
                            Some(pp.to_path_buf())
                        } else {
                            None
                        }
                    })
                }
            })
            .unwrap_or_else(|| std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")))
    };

    let bin_dir = base_dir.join("bin");
    let sword_dir = bin_dir.join("SWORD/bin/SWORD");
    let sword_bin = sword_dir.join("SWORD");
    let display_script = bin_dir.join("display_SWORD2_output.pl");

    // Ensure output directory exists
    let output_dir = std::fs::canonicalize(&cli.output_dir).unwrap_or_else(|_| {
        std::fs::create_dir_all(&cli.output_dir).ok();
        cli.output_dir.clone()
    });
    std::fs::create_dir_all(&output_dir)?;

    // Step 1: Obtain the structure file
    let (input_path, pdb_id_base) = resolve_input(&cli, &output_dir)?;

    // Step 2: Parse the structure
    tracing::info!("Parsing structure from {}", input_path.display());
    let structure = pdb::parse_pdb(&input_path)
        .with_context(|| format!("Failed to parse {}", input_path.display()))?;

    // Determine chain
    let chain_id = if let Some(ref chain_str) = cli.chain {
        chain_str.chars().next().unwrap_or('A')
    } else if let Some(model) = structure.first_model() {
        // Find the first chain that contains at least one standard amino acid residue
        // (skip nucleic acid chains like DNA/RNA)
        if let Some(chain) = model.chains.iter().find(|c| {
            c.residues.iter().any(|r| pdb::amino_acids::is_standard(&r.name))
        }) {
            tracing::info!("No chain specified. Using first protein chain '{}'", chain.id);
            chain.id
        } else if let Some(chain) = model.chains.first() {
            tracing::info!("No chain specified. Using first chain '{}'", chain.id);
            chain.id
        } else {
            anyhow::bail!("No chains found in the structure");
        }
    } else {
        anyhow::bail!("No models found in the structure");
    };

    // Verify chain exists
    let model = structure
        .first_model()
        .ok_or_else(|| anyhow::anyhow!("No models found"))?;
    let chain = model
        .get_chain(chain_id)
        .ok_or_else(|| {
            let available: Vec<String> = model.chains.iter().map(|c| c.id.to_string()).collect();
            anyhow::anyhow!(
                "Chain '{}' not found. Available chains: {}",
                chain_id,
                available.join(", ")
            )
        })?;

    let pdb_id_chain = format!("{}_{}", pdb_id_base, chain_id);
    let results_dir = output_dir.join(&pdb_id_chain);
    std::fs::create_dir_all(&results_dir)?;

    tracing::info!(">>> {} ({} residues)", pdb_id_chain, chain.len());
    tracing::info!(">>> Using {} cpus", num_threads);

    // Step 3: Clean PDB - remove non-standard residues, insertion codes, renumber from 1
    // This replicates the Python: prot.select("protein and not nonstdaa and not hetatm")
    tracing::info!("Write a clean version of the PDB: remove non standard residues");
    let (cleaned_chain, original_resnums) = pdb::writer::clean_chain_for_sword(chain);

    if cleaned_chain.is_empty() {
        anyhow::bail!(
            "No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file."
        );
    }

    let prot_len = cleaned_chain.len();
    tracing::info!("Clean chain: {} residues, sequence: {}", prot_len, cleaned_chain.get_sequence());

    // Write clean PDB file
    let pdb_chain_file = results_dir.join(format!("{}.pdb", pdb_id_chain));
    pdb::write_pdb(&cleaned_chain, &pdb_chain_file)?;

    // Remove the .pdb extension for SWORD (SWORD expects file without extension)
    let pdb_no_ext = results_dir.join(&pdb_id_chain);
    std::fs::rename(&pdb_chain_file, &pdb_no_ext)?;

    // Step 4: Compile DSSP if needed (first run of SWORD)
    let dssp_path = sword_dir.join("bin/Dssp/dsspcmbi");
    if !dssp_path.exists() {
        tracing::info!("Compiling DSSP dependency (first run)");
        let _ = std::process::Command::new(&sword_bin)
            .output();
    }

    // Step 5: Run SWORD binary
    tracing::info!("Launch SWORD");
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

    let sword_output = sword::run_sword_binary(
        &pdb_no_ext, // Pass file without .pdb extension
        &config,
    ).context("Failed to run SWORD binary")?;

    // Save raw SWORD output
    std::fs::write(results_dir.join("sword.txt"), sword_output.join("\n") + "\n")?;

    // Step 6: Parse SWORD output
    tracing::info!("Parse SWORD output");
    let sword_results = sword::parse_sword_output(&sword_output)
        .context("Failed to parse SWORD output")?;

    tracing::info!(
        "Found {} partitioning(s), ambiguity: {}",
        sword_results.domains.len(),
        sword_results.ambiguity
    );

    // Step 7: Calculate energies
    let energies: HashMap<energy::EnergyKey, energy::EnergyResult> =
        if !cli.disable_energies {
            tracing::info!("Calculate pseudo-energies of Domains");
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

    // Step 8: Write SWORD partitionings (text + JSON)
    tracing::info!("Write the SWORD results");
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

    // Step 9: Write Peeling results
    tracing::info!("Write Peeling results");
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
            Vec::new()
        };

        if !ori_resnums.is_empty() {
            let peeling_levels = peeling::parse_peeling_log(&peeling_log, &ori_resnums)?;

            // Calculate peeling energies if enabled
            let peeling_energies = if !cli.disable_energies {
                let energy_config =
                    energy::EnergyConfig::from_bin_dir(&bin_dir.to_string_lossy());
                let pdb_path_str = results_dir.join(&pdb_id_chain).to_string_lossy().to_string();
                let mut pe: HashMap<(i32, i32), (Option<f64>, Option<f64>)> = HashMap::new();
                for level in &peeling_levels {
                    for &(start, end) in &level.pus {
                        if pe.contains_key(&(start, end)) {
                            continue;
                        }
                        let pu_res_list = energy::build_residue_list((start, end), &chain_id.to_string());
                        if let Ok(result) = energy::get_energy_and_z_score(
                            &energy_config,
                            &pdb_path_str,
                            Some(&pu_res_list),
                        ) {
                            pe.insert((start, end), (result.energy, result.z_score));
                        }
                    }
                }
                Some(pe)
            } else {
                None
            };

            peeling::write_peeling_summary(
                &peeling_levels,
                &results_dir.join("PEELING_summary.txt"),
                peeling_energies.as_ref(),
            )?;
            tracing::info!("Wrote peeling summary");
        }
    }

    // Step 10: Generate plots (skip for now - user said don't worry about images)
    if !cli.disable_plots {
        let contact_matrix_dir = results_dir.join("Contact_Probability_Matrix");
        std::fs::create_dir_all(&contact_matrix_dir)?;

        let proba_mat_file = results_dir
            .join("PDBs_Clean")
            .join(&pdb_id_chain)
            .join("file_proba_contact.mat");

        if proba_mat_file.exists() {
            tracing::info!("Skipping contact probability matrix plots (use Python version for plots)");
        }
    }

    // Step 11: Calculate junction consistencies
    let stat_script = bin_dir.join("stat_pu_domains_from_SWORD.pl");
    if stat_script.exists() {
        tracing::info!("Calculate junctions consistencies");
        let junctions_output = std::process::Command::new(&stat_script)
            .arg(results_dir.join("sword.txt"))
            .output()
            .context("Failed to run junctions script")?;

        if junctions_output.status.success() {
            let stdout = String::from_utf8_lossy(&junctions_output.stdout);
            // Python writes each line + "\n", including trailing empty lines
            let mut content = String::new();
            for line in stdout.lines() {
                content.push_str(line);
                content.push('\n');
            }
            content.push('\n'); // Match Python's trailing newline
            std::fs::write(
                results_dir.join("junctions_consistencies.txt"),
                &content,
            )?;
        }
    }

    // Step 12: Write mapping file
    pdb::writer::write_mapping_file(
        &original_resnums,
        &results_dir.join("mapping_auth_resnums.txt"),
    )?;

    // Step 13: Clean and prepare results (same as Python)
    tracing::info!("Clean and prepare results");
    let pdbs_stand = results_dir.join("PDBs_Stand");
    if pdbs_stand.exists() {
        let _ = std::fs::remove_dir_all(&pdbs_stand);
    }

    let pdbs_clean = results_dir.join("PDBs_Clean");
    let sword_dir_dest = results_dir.join("SWORD");
    if pdbs_clean.exists() {
        let _ = std::fs::rename(&pdbs_clean, &sword_dir_dest);
    }

    // Move junctions file into Junctions/ directory
    let junctions_dir = results_dir.join("Junctions");
    let junctions_file = results_dir.join("junctions_consistencies.txt");
    if junctions_file.exists() {
        std::fs::create_dir_all(&junctions_dir)?;
        let _ = std::fs::rename(&junctions_file, junctions_dir.join("junctions_consistencies.txt"));
    }

    // Move Peeling directory
    let peeling_dir_glob = results_dir.join("SWORD").join(&pdb_id_chain).join("Peeling");
    if peeling_dir_glob.exists() {
        let _ = std::fs::rename(&peeling_dir_glob, results_dir.join("Protein_Units"));
    }

    let elapsed = start.elapsed();
    tracing::info!("Results can be found here: {}", results_dir.display());
    tracing::info!("Total runtime: {} seconds", elapsed.as_secs());

    Ok(())
}

/// Resolve the input source to a file path and base name.
fn resolve_input(cli: &Cli, output_dir: &PathBuf) -> Result<(PathBuf, String)> {
    if let Some(ref input_file) = cli.input_file {
        let base = input_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        Ok((input_file.clone(), base))
    } else if let Some(ref uniprot_id) = cli.uniprot_id {
        let path = fetch::fetch_alphafold(uniprot_id, output_dir)?;
        Ok((path, uniprot_id.clone()))
    } else if let Some(ref mgnify_id) = cli.mgnify_id {
        let path = fetch::fetch_esm(mgnify_id, output_dir)?;
        Ok((path, mgnify_id.clone()))
    } else if let Some(ref pdb_id) = cli.pdb_id {
        let path = fetch::fetch_pdb(pdb_id, output_dir)?;
        Ok((path, pdb_id.to_uppercase()))
    } else {
        anyhow::bail!("No input source specified");
    }
}
