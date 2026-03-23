//! SWORD2 CLI: SWift and Optimized Recognition of protein Domains.
//!
//! Command-line interface for running the SWORD2 protein domain recognition pipeline.
//! Fully in Rust — no Perl dependencies.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};
use clap::{ArgAction, Parser};
use console::Style;
use indicatif::{ProgressBar, ProgressStyle};

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

    /// Path to SWORD2 base directory (defaults to parent of binary location)
    #[arg(long)]
    base_dir: Option<PathBuf>,

    /// Number of random shuffles for Z-score calculation (default: 2000, lower = faster)
    #[arg(short = 's', long, default_value = "2000")]
    num_shuffles: usize,

    /// Increase verbosity (-v steps, -vv debug, -vvv trace)
    #[arg(short = 'v', long = "verbose", action = ArgAction::Count)]
    verbosity: u8,

    /// Suppress all output except errors
    #[arg(short = 'q', long)]
    quiet: bool,
}

// ---------------------------------------------------------------------------
// Reporter — manages all user-visible terminal output
// ---------------------------------------------------------------------------

struct Reporter {
    verbosity: u8,
    quiet: bool,
    spinner: Option<ProgressBar>,
    step_start: Option<Instant>,
    // Styles
    s_ok: Style,
    s_header: Style,
    s_step: Style,
    s_dim: Style,
    s_warn: Style,
    s_err: Style,
}

impl Reporter {
    fn new(verbosity: u8, quiet: bool) -> Self {
        let spinner = if verbosity == 0 && !quiet {
            let pb = ProgressBar::new_spinner();
            pb.set_style(
                ProgressStyle::with_template(" {spinner:.cyan}  {msg}")
                    .unwrap()
                    .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏", " "]),
            );
            pb.enable_steady_tick(Duration::from_millis(80));
            Some(pb)
        } else {
            None
        };

        Self {
            verbosity,
            quiet,
            spinner,
            step_start: None,
            s_ok: Style::new().green().bold(),
            s_header: Style::new().bold(),
            s_step: Style::new().cyan(),
            s_dim: Style::new().dim(),
            s_warn: Style::new().yellow().bold(),
            s_err: Style::new().red().bold(),
        }
    }

    /// Print the header block (only in -v mode).
    fn begin(&self, id: &str, chain: char, residues: usize, threads: usize) {
        if self.quiet {
            return;
        }
        if self.verbosity >= 1 {
            eprintln!(
                "\n {}  {}",
                self.s_header.apply_to("▸"),
                self.s_header.apply_to(format!(
                    "{}  ·  chain {}  ·  {} residues  ·  {} threads",
                    id, chain, residues, threads
                )),
            );
            eprintln!();
        } else if let Some(ref pb) = self.spinner {
            pb.set_message(format!("{} · Initializing…", id));
        }
    }

    /// Begin a step — show spinner text (v=0) or a running line (v≥1).
    fn step(&mut self, label: &str) {
        if self.quiet {
            return;
        }
        self.step_start = Some(Instant::now());
        if self.verbosity >= 1 {
            // Don't print running indicator — step_done will print the completed line
        } else if let Some(ref pb) = self.spinner {
            // Update the spinner message with the current protein ID prefix
            let current = pb.message();
            // Keep prefix (before first ·) and replace suffix
            let prefix = current.split('·').next().unwrap_or("").trim();
            pb.set_message(format!("{} · {}…", prefix, label));
        }
    }

    /// Complete the current step — print ✓ line (v≥1) or do nothing (v=0, spinner updates).
    fn step_done(&mut self, label: &str, detail: Option<&str>) {
        if self.quiet {
            return;
        }
        if self.verbosity >= 1 {
            let elapsed = self.step_start.map(|s| s.elapsed());
            let mut line = format!(
                " {}  {:<40}",
                self.s_ok.apply_to("✓"),
                label,
            );
            if let Some(d) = detail {
                write!(line, " {}", self.s_dim.apply_to(d)).ok();
            }
            if self.verbosity >= 2 {
                if let Some(el) = elapsed {
                    write!(line, "  {}", self.s_dim.apply_to(format_duration(el))).ok();
                }
            }
            eprintln!("{}", line.trim_end());
        }
        self.step_start = None;
    }

    /// Print a warning.
    #[allow(dead_code)]
    fn warn(&self, msg: &str) {
        if self.quiet {
            return;
        }
        if self.verbosity >= 1 {
            eprintln!(" {}  {}", self.s_warn.apply_to("⚠"), msg);
        }
        // In v=0, warnings are suppressed (tracing handles them for -vv+)
    }

    /// Print the final summary line.
    fn finish(&self, id: &str, n_domains: usize, residues: usize, elapsed: Duration, path: &std::path::Path) {
        if self.quiet {
            return;
        }
        if let Some(ref pb) = self.spinner {
            pb.finish_and_clear();
        }

        let dur = format_duration(elapsed);
        let rel_path = path
            .strip_prefix(std::env::current_dir().unwrap_or_default())
            .unwrap_or(path);

        if self.verbosity >= 1 {
            eprintln!();
            eprintln!(
                " {}  {}  ·  {} domain{}  ·  {}  →  {}",
                self.s_ok.apply_to("✓"),
                self.s_header.apply_to("Done"),
                n_domains,
                if n_domains == 1 { "" } else { "s" },
                self.s_dim.apply_to(&dur),
                self.s_step.apply_to(rel_path.display()),
            );
            eprintln!();
        } else {
            eprintln!(
                " {}  {}  ·  {} residues  ·  {} domain{}  ·  {}  →  {}",
                self.s_ok.apply_to("✓"),
                self.s_header.apply_to(id),
                residues,
                n_domains,
                if n_domains == 1 { "" } else { "s" },
                self.s_dim.apply_to(&dur),
                self.s_step.apply_to(rel_path.display()),
            );
        }
    }

    /// Print a fatal error before exit.
    #[allow(dead_code)]
    fn error(&self, msg: &str) {
        eprintln!(" {}  {}", self.s_err.apply_to("✗"), msg);
    }
}

/// Format a Duration into a human-readable string.
fn format_duration(d: Duration) -> String {
    let ms = d.as_millis();
    if ms < 1000 {
        format!("{}ms", ms)
    } else {
        format!("{:.2}s", d.as_secs_f64())
    }
}

/// Configure tracing-subscriber based on verbosity.
fn setup_logging(verbosity: u8, quiet: bool) {
    use tracing_subscriber::fmt::format::FmtSpan;

    let level = if quiet {
        "error"
    } else {
        match verbosity {
            0 | 1 => "warn",
            2 => "debug",
            _ => "trace",
        }
    };

    let directive = format!("sword2={}", level);
    let env_filter = tracing_subscriber::EnvFilter::from_default_env()
        .add_directive(directive.parse().unwrap());

    if verbosity >= 3 {
        // Full trace: timestamps, spans, level, target
        tracing_subscriber::fmt()
            .with_env_filter(env_filter)
            .with_span_events(FmtSpan::CLOSE)
            .init();
    } else {
        // Compact: no timestamps, colored level prefix
        tracing_subscriber::fmt()
            .with_env_filter(env_filter)
            .without_time()
            .with_target(verbosity >= 2)
            .init();
    }
}

fn main() -> Result<()> {
    let start = Instant::now();

    let cli = Cli::parse();

    // Initialize logging based on verbosity
    setup_logging(cli.verbosity, cli.quiet);

    let mut reporter = Reporter::new(cli.verbosity, cli.quiet);

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

    // Ensure output directory exists
    let output_dir = std::fs::canonicalize(&cli.output_dir).unwrap_or_else(|_| {
        std::fs::create_dir_all(&cli.output_dir).ok();
        cli.output_dir.clone()
    });
    std::fs::create_dir_all(&output_dir)?;

    // Step 1: Obtain the structure file
    reporter.step("Fetch structure");
    let (input_path, pdb_id_base, is_fetched) = resolve_input(&cli, &output_dir)?;
    reporter.step_done("Fetch structure", None);

    // Step 2: Parse the structure
    reporter.step("Parse structure");
    tracing::debug!("Parsing structure from {}", input_path.display());
    let structure = pdb::parse_pdb(&input_path)
        .with_context(|| format!("Failed to parse {}", input_path.display()))?;

    // Determine chain
    let chain_id = if let Some(ref chain_str) = cli.chain {
        chain_str.chars().next().unwrap_or('A')
    } else if let Some(model) = structure.first_model() {
        // Find the first chain that contains at least one standard amino acid residue
        if let Some(chain) = model.chains.iter().find(|c| {
            c.residues.iter().any(|r| pdb::amino_acids::is_standard(&r.name))
        }) {
            tracing::debug!("No chain specified. Using first protein chain '{}'", chain.id);
            chain.id
        } else if let Some(chain) = model.chains.first() {
            tracing::debug!("No chain specified. Using first chain '{}'", chain.id);
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

    // Delete downloaded PDB from output root (it was only needed for parsing)
    if is_fetched {
        let _ = std::fs::remove_file(&input_path);
    }

    tracing::debug!(">>> {} ({} residues)", pdb_id_chain, chain.len());
    tracing::debug!(">>> Using {} cpus", num_threads);

    // Step 3: Clean PDB - remove non-standard residues, insertion codes, renumber from 1
    tracing::debug!("Write a clean version of the PDB: remove non standard residues");
    let (cleaned_chain, original_resnums) = pdb::writer::clean_chain_for_sword(chain);

    if cleaned_chain.is_empty() {
        anyhow::bail!(
            "No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file."
        );
    }

    let prot_len = cleaned_chain.len();
    tracing::debug!("Clean chain: {} residues, sequence: {}", prot_len, cleaned_chain.get_sequence());
    reporter.step_done("Parse & clean PDB", Some(&format!("{} residues", prot_len)));

    // Print header now that we know all details
    reporter.begin(&pdb_id_chain, chain_id, prot_len, num_threads);

    // Write clean PDB as input.pdb
    let input_pdb = results_dir.join("input.pdb");
    pdb::write_pdb(&cleaned_chain, &input_pdb)?;

    // Write .num file with sequential 1-based numbering matching the clean PDB.
    // Must exist before run_pipeline.
    {
        let intermediate_dir = results_dir.join("intermediate");
        std::fs::create_dir_all(&intermediate_dir)?;
        let num_file = intermediate_dir.join(format!("{}.num", pdb_id_chain));
        let n = original_resnums.len();
        let num_content: Vec<String> = (1..=n).map(|i| i.to_string()).collect();
        std::fs::write(&num_file, num_content.join(" "))?;
    }

    // Step 4: Run the SWORD pipeline (DSSP is pure Rust, no compilation needed)
    reporter.step("SWORD pipeline");
    tracing::debug!("Launch SWORD pipeline");
    let config = sword::SwordConfig {
        compute_energies: !cli.disable_energies,
        generate_plots: !cli.disable_plots,
        num_threads,
        output_dir: results_dir.to_string_lossy().to_string(),
        max_alternatives: 9,
    };

    let (sword_output, sword_results) = sword::run_pipeline(
        &input_pdb,
        &pdb_id_chain,
        &config,
    ).context("Failed to run SWORD pipeline")?;

    // Step 6: Parse SWORD output (already done in run_pipeline)
    let n_domains = sword_results.domains.first().map_or(0, |p| p.nb_domains);
    tracing::debug!(
        "Found {} partitioning(s), ambiguity: {}",
        sword_results.domains.len(),
        sword_results.ambiguity
    );
    reporter.step_done(
        "SWORD pipeline",
        Some(&format!(
            "{} domain{}, {} partitioning{}",
            n_domains,
            if n_domains == 1 { "" } else { "s" },
            sword_results.domains.len(),
            if sword_results.domains.len() == 1 { "" } else { "s" },
        )),
    );

    // Step 7: Calculate energies
    // First, collect all unique PU ranges across SWORD partitions (for batch dedup)
    let energy_config = if !cli.disable_energies {
        let mut ec = energy::EnergyConfig::from_bin_dir(&bin_dir.to_string_lossy());
        ec.num_shuffles = cli.num_shuffles;
        Some(ec)
    } else {
        None
    };
    let pdb_path_str = input_pdb.to_string_lossy().to_string();
    let chain_str = chain_id.to_string();

    // Pre-compute all unique PU energies in one parallel batch
    let pu_energy_cache: HashMap<(i32, i32), energy::EnergyResult> =
        if let Some(ref ec) = energy_config {
            reporter.step("Pseudo-energies (PUs)");
            let mut all_pu_ranges: std::collections::HashSet<(i32, i32)> = std::collections::HashSet::new();
            for part in &sword_results.domains {
                for domain in &part.boundaries {
                    for &(s, e) in domain {
                        all_pu_ranges.insert((s, e));
                    }
                }
            }
            let unique_ranges: Vec<(i32, i32)> = all_pu_ranges.into_iter().collect();
            let cache = energy::compute_pu_energies_batch(ec, &pdb_path_str, &chain_str, &unique_ranges);
            reporter.step_done(
                "Pseudo-energies (PUs)",
                Some(&format!("{} unique ranges", unique_ranges.len())),
            );
            cache
        } else {
            HashMap::new()
        };

    // Now compute domain-level energies using cached PU results
    let energies: HashMap<energy::EnergyKey, energy::EnergyResult> =
        if let Some(ref ec) = energy_config {
            reporter.step("Pseudo-energies (domains)");
            let e = energy::calculate_all_energies_with_cache(
                ec,
                &pdb_path_str,
                &chain_str,
                &sword_results.domains,
                &pu_energy_cache,
            );
            reporter.step_done("Pseudo-energies (domains)", None);
            e
        } else {
            HashMap::new()
        };

    // Step 8: Write SWORD partitionings (text + JSON)
    reporter.step("Write results");
    tracing::debug!("Write the SWORD results");
    output::write_sword_summary(
        &sword_results,
        &energies,
        cli.disable_energies,
        &results_dir.join("summary.txt"),
    )?;
    output::write_sword_summary_json(
        &sword_results,
        &energies,
        cli.disable_energies,
        &results_dir.join("summary.json"),
    )?;

    // Step 9: Write Peeling results
    tracing::debug!("Write Peeling results");
    let peeling_num = results_dir
        .join("intermediate")
        .join(format!("{}.num", pdb_id_chain));
    let peeling_log = results_dir
        .join("intermediate")
        .join("peeling.log");

    if peeling_log.exists() {
        if peeling_num.exists() {
            let peeling_results = peeling::load_legacy_results(&peeling_log, &peeling_num)?;

            // Calculate peeling energies if enabled — reuse global PU cache + compute missing
            let peeling_energies = if let Some(ref ec) = energy_config {
                // Collect unique PU ranges not already in the cache
                let mut missing_ranges: Vec<(i32, i32)> = Vec::new();
                let mut seen_pus: std::collections::HashSet<(i32, i32)> = std::collections::HashSet::new();
                for level in &peeling_results.levels {
                    for range in &level.pus {
                        if seen_pus.insert((range.start, range.end))
                            && !pu_energy_cache.contains_key(&(range.start, range.end))
                        {
                            missing_ranges.push((range.start, range.end));
                        }
                    }
                }
                // Compute only the missing ones in parallel
                let extra = if !missing_ranges.is_empty() {
                    energy::compute_pu_energies_batch(ec, &pdb_path_str, &chain_str, &missing_ranges)
                } else {
                    HashMap::new()
                };
                // Build peeling energy map from both caches
                let mut pe: HashMap<(i32, i32), (Option<f64>, Option<f64>)> = HashMap::new();
                for level in &peeling_results.levels {
                    for range in &level.pus {
                        if pe.contains_key(&(range.start, range.end)) {
                            continue;
                        }
                        if let Some(r) = pu_energy_cache
                            .get(&(range.start, range.end))
                            .or_else(|| extra.get(&(range.start, range.end)))
                        {
                            pe.insert((range.start, range.end), (r.energy, r.z_score));
                        }
                    }
                }
                Some(pe)
            } else {
                None
            };

            peeling::write_peeling_summary(
                &peeling_results,
                &results_dir.join("peeling.txt"),
                peeling_energies.as_ref(),
            )?;
            peeling::write_peeling_summary_json(
                &peeling_results,
                &results_dir.join("peeling.json"),
                peeling_energies.as_ref(),
            )?;
            tracing::debug!("Wrote peeling summary");
        }
    }
    reporter.step_done("Write results", None);

    // Step 10: Generate plots
    if config.generate_plots {
        reporter.step("Generate plots");
        let plots_dir = results_dir.join("plots");
        std::fs::create_dir_all(&plots_dir)?;

        let proba_mat_file = results_dir
            .join("intermediate")
            .join("contact_matrix.mat");

        // Domain consistency histogram (SVG)
        let histogram_output = plots_dir.join("domain_histogram.svg");
        let domain_counts = plot::count_domains(&sword_results.domains);
        if let Err(err) = plot::write_domain_histogram(
            &domain_counts,
            &histogram_output.to_string_lossy(),
        ) {
            tracing::warn!(error = %err, "Failed to write domain consistency histogram");
            reporter.warn("Could not generate the domain consistency histogram");
        }

        // 3-level contact probability matrix plots (PNG)
        if !proba_mat_file.exists() {
            tracing::warn!(
                matrix = %proba_mat_file.display(),
                "Skipping contact probability matrix plots because the matrix file is missing"
            );
            reporter.warn("Skipping contact matrix plots because the matrix file is missing");
        } else {
            match plot::load_contact_matrix(&proba_mat_file) {
                Ok(matrix) => {
                    let pu_colors = plot::assign_pu_colors(&sword_results.domains);

                    for (i, partition) in sword_results.domains.iter().enumerate() {
                        if let Err(err) = plot::generate_alternative_plots(
                            &matrix,
                            i,
                            partition,
                            &pu_colors,
                            &plots_dir,
                        ) {
                            tracing::warn!(
                                error = %err,
                                alt = i,
                                "Failed to write contact probability matrix plots for alternative {}",
                                i
                            );
                        }
                    }
                }
                Err(err) => {
                    tracing::warn!(error = %err, "Failed to load contact probability matrix input");
                    reporter.warn("Skipping contact matrix plots because the matrix input could not be loaded");
                }
            }
        }

        reporter.step_done("Generate plots", Some("PNG outputs"));
    }

    // Step 11: Calculate junction consistencies (pure Rust — no Perl!)
    reporter.step("Junctions & cleanup");
    tracing::debug!("Calculate junctions consistencies");
    let junctions_content =
        sword::junctions::calculate_junction_consistencies(&sword_output);
    if !junctions_content.is_empty() {
        std::fs::write(
            results_dir.join("junctions.txt"),
            &junctions_content,
        )?;
    }

    // Step 12: Write mapping file
    pdb::writer::write_mapping_file(
        &original_resnums,
        &results_dir.join("residue_mapping.txt"),
    )?;

    // Step 13: Clean up legacy artifacts
    tracing::debug!("Clean up results");
    let pdbs_stand = results_dir.join("PDBs_Stand");
    if pdbs_stand.exists() {
        let _ = std::fs::remove_dir_all(&pdbs_stand);
    }
    reporter.step_done("Junctions & cleanup", None);

    let elapsed = start.elapsed();
    reporter.finish(&pdb_id_chain, n_domains, prot_len, elapsed, &results_dir);

    Ok(())
}

/// Resolve the input source to a file path, base name, and whether it was fetched.
fn resolve_input(cli: &Cli, output_dir: &PathBuf) -> Result<(PathBuf, String, bool)> {
    if let Some(ref input_file) = cli.input_file {
        let base = input_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        Ok((input_file.clone(), base, false))
    } else if let Some(ref uniprot_id) = cli.uniprot_id {
        let path = fetch::fetch_alphafold(uniprot_id, output_dir)?;
        Ok((path, uniprot_id.clone(), true))
    } else if let Some(ref mgnify_id) = cli.mgnify_id {
        let path = fetch::fetch_esm(mgnify_id, output_dir)?;
        Ok((path, mgnify_id.clone(), true))
    } else if let Some(ref pdb_id) = cli.pdb_id {
        let path = fetch::fetch_pdb(pdb_id, output_dir)?;
        Ok((path, pdb_id.to_uppercase(), true))
    } else {
        anyhow::bail!("No input source specified");
    }
}
