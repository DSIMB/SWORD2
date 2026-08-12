//! SWORD2 CLI: SWift and Optimized Recognition of protein Domains.
//!
//! Command-line interface for running the SWORD2 protein domain recognition pipeline.
//! Fully in Rust — no Perl dependencies.

use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};
use clap::{ArgAction, Args, Parser, Subcommand};
use console::Style;
use indicatif::{ProgressBar, ProgressStyle};

use sword2_lib::{energy, fetch, output, pdb, peeling, plot, sword};

/// SWORD2: SWift and Optimized Recognition of protein Domains
#[derive(Parser, Debug)]
#[command(name = "sword2", version, about)]
struct Cli {
    /// Optional subcommand. When omitted, the full domain-partitioning pipeline runs.
    #[command(subcommand)]
    command: Option<Commands>,

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
    nmr_model: i32,

    /// Enable pseudo-energy calculations
    #[arg(short = 'E', long)]
    energies: bool,

    /// Use the pairwise-trained reranker (energy Z-score + boundary secondary
    /// structure) to pick the winning partition instead of the legacy
    /// distance_model selection. Off by default; benchmark before enabling.
    #[arg(long)]
    use_pairwise_reranker: bool,

    /// Select domain count and partition with the embedded factorized structural
    /// ranker. Experimental and off by default; incomplete evidence falls back
    /// to the legacy selector.
    #[arg(long, conflicts_with_all = [
        "use_pairwise_reranker",
        "use_count_calibration",
        "count_lambda",
        "use_geometry_metrics",
        "geometry_lambda",
    ])]
    use_factorized_ranker: bool,

    /// Bias domain-count selection toward a length-predicted count (experimental)
    #[arg(long)]
    use_count_calibration: bool,

    /// Penalty weight for count calibration (default: fitted constant)
    #[arg(long)]
    count_lambda: Option<f64>,

    /// Reorder the winning candidate using analytical "ideal sphere" geometry
    /// criteria (sphericity, density, inter-domain interface fraction).
    /// These criteria are always computed and shown in the output regardless
    /// of this flag; it only changes which candidate gets picked.
    #[arg(long)]
    use_geometry_metrics: bool,

    /// Penalty weight for geometry-based reordering (default: provisional constant)
    #[arg(long)]
    geometry_lambda: Option<f64>,

    /// Enable generation of contact probability matrix plots
    #[arg(short = 'P', long)]
    plots: bool,

    /// Number of threads for parallel computation (0 = all CPUs)
    #[arg(short = 'j', long, default_value = "0")]
    threads: usize,

    /// Path to SWORD2 installation directory (defaults to parent of binary location)
    #[arg(long)]
    install_dir: Option<PathBuf>,

    /// Number of random shuffles for Z-score calculation (default: 2000, lower = faster)
    #[arg(short = 'z', long, default_value = "2000")]
    zscore_shuffles: usize,

    /// Skip if output already exists (summary.json present in output directory)
    #[arg(long)]
    skip_existing: bool,

    /// Write one PDB file per domain of the optimal partition into <output>/domains_optimal/
    #[arg(long)]
    extract_domains: bool,

    /// Drop residues with pLDDT (AlphaFold/ESM confidence score) below this value [0–100].
    /// Use e.g. --min-plddt 70 to discard low-confidence disordered regions before analysis.
    /// Has no effect on experimental PDB structures (their B-factors mean something different).
    #[arg(long, value_name = "0-100")]
    min_plddt: Option<f64>,

    /// Output format for stdout summary: text (default), tsv, json
    #[arg(long, default_value = "text")]
    format: output::OutputFormat,

    /// Fetch legacy PDB format (.pdb) from RCSB instead of mmCIF (.cif).
    /// Use this only if you specifically need PDB format; some entries are
    /// unavailable in PDB format and will error.
    #[arg(long)]
    legacy_pdb: bool,

    /// Batch file: one structure per line (PDB ID, af:UNIPROT, esm:MGNIFY, or file path)
    #[arg(long)]
    batch: Option<PathBuf>,

    /// Increase verbosity (-v steps, -vv debug, -vvv trace)
    #[arg(short = 'v', long = "verbose", action = ArgAction::Count)]
    verbosity: u8,

    /// Suppress all output except errors
    #[arg(short = 'q', long)]
    quiet: bool,
}

/// Subcommands (the default no-subcommand form runs the full pipeline).
#[derive(Subcommand, Debug)]
enum Commands {
    /// Score one structure's pseudo-energy (and Z-score) for a whole structure
    /// or residue subset. Prints `Pseudo-energy = …` and `Z-score = …`, one per
    /// line. Useful for scoring a single structure or debugging.
    Score(ScoreArgs),
}

/// Arguments for the `score` subcommand.
#[derive(Args, Debug)]
struct ScoreArgs {
    /// Path to the (cleaned) structure file to score.
    #[arg(long)]
    pdb: PathBuf,

    /// Path to the potential directory (e.g. bin/mypmfs-master/025_30_100_potential).
    #[arg(long)]
    potential_dir: String,

    /// Skip the Z-score computation (energy only).
    #[arg(long)]
    no_zscore: bool,

    /// Number of random shuffles for the Z-score.
    #[arg(long, default_value = "2000")]
    shuffles: usize,

    /// Threads for parallel decoy scoring (0 = all CPUs).
    #[arg(long, default_value = "0")]
    cpu: usize,

    /// Restrict scoring to a residue subset (comma-separated `numchain` tokens, e.g. "1A,2A").
    #[arg(long)]
    residues: Option<String>,
}

/// Run the `score` subcommand: load potentials, score the structure, print results.
fn run_score(args: &ScoreArgs) -> Result<()> {
    if args.cpu > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.cpu)
            .build_global()
            .ok();
    }
    let pot = energy::Potentials::load(&args.potential_dir)?;
    let pdb_str = args.pdb.to_string_lossy();
    let result = energy::score_structure(
        &pot,
        &pdb_str,
        args.residues.as_deref(),
        args.shuffles,
        !args.no_zscore,
    )?;
    if let Some(e) = result.energy {
        println!("Pseudo-energy = {e}");
    }
    if !args.no_zscore {
        if let Some(z) = result.z_score {
            println!("Z-score = {z}");
        }
    }
    Ok(())
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
            let mut line = format!(" {}  {:<40}", self.s_ok.apply_to("✓"), label,);
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

    /// Print a plain informational note (always visible unless quiet).
    fn note(&self, msg: &str) {
        if self.quiet {
            return;
        }
        if let Some(ref pb) = self.spinner {
            pb.println(format!(" {}  {}", self.s_ok.apply_to("·"), msg));
        } else {
            eprintln!(" {}  {}", self.s_ok.apply_to("·"), msg);
        }
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
    fn finish(
        &self,
        id: &str,
        n_domains: usize,
        residues: usize,
        elapsed: Duration,
        path: &std::path::Path,
    ) {
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
    let env_filter =
        tracing_subscriber::EnvFilter::from_default_env().add_directive(directive.parse().unwrap());

    if verbosity >= 3 {
        // Full trace: timestamps, spans, level, target
        tracing_subscriber::fmt()
            .with_env_filter(env_filter)
            .with_writer(std::io::stderr)
            .with_span_events(FmtSpan::CLOSE)
            .init();
    } else {
        // Compact: no timestamps, colored level prefix
        tracing_subscriber::fmt()
            .with_env_filter(env_filter)
            .with_writer(std::io::stderr)
            .without_time()
            .with_target(verbosity >= 2)
            .init();
    }
}

// ---------------------------------------------------------------------------
// Batch support — entry descriptor and line parser
// ---------------------------------------------------------------------------

/// One structure to process (populated from CLI flags or a batch file line).
struct EntryArgs {
    pdb_id: Option<String>,
    uniprot_id: Option<String>,
    mgnify_id: Option<String>,
    input_file: Option<PathBuf>,
    /// Explicit chain override (overrides auto-detection).
    chain: Option<char>,
}

impl EntryArgs {
    fn label(&self) -> &str {
        self.pdb_id
            .as_deref()
            .or(self.uniprot_id.as_deref())
            .or(self.mgnify_id.as_deref())
            .or_else(|| self.input_file.as_ref().and_then(|p| p.to_str()))
            .unwrap_or("?")
    }

    fn is_predicted(&self) -> bool {
        self.uniprot_id.is_some() || self.mgnify_id.is_some()
    }
}

/// Parse one batch-file line into an `EntryArgs`.
///
/// Prefixes: `af:` → AlphaFold, `esm:` → ESM Atlas.
/// Paths: any value containing `/` or starting with `.`.
/// Optional `:A` suffix on any value selects a chain.
/// Lines starting with `#` or empty are skipped (return `None`).
fn parse_batch_line(line: &str) -> Option<EntryArgs> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    if let Some(rest) = line.strip_prefix("af:") {
        let (value, chain) = split_chain_suffix(rest);
        return Some(EntryArgs {
            pdb_id: None,
            uniprot_id: Some(value.to_string()),
            mgnify_id: None,
            input_file: None,
            chain,
        });
    }
    if let Some(rest) = line.strip_prefix("esm:") {
        let (value, chain) = split_chain_suffix(rest);
        return Some(EntryArgs {
            pdb_id: None,
            uniprot_id: None,
            mgnify_id: Some(value.to_string()),
            input_file: None,
            chain,
        });
    }
    if line.contains('/') || line.starts_with('.') {
        let (value, chain) = split_chain_suffix(line);
        return Some(EntryArgs {
            pdb_id: None,
            uniprot_id: None,
            mgnify_id: None,
            input_file: Some(PathBuf::from(value)),
            chain,
        });
    }
    let (value, chain) = split_chain_suffix(line);
    Some(EntryArgs {
        pdb_id: Some(value.to_string()),
        uniprot_id: None,
        mgnify_id: None,
        input_file: None,
        chain,
    })
}

/// Split `"VALUE:A"` → `("VALUE", Some('A'))`, or `("VALUE", None)` if no single-char suffix.
fn split_chain_suffix(s: &str) -> (&str, Option<char>) {
    if let Some(pos) = s.rfind(':') {
        let suffix = &s[pos + 1..];
        if suffix.len() == 1 {
            if let Some(c) = suffix.chars().next().filter(|c| c.is_ascii_alphabetic()) {
                return (&s[..pos], Some(c));
            }
        }
    }
    (s, None)
}

// ---------------------------------------------------------------------------
// Per-entry pipeline
// ---------------------------------------------------------------------------

/// Run the full SWORD2 pipeline for a single structure entry.
fn process_entry(
    entry: &EntryArgs,
    cli: &Cli,
    reporter: &mut Reporter,
    bin_dir: &std::path::Path,
    output_dir: &PathBuf,
    num_threads: usize,
) -> Result<()> {
    let start = Instant::now();

    // Step 1: Obtain the structure file
    reporter.step("Fetch structure");
    let (input_path, pdb_id_base, is_fetched) = resolve_input(entry, output_dir, cli.legacy_pdb)?;
    reporter.step_done("Fetch structure", None);

    // Step 2: Parse the structure
    reporter.step("Parse structure");
    tracing::debug!("Parsing structure from {}", input_path.display());
    let structure = pdb::parse_pdb(&input_path)
        .with_context(|| format!("Failed to parse {}", input_path.display()))?;

    // Determine chain
    let chain_id = if let Some(c) = entry.chain {
        c
    } else if let Some(model) = structure.first_model() {
        // Find the first chain that contains at least one standard amino acid residue
        if let Some(chain) = model.chains.iter().find(|c| {
            c.residues
                .iter()
                .any(|r| pdb::amino_acids::is_standard(&r.name))
        }) {
            tracing::debug!(
                "No chain specified. Using first protein chain '{}'",
                chain.id
            );
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
    let chain = model.get_chain(chain_id).ok_or_else(|| {
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

    if cli.skip_existing && results_dir.join("summary.json").exists() {
        if let Some(ref pb) = reporter.spinner {
            pb.finish_and_clear();
        }
        if !cli.quiet {
            eprintln!(" ⏭  {}  · already exists, skipping", pdb_id_chain);
        }
        return Ok(());
    }

    // Delete downloaded PDB from output root (it was only needed for parsing)
    if is_fetched {
        let _ = std::fs::remove_file(&input_path);
    }

    tracing::debug!(">>> {} ({} residues)", pdb_id_chain, chain.len());
    tracing::debug!(">>> Using {} cpus", num_threads);

    // Step 3: Clean PDB - remove non-standard residues, insertion codes, renumber from 1
    tracing::debug!("Write a clean version of the PDB: remove non standard residues");
    // Optional pLDDT filter (AlphaFold/ESM structures store confidence in B-factor column)
    let plddt_filtered;
    let chain = if let Some(min_plddt) = cli.min_plddt {
        if !entry.is_predicted() {
            tracing::warn!(
                "--min-plddt is intended for AlphaFold/ESM structures; \
                 B-factors in experimental PDB structures have different semantics"
            );
        }
        plddt_filtered = pdb::writer::filter_by_plddt(chain, min_plddt);
        &plddt_filtered
    } else {
        chain
    };
    let (cleaned_chain, original_resnums) = pdb::writer::clean_chain_for_sword(chain);

    if cleaned_chain.is_empty() {
        anyhow::bail!(
            "No atomic data is left after trying to keep the 20 classical residues. Please check your PDB file."
        );
    }

    let prot_len = cleaned_chain.len();
    tracing::debug!(
        "Clean chain: {} residues, sequence: {}",
        prot_len,
        cleaned_chain.get_sequence()
    );
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

    // Reduced-shuffle energy config for candidate rescoring/training dump —
    // separate from the `-E` display config below (which uses cli.zscore_shuffles,
    // default 2000). Construction is free (potentials load lazily on first use),
    // so this is built unconditionally without affecting default runtime.
    let mut rerank_energy_config = energy::EnergyConfig::from_bin_dir(&bin_dir.to_string_lossy());
    rerank_energy_config.num_shuffles = 200;

    let config = sword::SwordConfig {
        compute_energies: cli.energies,
        generate_plots: cli.plots,
        num_threads,
        output_dir: results_dir.to_string_lossy().to_string(),
        max_alternatives: 9,
        energy_config: Some(rerank_energy_config),
        chain_id: chain_id.to_string(),
        use_pairwise_reranker: cli.use_pairwise_reranker,
        use_factorized_ranker: cli.use_factorized_ranker,
        use_count_calibration: cli.use_count_calibration,
        count_lambda: cli.count_lambda,
        use_geometry_metrics: cli.use_geometry_metrics,
        geometry_lambda: cli.geometry_lambda,
    };

    let (sword_output, sword_results) = sword::run_pipeline(&input_pdb, &pdb_id_chain, &config)
        .context("Failed to run SWORD pipeline")?;

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
            if sword_results.domains.len() == 1 {
                ""
            } else {
                "s"
            },
        )),
    );

    // Step 7: Calculate energies
    // First, collect all unique PU ranges across SWORD partitions (for batch dedup)
    let energy_config = if cli.energies {
        let mut ec = energy::EnergyConfig::from_bin_dir(&bin_dir.to_string_lossy());
        ec.num_shuffles = cli.zscore_shuffles;
        // Load potentials once up front, before the parallel scoring batch.
        ec.preload()?;
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
            let mut all_pu_ranges: std::collections::HashSet<(i32, i32)> =
                std::collections::HashSet::new();
            for part in &sword_results.domains {
                for domain in &part.boundaries {
                    for &(s, e) in domain {
                        all_pu_ranges.insert((s, e));
                    }
                }
            }
            let unique_ranges: Vec<(i32, i32)> = all_pu_ranges.into_iter().collect();
            let cache =
                energy::compute_pu_energies_batch(ec, &pdb_path_str, &chain_str, &unique_ranges);
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
        !cli.energies,
        &results_dir.join("summary.txt"),
    )?;
    output::write_sword_summary_json(
        &sword_results,
        &energies,
        !cli.energies,
        &results_dir.join("summary.json"),
    )?;

    // Step 9: Write Peeling results
    tracing::debug!("Write Peeling results");
    let peeling_num = results_dir
        .join("intermediate")
        .join(format!("{}.num", pdb_id_chain));
    let peeling_log = results_dir.join("intermediate").join("peeling.log");

    if peeling_log.exists() && peeling_num.exists() {
        let peeling_results = peeling::load_legacy_results(&peeling_log, &peeling_num)?;

        // Calculate peeling energies if enabled — reuse global PU cache + compute missing
        let peeling_energies = if let Some(ref ec) = energy_config {
            // Collect unique PU ranges not already in the cache
            let mut missing_ranges: Vec<(i32, i32)> = Vec::new();
            let mut seen_pus: std::collections::HashSet<(i32, i32)> =
                std::collections::HashSet::new();
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
    reporter.step_done("Write results", None);

    // Domain PDB extraction (optimal partition only)
    if cli.extract_domains {
        if let Some(best) = sword_results.domains.first() {
            let domains_dir = results_dir.join("domains_optimal");
            pdb::writer::write_domain_pdbs(&cleaned_chain, &best.boundaries, &domains_dir)
                .context("Failed to write domain PDBs")?;
            // The directory name (domains_optimal) already conveys this; only
            // surface the note when the user asks for more verbose output.
            if cli.verbosity >= 1 {
                reporter.note(&format!(
                    "Extracted {} domain PDB(s) (optimal partition) → {}",
                    best.boundaries.len(),
                    domains_dir.display()
                ));
            }
        }
    }

    // Step 10: Generate plots
    if config.generate_plots {
        reporter.step("Generate plots");
        let plots_dir = results_dir.join("plots");
        std::fs::create_dir_all(&plots_dir)?;

        let proba_mat_file = results_dir.join("intermediate").join("contact_matrix.mat");

        // Domain consistency histogram (SVG)
        let histogram_output = plots_dir.join("domain_histogram.svg");
        let domain_counts = plot::count_domains(&sword_results.domains);
        if let Err(err) =
            plot::write_domain_histogram(&domain_counts, &histogram_output.to_string_lossy())
        {
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
                            &matrix, i, partition, &pu_colors, &plots_dir,
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
    let junctions_content = sword::junctions::calculate_junction_consistencies(&sword_output);
    if !junctions_content.is_empty() {
        std::fs::write(results_dir.join("junctions.txt"), &junctions_content)?;
    }

    // Step 12: Write mapping file
    pdb::writer::write_mapping_file(&original_resnums, &results_dir.join("residue_mapping.txt"))?;

    // Step 13: Clean up legacy artifacts
    tracing::debug!("Clean up results");
    let pdbs_stand = results_dir.join("PDBs_Stand");
    if pdbs_stand.exists() {
        let _ = std::fs::remove_dir_all(&pdbs_stand);
    }
    reporter.step_done("Junctions & cleanup", None);

    output::write_stdout_summary(
        &pdb_id_chain,
        chain_id,
        &sword_results,
        &energies,
        cli.format,
    );

    let elapsed = start.elapsed();
    reporter.finish(&pdb_id_chain, n_domains, prot_len, elapsed, &results_dir);

    Ok(())
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Subcommands short-circuit the full pipeline.
    if let Some(Commands::Score(args)) = &cli.command {
        return run_score(args);
    }

    setup_logging(cli.verbosity, cli.quiet);

    if cli.format == output::OutputFormat::Tsv {
        output::write_tsv_header();
    }

    let num_threads = if cli.threads == 0 {
        num_cpus::get()
    } else {
        cli.threads
    };

    let base_dir = if let Some(ref bd) = cli.install_dir {
        bd.clone()
    } else {
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

    let output_dir = std::fs::canonicalize(&cli.output_dir).unwrap_or_else(|_| {
        std::fs::create_dir_all(&cli.output_dir).ok();
        cli.output_dir.clone()
    });
    std::fs::create_dir_all(&output_dir)?;

    if let Some(ref batch_file) = cli.batch {
        // Batch mode: read entries from file, process each in turn
        let content = std::fs::read_to_string(batch_file)
            .with_context(|| format!("Failed to read batch file {}", batch_file.display()))?;
        let entries: Vec<EntryArgs> = content.lines().filter_map(parse_batch_line).collect();

        if entries.is_empty() {
            anyhow::bail!("Batch file contains no valid entries");
        }

        let mut n_ok = 0usize;
        let mut n_err = 0usize;
        for entry in &entries {
            let mut reporter = Reporter::new(cli.verbosity, cli.quiet);
            match process_entry(
                &entry,
                &cli,
                &mut reporter,
                &bin_dir,
                &output_dir,
                num_threads,
            ) {
                Ok(()) => n_ok += 1,
                Err(e) => {
                    n_err += 1;
                    if !cli.quiet {
                        eprintln!(" ✗  {}  · {}", entry.label(), e);
                    }
                }
            }
        }
        if !cli.quiet {
            eprintln!("\nBatch complete: {} ok, {} failed", n_ok, n_err);
        }
    } else {
        // Single mode: validate that exactly one input source is given
        if cli.pdb_id.is_none()
            && cli.uniprot_id.is_none()
            && cli.mgnify_id.is_none()
            && cli.input_file.is_none()
        {
            anyhow::bail!(
                "Please provide an input source: --pdb-id, --uniprot-id, --mgnify-id, \
                 --input-file, or --batch"
            );
        }
        let entry = EntryArgs {
            pdb_id: cli.pdb_id.clone(),
            uniprot_id: cli.uniprot_id.clone(),
            mgnify_id: cli.mgnify_id.clone(),
            input_file: cli.input_file.clone(),
            chain: cli.chain.as_ref().and_then(|s| s.chars().next()),
        };
        let mut reporter = Reporter::new(cli.verbosity, cli.quiet);
        process_entry(
            &entry,
            &cli,
            &mut reporter,
            &bin_dir,
            &output_dir,
            num_threads,
        )?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Input resolution
// ---------------------------------------------------------------------------

/// Resolve an entry to a local file path, a base name for output, and a fetch flag.
fn resolve_input(
    entry: &EntryArgs,
    output_dir: &PathBuf,
    legacy_pdb: bool,
) -> Result<(PathBuf, String, bool)> {
    if let Some(ref input_file) = entry.input_file {
        let base = input_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        Ok((input_file.clone(), base, false))
    } else if let Some(ref uniprot_id) = entry.uniprot_id {
        let path = fetch::fetch_alphafold(uniprot_id, output_dir)?;
        Ok((path, uniprot_id.clone(), true))
    } else if let Some(ref mgnify_id) = entry.mgnify_id {
        let path = fetch::fetch_esm(mgnify_id, output_dir)?;
        Ok((path, mgnify_id.clone(), true))
    } else if let Some(ref pdb_id) = entry.pdb_id {
        let path = fetch::fetch_pdb(pdb_id, output_dir, legacy_pdb)?;
        Ok((path, pdb_id.to_uppercase(), true))
    } else {
        anyhow::bail!("No input source specified");
    }
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::{error::ErrorKind, Parser};

    #[test]
    fn factorized_cli_defaults_off_and_explicit_flag_enables_it() {
        let default = Cli::try_parse_from(["sword2"]).unwrap();
        assert!(!default.use_factorized_ranker);

        let enabled = Cli::try_parse_from(["sword2", "--use-factorized-ranker"]).unwrap();
        assert!(enabled.use_factorized_ranker);
    }

    #[test]
    fn factorized_cli_conflicts_with_every_legacy_experimental_selector() {
        for conflicting in [
            vec!["--use-pairwise-reranker"],
            vec!["--use-count-calibration"],
            vec!["--count-lambda", "0.1"],
            vec!["--use-geometry-metrics"],
            vec!["--geometry-lambda", "0.1"],
        ] {
            let mut args = vec!["sword2", "--use-factorized-ranker"];
            args.extend(conflicting);
            let error = Cli::try_parse_from(args).unwrap_err();
            assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
        }
    }

    #[test]
    fn factorized_cli_does_not_expose_future_legacy_selector_flag() {
        let error = Cli::try_parse_from(["sword2", "--legacy-selector"]).unwrap_err();
        assert_eq!(error.kind(), ErrorKind::UnknownArgument);
    }
}
