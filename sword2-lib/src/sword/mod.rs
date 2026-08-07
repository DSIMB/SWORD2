//! SWORD algorithm orchestration.
//!
//! This module coordinates the full SWORD2 pipeline entirely in Rust:
//! 1. DSSP secondary structure assignment (pure Rust)
//! 2. Protein Peeling (pure Rust)
//! 3. ComputeMeasure (PU merging — pure Rust)
//! 4. ParseMeasure + prediction model (domain selection — pure Rust)
//! 5. Quality scoring and display (pure Rust)

use std::path::{Path, PathBuf};
use std::sync::LazyLock;

use anyhow::{Context, Result};
use regex::Regex;

static DIGITS_RE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d+").unwrap());
static AMB_RE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"^A-index = (\++)$").unwrap());
static ASS_RE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r"\d{1,}\s+\|").unwrap());

use crate::output::Partitioning;

pub(crate) mod candidate_features;
pub mod compute_jones;
pub mod compute_measure;
pub mod count_calibration;
pub mod distance_model;
pub(crate) mod factorized_ranker;
pub mod geometry_metrics;
pub mod junctions;
pub mod parse_measure;
pub(crate) mod reranker;

/// Configuration for a SWORD2 run.
#[derive(Debug, Clone)]
pub struct SwordConfig {
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
    /// Reduced-shuffle energy config used for candidate rescoring and the
    /// training dump — separate from the `-E` display-time energy config,
    /// which main.rs builds independently at the full shuffle count. `None`
    /// disables energy-based rescoring entirely (dump rows omit `energy_z`,
    /// `use_pairwise_reranker` falls back to distance_model-only scoring).
    pub energy_config: Option<crate::energy::EnergyConfig>,
    /// PDB chain letter (e.g. "A"), needed to build energy residue lists.
    pub chain_id: String,
    /// Use the pairwise-trained reranker to pick the winning candidate
    /// instead of the legacy distance_model-based selection. Off by default.
    pub use_pairwise_reranker: bool,
    /// Bias N_dom selection toward a length-predicted domain count.
    pub use_count_calibration: bool,
    /// Override the calibration penalty weight (for A/B sweeps); None = fitted default.
    pub count_lambda: Option<f64>,
    /// Reorder the winning candidate using the analytical "ideal sphere"
    /// geometry criteria (sphericity/density/interface-fraction z-scores).
    /// Off by default — geometry criteria are always computed and displayed
    /// regardless of this flag, but only change the pick when set.
    pub use_geometry_metrics: bool,
    /// Override the geometry selection penalty weight (for A/B sweeps);
    /// None = `DEFAULT_GEOMETRY_LAMBDA`.
    pub geometry_lambda: Option<f64>,
}

impl Default for SwordConfig {
    fn default() -> Self {
        Self {
            compute_energies: true,
            generate_plots: true,
            num_threads: num_cpus::get(),
            output_dir: ".".to_string(),
            max_alternatives: 9,
            energy_config: None,
            chain_id: "A".to_string(),
            use_pairwise_reranker: false,
            use_count_calibration: false,
            count_lambda: None,
            use_geometry_metrics: false,
            geometry_lambda: None,
        }
    }
}

/// Default weight for the geometry-score penalty in opt-in candidate
/// reordering (`score = dist_model - lambda * geometry_score`). Provisional —
/// unlike `count_calibration`'s fitted lambda, this hasn't had an A/B sweep
/// against a benchmark yet; override via `SwordConfig::geometry_lambda` or
/// `--geometry-lambda` while tuning.
const DEFAULT_GEOMETRY_LAMBDA: f64 = 0.1;

fn factorized_first_pass_indices(
    measures: &[compute_measure::MeasureLine],
    dir_data: &str,
    pdb_name: &str,
) -> Vec<usize> {
    parse_measure::parse_measure_indices(measures, dir_data, pdb_name, false, 0, 3, 3, true)
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
    /// Upper-tail p-value of the worst domain's sphericity (kappa^2) vs. the
    /// CATH reference; small = anomalously elongated. Always computed.
    pub sphericity_p: Option<f64>,
    /// Lower-tail p-value of the worst domain's Ca density vs. the CATH
    /// reference; small = anomalously loose/non-compact. Always computed.
    pub density_p: Option<f64>,
    /// Upper-tail p-value of the widest adjacent-domain interface fraction
    /// vs. the CATH reference; small = anomalously wide (likely over-split).
    /// `None` for single-domain partitions.
    pub interface_p: Option<f64>,
    /// Combined analytical geometry penalty (see `geometry_metrics::GeometryReport`).
    pub geometry_score: Option<f64>,
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
    pdb_name: &str,
    config: &SwordConfig,
) -> Result<(Vec<String>, SwordResults)> {
    let results_dir = PathBuf::from(&config.output_dir);
    let intermediate_dir = results_dir.join("intermediate");
    std::fs::create_dir_all(&intermediate_dir)?;

    // Copy PDB file to intermediate directory
    let pdb_file_src = pdb_path.to_path_buf();
    let pdb_file_dst = intermediate_dir.join(format!("{}.pdb", pdb_name));
    if pdb_file_src.exists() && !pdb_file_dst.exists() {
        std::fs::copy(&pdb_file_src, &pdb_file_dst)?;
    }

    // Step 1: Run DSSP (pure Rust)
    let dssp_file = intermediate_dir.join(format!("{}.dssp", pdb_name));
    let dssp_result = if !dssp_file.exists() {
        tracing::debug!("Running DSSP on {}", pdb_file_dst.display());
        let s2d_file = intermediate_dir.join(format!("{}.s2d", pdb_name));
        Some(crate::dssp::run_dssp(
            &pdb_file_dst,
            &dssp_file,
            &s2d_file,
            pdb_name,
        )?)
    } else {
        None
    };

    // Extract CA coordinates from the clean PDB. Needed for Peeling (when it
    // runs) and, regardless of whether Peeling was cached, for the
    // geometry_metrics criteria computed later — so this is hoisted out of
    // the cache-hit guard below rather than only computed on a Peeling run.
    let pdb_struct = crate::pdb::parse_pdb(&pdb_file_dst)
        .with_context(|| format!("Failed to parse clean PDB: {}", pdb_file_dst.display()))?;
    let ca_coords: Vec<[f64; 3]> = pdb_struct
        .first_model()
        .map(|m| {
            m.chains
                .iter()
                .flat_map(|c| c.residues.iter())
                .filter_map(|r| r.get_ca())
                .map(|a| [a.coord.x, a.coord.y, a.coord.z])
                .collect()
        })
        .unwrap_or_default();

    // Step 2: Run Peeling (native Rust)
    let pu_delineation_file = intermediate_dir.join("pu_delineation.mtx");
    let peeling_output = if !pu_delineation_file.exists() {
        tracing::debug!("Running Peeling on {}", pdb_file_dst.display());

        let peeling_config = crate::peeling::PeelingConfig::default();
        let output = crate::peeling::run_peeling(&ca_coords, &dssp_file, &peeling_config)?;

        // Write files for downstream compatibility (plots, debugging)
        output.write_peeling_log(&intermediate_dir.join("peeling.log"))?;
        output.write_pu_contact_matrix(&intermediate_dir.join("pu_contact.mtx"))?;
        output.write_pu_delineation(&pu_delineation_file)?;
        output
            .contact_matrix
            .write_matrix_file(&intermediate_dir.join("contact_matrix.mat"))?;

        Some(output)
    } else {
        None
    };

    // Step 3: Prepare .num file (residue number mapping)
    // The .num file should already exist from main.rs clean_chain_for_sword
    // but if not, create a default one from the PDB
    let num_file = intermediate_dir.join(format!("{}.num", pdb_name));
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
    let has_peeling = peeling_output.is_some() || pu_delineation_file.exists();
    if !has_peeling {
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
                sphericity_p: None,
                density_p: None,
                interface_p: None,
                geometry_score: None,
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

    // Run ComputeMeasure (in-memory if peeling output available, file-based if cached)
    tracing::debug!("Computing criteria for PUs merging");
    let measure_corpus = peeling_output.as_ref().map(|po| {
        compute_measure::compute_measure_from_data_with_provenance(
            &po.final_pu_contacts,
            &po.final_pu_delineation,
        )
    });
    let cached_measure_lines = measure_corpus.is_none().then(|| {
        let contact_matrix_file = intermediate_dir.join("pu_contact.mtx");
        compute_measure::compute_measure(&contact_matrix_file, &pu_delineation_file, 0.0001)
    });
    let measure_lines = measure_corpus
        .as_ref()
        .map(|corpus| corpus.lines.as_slice())
        .or_else(|| cached_measure_lines.as_deref())
        .expect("a fresh or cached measure corpus is always available");

    // DEBUG: count measure lines per domain count
    {
        let mut counts: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        for ml in measure_lines {
            *counts.entry(ml.num_domains).or_insert(0) += 1;
        }
        tracing::info!("ComputeMeasure lines per domain count: {:?}", counts);
    }

    // Legacy first-pass selection follows the configured alternative settings.
    let legacy_first_pass_indices = parse_measure::parse_measure_indices(
        &measure_lines,
        &results_dir.join("intermediate").to_string_lossy(),
        pdb_name,
        false,
        0,
        alt_b,
        alt_l,
        true,
    );
    let relevant_measure: Vec<String> = legacy_first_pass_indices
        .iter()
        .filter_map(|&index| measure_lines.get(index))
        .map(compute_measure::MeasureLine::to_line)
        .collect();

    // The factorized lattice is always built from its fixed (3, 3) first
    // pass, independently of the legacy output-alternative configuration.
    let factorized_indices = factorized_first_pass_indices(
        &measure_lines,
        &results_dir.join("intermediate").to_string_lossy(),
        pdb_name,
    );
    let mut candidate_lattice = factorized_ranker::lattice::CandidateLattice::from_first_pass(
        &measure_lines,
        &factorized_indices,
        tab_num.len(),
    )
    .map_err(|error| anyhow::anyhow!(error))?;
    let _factorized_context = match (
        dssp_result.as_ref(),
        peeling_output.as_ref(),
        measure_corpus.as_ref(),
    ) {
        (Some(dssp), Some(peeling), Some(corpus)) => Some(
            factorized_ranker::prepare_factorized_context(
                Some(&ca_coords),
                Some(&dssp.chain),
                &peeling.iterations,
                Some((&peeling.contact_matrix, &corpus.provenance)),
            )
            .and_then(|context| {
                candidate_lattice
                    .attach_hierarchy(context.measure_provenance, context.iterations)?;
                Ok(context)
            }),
        ),
        _ => None,
    };

    let first_pass_measures: Vec<compute_measure::MeasureLine> = legacy_first_pass_indices
        .iter()
        .filter_map(|&index| measure_lines.get(index).cloned())
        .collect();
    let legacy_selection = factorized_ranker::lattice::select_legacy(
        &first_pass_measures,
        tab_num.len(),
        config.use_count_calibration,
        config.count_lambda,
    );
    let n_dom = legacy_selection.num_domains;
    let to_print_first_pass = legacy_selection.measure_line;
    tracing::info!("Distance-model n_dom selection: n_dom={}", n_dom);

    // Second ParseMeasure pass — select assignments around predicted N_dom
    tracing::debug!("Selecting domain assignments around N_dom={}", n_dom);

    // DEBUG: dump relevant_measure from first pass
    {
        let mut counts: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        for rm in &relevant_measure {
            let nd: usize = rm
                .split('|')
                .next()
                .unwrap_or("0")
                .trim()
                .parse()
                .unwrap_or(0);
            *counts.entry(nd).or_insert(0) += 1;
        }
        tracing::info!(
            "First ParseMeasure: {} lines, per level: {:?}",
            relevant_measure.len(),
            counts
        );
    }

    let relevant_measure2 = parse_measure::parse_measure(
        &measure_lines,
        &results_dir.join("intermediate").to_string_lossy(),
        pdb_name,
        true,
        n_dom,
        alt_b,
        alt_l,
        true,
    );

    // Analytical "ideal sphere" geometry criteria (sphericity, density,
    // inter-domain interface fraction): always computable from ca_coords
    // alone, used purely as a decision aid for display unless
    // `use_geometry_metrics` opts into using them to reorder the pick.
    // Independent of, and untouched by, the pairwise reranker branch below.
    let geometry_reference = geometry_metrics::ReferenceDistributions::embedded();

    // Training dump: set SWORD2_DUMP_CANDIDATES=/path/to/output.csv to record the
    // same alt_b/alt_l shortlist the pairwise reranker scores at inference time,
    // with the same features (boundary_coil_fraction, energy_z, modal_count_distance)
    // for offline training. Runs on relevant_measure2 (post-shortlist) rather than
    // the raw measure_lines so training and inference see identical feature
    // distributions.
    if let Ok(dump_path) = std::env::var("SWORD2_DUMP_CANDIDATES") {
        use std::io::Write as _;

        let ss_types: &[crate::peeling::algorithm::SsType] = peeling_output
            .as_ref()
            .map(|po| po.ss_types.as_slice())
            .unwrap_or(&[]);
        let pdb_path_str = pdb_file_dst.to_string_lossy().to_string();

        let mut parsed_rows: Vec<(usize, usize, f64, f64, f64, String, String)> = Vec::new();
        for rm in &relevant_measure2 {
            let fields: Vec<&str> = rm.split('|').collect();
            if fields.len() < 7 {
                continue;
            }
            let nd: usize = fields[0].trim().parse().unwrap_or(0);
            let min_size: usize = fields[1].trim().parse().unwrap_or(0);
            let raw_del = fields[2].trim().to_string();
            let max_cr: f64 = fields[3].trim().parse().unwrap_or(0.0);
            let density_min: f64 = fields[5].trim().parse().unwrap_or(0.0);
            let mean_density: f64 = fields[6].trim().parse().unwrap_or(0.0);
            let remapped_del = remap_residue_numbers(&raw_del, &tab_num);
            parsed_rows.push((
                nd,
                min_size,
                max_cr,
                density_min,
                mean_density,
                raw_del,
                remapped_del,
            ));
        }

        let modal = candidate_features::modal_num_domains(
            &parsed_rows.iter().map(|r| r.0).collect::<Vec<_>>(),
        );

        let file_existed = std::path::Path::new(&dump_path).exists();
        if let Ok(mut f) = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&dump_path)
        {
            if !file_existed {
                let _ = writeln!(
                    f,
                    "chain_id,output_dir,num_domains,min_size,max_cr,density_min,mean_density,delineation,boundary_coil_fraction,energy_z,modal_count_distance"
                );
            }
            for (nd, min_size, max_cr, density_min, mean_density, raw_del, remapped_del) in
                &parsed_rows
            {
                let row = candidate_features::build_dump_row(
                    *nd,
                    *min_size,
                    *max_cr,
                    *density_min,
                    *mean_density,
                    raw_del,
                    remapped_del,
                    ss_types,
                    config.energy_config.as_ref(),
                    &pdb_path_str,
                    &config.chain_id,
                    modal,
                );
                let energy_z_str = row.energy_z.map(|z| z.to_string()).unwrap_or_default();
                let _ = writeln!(
                    f,
                    "{},{},{},{},{:.6},{:.6},{:.6},\"{}\",{:.6},{},{:.1}",
                    pdb_name,
                    results_dir.display(),
                    row.num_domains,
                    row.min_size,
                    row.max_cr,
                    row.density_min,
                    row.mean_density,
                    row.delineation,
                    row.boundary_coil_fraction,
                    energy_z_str,
                    row.modal_count_distance,
                );
            }
        }
    }

    // Find the measure line for the predicted N_dom
    // DEBUG: dump relevant_measure2 from second pass
    {
        tracing::info!(
            "Second ParseMeasure: {} lines total",
            relevant_measure2.len()
        );
        for (i, rm) in relevant_measure2.iter().enumerate() {
            let fields: Vec<&str> = rm.split('|').collect();
            let nd: usize = fields[0].trim().parse().unwrap_or(0);
            let del = if fields.len() > 2 {
                fields[2].trim()
            } else {
                "?"
            };
            let cr: f64 = if fields.len() > 3 {
                fields[3].trim().parse().unwrap_or(0.0)
            } else {
                0.0
            };
            let cpd: f64 = if fields.len() > 5 {
                fields[5].trim().parse().unwrap_or(0.0)
            } else {
                0.0
            };
            let dist_signed = crate::sword::distance_model::distance_model(cr, cpd, 1);
            tracing::info!(
                "  [{}] nd={} dist_signed={:.4} abs={:.4} del={}",
                i,
                nd,
                dist_signed,
                dist_signed.abs(),
                del
            );
        }
    }

    let mut to_print = if config.use_pairwise_reranker && relevant_measure2.len() > 1 {
        let ss_types: &[crate::peeling::algorithm::SsType] = peeling_output
            .as_ref()
            .map(|po| po.ss_types.as_slice())
            .unwrap_or(&[]);
        let pdb_path_str = pdb_file_dst.to_string_lossy().to_string();

        let mut raw_dels: Vec<String> = Vec::new();
        let mut remapped_dels: Vec<String> = Vec::new();
        let mut parsed: Vec<(usize, usize, f64, f64, f64)> = Vec::new();
        for rm in &relevant_measure2 {
            let fields: Vec<&str> = rm.split('|').collect();
            if fields.len() < 7 {
                continue;
            }
            let nd: usize = fields[0].trim().parse().unwrap_or(0);
            let min_size: usize = fields[1].trim().parse().unwrap_or(0);
            let raw_del = fields[2].trim().to_string();
            let max_cr: f64 = fields[3].trim().parse().unwrap_or(0.0);
            let density_min: f64 = fields[5].trim().parse().unwrap_or(0.0);
            let mean_density: f64 = fields[6].trim().parse().unwrap_or(0.0);
            remapped_dels.push(remap_residue_numbers(&raw_del, &tab_num));
            raw_dels.push(raw_del);
            parsed.push((nd, min_size, max_cr, density_min, mean_density));
        }

        let inputs: Vec<reranker::CandidateInput> = parsed
            .iter()
            .zip(raw_dels.iter())
            .zip(remapped_dels.iter())
            .map(
                |((&(nd, min_size, max_cr, density_min, mean_density), raw), remapped)| {
                    reranker::CandidateInput {
                        num_domains: nd,
                        min_size,
                        max_cr,
                        density_min,
                        mean_density,
                        raw_delineation: raw,
                        remapped_delineation: remapped,
                    }
                },
            )
            .collect();

        if inputs.is_empty() {
            String::new()
        } else {
            let winner = reranker::rerank(
                &inputs,
                ss_types,
                config.energy_config.as_ref(),
                &pdb_path_str,
                &config.chain_id,
            );
            relevant_measure2[winner].clone()
        }
    } else {
        String::new()
    };

    if to_print.is_empty() && config.use_geometry_metrics && relevant_measure2.len() > 1 {
        // Analytical opt-in reordering: score = dist_model - lambda * G,
        // independent of (and evaluated after) the pairwise reranker above.
        let lambda_geo = config.geometry_lambda.unwrap_or(DEFAULT_GEOMETRY_LAMBDA);
        let mut best_idx: Option<usize> = None;
        let mut best_score = f64::NEG_INFINITY;
        for (i, rm) in relevant_measure2.iter().enumerate() {
            let fields: Vec<&str> = rm.split('|').collect();
            if fields.len() < 7 {
                continue;
            }
            let cr: f64 = fields[3].trim().parse().unwrap_or(0.0);
            let cpd: f64 = fields[5].trim().parse().unwrap_or(0.0);
            let dist = distance_model::distance_model(cr, cpd, 1);
            let penalty = geometry_report_for_raw_delineation(fields[2].trim(), &ca_coords, geometry_reference)
                .map(|r| r.geometry_score)
                .unwrap_or(0.0);
            let score = dist - lambda_geo * penalty;
            if score > best_score {
                best_score = score;
                best_idx = Some(i);
            }
        }
        if let Some(idx) = best_idx {
            to_print = relevant_measure2[idx].clone();
        }
    }

    if to_print.is_empty() {
        // Legacy distance_model-based selection (also the fallback when the
        // reranker is disabled, has <2 candidates, or all candidates failed
        // to parse above).
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
    }
    // Fallback: the second-pass distance filter (dist < 0.2) excludes deeply good-zone
    // candidates. If the optimal was filtered out, use the first-pass representative.
    if to_print.is_empty() {
        tracing::debug!(
            "to_print empty after second pass; using first-pass representative for nd={}",
            n_dom
        );
        to_print = to_print_first_pass;
    }

    // Quality and display
    let geometry_ctx = GeometryContext { ca_coords: &ca_coords, reference: geometry_reference };
    let output_lines = quality_and_display(
        &tab_num,
        &to_print,
        n_dom,
        &relevant_measure2,
        alt_b,
        alt_l,
        &geometry_ctx,
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

    let last = relevant_measure.len().saturating_sub(1);
    for line_str in &relevant_measure[..last] {
        let fields: Vec<&str> = line_str.split('|').collect();
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

/// Bundles what the geometry criteria need to be computed during display:
/// Cα coordinates plus the reference distributions to score against.
struct GeometryContext<'a> {
    ca_coords: &'a [[f64; 3]],
    reference: &'a geometry_metrics::ReferenceDistributions,
}

/// Generate quality and display output lines.
///
/// Port of `quality_and_display()` from SWORD Perl script.
fn quality_and_display(
    tab_num: &[i32],
    to_print: &str,
    n_dom: usize,
    relevant_measure: &[String],
    alt_b: usize,
    alt_l: usize,
    geometry_ctx: &GeometryContext,
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

    let geometry = geometry_report_for_raw_delineation(
        to_print_fields[2].trim(),
        geometry_ctx.ca_coords,
        geometry_ctx.reference,
    );

    // Output: ambiguity will be computed at the end
    // Header
    output.push(format!(
        "{:<2}|{:<3}|{:>60}|{:>12.6}|{:>10}|{}|",
        to_print_fields[0].trim(),
        to_print_fields[1].trim(),
        delineation,
        avg_k,
        calc_print,
        format_geometry_fields(geometry.as_ref()),
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
                    let alt_geometry = geometry_report_for_raw_delineation(
                        fields[2].trim(),
                        geometry_ctx.ca_coords,
                        geometry_ctx.reference,
                    );
                    output.push(format!(
                        "{:<2}|{:<3}|{:>60}|{:>12.6}|{:>10}|{}|",
                        fields[0].trim(),
                        fields[1].trim(),
                        alt_del,
                        this_k,
                        alt_quality,
                        format_geometry_fields(alt_geometry.as_ref()),
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

/// Parse a raw (0-based) delineation string and compute its geometry report,
/// or `None` if there are no usable Cα coordinates or domains to score.
fn geometry_report_for_raw_delineation(
    raw_delineation: &str,
    ca_coords: &[[f64; 3]],
    reference: &geometry_metrics::ReferenceDistributions,
) -> Option<geometry_metrics::GeometryReport> {
    if ca_coords.is_empty() {
        return None;
    }
    let domains = geometry_metrics::parse_domain_indices(raw_delineation);
    if domains.is_empty() {
        return None;
    }
    Some(geometry_metrics::geometry_report(ca_coords, &domains, reference))
}

/// Format a geometry report as four `|`-delimited fields (sphericity_p,
/// density_p, interface_p, geometry_score), `n/a` where not applicable —
/// appended to display lines alongside the existing quality star rating.
fn format_geometry_fields(report: Option<&geometry_metrics::GeometryReport>) -> String {
    match report {
        Some(r) => format!(
            "{:.4}|{:.4}|{}|{:.4}",
            r.sphericity_p(),
            r.density_p(),
            r.interface_p()
                .map(|p| format!("{:.4}", p))
                .unwrap_or_else(|| "n/a".to_string()),
            r.geometry_score,
        ),
        None => "n/a|n/a|n/a|n/a".to_string(),
    }
}

/// Remap renumbered residue positions in a delineation string to original numbering.
///
/// The delineation contains 0-based indices from file_pu_delineation.mtx.
/// `tab_num` is 0-indexed: tab_num[i] = original residue number for index i.
/// This matches the Perl original: `$$ref_tab_num[$1]` (direct 0-based access).
fn remap_residue_numbers(delineation: &str, tab_num: &[i32]) -> String {
    DIGITS_RE
        .replace_all(delineation, |caps: &regex::Captures| {
            let idx: usize = caps[0].parse().unwrap_or(0);
            if idx < tab_num.len() {
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
    if c5 >= 5 {
        "+++++".to_string()
    } else if c4 >= 4 {
        "++++".to_string()
    } else if c3 >= 3 {
        "+++".to_string()
    } else if c2 >= 2 {
        "++".to_string()
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
            // Geometry fields (sphericity_p, density_p, interface_p, geometry_score) are
            // appended after quality by quality_and_display; absent on older-format lines.
            let sphericity_p = parts.get(5).and_then(|s| s.parse::<f64>().ok());
            let density_p = parts.get(6).and_then(|s| s.parse::<f64>().ok());
            let interface_p = parts.get(7).and_then(|s| s.parse::<f64>().ok());
            let geometry_score = parts.get(8).and_then(|s| s.parse::<f64>().ok());

            domains.push(SwordPartition {
                nb_domains,
                min_size,
                boundaries,
                average_k,
                quality,
                sphericity_p,
                density_p,
                interface_p,
                geometry_score,
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
                sphericity_p: part.sphericity_p,
                density_p: part.density_p,
                interface_p: part.interface_p,
                geometry_score: part.geometry_score,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn measure_line(num_domains: usize) -> compute_measure::MeasureLine {
        compute_measure::MeasureLine {
            num_domains,
            min_size: 1,
            delineation: "0-1".to_string(),
            max_cr: 0.1,
            mean_cr: 0.0,
            density_min: 3.0,
            mean_density: 3.0,
        }
    }

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
        let output = vec!["  2 | 30 | 1-50;151-200 51-150 | 3.5 | ****".to_string()];
        let results = parse_sword_output(&output).unwrap();
        assert_eq!(results.domains.len(), 1);
        assert_eq!(results.domains[0].boundaries[0], vec![(1, 50), (151, 200)]);
        assert_eq!(results.domains[0].boundaries[1], vec![(51, 150)]);
    }

    #[test]
    fn test_parse_sword_no_ambiguity() {
        let output = vec!["  2 | 30 | 1-100 101-200 | 3.5 | *****".to_string()];
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

    #[test]
    fn count_calibration_defaults_to_off() {
        let config = SwordConfig::default();
        assert!(!config.use_count_calibration);
        assert_eq!(config.count_lambda, None);
    }

    #[test]
    fn factorized_first_pass_uses_fixed_shortlist_parameters() {
        let measures: Vec<compute_measure::MeasureLine> = [7, 6, 5, 4, 3, 2, 1]
            .into_iter()
            .flat_map(|num_domains| (0..4).map(move |_| measure_line(num_domains)))
            .collect();
        let fixed = factorized_first_pass_indices(&measures, "", "test");
        let configured = parse_measure::parse_measure_indices(
            &measures, "", "test", false, 0, 1, 3, true,
        );
        assert_ne!(fixed, configured);
    }

    #[test]
    fn calibration_off_matches_raw_distance_argmax() {
        use crate::sword::count_calibration::CountCalibration;
        use crate::sword::distance_model::distance_model;

        // three level representatives: (nd, cr, cpd)
        let reps = [(2usize, 0.15, 3.3), (3, 0.21, 3.1), (5, 0.46, 2.7)];

        // raw argmax (what ships today)
        let raw_best = reps
            .iter()
            .max_by(|a, b| {
                distance_model(a.1, a.2, 1)
                    .partial_cmp(&distance_model(b.1, b.2, 1))
                    .unwrap()
            })
            .unwrap()
            .0;

        // lambda = 0 must reproduce the raw argmax exactly
        let c = CountCalibration {
            intercept: 1.2,
            len_coef: 0.0035,
            lambda: 0.0,
        };
        let expected = c.expected_num_domains(300);
        let calib_best = reps
            .iter()
            .max_by(|a, b| {
                c.adjusted_score(distance_model(a.1, a.2, 1), a.0, expected)
                    .partial_cmp(&c.adjusted_score(distance_model(b.1, b.2, 1), b.0, expected))
                    .unwrap()
            })
            .unwrap()
            .0;

        assert_eq!(raw_best, calib_best);
    }
}
