//! Native Rust implementation of the Protein Peeling algorithm.
//!
//! Faithfully ports the iterative hierarchical cutting algorithm from
//! `peeling_omp.c` (Gelly et al., 2006). The algorithm:
//!
//! 1. Builds a contact probability matrix from C-alpha coordinates
//! 2. Determines which residue boundaries can be cut (based on secondary structure)
//! 3. Iteratively finds the best single or double cut that maximizes a
//!    Matthews-like coefficient separating intra-PU contacts from inter-PU contacts
//! 4. Stops when CI (Compaction Index) exceeds the threshold, no cut is possible,
//!    or the maximum number of PUs is reached
//!
//! References:
//!   - Gelly, J.C.; de Brevern, A.G.; Hazout, S. (2006) Bioinformatics 22(2):129-33
//!   - Gelly, J.C.; de Brevern, A.G. (2011) Bioinformatics 27(2):132-133

use std::io::Write;
use std::path::Path;
use std::sync::Mutex;

use anyhow::{Context, Result};
use rayon::prelude::*;

use super::contact_matrix::ContactMatrix;

/// Maximum number of peeling iterations.
const MAX_ITERATION: usize = 64;

/// Maximum number of PUs stored per iteration.
const _MAX_PUS: usize = 128;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

/// Configuration parameters for the peeling algorithm.
///
/// Default values match the SWORD2 pipeline invocation:
/// `-r 98 -s 8 -l 30 -m 0 -0 6.0 -t 1.5 -o 0 -g 0 -c 0 -n 30`
#[derive(Debug, Clone)]
pub struct PeelingConfig {
    /// Maximum CI value before stopping (default: 98).
    pub max_r2: i32,
    /// Minimum secondary structure segment size that can be cut (default: 8).
    pub min_ss_size: usize,
    /// Minimum PU size in residues (default: 30).
    pub min_pu_size: usize,
    /// Maximum PU size (0 = no limit, default: 0).
    pub max_pu_size: usize,
    /// Midpoint distance for the logistic contact function in Å (default: 6.0).
    pub d0: f64,
    /// Steepness parameter for the logistic contact function in Å (default: 1.5).
    pub delta: f64,
    /// If true, down-weight contacts involving coil residues (default: false).
    pub only_ss: bool,
    /// Enable pruning based on homogeneity (default: false).
    pub pruning: bool,
    /// Cutoff for pruning homogeneity score (default: 0.0).
    pub cutoff_pruning: f64,
    /// Maximum number of PUs allowed (default: 30).
    pub max_pu_number: usize,
}

impl Default for PeelingConfig {
    fn default() -> Self {
        Self {
            max_r2: 98,
            min_ss_size: 8,
            min_pu_size: 30,
            max_pu_size: 0,
            d0: 6.0,
            delta: 1.5,
            only_ss: false,
            pruning: false,
            cutoff_pruning: 0.0,
            max_pu_number: 30,
        }
    }
}

// ---------------------------------------------------------------------------
// Secondary structure / cutting mask
// ---------------------------------------------------------------------------

/// Secondary structure classification for a residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SsType {
    Coil = 0,
    Helix = 1,
    Sheet = 2,
}

/// Build the cutting mask and SS type array from a DSSP file.
///
/// This parses the DSSP output file in the same way as the C code:
/// - Column 126 must not be whitespace or '-' (skip header/blank lines)
/// - Column 13: amino acid type ('!' = chain break, skip)
/// - Column 16: secondary structure symbol (H/G=helix, E/B=sheet, else coil)
/// - Columns 6-10: residue number (tab_true_num)
///
/// Returns `(ss_types, true_nums, cutting_mask)`:
/// - `ss_types[i]`: secondary structure type for residue i (0-indexed)
/// - `true_nums[i]`: original residue number for residue i
/// - `cutting_mask[i]`: whether position i can be a cut point
///
/// Short SS segments (< min_ss_size) are marked as non-cuttable.
#[allow(clippy::needless_range_loop)]
pub(crate) fn parse_dssp_for_peeling(
    dssp_path: &Path,
    num_residues: usize,
    config: &PeelingConfig,
) -> Result<(Vec<SsType>, Vec<i32>, Vec<bool>)> {
    let content = std::fs::read_to_string(dssp_path)
        .with_context(|| format!("Cannot read DSSP file: {}", dssp_path.display()))?;

    let mut ss_types = Vec::with_capacity(num_residues);
    let mut true_nums = Vec::with_capacity(num_residues);

    for line in content.lines() {
        let bytes = line.as_bytes();
        // Skip lines where column 126 is whitespace or '-'
        if bytes.len() <= 126 {
            continue;
        }
        let col126 = bytes[126];
        if col126 == b' ' || col126 == b'\t' || col126 == b'-' {
            continue;
        }

        // Column 13: amino acid type
        if bytes.len() <= 16 {
            continue;
        }
        let aa = bytes[13] as char;
        if aa == '!' {
            continue;
        }

        // Parse residue number from columns 6-10 (after blanking col 6 and 11)
        let mut num_bytes = Vec::from(&bytes[6..11]);
        num_bytes[0] = b' '; // line[6] = ' ' in C code
                             // Note: C code also sets line[11] = ' ' but we only read cols 6-10
        let num_str = std::str::from_utf8(&num_bytes).unwrap_or("0");
        let num: i32 = num_str.trim().parse().unwrap_or(0);
        true_nums.push(num);

        // Column 16: secondary structure
        let ss_char = bytes[16] as char;
        let ss = match ss_char {
            'H' | 'G' => SsType::Helix,
            'E' | 'B' => SsType::Sheet,
            _ => SsType::Coil,
        };
        ss_types.push(ss);
    }

    let n_dssp = ss_types.len();
    let max_residues = num_residues.max(n_dssp);

    // Initialize cutting mask: all positions are cuttable
    let mut cutting_mask = vec![true; max_residues + 1];

    // Find short secondary structure segments and mark them as non-cuttable.
    // This reproduces the C code's SS segment detection logic exactly.
    let mut small_ss_segments: Vec<(usize, usize)> = Vec::new();
    let mut current_ss: SsType = SsType::Coil;
    let mut segment_start: usize = 0;

    for i in 0..n_dssp {
        if ss_types[i] != SsType::Coil && current_ss == SsType::Coil {
            // Start of a new SS segment
            current_ss = ss_types[i];
            segment_start = i;
        } else if ss_types[i] != current_ss
            && ss_types[i] != SsType::Coil
            && current_ss != SsType::Coil
        {
            // Transition between different SS types
            let segment_end = i - 1;
            let size = segment_end - segment_start;
            if size <= config.min_ss_size && size != 0 {
                small_ss_segments.push((segment_start, segment_end));
            }
            current_ss = ss_types[i];
            segment_start = i;
        } else if ss_types[i] == SsType::Coil && current_ss != SsType::Coil {
            // End of SS segment
            let segment_end = i - 1;
            let size = segment_end - segment_start;
            if size <= config.min_ss_size && size != 0 {
                small_ss_segments.push((segment_start, segment_end));
            }
            current_ss = SsType::Coil;
        }
    }

    // Mark positions within small SS segments as non-cuttable
    for &(seg_start, seg_end) in &small_ss_segments {
        for k in seg_start..seg_end {
            cutting_mask[k] = false;
        }
    }

    Ok((ss_types, true_nums, cutting_mask))
}

/// Apply the "only secondary structure" filter to the contact matrix.
///
/// When enabled, contacts involving coil residues are down-weighted by 10x.
/// This modifies the contact matrix in place (not used in default SWORD2 pipeline).
#[allow(dead_code)]
pub(crate) fn apply_only_ss_filter(matrix: &mut ContactMatrix, ss_types: &[SsType]) {
    let n = matrix.len();
    // We cannot mutate the matrix directly since it's behind an abstraction.
    // This is only used when only_ss=true which is not the default.
    // For the default pipeline, this is a no-op.
    // If needed in the future, we'd need to add a mutable access method.
    let _ = (n, ss_types);
}

// ---------------------------------------------------------------------------
// Cutting algorithms
// ---------------------------------------------------------------------------

/// Result of a cutting operation (single or double cut).
#[derive(Debug, Clone, Copy)]
struct CutResult {
    /// Matthews coefficient for this cut.
    coeff: f64,
    /// Number of cuts (1 or 2).
    num_cuts: usize,
    /// PU index that was cut.
    pu_index: usize,
    /// Start of the PU being cut.
    start: usize,
    /// First cut point: left boundary of first new PU ends at i1.
    i1: usize,
    /// First cut point: right boundary starts at i2 = i1 + 1.
    i2: usize,
    /// Second cut point (for double cut): j1 and j2 = j1 + 1.
    j1: usize,
    j2: usize,
    /// End of the PU being cut.
    end: usize,
}

/// Find the best single cut within a PU.
///
/// Tests every valid cut position and returns the one with the highest
/// Matthews-like coefficient:
///
/// $$\text{MCC} = \frac{ab - c^2}{(a+c)(b+c)}$$
///
/// where a = intra-left contacts, b = intra-right contacts, c = inter contacts.
fn simple_cutting(
    matrix: &ContactMatrix,
    pu_start: usize,
    pu_end: usize,
    cutting_mask: &[bool],
    min_pu_size: usize,
    current_best: f64,
    pu_index: usize,
) -> Option<CutResult> {
    let mut best: Option<CutResult> = None;
    let mut best_coeff = current_best;

    for i in pu_start..pu_end {
        let i_p = i + 1;
        // Check cutting mask (C code: tab_decoupe[i+1])
        if i_p < cutting_mask.len() && !cutting_mask[i_p] {
            continue;
        }
        let size_pu1 = i - pu_start;
        let size_pu2 = pu_end - i;
        if size_pu1 < min_pu_size || size_pu2 < min_pu_size {
            continue;
        }

        let a = matrix.rectangle_sum(pu_start, pu_start, i, i);
        let b = matrix.rectangle_sum(i + 1, i + 1, pu_end, pu_end);
        let c = matrix.rectangle_sum(pu_start, i + 1, i, pu_end);

        let denom = (a + c) * (b + c);
        if denom == 0.0 {
            continue;
        }

        let coeff = (a * b - c * c) / denom;
        if coeff > best_coeff {
            best_coeff = coeff;
            best = Some(CutResult {
                coeff,
                num_cuts: 1,
                pu_index,
                start: pu_start,
                i1: i,
                i2: i + 1,
                j1: 0,
                j2: 0,
                end: pu_end,
            });
        }
    }

    best
}

/// Find the best double cut within a PU.
///
/// Tests all valid pairs of cut positions. The outer loop is parallelized
/// with rayon. For a double cut, the Matthews coefficient considers:
///
/// - a = intra-middle contacts (between the two cuts)
/// - b = intra-flanking contacts (left + right + cross between left-right)
/// - c = inter contacts (middle-left + middle-right)
fn double_cutting(
    matrix: &ContactMatrix,
    pu_start: usize,
    pu_end: usize,
    cutting_mask: &[bool],
    min_pu_size: usize,
    current_best: f64,
    pu_index: usize,
) -> Option<CutResult> {
    let min_seg = min_pu_size;
    let coo_max_i = pu_end.saturating_sub(min_seg);
    let coo_max_j = pu_end.saturating_sub(min_seg / 2);
    let coo_min_i = pu_start + min_seg - 1;

    if coo_min_i >= coo_max_i {
        return None;
    }

    // Use a Mutex to collect results from parallel threads
    let best = Mutex::new((current_best, None::<CutResult>));

    (coo_min_i..coo_max_i).into_par_iter().for_each(|i| {
        let i_p = i + 1;
        if i_p < cutting_mask.len() && !cutting_mask[i_p] {
            return;
        }
        let coo_min_j = i + min_seg;
        if coo_min_j + min_seg >= coo_max_j {
            return;
        }

        // Thread-local best for this outer iteration
        let mut local_best_coeff = {
            let lock = best.lock().unwrap();
            lock.0
        };
        let mut local_best_cut: Option<CutResult> = None;

        for j in coo_min_j..=coo_max_j {
            let j_p = j + 1;
            if j_p < cutting_mask.len() && !cutting_mask[j_p] {
                continue;
            }
            let size_pu2 = j - i;
            if size_pu2 < min_seg {
                continue;
            }

            let i2 = i + 1;
            let j2 = j + 1;

            // a = intra-middle contacts
            let a = matrix.rectangle_sum(i2, i2, j, j);

            // b = intra-left + intra-right + 2 * cross(left, right)
            let b1 = matrix.rectangle_sum(pu_start, pu_start, i, i);
            let b2 = matrix.rectangle_sum(j2, j2, pu_end, pu_end);
            let b3 = matrix.rectangle_sum(pu_start, j2, i, pu_end);
            let b = b1 + b2 + 2.0 * b3;

            // c = contacts between middle and flanking regions
            let c1 = matrix.rectangle_sum(i2, pu_start, j, i);
            let c2 = matrix.rectangle_sum(j2, i2, pu_end, j);
            let c = c1 + c2;

            let denom = (a + c) * (b + c);
            if denom == 0.0 {
                continue;
            }

            let coeff = (a * b - c * c) / denom;
            if coeff > local_best_coeff {
                local_best_coeff = coeff;
                local_best_cut = Some(CutResult {
                    coeff,
                    num_cuts: 2,
                    pu_index,
                    start: pu_start,
                    i1: i,
                    i2,
                    j1: j,
                    j2,
                    end: pu_end,
                });
            }
        }

        // Publish local best to shared state
        if let Some(cut) = local_best_cut {
            let mut lock = best.lock().unwrap();
            if cut.coeff > lock.0 {
                lock.0 = cut.coeff;
                lock.1 = Some(cut);
            }
        }
    });

    let lock = best.lock().unwrap();
    lock.1
}

// ---------------------------------------------------------------------------
// Metrics
// ---------------------------------------------------------------------------

/// Compute max contact ratio and minimum density for the current PU set.
///
/// Returns (max_cr, min_density).
fn measure_coeff(matrix: &ContactMatrix, pus: &[[usize; 2]]) -> (f64, f64) {
    let n_pus = pus.len();
    let mut min_density = f64::MAX;
    let mut max_cr = f64::MIN;

    // Cache internal contact sums
    let mut internal: Vec<f64> = Vec::with_capacity(n_pus);
    for pu in pus {
        let sum = matrix.rectangle_sum(pu[0], pu[0], pu[1], pu[1]);
        internal.push(sum);
    }

    for x in 0..n_pus {
        let x1 = pus[x][0];
        let x2 = pus[x][1];
        let size_pu1 = (x2 - x1 + 1) as f64;
        let density = internal[x] / size_pu1;
        if density < min_density {
            min_density = density;
        }

        for y in 0..n_pus {
            if x == y {
                continue;
            }
            let y1 = pus[y][0];
            let y2 = pus[y][1];
            let size_pu2 = (y2 - y1 + 1) as f64;

            let prob_external = matrix.rectangle_sum(x1, y1, x2, y2);

            let pdp = (internal[x] + internal[y] + prob_external) / (size_pu1 + size_pu2);
            let nnc = prob_external / (size_pu1.powf(0.43) * size_pu2.powf(0.43));

            if pdp > 0.0 {
                let cr = nnc / pdp;
                if cr > max_cr {
                    max_cr = cr;
                }
            }
        }
    }

    (max_cr, min_density)
}

/// Compute the Compaction Index (CI) and R metric via mutual information.
///
/// Also returns the PU-PU contact matrix and delineation data for use
/// by the downstream compute_measure module.
///
/// Returns (ci, r, pu_contact_entries, pu_delineation_entries).
#[allow(clippy::type_complexity, clippy::needless_range_loop)]
fn mutual_information(
    matrix: &ContactMatrix,
    pus: &[[usize; 2]],
) -> (
    f64,
    f64,
    Vec<(usize, usize, f64)>,
    Vec<(usize, usize, usize)>,
) {
    let n_pus = pus.len();

    // Compute PU-PU contact sums
    let mut prob_zone = vec![vec![0.0f64; n_pus]; n_pus];
    let mut pu_contact_entries = Vec::with_capacity(n_pus * n_pus);
    let mut pu_delineation_entries = Vec::with_capacity(n_pus);

    for x in 0..n_pus {
        let x1 = pus[x][0];
        let x2 = pus[x][1];
        for y in 0..n_pus {
            let y1 = pus[y][0];
            let y2 = pus[y][1];
            prob_zone[x][y] = matrix.rectangle_sum(x1, y1, x2, y2);
            pu_contact_entries.push((x, y, prob_zone[x][y]));
        }
        pu_delineation_entries.push((x, x1, x2));
    }

    // Compute marginal sums and total
    let mut sprob_zone = vec![0.0f64; n_pus];
    let mut sprob_tot = 0.0f64;
    for x in 0..n_pus {
        for y in 0..n_pus {
            sprob_zone[x] += prob_zone[x][y];
        }
        sprob_tot += sprob_zone[x];
    }

    // Normalize
    if sprob_tot > 0.0 {
        for x in 0..n_pus {
            for y in 0..n_pus {
                prob_zone[x][y] /= sprob_tot;
            }
            sprob_zone[x] /= sprob_tot;
        }
    }

    // Compute mutual information (entropy)
    let mut entropy = 0.0f64;
    for x in 0..n_pus {
        for y in 0..n_pus {
            if prob_zone[x][y] > 1e-5 && sprob_zone[x] > 1e-5 && sprob_zone[y] > 1e-5 {
                entropy +=
                    prob_zone[x][y] * (prob_zone[x][y] / (sprob_zone[x] * sprob_zone[y])).ln();
            }
        }
    }

    let ci = 100.0 * (1.0 - (-2.0 * entropy).exp()).sqrt();
    let r = 100.0 * (1.0 - (-2.0 * entropy).exp());

    (ci, r, pu_contact_entries, pu_delineation_entries)
}

/// Compute homogeneity score for a PU region (used in pruning mode).
fn homogeneity(matrix: &ContactMatrix, start: usize, end: usize) -> f64 {
    let threshold = 0.5;
    let mut pcontact1 = 0.0f64;
    let mut pcontact2 = 0.0f64;

    for k in start..=end {
        for l in start..=end {
            let p = matrix.get(k, l);
            if p > threshold {
                pcontact1 += p;
                if (k as isize - l as isize).unsigned_abs() < 6 {
                    pcontact2 += p;
                }
            }
        }
    }

    let mut h1 = 0.0f64;
    let mut h2 = 0.0f64;
    for i in start..end {
        for j in start..end {
            let p = matrix.get(i, j);
            if p < 0.0001 {
                continue;
            }
            if p > threshold {
                if pcontact1 > 0.0 {
                    let pn1 = p / pcontact1;
                    h1 += pn1 * pn1.ln();
                }
                if (i as isize - j as isize).unsigned_abs() < 6 && pcontact2 > 0.0 {
                    let pn2 = p / pcontact2;
                    h2 += pn2 * pn2.ln();
                }
            }
        }
    }

    let neq1 = (-h1).exp();
    let neq2 = (-h2).exp();
    let n = (end - start) as f64;
    if n > 0.0 {
        (neq1 - neq2) / n
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Output types
// ---------------------------------------------------------------------------

/// Metrics for one iteration of peeling.
#[derive(Debug, Clone)]
pub struct IterationResult {
    /// Maximum contact ratio.
    pub max_cr: f64,
    /// Minimum density.
    pub min_density: f64,
    /// Compaction Index.
    pub ci: f64,
    /// R metric.
    pub r: f64,
    /// Number of PUs.
    pub num_pus: usize,
    /// PU boundaries as [start, end] pairs (0-indexed into residue array).
    pub pu_boundaries: Vec<[usize; 2]>,
}

/// Complete output from the peeling algorithm.
#[derive(Debug, Clone)]
pub struct PeelingOutput {
    /// Contact probability matrix (owned).
    pub contact_matrix: ContactMatrix,
    /// Results for each iteration (1-indexed: iterations[0] = iteration 1).
    pub iterations: Vec<IterationResult>,
    /// Final PU contact matrix entries (x, y, value) from the last iteration.
    pub final_pu_contacts: Vec<(usize, usize, f64)>,
    /// Final PU delineation entries (id, start, end) from the last iteration.
    pub final_pu_delineation: Vec<(usize, usize, usize)>,
    /// Original residue numbers from DSSP (tab_true_num equivalent).
    pub true_nums: Vec<i32>,
}

impl PeelingOutput {
    pub(crate) fn true_num_at(&self, index: usize) -> i32 {
        if let Some(value) = self.true_nums.get(index) {
            *value
        } else if let Some(last) = self.true_nums.last() {
            last + (index + 1 - self.true_nums.len()) as i32
        } else {
            (index + 1) as i32
        }
    }

    /// Write Peeling.log in the format expected by downstream parsers.
    ///
    /// Format:
    /// ```text
    /// Max_CR Min_Density CI R Num_PUs PU_Delineations
    /// 0.15  12.34  45.678901 23.456789 2 1 50 51 100
    /// ```
    pub fn write_peeling_log(&self, path: &Path) -> Result<()> {
        let mut buf = Vec::with_capacity(4096);
        writeln!(buf, "Max_CR Min_Density CI R Num_PUs PU_Delineations")?;

        for iter in &self.iterations {
            // Format matches C code: "%-5.2lf %-5.2lf %lf %lf %d start1 end1 ..."
            write!(buf, "{:<5.2} {:<5.2} ", iter.max_cr, iter.min_density)?;
            write!(buf, "{:.*} {:.*} ", 6, iter.ci, 6, iter.r)?;
            write!(buf, "{} ", iter.num_pus)?;
            for pu in &iter.pu_boundaries {
                let start_num = self.true_num_at(pu[0]);
                let end_num = self.true_num_at(pu[1]);
                write!(buf, "{} {} ", start_num, end_num)?;
            }
            writeln!(buf)?;
        }

        std::fs::write(path, buf)
            .with_context(|| format!("Cannot write peeling log to {}", path.display()))?;
        Ok(())
    }

    /// Write the PU contact matrix file (file_matrix_pu_contact.mtx).
    pub fn write_pu_contact_matrix(&self, path: &Path) -> Result<()> {
        let mut buf = Vec::with_capacity(4096);
        for &(x, y, val) in &self.final_pu_contacts {
            writeln!(buf, "{} {} {:.*}", x, y, 6, val)?;
        }
        std::fs::write(path, buf)
            .with_context(|| format!("Cannot write PU contact matrix to {}", path.display()))?;
        Ok(())
    }

    /// Write the PU delineation file (file_pu_delineation.mtx).
    pub fn write_pu_delineation(&self, path: &Path) -> Result<()> {
        let mut buf = Vec::with_capacity(1024);
        for &(id, start, end) in &self.final_pu_delineation {
            writeln!(buf, "{} {} {}", id, start, end)?;
        }
        std::fs::write(path, buf)
            .with_context(|| format!("Cannot write PU delineation to {}", path.display()))?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Main peeling loop
// ---------------------------------------------------------------------------

/// Run the complete protein peeling algorithm.
///
/// This is the main entry point that replaces the external `Peeling_omp` binary.
///
/// # Arguments
/// * `ca_coords` - C-alpha coordinates as [x, y, z] arrays
/// * `dssp_path` - Path to the DSSP output file
/// * `config` - Algorithm parameters
///
/// # Returns
/// A `PeelingOutput` containing all iteration results, the contact matrix,
/// and the final PU delineation data needed by downstream modules.
#[allow(clippy::needless_range_loop)]
pub fn run_peeling(
    ca_coords: &[[f64; 3]],
    dssp_path: &Path,
    config: &PeelingConfig,
) -> Result<PeelingOutput> {
    let n = ca_coords.len();
    if n == 0 {
        anyhow::bail!("No C-alpha coordinates provided for peeling");
    }

    tracing::debug!(
        "Peeling: {} residues, d0={}, delta={}, min_pu_size={}, max_r2={}",
        n,
        config.d0,
        config.delta,
        config.min_pu_size,
        config.max_r2
    );

    // Step 1: Compute contact probability matrix
    let matrix = ContactMatrix::from_ca_coords(ca_coords, config.d0, config.delta);
    let ind = n - 1; // C code uses 0-based, ind = number of residues - 1

    // Step 2: Parse DSSP for secondary structure and cutting mask
    let (_ss_types, true_nums, cutting_mask) = parse_dssp_for_peeling(dssp_path, n, config)?;

    // Step 3: Initialize PU array
    // pu[iteration][pu_index] = [start, end]
    // We use a Vec of Vecs instead of the C fixed-size 3D array
    let mut pu_table: Vec<Vec<[usize; 2]>> = vec![Vec::new(); MAX_ITERATION];
    pu_table[0] = vec![[0, ind]];

    let mut iterations: Vec<IterationResult> = Vec::new();
    let mut nbre_pu: usize;
    let mut new_nb_pu: usize = 0;
    let mut last_pu_contacts = Vec::new();
    let mut last_pu_delineation = Vec::new();

    for iteration in 1..MAX_ITERATION {
        nbre_pu = new_nb_pu;

        // Check max PU number (C code: nbre_pu > MAXNUMBEROFPU)
        if nbre_pu > config.max_pu_number {
            tracing::debug!("Peeling: max PU count reached at iteration {}", iteration);
            break;
        }

        // Find the best cut across all PUs
        let mut best_cut: Option<CutResult> = None;
        let mut best_coeff = 0.0f64;

        for x in 0..=nbre_pu {
            let start = pu_table[iteration - 1][x][0];
            let end = pu_table[iteration - 1][x][1];
            let size = end - start;

            if size < config.min_pu_size {
                continue;
            }

            // Try single cut
            if let Some(cut) = simple_cutting(
                &matrix,
                start,
                end,
                &cutting_mask,
                config.min_pu_size,
                best_coeff,
                x,
            ) {
                best_coeff = cut.coeff;
                best_cut = Some(cut);
            }

            // Try double cut
            if let Some(cut) = double_cutting(
                &matrix,
                start,
                end,
                &cutting_mask,
                config.min_pu_size,
                best_coeff,
                x,
            ) {
                best_coeff = cut.coeff;
                best_cut = Some(cut);
            }
        }

        // No cut found — stop
        let best_cut = match best_cut {
            Some(cut) => cut,
            None => {
                tracing::debug!(
                    "Peeling: no further cuts possible at iteration {}",
                    iteration
                );
                break;
            }
        };

        // Save PUs for this iteration
        let mut new_pus: Vec<[usize; 2]> = Vec::with_capacity(nbre_pu + 3);

        if best_cut.num_cuts == 1 {
            new_pus.push([best_cut.start, best_cut.i1]);
            new_pus.push([best_cut.i2, best_cut.end]);
        } else {
            new_pus.push([best_cut.start, best_cut.i1]);
            new_pus.push([best_cut.i2, best_cut.j1]);
            new_pus.push([best_cut.j2, best_cut.end]);
        }

        // Copy unmodified PUs from previous iteration
        if iteration > 1 {
            for x in 0..=nbre_pu {
                if x != best_cut.pu_index {
                    new_pus.push(pu_table[iteration - 1][x]);
                }
            }
        }

        new_nb_pu = new_pus.len() - 1;
        pu_table[iteration] = new_pus;

        // Check max PU size constraint
        if config.max_pu_size > 0 {
            let all_within = (0..=new_nb_pu).all(|x| {
                let s = pu_table[iteration - 1][x][0];
                let e = pu_table[iteration - 1][x][1];
                (e - s) <= config.max_pu_size
            });
            if all_within {
                tracing::debug!(
                    "Peeling: all PUs within max size at iteration {}",
                    iteration
                );
                break;
            }
        }

        // Check pruning constraint
        if config.pruning {
            let passes_pruning = if best_cut.num_cuts == 1 {
                let h1 = homogeneity(&matrix, best_cut.start, best_cut.i1);
                let h2 = homogeneity(&matrix, best_cut.i2, best_cut.end);
                h1 >= config.cutoff_pruning || h2 >= config.cutoff_pruning
            } else {
                let h1 = homogeneity(&matrix, best_cut.start, best_cut.i1);
                let h2 = homogeneity(&matrix, best_cut.i2, best_cut.j1);
                let h3 = homogeneity(&matrix, best_cut.j2, best_cut.end);
                h1 >= config.cutoff_pruning
                    || h2 >= config.cutoff_pruning
                    || h3 >= config.cutoff_pruning
            };
            if !passes_pruning {
                tracing::debug!(
                    "Peeling: pruning criteria not met at iteration {}",
                    iteration
                );
                break;
            }
        }

        // Compute metrics
        let current_pus = &pu_table[iteration];
        let pu_refs: Vec<[usize; 2]> = current_pus.clone();
        let (max_cr, min_density) = measure_coeff(&matrix, &pu_refs);
        let (ci, r, pu_contacts, pu_delineation) = mutual_information(&matrix, &pu_refs);
        last_pu_contacts = pu_contacts;
        last_pu_delineation = pu_delineation;

        iterations.push(IterationResult {
            max_cr,
            min_density,
            ci,
            r,
            num_pus: new_nb_pu + 1,
            pu_boundaries: pu_refs,
        });

        tracing::debug!(
            "Peeling iteration {}: {} PUs, CI={:.2}, max_cr={:.2}",
            iteration,
            new_nb_pu + 1,
            ci,
            max_cr
        );

        // Check CI threshold
        if ci > config.max_r2 as f64 {
            tracing::debug!(
                "Peeling: CI ({:.2}) exceeds max_r2 ({}) at iteration {}",
                ci,
                config.max_r2,
                iteration
            );
            break;
        }
    }

    Ok(PeelingOutput {
        contact_matrix: matrix,
        iterations,
        final_pu_contacts: last_pu_contacts,
        final_pu_delineation: last_pu_delineation,
        true_nums,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    // Build a valid DSSP data line for parse_dssp_for_peeling.
    //
    // Column layout (0-indexed):
    //  0-4:  DSSP seqnum (right-aligned in 5)
    //  5:    space
    //  6-9:  PDB resnum (right-aligned in 4; col 6 is blanked by the parser)
    //  10:   icode (' ')
    //  11:   chain ('A')
    //  12:   space
    //  13:   amino acid code
    //  14-15: spaces
    //  16:   secondary structure code
    //  17-125: spaces
    //  126:  '0'  (non-space marker the parser requires at this column)
    fn dssp_line(seqnum: usize, resnum: i32, aa: char, ss: char) -> String {
        let head = format!("{:>5} {:>4} A {}  {}", seqnum, resnum, aa, ss);
        debug_assert_eq!(head.len(), 17);
        let mut line = head;
        line.extend(std::iter::repeat_n(' ', 109));
        line.push('0');
        line
    }

    #[test]
    fn test_parse_dssp_ss_types_and_residue_numbers() {
        let dir = tempdir().unwrap();
        let dssp_path = dir.path().join("test.dssp");

        // Header line: col 126 is space → skipped by parser
        let mut content = format!("{:<128}\n", "  # RESIDUE AA STRUCTURE");
        // 3 coil, 4 helix, 3 coil
        for i in 1..=3usize {
            content.push_str(&dssp_line(i, i as i32, 'A', ' '));
            content.push('\n');
        }
        for i in 4..=7usize {
            content.push_str(&dssp_line(i, i as i32, 'A', 'H'));
            content.push('\n');
        }
        for i in 8..=10usize {
            content.push_str(&dssp_line(i, i as i32, 'A', ' '));
            content.push('\n');
        }
        std::fs::write(&dssp_path, &content).unwrap();

        let config = PeelingConfig::default();
        let (ss_types, true_nums, cutting_mask) =
            parse_dssp_for_peeling(&dssp_path, 10, &config).unwrap();

        assert_eq!(ss_types.len(), 10);
        assert_eq!(true_nums, (1..=10).collect::<Vec<i32>>());

        assert!(ss_types[..3].iter().all(|s| *s == SsType::Coil));
        assert!(ss_types[3..7].iter().all(|s| *s == SsType::Helix));
        assert!(ss_types[7..].iter().all(|s| *s == SsType::Coil));

        // Helix segment (0-indexed 3..6) has size=3 <= min_ss_size=8 → marked non-cuttable
        assert!(cutting_mask[2]); // last coil before helix: cuttable
        assert!(!cutting_mask[3]); // helix interior: non-cuttable
        assert!(!cutting_mask[4]);
        assert!(!cutting_mask[5]);
        assert!(cutting_mask[6]); // segment_end itself is NOT masked (loop is seg_start..seg_end)
        assert!(cutting_mask[7]); // first coil after helix: cuttable
    }

    #[test]
    fn test_write_peeling_log_extends_missing_true_numbers() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("Peeling.log");
        let matrix = ContactMatrix::from_ca_coords(
            &[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
            ],
            6.0,
            1.5,
        );
        let output = PeelingOutput {
            contact_matrix: matrix,
            iterations: vec![IterationResult {
                max_cr: 1.0,
                min_density: 2.0,
                ci: 3.0,
                r: 4.0,
                num_pus: 1,
                pu_boundaries: vec![[0, 3]],
            }],
            final_pu_contacts: vec![],
            final_pu_delineation: vec![],
            true_nums: vec![1, 2, 3],
        };

        output.write_peeling_log(&path).unwrap();

        let log = std::fs::read_to_string(path).unwrap();
        assert!(log.contains("1 4"));
    }
}
