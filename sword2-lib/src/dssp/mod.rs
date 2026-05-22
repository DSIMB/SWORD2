//! DSSP secondary structure assignment — pure Rust implementation.
//!
//! This module provides a complete DSSP (Define Secondary Structure of Proteins)
//! implementation ported from the classic Kabsch & Sander algorithm (1983).
//!
//! Key optimization: spatial grid for H-bond detection (O(N·k) vs original O(N²)).
//!
//! Reference: Kabsch, W. and Sander, C. (1983) Biopolymers 22, 2577-2637.

pub mod angles;
pub mod backbone;
pub mod bridge;
pub mod format;
pub mod hbond;
pub mod helix;
pub mod types;

use std::path::Path;

use anyhow::Result;

pub use types::DsspChain;

/// Secondary structure types from DSSP.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SecondaryStructure {
    Helix,
    Sheet,
    Turn,
    Coil,
}

impl SecondaryStructure {
    /// Parse a single DSSP character code.
    pub fn from_dssp_char(c: char) -> Self {
        match c {
            'H' | 'G' | 'I' => SecondaryStructure::Helix,
            'E' | 'B' => SecondaryStructure::Sheet,
            'T' | 'S' => SecondaryStructure::Turn,
            _ => SecondaryStructure::Coil,
        }
    }

    /// Convert to single-character code.
    pub fn to_char(self) -> char {
        match self {
            SecondaryStructure::Helix => 'H',
            SecondaryStructure::Sheet => 'E',
            SecondaryStructure::Turn => 'T',
            SecondaryStructure::Coil => 'C',
        }
    }
}

/// Result of running DSSP on a structure.
pub struct DsspResult {
    /// The DSSP chain with all assignments.
    pub chain: DsspChain,
    /// Path to the written .dssp file.
    pub dssp_path: std::path::PathBuf,
    /// Path to the written .s2d file.
    pub s2d_path: std::path::PathBuf,
}

/// Run the full DSSP algorithm on a cleaned PDB file.
///
/// This replaces the external `dsspcmbi` binary. Reads the PDB file, extracts
/// backbone atoms, computes H-bonds, assigns secondary structure, and writes
/// both a .dssp file (for Peeling) and a .s2d file (for SWORD pipeline).
///
/// # Arguments
/// * `pdb_path` - Path to the cleaned PDB file
/// * `dssp_path` - Where to write the .dssp output
/// * `s2d_path` - Where to write the .s2d output
/// * `pdb_name` - Name used in DSSP header
pub fn run_dssp(
    pdb_path: &Path,
    dssp_path: &Path,
    s2d_path: &Path,
    pdb_name: &str,
) -> Result<DsspResult> {
    // Step 1: Extract backbone atoms and synthesize H positions
    let mut chain = backbone::extract_backbone(pdb_path)?;

    tracing::debug!(
        "DSSP: {} residues extracted (including {} chain breaks)",
        chain.len,
        (1..=chain.len).filter(|&i| chain.get(i).aa == '!').count()
    );

    // Step 2: Calculate backbone angles (kappa, alpha, chirality)
    angles::calculate_angles(&mut chain);

    // Step 3: Detect H-bonds (spatial grid optimization)
    hbond::flag_hydrogen_bonds(&mut chain);

    // Step 4: Detect beta-bridges, build ladders, assemble sheets
    bridge::flag_bridges(&mut chain);

    // Step 5: Detect turns, assign helix/bend/SS symbols
    helix::flag_turns(&mut chain);

    // Step 6: Write output files
    format::write_dssp(&chain, pdb_name, dssp_path)?;
    format::write_s2d(&chain, pdb_name, s2d_path)?;

    tracing::debug!(
        "DSSP: wrote {} and {}",
        dssp_path.display(),
        s2d_path.display()
    );

    Ok(DsspResult {
        chain,
        dssp_path: dssp_path.to_path_buf(),
        s2d_path: s2d_path.to_path_buf(),
    })
}
