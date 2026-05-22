//! Core data structures for DSSP secondary structure assignment.

use crate::pdb::types::Point3D;

/// Physical constants from the original DSSP algorithm.
pub const BREAKDIST: f64 = 2.5; // Max peptide bond C-N distance (Angstrom)
pub const CADIST: f64 = 9.0; // CA distance cutoff for H-bond candidates
pub const Q: f64 = -27888.0; // Electrostatic coupling constant (cal/mol)
pub const HBLOW: i64 = -9900; // Min H-bond energy (cal/mol)
pub const HBHIGH: i64 = -500; // Max H-bond energy (cal/mol)
pub const DIST_MIN: f64 = 0.5; // Smallest allowed distance between atoms
pub const RADIAN: f64 = 57.29578; // 180/PI

/// Maximum number of bridges.
pub const MAXBRIDGE: usize = 2000;

/// A backbone residue with coordinates and DSSP state.
#[derive(Debug, Clone)]
pub struct BackboneResidue {
    /// One-letter amino acid code ('!' for chain break).
    pub aa: char,
    /// PDB residue identifier string (6 chars: resnum + icode + chain).
    pub aaident: [u8; 6],
    /// Three-letter residue code.
    pub three_letter: [u8; 4],
    /// Backbone atom coordinates.
    pub n: Point3D,
    pub ca: Point3D,
    pub c: Point3D,
    pub o: Point3D,
    /// Synthesized backbone H position.
    pub h: Point3D,
    /// Whether this residue has a valid amide H (false for Pro, chain start).
    pub has_h: bool,
    /// Secondary structure columns:
    /// [0]=symbol, [1]=turn3, [2]=turn4, [3]=turn5, [4]=bend, [5]=chirality, [6]=beta1, [7]=beta2
    pub ss: [char; 8],
    /// Bridge partners for beta1 and beta2.
    pub partner: [usize; 2],
    /// Sheet label character.
    pub sheet_label: char,
    /// Best 2 acceptor H-bonds (this residue's NH donates to acceptor's CO).
    pub acceptor: [HydrogenBond; 2],
    /// Best 2 donor H-bonds (donor's NH donates to this residue's CO).
    pub donor: [HydrogenBond; 2],
    /// Solvent accessibility (not computed when -na flag is used).
    pub access: i64,
    /// Virtual bend angle kappa.
    pub kappa: f64,
    /// Virtual dihedral alpha.
    pub alpha: f64,
}

impl BackboneResidue {
    /// Create a chain break marker residue.
    pub fn chain_break() -> Self {
        Self {
            aa: '!',
            ..Self::default()
        }
    }
}

impl Default for BackboneResidue {
    fn default() -> Self {
        Self {
            aa: ' ',
            aaident: [b' '; 6],
            three_letter: [b' '; 4],
            n: Point3D::new(0.0, 0.0, 0.0),
            ca: Point3D::new(0.0, 0.0, 0.0),
            c: Point3D::new(0.0, 0.0, 0.0),
            o: Point3D::new(0.0, 0.0, 0.0),
            h: Point3D::new(0.0, 0.0, 0.0),
            has_h: false,
            ss: [' '; 8],
            partner: [0; 2],
            sheet_label: ' ',
            acceptor: [HydrogenBond::default(); 2],
            donor: [HydrogenBond::default(); 2],
            access: 0,
            kappa: 360.0,
            alpha: 360.0,
        }
    }
}

/// A hydrogen bond with residue index and energy.
#[derive(Debug, Clone, Copy, Default)]
pub struct HydrogenBond {
    /// 1-based residue index (0 = no bond).
    pub residue: usize,
    /// Energy in cal/mol.
    pub energy: i64,
}

/// Bridge type between two residues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BridgeType {
    Parallel,
    Antiparallel,
    NoBridge,
}

/// A beta-bridge/ladder entry in the bridge table.
#[derive(Debug, Clone)]
pub struct Bridge {
    pub sheet_name: char,
    pub ladder_name: char,
    pub btype: BridgeType,
    /// Set of linked ladder indices.
    pub link_set: Vec<bool>,
    /// i-strand begin/end (1-based).
    pub ib: usize,
    pub ie: usize,
    /// j-strand begin/end (1-based).
    pub jb: usize,
    pub je: usize,
    /// Link to previous/next ladder in the chain.
    pub from: usize,
    pub towards: usize,
}

impl Bridge {
    pub fn new(max_bridges: usize) -> Self {
        Self {
            sheet_name: ' ',
            ladder_name: ' ',
            btype: BridgeType::NoBridge,
            link_set: vec![false; max_bridges + 1],
            ib: 0,
            ie: 0,
            jb: 0,
            je: 0,
            from: 0,
            towards: 0,
        }
    }
}

/// The DSSP chain: 1-indexed array of backbone residues.
pub struct DsspChain {
    /// Residues, 1-indexed. Index 0 is unused (sentinel).
    pub residues: Vec<BackboneResidue>,
    /// Number of actual residues (length of chain, excluding index 0).
    pub len: usize,
}

impl Default for DsspChain {
    fn default() -> Self {
        Self::new()
    }
}

impl DsspChain {
    pub fn new() -> Self {
        Self {
            residues: vec![BackboneResidue::default()], // index 0 sentinel
            len: 0,
        }
    }

    /// Get residue at 1-based index.
    pub fn get(&self, i: usize) -> &BackboneResidue {
        &self.residues[i]
    }

    /// Get mutable residue at 1-based index.
    pub fn get_mut(&mut self, i: usize) -> &mut BackboneResidue {
        &mut self.residues[i]
    }

    /// Push a residue onto the chain. Returns its 1-based index.
    pub fn push(&mut self, res: BackboneResidue) -> usize {
        self.residues.push(res);
        self.len = self.residues.len() - 1;
        self.len
    }

    /// Check if there's no chain break between residues i and j (inclusive, 1-based).
    pub fn no_chain_break(&self, i: usize, j: usize) -> bool {
        if i < 1 || j > self.len || i > j {
            return false;
        }
        for k in i..=j {
            if self.residues[k].aa == '!' {
                return false;
            }
        }
        true
    }
}
