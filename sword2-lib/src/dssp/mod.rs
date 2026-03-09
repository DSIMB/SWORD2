//! DSSP secondary structure assignment.
//!
//! This module will provide functionality for parsing DSSP output
//! and assigning secondary structure to residues.

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
