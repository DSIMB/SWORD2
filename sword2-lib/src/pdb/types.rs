// Core PDB data types for representing protein structures.
// Ported from Python PDB.py and Perl ParsePDB.pm.

use std::fmt;

use super::amino_acids;

// ---------------------------------------------------------------------------
// Point3D
// ---------------------------------------------------------------------------

/// A 3D coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3D {
    /// Create a new 3D point.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Euclidean distance to another point.
    pub fn distance_to(&self, other: &Point3D) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

impl fmt::Display for Point3D {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({:.3}, {:.3}, {:.3})", self.x, self.y, self.z)
    }
}

impl std::ops::Sub for Point3D {
    type Output = Point3D;
    fn sub(self, rhs: Point3D) -> Point3D {
        Point3D::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl std::ops::Add for Point3D {
    type Output = Point3D;
    fn add(self, rhs: Point3D) -> Point3D {
        Point3D::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

// ---------------------------------------------------------------------------
// Atom
// ---------------------------------------------------------------------------

/// Represents a single ATOM / HETATM record in PDB format.
#[derive(Debug, Clone)]
pub struct Atom {
    pub serial: i32,
    pub name: String,
    pub alt_loc: char,
    pub res_name: String,
    pub chain_id: char,
    pub res_seq: i32,
    pub icode: char,
    pub coord: Point3D,
    pub occupancy: f64,
    pub temp_factor: f64,
    pub element: String,
    pub charge: String,
    pub is_hetatm: bool,
}

impl Atom {
    /// Create a new Atom from all fields.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        serial: i32,
        name: &str,
        alt_loc: char,
        res_name: &str,
        chain_id: char,
        res_seq: i32,
        icode: char,
        coord: Point3D,
        occupancy: f64,
        temp_factor: f64,
        element: &str,
        charge: &str,
        is_hetatm: bool,
    ) -> Self {
        Self {
            serial,
            name: name.to_string(),
            alt_loc,
            res_name: res_name.to_string(),
            chain_id,
            res_seq,
            icode,
            coord,
            occupancy,
            temp_factor,
            element: element.to_string(),
            charge: charge.to_string(),
            is_hetatm,
        }
    }

    /// Returns `true` if this is a C-alpha atom.
    pub fn is_ca(&self) -> bool {
        self.name.trim() == "CA"
    }

    /// Returns `true` if this is a backbone atom (N, CA, C, O).
    pub fn is_backbone(&self) -> bool {
        matches!(self.name.trim(), "N" | "CA" | "C" | "O")
    }
}

impl fmt::Display for Atom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {} {} {}{}",
            if self.is_hetatm { "HETATM" } else { "ATOM" },
            self.serial,
            self.name,
            self.res_name,
            self.res_seq,
        )
    }
}

// ---------------------------------------------------------------------------
// Residue
// ---------------------------------------------------------------------------

/// A residue is a group of atoms sharing the same (chain_id, res_seq, icode).
#[derive(Debug, Clone)]
pub struct Residue {
    pub name: String,
    pub seq_num: i32,
    pub icode: char,
    pub chain_id: char,
    pub atoms: Vec<Atom>,
}

impl Residue {
    /// Create a new residue.
    pub fn new(name: &str, seq_num: i32, icode: char, chain_id: char) -> Self {
        Self {
            name: name.to_string(),
            seq_num,
            icode,
            chain_id,
            atoms: Vec::new(),
        }
    }

    /// Find the C-alpha atom in this residue, if present.
    pub fn get_ca(&self) -> Option<&Atom> {
        self.atoms.iter().find(|a| a.is_ca())
    }

    /// The residue number (seq_num).
    pub fn residue_number(&self) -> i32 {
        self.seq_num
    }

    /// One-letter amino acid code for this residue, or 'X' if unknown.
    pub fn one_letter_code(&self) -> char {
        amino_acids::three_to_one(&self.name).unwrap_or('X')
    }
}

impl fmt::Display for Residue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}", self.name, self.seq_num, self.icode)
    }
}

// ---------------------------------------------------------------------------
// Chain
// ---------------------------------------------------------------------------

/// A chain is a sequence of residues identified by a single-character chain ID.
#[derive(Debug, Clone)]
pub struct Chain {
    pub id: char,
    pub residues: Vec<Residue>,
}

impl Chain {
    /// Create a new empty chain.
    pub fn new(id: char) -> Self {
        Self {
            id,
            residues: Vec::new(),
        }
    }

    /// Number of residues.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// Whether the chain has no residues.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// Get the amino acid sequence as a String of one-letter codes.
    pub fn get_sequence(&self) -> String {
        self.residues.iter().map(|r| r.one_letter_code()).collect()
    }

    /// Find a residue by its sequence number. Returns the first match.
    pub fn get_residue_by_num(&self, seq_num: i32) -> Option<&Residue> {
        self.residues.iter().find(|r| r.seq_num == seq_num)
    }

    /// Find a residue by its sequence number (mutable).
    pub fn get_residue_by_num_mut(&mut self, seq_num: i32) -> Option<&mut Residue> {
        self.residues.iter_mut().find(|r| r.seq_num == seq_num)
    }

    /// Collect all C-alpha atoms.
    pub fn ca_atoms(&self) -> Vec<&Atom> {
        self.residues.iter().filter_map(|r| r.get_ca()).collect()
    }
}

impl fmt::Display for Chain {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Chain {} ({} residues)", self.id, self.len())
    }
}

// ---------------------------------------------------------------------------
// Model
// ---------------------------------------------------------------------------

/// A model contains chains (NMR structures can have multiple models).
#[derive(Debug, Clone)]
pub struct Model {
    pub serial: i32,
    pub chains: Vec<Chain>,
}

impl Model {
    /// Create a new empty model.
    pub fn new(serial: i32) -> Self {
        Self {
            serial,
            chains: Vec::new(),
        }
    }

    /// Get a chain by its ID.
    pub fn get_chain(&self, chain_id: char) -> Option<&Chain> {
        self.chains.iter().find(|c| c.id == chain_id)
    }

    /// Get a mutable reference to a chain by its ID.
    pub fn get_chain_mut(&mut self, chain_id: char) -> Option<&mut Chain> {
        self.chains.iter_mut().find(|c| c.id == chain_id)
    }

    /// List all chain IDs present in this model.
    pub fn chain_ids(&self) -> Vec<char> {
        self.chains.iter().map(|c| c.id).collect()
    }
}

impl fmt::Display for Model {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Model {} ({} chains)",
            self.serial,
            self.chains.len()
        )
    }
}

// ---------------------------------------------------------------------------
// Structure
// ---------------------------------------------------------------------------

/// Top-level structure parsed from a PDB or mmCIF file.
#[derive(Debug, Clone)]
pub struct Structure {
    pub name: String,
    pub models: Vec<Model>,
    pub header_lines: Vec<String>,
}

impl Structure {
    /// Create a new empty structure.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            models: Vec::new(),
            header_lines: Vec::new(),
        }
    }

    /// Get a model by its serial number.
    pub fn get_model(&self, serial: i32) -> Option<&Model> {
        self.models.iter().find(|m| m.serial == serial)
    }

    /// Get the first model (most common use case).
    pub fn first_model(&self) -> Option<&Model> {
        self.models.first()
    }

    /// Get a mutable reference to the first model.
    pub fn first_model_mut(&mut self) -> Option<&mut Model> {
        self.models.first_mut()
    }

    /// List all chain IDs across all models (typically from first model).
    pub fn chain_ids(&self) -> Vec<char> {
        match self.first_model() {
            Some(model) => model.chain_ids(),
            None => Vec::new(),
        }
    }

    /// Total number of atoms across all models.
    pub fn atom_count(&self) -> usize {
        self.models
            .iter()
            .flat_map(|m| &m.chains)
            .flat_map(|c| &c.residues)
            .map(|r| r.atoms.len())
            .sum()
    }
}

impl fmt::Display for Structure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Structure '{}' ({} models, {} atoms)",
            self.name,
            self.models.len(),
            self.atom_count()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point3d_distance() {
        let a = Point3D::new(0.0, 0.0, 0.0);
        let b = Point3D::new(3.0, 4.0, 0.0);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_atom_is_ca() {
        let atom = Atom::new(
            1, " CA ", ' ', "ALA", 'A', 1, ' ',
            Point3D::new(0.0, 0.0, 0.0),
            1.0, 0.0, "C", "", false,
        );
        assert!(atom.is_ca());
    }

    #[test]
    fn test_chain_sequence() {
        let mut chain = Chain::new('A');
        let mut r1 = Residue::new("ALA", 1, ' ', 'A');
        let mut r2 = Residue::new("GLY", 2, ' ', 'A');
        // Add a dummy CA to each residue so they are valid
        r1.atoms.push(Atom::new(1, " CA ", ' ', "ALA", 'A', 1, ' ', Point3D::new(0.0,0.0,0.0), 1.0, 0.0, "C", "", false));
        r2.atoms.push(Atom::new(2, " CA ", ' ', "GLY", 'A', 2, ' ', Point3D::new(1.0,0.0,0.0), 1.0, 0.0, "C", "", false));
        chain.residues.push(r1);
        chain.residues.push(r2);
        assert_eq!(chain.get_sequence(), "AG");
        assert_eq!(chain.len(), 2);
    }
}
