//! PDB file parsing and protein structure representation.

pub mod amino_acids;
pub mod parser;
pub mod types;
pub mod writer;

pub use amino_acids::three_to_one;
pub use parser::parse_pdb;
pub use types::{Atom, Chain, Model, Point3D, Residue, Structure};
pub use writer::write_pdb;
