//! PDB file writer.
//!
//! Writes a cleaned protein structure to PDB format,
//! matching the output format of ProDy's writePDB function.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use super::types::Chain;

/// Clean a chain for SWORD processing:
/// - Keep only ATOM records (not HETATM)
/// - Keep only standard amino acids (the 20 classical residues)
/// - Remove residues with insertion codes
/// - Renumber residues sequentially from 1
///
/// Returns (cleaned_chain, original_residue_numbers) where original_residue_numbers
/// maps new index to original residue number.
pub fn clean_chain_for_sword(chain: &Chain) -> (Chain, Vec<i32>) {
    let mut cleaned = Chain::new(chain.id);
    let mut original_resnums = Vec::new();

    for residue in &chain.residues {
        if !super::structural_quality::is_sword_candidate_residue(residue) {
            continue;
        }

        original_resnums.push(residue.seq_num);
    }

    // Now create the cleaned chain with renumbered residues
    let mut new_resnum = 1i32;
    for residue in &chain.residues {
        if !super::structural_quality::is_sword_candidate_residue(residue) {
            continue;
        }

        let mut new_residue = super::types::Residue::new(&residue.name, new_resnum, ' ', chain.id);

        // Copy non-HETATM atoms with updated residue numbers
        for atom in &residue.atoms {
            if atom.is_hetatm {
                continue;
            }
            let mut new_atom = atom.clone();
            new_atom.res_seq = new_resnum;
            new_atom.icode = ' ';
            new_residue.atoms.push(new_atom);
        }

        cleaned.residues.push(new_residue);
        new_resnum += 1;
    }

    (cleaned, original_resnums)
}

/// Write a chain as a PDB file in standard PDB format.
///
/// This matches the output format of ProDy's writePDB so SWORD can read it.
/// PDB ATOM record format (1-based columns):
///  1- 6  Record name ("ATOM  " or "HETATM")
///  7-11  Atom serial number (right-justified)
/// 12     Space
/// 13-16  Atom name (4 chars, preserved as-is from parsing)
/// 17     Alternate location indicator
/// 18-20  Residue name (right-justified in 3 chars)
/// 21     Space
/// 22     Chain ID
/// 23-26  Residue sequence number (right-justified)
/// 27     Code for insertion of residues
/// 28-30  Spaces
/// 31-38  X coordinate (8.3)
/// 39-46  Y coordinate (8.3)
/// 47-54  Z coordinate (8.3)
/// 55-60  Occupancy (6.2)
/// 61-66  Temperature factor (6.2)
/// 67-76  Spaces
/// 77-78  Element symbol (right-justified)
/// 79-80  Charge
pub fn write_pdb(chain: &Chain, output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Cannot create PDB file: {}", output_path.display()))?;

    let mut serial = 1i32;
    for residue in &chain.residues {
        for atom in &residue.atoms {
            let record = if atom.is_hetatm { "HETATM" } else { "ATOM  " };

            // Atom name: must be exactly 4 chars following PDB convention.
            // PDB format: 1-char element names (N, C, O, S) start at col 14,
            // so they are right-justified: " N  ", " C  ", " O  ", " S  "
            // 2-char names start at col 13: " CA ", " CB ", " CG ", " NZ "
            // 3-char names: " CG1", " OD1", " NE2"
            // 4-char names: "1HG2", "2HD2" (hydrogen naming)
            let trimmed = atom.name.trim();
            let atom_name = if trimmed.len() == 4 {
                trimmed.to_string()
            } else if trimmed.len() == 3 {
                format!(" {}", trimmed)
            } else if trimmed.len() == 2 {
                format!(" {} ", trimmed)
            } else if trimmed.len() == 1 {
                format!(" {}  ", trimmed)
            } else {
                format!("{:<4}", trimmed)
            };

            let alt_loc = atom.alt_loc;

            // Build the line using exact column positions
            let line = format!(
                "{}{:>5} {}{}{:>3} {}{:>4}{}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {:>2}{:>2}",
                record,          // 1-6
                serial,          // 7-11
                atom_name,       // 13-16 (space at col 12 from format)
                alt_loc,         // 17
                residue.name,    // 18-20
                chain.id,        // 22 (space at 21 from format)
                residue.seq_num, // 23-26
                ' ',             // 27 (icode)
                atom.coord.x,
                atom.coord.y,
                atom.coord.z,
                atom.occupancy,
                atom.temp_factor,
                atom.element,
                atom.charge,
            );

            writeln!(f, "{}", line)?;
            serial += 1;
        }
    }
    writeln!(f, "END")?;

    Ok(())
}

/// Remove residues whose mean B-factor (pLDDT in AlphaFold/ESM) is below `min_plddt`.
///
/// Returns a new chain without the low-confidence residues. Renumbering is left
/// to `clean_chain_for_sword`, which always runs immediately after.
pub fn filter_by_plddt(chain: &Chain, min_plddt: f64) -> Chain {
    let mut filtered = Chain::new(chain.id);
    for residue in &chain.residues {
        if residue.atoms.is_empty() {
            continue;
        }
        let mean_b: f64 = residue.atoms.iter().map(|a| a.temp_factor).sum::<f64>()
            / residue.atoms.len() as f64;
        if mean_b >= min_plddt {
            filtered.residues.push(residue.clone());
        }
    }
    filtered
}

/// Write one PDB file per domain in a partition.
///
/// Each file contains only the residues whose seq_num falls within any
/// (start, end) segment of that domain. Files are written to `output_dir`
/// as `domain_1.pdb`, `domain_2.pdb`, etc.
pub fn write_domain_pdbs(
    chain: &Chain,
    boundaries: &[Vec<(i32, i32)>],
    output_dir: &Path,
) -> Result<()> {
    std::fs::create_dir_all(output_dir)?;

    for (domain_idx, segments) in boundaries.iter().enumerate() {
        let mut domain_chain = Chain::new(chain.id);
        for residue in &chain.residues {
            let in_domain = segments
                .iter()
                .any(|&(start, end)| residue.seq_num >= start && residue.seq_num <= end);
            if in_domain {
                domain_chain.residues.push(residue.clone());
            }
        }
        let out_path = output_dir.join(format!("domain_{}.pdb", domain_idx + 1));
        write_pdb(&domain_chain, &out_path)?;
    }
    Ok(())
}

/// Write the residue number mapping file.
///
/// Format matches Python: "# Mapping of authors PDB residues numbers\n# with the new one...\nORIGINAL RENUM\n{orig} {new}\n..."
pub fn write_mapping_file(original_resnums: &[i32], output_path: &Path) -> Result<()> {
    let mut f = std::fs::File::create(output_path)
        .with_context(|| format!("Cannot create mapping file: {}", output_path.display()))?;

    writeln!(f, "# Mapping of authors PDB residues numbers ")?;
    writeln!(f, "# with the new one presented on the server")?;
    writeln!(f, "ORIGINAL RENUM")?;
    for (i, &orig) in original_resnums.iter().enumerate() {
        writeln!(f, "{} {}", orig, i + 1)?;
    }

    Ok(())
}
