//! PDB file parser.
//!
//! Parses ATOM/HETATM records from PDB-format files into a [`Structure`].
//! Ported from the Python PDB.py parser.

use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::GzDecoder;
use pdbtbx::StrictnessLevel;

use super::types::{Atom, Chain, Model, Point3D, Residue, Structure};

/// Parse a protein structure file from the given path.
///
/// Supports plain `.pdb`, `.cif` (mmCIF) files and gzip-compressed versions.
pub fn parse_pdb(path: &Path) -> Result<Structure> {
    let name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");
    
    // Check if it's mmCIF (possibly gzipped)
    if ext == "cif" || (ext == "gz" && path.to_str().unwrap_or("").ends_with(".cif.gz")) {
        return parse_mmcif(path);
    }

    let file = fs::File::open(path).with_context(|| format!("Cannot open {}", path.display()))?;

    let reader: Box<dyn BufRead> = if ext == "gz" {
        Box::new(BufReader::new(GzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    parse_pdb_reader(reader, &name)
}

/// Parse an mmCIF file using pdbtbx and convert it to our internal Structure type.
pub fn parse_mmcif(path: &Path) -> Result<Structure> {
    let path_str = path.to_str().ok_or_else(|| anyhow::anyhow!("Invalid path"))?;
    let (pdbtbx_struct, _warnings) = pdbtbx::ReadOptions::new()
        .set_level(StrictnessLevel::Loose)
        .read(path_str)
        .map_err(|e| anyhow::anyhow!("pdbtbx error: {:?}", e))?;

    let name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string();

    let mut structure = Structure::new(&name);

    for pdbtbx_model in pdbtbx_struct.models() {
        let mut model = Model::new(pdbtbx_model.serial_number() as i32);
        for pdbtbx_chain in pdbtbx_model.chains() {
            let chain_id = pdbtbx_chain.id().chars().next().unwrap_or(' ');
            let mut chain = Chain::new(chain_id);

            for pdbtbx_residue in pdbtbx_chain.residues() {
                let res_name = pdbtbx_residue.name().unwrap_or("UNK");
                let res_serial = pdbtbx_residue.serial_number() as i32;
                let icode = pdbtbx_residue.insertion_code().and_then(|s| s.chars().next()).unwrap_or(' ');
                
                let mut residue = Residue::new(
                    res_name,
                    res_serial,
                    icode,
                    chain_id,
                );

                for pdbtbx_conformer in pdbtbx_residue.conformers() {
                    let alt_loc = pdbtbx_conformer.alternative_location()
                        .and_then(|s| s.chars().next())
                        .unwrap_or(' ');

                    for pdbtbx_atom in pdbtbx_conformer.atoms() {
                        let atom = Atom::new(
                            pdbtbx_atom.serial_number() as i32,
                            pdbtbx_atom.name(),
                            alt_loc,
                            res_name,
                            chain_id,
                            res_serial,
                            icode,
                            Point3D::new(pdbtbx_atom.x(), pdbtbx_atom.y(), pdbtbx_atom.z()),
                            pdbtbx_atom.occupancy(),
                            pdbtbx_atom.b_factor(),
                            pdbtbx_atom.element().map(|e| e.to_string()).unwrap_or_default().as_str(),
                            "", // Charge
                            pdbtbx_atom.hetero(),
                        );
                        residue.atoms.push(atom);
                    }
                }
                chain.residues.push(residue);
            }
            model.chains.push(chain);
        }
        structure.models.push(model);
    }

    Ok(structure)
}

/// Parse PDB content from a string.
pub fn parse_pdb_str(content: &str, name: &str) -> Result<Structure> {
    let reader = BufReader::new(content.as_bytes());
    parse_pdb_reader(reader, name)
}

/// Core parser that reads PDB records from any `BufRead` source.
fn parse_pdb_reader<R: BufRead>(reader: R, name: &str) -> Result<Structure> {
    let mut structure = Structure::new(name);
    let mut current_model = Model::new(1);
    let mut in_model = false;
    let mut model_count = 0;

    for line_result in reader.lines() {
        let line = line_result?;
        let record_type = if line.len() >= 6 {
            line[..6].trim()
        } else {
            line.trim()
        };

        match record_type {
            "MODEL" => {
                in_model = true;
                let serial = line
                    .get(10..14)
                    .and_then(|s| s.trim().parse::<i32>().ok())
                    .unwrap_or(model_count + 1);
                current_model = Model::new(serial);
            }
            "ENDMDL" => {
                structure.models.push(current_model);
                current_model = Model::new(0);
                model_count += 1;
            }
            "ATOM" | "HETATM" => {
                if let Some(atom) = parse_atom_line(&line) {
                    add_atom_to_model(&mut current_model, atom);
                }
            }
            "HEADER" | "TITLE" | "REMARK" | "SEQRES" | "DBREF" | "COMPND" | "SOURCE"
            | "KEYWDS" | "EXPDTA" | "AUTHOR" | "REVDAT" | "JRNL" => {
                structure.header_lines.push(line.clone());
            }
            "END" => {
                break;
            }
            _ => {}
        }
    }

    // If no MODEL/ENDMDL records were found, use the accumulated atoms as model 1
    if !in_model && !current_model.chains.is_empty() {
        structure.models.push(current_model);
    }

    Ok(structure)
}

/// Parse a single ATOM or HETATM line according to the PDB fixed-column format.
///
/// PDB ATOM record format (1-based columns):
///  1- 6  Record name ("ATOM  " or "HETATM")
///  7-11  Atom serial number
/// 13-16  Atom name
/// 17     Alternate location indicator
/// 18-20  Residue name
/// 22     Chain ID
/// 23-26  Residue sequence number
/// 27     Code for insertion of residues
/// 31-38  X coordinate (8.3)
/// 39-46  Y coordinate (8.3)
/// 47-54  Z coordinate (8.3)
/// 55-60  Occupancy (6.2)
/// 61-66  Temperature factor (6.2)
/// 77-78  Element symbol
/// 79-80  Charge
fn parse_atom_line(line: &str) -> Option<Atom> {
    if line.len() < 54 {
        return None;
    }

    let is_hetatm = line.starts_with("HETATM");

    let serial = line.get(6..11)?.trim().parse::<i32>().ok()?;
    let name = line.get(12..16).unwrap_or("    ").to_string();
    let alt_loc = line.as_bytes().get(16).map(|&b| b as char).unwrap_or(' ');
    let res_name = line.get(17..20).unwrap_or("   ").trim().to_string();
    let chain_id = line.as_bytes().get(21).map(|&b| b as char).unwrap_or(' ');
    let res_seq = line.get(22..26)?.trim().parse::<i32>().ok()?;
    let icode = line.as_bytes().get(26).map(|&b| b as char).unwrap_or(' ');

    let x = line.get(30..38)?.trim().parse::<f64>().ok()?;
    let y = line.get(38..46)?.trim().parse::<f64>().ok()?;
    let z = line.get(46..54)?.trim().parse::<f64>().ok()?;

    let occupancy = line
        .get(54..60)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(1.0);
    let temp_factor = line
        .get(60..66)
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(0.0);
    let element = line
        .get(76..78)
        .map(|s| s.trim().to_string())
        .unwrap_or_default();
    let charge = line
        .get(78..80)
        .map(|s| s.trim().to_string())
        .unwrap_or_default();

    Some(Atom::new(
        serial, &name, alt_loc, &res_name, chain_id, res_seq, icode,
        Point3D::new(x, y, z),
        occupancy, temp_factor, &element, &charge, is_hetatm,
    ))
}

/// Insert an atom into the correct chain and residue of the model.
fn add_atom_to_model(model: &mut Model, atom: Atom) {
    let chain_id = atom.chain_id;
    let res_seq = atom.res_seq;
    let icode = atom.icode;
    let res_name = atom.res_name.clone();

    // Find or create chain
    let chain = if let Some(pos) = model.chains.iter().position(|c| c.id == chain_id) {
        &mut model.chains[pos]
    } else {
        model.chains.push(Chain::new(chain_id));
        model.chains.last_mut().unwrap()
    };

    // Find or create residue
    let residue = if let Some(pos) = chain
        .residues
        .iter()
        .position(|r| r.seq_num == res_seq && r.icode == icode)
    {
        &mut chain.residues[pos]
    } else {
        chain
            .residues
            .push(Residue::new(&res_name, res_seq, icode, chain_id));
        chain.residues.last_mut().unwrap()
    };

    residue.atoms.push(atom);
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_PDB: &str = "\
HEADER    TEST STRUCTURE
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00  0.00           C
ATOM      4  O   ALA A   1       4.000   5.000   6.000  1.00  0.00           O
ATOM      5  N   GLY A   2       5.000   6.000   7.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       6.000   7.000   8.000  1.00  0.00           C
ATOM      7  C   GLY A   2       7.000   8.000   9.000  1.00  0.00           C
ATOM      8  O   GLY A   2       8.000   9.000  10.000  1.00  0.00           O
END
";

    #[test]
    fn test_parse_pdb_str() {
        let structure = parse_pdb_str(SAMPLE_PDB, "test").unwrap();
        assert_eq!(structure.name, "test");
        assert_eq!(structure.models.len(), 1);

        let model = structure.first_model().unwrap();
        assert_eq!(model.chains.len(), 1);

        let chain = model.get_chain('A').unwrap();
        assert_eq!(chain.len(), 2);
        assert_eq!(chain.get_sequence(), "AG");

        let res1 = &chain.residues[0];
        assert_eq!(res1.name, "ALA");
        assert_eq!(res1.atoms.len(), 4);

        let ca = res1.get_ca().unwrap();
        assert!((ca.coord.x - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_parse_multimodel() {
        let pdb = "\
MODEL        1
ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C
ENDMDL
MODEL        2
ATOM      1  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00           C
ENDMDL
END
";
        let structure = parse_pdb_str(pdb, "multi").unwrap();
        assert_eq!(structure.models.len(), 2);
        assert_eq!(structure.models[0].serial, 1);
        assert_eq!(structure.models[1].serial, 2);
    }
}
