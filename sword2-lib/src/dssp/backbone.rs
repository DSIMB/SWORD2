//! Backbone extraction from PDB and H-atom synthesis.
//!
//! Reads a cleaned PDB file and builds a DsspChain with backbone atom coordinates.
//! Synthesizes backbone H atoms using the C(i-1)->O(i-1) direction.

use std::path::Path;

use anyhow::{Context, Result};

use crate::pdb::types::Point3D;

use super::types::{BackboneResidue, DsspChain, BREAKDIST};

/// One-letter code lookup from three-letter residue name.
fn one_letter_code(three: &str) -> char {
    match three.trim() {
        "ALA" => 'A', "ARG" => 'R', "ASN" => 'N', "ASP" => 'D',
        "CYS" => 'C', "GLU" => 'E', "GLN" => 'Q', "GLY" => 'G',
        "HIS" => 'H', "ILE" => 'I', "LEU" => 'L', "LYS" => 'K',
        "MET" => 'M', "PHE" => 'F', "PRO" => 'P', "SER" => 'S',
        "THR" => 'T', "TRP" => 'W', "TYR" => 'Y', "VAL" => 'V',
        "ASX" => 'B', "GLX" => 'Z', "CSH" | "CSS" | "CYX" => 'C',
        _ => 'X',
    }
}

/// Parse a cleaned PDB file and build a DsspChain with backbone atoms and synthesized H.
pub fn extract_backbone(pdb_path: &Path) -> Result<DsspChain> {
    let content = std::fs::read_to_string(pdb_path)
        .with_context(|| format!("Failed to read PDB file: {}", pdb_path.display()))?;

    let mut chain = DsspChain::new();

    // Collect ATOM records grouped by residue
    let mut current_res: Option<ResidueBuilder> = None;

    for line in content.lines() {
        if !line.starts_with("ATOM") || line.len() < 54 {
            continue;
        }

        let atom_name = &line[12..16];
        let res_name = line[17..20].trim();
        let chain_id = line.as_bytes().get(21).map(|&b| b as char).unwrap_or(' ');
        let res_seq: i32 = line[22..26].trim().parse().unwrap_or(0);
        let icode = line.as_bytes().get(26).map(|&b| b as char).unwrap_or(' ');

        let x: f64 = line[30..38].trim().parse().unwrap_or(0.0);
        let y: f64 = line[38..46].trim().parse().unwrap_or(0.0);
        let z: f64 = line[46..54].trim().parse().unwrap_or(0.0);
        let coord = Point3D::new(x, y, z);

        // Check if we've moved to a new residue
        let new_residue = match &current_res {
            Some(r) => r.res_seq != res_seq || r.icode != icode,
            None => true,
        };

        if new_residue {
            // Finalize previous residue
            if let Some(builder) = current_res.take() {
                if let Some(res) = builder.finalize() {
                    add_residue_to_chain(&mut chain, res);
                }
            }
            current_res = Some(ResidueBuilder::new(res_name, chain_id, res_seq, icode));
        }

        // Store backbone atom coordinates
        if let Some(ref mut builder) = current_res {
            let name = atom_name.trim();
            // Only take first alternate location
            let altloc = line.as_bytes().get(16).map(|&b| b as char).unwrap_or(' ');
            if altloc != ' ' && altloc != 'A' {
                continue;
            }
            match name {
                "N" => builder.n = Some(coord),
                "CA" => builder.ca = Some(coord),
                "C" => builder.c = Some(coord),
                "O" => builder.o = Some(coord),
                _ => {}
            }
        }
    }

    // Finalize last residue
    if let Some(builder) = current_res.take() {
        if let Some(res) = builder.finalize() {
            add_residue_to_chain(&mut chain, res);
        }
    }

    Ok(chain)
}

/// Add a residue to the chain, inserting chain breaks as needed and synthesizing H.
fn add_residue_to_chain(chain: &mut DsspChain, mut res: BackboneResidue) {
    // Default: H at N position (no amide H for first residue or after break)
    res.h = res.n;
    res.has_h = false;

    if chain.len > 0 {
        let prev = chain.get(chain.len);

        // Check for chain break
        if prev.aa != '!' {
            let cn_dist = prev.c.distance_to(&res.n);
            if cn_dist > BREAKDIST {
                // Insert chain break marker
                chain.push(BackboneResidue::chain_break());
            }
        }

        // Synthesize backbone H: H(i) = N(i) + normalize(C(i-1) - O(i-1))
        let prev_idx = chain.len;
        let prev = chain.get(prev_idx);
        if prev.aa != '!' && res.aa != 'P' {
            let co = prev.c - prev.o;
            let co_len = (co.x * co.x + co.y * co.y + co.z * co.z).sqrt();
            if co_len > 0.0 {
                res.h = Point3D::new(
                    res.n.x + co.x / co_len,
                    res.n.y + co.y / co_len,
                    res.n.z + co.z / co_len,
                );
                res.has_h = true;
            }
        }
    }

    chain.push(res);
}

/// Builder for collecting backbone atoms of a single residue.
struct ResidueBuilder {
    res_name: String,
    chain_id: char,
    res_seq: i32,
    icode: char,
    n: Option<Point3D>,
    ca: Option<Point3D>,
    c: Option<Point3D>,
    o: Option<Point3D>,
}

impl ResidueBuilder {
    fn new(res_name: &str, chain_id: char, res_seq: i32, icode: char) -> Self {
        Self {
            res_name: res_name.to_string(),
            chain_id,
            res_seq,
            icode,
            n: None,
            ca: None,
            c: None,
            o: None,
        }
    }

    /// Finalize into a BackboneResidue if all backbone atoms present.
    fn finalize(self) -> Option<BackboneResidue> {
        let n = self.n?;
        let ca = self.ca?;
        let c = self.c?;
        let o = self.o?;

        let aa = one_letter_code(&self.res_name);
        if aa == '-' {
            return None;
        }

        // Build aaident: 4 chars resnum + icode + chain_id (matching DSSP format)
        let mut aaident = [b' '; 6];
        let resnum_str = format!("{:>4}", self.res_seq);
        for (i, b) in resnum_str.bytes().take(4).enumerate() {
            aaident[i] = b;
        }
        aaident[4] = self.icode as u8;
        aaident[5] = self.chain_id as u8;

        let mut three_letter = [b' '; 4];
        for (i, b) in self.res_name.bytes().take(3).enumerate() {
            three_letter[i] = b;
        }

        Some(BackboneResidue {
            aa,
            aaident,
            three_letter,
            n,
            ca,
            c,
            o,
            h: n, // Will be overwritten by add_residue_to_chain
            has_h: false,
            ..BackboneResidue::default()
        })
    }
}
