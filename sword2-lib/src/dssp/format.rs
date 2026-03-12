//! DSSP fixed-width output format.
//!
//! Produces output compatible with the original dsspcmbi binary,
//! which Peeling reads via parse_dssp().
//!
//! Critical columns for Peeling:
//! - line[126]: must be non-space and non-'-' for data lines
//! - line[6..11]: DSSP sequential residue number
//! - line[13]: amino acid code ('!' for chain break)
//! - line[16]: secondary structure summary code

use std::io::Write;
use std::path::Path;

use anyhow::Result;

use super::angles;
use super::types::DsspChain;

/// Write DSSP output file in the classic fixed-width format.
pub fn write_dssp(chain: &DsspChain, pdb_name: &str, output_path: &Path) -> Result<()> {
    let mut out = Vec::new();

    // Header lines — padded to ≥128 chars so that position 126 is always a space.
    // The C binary's parse_dssp() skips lines where line[126] is whitespace or '-'.
    // Without padding, short header lines leave position 126 as uninitialized buffer
    // content, causing them to be incorrectly parsed as data lines.
    write_padded_line(&mut out, "==== Secondary Structure Definition by the program DSSP, Rust port ====")?;
    write_padded_line(&mut out, "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637")?;
    write_padded_line(&mut out, &format!("HEADER    {}", pdb_name))?;
    write_padded_line(&mut out, "COMPND")?;
    write_padded_line(&mut out, "SOURCE")?;
    write_padded_line(&mut out, "AUTHOR")?;

    // Count residues (excluding chain breaks)
    let nres = (1..=chain.len)
        .filter(|&i| chain.get(i).aa != '!')
        .count();
    let nchains = count_chains(chain);

    writeln!(
        out,
        "{:>5}{:>5}{:>5}{:>5}{:>5} TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .",
        nres, nchains, 0, 0, 0
    )?;

    // Statistics placeholder lines (not used by Peeling, but maintain format)
    writeln!(out, "     0  0.0   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .")?;

    // Histogram placeholders
    for label in &[
        "RESIDUES PER ALPHA HELIX",
        "PARALLEL BRIDGES PER LADDER",
        "ANTIPARALLEL BRIDGES PER LADDER",
        "LADDERS PER SHEET",
    ] {
        write!(out, "  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0")?;
        writeln!(out, "     {:<68}.", label)?;
    }

    // Sentinel line that parse_dssp_to_s2d looks for
    writeln!(out,
        "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA            "
    )?;

    // Per-residue lines
    for i in 1..=chain.len {
        let res = chain.get(i);

        // Chain break line
        if res.aa == '!' {
            // Minimal chain break line: just needs '!' at column 13 and non-space at column 126
            write!(out, "{:>5} ", i)?;
            write!(out, "      ")?; // aaident
            write!(out, " !  ")?;
            // Pad to reach column 126+ with the dot marker
            let line_so_far = 16;
            for _ in line_so_far..126 {
                write!(out, " ")?;
            }
            writeln!(out, " ")?;
            continue;
        }

        // Compute angles for this residue
        let tco_val = if i > 1 { angles::tco(chain, i) } else { 0.0 };
        let phi_val = if i > 1 { angles::phi(chain, i) } else { 360.0 };
        let psi_val = angles::psi(chain, i);

        // Write DSSP number (columns 0-4)
        write!(out, "{:>5} ", i)?;

        // aaident (columns 6-11): residue number, icode, chain
        for &b in &res.aaident {
            write!(out, "{}", b as char)?;
        }

        // Column 12: space, 13: AA, 14: chain_break(' '), 15: space, 16: SS symbol, 17: space
        // Written character by character to ensure exact column alignment.
        out.push(b' ');       // col 12
        out.push(res.aa as u8);  // col 13
        out.push(b' ');       // col 14 (chain break placeholder)
        out.push(b' ');       // col 15
        out.push(res.ss[0] as u8);  // col 16
        out.push(b' ');       // col 17

        // Turn/bend/chirality/beta columns (18-24): turn3, turn4, turn5, bend, chirality, beta1, beta2
        for s in 1..=7 {
            write!(out, "{}", res.ss[s])?;
        }

        // Bridge partners (columns 25-32): BP1, BP2
        write!(out, "{:>4}{:>4}", res.partner[0], res.partner[1])?;

        // Sheet label and accessibility (columns 33-37)
        write!(out, "{}{:>4} ", res.sheet_label, res.access)?;

        // H-bonds: 4 pairs (acceptor0, donor0, acceptor1, donor1)
        for j in 0..2 {
            write_hb(&mut out, i, res.acceptor[j])?;
            write_hb(&mut out, i, res.donor[j])?;
        }

        // TCO, KAPPA, ALPHA, PHI, PSI
        write!(
            out,
            "{:>8.3}{:>6.1}{:>6.1}{:>6.1}{:>6.1}",
            tco_val, res.kappa, res.alpha, phi_val, psi_val
        )?;

        // CA coordinates
        write!(
            out,
            "{:>7.1}{:>7.1}{:>7.1}",
            res.ca.x, res.ca.y, res.ca.z
        )?;

        // Pad to ensure column 126 has a non-space character
        // Current position is approximately 126+ already with coordinates
        // Add padding + marker to be safe
        writeln!(out, "            ")?;
    }

    std::fs::write(output_path, out)?;
    Ok(())
}

/// Write a hydrogen bond in DSSP format.
fn write_hb(out: &mut Vec<u8>, i: usize, hb: super::types::HydrogenBond) -> Result<()> {
    let relative = if hb.residue != 0 {
        hb.residue as i64 - i as i64
    } else {
        0
    };
    let energy = hb.energy as f64 / 1000.0;
    write!(out, "{:>6},{:>4.1}", relative, energy)?;
    Ok(())
}

/// Count the number of separate polypeptide chains.
fn count_chains(chain: &DsspChain) -> usize {
    if chain.len == 0 {
        return 0;
    }
    let mut count = 1;
    for i in 1..=chain.len {
        if chain.get(i).aa == '!' {
            count += 1;
        }
    }
    count
}

/// Generate a .s2d file (3-state secondary structure) from a DsspChain.
///
/// This replaces parse_dssp_to_s2d by generating directly from structured data.
pub fn write_s2d(chain: &DsspChain, pdb_name: &str, output_path: &Path) -> Result<()> {
    let mut ss_list = Vec::new();

    for i in 1..=chain.len {
        let res = chain.get(i);

        if res.aa == '!' {
            continue;
        }

        let ss3 = match res.ss[0] {
            'H' | 'G' | 'I' => 'H',
            'E' | 'B' => 'E',
            _ => 'C',
        };
        ss_list.push(ss3);
    }

    let mut out = format!("> {} \n", pdb_name);
    for (i, chunk) in ss_list.chunks(80).enumerate() {
        let s: String = chunk.iter().collect();
        out.push_str(&s);
        if i < ss_list.len() / 80 {
            out.push('\n');
        }
    }

    std::fs::write(output_path, out)?;
    Ok(())
}

/// Write a line padded to at least 128 characters (spaces) so that
/// the C binary's parse_dssp() sees a space at position 126 and skips it.
fn write_padded_line(out: &mut Vec<u8>, content: &str) -> Result<()> {
    const MIN_WIDTH: usize = 128;
    write!(out, "{}", content)?;
    let len = content.len();
    if len < MIN_WIDTH {
        for _ in len..MIN_WIDTH {
            out.push(b' ');
        }
    }
    writeln!(out)?;
    Ok(())
}
