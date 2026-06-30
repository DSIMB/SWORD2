//! Beta-bridge detection, ladder extension, and sheet assembly.
//!
//! Ported from Dssp.c: Testbridge, Ladder, Extendladder, Sheet, Markstrands.

use super::hbond::test_bond;
use super::types::{Bridge, BridgeType, DsspChain, MAXBRIDGE};

/// Detect all beta-bridges and assemble into ladders and sheets.
pub fn flag_bridges(chain: &mut DsspChain) {
    let mut bridge_table: Vec<Bridge> = Vec::new();
    let mut nbridge: usize = 0;

    // Test all residue pairs for bridges
    let len = chain.len;
    for i in 1..=len {
        if !chain.no_chain_break(i.saturating_sub(1).max(1), (i + 1).min(len)) {
            continue;
        }

        let mut j1: usize = 0;
        let mut j2: usize = 0;
        let mut j = i + 3;

        while j2 == 0 && j < len {
            if chain.no_chain_break(j.saturating_sub(1).max(1), (j + 1).min(len)) {
                let btype = test_bridge(chain, i, j);
                if btype != BridgeType::NoBridge {
                    if j1 == 0 {
                        j1 = j;
                        ladder(&mut bridge_table, &mut nbridge, chain, i, j, btype);
                    } else if j != j1 {
                        j2 = j;
                        ladder(&mut bridge_table, &mut nbridge, chain, i, j, btype);
                    }
                }
            }
            j += 1;
        }
    }

    if nbridge == 0 {
        return;
    }

    // Extend ladders with bulges
    extend_ladder(&mut bridge_table, nbridge, chain);

    // Assemble sheets
    sheet(&mut bridge_table, nbridge);

    // Mark strand residues in chain
    mark_strands(chain, &bridge_table, nbridge);
}

/// Test for a beta-bridge between residues i and j.
fn test_bridge(chain: &DsspChain, i: usize, j: usize) -> BridgeType {
    // Need i-1..i+1 and j-1..j+1 to all be valid
    if i < 2 || j < 2 || i + 1 > chain.len || j + 1 > chain.len {
        return BridgeType::NoBridge;
    }

    // Parallel: (i+1→j) & (j→i-1) OR (j+1→i) & (i→j-1)
    if (test_bond(chain, i + 1, j) && test_bond(chain, j, i - 1))
        || (test_bond(chain, j + 1, i) && test_bond(chain, i, j - 1))
    {
        return BridgeType::Parallel;
    }

    // Antiparallel: (i+1→j-1) & (j+1→i-1) OR (j→i) & (i→j)
    if (test_bond(chain, i + 1, j - 1) && test_bond(chain, j + 1, i - 1))
        || (test_bond(chain, j, i) && test_bond(chain, i, j))
    {
        return BridgeType::Antiparallel;
    }

    BridgeType::NoBridge
}

/// Add a bridge to an existing ladder or create a new one.
#[allow(clippy::needless_range_loop)]
fn ladder(
    table: &mut Vec<Bridge>,
    nbridge: &mut usize,
    chain: &DsspChain,
    i: usize,
    j: usize,
    btype: BridgeType,
) {
    if i >= j {
        return;
    }

    // Try to extend an existing ladder
    for k in 0..*nbridge {
        let b = &table[k];
        let extends = b.btype == btype
            && i == b.ie + 1
            && chain.no_chain_break(b.ie, i)
            && ((j == b.je + 1 && btype == BridgeType::Parallel && chain.no_chain_break(b.je, j))
                || (j == b.jb - 1
                    && btype == BridgeType::Antiparallel
                    && chain.no_chain_break(j, b.jb)));

        if extends {
            let b = &mut table[k];
            b.ie = i;
            if btype == BridgeType::Parallel {
                b.je = j;
            } else {
                b.jb = j;
            }
            return;
        }
    }

    // Create new ladder entry
    if *nbridge >= MAXBRIDGE {
        tracing::warn!("Bridge table overflow");
        return;
    }

    let mut bridge = Bridge::new(MAXBRIDGE);
    bridge.ib = i;
    bridge.ie = i;
    bridge.jb = j;
    bridge.je = j;
    bridge.btype = btype;
    bridge.from = 0;
    bridge.towards = 0;

    table.push(bridge);
    *nbridge += 1;
}

/// Extend ladders by detecting bulges (gaps < 6 residues).
fn extend_ladder(table: &mut [Bridge], nbridge: usize, chain: &DsspChain) {
    // Link ladders that can be joined via bulges
    for i in 0..nbridge {
        let mut j = i + 1;
        while j < nbridge && table[i].towards == 0 {
            let ie = table[i].ie;
            let ib1 = table[j].ib;
            let jb1 = table[j].jb;
            let je1 = table[j].je;

            let mut bulge = chain.no_chain_break(ie, ib1)
                && ib1.saturating_sub(ie) < 6
                && table[j].btype == table[i].btype
                && table[j].from == 0;

            if bulge {
                match table[i].btype {
                    BridgeType::Parallel => {
                        let je_i = table[i].je;
                        bulge = (jb1.saturating_sub(je_i) < 6 && ib1.saturating_sub(ie) < 3
                            || jb1.saturating_sub(je_i) < 3)
                            && chain.no_chain_break(je_i, jb1);
                    }
                    BridgeType::Antiparallel => {
                        let jb_i = table[i].jb;
                        bulge = (jb_i.saturating_sub(je1) < 6 && ib1.saturating_sub(ie) < 3
                            || jb_i.saturating_sub(je1) < 3)
                            && chain.no_chain_break(je1, jb_i);
                    }
                    BridgeType::NoBridge => {}
                }
            }

            if bulge {
                table[i].towards = j + 1; // 1-based index
                table[j].from = i + 1;
            }
            j += 1;
        }
    }

    // Build link sets: each ladder head gets the set of all ladders in its chain
    for i in 0..nbridge {
        if table[i].from == 0 {
            // Head of a ladder chain — collect all linked ladders
            let mut link_set = vec![false; MAXBRIDGE + 1];
            link_set[i + 1] = true; // 1-based
            let mut j = table[i].towards;
            while j != 0 {
                link_set[j] = true;
                j = table[j - 1].towards;
            }
            table[i].link_set = link_set.clone();

            // Propagate to all members
            j = table[i].towards;
            while j != 0 {
                table[j - 1].link_set = link_set.clone();
                j = table[j - 1].towards;
            }
        }
    }
}

/// Assemble ladders sharing residues into sheets, assign sheet/ladder labels.
#[allow(clippy::needless_range_loop)]
fn sheet(table: &mut [Bridge], nbridge: usize) {
    // ladderset: which ladders still need sheet assignment (1-based)
    let mut ladder_set = vec![true; nbridge + 1]; // [1..nbridge]
    ladder_set[0] = false;

    let mut sheet_label = b'A' - 1;
    let mut strand_label: u8 = b'A' - 1; // 'A' through 'Z' for strand labels

    while ladder_set[1..=nbridge].iter().any(|&x| x) {
        sheet_label += 1;
        if sheet_label == b'Z' + 1 {
            sheet_label = b'a';
        }
        if sheet_label > b'z' {
            sheet_label = b'A';
        }

        // Find first unassigned ladder
        let first = match (1..=nbridge).find(|&l| ladder_set[l]) {
            Some(l) => l,
            None => break,
        };

        // Start sheet with this ladder's link set
        let mut sheet_set = table[first - 1].link_set.clone();
        // Remove from ladder_set
        for k in 1..=nbridge {
            if sheet_set.get(k).copied().unwrap_or(false) {
                ladder_set[k] = false;
            }
        }

        // Expand: find any remaining ladder that shares residues with the sheet
        loop {
            let mut expanded = false;
            for l1 in 1..=nbridge {
                if !sheet_set.get(l1).copied().unwrap_or(false) {
                    continue;
                }
                for l2 in 1..=nbridge {
                    if !ladder_set[l2] {
                        continue;
                    }
                    if link(table, l1 - 1, l2 - 1) {
                        // Add l2's link set to sheet
                        let ls = table[l2 - 1].link_set.clone();
                        for k in 0..ls.len().min(sheet_set.len()) {
                            if ls[k] {
                                sheet_set[k] = true;
                                ladder_set[k] = false;
                            }
                        }
                        expanded = true;
                    }
                }
            }
            if !expanded {
                break;
            }
        }

        // Assign labels
        for i in 1..=nbridge {
            if !sheet_set.get(i).copied().unwrap_or(false) {
                continue;
            }
            if table[i - 1].from != 0 {
                continue;
            }
            // Head of a ladder chain
            strand_label += 1;
            if strand_label > b'Z' {
                strand_label = b'A';
            }

            let label = if table[i - 1].btype == BridgeType::Parallel {
                (strand_label + 32) as char // lowercase
            } else {
                strand_label as char
            };

            table[i - 1].sheet_name = sheet_label as char;
            table[i - 1].ladder_name = label;

            // Propagate to linked ladders
            let mut j = table[i - 1].towards;
            while j != 0 {
                table[j - 1].sheet_name = sheet_label as char;
                table[j - 1].ladder_name = label;
                table[j - 1].link_set = sheet_set.clone();
                j = table[j - 1].towards;
            }
            table[i - 1].link_set = sheet_set.clone();
        }
    }
}

/// Check if two ladders share any residues.
fn link(table: &[Bridge], l1: usize, l2: usize) -> bool {
    let a = &table[l1];
    let b = &table[l2];
    (a.ie >= b.ib && a.ib <= b.ie)
        || (a.ie >= b.jb && a.ib <= b.je)
        || (a.je >= b.ib && a.jb <= b.ie)
        || (a.je >= b.jb && a.jb <= b.je)
}

fn partner_index(base: usize, add: usize, sub: usize) -> usize {
    base.checked_add(add)
        .and_then(|value| value.checked_sub(sub))
        .unwrap_or(0)
}

/// Mark strand residues in the chain based on bridge table.
fn mark_strands(chain: &mut DsspChain, table: &[Bridge], nbridge: usize) {
    for i in 0..nbridge {
        if table[i].from != 0 {
            continue; // Only process ladder heads
        }

        // Determine which beta column (beta1=6, beta2=7) to use for i-strand and j-strand
        // Check if beta1 column is free for all residues in this ladder chain
        let mut ib0 = chain.len;
        let mut ie0 = 0usize;
        let mut jb0 = chain.len;
        let mut je0 = 0usize;

        // Check occupancy of beta columns
        let mut i_col_free = [true, true]; // beta1, beta2
        let mut j_col_free = [true, true];

        let mut j_idx = i;
        loop {
            let b = &table[j_idx];
            for l in b.ib..=b.ie {
                if chain.get(l).ss[6] != ' ' {
                    i_col_free[0] = false;
                }
                if chain.get(l).ss[7] != ' ' {
                    i_col_free[1] = false;
                }
            }
            for l in b.jb..=b.je {
                if chain.get(l).ss[6] != ' ' {
                    j_col_free[0] = false;
                }
                if chain.get(l).ss[7] != ' ' {
                    j_col_free[1] = false;
                }
            }
            if b.ib < ib0 {
                ib0 = b.ib;
            }
            if b.ie > ie0 {
                ie0 = b.ie;
            }
            if b.jb < jb0 {
                jb0 = b.jb;
            }
            if b.je > je0 {
                je0 = b.je;
            }

            if b.towards == 0 {
                break;
            }
            j_idx = b.towards - 1;
        }

        let betai: usize = if i_col_free[0] { 6 } else { 7 }; // ss index
        let betaj: usize = if j_col_free[0] { 6 } else { 7 };

        // Assign ladder labels and partners
        j_idx = i;
        loop {
            let ladder_name = table[j_idx].ladder_name;
            let btype = table[j_idx].btype;
            let ib = table[j_idx].ib;
            let ie = table[j_idx].ie;
            let jb = table[j_idx].jb;
            let je = table[j_idx].je;

            for l in ib..=ie {
                chain.get_mut(l).ss[betai] = ladder_name;
                chain.get_mut(l).partner[betai - 6] = match btype {
                    BridgeType::Parallel => partner_index(jb, l, ib),
                    BridgeType::Antiparallel => partner_index(je, ib, l),
                    BridgeType::NoBridge => 0,
                };
            }
            for l in jb..=je {
                chain.get_mut(l).ss[betaj] = ladder_name;
                chain.get_mut(l).partner[betaj - 6] = match btype {
                    BridgeType::Parallel => partner_index(ib, l, jb),
                    BridgeType::Antiparallel => partner_index(ie, jb, l),
                    BridgeType::NoBridge => 0,
                };
            }

            let towards = table[j_idx].towards;
            if towards == 0 {
                break;
            }
            j_idx = towards - 1;
        }

        // Mark symbol column: 'E' if ladder spans >1 residue, 'B' for isolated bridge
        let cc = if ib0 == ie0 { 'B' } else { 'E' };
        for l in ib0..=ie0 {
            if chain.get(l).ss[0] != 'E' {
                chain.get_mut(l).ss[0] = cc;
            }
        }
        let cc = if jb0 == je0 { 'B' } else { 'E' };
        for l in jb0..=je0 {
            if chain.get(l).ss[0] != 'E' {
                chain.get_mut(l).ss[0] = cc;
            }
        }
    }

    // Assign sheet labels to all residues in bridges
    for b in &table[..nbridge] {
        for l in b.ib..=b.ie {
            chain.get_mut(l).sheet_label = b.sheet_name;
        }
        for l in b.jb..=b.je {
            chain.get_mut(l).sheet_label = b.sheet_name;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dssp::types::{BackboneResidue, MAXBRIDGE};

    #[test]
    fn mark_strands_does_not_underflow_antiparallel_partner_indices() {
        let mut chain = DsspChain::new();
        for _ in 0..6 {
            chain.push(BackboneResidue::default());
        }

        let mut bridge = Bridge::new(MAXBRIDGE);
        bridge.sheet_name = 'A';
        bridge.ladder_name = 'A';
        bridge.btype = BridgeType::Antiparallel;
        bridge.ib = 1;
        bridge.ie = 1;
        bridge.jb = 3;
        bridge.je = 5;

        mark_strands(&mut chain, &[bridge], 1);

        assert_eq!(chain.get(5).partner[0], 0);
    }
}
