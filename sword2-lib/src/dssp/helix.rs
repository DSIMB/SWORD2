//! Turn detection (3/4/5-turns) and helix/bend/symbol assignment.
//!
//! Ported from Dssp.c: Flagturn and Flagsymbol.

use super::hbond::test_bond;
use super::types::DsspChain;

// SS column indices: 0=symbol, 1=turn3, 2=turn4, 3=turn5, 4=bend, 5=chirality, 6=beta1, 7=beta2

/// Detect turns and assign helix/bend symbols.
pub fn flag_turns(chain: &mut DsspChain) {
    let len = chain.len;

    // Detect 3-, 4-, 5-turns
    // turn3 -> ss[1], k=3, marker='3'
    // turn4 -> ss[2], k=4, marker='4'
    // turn5 -> ss[3], k=5, marker='5'
    for (turn_idx, k, cc) in [(1usize, 3usize, '3'), (2, 4, '4'), (3, 5, '5')] {
        if len < k + 1 {
            continue;
        }
        for i in 1..=len - k {
            if chain.no_chain_break(i, i + k) && test_bond(chain, i + k, i) {
                // Mark the end of the turn
                chain.get_mut(i + k).ss[turn_idx] = '<';

                // Mark interior residues
                for j in 1..k {
                    if chain.get(i + j).ss[turn_idx] == ' ' {
                        chain.get_mut(i + j).ss[turn_idx] = cc;
                    }
                }

                // Mark the start
                let start_ss = chain.get(i).ss[turn_idx];
                if start_ss == '<' {
                    chain.get_mut(i).ss[turn_idx] = 'X';
                } else {
                    chain.get_mut(i).ss[turn_idx] = '>';
                }
            }
        }
    }

    // Mark bends (kappa > 70 degrees)
    for i in 1..=len {
        let kappa = chain.get(i).kappa;
        if kappa != 360.0 && kappa > 70.0 {
            chain.get_mut(i).ss[4] = 'S';
        }
    }

    // Assign final SS symbols
    flag_symbol(chain);
}

/// Assign the summary SS symbol (column 0) based on turns and bridges.
fn flag_symbol(chain: &mut DsspChain) {
    let len = chain.len;

    let is_nh = |c: char| -> bool { c == '>' || c == 'X' };

    // Alpha helix (H): two consecutive 4-turn starts
    if len >= 5 {
        for i in 2..=len - 4 {
            if is_nh(chain.get(i - 1).ss[2]) && is_nh(chain.get(i).ss[2]) {
                for j in i..=i + 3 {
                    chain.get_mut(j).ss[0] = 'H';
                }
            }
        }
    }

    // 3₁₀ helix (G): two consecutive 3-turn starts, only if not already H
    if len >= 4 {
        for i in 2..=len - 3 {
            if is_nh(chain.get(i - 1).ss[1]) && is_nh(chain.get(i).ss[1]) {
                let empty = (i..=i + 2).all(|j| {
                    let s = chain.get(j).ss[0];
                    s == 'G' || s == ' '
                });
                if empty {
                    for j in i..=i + 2 {
                        chain.get_mut(j).ss[0] = 'G';
                    }
                }
            }
        }
    }

    // Pi helix (I): two consecutive 5-turn starts, only if not already H/G
    if len >= 6 {
        for i in 2..=len - 5 {
            if is_nh(chain.get(i - 1).ss[3]) && is_nh(chain.get(i).ss[3]) {
                let empty = (i..=i + 4).all(|j| {
                    let s = chain.get(j).ss[0];
                    s == 'I' || s == ' '
                });
                if empty {
                    for j in i..=i + 4 {
                        chain.get_mut(j).ss[0] = 'I';
                    }
                }
            }
        }
    }

    // Remaining: T for turn residues, S for bend
    for i in 2..len {
        if chain.get(i).ss[0] != ' ' {
            continue;
        }

        let mut cc = ' ';
        // Check if this residue is inside any n-turn
        // turn3 (k=3): check if any of i-1, i-2 starts a 3-turn
        // turn4 (k=4): check if any of i-1, i-2, i-3 starts a 4-turn
        // turn5 (k=5): check if any of i-1, i-2, i-3, i-4 starts a 5-turn
        let mut j = 1usize; // number of residues to look back
        for turn_idx in 1..=3usize {
            j += 1; // j = 2, 3, 4 for turn3, turn4, turn5
            for k in 1..=j {
                if i > k {
                    let s = chain.get(i - k).ss[turn_idx];
                    if is_nh(s) {
                        cc = 'T';
                    }
                }
            }
        }

        if cc == ' ' {
            cc = chain.get(i).ss[4]; // bend 'S' or ' '
        }

        chain.get_mut(i).ss[0] = cc;
    }
}
