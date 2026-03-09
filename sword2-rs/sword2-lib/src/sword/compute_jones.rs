//! Jones overlap scoring between domain assignments.
//!
//! Port of `ComputeJones.pl`. Determines whether two domain decompositions
//! overlap by ≥85% (using permutation-based comparison of domain assignments).

use std::path::Path;
use std::fs;

/// Compute the Jones overlap criterion between two domain delineations.
///
/// Returns `(criterion, percentage)`:
///   - `criterion`: 1 if overlap ≥ 85%, else 0.
///   - `percentage`: the best overlap percentage found.
///
/// `pdb_name` is the PDB identifier (with chain suffix and trailing "_1_").
/// `delineation1` and `delineation2` are space-separated domain boundary strings
/// like "1-50 51-100" or "1-30;60-100 31-59".
/// `dir_data` is the path to PDBs_Clean/ directory.
pub fn compute_jones(
    pdb_name: &str,
    delineation1: &str,
    delineation2: &str,
    dir_data: &str,
) -> (i32, f64) {
    let tab_del1: Vec<&str> = delineation1.split_whitespace().collect();
    let tab_del2: Vec<&str> = delineation2.split_whitespace().collect();

    // The original Perl script expects first two tokens to be pdb and level info.
    // When called from parse_measure, the input may not have that prefix.
    // Handle both cases gracefully.
    let skip1 = if tab_del1.len() > 2 && tab_del1[0].chars().any(|c| c.is_alphabetic()) { 2 } else { 0 };
    let skip2 = if tab_del2.len() > 2 && tab_del2[0].chars().any(|c| c.is_alphabetic()) { 2 } else { 0 };

    let del1_domains = &tab_del1[skip1..];
    let del2_domains = &tab_del2[skip2..];

    let (dom_peel_tokens, dom_auth_tokens) = if del1_domains.len() < del2_domains.len() {
        (del2_domains, del1_domains)
    } else {
        (del1_domains, del2_domains)
    };

    let dom_auth = dom_auth_tokens.join(" ");
    let dom_peel = dom_peel_tokens.join(" ");

    // Read residue numbers from the PDB file
    let pdb_base = &pdb_name[..pdb_name.len().saturating_sub(3)]; // strip "_1_"
    let pdb_file_path = Path::new(dir_data).join(pdb_base).join(format!("{}.pdb", pdb_base));
    let tab_num_aa = read_ca_residue_numbers(&pdb_file_path);

    if tab_num_aa.is_empty() {
        return (0, 0.0);
    }

    let length_authors = tab_num_aa.len();

    // Parse auth domain delineation and assign domain IDs to each residue position
    let auth_domains: Vec<&str> = dom_auth.split_whitespace().collect();
    let mut tab_protein_domain_authors = vec![-1i32; length_authors];

    for pos in 0..length_authors {
        for (num_dom, domain_str) in auth_domains.iter().enumerate() {
            let pu_parts: Vec<&str> = domain_str.split(';').collect();
            for pu in &pu_parts {
                if let Some((start, end)) = parse_range(pu) {
                    if tab_num_aa[pos] >= start && tab_num_aa[pos] <= end {
                        tab_protein_domain_authors[pos] = num_dom as i32;
                    }
                }
            }
        }
    }

    // Parse peel domain delineation
    let peel_domains: Vec<&str> = dom_peel.split_whitespace().collect();
    let n_peel = peel_domains.len();

    // Find peel range to compute length
    let mut peel_start = i32::MAX;
    let mut peel_end = i32::MIN;
    for domain_str in &peel_domains {
        for seg in domain_str.split(';') {
            if let Some((s, e)) = parse_range(seg) {
                peel_start = peel_start.min(s);
                peel_end = peel_end.max(e);
            }
        }
    }

    if peel_start == i32::MAX {
        return (0, 0.0);
    }

    let peel_length = (peel_end - peel_start) as usize + 1;

    // Assign domain IDs to each position in the peel decomposition
    let mut tab_protein_domain_peel = vec![-1i32; peel_length];
    for pos in 0..peel_length {
        for (num_dom, domain_str) in peel_domains.iter().enumerate() {
            let pu_parts: Vec<&str> = domain_str.split(';').collect();
            for pu in &pu_parts {
                if let Some((s, e)) = parse_range(pu) {
                    let adj_s = (s - peel_start) as usize;
                    let adj_e = (e - peel_start) as usize;
                    if pos >= adj_s && pos <= adj_e {
                        tab_protein_domain_peel[pos] = num_dom as i32;
                    }
                }
            }
        }
    }

    // Skip if too many permutations (> 7 domains peel vs fewer auth)
    if n_peel > 7 && auth_domains.len() < n_peel {
        return (0, 0.0);
    }

    // Try all permutations to find the best Jones overlap
    let domain_ids: Vec<usize> = (0..n_peel).collect();
    let mut best_jones: i32 = 0;
    let threshold = 85;

    for_each_permutation(&domain_ids, &mut |perm| {
        let mut jones = 0i32;
        let mut total = 0i32;

        for pos in 0..tab_protein_domain_peel.len() {
            if pos >= tab_protein_domain_authors.len() {
                break;
            }
            let peel_dom = tab_protein_domain_peel[pos];
            if peel_dom < 0 {
                continue;
            }
            let remapped = perm[peel_dom as usize] as i32;
            if tab_protein_domain_authors[pos] == -1 {
                continue;
            }
            total += 1;
            if remapped == tab_protein_domain_authors[pos] {
                jones += 1;
            }
        }

        if jones > best_jones {
            best_jones = jones;
        }
        // Early exit if we already exceed threshold
        best_jones >= threshold as i32
    });

    let total = tab_protein_domain_authors
        .iter()
        .filter(|&&d| d != -1)
        .count() as i32;

    if total <= 0 {
        return (0, 0.0);
    }

    let percentage = best_jones as f64 / total as f64 * 100.0;
    let criterion = if percentage >= 85.0 { 1 } else { 0 };

    (criterion, percentage)
}

/// Parse a range string like "10-50" into (10, 50).
fn parse_range(s: &str) -> Option<(i32, i32)> {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() == 2 {
        let start = parts[0].parse::<i32>().ok()?;
        let end = parts[1].parse::<i32>().ok()?;
        Some((start, end))
    } else {
        None
    }
}

/// Read CA residue numbers from a PDB file.
fn read_ca_residue_numbers(pdb_path: &Path) -> Vec<i32> {
    let content = match fs::read_to_string(pdb_path) {
        Ok(c) => c,
        Err(_) => return Vec::new(),
    };

    let mut nums = Vec::new();
    for line in content.lines() {
        if !line.starts_with("ATOM") {
            if line.starts_with("TER") || line.starts_with("ENDMDL") {
                break;
            }
            continue;
        }
        // Check atom name (columns 13-16)
        if line.len() < 54 {
            continue;
        }
        let atom_name = line[12..16].trim();
        if atom_name != "CA" {
            continue;
        }
        // Check for standard amino acids
        let aa = line[17..20].trim();
        if !is_standard_aa(aa) {
            continue;
        }
        // Residue number (columns 23-26)
        if let Ok(num) = line[22..26].trim().parse::<i32>() {
            nums.push(num);
        }
    }
    nums
}

fn is_standard_aa(aa: &str) -> bool {
    matches!(
        aa,
        "ALA" | "CYS" | "ASP" | "GLU" | "PHE" | "GLY" | "HIS" | "ILE" | "LYS" | "LEU"
            | "MET" | "ASN" | "PRO" | "GLN" | "ARG" | "SER" | "THR" | "VAL" | "TRP" | "TYR"
            | "UNK"
    )
}

/// Iterate over all permutations of `items`, calling `f` for each.
/// If `f` returns `true`, stop early.
fn for_each_permutation<T: Clone>(items: &[T], f: &mut dyn FnMut(&[T]) -> bool) {
    let n = items.len();
    if n == 0 {
        return;
    }
    let mut perm: Vec<T> = items.to_vec();
    let mut c = vec![0usize; n];
    if f(&perm) {
        return;
    }

    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 {
                perm.swap(0, i);
            } else {
                perm.swap(c[i], i);
            }
            if f(&perm) {
                return;
            }
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_range() {
        assert_eq!(parse_range("1-100"), Some((1, 100)));
        assert_eq!(parse_range("50-200"), Some((50, 200)));
        assert_eq!(parse_range("abc"), None);
    }

    #[test]
    fn test_permutation_count() {
        let items = vec![0, 1, 2];
        let mut count = 0;
        for_each_permutation(&items, &mut |_| {
            count += 1;
            false
        });
        assert_eq!(count, 6); // 3! = 6
    }
}
