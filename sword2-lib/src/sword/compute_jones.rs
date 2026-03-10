//! Jones overlap scoring between domain assignments.
//!
//! Determines whether two domain decompositions overlap by ≥85%.
//! Uses the Hungarian (Kuhn-Munkres) algorithm in O(n³) to find the optimal
//! domain-to-domain mapping, replacing the previous O(n!) permutation search.

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
    compute_jones_with_cache(pdb_name, delineation1, delineation2, dir_data, None)
}

/// Compute the Jones overlap, optionally using pre-loaded CA residue numbers.
pub fn compute_jones_with_cache(
    pdb_name: &str,
    delineation1: &str,
    delineation2: &str,
    dir_data: &str,
    cached_residue_nums: Option<&[i32]>,
) -> (i32, f64) {
    let tab_del1: Vec<&str> = delineation1.split_whitespace().collect();
    let tab_del2: Vec<&str> = delineation2.split_whitespace().collect();

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

    // Read residue numbers from cache or PDB file
    let tab_num_aa: Vec<i32>;
    let residue_nums = if let Some(cached) = cached_residue_nums {
        cached
    } else {
        let pdb_base = &pdb_name[..pdb_name.len().saturating_sub(3)];
        let pdb_file_path = Path::new(dir_data).join(pdb_base).join(format!("{}.pdb", pdb_base));
        tab_num_aa = read_ca_residue_numbers(&pdb_file_path);
        &tab_num_aa
    };

    if residue_nums.is_empty() {
        return (0, 0.0);
    }

    let length_authors = residue_nums.len();

    // Parse auth domain delineation and assign domain IDs to each residue position
    let auth_domains: Vec<&str> = dom_auth.split_whitespace().collect();
    let n_auth = auth_domains.len();
    let mut tab_protein_domain_authors = vec![-1i32; length_authors];

    for pos in 0..length_authors {
        for (num_dom, domain_str) in auth_domains.iter().enumerate() {
            let pu_parts: Vec<&str> = domain_str.split(';').collect();
            for pu in &pu_parts {
                if let Some((start, end)) = parse_range(pu) {
                    if residue_nums[pos] >= start && residue_nums[pos] <= end {
                        tab_protein_domain_authors[pos] = num_dom as i32;
                    }
                }
            }
        }
    }

    // Parse peel domain delineation
    let peel_domains: Vec<&str> = dom_peel.split_whitespace().collect();
    let n_peel = peel_domains.len();

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

    // Build overlap cost matrix for Hungarian algorithm
    // overlap[i][j] = number of residues assigned to peel domain i AND auth domain j
    let n = n_peel.max(n_auth);
    let mut overlap = vec![vec![0i32; n]; n];

    for pos in 0..tab_protein_domain_peel.len() {
        if pos >= tab_protein_domain_authors.len() {
            break;
        }
        let peel_dom = tab_protein_domain_peel[pos];
        let auth_dom = tab_protein_domain_authors[pos];
        if peel_dom >= 0 && auth_dom >= 0 {
            overlap[peel_dom as usize][auth_dom as usize] += 1;
        }
    }

    // Find optimal assignment using Hungarian algorithm (maximize overlap)
    // Convert to cost minimization: cost = max_val - overlap
    let max_val = *overlap.iter().flat_map(|row| row.iter()).max().unwrap_or(&0);
    let cost: Vec<Vec<i32>> = overlap
        .iter()
        .map(|row| row.iter().map(|&v| max_val - v).collect())
        .collect();

    let assignment = hungarian_algorithm(&cost);

    // Sum up the overlaps for the optimal assignment
    let best_jones: i32 = assignment
        .iter()
        .enumerate()
        .map(|(i, &j)| overlap[i][j])
        .sum();

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

/// Hungarian (Kuhn-Munkres) algorithm for the assignment problem.
///
/// Given an n×n cost matrix, finds an assignment of rows to columns that
/// minimizes total cost. Returns a vector where result[i] = column assigned to row i.
///
/// Time complexity: O(n³)
fn hungarian_algorithm(cost: &[Vec<i32>]) -> Vec<usize> {
    let n = cost.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![0];
    }

    // u[i] and v[j] are potentials for rows and columns (1-indexed, 0 is dummy)
    let mut u = vec![0i64; n + 1];
    let mut v = vec![0i64; n + 1];
    // p[j] = row assigned to column j
    let mut p = vec![0usize; n + 1];
    // way[j] = previous column in alternating path
    let mut way = vec![0usize; n + 1];

    for i in 1..=n {
        // Start augmenting path from row i
        p[0] = i;
        let mut j0 = 0usize;
        let mut minv = vec![i64::MAX; n + 1];
        let mut used = vec![false; n + 1];

        loop {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = i64::MAX;
            let mut j1 = 0usize;

            for j in 1..=n {
                if used[j] {
                    continue;
                }
                let cur = cost[i0 - 1][j - 1] as i64 - u[i0] - v[j];
                if cur < minv[j] {
                    minv[j] = cur;
                    way[j] = j0;
                }
                if minv[j] < delta {
                    delta = minv[j];
                    j1 = j;
                }
            }

            for j in 0..=n {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }

            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }

        // Update assignment along the alternating path
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }

    // Convert from column→row to row→column
    let mut result = vec![0usize; n];
    for j in 1..=n {
        if p[j] > 0 {
            result[p[j] - 1] = j - 1;
        }
    }
    result
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
pub fn read_ca_residue_numbers(pdb_path: &Path) -> Vec<i32> {
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
        if line.len() < 54 {
            continue;
        }
        let atom_name = line[12..16].trim();
        if atom_name != "CA" {
            continue;
        }
        let aa = line[17..20].trim();
        if !is_standard_aa(aa) {
            continue;
        }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hungarian_identity() {
        // Cost matrix where optimal is diagonal assignment
        let cost = vec![
            vec![0, 10, 10],
            vec![10, 0, 10],
            vec![10, 10, 0],
        ];
        let result = hungarian_algorithm(&cost);
        assert_eq!(result, vec![0, 1, 2]);
    }

    #[test]
    fn test_hungarian_swap() {
        // Cost matrix where optimal swaps 0 and 1
        let cost = vec![
            vec![10, 0],
            vec![0, 10],
        ];
        let result = hungarian_algorithm(&cost);
        assert_eq!(result, vec![1, 0]);
    }

    #[test]
    fn test_hungarian_3x3() {
        let cost = vec![
            vec![1, 2, 3],
            vec![2, 4, 6],
            vec![3, 6, 9],
        ];
        let result = hungarian_algorithm(&cost);
        // Total cost should be minimized
        let total: i32 = result.iter().enumerate().map(|(i, &j)| cost[i][j]).sum();
        // Optimal: 0→2(3), 1→1(4), 2→0(3) = 10
        assert_eq!(total, 10);
    }

    #[test]
    fn test_hungarian_single() {
        let cost = vec![vec![5]];
        let result = hungarian_algorithm(&cost);
        assert_eq!(result, vec![0]);
    }

    #[test]
    fn test_parse_range() {
        assert_eq!(parse_range("10-50"), Some((10, 50)));
        assert_eq!(parse_range("abc"), None);
    }
}
