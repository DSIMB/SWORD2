//! Junction/hinge consistency calculation.
//!
//! Port of `stat_pu_domains_from_SWORD.pl`.
//! Computes how consistently each residue position appears as a domain
//! boundary (junction) across the set of SWORD partitionings.

use std::collections::BTreeMap;

/// Compute junction consistency statistics from SWORD output lines.
///
/// Each SWORD output line has the format (pipe-delimited):
/// ```text
///   ndomains | min | boundaries | avg_kappa | quality |
/// ```
///
/// Returns a formatted string matching the original Perl script output.
pub fn calculate_junction_consistencies(sword_lines: &[String]) -> String {
    let mut hash_junctions: BTreeMap<i32, i32> = BTreeMap::new();
    let mut hash_weighted: BTreeMap<i32, i32> = BTreeMap::new();
    let mut total: i32 = 0;

    for line in sword_lines {
        let parts: Vec<&str> = line.split('|').collect();
        if parts.len() < 5 {
            continue;
        }

        // parts[0] = ndomains, parts[1] = min, parts[2] = boundaries,
        // parts[3] = avg_kappa, parts[4] = quality
        let ndom_str = parts[0].trim();
        if ndom_str.parse::<usize>().is_err() {
            continue;
        }

        let quality = parts[4].trim();
        // Skip lines where quality contains 'n' (e.g., "n/a")
        if quality.contains('n') {
            continue;
        }
        let cnt = quality.matches('*').count() as i32;
        total += 1;

        let delineation = parts[2].trim();
        let tab_domains: Vec<&str> = delineation.split_whitespace().collect();

        for domain in &tab_domains {
            // Extract start junction: number at the beginning before '-'
            if let Some(pos) = domain.find('-') {
                if let Ok(junction) = domain[..pos].parse::<i32>() {
                    *hash_junctions.entry(junction).or_insert(0) += 1;
                    *hash_weighted.entry(junction).or_insert(0) += cnt;
                }
            }
            // Extract end junction: number at the end after last '-'
            if let Some(pos) = domain.rfind('-') {
                if let Ok(junction) = domain[pos + 1..].parse::<i32>() {
                    *hash_junctions.entry(junction).or_insert(0) += 1;
                    *hash_weighted.entry(junction).or_insert(0) += cnt;
                }
            }
        }
    }

    let mut output = String::new();
    output.push_str("#Hinge/Junction consistency:\n");
    output.push_str("#Jnct  Cnt  Raw  Wei\n");

    if total > 0 {
        for (&junction, &count) in &hash_junctions {
            let weighted = hash_weighted.get(&junction).copied().unwrap_or(0);
            output.push_str(&format!(
                "{:<6} {:<4} {:4.2} {:4.2}\n",
                junction,
                count,
                count as f64 / total as f64,
                weighted as f64 / total as f64
            ));
        }
    }

    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_junction() {
        let lines = vec![
            "  2 | 30 | 1-100 101-200 | 3.5 | ***** |".to_string(),
            "  3 | 20 | 1-50 51-150 151-200 | 2.1 | *** |".to_string(),
        ];
        let result = calculate_junction_consistencies(&lines);
        assert!(result.contains("#Hinge/Junction consistency:"));
        // Should contain junction residue numbers from the boundaries
        assert!(result.contains("100"), "Expected junction at 100, got:\n{}", result);
        assert!(result.contains("200"), "Expected junction at 200, got:\n{}", result);
    }

    #[test]
    fn test_skip_na_quality() {
        let lines = vec![
            "  2 | 30 | 1-100 101-200 | 3.5 | n/a |".to_string(),
        ];
        let result = calculate_junction_consistencies(&lines);
        // Should only have header lines, no junction data
        assert!(!result.contains("100"));
    }
}
