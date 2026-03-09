//! Hierarchical PU merging / domain reconstruction algorithm.
//!
//! Port of `ComputeMeasure.pl` — the core algorithm that takes PU delineations
//! and a contact probability matrix and iteratively merges PUs into domains,
//! computing quality measures at each level.

use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

/// A single measure line from the hierarchical merging process.
#[derive(Debug, Clone)]
pub struct MeasureLine {
    /// Number of domains at this level.
    pub num_domains: usize,
    /// Minimum domain size.
    pub min_size: usize,
    /// Domain delineation string (space-separated, semicolons for discontinuous).
    pub delineation: String,
    /// Maximum contact ratio.
    pub max_cr: f64,
    /// Mean contact ratio (always 0 in current implementation).
    pub mean_cr: f64,
    /// Minimum density.
    pub density_min: f64,
    /// Mean density.
    pub mean_density: f64,
}

impl MeasureLine {
    /// Format as the pipe-delimited string that downstream code parses.
    pub fn to_line(&self) -> String {
        format!(
            "{:<2}|{:<3}|{:>50}|{:>10.6}|{:>10.6}|{:>10.6}|{:>10.6}",
            self.num_domains,
            self.min_size,
            self.delineation,
            self.max_cr,
            self.mean_cr,
            self.density_min,
            self.mean_density,
        )
    }
}

/// PU information.
#[derive(Debug, Clone)]
struct PuInfo {
    size: i32,
}

/// Run the full ComputeMeasure algorithm.
///
/// Reads the contact matrix and PU delineation files produced by Peeling_omp,
/// then hierarchically merges PUs into domains.
///
/// Returns a vector of `MeasureLine` items — one per domain decomposition
/// explored during the merging process.
pub fn compute_measure(
    file_contact: &Path,
    file_pu: &Path,
    _cutoff_pdp: f64,
) -> Vec<MeasureLine> {
    let max_number_results: usize = 500;
    let cutoff_size_domain: usize = 30;

    // 1) Read PU delineation
    let pu_content = fs::read_to_string(file_pu).unwrap_or_default();
    let mut hash_pu: BTreeMap<i32, (i32, i32, i32)> = BTreeMap::new(); // start -> (id, start, end)

    for line in pu_content.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            continue;
        }
        let id_pu: i32 = fields[0].parse().unwrap_or(0);
        let start: i32 = fields[1].parse().unwrap_or(0);
        let end: i32 = fields[2].parse().unwrap_or(0);
        hash_pu.insert(start, (id_pu, start, end));
    }

    if hash_pu.is_empty() {
        return Vec::new();
    }

    // Sort PUs by start position and build index mapping
    let mut pu_list: Vec<PuInfo> = Vec::new();
    let mut id_to_idx: BTreeMap<i32, usize> = BTreeMap::new();
    let mut pu_start_end: BTreeMap<usize, (i32, i32)> = BTreeMap::new();
    let mut total_size: usize = 0;

    for (pu_idx, (_start, (id_pu, s, e))) in hash_pu.iter().enumerate() {
        let size = e - s + 1;
        total_size += size as usize;
        id_to_idx.insert(*id_pu, pu_idx);
        pu_start_end.insert(pu_idx + 1, (*s, *e));
        pu_list.push(PuInfo {
            size,
        });
    }

    // 2) Read contact matrix
    let contact_content = fs::read_to_string(file_contact).unwrap_or_default();
    let mut tab_matrix: Vec<Vec<f64>> = vec![vec![0.0; pu_list.len()]; pu_list.len()];

    for line in contact_content.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            continue;
        }
        let x: i32 = fields[0].parse().unwrap_or(-1);
        let y: i32 = fields[1].parse().unwrap_or(-1);
        let value: f64 = fields[2].parse().unwrap_or(0.0);

        if let (Some(&a), Some(&b)) = (id_to_idx.get(&x), id_to_idx.get(&y)) {
            tab_matrix[a][b] = value;
        }
    }

    let pu_sizes: Vec<f64> = pu_list.iter().map(|p| p.size as f64).collect();

    // 3) Initial domain assignment: each PU is its own domain
    let initial_domains: Vec<String> = (1..=pu_list.len())
        .map(|i| i.to_string())
        .collect();

    // Compute initial measure
    let mut output_lines: Vec<MeasureLine> = Vec::new();

    let (min_size_1, max_cr_1, _mean_cr_1, density_min_1, mean_density_1) =
        measure_domain(&initial_domains, &tab_matrix, &pu_sizes);
    let delineation_1 = print_domain(&initial_domains, &pu_start_end);
    output_lines.push(MeasureLine {
        num_domains: initial_domains.len(),
        min_size: min_size_1,
        delineation: delineation_1,
        max_cr: max_cr_1,
        mean_cr: 0.0,
        density_min: density_min_1,
        mean_density: mean_density_1,
    });

    // 4) Iterative merging
    let mut all_tab_domains: Vec<Vec<String>> = vec![initial_domains];
    let mut number_domains = pu_list.len();

    while number_domains > 0 {
        // Remove duplicates
        let mut unique_domains: Vec<Vec<String>> = Vec::new();
        let mut seen: HashSet<String> = HashSet::new();
        for doms in &all_tab_domains {
            let key = doms.join(" ");
            if seen.insert(key) {
                unique_domains.push(doms.clone());
            }
        }
        all_tab_domains = unique_domains;

        // Try merging PU pairs
        let mut all_results: Vec<MergeResult> = Vec::new();
        for doms in &all_tab_domains {
            let mut results = compute_merge_pu(doms, &tab_matrix, &pu_sizes);
            all_results.append(&mut results);
        }

        if all_results.is_empty() {
            break;
        }

        // Sort fragments within each result
        for result in &mut all_results {
            sort_domain_fragments(&mut result.new_domains);
        }

        // Sort results by ratio_pdp (descending) and take top N
        all_results.sort_by(|a, b| b.ratio_pdp.partial_cmp(&a.ratio_pdp).unwrap_or(std::cmp::Ordering::Equal));
        all_results.truncate(max_number_results + 1);

        // Compute measures and produce output for this merging level
        let mut new_all_tab_domains: Vec<Vec<String>> = Vec::new();
        let mut seen_delineations: HashSet<String> = HashSet::new();
        let mut level_lines: Vec<(f64, MeasureLine, Vec<String>)> = Vec::new();

        for result in &all_results {
            let print_dom = merge_dom(&print_domain(&result.new_domains, &pu_start_end));

            if seen_delineations.contains(&print_dom) {
                continue;
            }
            seen_delineations.insert(print_dom.clone());

            let (min_sz, max_cr, _mean_cr, den_min, mean_den) =
                measure_domain(&result.new_domains, &tab_matrix, &pu_sizes);

            level_lines.push((mean_den, MeasureLine {
                num_domains: number_domains - 1,
                min_size: min_sz,
                delineation: print_dom,
                max_cr,
                mean_cr: 0.0,
                density_min: den_min,
                mean_density: mean_den,
            }, result.new_domains.clone()));
        }

        // Sort by mean_density (ascending)
        level_lines.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

        let limit_number_of_domains = (total_size + cutoff_size_domain - 1) / cutoff_size_domain;
        let mut count = 0;

        for (_, ml, doms) in &level_lines {
            new_all_tab_domains.push(doms.clone());
            if number_domains <= limit_number_of_domains {
                output_lines.push(ml.clone());
                count += 1;
                if count > max_number_results {
                    break;
                }
            }
        }

        number_domains = number_domains.saturating_sub(1);
        all_tab_domains = new_all_tab_domains;

        if number_domains == 0 || all_tab_domains.is_empty() {
            break;
        }
        // Safety: also stop if first domain set is a single domain
        if all_tab_domains.iter().all(|d| d.len() <= 1) {
            break;
        }
    }

    output_lines
}

/// Result of attempting to merge two PUs.
#[derive(Debug, Clone)]
struct MergeResult {
    new_domains: Vec<String>,
    ratio_pdp: f64,
}

/// Try all pair-wise merges of domains and return candidate results.
fn compute_merge_pu(
    domains: &[String],
    tab_matrix: &[Vec<f64>],
    pu_sizes: &[f64],
) -> Vec<MergeResult> {
    let mut results = Vec::new();

    for j in 0..domains.len() {
        let sub_doms_1 = parse_sub_domains(&domains[j]);
        let (size_int_1, contact_int_1) = measure_internal(&sub_doms_1, tab_matrix, pu_sizes);

        for k in (j + 1)..domains.len() {
            let sub_doms_2 = parse_sub_domains(&domains[k]);
            let (size_int_2, contact_int_2) = measure_internal(&sub_doms_2, tab_matrix, pu_sizes);
            let (_size_ext, contact_ext) =
                measure_external(&sub_doms_1, &sub_doms_2, tab_matrix, pu_sizes);

            let contact_ext_safe = if contact_ext <= 0.001 { 0.001 } else { contact_ext };
            let nnc = contact_ext_safe
                / (size_int_1.powf(0.43) * size_int_2.powf(0.43));

            let denom_int = (contact_int_1 + contact_int_2) / (size_int_1 * size_int_2);
            let ratio_pdp = if denom_int.abs() < 1e-15 { 0.0 } else { nnc / denom_int };

            if results.is_empty() || nnc > 0.01 {
                let new_pu = format!("{};{}", domains[j], domains[k]);
                let mut new_domains = domains.to_vec();
                new_domains[j] = new_pu;
                new_domains.remove(k);

                results.push(MergeResult {
                    new_domains,
                    ratio_pdp,
                });
            }
        }
    }

    results
}

/// Measure internal contacts and total size for a (possibly multi-PU) domain.
fn measure_internal(
    sub_domains: &[usize],
    tab_matrix: &[Vec<f64>],
    pu_sizes: &[f64],
) -> (f64, f64) {
    let mut contact = 0.0;
    let mut size = 0.0;

    for &l in sub_domains {
        for &m in sub_domains {
            if l < tab_matrix.len() && m < tab_matrix.len() {
                contact += tab_matrix[l][m];
            }
            if l == m && l < pu_sizes.len() {
                size += pu_sizes[l];
            }
        }
    }
    contact /= 2.0;
    (size, contact)
}

/// Measure external contacts between two domains.
fn measure_external(
    sub_doms_1: &[usize],
    sub_doms_2: &[usize],
    tab_matrix: &[Vec<f64>],
    _pu_sizes: &[f64],
) -> (f64, f64) {
    let mut contact = 0.0;

    for &l in sub_doms_1 {
        for &m in sub_doms_2 {
            if l < tab_matrix.len() && m < tab_matrix.len() {
                contact += tab_matrix[l][m];
            }
        }
    }
    (0.0, contact)
}

/// Measure domain partition quality: min size, max CR, mean CR, min density, mean density.
fn measure_domain(
    domains: &[String],
    tab_matrix: &[Vec<f64>],
    pu_sizes: &[f64],
) -> (usize, f64, f64, f64, f64) {
    let n = domains.len();
    if n <= 1 {
        let sub = parse_sub_domains(&domains[0]);
        let (sz, _) = measure_internal(&sub, tab_matrix, pu_sizes);
        return (sz as usize, 0.0, 0.0, 0.0, 0.0);
    }

    let mut min_size = usize::MAX;
    let mut max_cr: f64 = 0.0;
    let mut min_density: f64 = f64::MAX;
    let mut density_tot: f64 = 0.0;

    for j in 0..n {
        let sub1 = parse_sub_domains(&domains[j]);
        let (size1, contact1) = measure_internal(&sub1, tab_matrix, pu_sizes);
        if (size1 as usize) < min_size {
            min_size = size1 as usize;
        }
        let density1 = if size1 > 0.0 { contact1 / size1 } else { 0.0 };
        density_tot += density1;

        for k in (j + 1)..n {
            let sub2 = parse_sub_domains(&domains[k]);
            let (size2, contact2) = measure_internal(&sub2, tab_matrix, pu_sizes);
            let (_, contact_ext) = measure_external(&sub1, &sub2, tab_matrix, pu_sizes);

            if (size2 as usize) < min_size {
                min_size = size2 as usize;
            }

            let density2 = if size2 > 0.0 { contact2 / size2 } else { 0.0 };
            if density1 < min_density {
                min_density = density1;
            }
            if density2 < min_density {
                min_density = density2;
            }

            let current_nnc = if size1 > 0.0 && size2 > 0.0 {
                contact_ext / (size1.powf(0.43) * size2.powf(0.43))
            } else {
                0.0
            };

            let total_contact = contact_ext + contact1 + contact2;
            let total_size = size1 + size2;
            let criterion = if total_size > 0.0 { total_contact / total_size } else { 1.0 };
            let current_cr = if criterion > 0.0 { current_nnc / criterion } else { 0.0 };

            if current_cr > max_cr {
                max_cr = current_cr;
            }
        }
    }

    let mean_density = density_tot / n as f64;
    // mean_cr is always set to 0 in the original Perl code
    (min_size, max_cr, 0.0, min_density, mean_density)
}

/// Parse a domain string like "1;3" into PU indices (0-based).
fn parse_sub_domains(dom: &str) -> Vec<usize> {
    dom.split(';')
        .filter_map(|s| s.trim().parse::<usize>().ok().map(|v| v.saturating_sub(1)))
        .collect()
}

/// Sort domain fragments lexically by first PU number.
fn sort_domain_fragments(domains: &mut Vec<String>) {
    for dom in domains.iter_mut() {
        if dom.contains(';') {
            let mut parts: Vec<usize> = dom
                .split(';')
                .filter_map(|s| s.parse::<usize>().ok())
                .collect();
            parts.sort();
            *dom = parts.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(";");
        }
    }
    domains.sort_by(|a, b| {
        let first_a = a.split(';').next().and_then(|s| s.parse::<usize>().ok()).unwrap_or(0);
        let first_b = b.split(';').next().and_then(|s| s.parse::<usize>().ok()).unwrap_or(0);
        first_a.cmp(&first_b)
    });
}

/// Convert domain PU-id list to boundary delineation string.
fn print_domain(domains: &[String], pu_start_end: &BTreeMap<usize, (i32, i32)>) -> String {
    let mut parts = Vec::new();

    for dom in domains {
        let sub_pus: Vec<usize> = dom
            .split(';')
            .filter_map(|s| s.parse::<usize>().ok())
            .collect();

        let mut delineation = String::new();
        let mut old_pu = 0usize;

        for (bi, &pu_id) in sub_pus.iter().enumerate() {
            if let Some(&(s, e)) = pu_start_end.get(&pu_id) {
                if bi > 0 && pu_id == old_pu + 1 {
                    // Consecutive PU: extend previous range
                    // Replace the last end number
                    if let Some(pos) = delineation.rfind('-') {
                        // Find the last end value and replace it
                        let new_delineation = format!("{}{};", &delineation[..pos + 1], e);
                        delineation = new_delineation;
                    }
                } else {
                    delineation.push_str(&format!("{}-{};", s, e));
                }
                old_pu = pu_id;
            }
        }

        // Remove trailing semicolon
        if delineation.ends_with(';') {
            delineation.pop();
        }
        parts.push(delineation);
    }

    parts.join(" ")
}

/// Merge adjacent PU segments in a domain delineation string.
///
/// E.g., "1-50;51-100" → "1-100" when the segments are contiguous.
fn merge_dom(dom_peel: &str) -> String {
    let tab_dom_peel: Vec<&str> = dom_peel.split_whitespace().collect();
    let mut result_parts = Vec::new();

    for dom in &tab_dom_peel {
        if dom.contains(';') {
            let mut segments: Vec<(i32, i32)> = dom
                .split(';')
                .filter(|s| !s.is_empty())
                .filter_map(|seg| {
                    let parts: Vec<&str> = seg.split('-').collect();
                    if parts.len() == 2 {
                        Some((
                            parts[0].parse::<i32>().unwrap_or(0),
                            parts[1].parse::<i32>().unwrap_or(0),
                        ))
                    } else {
                        None
                    }
                })
                .collect();

            // Merge contiguous segments
            let mut merged = true;
            while merged {
                merged = false;
                for i in 0..segments.len() {
                    for j in (i + 1)..segments.len() {
                        if segments[i].1 + 1 == segments[j].0 {
                            segments[i].1 = segments[j].1;
                            segments.remove(j);
                            merged = true;
                            break;
                        } else if segments[j].1 + 1 == segments[i].0 {
                            segments[i].0 = segments[j].0;
                            segments.remove(j);
                            merged = true;
                            break;
                        }
                    }
                    if merged {
                        break;
                    }
                }
            }

            let seg_str: Vec<String> = segments
                .iter()
                .map(|(s, e)| format!("{}-{}", s, e))
                .collect();
            result_parts.push(seg_str.join(";"));
        } else {
            result_parts.push(dom.to_string());
        }
    }

    result_parts.join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_dom_contiguous() {
        assert_eq!(merge_dom("1-50;51-100"), "1-100");
    }

    #[test]
    fn test_merge_dom_non_contiguous() {
        assert_eq!(merge_dom("1-50;80-100"), "1-50;80-100");
    }

    #[test]
    fn test_parse_sub_domains() {
        assert_eq!(parse_sub_domains("1;3;5"), vec![0, 2, 4]);
        assert_eq!(parse_sub_domains("2"), vec![1]);
    }
}
