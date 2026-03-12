//! Domain assignment selection and filtering.
//!
//! Port of `ParseMeasure.pm`. Filters the output of ComputeMeasure to select
//! relevant domain assignments with alternative positions and levels around
//! the predicted optimal number of domains.

use super::compute_jones;
use super::distance_model;

/// Filter measure lines to select relevant domain assignments.
///
/// `measure_lines` - output from ComputeMeasure (pipe-separated strings).
/// `dir_data` - path to PDBs_Clean/ directory.
/// `pdb` - PDB identifier (will be suffixed with "_1_").
/// `option_alt_dist` - whether to filter by DistanceModel distance.
/// `n_dom` - predicted optimal number of domains (0 on first pass).
/// `alt_b` - number of boundary alternatives.
/// `alt_l` - number of level alternatives.
/// `option_alt_diff` - whether to check Jones overlap for uniqueness.
///
/// Returns the filtered measure lines.
pub fn parse_measure(
    measure_lines: &[String],
    dir_data: &str,
    pdb: &str,
    option_alt_dist: bool,
    n_dom: usize,
    alt_b: usize,
    alt_l: usize,
    option_alt_diff: bool,
) -> Vec<String> {
    let pdb_id = format!("{}_1_", pdb);
    let max_dist: f64 = 0.2;

    if option_alt_diff {
        parse_measure_with_diff(
            measure_lines,
            dir_data,
            &pdb_id,
            option_alt_dist,
            n_dom,
            alt_b,
            alt_l,
            max_dist,
        )
    } else {
        parse_measure_simple(
            measure_lines,
            dir_data,
            option_alt_dist,
            n_dom,
            alt_b,
            alt_l,
            max_dist,
        )
    }
}

fn get_num_domains(line: &str) -> usize {
    line.chars()
        .take(2)
        .collect::<String>()
        .trim()
        .parse::<usize>()
        .unwrap_or(0)
}

fn is_comment(line: &str) -> bool {
    line.starts_with('#')
}

fn get_fields(line: &str) -> Vec<String> {
    line.split('|').map(|s| s.trim().to_string()).collect()
}

/// Parse measure with diff (option_alt_diff == true).
fn parse_measure_with_diff(
    measure_lines: &[String],
    dir_data: &str,
    pdb_id: &str,
    option_alt_dist: bool,
    n_dom: usize,
    alt_b: usize,
    alt_l: usize,
    max_dist: f64,
) -> Vec<String> {
    let mut clean_measure: Vec<String> = Vec::new();
    let mut temp_measure: Vec<String> = Vec::new();
    let mut max_dom: usize = 0;

    // Pre-load CA residue numbers once for all Jones overlap checks
    let pdb_base = &pdb_id[..pdb_id.len().saturating_sub(3)];
    let pdb_file_path = std::path::Path::new(dir_data)
        .join(pdb_base)
        .join(format!("{}.pdb", pdb_base));
    let cached_residue_nums = compute_jones::read_ca_residue_numbers(&pdb_file_path);

    // First pass: group and filter by domain count levels
    for line in measure_lines {
        if is_comment(line) {
            continue;
        }

        let nd = get_num_domains(line);

        if max_dom == 0 {
            if nd > (n_dom + alt_l) && nd > 6 {
                clean_measure.push(line.clone());
            } else {
                max_dom = nd;
            }
        }

        if nd == max_dom {
            if option_alt_dist {
                let fields = get_fields(line);
                if fields.len() > 5 {
                    let cr: f64 = fields[3].parse().unwrap_or(0.0);
                    let cpd: f64 = fields[5].parse().unwrap_or(0.0);
                    // EXPERIMENT: signed like Perl
                    if distance_model::distance_model(cr, cpd, 1) < max_dist {
                        temp_measure.push(line.clone());
                    }
                }
            } else {
                temp_measure.push(line.clone());
            }
        }

        if nd < max_dom {
            max_dom = nd;
            temp_measure.reverse();

            if temp_measure.len() > 1 && alt_b > 1 {
                let mut temp_measure2: Vec<String> = Vec::new();
                let mut i = 0;
                while i < alt_b && i <= temp_measure.len().saturating_sub(1) {
                    temp_measure2.insert(0, temp_measure[i].clone());

                    if i == temp_measure.len() - 1 || i == alt_b - 1 {
                        break;
                    }

                    // Check Jones overlap with subsequent entries
                    let mut cpt = 0;
                    let mut j = i + 1;
                    while j < temp_measure.len() {
                        cpt += 1;
                        if cpt > 20 {
                            break;
                        }

                        let fields_i = get_fields(&temp_measure[i]);
                        let fields_j = get_fields(&temp_measure[j]);

                        if fields_i.len() > 2 && fields_j.len() > 2 {
                            let del1 = fields_i[2].replace(' ', "_").replace(';', "\\;");
                            let del2 = fields_j[2].replace(' ', "_").replace(';', "\\;");
                            let del_str1 = format!("{}{}", pdb_id, del1);
                            let del_str2 = format!("{}{}", pdb_id, del2);

                            let (criterion, pct) = compute_jones::compute_jones_with_cache(
                                pdb_id,
                                &del_str1,
                                &del_str2,
                                dir_data,
                                Some(&cached_residue_nums),
                            );
                            tracing::info!("  Jones[i={},j={}] criterion={} pct={:.1}% del_i={} del_j={}",
                                i, j, criterion, pct,
                                fields_i[2].trim(), fields_j[2].trim());
                            if criterion == 1 {
                                // Too similar — remove
                                tracing::info!("    REMOVED j={}", j);
                                temp_measure.remove(j);
                                continue; // don't increment j
                            }
                        }

                        if i == alt_b.saturating_sub(2) {
                            break;
                        }
                        j += 1;
                    }
                    i += 1;
                }
                clean_measure.extend(temp_measure2.clone());
                
            } else if !temp_measure.is_empty() {
                clean_measure.push(temp_measure[0].clone());
            }

            temp_measure.clear();

            // Check if current line satisfies distance criteria
            let fields = get_fields(line);
            if fields.len() > 5 {
                let cr: f64 = fields[3].parse().unwrap_or(0.0);
                let cpd: f64 = fields[5].parse().unwrap_or(0.0);
                // EXPERIMENT: signed like Perl
                if distance_model::distance_model(cr, cpd, 1) < max_dist {
                    temp_measure.push(line.clone());
                }
            }
        }
    }

    // Add last measure line
    if let Some(last) = measure_lines.last() {
        clean_measure.push(last.clone());
    }

    if alt_b == 1 {
        return clean_measure;
    }

    // Second pass: select alt_b alternatives per level
    let mut relevant_measure: Vec<String> = Vec::new();
    let mut max_dom2: usize = 0;

    for (id_line, line) in clean_measure.iter().enumerate() {
        let nd = get_num_domains(line);
        if max_dom2 == 0 {
            max_dom2 = nd;
        }

        if nd < max_dom2 {
            for i in 1..=alt_b {
                if id_line >= i {
                    let prev_nd = get_num_domains(&clean_measure[id_line - i]);
                    if prev_nd == max_dom2 {
                        relevant_measure.push(clean_measure[id_line - i].clone());
                    }
                }
            }
            max_dom2 = nd;
        }
    }

    // Add last line
    if let Some(last) = measure_lines.last() {
        relevant_measure.push(last.clone());
    }

    relevant_measure
}

/// Simpler parse without diff checking.
fn parse_measure_simple(
    measure_lines: &[String],
    _dir_data: &str,
    option_alt_dist: bool,
    _n_dom: usize,
    alt_b: usize,
    _alt_l: usize,
    max_dist: f64,
) -> Vec<String> {
    let mut relevant_measure: Vec<String> = Vec::new();
    let mut max_dom: usize = 0;

    for (id_line, line) in measure_lines.iter().enumerate() {
        if is_comment(line) {
            continue;
        }

        let nd = get_num_domains(line);
        if max_dom == 0 {
            max_dom = nd;
        }

        if nd < max_dom {
            for i in 1..=(alt_b + 1) {
                if id_line >= i {
                    let prev_nd = get_num_domains(&measure_lines[id_line - i]);
                    if prev_nd == max_dom {
                        if option_alt_dist {
                            let fields = get_fields(&measure_lines[id_line - i]);
                            if fields.len() > 5 {
                                let cr: f64 = fields[3].parse().unwrap_or(0.0);
                                let cpd: f64 = fields[5].parse().unwrap_or(0.0);
                                if distance_model::distance_model(cr, cpd, 1).abs() < max_dist {
                                    relevant_measure
                                        .push(measure_lines[id_line - i].clone());
                                }
                            }
                        } else {
                            relevant_measure.push(measure_lines[id_line - i].clone());
                        }
                    }
                }
            }
            max_dom = nd;
        }
    }

    if let Some(last) = measure_lines.last() {
        relevant_measure.push(last.clone());
    }

    relevant_measure
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_num_domains() {
        assert_eq!(get_num_domains(" 5|30|..."), 5);
        assert_eq!(get_num_domains("12|30|..."), 12);
    }
}
