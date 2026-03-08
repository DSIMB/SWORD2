//! Visualization and plotting utilities.
//!
//! Generates contact probability matrix plots and domain consistency histograms
//! as SVG files using the `plotters` crate.

use std::path::Path;

use anyhow::{Context, Result};
use ndarray::Array2;
use plotters::prelude::*;

/// Color palette for Protein Units (pastel colors).
const PU_COLORS: &[(u8, u8, u8)] = &[
    (186, 234, 229), // #baeae5
    (225, 198, 91),  // #e1c65b
    (180, 188, 247), // #b4bcf7
    (208, 228, 123), // #d0e47b
    (240, 168, 229), // #f0a8e5
    (109, 228, 172), // #6de4ac
    (216, 192, 228), // #d8c0e4
    (165, 225, 141), // #a5e18d
    (104, 209, 241), // #68d1f1
    (243, 177, 117), // #f3b175
    (99, 227, 216),  // #63e3d8
    (235, 186, 186), // #ebbaba
    (195, 210, 140), // #c3d28c
    (170, 197, 226), // #aac5e2
    (232, 218, 146), // #e8da92
    (188, 219, 236), // #bcdbec
    (225, 194, 152), // #e1c298
    (152, 199, 198), // #98c7c6
    (171, 221, 180), // #abddb4
    (212, 216, 187), // #d4d8bb
];

/// Color palette for domains (saturated colors).
const DOMAIN_COLORS: &[(u8, u8, u8)] = &[
    (39, 163, 180),  // #27a3b4
    (192, 132, 35),  // #c08423
    (216, 62, 124),  // #d83e7c
    (152, 106, 53),  // #986a35
    (104, 111, 223), // #686fdf
    (85, 154, 59),   // #559a3b
    (118, 61, 166),  // #763da6
    (141, 141, 54),  // #8d8d36
    (206, 97, 199),  // #ce61c7
    (64, 96, 33),    // #406021
    (164, 44, 136),  // #a42c88
    (61, 149, 107),  // #3d956b
    (52, 29, 121),   // #341d79
    (207, 58, 68),   // #cf3a44
    (60, 140, 201),  // #3c8cc9
    (207, 100, 48),  // #cf6430
    (77, 75, 146),   // #4d4b92
    (125, 49, 25),   // #7d3119
    (123, 129, 205), // #7b81cd
    (207, 108, 97),  // #cf6c61
    (64, 29, 86),    // #401d56
    (201, 95, 122),  // #c95f7a
    (121, 46, 101),  // #792e65
    (130, 38, 58),   // #82263a
    (184, 109, 168), // #b86da8
];

/// Load a contact probability matrix from a whitespace-delimited text file.
///
/// Each line is a row of space-separated float values.
pub fn load_contact_matrix(path: &Path) -> Result<Array2<f64>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Cannot read matrix file: {}", path.display()))?;

    let rows: Vec<Vec<f64>> = content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            line.split_whitespace()
                .filter_map(|s| s.parse::<f64>().ok())
                .collect()
        })
        .collect();

    if rows.is_empty() {
        anyhow::bail!("Empty matrix file: {}", path.display());
    }

    let n_rows = rows.len();
    let n_cols = rows[0].len();

    let flat: Vec<f64> = rows.into_iter().flatten().collect();
    Array2::from_shape_vec((n_rows, n_cols), flat)
        .with_context(|| "Failed to construct matrix from file data")
}

/// Generate a contact probability matrix SVG plot.
///
/// Draws the matrix as a heatmap with PU boundary rectangles overlaid.
///
/// # Arguments
/// * `matrix` - The contact probability matrix (NxN)
/// * `output_path` - Path to write the output SVG file
/// * `title` - Plot title
/// * `pu_boundaries` - List of (start, end) PU boundary pairs to draw rectangles for
pub fn write_contact_matrix(
    matrix: &Array2<f64>,
    output_path: &str,
    title: &str,
    pu_boundaries: &[(i32, i32)],
) -> Result<()> {
    let (n_rows, n_cols) = matrix.dim();
    let width = 600u32;
    let height = 700u32;

    let root = SVGBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 14))
        .margin(10)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0i32..n_cols as i32, 0i32..n_rows as i32)?;

    chart
        .configure_mesh()
        .x_desc("Residues")
        .y_desc("Residues")
        .draw()?;

    // Find max value for color scaling
    let max_val = matrix.iter().cloned().fold(0.0f64, f64::max).max(1e-10);

    // Draw the heatmap as colored rectangles
    for i in 0..n_rows {
        for j in 0..n_cols {
            let val = matrix[[i, j]];
            let intensity = (val / max_val).min(1.0);
            // RdPu-like color map: white -> pink -> purple
            let r = (255.0 - intensity * 128.0) as u8;
            let g = (255.0 - intensity * 200.0) as u8;
            let b = (255.0 - intensity * 100.0) as u8;
            let color = RGBColor(r, g, b);

            chart.draw_series(std::iter::once(Rectangle::new(
                [(j as i32, (n_rows - 1 - i) as i32), (j as i32 + 1, (n_rows - i) as i32)],
                color.filled(),
            )))?;
        }
    }

    // Draw PU boundary rectangles
    for (idx, &(start, end)) in pu_boundaries.iter().enumerate() {
        let color_idx = idx % PU_COLORS.len();
        let (r, g, b) = PU_COLORS[color_idx];
        let color = RGBColor(r, g, b);
        let s = start - 1; // 0-based
        let e = end; // exclusive
        let y_s = (n_rows as i32) - e;
        let y_e = (n_rows as i32) - s;

        // Draw rectangle outline
        chart.draw_series(std::iter::once(Rectangle::new(
            [(s, y_s), (e, y_e)],
            ShapeStyle::from(color).stroke_width(2),
        )))?;
    }

    root.present()?;
    Ok(())
}

/// Domain count info for histogram.
#[derive(Debug, Clone)]
pub struct DomainCount {
    /// Domain label (boundary description).
    pub label: String,
    /// How many partitionings include this domain.
    pub count: usize,
    /// Color for this domain.
    pub color: (u8, u8, u8),
}

/// Generate a domain consistency histogram as SVG.
///
/// Shows how consistently each unique domain appears across
/// alternative partitionings.
///
/// Ported from Python `write_domains_histogram()`.
pub fn write_domain_histogram(domain_counts: &[DomainCount], output_path: &str) -> Result<()> {
    if domain_counts.is_empty() {
        return Ok(());
    }

    let width = 800u32;
    let height = 500u32;

    let root = SVGBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_count = domain_counts.iter().map(|d| d.count).max().unwrap_or(1);
    let n_bars = domain_counts.len();

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Consistency of domains determined by SWORD",
            ("sans-serif", 16),
        )
        .margin(20)
        .x_label_area_size(80)
        .y_label_area_size(40)
        .build_cartesian_2d(0i32..n_bars as i32, 0i32..(max_count as i32 + 1))?;

    chart
        .configure_mesh()
        .x_desc("SWORD Domains")
        .y_desc("Count")
        .x_labels(n_bars)
        .x_label_formatter(&|idx| {
            domain_counts
                .get(*idx as usize)
                .map(|d| d.label.clone())
                .unwrap_or_default()
        })
        .draw()?;

    // Draw bars
    for (i, dc) in domain_counts.iter().enumerate() {
        let (r, g, b) = dc.color;
        let color = RGBColor(r, g, b);
        chart.draw_series(std::iter::once(Rectangle::new(
            [(i as i32, 0), (i as i32 + 1, dc.count as i32)],
            color.filled(),
        )))?;
    }

    root.present()?;
    Ok(())
}

/// Count unique domains across all partitionings and assign colors.
///
/// Returns domain counts sorted by frequency (descending).
///
/// Ported from Python `write_domains_histogram()` + `define_colors()`.
pub fn count_domains(
    partitions: &[crate::sword::SwordPartition],
) -> Vec<DomainCount> {
    use std::collections::HashMap;

    let mut domain_freq: HashMap<Vec<(i32, i32)>, usize> = HashMap::new();

    for part in partitions {
        for domain in &part.boundaries {
            let mut sorted = domain.clone();
            sorted.sort_by_key(|&(s, _)| s);
            *domain_freq.entry(sorted).or_insert(0) += 1;
        }
    }

    // Sort by count descending
    let mut sorted: Vec<_> = domain_freq.into_iter().collect();
    sorted.sort_by(|a, b| b.1.cmp(&a.1));

    sorted
        .into_iter()
        .enumerate()
        .map(|(i, (boundaries, count))| {
            let label = boundaries
                .iter()
                .map(|(s, e)| format!("({}, {})", s, e))
                .collect::<Vec<_>>()
                .join(", ");
            let color_idx = i % DOMAIN_COLORS.len();
            DomainCount {
                label,
                count,
                color: DOMAIN_COLORS[color_idx],
            }
        })
        .collect()
}

/// Get the PU color for a given index.
pub fn get_pu_color(index: usize) -> (u8, u8, u8) {
    PU_COLORS[index % PU_COLORS.len()]
}

/// Get the domain color for a given index.
pub fn get_domain_color(index: usize) -> (u8, u8, u8) {
    DOMAIN_COLORS[index % DOMAIN_COLORS.len()]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sword::SwordPartition;

    #[test]
    fn test_count_domains() {
        let partitions = vec![
            SwordPartition {
                nb_domains: 2,
                min_size: 30,
                boundaries: vec![vec![(1, 100)], vec![(101, 200)]],
                average_k: 3.5,
                quality: "*****".to_string(),
            },
            SwordPartition {
                nb_domains: 2,
                min_size: 30,
                boundaries: vec![vec![(1, 100)], vec![(101, 200)]],
                average_k: 3.0,
                quality: "****".to_string(),
            },
            SwordPartition {
                nb_domains: 3,
                min_size: 20,
                boundaries: vec![vec![(1, 50)], vec![(51, 150)], vec![(151, 200)]],
                average_k: 2.0,
                quality: "***".to_string(),
            },
        ];

        let counts = count_domains(&partitions);
        // (1,100) and (101,200) appear in 2 partitions each
        assert!(counts[0].count >= 2);
    }
}
