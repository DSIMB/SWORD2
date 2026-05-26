//! Visualization and plotting utilities.
//!
//! Generates contact probability matrix plots (PNG) and domain consistency
//! histograms (SVG) matching the original Python SWORD2 output.
//!
//! Three levels of contact probability matrix plots are produced per alternative:
//! 1. Alternative-level: all PU rectangles
//! 2. Domain-level: one domain's PU rectangles
//! 3. PU-level: single PU rectangle

use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::path::Path;
use std::sync::Once;

use anyhow::{Context, Result};
use ndarray::Array2;
use plotters::prelude::*;

const EMBEDDED_SANS_FONT: &[u8] = include_bytes!("../../assets/fonts/Abel-Regular.ttf");
static REGISTER_EMBEDDED_FONT: Once = Once::new();

fn register_embedded_font() {
    REGISTER_EMBEDDED_FONT.call_once(|| {
        assert!(
            plotters::style::register_font("sans-serif", FontStyle::Normal, EMBEDDED_SANS_FONT)
                .is_ok(),
            "embedded plot font must be a valid TrueType font"
        );
    });
}

/// Color palette for Protein Units (pastel colors).
pub const PU_COLORS: &[(u8, u8, u8)] = &[
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
pub const DOMAIN_COLORS: &[(u8, u8, u8)] = &[
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

/// RdPu colormap approximation (11 stops interpolated from matplotlib).
const RDPU_STOPS: &[(f64, u8, u8, u8)] = &[
    (0.0, 255, 247, 243),
    (0.1, 253, 224, 221),
    (0.2, 252, 197, 192),
    (0.3, 250, 159, 181),
    (0.4, 247, 104, 161),
    (0.5, 221, 52, 151),
    (0.6, 174, 1, 126),
    (0.7, 122, 1, 119),
    (0.8, 73, 0, 106),
    (0.9, 47, 0, 79),
    (1.0, 47, 0, 79),
];

/// Interpolate a color from the RdPu colormap.
fn rdpu_color(t: f64) -> (u8, u8, u8) {
    let t = t.clamp(0.0, 1.0);
    for window in RDPU_STOPS.windows(2) {
        let (t0, r0, g0, b0) = window[0];
        let (t1, r1, g1, b1) = window[1];
        if t >= t0 && t <= t1 {
            let frac = if (t1 - t0).abs() < 1e-12 {
                0.0
            } else {
                (t - t0) / (t1 - t0)
            };
            let r = (r0 as f64 + frac * (r1 as f64 - r0 as f64)) as u8;
            let g = (g0 as f64 + frac * (g1 as f64 - g0 as f64)) as u8;
            let b = (b0 as f64 + frac * (b1 as f64 - b0 as f64)) as u8;
            return (r, g, b);
        }
    }
    let (_, r, g, b) = RDPU_STOPS[RDPU_STOPS.len() - 1];
    (r, g, b)
}

/// Load a contact probability matrix from a whitespace-delimited text file.
pub fn load_contact_matrix(path: &Path) -> Result<Array2<f64>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Cannot read matrix file: {}", path.display()))?;

    let rows: Vec<Vec<f64>> = content
        .lines()
        .filter_map(|line| {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                return None;
            }
            let row: Vec<f64> = trimmed
                .split_whitespace()
                .map(str::parse::<f64>)
                .collect::<std::result::Result<_, _>>()
                .ok()?;
            if row.is_empty() {
                None
            } else {
                Some(row)
            }
        })
        .collect();

    if rows.is_empty() {
        anyhow::bail!("Empty matrix file: {}", path.display());
    }

    let n_rows = rows.len();
    let n_cols = rows[0].len();

    if rows.iter().any(|row| row.len() != n_cols) {
        anyhow::bail!("Non-rectangular matrix file: {}", path.display());
    }

    let flat: Vec<f64> = rows.into_iter().flatten().collect();
    Array2::from_shape_vec((n_rows, n_cols), flat)
        .with_context(|| "Failed to construct matrix from file data")
}

/// A PU boundary with its assigned color.
#[derive(Debug, Clone)]
pub struct ColoredPu {
    pub start: i32,
    pub end: i32,
    pub color: (u8, u8, u8),
}

/// Assign unique colors to PUs across all partitionings (matching original Python logic).
///
/// Returns a map from (start, end) to color.
pub fn assign_pu_colors(
    partitions: &[crate::sword::SwordPartition],
) -> HashMap<(i32, i32), (u8, u8, u8)> {
    let mut pu_colors = HashMap::new();
    let mut color_idx = 0usize;

    for part in partitions {
        for domain in &part.boundaries {
            for &(start, end) in domain {
                if let Entry::Vacant(e) = pu_colors.entry((start, end)) {
                    e.insert(PU_COLORS[color_idx % PU_COLORS.len()]);
                    color_idx += 1;
                }
            }
        }
    }

    pu_colors
}

/// Assign unique colors to domains across all partitionings (matching original Python logic).
///
/// Returns a map from sorted domain boundaries to color.
pub fn assign_domain_colors(
    partitions: &[crate::sword::SwordPartition],
) -> HashMap<Vec<(i32, i32)>, (u8, u8, u8)> {
    let mut dom_colors = HashMap::new();
    let mut color_idx = 0usize;

    for part in partitions {
        for domain in &part.boundaries {
            let mut sorted = domain.clone();
            sorted.sort_by_key(|&(s, _)| s);
            if let Entry::Vacant(e) = dom_colors.entry(sorted) {
                e.insert(DOMAIN_COLORS[color_idx % DOMAIN_COLORS.len()]);
                color_idx += 1;
            }
        }
    }

    dom_colors
}

/// Generate all 3-level PNG plots for a single alternative partitioning.
///
/// Produces:
/// - `contact_probability_matrix_alternative_{alt_idx}.png` (all PUs)
/// - `contact_probability_matrix_alternative_{alt_idx}_domain_{dom_idx}.png` (per domain)
/// - `contact_probability_matrix_alternative_{alt_idx}_domain_{dom_idx}_pu_{start}_{end}.png` (per PU)
pub fn generate_alternative_plots(
    matrix: &Array2<f64>,
    alt_idx: usize,
    partition: &crate::sword::SwordPartition,
    pu_colors: &HashMap<(i32, i32), (u8, u8, u8)>,
    output_dir: &Path,
) -> Result<()> {
    let n = matrix.nrows();

    // Collect all PUs for this alternative
    let mut all_pus: Vec<ColoredPu> = Vec::new();
    for domain in &partition.boundaries {
        for &(start, end) in domain {
            let color = pu_colors
                .get(&(start, end))
                .copied()
                .unwrap_or((200, 200, 200));
            all_pus.push(ColoredPu { start, end, color });
        }
    }

    // 1) Alternative-level plot: all PU rectangles
    let alt_title = if alt_idx == 0 {
        "Contact Probability Map of the\noptimal partition (all Protein Units)".to_string()
    } else {
        format!(
            "Contact Probability Map of the alternative\npartition n\u{b0}{} (all Protein Units)",
            alt_idx
        )
    };
    let alt_path = output_dir.join(format!("alt{}.png", alt_idx));
    write_contact_matrix_png(matrix, n, &alt_title, &all_pus, &alt_path, true)?;

    // 2) Domain-level plots
    for (j, domain) in partition.boundaries.iter().enumerate() {
        let dom_pus: Vec<ColoredPu> = domain
            .iter()
            .map(|&(start, end)| {
                let color = pu_colors
                    .get(&(start, end))
                    .copied()
                    .unwrap_or((200, 200, 200));
                ColoredPu { start, end, color }
            })
            .collect();

        let dom_title = if alt_idx == 0 {
            format!(
                "Contact Probability Map of the domain {}\nof the optimal partition",
                j + 1
            )
        } else {
            format!(
                "Contact Probability Map of the domain {}\nof the alternative partition n\u{b0}{}",
                j + 1,
                alt_idx
            )
        };
        let dom_path = output_dir.join(format!("alt{}_dom{}.png", alt_idx, j));
        write_contact_matrix_png(matrix, n, &dom_title, &dom_pus, &dom_path, true)?;

        // 3) PU-level plots
        for &(start, end) in domain {
            let color = pu_colors
                .get(&(start, end))
                .copied()
                .unwrap_or((200, 200, 200));
            let pu_vec = vec![ColoredPu { start, end, color }];

            let pu_title = if alt_idx == 0 {
                format!(
                    "Contact Probability Map of PU {}-{} of the domain {}\nof the optimal partition",
                    start, end, j + 1
                )
            } else {
                format!(
                    "Contact Probability Map of PU {}-{} of the domain {}\nof the alternative partition n\u{b0}{}",
                    start, end, j + 1, alt_idx
                )
            };
            let pu_path =
                output_dir.join(format!("alt{}_dom{}_pu_{}_{}.png", alt_idx, j, start, end));
            write_contact_matrix_png(matrix, n, &pu_title, &pu_vec, &pu_path, false)?;
        }
    }

    Ok(())
}

/// Write a contact probability matrix PNG with PU boundary rectangles and legend.
fn write_contact_matrix_png(
    matrix: &Array2<f64>,
    n: usize,
    title: &str,
    pus: &[ColoredPu],
    output_path: &Path,
    large_format: bool,
) -> Result<()> {
    register_embedded_font();

    // Matching original matplotlib: figsize=(6,9) @ 150 dpi = 900×1350 for large,
    // figsize=(5,6.5) @ 150 dpi = 750×975 for PU-level
    let (width, height) = if large_format {
        (900u32, 1350u32)
    } else {
        (750u32, 975u32)
    };
    let title_font_size = if large_format { 24 } else { 22 };
    let axis_font_size = if large_format { 18 } else { 16 };
    let axis_desc_font_size = if large_format { 20 } else { 18 };
    let legend_title_font_size = if large_format { 18 } else { 17 };
    let legend_label_font_size = if large_format { 16 } else { 15 };

    let path_str = output_path.to_string_lossy().to_string();
    let root = BitMapBackend::new(&path_str, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    // Layout: title at top, plot area, legend at bottom
    let title_height = if large_format { 84u32 } else { 76u32 };
    let legend_height = if large_format { 140u32 } else { 108u32 };
    let margin = if large_format { 24u32 } else { 20u32 };
    let label_area = if large_format { 60u32 } else { 56u32 };

    // Title area
    let (title_area, rest) = root.split_vertically(title_height);
    let (plot_area, legend_area) = rest.split_vertically(height - title_height - legend_height);
    let plot_height = height - title_height - legend_height;
    let square_side = width.min(plot_height);
    let horizontal_padding = (width - square_side) / 2;
    let vertical_padding = (plot_height - square_side) / 2;

    let (_, plot_area) = plot_area.split_vertically(vertical_padding);
    let (plot_area, _) = plot_area.split_vertically(square_side);
    let (_, plot_area) = plot_area.split_horizontally(horizontal_padding);
    let (plot_area, _) = plot_area.split_horizontally(square_side);

    // Draw title
    for (line_idx, line) in title.lines().enumerate() {
        title_area.draw_text(
            line,
            &TextStyle::from(("sans-serif", title_font_size).into_font()).color(&BLACK),
            (
                width as i32 / 2 - (line.len() as i32 * (title_font_size / 4)),
                10 + line_idx as i32 * (title_font_size + 6),
            ),
        )?;
    }

    // Build chart in plot area
    let mut chart = ChartBuilder::on(&plot_area)
        .margin(margin)
        .x_label_area_size(label_area)
        .y_label_area_size(label_area)
        .build_cartesian_2d(0i32..n as i32, 0i32..n as i32)?;

    chart
        .configure_mesh()
        .x_desc("Residues")
        .y_desc("Residues")
        .label_style(("sans-serif", axis_font_size))
        .axis_desc_style(("sans-serif", axis_desc_font_size))
        .draw()?;

    // Find max value for color scaling
    let max_val = matrix.iter().cloned().fold(0.0f64, f64::max).max(1e-10);

    // Draw the heatmap (Y-axis inverted like matplotlib's imshow)
    for i in 0..n {
        for j in 0..n {
            let val = matrix[[i, j]];
            let t = (val / max_val).min(1.0);
            let (r, g, b) = rdpu_color(t);
            let color = RGBColor(r, g, b);

            // imshow convention: row 0 at top → invert Y
            let yi = (n - 1 - i) as i32;
            chart.draw_series(std::iter::once(Rectangle::new(
                [(j as i32, yi), (j as i32 + 1, yi + 1)],
                color.filled(),
            )))?;
        }
    }

    // Draw PU boundary rectangles (linewidth=~2px to approximate matplotlib 1.5pt @150dpi)
    for pu in pus {
        let (r, g, b) = pu.color;
        let color = RGBColor(r, g, b);
        let s = pu.start - 1; // 0-based
        let l = pu.end - pu.start; // length (matplotlib uses l = end - start, NOT +1)
        let e = s + l;

        // Inverted Y: top of rect = n - s - l, bottom = n - s
        let y_top = (n as i32) - s - l;
        let y_bot = (n as i32) - s;

        chart.draw_series(std::iter::once(Rectangle::new(
            [(s, y_top), (e, y_bot)],
            ShapeStyle::from(color).stroke_width(2),
        )))?;
    }

    // Draw legend in legend area
    let legend_x_start = margin as i32 + 30;
    let legend_y_start = 12i32;
    let cols = 3;
    let col_width = ((width as i32 - 2 * legend_x_start) / cols).max(150);
    let row_height = 28;

    // Legend title
    let legend_title = if pus.len() == 1 {
        "Protein Unit"
    } else {
        "Protein Units"
    };
    legend_area.draw_text(
        legend_title,
        &TextStyle::from(("sans-serif", legend_title_font_size).into_font()).color(&BLACK),
        (
            width as i32 / 2 - (legend_title.len() as i32 * 4),
            legend_y_start,
        ),
    )?;

    for (i, pu) in pus.iter().enumerate() {
        let col = (i % cols as usize) as i32;
        let row = (i / cols as usize) as i32;
        let x = legend_x_start + col * col_width;
        let y = legend_y_start + 24 + row * row_height;

        // Color swatch
        let (r, g, b) = pu.color;
        let color = RGBColor(r, g, b);
        legend_area.draw(&plotters::element::Rectangle::new(
            [(x, y), (x + 18, y + 16)],
            color.filled(),
        ))?;

        // Label
        let label = format!("{}-{}", pu.start, pu.end);
        legend_area.draw_text(
            &label,
            &TextStyle::from(("sans-serif", legend_label_font_size).into_font()).color(&BLACK),
            (x + 24, y - 1),
        )?;
    }

    root.present()
        .map_err(|e| anyhow::anyhow!("Failed to write PNG: {}", e))?;
    Ok(())
}

/// Domain count info for histogram.
#[derive(Debug, Clone)]
pub struct DomainCount {
    pub label: String,
    pub count: usize,
    pub color: (u8, u8, u8),
}

/// Generate a domain consistency histogram as SVG.
pub fn write_domain_histogram(domain_counts: &[DomainCount], output_path: &str) -> Result<()> {
    if domain_counts.is_empty() {
        return Ok(());
    }
    register_embedded_font();

    let width = 960u32;
    let height = 620u32;
    let caption_font_size = 22;
    let label_font_size = 16;
    let axis_desc_font_size = 18;

    let root = SVGBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_count = domain_counts.iter().map(|d| d.count).max().unwrap_or(1);
    let n_bars = domain_counts.len();

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Consistency of domains determined by SWORD",
            ("sans-serif", caption_font_size),
        )
        .margin(30)
        .x_label_area_size(130)
        .y_label_area_size(70)
        .build_cartesian_2d(0i32..n_bars as i32, 0i32..(max_count as i32 + 1))?;

    chart
        .configure_mesh()
        .x_desc("SWORD Domains")
        .y_desc("Count")
        .label_style(("sans-serif", label_font_size))
        .axis_desc_style(("sans-serif", axis_desc_font_size))
        .x_labels(n_bars)
        .x_label_formatter(&|idx| {
            domain_counts
                .get(*idx as usize)
                .map(|d| d.label.clone())
                .unwrap_or_default()
        })
        .draw()?;

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
pub fn count_domains(partitions: &[crate::sword::SwordPartition]) -> Vec<DomainCount> {
    let mut domain_freq: HashMap<Vec<(i32, i32)>, usize> = HashMap::new();

    for part in partitions {
        for domain in &part.boundaries {
            let mut sorted = domain.clone();
            sorted.sort_by_key(|&(s, _)| s);
            *domain_freq.entry(sorted).or_insert(0) += 1;
        }
    }

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
    fn png_plot_renders_labels_without_system_fonts() {
        let dir = tempfile::tempdir().expect("temporary output directory");
        let output = dir.path().join("plot.png");
        let matrix = Array2::from_elem((1, 1), 1.0);

        write_contact_matrix_png(&matrix, 1, "Plot title", &[], &output, false)
            .expect("write labeled PNG");

        assert!(output.metadata().expect("PNG metadata").len() > 0);
    }

    #[test]
    fn svg_histogram_renders_labels_without_system_fonts() {
        let dir = tempfile::tempdir().expect("temporary output directory");
        let output = dir.path().join("histogram.svg");
        let counts = [DomainCount {
            label: "(1, 20)".to_string(),
            count: 1,
            color: DOMAIN_COLORS[0],
        }];

        write_domain_histogram(&counts, output.to_str().expect("UTF-8 path"))
            .expect("write labeled SVG");

        assert!(output.metadata().expect("SVG metadata").len() > 0);
    }

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
