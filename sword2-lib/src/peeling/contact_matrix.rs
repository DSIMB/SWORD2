//! Contact probability matrix with 2D cumulative sum table.
//!
//! The contact probability between two C-alpha atoms is defined by the
//! logistic function:
//!
//! $$p(d) = \frac{1}{1 + e^{(d - D_0) / \Delta}}$$
//!
//! where $d$ is the Euclidean distance between the two C-alpha atoms,
//! $D_0$ is the midpoint distance (default 6.0 Å), and $\Delta$ controls
//! the steepness (default 1.5 Å).
//!
//! A precomputed 2D prefix-sum table enables O(1) rectangle sum queries,
//! which is critical for the iterative cutting algorithm.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use rayon::prelude::*;

/// Contact probability matrix with precomputed cumulative sums.
///
/// Stores an N×N symmetric matrix of contact probabilities and an
/// (N+1)×(N+1) cumulative sum table for O(1) rectangle queries.
/// Both are stored as flat `Vec<f64>` for cache-friendly access.
#[derive(Debug, Clone)]
pub struct ContactMatrix {
    /// Number of residues.
    n: usize,
    /// Contact probabilities, row-major N×N.
    data: Vec<f64>,
    /// Cumulative sum table, row-major (N+1)×(N+1).
    /// cum[i][j] = sum of data[0..i-1][0..j-1].
    cum: Vec<f64>,
    /// Total contact probability (sum of all entries).
    pub total_contact: f64,
}

impl ContactMatrix {
    /// Compute the contact probability matrix from C-alpha coordinates.
    ///
    /// Uses the logistic function p = 1 / (1 + exp((d - d0) / delta)).
    /// The matrix computation is parallelized with rayon.
    pub fn from_ca_coords(ca_coords: &[[f64; 3]], d0: f64, delta: f64) -> Self {
        let n = ca_coords.len();
        let mut data = vec![0.0f64; n * n];

        // Compute upper triangle (including diagonal) in parallel.
        // Each row is independent so we parallelize over rows.
        let row_sums: Vec<(Vec<f64>, f64)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut row = vec![0.0f64; n];
                let mut row_total = 0.0f64;
                let (xi, yi, zi) = (ca_coords[i][0], ca_coords[i][1], ca_coords[i][2]);
                for j in i..n {
                    let dx = xi - ca_coords[j][0];
                    let dy = yi - ca_coords[j][1];
                    let dz = zi - ca_coords[j][2];
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                    let p = 1.0 / (1.0 + ((dist - d0) / delta).exp());
                    row[j] = p;
                    if i != j {
                        row_total += 2.0 * p;
                    } else {
                        row_total += p;
                    }
                }
                (row, row_total)
            })
            .collect();

        let mut total_contact = 0.0f64;
        for (i, (row, row_total)) in row_sums.into_iter().enumerate() {
            total_contact += row_total;
            for j in i..n {
                let p = row[j];
                data[i * n + j] = p;
                data[j * n + i] = p;
            }
        }

        // Build cumulative sum table, size (n+1) × (n+1).
        // cum[i][j] = sum of data[0..i-1][0..j-1]
        let cum_n = n + 1;
        let mut cum = vec![0.0f64; cum_n * cum_n];
        for i in 1..cum_n {
            for j in 1..cum_n {
                cum[i * cum_n + j] = data[(i - 1) * n + (j - 1)]
                    + cum[(i - 1) * cum_n + j]
                    + cum[i * cum_n + (j - 1)]
                    - cum[(i - 1) * cum_n + (j - 1)];
            }
        }

        Self {
            n,
            data,
            cum,
            total_contact,
        }
    }

    /// Number of residues.
    #[inline]
    pub fn len(&self) -> usize {
        self.n
    }

    /// Whether the matrix is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Get contact probability between residues i and j (0-indexed).
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n + j]
    }

    /// Sum of contact probabilities over a rectangle [row_start..=row_end, col_start..=col_end].
    ///
    /// All indices are 0-based. Uses the cumulative sum table for O(1) computation.
    /// Matches the C `get_rectangle_sum()` function exactly.
    #[inline]
    pub fn rectangle_sum(
        &self,
        row_start: usize,
        col_start: usize,
        row_end: usize,
        col_end: usize,
    ) -> f64 {
        // Convert to 1-based indices into the cumulative sum table
        let r1 = row_start; // cum index: row_start (before +1, this is the "start-1")
        let c1 = col_start;
        let r2 = row_end + 1;
        let c2 = col_end + 1;
        let cn = self.n + 1;
        self.cum[r2 * cn + c2] - self.cum[r1 * cn + c2] - self.cum[r2 * cn + c1]
            + self.cum[r1 * cn + c1]
    }

    /// Write the contact probability matrix to a file in the format
    /// expected by the plot module and compatible with the C output.
    pub fn write_matrix_file(&self, path: &Path) -> Result<()> {
        let mut buf = Vec::with_capacity(self.n * self.n * 9);
        writeln!(buf, "# Contact probability matrix")?;
        for i in 0..self.n {
            for j in 0..self.n {
                if j > 0 {
                    write!(buf, " ")?;
                }
                write!(buf, "{:7.5}", self.data[i * self.n + j])?;
            }
            writeln!(buf)?;
        }
        std::fs::write(path, buf)
            .with_context(|| format!("Cannot write contact matrix to {}", path.display()))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contact_matrix_small() {
        // Three residues at known positions
        let coords = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let mat = ContactMatrix::from_ca_coords(&coords, 6.0, 1.5);
        assert_eq!(mat.len(), 3);

        // Diagonal entries: distance=0, p = 1/(1+exp(-4)) ≈ 0.982
        let diag = mat.get(0, 0);
        assert!((diag - 1.0 / (1.0 + (-6.0f64 / 1.5).exp())).abs() < 1e-10);

        // Symmetry
        assert!((mat.get(0, 1) - mat.get(1, 0)).abs() < 1e-15);
        assert!((mat.get(0, 2) - mat.get(2, 0)).abs() < 1e-15);

        // Rectangle sum of full matrix should equal total_contact
        let full_sum = mat.rectangle_sum(0, 0, 2, 2);
        assert!((full_sum - mat.total_contact).abs() < 1e-10);

        // Rectangle sum of single element
        let single = mat.rectangle_sum(1, 1, 1, 1);
        assert!((single - mat.get(1, 1)).abs() < 1e-15);
    }

    #[test]
    fn test_rectangle_sum_consistency() {
        // 4 residues, verify rectangle sums against brute force
        let coords = vec![
            [0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0],
            [0.0, 5.0, 0.0],
            [5.0, 5.0, 0.0],
        ];
        let mat = ContactMatrix::from_ca_coords(&coords, 6.0, 1.5);

        // Brute-force sum over [1..=2, 0..=1]
        let brute = mat.get(1, 0) + mat.get(1, 1) + mat.get(2, 0) + mat.get(2, 1);
        let fast = mat.rectangle_sum(1, 0, 2, 1);
        assert!((brute - fast).abs() < 1e-12);
    }
}
