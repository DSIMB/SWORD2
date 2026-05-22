//! H-bond detection with spatial grid optimization.
//!
//! The original DSSP uses O(N^2) pairwise search with a 9.0A CA distance cutoff.
//! We use a spatial grid (cell size = CADIST) to only check 27 neighbor cells,
//! reducing to O(N*k) where k is the average number of neighbors.

use super::types::{DsspChain, HydrogenBond, CADIST, DIST_MIN, HBHIGH, HBLOW, Q};

/// Detect all backbone H-bonds in the chain using a spatial grid.
pub fn flag_hydrogen_bonds(chain: &mut DsspChain) {
    let len = chain.len;
    if len < 2 {
        return;
    }

    // Build spatial grid for CA atoms
    let grid = SpatialGrid::build(chain);

    // For each residue, find neighbors via grid and compute H-bonds
    for i in 1..=len {
        if !chain.no_chain_break(i, i) {
            continue;
        }

        let neighbors = grid.neighbors(chain.get(i).ca.x, chain.get(i).ca.y, chain.get(i).ca.z);

        for &j in &neighbors {
            if j <= i {
                continue;
            }
            if !chain.no_chain_break(j, j) {
                continue;
            }

            let dist_sq = {
                let ci = &chain.get(i).ca;
                let cj = &chain.get(j).ca;
                let dx = ci.x - cj.x;
                let dy = ci.y - cj.y;
                let dz = ci.z - cj.z;
                dx * dx + dy * dy + dz * dz
            };

            if dist_sq >= CADIST * CADIST {
                continue;
            }

            // i is donor (NH), j is acceptor (CO)
            set_bonds(chain, i, j);
            // j is donor (NH), i is acceptor (CO) — skip i,i+1 pair
            if j != i + 1 {
                set_bonds(chain, j, i);
            }
        }
    }
}

/// Compute H-bond energy between donor i (NH) and acceptor j (CO).
/// Returns energy in cal/mol.
fn bond_energy(chain: &DsspChain, donor: usize, acceptor: usize) -> i64 {
    let d = chain.get(donor);
    let a = chain.get(acceptor);

    // Proline has no amide H
    if d.aa == 'P' {
        return 0;
    }

    let dho = d.h.distance_to(&a.o);
    let dhc = d.h.distance_to(&a.c);
    let dnc = d.n.distance_to(&a.c);
    let dno = d.n.distance_to(&a.o);

    if dho < DIST_MIN || dhc < DIST_MIN || dnc < DIST_MIN || dno < DIST_MIN {
        return HBLOW;
    }

    let e = (Q / dho - Q / dhc + Q / dnc - Q / dno + 0.5) as i64;
    e.max(HBLOW)
}

/// Update H-bond arrays for donor i and acceptor j.
fn set_bonds(chain: &mut DsspChain, donor: usize, acceptor: usize) {
    let energy = bond_energy(chain, donor, acceptor);

    // Update donor's acceptor list (CO(j) is acceptor of NH(i))
    let hb = HydrogenBond {
        residue: acceptor,
        energy,
    };
    update_bonds(&mut chain.get_mut(donor).acceptor, hb);

    // Update acceptor's donor list (NH(i) donates to CO(j))
    let hb = HydrogenBond {
        residue: donor,
        energy,
    };
    update_bonds(&mut chain.get_mut(acceptor).donor, hb);
}

/// Keep the two lowest-energy H-bonds.
fn update_bonds(bonds: &mut [HydrogenBond; 2], hb: HydrogenBond) {
    if hb.energy < bonds[0].energy {
        bonds[1] = bonds[0];
        bonds[0] = hb;
    } else if hb.energy < bonds[1].energy {
        bonds[1] = hb;
    }
}

/// Test if residue i's NH donates to residue j's CO.
pub fn test_bond(chain: &DsspChain, i: usize, j: usize) -> bool {
    let acc = &chain.get(i).acceptor;
    (acc[0].residue == j && acc[0].energy < HBHIGH)
        || (acc[1].residue == j && acc[1].energy < HBHIGH)
}

// ---------------------------------------------------------------------------
// Spatial grid for O(N*k) neighbor lookups
// ---------------------------------------------------------------------------

struct SpatialGrid {
    cell_size: f64,
    // Grid dimensions
    nx: usize,
    ny: usize,
    nz: usize,
    // Origin
    ox: f64,
    oy: f64,
    oz: f64,
    // cells[flat_index] = list of 1-based residue indices
    cells: Vec<Vec<usize>>,
}

impl SpatialGrid {
    fn build(chain: &DsspChain) -> Self {
        let cell_size = CADIST;

        if chain.len == 0 {
            return Self {
                cell_size,
                nx: 0,
                ny: 0,
                nz: 0,
                ox: 0.0,
                oy: 0.0,
                oz: 0.0,
                cells: Vec::new(),
            };
        }

        // Find bounding box of CA atoms (only non-break residues)
        let mut min_x = f64::MAX;
        let mut min_y = f64::MAX;
        let mut min_z = f64::MAX;
        let mut max_x = f64::MIN;
        let mut max_y = f64::MIN;
        let mut max_z = f64::MIN;

        for i in 1..=chain.len {
            if chain.get(i).aa == '!' {
                continue;
            }
            let ca = &chain.get(i).ca;
            min_x = min_x.min(ca.x);
            min_y = min_y.min(ca.y);
            min_z = min_z.min(ca.z);
            max_x = max_x.max(ca.x);
            max_y = max_y.max(ca.y);
            max_z = max_z.max(ca.z);
        }

        let ox = min_x - cell_size;
        let oy = min_y - cell_size;
        let oz = min_z - cell_size;
        let nx = ((max_x - ox) / cell_size) as usize + 2;
        let ny = ((max_y - oy) / cell_size) as usize + 2;
        let nz = ((max_z - oz) / cell_size) as usize + 2;

        let total = nx * ny * nz;
        let mut cells = vec![Vec::new(); total];

        for i in 1..=chain.len {
            if chain.get(i).aa == '!' {
                continue;
            }
            let ca = &chain.get(i).ca;
            let cx = ((ca.x - ox) / cell_size) as usize;
            let cy = ((ca.y - oy) / cell_size) as usize;
            let cz = ((ca.z - oz) / cell_size) as usize;
            let idx = cx * ny * nz + cy * nz + cz;
            if idx < total {
                cells[idx].push(i);
            }
        }

        Self {
            cell_size,
            nx,
            ny,
            nz,
            ox,
            oy,
            oz,
            cells,
        }
    }

    /// Return all residue indices within CADIST of the given point.
    fn neighbors(&self, x: f64, y: f64, z: f64) -> Vec<usize> {
        let cx = ((x - self.ox) / self.cell_size) as isize;
        let cy = ((y - self.oy) / self.cell_size) as isize;
        let cz = ((z - self.oz) / self.cell_size) as isize;

        let mut result = Vec::new();

        for dx in -1..=1isize {
            let gx = cx + dx;
            if gx < 0 || gx >= self.nx as isize {
                continue;
            }
            for dy in -1..=1isize {
                let gy = cy + dy;
                if gy < 0 || gy >= self.ny as isize {
                    continue;
                }
                for dz in -1..=1isize {
                    let gz = cz + dz;
                    if gz < 0 || gz >= self.nz as isize {
                        continue;
                    }
                    let idx = gx as usize * self.ny * self.nz + gy as usize * self.nz + gz as usize;
                    if idx < self.cells.len() {
                        result.extend_from_slice(&self.cells[idx]);
                    }
                }
            }
        }

        result
    }
}
