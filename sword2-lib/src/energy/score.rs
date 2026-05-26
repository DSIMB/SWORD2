//! Pure-Rust pseudo-energy scorer using the mypmfs statistical potentials
//! (CA representation, linear interpolation).
//!
//! Algorithm: load `.nrg` potentials + `parameters.log` + `xvector.dat` → extract
//! CA atoms → pairwise distances within `[distmin, distmax]` and sequence
//! separation `(diffmin, diffmax)` → per atom-pair linear-interpolate the
//! potential by distance → sum = pseudo-energy. Z-score: shuffle the sequence
//! `num_shuffles` times (similarity-constrained), recompute, `z = (E - μ) / σ`.
//!
//! Scope: linear interpolation only (no cubic spline); CA/BB representations only.
//!
//! Implementation note: at distances in `[last_bin, distmax)` we clamp to the
//! last valid bin interval, avoiding an invalid access at the upper edge.
//! See `KNOWN_DIVERGENCES.md`.

/// Loaded pseudo-energy potentials and scoring parameters for one representation.
#[derive(Debug)]
pub struct Potentials {
    /// Distance bin centers (Å), ascending, evenly spaced.
    xvector: Vec<f64>,
    /// Energy table: key = sorted-concat of two atom-type labels (e.g. "ACARCA"),
    /// value = one energy per distance bin (aligned with `xvector`).
    table: std::collections::HashMap<String, Vec<f64>>,
    /// Minimum interatomic distance (Å) for a pair to count.
    distmin: f64,
    /// Maximum interatomic distance (Å) for a pair to count.
    distmax: f64,
    /// Minimum sequence separation (exclusive) for an intra-chain pair.
    diffmin: i32,
    /// Maximum sequence separation (exclusive) for an intra-chain pair.
    diffmax: i32,
    /// Atom-type labels for the active representation (e.g. the 20 CA types).
    atypes: std::collections::BTreeSet<String>,
}

/// One-letter codes for the 20 standard amino acids.
const AA_ONE: [&str; 20] = [
    "A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
    "V",
];

impl Potentials {
    /// Load potentials and parameters from a mypmfs potential directory.
    ///
    /// Reads `parameters.log` (DIFFMIN/DIFFMAX/DISTMAX/DISTMIN/REPRES),
    /// `xvector.dat` (distance bins), and one `<pair>.nrg` file per unordered
    /// atom-type pair. Only the CA and BB representations are supported (the only
    /// ones the SWORD2 pipeline uses).
    pub fn load(potential_dir: &str) -> anyhow::Result<Self> {
        use anyhow::Context;
        let dir = std::path::Path::new(potential_dir);

        // parameters.log — parse by key (KEY=VALUE per line).
        let params = std::fs::read_to_string(dir.join("parameters.log"))
            .with_context(|| format!("reading {potential_dir}/parameters.log"))?;
        let mut diffmin = None;
        let mut diffmax = None;
        let mut distmax = None;
        let mut distmin = None;
        let mut repres = None;
        for line in params.lines() {
            if let Some((k, v)) = line.split_once('=') {
                let v = v.trim();
                match k.trim() {
                    "DIFFMIN" => diffmin = v.parse().ok(),
                    "DIFFMAX" => diffmax = v.parse().ok(),
                    "DISTMAX" => distmax = v.parse().ok(),
                    "DISTMIN" => distmin = v.parse().ok(),
                    "REPRES" => repres = Some(v.to_string()),
                    _ => {}
                }
            }
        }
        let diffmin = diffmin.context("DIFFMIN missing in parameters.log")?;
        let diffmax = diffmax.context("DIFFMAX missing in parameters.log")?;
        let distmax = distmax.context("DISTMAX missing in parameters.log")?;
        let distmin = distmin.context("DISTMIN missing in parameters.log")?;
        let repres = repres.context("REPRES missing in parameters.log")?;

        let suffix = match repres.as_str() {
            "CA" => "CA",
            "BB" => "BB",
            other => {
                anyhow::bail!("unsupported representation {other:?}; only CA and BB are supported")
            }
        };
        let atypes: std::collections::BTreeSet<String> =
            AA_ONE.iter().map(|aa| format!("{aa}{suffix}")).collect();

        // xvector.dat — distance bins.
        let xvec_raw = std::fs::read_to_string(dir.join("xvector.dat"))
            .with_context(|| format!("reading {potential_dir}/xvector.dat"))?;
        let xvector: Vec<f64> = xvec_raw
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.trim().parse::<f64>())
            .collect::<Result<_, _>>()
            .context("parsing xvector.dat")?;

        // One .nrg file per unordered pair (sorted-concat key).
        let sorted: Vec<&String> = atypes.iter().collect();
        let mut table = std::collections::HashMap::new();
        for i in 0..sorted.len() {
            for j in i..sorted.len() {
                let key = format!("{}{}", sorted[i], sorted[j]);
                let path = dir.join(format!("{key}.nrg"));
                let raw = std::fs::read_to_string(&path)
                    .with_context(|| format!("reading {}", path.display()))?;
                let energies: Vec<f64> = raw
                    .lines()
                    .filter(|l| !l.trim().is_empty())
                    .map(|l| l.trim().parse::<f64>())
                    .collect::<Result<_, _>>()
                    .with_context(|| format!("parsing {}", path.display()))?;
                table.insert(key, energies);
            }
        }

        Ok(Potentials {
            xvector,
            table,
            distmin,
            distmax,
            diffmin,
            diffmax,
            atypes,
        })
    }

    /// Bin index and interpolation fraction for distance `x`.
    ///
    /// Returns `(i, frac)` such that the interpolated value is
    /// `table[i] * (1 - frac) + table[i+1] * frac`. The index is clamped to a
    /// valid interval `[i, i+1]` at the upper boundary.
    fn bin_and_frac(&self, x: f64) -> (usize, f64) {
        let last = self.xvector.len() - 1;
        let mut i = last;
        while i > 0 && self.xvector[i] > x {
            i -= 1;
        }
        if i > last - 1 {
            i = last - 1;
        }
        let (x0, x1) = (self.xvector[i], self.xvector[i + 1]);
        (i, (x - x0) / (x1 - x0))
    }

    /// Linear-interpolate the potential for `atom_pair` at distance `x`.
    fn interpolate(&self, x: f64, atom_pair: &str) -> f64 {
        let energies = &self.table[atom_pair];
        let (i, frac) = self.bin_and_frac(x);
        energies[i] * (1.0 - frac) + energies[i + 1] * frac
    }
}

/// A scoring atom (one per residue in the CA representation).
#[derive(Debug, Clone)]
struct ScoreAtom {
    resnum: i32,
    /// Atom-type label: one-letter residue code + atom name (e.g. "ACA").
    label: String,
    x: f64,
    y: f64,
    z: f64,
}

/// A retained interatomic interaction: the distance and the two atom indices.
#[derive(Debug, Clone)]
struct Interaction {
    dist: f64,
    a: usize,
    b: usize,
}

/// Build the retained intra-chain interactions (ports `distances()` for a single
/// chain): all i<j atom pairs whose sequence separation is in `(diffmin, diffmax)`
/// and whose Euclidean distance is in `(distmin, distmax)`.
fn build_interactions(
    atoms: &[ScoreAtom],
    distmin: f64,
    distmax: f64,
    diffmin: i32,
    diffmax: i32,
) -> Vec<Interaction> {
    let mut interactions = Vec::new();
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let numdiff = atoms[j].resnum - atoms[i].resnum;
            if numdiff > diffmin && numdiff < diffmax {
                let (dx, dy, dz) = (
                    atoms[i].x - atoms[j].x,
                    atoms[i].y - atoms[j].y,
                    atoms[i].z - atoms[j].z,
                );
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < distmax && dist > distmin {
                    interactions.push(Interaction { dist, a: i, b: j });
                }
            }
        }
    }
    interactions
}

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

/// Max sequence identity with the query allowed for a shuffled decoy.
const SIM_THRES: f64 = 0.5;

/// Permute a per-residue sequence until it is at most `simthres` identical to the
/// original (ports the decoy-shuffle inner loop).
///
/// Returns the new value assigned to each residue index. A safety cap prevents a
/// hang on degenerate inputs (n < 2, or too few distinct values to drop below
/// `simthres`).
fn shuffle_seq<T: Clone + PartialEq + Default>(
    seq: &[T],
    simthres: f64,
    rng: &mut StdRng,
) -> Vec<T> {
    let n = seq.len();
    if n < 2 {
        return seq.to_vec();
    }
    // perm[i] = original residue index currently sitting at sequence position i.
    let mut perm: Vec<usize> = (0..n).collect();
    let mut similarity = 1.0;
    let cap = 100 * n;
    let mut iters = 0;
    while similarity > simthres && iters < cap {
        iters += 1;
        let i = rng.gen_range(0..n);
        let j = rng.gen_range(0..n);
        if i != j {
            perm.swap(i, j);
            let count = (0..n).filter(|&k| seq[perm[k]] == seq[k]).count();
            similarity = count as f64 / n as f64;
        }
    }
    // The residue that landed at position i takes the value originally at i.
    let mut new_seq = vec![T::default(); n];
    for i in 0..n {
        new_seq[perm[i]] = seq[i].clone();
    }
    new_seq
}

/// Compute the Z-score: `z = (E - μ) / σ` over `num_shuffles` decoy energies.
/// Returns `Some(0.0)` when the query energy is zero and `None` when sigma is
/// zero (undefined Z-score). Decoys are seeded for reproducibility.
///
/// Hot path is allocation-free: residue labels are mapped to small integer type
/// indices, potential curves are gathered into a dense `K×K` table, and each
/// interaction's `(bin, frac)` is precomputed — so each decoy energy is a sum of
/// array lookups and two multiply-adds, with no string hashing.
fn z_score(
    pot: &Potentials,
    atoms: &[ScoreAtom],
    interactions: &[Interaction],
    num_shuffles: usize,
    base_seed: u64,
) -> Option<f64> {
    let total = raw_energy(pot, atoms, interactions);
    if total == 0.0 || interactions.is_empty() {
        return Some(0.0);
    }

    // Map atom-type labels → indices 0..k over the representation's atom types.
    let sorted_types: Vec<&String> = pot.atypes.iter().collect();
    let k = sorted_types.len();
    let type_index: std::collections::HashMap<&str, usize> = sorted_types
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    // curve_of[a * k + b] = potential curve for the (a, b) type pair.
    let curve_of: Vec<&Vec<f64>> = (0..k)
        .flat_map(|a| (0..k).map(move |b| (a, b)))
        .map(|(a, b)| &pot.table[&pair_key(sorted_types[a], sorted_types[b])])
        .collect();
    let types: Vec<usize> = atoms.iter().map(|a| type_index[a.label.as_str()]).collect();

    // Per-interaction precomputation: endpoints + interpolation weights.
    struct Pre {
        a: usize,
        b: usize,
        bin: usize,
        w0: f64,
        w1: f64,
    }
    let pre: Vec<Pre> = interactions
        .iter()
        .map(|it| {
            let (bin, frac) = pot.bin_and_frac(it.dist);
            Pre {
                a: it.a,
                b: it.b,
                bin,
                w0: 1.0 - frac,
                w1: frac,
            }
        })
        .collect();

    let rand_energies: Vec<f64> = (0..num_shuffles)
        .into_par_iter()
        .map(|d| {
            // Distinct, reproducible seed per decoy.
            let seed = base_seed ^ (d as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
            let mut rng = StdRng::seed_from_u64(seed);
            let new_types = shuffle_seq(&types, SIM_THRES, &mut rng);
            pre.iter()
                .map(|p| {
                    let c = curve_of[new_types[p.a] * k + new_types[p.b]];
                    c[p.bin] * p.w0 + c[p.bin + 1] * p.w1
                })
                .sum::<f64>()
        })
        .collect();

    let n = rand_energies.len() as f64;
    let mean = rand_energies.iter().sum::<f64>() / n;
    let variance = rand_energies
        .iter()
        .map(|e| (e - mean).powi(2))
        .sum::<f64>()
        / n;
    let stdev = variance.sqrt();
    if stdev == 0.0 {
        return None;
    }
    Some((total - mean) / stdev)
}

/// Fixed RNG seed so Z-scores are reproducible across runs.
const DEFAULT_SEED: u64 = 0x00C0_FFEE_0000_0001;

/// Score a structure: pseudo-energy and (optionally) Z-score, for the whole
/// structure or a residue subset (`residue_list`, comma-separated `numchain`
/// tokens like `"12A,13A"`).
pub fn score(
    pot: &Potentials,
    pdb_path: &str,
    residue_list: Option<&str>,
    num_shuffles: usize,
    compute_z: bool,
) -> anyhow::Result<super::EnergyResult> {
    use anyhow::Context;
    let content = std::fs::read_to_string(pdb_path)
        .with_context(|| format!("reading structure {pdb_path}"))?;

    let filter: Option<std::collections::HashSet<String>> = residue_list.map(|list| {
        list.split(',')
            .map(|t| t.trim().to_string())
            .filter(|t| !t.is_empty())
            .collect()
    });
    let atoms = parse_atoms(&content, &pot.atypes, filter.as_ref());

    // No atoms in the selected residue subset means there is no score.
    if atoms.is_empty() {
        return Ok(super::EnergyResult {
            energy: None,
            z_score: None,
        });
    }

    let interactions =
        build_interactions(&atoms, pot.distmin, pot.distmax, pot.diffmin, pot.diffmax);
    let energy = raw_energy(pot, &atoms, &interactions);
    let z = if compute_z {
        z_score(pot, &atoms, &interactions, num_shuffles, DEFAULT_SEED)
    } else {
        None
    };
    Ok(super::EnergyResult {
        energy: Some(energy),
        z_score: z,
    })
}

/// Map a three-letter residue name to its one-letter code (20 standard AAs only).
fn three_to_one(resname: &str) -> Option<&'static str> {
    Some(match resname {
        "ALA" => "A",
        "ARG" => "R",
        "ASN" => "N",
        "ASP" => "D",
        "CYS" => "C",
        "GLU" => "E",
        "GLN" => "Q",
        "GLY" => "G",
        "HIS" => "H",
        "ILE" => "I",
        "LEU" => "L",
        "LYS" => "K",
        "MET" => "M",
        "PHE" => "F",
        "PRO" => "P",
        "SER" => "S",
        "THR" => "T",
        "TRP" => "W",
        "TYR" => "Y",
        "VAL" => "V",
        _ => return None,
    })
}

/// Parse the CA (or other-representation) atoms from PDB text: first model only;
/// `ATOM` records with a known
/// standard residue whose `<one-letter><atom-name>` label is in `atypes`; altloc
/// blank or "A"; optionally restricted to residues in `residue_filter` (tokens
/// like "12A" = resnum+chain).
fn parse_atoms(
    content: &str,
    atypes: &std::collections::BTreeSet<String>,
    residue_filter: Option<&std::collections::HashSet<String>>,
) -> Vec<ScoreAtom> {
    let mut atoms = Vec::new();
    for line in content.lines() {
        // Score only the first structure model.
        if line.starts_with("ENDMDL") {
            break;
        }
        // Coordinates require all columns through the z value at [46:54].
        if !line.starts_with("ATOM") || line.len() <= 53 {
            continue;
        }
        let bytes = line.as_bytes();
        let altloc = bytes[16] as char;
        if altloc != ' ' && altloc != 'A' {
            continue;
        }
        let resname = line[17..20].trim();
        let Some(one) = three_to_one(resname) else {
            continue;
        };
        let atomname = line[12..16].trim();
        let label = format!("{one}{atomname}");
        if !atypes.contains(&label) {
            continue;
        }
        let chain = bytes[21] as char;
        let resnum: i32 = match line[22..26].trim().parse() {
            Ok(n) => n,
            Err(_) => continue,
        };
        if let Some(filter) = residue_filter {
            if !filter.contains(&format!("{resnum}{chain}")) {
                continue;
            }
        }
        let (Ok(x), Ok(y), Ok(z)) = (
            line[30..38].trim().parse::<f64>(),
            line[38..46].trim().parse::<f64>(),
            line[46..54].trim().parse::<f64>(),
        ) else {
            continue;
        };
        atoms.push(ScoreAtom {
            resnum,
            label,
            x,
            y,
            z,
        });
    }
    atoms
}

/// Sorted-concat key for an atom pair, matching the `.nrg` filename convention
/// (lexicographically smaller label first).
fn pair_key(la: &str, lb: &str) -> String {
    if la <= lb {
        format!("{la}{lb}")
    } else {
        format!("{lb}{la}")
    }
}

/// Total pseudo-energy: sum over interactions of the interpolated potential at
/// each interaction's distance, keyed by the sorted pair of the two atom labels.
fn raw_energy(pot: &Potentials, atoms: &[ScoreAtom], interactions: &[Interaction]) -> f64 {
    interactions
        .iter()
        .map(|it| {
            let key = pair_key(&atoms[it.a].label, &atoms[it.b].label);
            pot.interpolate(it.dist, &key)
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{BTreeSet, HashMap, HashSet};

    fn atom(resnum: i32, label: &str, x: f64, y: f64, z: f64) -> ScoreAtom {
        ScoreAtom {
            resnum,
            label: label.to_string(),
            x,
            y,
            z,
        }
    }

    /// Format a PDB ATOM line with correct column placement.
    #[allow(clippy::too_many_arguments)]
    fn pdb_line(
        serial: i32,
        atom: &str,
        altloc: &str,
        resname: &str,
        chain: &str,
        resnum: i32,
        x: f64,
        y: f64,
        z: f64,
    ) -> String {
        // Atom name occupies cols 13-16; names <4 chars are left-padded by one space.
        let atom_field = if atom.len() >= 4 {
            atom.to_string()
        } else {
            format!(" {atom:<3}")
        };
        format!(
            "ATOM  {serial:>5} {atom_field}{altloc:<1}{resname:>3} {chain:<1}{resnum:>4}    {x:>8.3}{y:>8.3}{z:>8.3}  1.00  0.00"
        )
    }

    fn ca_types(labels: &[&str]) -> BTreeSet<String> {
        labels.iter().map(|s| s.to_string()).collect()
    }

    fn synthetic_potentials() -> Potentials {
        let mut table = HashMap::new();
        table.insert("XY".to_string(), vec![10.0, 20.0, 40.0, 80.0]);
        Potentials {
            xvector: vec![0.0, 1.0, 2.0, 3.0],
            table,
            distmin: 0.0,
            distmax: 15.0,
            diffmin: 3,
            diffmax: 5100,
            atypes: BTreeSet::new(),
        }
    }

    #[test]
    fn interpolate_midpoint_of_first_interval() {
        let p = synthetic_potentials();
        // x=0.5 → interval [0,1], 10 + 0.5*(20-10) = 15
        assert!((p.interpolate(0.5, "XY") - 15.0).abs() < 1e-12);
    }

    #[test]
    fn interpolate_midpoint_of_last_interval() {
        let p = synthetic_potentials();
        // x=2.5 → interval [2,3], 40 + 0.5*(80-40) = 60
        assert!((p.interpolate(2.5, "XY") - 60.0).abs() < 1e-12);
    }

    #[test]
    fn interpolate_at_and_beyond_last_bin_clamps() {
        let p = synthetic_potentials();
        // x=3.0 (== last bin): clamp to interval [2,3] → exactly table[3] = 80
        assert!((p.interpolate(3.0, "XY") - 80.0).abs() < 1e-12);
    }

    #[test]
    fn build_interactions_filters_by_separation_and_distance() {
        // resnum 1 @origin; 2 too close in sequence; 5 @4Å kept; 10 @30Å too far.
        let atoms = vec![
            atom(1, "ACA", 0.0, 0.0, 0.0),
            atom(2, "RCA", 3.0, 0.0, 0.0), // sep 1 (not > 3) → excluded
            atom(5, "ACA", 0.0, 4.0, 0.0), // sep 4, dist 4 → kept
            atom(10, "RCA", 0.0, 0.0, 30.0), // sep 9 but dist 30 (> 15) → excluded
        ];
        let inter = build_interactions(&atoms, 0.0, 15.0, 3, 5100);
        assert_eq!(inter.len(), 1);
        assert_eq!((inter[0].a, inter[0].b), (0, 2));
        assert!((inter[0].dist - 4.0).abs() < 1e-12);
    }

    #[test]
    fn raw_energy_sums_interpolated_pairs() {
        let p = synthetic_potentials(); // table {"XY": [10,20,40,80]}, bins [0,1,2,3]
        let atoms = vec![
            atom(1, "X", 0.0, 0.0, 0.0),
            atom(2, "Y", 0.0, 0.0, 0.0),
            atom(3, "X", 0.0, 0.0, 0.0),
        ];
        let inter = vec![
            Interaction {
                dist: 0.5,
                a: 0,
                b: 1,
            }, // XY @0.5 → 15
            Interaction {
                dist: 2.5,
                a: 1,
                b: 2,
            }, // XY @2.5 → 60
        ];
        assert!((raw_energy(&p, &atoms, &inter) - 75.0).abs() < 1e-12);
    }

    #[test]
    fn parse_atoms_keeps_ca_skips_altloc_and_nonstandard() {
        let atypes = ca_types(&["ACA", "RCA"]);
        let content = [
            pdb_line(1, "N", " ", "ALA", "A", 1, 0.0, 0.0, 0.0), // not CA → skip
            pdb_line(2, "CA", " ", "ALA", "A", 1, 1.0, 0.0, 0.0), // keep → ACA
            pdb_line(3, "CA", "B", "ARG", "A", 2, 2.0, 0.0, 0.0), // altloc B → skip
            pdb_line(4, "CA", " ", "ARG", "A", 2, 3.0, 0.0, 0.0), // keep → RCA
            pdb_line(5, "CA", " ", "UNK", "A", 3, 4.0, 0.0, 0.0), // non-standard → skip
        ]
        .join("\n");
        let atoms = parse_atoms(&content, &atypes, None);
        assert_eq!(atoms.len(), 2);
        assert_eq!(atoms[0].label, "ACA");
        assert_eq!((atoms[0].resnum, atoms[0].x), (1, 1.0));
        assert_eq!(atoms[1].label, "RCA");
        assert_eq!((atoms[1].resnum, atoms[1].x), (2, 3.0));
    }

    #[test]
    fn parse_atoms_applies_residue_filter() {
        let atypes = ca_types(&["ACA", "RCA"]);
        let content = [
            pdb_line(2, "CA", " ", "ALA", "A", 1, 1.0, 0.0, 0.0),
            pdb_line(4, "CA", " ", "ARG", "A", 2, 3.0, 0.0, 0.0),
        ]
        .join("\n");
        let filter: HashSet<String> = ["1A".to_string()].into_iter().collect();
        let atoms = parse_atoms(&content, &atypes, Some(&filter));
        assert_eq!(atoms.len(), 1);
        assert_eq!(atoms[0].resnum, 1);
    }

    fn potentials_xy() -> Potentials {
        // All four ordered label pairs over {X, Y}; sorted keys are XX, XY, YY.
        let mut table = HashMap::new();
        table.insert("XX".to_string(), vec![1.0, 2.0, 3.0, 4.0]);
        table.insert("XY".to_string(), vec![10.0, 20.0, 40.0, 80.0]);
        table.insert("YY".to_string(), vec![100.0, 200.0, 300.0, 400.0]);
        Potentials {
            xvector: vec![0.0, 1.0, 2.0, 3.0],
            table,
            distmin: 0.0,
            distmax: 15.0,
            diffmin: 3,
            diffmax: 5100,
            atypes: ["X".to_string(), "Y".to_string()].into_iter().collect(),
        }
    }

    #[test]
    fn shuffle_labels_is_permutation_below_similarity() {
        let labels: Vec<String> = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut rng = StdRng::seed_from_u64(123);
        let shuffled = shuffle_seq(&labels, 0.5, &mut rng);
        // Same multiset of labels.
        let mut a = labels.clone();
        a.sort();
        let mut b = shuffled.clone();
        b.sort();
        assert_eq!(a, b);
        // At most 50% of positions retain their original label.
        let same = labels.iter().zip(&shuffled).filter(|(o, n)| o == n).count();
        assert!(same as f64 / labels.len() as f64 <= 0.5);
    }

    #[test]
    fn shuffle_labels_is_deterministic_for_a_seed() {
        let labels: Vec<String> = ["X", "Y", "X", "Y", "X", "Y"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let s1 = shuffle_seq(&labels, 0.5, &mut StdRng::seed_from_u64(7));
        let s2 = shuffle_seq(&labels, 0.5, &mut StdRng::seed_from_u64(7));
        assert_eq!(s1, s2);
    }

    #[test]
    fn z_score_is_zero_without_interactions() {
        let p = potentials_xy();
        let atoms = vec![atom(1, "X", 0.0, 0.0, 0.0)];
        assert_eq!(z_score(&p, &atoms, &[], 100, 42), Some(0.0));
    }

    #[test]
    fn z_score_is_deterministic_for_a_seed() {
        let p = potentials_xy();
        let atoms: Vec<ScoreAtom> = (1..=8)
            .map(|i| atom(i, if i % 2 == 0 { "X" } else { "Y" }, i as f64, 0.0, 0.0))
            .collect();
        let inter = build_interactions(&atoms, 0.0, 15.0, 3, 5100);
        assert!(!inter.is_empty());
        let z1 = z_score(&p, &atoms, &inter, 200, 42);
        let z2 = z_score(&p, &atoms, &inter, 200, 42);
        assert_eq!(z1, z2);
    }

    /// Path to the committed potential directory used by the pipeline.
    fn real_potential_dir() -> String {
        format!(
            "{}/../bin/mypmfs-master/025_30_100_potential",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    #[test]
    fn load_real_potentials() {
        let p = Potentials::load(&real_potential_dir()).expect("load potentials");
        assert_eq!(p.diffmin, 3);
        assert_eq!(p.diffmax, 5100);
        assert_eq!(p.distmax, 15.0);
        assert_eq!(p.distmin, 0.0);
        assert_eq!(p.atypes.len(), 20);
        assert!(p.atypes.contains("ACA"));
        assert_eq!(p.xvector.len(), 150);
        // 20 CA types → 20*21/2 = 210 unordered pairs.
        assert_eq!(p.table.len(), 210);
        assert_eq!(p.table["ACAACA"].len(), 150);
    }
}
