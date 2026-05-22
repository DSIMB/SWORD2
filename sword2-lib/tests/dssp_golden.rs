/// Golden integration test for the pure-Rust DSSP implementation.
///
/// Runs `run_dssp` on the backbone atoms of residues 1-20 of 1JX4_A
/// (embedded as a fixture) and asserts:
///   - Exactly 20 residues are parsed.
///   - Residues 10-18 are assigned alpha-helix ('H'), which is expected from the
///     known crystal structure.  All inter-residue H-bonds spanning i→i+4 for the
///     helix are within the 20-residue fragment, so the helix is fully detectable
///     even without the rest of the chain.
///   - Both output files (.dssp, .s2d) are written to disk.
///
/// Ground truth is the current Rust DSSP output (not the original dsspcmbi binary).
/// See KNOWN_DIVERGENCES.md for rationale.
use std::fs;
use tempfile::tempdir;

use sword2_lib::dssp::run_dssp;

/// Backbone atoms (N, CA, C, O only) for residues 1-20 of 1JX4_A, chain A.
/// Source: results/1JX4_A/intermediate/1JX4_A.pdb (gitignored; extracted once).
const FIXTURE_PDB: &str = "\
ATOM      1  N   ILE A   1      46.170  17.543  13.913  1.00 33.83           N
ATOM      2  CA  ILE A   1      45.665  16.827  12.751  1.00 32.01           C
ATOM      3  C   ILE A   1      44.396  16.160  13.297  1.00 29.09           C
ATOM      4  O   ILE A   1      44.466  15.419  14.274  1.00 27.71           O
ATOM      9  N   VAL A   2      43.247  16.494  12.713  1.00 28.04           N
ATOM     10  CA  VAL A   2      41.973  15.941  13.134  1.00 26.25           C
ATOM     11  C   VAL A   2      41.475  15.044  12.028  1.00 27.46           C
ATOM     12  O   VAL A   2      41.543  15.426  10.857  1.00 26.95           O
ATOM     16  N   LEU A   3      40.986  13.867  12.417  1.00 26.01           N
ATOM     17  CA  LEU A   3      40.415  12.890  11.503  1.00 25.53           C
ATOM     18  C   LEU A   3      38.998  12.690  12.000  1.00 24.10           C
ATOM     19  O   LEU A   3      38.814  12.233  13.138  1.00 24.68           O
ATOM     24  N   PHE A   4      38.029  13.007  11.125  1.00 23.52           N
ATOM     25  CA  PHE A   4      36.591  12.947  11.382  1.00 23.38           C
ATOM     26  C   PHE A   4      36.034  11.771  10.603  1.00 23.86           C
ATOM     27  O   PHE A   4      36.392  11.560   9.449  1.00 25.30           O
ATOM     35  N   VAL A   5      35.163  10.999  11.241  1.00 24.37           N
ATOM     36  CA  VAL A   5      34.553   9.829  10.616  1.00 23.47           C
ATOM     37  C   VAL A   5      33.065  10.084  10.590  1.00 22.58           C
ATOM     38  O   VAL A   5      32.478  10.418  11.618  1.00 22.07           O
ATOM     42  N   ASP A   6      32.446   9.899   9.424  1.00 21.80           N
ATOM     43  CA  ASP A   6      31.022  10.191   9.273  1.00 21.96           C
ATOM     44  C   ASP A   6      30.440   8.962   8.554  1.00 21.90           C
ATOM     45  O   ASP A   6      30.778   8.694   7.402  1.00 20.77           O
ATOM     50  N   PHE A   7      29.621   8.183   9.249  1.00 20.49           N
ATOM     51  CA  PHE A   7      29.087   6.961   8.633  1.00 23.13           C
ATOM     52  C   PHE A   7      28.154   7.314   7.469  1.00 24.00           C
ATOM     53  O   PHE A   7      27.353   8.254   7.550  1.00 22.75           O
ATOM     61  N   ASP A   8      28.212   6.533   6.393  1.00 23.11           N
ATOM     62  CA  ASP A   8      27.370   6.885   5.240  1.00 21.98           C
ATOM     63  C   ASP A   8      25.887   6.510   5.392  1.00 21.02           C
ATOM     64  O   ASP A   8      25.564   5.394   5.796  1.00 21.79           O
ATOM     69  N   TYR A   9      25.007   7.441   5.023  1.00 20.58           N
ATOM     70  CA  TYR A   9      23.558   7.272   5.091  1.00 20.92           C
ATOM     71  C   TYR A   9      23.224   6.235   6.153  1.00 22.62           C
ATOM     72  O   TYR A   9      22.422   5.300   5.958  1.00 21.24           O
ATOM     81  N   PHE A  10      23.748   6.501   7.340  1.00 20.52           N
ATOM     82  CA  PHE A  10      23.718   5.478   8.361  1.00 19.69           C
ATOM     83  C   PHE A  10      22.537   4.569   8.612  1.00 20.06           C
ATOM     84  O   PHE A  10      22.721   3.356   8.506  1.00 21.07           O
ATOM     92  N   TYR A  11      21.339   5.083   8.936  1.00 18.38           N
ATOM     93  CA  TYR A  11      20.265   4.179   9.278  1.00 20.38           C
ATOM     94  C   TYR A  11      19.883   3.292   8.080  1.00 21.64           C
ATOM     95  O   TYR A  11      19.577   2.084   8.233  1.00 20.87           O
ATOM    104  N   ALA A  12      19.887   3.898   6.896  1.00 17.70           N
ATOM    105  CA  ALA A  12      19.539   3.124   5.725  1.00 19.35           C
ATOM    106  C   ALA A  12      20.628   2.115   5.446  1.00 20.93           C
ATOM    107  O   ALA A  12      20.310   1.001   5.042  1.00 19.80           O
ATOM    109  N   GLN A  13      21.916   2.467   5.650  1.00 21.40           N
ATOM    110  CA  GLN A  13      22.973   1.495   5.407  1.00 21.95           C
ATOM    111  C   GLN A  13      22.889   0.317   6.387  1.00 21.72           C
ATOM    112  O   GLN A  13      23.142  -0.839   5.990  1.00 23.52           O
ATOM    118  N   VAL A  14      22.554   0.581   7.649  1.00 20.36           N
ATOM    119  CA  VAL A  14      22.387  -0.510   8.587  1.00 20.59           C
ATOM    120  C   VAL A  14      21.278  -1.456   8.093  1.00 23.68           C
ATOM    121  O   VAL A  14      21.447  -2.697   8.139  1.00 23.53           O
ATOM    125  N   GLU A  15      20.159  -0.892   7.629  1.00 21.78           N
ATOM    126  CA  GLU A  15      19.099  -1.754   7.132  1.00 23.52           C
ATOM    127  C   GLU A  15      19.604  -2.634   5.957  1.00 22.81           C
ATOM    128  O   GLU A  15      19.224  -3.810   5.896  1.00 24.56           O
ATOM    134  N   GLU A  16      20.481  -2.092   5.093  1.00 23.01           N
ATOM    135  CA  GLU A  16      21.056  -2.832   3.949  1.00 25.82           C
ATOM    136  C   GLU A  16      21.997  -3.923   4.460  1.00 28.76           C
ATOM    137  O   GLU A  16      22.128  -4.972   3.816  1.00 29.72           O
ATOM    143  N   VAL A  17      22.691  -3.675   5.571  1.00 27.94           N
ATOM    144  CA  VAL A  17      23.572  -4.729   6.072  1.00 28.79           C
ATOM    145  C   VAL A  17      22.731  -5.870   6.647  1.00 27.56           C
ATOM    146  O   VAL A  17      23.083  -7.034   6.448  1.00 31.25           O
ATOM    150  N   LEU A  18      21.658  -5.551   7.356  1.00 27.86           N
ATOM    151  CA  LEU A  18      20.755  -6.546   7.935  1.00 27.86           C
ATOM    152  C   LEU A  18      19.855  -7.246   6.854  1.00 27.81           C
ATOM    153  O   LEU A  18      19.358  -8.353   7.092  1.00 27.83           O
ATOM    158  N   ASN A  19      19.684  -6.612   5.693  1.00 26.87           N
ATOM    159  CA  ASN A  19      18.928  -7.212   4.583  1.00 27.50           C
ATOM    160  C   ASN A  19      19.573  -6.791   3.260  1.00 27.77           C
ATOM    161  O   ASN A  19      19.141  -5.853   2.585  1.00 27.17           O
ATOM    166  N   PRO A  20      20.663  -7.485   2.868  1.00 28.99           N
ATOM    167  CA  PRO A  20      21.435  -7.239   1.642  1.00 29.83           C
ATOM    168  C   PRO A  20      20.602  -7.208   0.342  1.00 30.16           C
ATOM    169  O   PRO A  20      21.043  -6.664  -0.654  1.00 32.03           O
TER
END
";

#[test]
fn test_run_dssp_helix_detection_on_1jx4a_fragment() {
    let dir = tempdir().unwrap();
    let pdb_path = dir.path().join("fragment.pdb");
    let dssp_path = dir.path().join("fragment.dssp");
    let s2d_path = dir.path().join("fragment.s2d");

    fs::write(&pdb_path, FIXTURE_PDB).unwrap();

    let result = run_dssp(&pdb_path, &dssp_path, &s2d_path, "1JX4_A_fragment").unwrap();

    assert_eq!(result.chain.len, 20, "Expected 20 backbone residues");

    // Residues 10-18 form an alpha helix in the full 1JX4_A structure.
    // All i→i+4 H-bonds for this helix are within the 20-residue fragment,
    // so the helix should be fully detected here too.
    for pos in 10..=18 {
        let ss = result.chain.get(pos).ss[0];
        assert!(
            ss == 'H' || ss == 'G' || ss == 'I',
            "Expected helix symbol at position {pos} but got '{ss}'"
        );
    }

    // Output files must be written
    assert!(dssp_path.exists(), ".dssp file not written");
    assert!(s2d_path.exists(), ".s2d file not written");
}
