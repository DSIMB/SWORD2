//! Validation: the pure-Rust scorer must agree with the C++ `scoring_omp` binary.
//!
//! Raw pseudo-energy is deterministic in both implementations → assert near-exact
//! agreement. The Z-score depends on random decoy shuffles (C++ seeds with
//! `time(NULL)`, Rust uses a fixed seed) → assert agreement within a statistical
//! tolerance only.
//!
//! Gated on the C++ binary being present: builds without it (energy is now pure
//! Rust) simply skip this test. Build the binary with `make -C bin/mypmfs-master`.

use std::path::Path;
use std::process::Command;

use sword2_lib::energy::{get_energy_and_z_score, EnergyConfig};

const SHUFFLES: usize = 500;

fn parse_field(stdout: &str, prefix: &str) -> Option<f64> {
    stdout
        .lines()
        .find_map(|l| l.strip_prefix(prefix)?.trim().parse::<f64>().ok())
}

#[test]
fn rust_scorer_agrees_with_cpp_scoring_omp() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let bin_dir = format!("{manifest}/../bin");
    let scoring_bin = format!("{bin_dir}/mypmfs-master/scoring_omp");
    if !Path::new(&scoring_bin).exists() {
        eprintln!("skip: {scoring_bin} not built (run `make -C bin/mypmfs-master`)");
        return;
    }
    let fixture = format!("{manifest}/tests/fixtures/1jx4_ca.pdb");

    let mut config = EnergyConfig::from_bin_dir(&bin_dir);
    config.num_shuffles = SHUFFLES;

    // Reference: C++ scoring_omp -z.
    let output = Command::new(&scoring_bin)
        .args(["-i", &fixture, "-d", &config.potential_dir, "-z", "-s"])
        .arg(SHUFFLES.to_string())
        .output()
        .expect("run scoring_omp");
    assert!(
        output.status.success(),
        "scoring_omp failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    let cpp_energy = parse_field(&stdout, "Pseudo-energy = ").expect("C++ energy");
    let cpp_z = parse_field(&stdout, "Z-score = ").expect("C++ z-score");

    // Rust scorer via the public API.
    let rust = get_energy_and_z_score(&config, &fixture, None).expect("rust score");
    let rust_energy = rust.energy.expect("rust energy");
    let rust_z = rust.z_score.expect("rust z-score");

    // Raw energy: deterministic → near-exact (a tiny relative tolerance absorbs
    // the [last_bin, distmax) interpolation-boundary divergence; see
    // KNOWN_DIVERGENCES.md).
    let rel = (rust_energy - cpp_energy).abs() / cpp_energy.abs().max(1.0);
    assert!(
        rel < 1e-3,
        "energy mismatch: rust={rust_energy}, cpp={cpp_energy} (rel={rel:.2e})"
    );

    // Z-score: different RNG → statistical agreement only.
    assert!(
        (rust_z - cpp_z).abs() < 0.5,
        "z-score mismatch beyond tolerance: rust={rust_z}, cpp={cpp_z}"
    );
}
