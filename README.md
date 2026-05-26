# SWORD2: SWift and Optimized Recognition of protein Domains

[![Docker Pulls](https://img.shields.io/docker/pulls/dsimb/sword2.svg)](https://hub.docker.com/r/dsimb/sword2)
![GitHub Release](https://img.shields.io/github/v/release/DSIMB/SWORD2)

SWORD2 produces multiple alternative domain assignments for a protein structure. Instead of forcing a single decomposition, it preserves ambiguity when the structure supports several plausible partitions.

This repository is the standalone Rust implementation of SWORD2. The current codebase is organized as a Rust workspace with a CLI crate and a shared library crate.

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://user-images.githubusercontent.com/25644865/181113256-4c2e9740-014f-4d57-91d0-f5beaf7d51d3.png" width="300">
  <img alt="SWORD2 logo" src="https://user-images.githubusercontent.com/25644865/181106191-e97f1ace-fb93-41d1-a4f0-6a84b7fcc2a1.png" width="300">
</picture>
</p>

## Highlights

- Pure Rust pipeline orchestration.
- Pure Rust DSSP implementation.
- Pure Rust Peeling implementation.
- Supports local structures, PDB IDs, AlphaFold UniProt accessions, and ESM Metagenomic Atlas MGnify IDs.
- Writes both human-readable summaries and JSON outputs.
- Optional pseudo-energy scoring from bundled `mypmfs` statistical potential data.

## Web Server

The public SWORD2 web server is available at [dsimb.inserm.fr/SWORD2](https://dsimb.inserm.fr/SWORD2).

## Publications

[Cretin, G., Galochkina, T., Vander Meersche, Y., de Brevern, A. G., Postic, G., & Gelly, J. C. (2022). SWORD2: hierarchical analysis of protein 3D structures. Nucleic Acids Research, 50(W1), W732-W738.](https://doi.org/10.1093/nar/gkac370)

[Postic, G., Ghouzam, Y., Chebrek, R., & Gelly, J. C. (2017). An ambiguity principle for assigning protein structural domains. Science Advances, 3(1), e1600552.](https://doi.org/10.1126/sciadv.1600552)

## Repository Layout

- `sword2-cli/`: CLI crate. Entry point for the standalone executable.
- `sword2-lib/`: core library with parsing, DSSP, peeling, energy, plotting, and output modules.
- `sword2-lib/assets/fonts/`: embedded OFL-licensed font used for plot labels.
- `bin/mypmfs-master/025_30_100_potential/`: pseudo-energy statistical potential data.
- `results/`: example outputs.

## Requirements

- Rust toolchain for building the workspace.

The current code uses Rust edition 2021 and builds with standard Cargo commands.

## Build

From the repository root:

```bash
cargo build --release
```

The executable will be written to:

```text
target/release/sword2
```

For convenience, the repository root also includes a launcher script so users can run:

```bash
./sword2 --help
```

The launcher forwards to `target/release/sword2` when available, falls back to `target/debug/sword2`, and prints a short build hint if the project has not been compiled yet.

For development builds:

```bash
cargo build
```

To run the test suite:

```bash
cargo test
```

## Energy Backend

Domain and PU pseudo-energies are computed in Rust (`sword2-lib/src/energy/score.rs`)
using the precomputed `mypmfs` statistical potentials.

The runtime expects the potential data to exist relative to the SWORD2 base directory:

```text
bin/mypmfs-master/025_30_100_potential/
```

If you do not need pseudo-energies, run with `--disable-energies`.

## Usage

The CLI accepts one input source per run:

- `--pdb-id` for a structure from the PDB.
- `--uniprot-id` for an AlphaFold model.
- `--mgnify-id` for an ESM Metagenomic Atlas model.
- `--input-file` for a local PDB or mmCIF file.

Basic examples:

```bash
# PDB structure
./sword2 -p 1jx4 -o results

# Local PDB or mmCIF file
./sword2 -i structure.pdb -o results

# AlphaFold model by UniProt accession
./sword2 -u Q5VSL9 -o results

# ESM Metagenomic Atlas model by MGnify ID
./sword2 -m MGYP000936678158 -o results
```

Useful variants:

```bash
# Analyze a specific chain
./sword2 -p 1jx4 -c A -o results

# Use a different model for NMR structures
./sword2 -i structure.pdb --model 2 -o results

# Skip pseudo-energy calculations
./sword2 -p 1jx4 -o results --disable-energies

# Skip plot generation
./sword2 -p 1jx4 -o results --disable-plots

# Reduce Z-score shuffle count for faster energy runs
./sword2 -p 1jx4 -o results --num-shuffles 500

# Increase verbosity
./sword2 -p 1jx4 -o results -v
./sword2 -p 1jx4 -o results -vv

# Quiet mode
./sword2 -p 1jx4 -o results -q
```

If you run the binary outside the repository root, point it to the project base directory so it can find `bin/`:

```bash
./sword2 -p 1jx4 -o results --base-dir /path/to/SWORD2
```

You can also run directly through Cargo:

```bash
cargo run -p sword --release -- -p 1jx4 -o results
```

## CLI Summary

Current top-level options:

- `-p, --pdb-id <PDB_ID>`: fetch a structure from the PDB.
- `-c, --chain <CHAIN>`: choose the chain to analyze.
- `-u, --uniprot-id <UNIPROT_ID>`: fetch an AlphaFold model.
- `-m, --mgnify-id <MGNIFY_ID>`: fetch an ESM Atlas model.
- `-i, --input-file <INPUT_FILE>`: analyze a local PDB or mmCIF file.
- `-o, --output-dir <OUTPUT_DIR>`: choose the output directory. Defaults to the current directory.
- `--model <MODEL>`: pick a structure model for NMR inputs. Defaults to `1`.
- `-e, --disable-energies`: skip pseudo-energy calculations.
- `-l, --disable-plots`: skip contact matrix plot generation.
- `-x, --cpu <CPU>`: set the number of worker threads. `0` means all CPUs.
- `--base-dir <BASE_DIR>`: locate the SWORD2 base directory containing `bin/`.
- `-s, --num-shuffles <NUM_SHUFFLES>`: adjust Z-score shuffle count. Defaults to `2000`.
- `-v, --verbose...`: increase logging verbosity.
- `-q, --quiet`: suppress non-error output.

For the exact CLI help of the binary you built:

```bash
./sword2 --help
```

## Output Structure

A run creates one directory per analyzed chain, named `<INPUT>_<CHAIN>` inside the selected output directory.

Example:

```bash
./sword2 -p 1jx4 -o results
```

This produces a directory like:

```text
results/1JX4_A/
  input.pdb
  summary.txt
  summary.json
  peeling.txt
  peeling.json
  junctions.txt
  residue_mapping.txt
  intermediate/
    1JX4_A.dssp
    1JX4_A.pdb
    contact_matrix.mat
    peeling.log
    ...
  plots/
    alt0.png
    alt0_dom0.png
    alt0_dom0_pu_1_163.png
    ...
    domain_histogram.svg
```

The `plots/` directory is only written when `--disable-plots` is not set. It contains one overview PNG per alternative partitioning, one PNG per domain, and one PNG per protein unit, plus a `domain_histogram.svg` consistency histogram.

## Output Files

- `summary.txt`: human-readable summary of the optimal and alternative partitions.
- `summary.json`: machine-readable JSON summary of the same partitionings.
- `peeling.txt`: text summary of peeling levels and protein units.
- `peeling.json`: JSON version of peeling results.
- `junctions.txt`: junction consistency report.
- `residue_mapping.txt`: mapping from cleaned residue numbering back to author residue numbering.
- `input.pdb`: cleaned copy of the analyzed structure.
- `intermediate/`: internal algorithm files (DSSP output, contact matrices, peeling logs).

## Performance Notes

- `--disable-energies` is the main fast mode and skips pseudo-energy calculations.
- `--num-shuffles` lets you trade Z-score precision for speed during energy calculations.
- `--cpu 0` uses all available CPUs.
- Memory usage depends strongly on protein size; large proteins will require substantially more memory than typical 100-300 residue inputs.

## Docker

A prebuilt image is available on Docker Hub:

```bash
docker pull dsimb/sword2:latest
```

Example run:

```bash
docker run --rm \
  -e USER_ID=$(id -u) \
  -e GROUP_ID=$(id -g) \
  -v "$(pwd)/results:/output" \
  dsimb/sword2:latest \
  -p 1jx4 -o /output
```

## Architecture Notes

The current implementation differs from the older standalone releases in a few important ways:

- DSSP is implemented in Rust.
- Peeling is implemented in Rust.
- The CLI is a Rust binary with Clap-based argument parsing.
- Pseudo-energy/Z-score scoring is implemented in Rust, using the precomputed `mypmfs` potentials.

## License

This project is distributed under the terms of the license in `LICENSE`.
