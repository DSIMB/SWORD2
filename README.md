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
- Optional pseudo-energy scoring through the bundled `mypmfs` backend.

## Web Server

The public SWORD2 web server is available at [dsimb.inserm.fr/SWORD2](https://dsimb.inserm.fr/SWORD2).

## Publications

[Cretin, G., Galochkina, T., Vander Meersche, Y., de Brevern, A. G., Postic, G., & Gelly, J. C. (2022). SWORD2: hierarchical analysis of protein 3D structures. Nucleic Acids Research, 50(W1), W732-W738.](https://doi.org/10.1093/nar/gkac370)

[Postic, G., Ghouzam, Y., Chebrek, R., & Gelly, J. C. (2017). An ambiguity principle for assigning protein structural domains. Science Advances, 3(1), e1600552.](https://doi.org/10.1126/sciadv.1600552)

## Repository Layout

- `sword2-cli/`: CLI crate. Entry point for the standalone executable.
- `sword2-lib/`: core library with parsing, DSSP, peeling, energy, plotting, and output modules.
- `bin/mypmfs-master/`: pseudo-energy scoring backend and potentials.
- `results/`: example outputs.

## Requirements

- Rust toolchain for building the workspace.
- A working C toolchain only if you need to rebuild the `mypmfs` scoring binary.

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

Domain and PU pseudo-energies are computed by the external `scoring_omp` binary used by the `mypmfs` backend.

The runtime expects these paths to exist relative to the SWORD2 base directory:

```text
bin/mypmfs-master/scoring_omp
bin/mypmfs-master/025_30_100_potential/
```

If you need to rebuild the scoring backend manually, use the bundled makefile:

```bash
make -C bin/mypmfs-master
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
  1JX4_A
  SWORD2_summary.txt
  SWORD2_summary.json
  PEELING_summary.txt
  PEELING_summary.json
  sword.txt
  mapping_auth_resnums.txt
  Junctions/
    junctions_consistencies.txt
  Protein_Units/
    Peeling.log
  SWORD/
    1JX4_A/
      ...
```

When plots are enabled, SWORD2 also writes a `Contact_Probability_Matrix/` directory containing:

- `domain_consistency_histogram.svg`
- `contact_probability_matrix_alternative_<N>.png`
- per-domain and per-PU contact matrix PNGs

## Output Files

- `SWORD2_summary.txt`: human-readable summary of the optimal and alternative partitions.
- `SWORD2_summary.json`: machine-readable JSON summary of the same partitionings.
- `PEELING_summary.txt`: text summary of peeling levels and PUs.
- `PEELING_summary.json`: JSON version of peeling results.
- `sword.txt`: raw SWORD pipeline output.
- `mapping_auth_resnums.txt`: mapping from cleaned residue numbering back to author residue numbering.
- `Junctions/junctions_consistencies.txt`: junction consistency report.

## Performance Notes

- `--disable-energies` is the main fast mode and avoids the external scoring backend.
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
- The only external runtime scoring component is `mypmfs/scoring_omp` for pseudo-energy calculations.

## License

This project is distributed under the terms of the license in `LICENSE`.
