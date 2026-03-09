# SWORD2: SWift and Optimized Recognition of protein Domains

[![Docker Pulls](https://img.shields.io/docker/pulls/dsimb/sword2.svg)](https://hub.docker.com/r/dsimb/sword2)
![GitHub Release](https://img.shields.io/github/v/release/DSIMB/SWORD2)

The SWORD2 partitioning algorithm produces multiple alternative domain assignments for a given protein structure. This unique approach handles ambiguous protein structure partitioning, admitting several solutions. The decomposition of the protein structure into domains is achieved through the hierarchical clustering of Protein Units, evolutionarily preserved structural descriptors at the interface between secondary structures and domains.

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://user-images.githubusercontent.com/25644865/181113256-4c2e9740-014f-4d57-91d0-f5beaf7d51d3.png" width="300">
  <img alt="" src="https://user-images.githubusercontent.com/25644865/181106191-e97f1ace-fb93-41d1-a4f0-6a84b7fcc2a1.png" width="300">
</picture>
</p>


## Webserver

This is the repository of the standalone version of the corresponding webserver:  
[dsimb.inserm.fr/SWORD2](https://dsimb.inserm.fr/SWORD2)


## Publications

[Cretin, G., Galochkina, T., Vander Meersche, Y., de Brevern, A. G., Postic, G., & Gelly, J. C. (2022).
SWORD2: hierarchical analysis of protein 3D structures. Nucleic acids research, gkac370.
50(W1), W732–W738 10.1093/nar/gkac370](https://doi.org/10.1093/nar/gkac370)

[Postic, G., Ghouzam, Y., Chebrek, R., & Gelly, J. C. (2017).
An ambiguity principle for assigning protein structural domains.
Science advances, 3(1), e1600552.10.1126/sciadv.1600552](https://doi.org/10.1126/sciadv.1600552)



## Install on Linux, macOS (Intel & Apple Silicon) and Windows

SWORD2 is implemented in Rust. It requires [Rust](https://www.rust-lang.org/tools/install) (1.70+).

On macOS, you also need the Xcode Command Line Tools with the license accepted:
```bash
xcode-select --install        # install if missing
sudo xcodebuild -license      # accept the license
```

```bash
# Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Build the Rust binary
cd sword2-rs
cargo build --release

# The binary is at sword2-rs/target/release/sword2
```

You still need to compile the native dependencies (SWORD, DSSP, MyPMFs, Peeling):
```bash
bash install.sh
```

### Docker

A Docker image is also available:
```bash
docker pull dsimb/sword2:latest
```

## How to use SWORD2

### Docker

```bash
# This will mount the local `results` directory inside the `/output` Docker image directory and run SWORD2 as you (local user)
# so that the results files are owned by you and not root.
(sudo) docker run --rm -e USER_ID=$(id -u) -e GROUP_ID=$(id -g) -v $(pwd)/results:/output dsimb/sword2:latest -p 1jx4 -o /output
```

### Local

Run from the repository root:

```bash
# On a PDB id:
./sword2-rs/target/release/sword2 -p 1jx4 -o results

# On an AlphaFold predicted structure using its Uniprot Accession Id:
./sword2-rs/target/release/sword2 -u Q76EI6 -o results

# On an ESMFold predicted structure using its MGnify Id from the ESM Metagenomic Atlas:
./sword2-rs/target/release/sword2 -m MGYP000936678158 -o results

# On your own PDB/mmCIF structure:
./sword2-rs/target/release/sword2 -i ./structure.pdb -o results

# On your own PDB/mmCIF of an NMR structure, model 2:
./sword2-rs/target/release/sword2 -i ./structure.pdb -d 2 -o results

# Specify the base directory explicitly (if not running from the repo root):
./sword2-rs/target/release/sword2 -p 1jx4 -o results --base-dir /path/to/SWORD2
```

Note: The Rust version currently skips plot generation (contact probability matrices, domain histograms).

## Fast mode 

#### Skip calculation of pseudo-energy (2x faster)
```bash
./sword2-rs/target/release/sword2 -p 1jx4 -o results --disable-energies
```

## Memory usage

SWORD2 will have a peak memory usage of 150-200 Mb for average protein sizes of 100-300 residues.
The largest proteins, for example 2700 residues can take up to 1.5 Gb memory.

## Help

To get the full help:
```console
$ ./sword2-rs/target/release/sword2 --help
Usage: sword2 [OPTIONS] --output <OUTPUT>

SWORD2: SWift and Optimized Recognition of protein Domains.
The SWORD2 partitioning algorithm produces multiple alternative
domain assignments for a given protein structure.
This unique approach handles ambiguous protein structure partitioning,
admitting several solutions. The decomposition of the protein structure
into domains is achieved through the hierarchical clustering of Protein Units,
evolutionarily preserved structural descriptors at the interface between
secondary structures and domains.

Options:
  -u, --uniprot-id <UNIPROT_ID>
          AlphaFold Uniprot Accession Id.
  -m, --mgnify-id <MGNIFY_ID>
          MGnify Id for the ESM Metagenomic Atlas.
  -p, --pdb-id <PDB_ID>
          PDB id to download from the PDB database.
  -i, --input-file <INPUT_FILE>
          Path to an input PDB or mmCIF file.
  -o, --output <OUTPUT>
          Output directory. Results will be generated inside in a dedicated directory named after OUTPUT/PDBID_CHAIN/
  -c, --pdb-chain <PDB_CHAIN>
          PDB chain. If not specified, the first chain in the PDB file will be used.
  -d, --model <MODEL>
          Model to parse. Especially useful for NMR files which contain several models. Default is 1. [default: 1]
  -x, --cpu <CPU>
          Number of CPUs to use. Default all (0).
  -e, --disable-energies
          Disable the calculation of pseudo-energy of domains and PUs.
      --disable-plots
          (Ignored in Rust port currently) Disable the generation of contact probability matrices plots.
  -b, --base-dir <BASE_DIR>
          Base directory string. Usually not needed if running from the source repository directory or installed with pip setup.
  -t, --tmp-dir <TMP_DIR>
          Change tmp directory name if you are running multiple SWORD instances in the same current working directory. [default: SWORD2_tmp]
  -h, --help
          Print help (see more with '--help')
  -V, --version
          Print version
```

## Output

Example:
```
$ ./sword2-rs/target/release/sword2 -p 1jx4 -o results
2024-09-19T11:14:27Z INFO  [sword2] Fetch PDB ID: 1jx4
2024-09-19T11:14:28Z INFO  [sword2] No chain specified. Using first chain 'A' in the PDB file.
2024-09-19T11:14:28Z INFO  [sword2] 
2024-09-19T11:14:28Z INFO  [sword2] >>>   1jx4_A (335 aa)
2024-09-19T11:14:28Z INFO  [sword2] >>>   Estimated runtime: 1 minutes and 6 seconds
2024-09-19T11:14:28Z INFO  [sword2] >>>   Using 10 cpus
2024-09-19T11:14:28Z INFO  [sword2] 
2024-09-19T11:14:28Z INFO  [sword2] Write a clean version of the PDB: remove non standard residues
2024-09-19T11:14:28Z INFO  [sword2] Launch SWORD
2024-09-19T11:14:34Z INFO  [sword2] Parse SWORD output
2024-09-19T11:14:34Z INFO  [sword2] Calculate pseudo-energies of Domains
2024-09-19T11:15:27Z INFO  [sword2] Write the SWORD results
2024-09-19T11:15:27Z INFO  [sword2] Write Peeling results
2024-09-19T11:15:54Z INFO  [sword2] Calculate junctions consistencies
2024-09-19T11:15:54Z INFO  [sword2] Clean and prepare results
2024-09-19T11:15:54Z INFO  [sword2] Results can be found here: results/1jx4_A
2024-09-19T11:15:54Z INFO  [sword2] Total runtime: 87 seconds
```

An easily parseable output in JSON format is generated for easier downstream tasks/analysis: `SWORD2_summary.json`
