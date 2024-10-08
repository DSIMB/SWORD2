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

(Easy & recommanded)  
Install the conda environment using the `environment.yml` file:
```
conda env create -f environment.yml
or
mamba env create -f environment.yml
```

Compile necessary dependencies:
```
bash install.sh
```

Otherwise, a docker is also available:
```
docker pull dsimb/sword2:latest
```

## How to use SWORD2

### Docker

```
# This will mount the local `results` directory inside the `/output` Docker image directory and run SWORD2 as you (local user)
# so that the results files are owned by you and not root.
(sudo) docker run --rm -e USER_ID=$(id -u) -e GROUP_ID=$(id -g) -v $(pwd)/results:/output dsimb/sword2:latest -p 1jx4 -o /output
```

### Conda
First, activate the working environment:
```
conda activate sword2
```

Then, launch SWORD2:

#### On a PDB id:
```
./SWORD2.py -p 1jx4 -o results
```

#### On an AlphaFold predicted structure using its Uniprot Accession Id:
```
./SWORD2.py -u Q76EI6 -o results
```

#### On an ESMFold predicted structure using its MGnify Id from the ESM Metagenomic Atlas:
```
./SWORD2.py -m MGYP000936678158 -o results
```

#### On your own PDB/mmCIF structure:
```
./SWORD2.py -i ./structure.pdb -o results
```

#### On your own PDB/mmCIF of an NMR structure, model 2:
```
./SWORD2.py -i ./structure.pdb -d 2 -o results
```

## Fast mode 

#### Skip calculation of pseudo-energy and plots of contact probability maps (2x faster)
```
./SWORD2.py -p 1jx4 -o results --disable-energy --disable-plots
```

## Memory usage

SWORD2 will have a peak memory usage of 150-200 Mb for average protein sizes of 100-300 residues.
The largest proteins, for example 2700 residues can take up to 1.5 Gb memory.

## Help

To get the full help:
```console
$ ./SWORD2.py --help
usage: SWORD2.py [-h] [--version] (-u UNIPROT_ID | -m MGNIFY_ID | -p PDB_ID | -i INPUT_FILE) [-c PDB_CHAIN] [-d MODEL] [-x CPU] [-e] [-l] -o OUTPUT

SWORD2: SWift and Optimized Recognition of protein Domains.
The SWORD2 partitioning algorithm produces multiple alternative
domain assignments for a given protein structure.
This unique approach handles ambiguous protein structure partitioning,
admitting several solutions. The decomposition of the protein structure
into domains is achieved through the hierarchical clustering of Protein Units,
evolutionarily preserved structural descriptors at the interface between
secondary structures and domains.

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -u UNIPROT_ID, --uniprot-id UNIPROT_ID
                        AlphaFold Uniprot Accession Id.
  -m MGNIFY_ID, --mgnify-id MGNIFY_ID
                        MGnify Id for the ESM Metagenomic Atlas.
  -p PDB_ID, --pdb-id PDB_ID
                        PDB id to download from the PDB database.
  -i INPUT_FILE, --input-file INPUT_FILE
                        Path to an input PDB or mmCIF file.

optional arguments:
  -c PDB_CHAIN, --pdb-chain PDB_CHAIN
                        PDB chain. If not specified, the first chain in the PDB file will be used.
  -d MODEL, --model MODEL
                        Model to parse. Especially useful for NMR files which contain several models. Default is 1.
  -x CPU, --cpu CPU     Number of CPUs to use. Default all (0). Max on this computer is: 10
  -e, --disable-energies
                        Disable the calculation of pseudo-energy of domains and PUs.
  -l, --disable-plots   Disable the generation of contact probability matrices plots.

required arguments:
  -o OUTPUT, --output OUTPUT
                        Output directory. Results will be generated inside in a dedicated directory named after OUTPUT/PDBID_CHAIN/
```


## Output

Example:
```
$ conda activate sword2
(sword2) $ ./SWORD2.py -p 1jx4 -o results
2024/09/19 11:14:27 INFO     Fetch PDB ID: 1jx4
2024/09/19 11:14:28 INFO     No chain specified. Using first chain 'A' in the PDB file.
2024/09/19 11:14:28 INFO
2024/09/19 11:14:28 INFO     >>>   1jx4_A (335 aa)
2024/09/19 11:14:28 INFO     >>>   Estimated runtime: 1 minutes and 6 seconds
2024/09/19 11:14:28 INFO     >>>   Using 10 cpus
2024/09/19 11:14:28 INFO
2024/09/19 11:14:28 INFO     Write a clean version of the PDB: remove non standard residues
2024/09/19 11:14:28 INFO     Launch SWORD
2024/09/19 11:14:34 INFO     Parse SWORD output
2024/09/19 11:14:34 INFO     Calculate pseudo-energies of Domains
2024/09/19 11:15:27 INFO     Write the SWORD results
2024/09/19 11:15:27 INFO     Write Peeling results
2024/09/19 11:15:48 INFO     Generate histogram of SWORD2 domains consistency
2024/09/19 11:15:49 INFO     Generate contact probability matrices
2024/09/19 11:15:54 INFO     Calculate junctions consistencies
2024/09/19 11:15:54 INFO     Clean and prepare results
2024/09/19 11:15:54 INFO     Results can be found here: results/1jx4_A
2024/09/19 11:15:54 INFO     Total runtime: 87 seconds
```

An easily parseable output in JSON format is generated for easier downstream tasks/analysis: `SWORD2_summary.json`



## Multithreading

Users can now run multiple SWORD2 jobs in parallel.
Here is an example script:
```python
import subprocess
from concurrent.futures import ThreadPoolExecutor

# List of 5 PDB codes
pdb_list = ['1jx4', '2c78', '1f5n', '1a8y', '1b89']

# Function to run a SWORD2 job
def run_sword2(pdb_code):
    command = f'./SWORD2.py -p {pdb_code} -o results'
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Job for {pdb_code} completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running job for {pdb_code}: {e}")

# Run 5 jobs using ThreadPoolExecutor
def run_sword2_jobs():
    with ThreadPoolExecutor(max_workers=5) as executor:  # Adjust max_workers based on your system's capabilities
        executor.map(run_sword2, pdb_list)

if __name__ == "__main__":
    run_sword2_jobs()
```
