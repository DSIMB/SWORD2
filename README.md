# SWORD2: SWift and Optimized Recognition of protein Domains

[![DOI](https://zenodo.org/badge/518046197.svg)](https://zenodo.org/badge/latestdoi/518046197)
[![Docker Pulls](https://img.shields.io/docker/pulls/dsimb/sword2.svg)](https://hub.docker.com/r/dsimb/sword2)


The SWORD2 partitioning algorithm produces multiple alternative domain assignments for a given protein structure. This unique approach handles ambiguous protein structure partitioning, admitting several solutions. The decomposition of the protein structure into domains is achieved through the hierarchical clustering of Protein Units, evolutionarily preserved structural descriptors at the interface between secondary structures and domains.

<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://user-images.githubusercontent.com/25644865/181113256-4c2e9740-014f-4d57-91d0-f5beaf7d51d3.png" width="300">
  <img alt="" src="https://user-images.githubusercontent.com/25644865/181106191-e97f1ace-fb93-41d1-a4f0-6a84b7fcc2a1.png" width="300">
</picture>
</p>

This is the repository of the standalone version of the corresponding webserver:  
https://www.dsimb.inserm.fr/SWORD2/index.html



## Publications

[Cretin, G., Galochkina, T., Vander Meersche, Y., de Brevern, A. G., Postic, G., & Gelly, J. C. (2022).
SWORD2: hierarchical analysis of protein 3D structures. Nucleic acids research, gkac370.
50(W1), W732–W738 10.1093/nar/gkac370](https://doi.org/10.1093/nar/gkac370)

[Postic, G., Ghouzam, Y., Chebrek, R., & Gelly, J. C. (2017).
An ambiguity principle for assigning protein structural domains.
Science advances, 3(1), e1600552.10.1126/sciadv.1600552](https://doi.org/10.1126/sciadv.1600552)



## Install (Linux, macOS <= 13 Intel, M1/M2, Windows)

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
docker pull dsimb/sword2
```

## How to use SWORD2

Launch the docker version:
```
# This will mount the local `results` directory inside the Docker image and run the program inside the image as you (local user)
# so that the results generated have good rights.
docker run -it -e LOCAL_UID=$(id -u $USER) -e LOCAL_GID=$(id -g $USER) -v $(pwd)/results:/app/results dsimb/sword2 -p 1jx4 -o results
```

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

## Help

To get the full help:
```
$ ./SWORD2.py --help                                                                                                                                                        ─╯
usage: SWORD2.py [-h] (-u UNIPROT_ID | -m MGNIFY_ID | -p PDB_ID | -i INPUT_FILE) [-c PDB_CHAIN] [-d MODEL] [-x CPU] -o OUTPUT

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
  -u UNIPROT_ID, --uniprot-id UNIPROT_ID
                        AlphaFold Uniprot Accession Id.
                        The corresponding predicted structure will be downloaded from the AlphaFold database.
  -m MGNIFY_ID, --mgnify-id MGNIFY_ID
                        MGnify Id.
                        The corresponding predicted structure will be downloaded from the ESM Metagenomic Atlas database.
  -p PDB_ID, --pdb-id PDB_ID
                        PDB id.
                        The corresponding structure will be downloaded from the PDB database.
  -i INPUT_FILE, --input-file INPUT_FILE
                        Path to an input PDB or mmCIF file.

optional arguments:
  -c PDB_CHAIN, --pdb-chain PDB_CHAIN
                        PDB chain. Default is A.
  -d MODEL, --model MODEL
                        Model to parse. Especially usefull for NMR files which contain several models. Default is 1.
  -x CPU, --cpu CPU     How many CPUs to use.
                        Default all (0).
                        Max on this computer is: 32

required arguments:
  -o OUTPUT, --output OUTPUT
                        Output directory.
                        Results will be generated inside in a dedicated directory
                        named after OUTPUT/PDBID_CHAIN/
```


## Output

Example:
```
$ ./SWORD2.py -p 1jx4 -o results                                                                                                                                            ─╯
2022/11/03 02:00:04 INFO     Fetch PDB ID: 1jx4
2022/11/03 02:00:04 DEBUG    Connecting wwPDB FTP server RCSB PDB (USA).
2022/11/03 02:00:08 DEBUG    1jx4 downloaded (results/1jx4_A/1jx4.pdb)
2022/11/03 02:00:08 DEBUG    PDB download via FTP completed (1 downloaded, 0 failed).
2022/11/03 02:00:08 DEBUG    3162 atoms and 1 coordinate set(s) were parsed in 0.03s.
2022/11/03 02:00:08 INFO
2022/11/03 02:00:08 INFO     >>>   Estimated runtime: 4 minutes
2022/11/03 02:00:08 INFO     >>>   Using 8 cpus
2022/11/03 02:00:08 INFO
2022/11/03 02:00:08 INFO     Write a clean version of the PDB: remove non standard residues
2022/11/03 02:00:08 INFO     Launch SWORD
2022/11/03 02:00:23 INFO     Parse SWORD output
2022/11/03 02:00:23 INFO     Calculate pseudo-energies of Domains
2022/11/03 02:01:57 INFO     Write the SWORD results
2022/11/03 02:01:57 INFO     Calculate pseudo-energies of PUs
2022/11/03 02:02:37 INFO     Write the Peeling results
2022/11/03 02:02:37 INFO     Generate histogram of SWORD2 domains consistency
2022/11/03 02:03:09 INFO     Calculate junctions consistencies
2022/11/03 02:03:09 INFO     Clean and prepare results
2022/11/03 02:03:09 INFO     Results can be found here: results/1jx4_A
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
