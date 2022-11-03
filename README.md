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

If necessary, a docker is also available:
```
docker pull dsimb/sword2
```

Compile necessary dependencies:
```
bash install.sh
```

## How to use SWORD2

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

#### On your own PDB/mmCIF structure:
```
./SWORD2.py -i ./structure.pdb -o results
```

## Help

To get the full help:
```
$ ./SWORD2.py --help                                                                                                                                                        ─╯
usage: SWORD2.py [-h] (-u UNIPROT_CODE | -p PDB_CODE | -i INPUT_FILE) [-c PDB_CHAIN] [-x CPU] -o OUTPUT

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
  -u UNIPROT_CODE, --uniprot-code UNIPROT_CODE
                        AlphaFold Uniprot Accession Id.
                        The corresponding structure will be downloaded from the AlphaFold database.
  -p PDB_CODE, --pdb-code PDB_CODE
                        PDB code.
                        The corresponding structure will be downloaded from the PDB database.
  -i INPUT_FILE, --input-file INPUT_FILE
                        Path to an input PDB or mmCIF file.

optional arguments:
  -c PDB_CHAIN, --pdb-chain PDB_CHAIN
                        PDB chain. Default is A.
  -x CPU, --cpu CPU     How many CPUs to use.
                        Default all (0).
                        Max on this computer is: 8

required arguments:
  -o OUTPUT, --output OUTPUT
                        Output directory.
                        Results will be generated inside in a dedicated directory
                        named after OUTPUT/PDBCODE_CHAIN/
```





