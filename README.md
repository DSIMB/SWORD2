# SWORD2

SWift and Optimized Recognition of protein Domains

The SWORD2 partitioning algorithm produces multiple alternative domain assignments for a given protein structure. This unique approach handles ambiguous protein structure partitioning, admitting several solutions. The decomposition of the protein structure into domains is achieved through the hierarchical clustering of Protein Units, evolutionarily preserved structural descriptors at the interface between secondary structures and domains.

This is the repository of the standalone version of SWORD2 proposed as a webserver: https://www.dsimb.inserm.fr/SWORD2/index.html

## Abstract

Understanding the functions and origins of proteins requires splitting these macromolecules into fragments that could be independent in terms of folding, activity, or evolution. For that purpose, structural domains are the typical level of analysis, but shorter segments, such as subdomains and supersecondary structures, are insightful as well. Here, we propose SWORD2, a web server for exploring how an input protein structure may be decomposed into ‘Protein Units’ that can be hierarchically assembled to delimit structural domains. For each partitioning solution, the relevance of the identified substructures is estimated through different measures. This multilevel analysis is achieved by integrating our previous work on domain delineation, ‘protein peeling’ and model quality assessment. We hope that SWORD2 will be useful to biologists searching for key regions in their proteins of interest and to bioinformaticians building datasets of protein structures.

## Install

Install the conda environment using the `environment.yml` file:
```
conda env create -f environment.yml
or
mamba env create -f environment.yml
```

## How to use SWORD2

First, activate the working environment:
```
conda activate sword2
```

Then, launch SWORD2:

### On a PDB id:
```
./SWORD2.py -p 1jx4 -o results
```

### On an AlphaFold predicted structure using its Uniprot Accession Id:
```
./SWORD2.py -u Q76EI6 -o results
```

### On your own PDB structure:
```
./SWORD2.py -f ./structure.pdb -o results
```

## Help

To get the full help:
```
$ ./SWORD2.py --help
usage: SWORD2.py [-h] (-u UNIPROT_CODE | -p PDB_CODE | -f PDB_FILE) [-c PDB_CHAIN] -o OUTPUT

options:
  -h, --help            show this help message and exit
  -u UNIPROT_CODE, --uniprot-code UNIPROT_CODE
                        AlphaFold Uniprot Accession Id
  -p PDB_CODE, --pdb-code PDB_CODE
                        PDB code
  -f PDB_FILE, --pdb-file PDB_FILE
                        PDB file path
  -c PDB_CHAIN, --pdb-chain PDB_CHAIN
                        PDB chain. Default is A.
  -o OUTPUT, --output OUTPUT
                        Output directory
```





