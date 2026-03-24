#!/bin/bash

#SBATCH -p DEEP_LEARNING
#SBATCH --time=7-00:00:00         # Walltime: max depends on the partition chosen
#SBATCH --cpus-per-task=90      # Nombre de coeurs
#SBATCH --mem=150G               # 20000 Mb = 20 Gb by default, you may ask for more if needed, up to maximum 60000 Mb.
#SBATCH --gres=gpu:1            # Nombre de GPU (aucune par défaut)
#SBATCH --nodes=1               # Number of nodes
#SBATCH --output=SWORD3_SLURM_OUT/%x.%j.out      # Standard output + error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=gabriel.cretin@inserm.fr
#SBATCH -J "SWORD3"             # Job name

# Move to the directory where the job was submitted
cd $SLURM_SUBMIT_DIR

ml purge

date;pwd

source .venv/bin/activate

python3 sword2_dl/train.py --config configs/default.yaml

date
