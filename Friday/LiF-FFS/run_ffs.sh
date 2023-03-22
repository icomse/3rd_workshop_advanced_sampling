#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --partition=RM-shared
#SBATCH --time=8:00:00
#SBATCH --account=see220002p
#SBATCH -e stderr-lif
#SBATCH -o stdout-lif
#SBATCH --mem=62g 

chmod +x scripts/setup.sh
chmod +x global-manager.py
chmod +x scripts/run.sh
chmod +x scripts/op.py

#use the singularity as below, or use python from an environment with required packages
singularity exec /ocean/projects/see220002p/shared/icomse_cpu.sif python global-manager.py    # run FFS
