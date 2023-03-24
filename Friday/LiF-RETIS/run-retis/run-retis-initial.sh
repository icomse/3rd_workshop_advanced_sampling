#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH -t 6:00:00
#SBATCH --job-name=LiFRETIS
#SBATCH -e stderr-LiF
#SBATCH -o stdout-LiF
#SBATCH --partition=RM-shared
#SBATCH --account=see220002p
#SBATCH --mem=20g

export OMP_NUM_THREADS=1
# Run RETIS with GROMACS 2021
singularity exec /ocean/projects/see220002p/shared/icomse_cpu.sif python global-manager-init.py

