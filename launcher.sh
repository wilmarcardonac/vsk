#!/bin/sh
#SBATCH --cpus-per-task=12
#SBATCH --job-name=vsk
#SBATCH --ntasks=1
#SBATCH --time=0-20:00:00
#SBATCH --mail-user=wilmar.cardona@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --partition=dpt
#SBATCH --output=slurm-%J.out

srun ./vsk