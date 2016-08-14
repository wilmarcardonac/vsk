#!/bin/sh 
#SBATCH --cpus-per-task=12
#SBATCH --job-name=vsk
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=Wilmar.Cardona@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --partition=dpt
#SBATCH --output=slurm-%J.out


srun ./vsk