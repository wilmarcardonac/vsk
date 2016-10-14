#!/bin/sh 
#SBATCH --cpus-per-task=12
#SBATCH --job-name=inpainted-smica-8
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-user=Wilmar.Cardona@unige.ch
#SBATCH --mail-type=ALL
#SBATCH --partition=dpt
#SBATCH --output=slurm-%J.out


srun ./vsk