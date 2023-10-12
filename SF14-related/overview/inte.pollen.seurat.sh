#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -J pollen
#SBATCH -p bigmem
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH -o pollen_out.txt
#SBATCH -e pollen_err.txt



srun -N 1 -n 1 Rscript inte.pollen.seurat.R




