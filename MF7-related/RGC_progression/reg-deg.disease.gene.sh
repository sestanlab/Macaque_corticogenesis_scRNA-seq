#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -J deg
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o deg_out.txt
#SBATCH -e deg_err.txt



srun -N 1 -n 1 Rscript reg-deg.disease.gene.R

