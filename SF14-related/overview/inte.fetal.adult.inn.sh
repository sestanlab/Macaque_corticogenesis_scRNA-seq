#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH -J seurat
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH -o seurat_out.txt
#SBATCH -e seurat_err.txt



srun -N 1 -n 1 Rscript inte.fetal.adult.inn.R




