#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J allinn
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=45G
#SBATCH -o allinn_out.txt
#SBATCH -e allinn_err.txt



Rscript inte.allinn.harmony.v1.R




