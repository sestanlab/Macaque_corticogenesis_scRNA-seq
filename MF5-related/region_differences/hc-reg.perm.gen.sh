#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -J perm
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o perm_out.txt
#SBATCH -e perm_err.txt



Rscript hc-reg.perm.gen.R


