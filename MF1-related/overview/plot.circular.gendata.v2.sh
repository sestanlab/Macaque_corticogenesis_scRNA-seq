#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -J GEN
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o GEN_out.txt
#SBATCH -e GEN_err.txt



Rscript plot.circular.gendata.v2.R



