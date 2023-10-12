#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -J pbavgsub
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o pbavgsub_out.txt
#SBATCH -e pbavgsub_err.txt



srun -N 1 -n 1 Rscript calc.pseudobulk.avg.subset.rgc.v2.R

