#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH -J chat
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o chat_out.txt
#SBATCH -e chat_err.txt


srun -N 1 -n 1 Rscript cellchat.calc.R


