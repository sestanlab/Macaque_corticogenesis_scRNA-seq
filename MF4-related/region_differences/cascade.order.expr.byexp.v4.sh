#!/bin/bash
#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 4
#SBATCH -J urdexp4
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=45G
#SBATCH -o urdexp4_out.txt
#SBATCH -e urdexp4_err.txt



srun -N 1 -n 1 cascade.Rscript order.expr.byexp.org.v4.R &
srun -N 1 -n 1 cascade.Rscript order.expr.byexp.trg.v4.R &
wait


