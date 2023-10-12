#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH -J scvelo
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o velo_out.txt
#SBATCH -e velo_err.txt


srun -N 1 -n 1 python lin.velo.stream.hem.py &
srun -N 1 -n 1 python lin.velo.stream.antihem.py &
srun -N 1 -n 1 python lin.velo.stream.fgf17.py &
srun -N 1 -n 1 python lin.velo.stream.av.py &
wait




