#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH -J ptime
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o ptime_out.txt
#SBATCH -e ptime_err.txt


srun -N 1 -n 1 python lin.velo.ptime.hem.py &
srun -N 1 -n 1 python lin.velo.ptime.antihem.py &
srun -N 1 -n 1 python lin.velo.ptime.fgf17.py &
srun -N 1 -n 1 python lin.velo.ptime.av.py &
wait




