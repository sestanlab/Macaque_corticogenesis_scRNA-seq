#!/bin/bash
#SBATCH -N 4
#SBATCH -n 4
#SBATCH -c 12
#SBATCH -J heatmap
#SBATCH -p scavenge
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH -o heatmap_out.txt
#SBATCH -e heatmap_err.txt


srun -N 1 -n 1 Rscript lin.plot.heat.R Hem &
srun -N 1 -n 1 Rscript lin.plot.heat.R AV &
srun -N 1 -n 1 Rscript lin.plot.heat.R Antihem &
srun -N 1 -n 1 Rscript lin.plot.heat.R FGF17 &
wait

