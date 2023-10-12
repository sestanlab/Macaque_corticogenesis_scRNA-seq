#!/bin/bash
#SBATCH --output dsq-net.perm.avg.task-%A_%1a-%N.out
#SBATCH --array 0-6%7
#SBATCH --job-name allg
#SBATCH -c 12 -p scavenge --mem-per-cpu=15G -N 1

# DO NOT EDIT LINE BELOW
/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/sestan.ycga/sm2726/NHPfetal/PAT_network/net.perm.avg.task.txt --status-dir /gpfs/gibbs/pi/sestan.ycga/sm2726/NHPfetal/PAT_network

