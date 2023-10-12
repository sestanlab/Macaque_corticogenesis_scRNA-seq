#!/bin/bash
#SBATCH --output dsq-augur.run.upper.task-%A_%2a-%N.out
#SBATCH --array 0-17%18
#SBATCH --job-name augur
#SBATCH -c 6 -p scavenge --mem-per-cpu=15G -N 1

# DO NOT EDIT LINE BELOW
/ycga-gpfs/apps/hpc/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/sestan.ycga/sm2726/NHPfetal/ExNv4/augur.run.upper.task.txt --status-dir /gpfs/gibbs/pi/sestan.ycga/sm2726/NHPfetal/ExNv4

