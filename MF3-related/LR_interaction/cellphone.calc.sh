#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -J cpdb
#SBATCH -p bigmem
#SBATCH -t 100:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH -o cpdb_out.txt
#SBATCH -e cpdb_err.txt


##Run the jobs 
niters=5000
ncores=24
outDir=/home/sm2726/project/NHPfetal/LR_PATRGC/load_files/CPDBraw
curDir=/home/sm2726/project/NHPfetal/LR_PATRGC/load_files/
source activate cellpheno


srun -N 1 -n 1 cellphonedb method statistical_analysis $curDir/PAT_RGC_meta.txt $curDir/PAT_RGC_count.txt --counts-data gene_name --threshold 0.04 --iterations=$niters --threads=$ncores --project-name=PATRGCv1 --output-path=$outDir --pvalue 0.05 > cpdb.out 2> cpdb.err

##plotDir=/home/sm2726/project/NHPfetal/MF_pat/load_files/raw/PATRGCv1
##cellphonedb plot dot_plot --means-path $plotDir/means.txt --pvalues-path $plotDir/pvalues.txt --output-path $plotDir --output-name PATRGC.dotplot.pdf &
##cellphonedb plot heatmap_plot $curDir/PAT_RGC_meta.txt --pvalues-path $plotDir/pvalues.txt --output-path $plotDir --count-name PATRGC.heatmap_count.pdf --log-name PATRGC.heatmap_log_count.pdf --count-network-name PATRGC.count_network.txt --interaction-count-name PATRGC.interaction_count.txt





