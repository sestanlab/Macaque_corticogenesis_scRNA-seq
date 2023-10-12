args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(dplyr)
library(tibble)
source("./net.fun.v2.R")


## system("cp ../MF_pat_v3/load_files/PAT_inte.organizer.inte.rds ./load_files/")
pat <- readRDS(file = "./load_files/PAT_inte.organizer.inte.rds")
pat@meta.data$subclass <- gsub("GE_NE_NKX2-1", "GE NKX2-1", pat@meta.data$subclass)


## Genes for permutation
##load(file = paste0("./load_files/", "PAT_markers_mres_wiGE.Rdata"))
sel_cls <- c("PC FGF17", "PC NKX2-1 NKX6-2", "PC NKX2-1 LMO1", "GE_RG_NKX2-1_DLK1", "GE_RG_NKX2-1_OLIG1", "PC RSPO3", "PC TTR", "PC SFRP2", "PC TCF7L2")
## system("cp ../MF_pat_v3/load_files/PAT_markers.Rdata ./load_files/PAT_markers.Rdata")
load(file = paste0("./load_files/PAT_markers.Rdata")
##cc_genes <- get_genes(input_genes = rownames(pat), gene_type = "cc", revised = TRUE)

genes <- lapply(sel_cls, function(cls) mar_res[[cls]]$gene) %>%
		unlist() %>% unique()


## Subset the seurat object
sel_cls <- c("PC FGF17", "PC NKX2-1", "GE NKX2-1", "PC RSPO3", "PC TTR", "PC SFRP2", "PC TCF7L2")
ctp <- sel_cls[as.numeric(args[1])]
pat <- pat[genes, ]
subseu <- subset(pat, subclass %in% ctp)


## Average expression function
perm_avgs <- CalcUMAPbin_PermAvg(object = subseu, dims = 1:30, nbins = 30, npermutations = 1000, ncores = 12)
perm_avgs <- lapply(perm_avgs, function(x) x[rownames(perm_avgs[[1]]), ])
saveRDS(perm_avgs, file = paste0("./load_files/", "Pseudobulk_permuted_along_UMAP_", ctp, ".rds"))


##codes <- paste0("Rscript net.perm.avg.R ", 1:7)
##writeLines(codes, con = "net.perm.avg.task.txt")
##dsq --job-file net.perm.avg.task.txt --batch-file net.perm.avg.sh -c 12 -p scavenge --mem-per-cpu=15G -J allg --max-jobs 7 -N 1



