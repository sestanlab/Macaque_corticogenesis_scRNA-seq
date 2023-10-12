library(Seurat)
library(dplyr)
library(Matrix)
source("./dis.fun.R")


## Generate pseudobulk average for a give cluster
rgc <- readRDS(file = "./load_files/RGC_data_09122022.rds")
rgc_pbavgs <-  ParallelAVG_Allcls(object = rgc, group.by = "subtype", nreps = 100, ncells = 50, nCores = 6, seed.use = 42)
saveRDS(rgc_pbavgs, file = paste0("./load_files/", "RGC_pseudobulk_avgs_v2.rds"))



## Load Disease gene sets
rgc_pbavgs <- readRDS(file = paste0("./load_files/", "RGC_pseudobulk_avgs_v2.rds"))
load("./load_files/Disease_genes_v3.Rdata") ## alltb, alllist
alllist <- lapply(alllist, function(x) intersect(x, rownames(rgc_pbavgs)))


auc <- GetModuleScore_NEW(assay.data = rgc_pbavgs, features = alllist, seed = 42, method = "aucell", input_dir = "./load_files/", file_name = paste0("AUCell_RGC_pseudobulk"), output_dir = "./report/", rethreshold_list = NULL, cellbin.size = 500, nCores = 4)






