source("../scripts/nhpf.fun.R")
source("./dis.fun.R")


## Generate data subset
if (!file.exists(paste0("./load_files/", "RGC_data_subset_09122022.rds"))){
	rgc <- readRDS(file = "./load_files/RGC_data_09122022.rds")
	rgc <- rgc[, rgc@meta.data$nCount_RNA >= 2000]

	new_ctm <- SubsampleData.count.equalUMI(counts = rgc$RNA@counts, numi = 2000, seed.use = 42, nCores = 12)
	rgc[["DS"]] <- CreateAssayObject(counts = new_ctm[, colnames(rgc)], min.cells = 0, min.features = 0)
	DefaultAssay(rgc) <- "DS"
	rgc <- NormalizeData(rgc, normalization.method = "LogNormalize")
	saveRDS(rgc, file = paste0("./load_files/", "RGC_data_subset_09122022.rds"))
}




## Generate pseudobulk average for a give cluster
rgc <- readRDS(file = paste0("./load_files/", "RGC_data_subset_09122022.rds"))
rgc_pbavgs <-  ParallelAVG_Allcls(object = rgc, group.by = "subtype", nreps = 100, ncells = 50, nCores = 6, seed.use = 42)
saveRDS(rgc_pbavgs, file = paste0("./load_files/", "RGC_pseudobulk_avgs_subset.rds"))



## Load Disease gene sets
rgc_pbavgs <- readRDS(file = paste0("./load_files/", "RGC_pseudobulk_avgs_subset.rds"))
load("./load_files/Disease_genes_v3.Rdata") ## alltb, alllist
alllist <- lapply(alllist, function(x) intersect(x, rownames(rgc_pbavgs)))


auc <- GetModuleScore_NEW(assay.data = rgc_pbavgs, features = alllist, seed = 42, method = "aucell", input_dir = "./load_files/", file_name = paste0("AUCell_RGC_pseudobulk_subset"), output_dir = "./report/", rethreshold_list = NULL, cellbin.size = 500, nCores = 4)













