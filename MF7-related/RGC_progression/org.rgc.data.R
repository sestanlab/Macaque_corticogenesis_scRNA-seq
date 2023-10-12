library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
inputdir <- "./load_files/"
outputdir <- "./report/"



##------------------------------------------------------------------------
## Load all data
fmeta <- readRDS(file = paste0("../../MF1/overview/load_files/intermediate/Reanno_E37-110.org.meta.rds"))
fmeta$subtype <- fmeta$subtype2


seu <- readRDS(file = paste0("../../MF1/overview/load_files/", "All.MNN.v1.mnn.rds"))
seu$RNA@scale.data <- matrix(0, nrow = 1, ncol = 1)
seu <- seu[, rownames(fmeta)]

## Update meta.data
sel_cls <- c("subtype", "subclass", "cbnage", "lobe", "region", "cell_origin", "samplename", "age")
for (ii in sel_cls){
	seu@meta.data[, ii] <- fmeta[colnames(seu), ii]
}


## Subset to NSC lineage
rgc <- subset(seu, subclass %in% c("dorsal NSC", "GE NSC"))


## Further prune the meta.data
rm_cols <- setdiff(colnames(rgc@meta.data), c(c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "seurat_clusters"), sel_cls))
for (ii in rm_cols){
	rgc@meta.data[, ii] <- NULL
}
saveRDS(rgc, file = "./load_files/RGC_data_09122022.rds")








