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
seu@meta.data$subtype <- fmeta[colnames(seu), "subtype"]
seu@meta.data$subclass <- fmeta[colnames(seu), "subclass"]
seu@meta.data$cbnage <- fmeta[colnames(seu), "cbnage"]
seu@meta.data$lobe <- fmeta[colnames(seu), "lobe"]
seu@meta.data$region <- fmeta[colnames(seu), "region"]
seu@meta.data$cell_origin <- fmeta[colnames(seu), "cell_origin"]
seu@meta.data$samplename <- fmeta[colnames(seu), "samplename"]
seu@meta.data$age <- fmeta[colnames(seu), "age"]



## Subset to ExN lineage
inn <- subset(seu, subclass %in% c("inIPC", "Inhibitory neurons"))


## Further prune the meta.data
rm_cols <- setdiff(colnames(inn@meta.data), c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "seurat_clusters", "subtype", "subclass", "cbnage", "lobe", "region", "cell_origin", "samplename", "age"))
for (ii in rm_cols){
	inn@meta.data[, ii] <- NULL
}


## Determine the integration batch
inn@meta.data$inte.batch.spread <- inn@meta.data$samplename
inn@meta.data$inte.batch.spread[inn@meta.data$inte.batch.spread %in% c("E37", "E42", "E43")] <- "E37-43"
inn@meta.data$inte.batch.spread[inn@meta.data$cbnage %in% "E77-78"] <- gsub("_.*_", "_", inn@meta.data$cell_origin)[inn@meta.data$cbnage %in% "E77-78"]


## Further remove some doublets for cleaner analysis
e93_ident <- readRDS(file = paste0("./load_files/intermediate/", "QC_E93_step1_ident.rds"))
e110_ident <- readRDS(file = paste0("./load_files/intermediate/", "QC_E110_step1_ident.rds"))
all_ident <- c(e93_ident, e110_ident)
db_cells <- names(all_ident)[all_ident %in% "doublets"]

inn <- inn[, setdiff(colnames(inn), db_cells)]
saveRDS(inn, file = "./load_files/InN_data_09122022.rds")








