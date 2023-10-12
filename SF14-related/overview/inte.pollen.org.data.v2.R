## R-4.1.0
library(Seurat)
library(SeuratData)
library(SeuratDisk)


Convert("./Pollen_data/GSE169122_MacaqueDevInhibitoryNeurons.h5ad", dest = "h5seurat", overwrite = TRUE)
pollen <- LoadH5Seurat("./Pollen_data/GSE169122_MacaqueDevInhibitoryNeurons.h5seurat", assays = "RNA")
saveRDS(pollen, file = paste0("./Pollen_data/Pollen_data_full_seuratv4.rds"))


library(Matrix)
pollen <- readRDS(file = paste0("./Pollen_data/Pollen_data_full_seuratv4.rds"))
meta <- pollen@meta.data
pca <- pollen$pca@cell.embeddings
umap <- pollen$umap@cell.embeddings
ctx <- pollen$RNA@counts
save(ctx, meta, pca, umap, file = paste0("./Pollen_data/Pollen_data_full_raw.Rdata"))



## R-3.6.1
## Create a new object (version-3 seuart)
source("../scripts/nhpf.fun.R")
load(file = paste0("./Pollen_data/Pollen_data_full_raw.Rdata"))
## ctx, meta, pca, umap, 


seu <- seu_prepare(counts = ctx, min.cells = 0, normalization.method = "LogNormalize", nfeatures = 2500, hvg.method = NULL, assay = "RNA")
seu[["pca"]] <- CreateDimReducObject(embeddings = pca[colnames(seu), ], key = "PC_", assay = DefaultAssay(seu))
seu[["umap"]] <- CreateDimReducObject(embeddings = umap[colnames(seu), ], key = "UMAP_", assay = DefaultAssay(seu))

seu@meta.data <- cbind(seu@meta.data, meta[colnames(seu), ])

seu@meta.data$age <- paste0("E", seu@meta.data$timepoint)
saveRDS(seu, file = paste0("./Pollen_data/Pollen_data_full_seuratv3.rds"))










