## R-4.1.0

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(dplyr)


loomDir <- "../../data/loom/merged/"
cbn_ldat <- readRDS(file = paste0(loomDir, "Dorsal_slim_velocyto_loom.rds"))


## Subset the loom matrices
inputdir <- "./load_files/"
#all_lins <- c("Hem", "FGF17", "NKX21", "Antihem")##[3]
#for (lin in all_lins){
#	message(paste0("Working on lineage: ", lin))
lin <- 'AV'
	data <- readRDS(file = paste0(inputdir, "Pat_lineage_", lin, ".rds"))
	slim_ldat <- lapply(cbn_ldat, function(x) x[, colnames(data)])


	## Transform to a seurat object
	seu <- as.Seurat(x = slim_ldat)
	seu[["RNA"]] <- seu[["spliced"]]
	seu[["pca"]] <- CreateDimReducObject(embeddings = data$mnn@cell.embeddings, key = "PC_", assay = "RNA")
	seu[["umap"]] <- CreateDimReducObject(embeddings = data$umap@cell.embeddings, key = "UMAP_", assay = "RNA")
	seu@meta.data$cluster <- data@meta.data$subtype
	Idents(seu) <- "cluster"


	DefaultAssay(seu) <- "RNA"
	SaveH5Seurat(seu, filename = paste0(inputdir, "velo.", lin, ".h5Seurat"), overwrite = TRUE)
	Convert(paste0(inputdir, "velo.", lin, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)
#}





