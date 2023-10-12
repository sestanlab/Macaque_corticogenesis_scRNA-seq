source("../scripts/nhpf.fun.R")
seu <- readRDS(file = paste0(inputdir, "PAT_inte.E37-43.final.03062022.rds"))


lin_list <- list(Hem = c("PC RSPO3", "IPC RSPO3 NEUROG1", "IPC RSPO3 NHLH1", "CR TP73"),
				FGF17 = c("PC FGF17", "IPC FGF17", "Neu TAGLN3 ONECUT2"),
				Antihem = c("PC SFRP2", "IPC SFRP2 ASCL1", "IPC SFRP2 MEIS2")
				) ## "PC NKX2-1 NKX6-2" doesn't seem to give rise to InNs and were not included.
lin_regs <- list(Hem = "OX",
				FGF17 = "FR",
				Antihem = "OX")


lin <- "Hem"
subseu <- subset(seu, subtype %in% lin_list[[lin]] & region %in% lin_regs[[lin]])


if (lin == "Antihem"){
	subseu@meta.data$subtype <- gsub("IPC SFRP2 MEIS2", "InN MEIS2", subseu@meta.data$subtype)
}


## do UMAP
subseu <- RunUMAP(subseu, dims = 1:30, umap.method = "umap-learn", metric = "correlation")

## rodo clustering & spot cycling cells
subseu <- FindNeighbors(subseu, dims = 1:30) %>%
                FindClusters(., resolution = 1)
DimFig(subseu, file_name = paste0("Lineage_", lin), plot.scale = 0.5, pt.size = 0.4, group.by = c("subtype", "seurat_clusters"))
FeatureFig(subseu, file_name = paste0("Lineage_", lin), plot.scale = 0.5, pt.size = 0.2, features = c("MKI67", "PCNA"))


## Remove cycling clusters
rm_cls <- switch(lin,
			FGF17 = c("0", "4", "6", "8"), 
			Antihem = c("2", "3", "6", "7"), 
			Hem = c("0", "5", "7", "11"))
fseu <- subseu[, !subseu@meta.data$seurat_clusters %in% rm_cls]

if (lin == "Hem"){
	fseu <- RunUMAP(fseu, dims = 1:30, umap.method = "umap-learn", metric = "correlation", min.dist = 0.5, n.neighbors = 35)
} else {
	fseu <- RunUMAP(fseu, dims = 1:30, umap.method = "umap-learn", metric = "correlation")
}


DimFig(fseu, file_name = paste0("Lineage_final_", lin), plot.scale = 0.5, pt.size = 0.4, group.by = c("subtype"))
saveRDS(fseu, file = paste0(inputdir, "Pat_lineage_", lin, ".rds"))







