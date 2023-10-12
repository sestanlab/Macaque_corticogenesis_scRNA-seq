source("../scripts/nhpf.fun.R")


##-----------------------------------------------------------------------------
## Prepare the data
pat <- readRDS(file = paste0("./load_files/", "PAT_inte.E37-43.final.03062022.rds"))

meta <- readRDS(file = paste0("../../MF1/overview/load_files/Reanno_E37-110.org.meta.10052022.rds")) ##"AntVen NKX2-1 NKX6-2", 


pat@meta.data$subtype <- "removed"
sh_cells <- intersect(rownames(meta), colnames(pat))
pat@meta.data[sh_cells, "subtype"] <- meta[sh_cells, "subtype"]

cls_ord <- c("AntVen NKX2-1 LMO1", "inIPC ASCL1 DLX1", "InN LHX8 ZIC1", "InN GNRH1", "InN HMX1") ## "PC NKX2-1 NKX6-2" doesn't seem to give rise to InNs and were not included.
pat <- subset(pat, subtype %in% cls_ord & region %in% "FR")



##-----------------------------------------------------------------------------
## Prepare the data
lin <- "AV"
pat <- RunUMAP(pat, dims = 1:30, umap.method = "umap-learn", metric = "correlation")


## rodo clustering & spot cycling cells
pat <- FindNeighbors(pat, dims = 1:30) %>%
                FindClusters(., resolution = 1)
DimFig(pat, file_name = paste0("Lineage_", lin), plot.scale = 0.5, pt.size = 0.4, group.by = c("subtype", "seurat_clusters"))
FeatureFig(pat, file_name = paste0("Lineage_", lin), plot.scale = 0.5, pt.size = 0.2, features = c("MKI67", "PCNA"))


## Remove cycling clusters
rm_cls <- c("8", "6", "7")
fseu <- pat[, !pat@meta.data$seurat_clusters %in% rm_cls]
fseu <- RunUMAP(fseu, dims = 1:30, umap.method = "umap-learn", metric = "correlation", min.dist = 0.4, n.neighbors = 40, spread = 0.6)

DimFig(fseu, file_name = paste0("Lineage_final_", lin), plot.scale = 0.5, pt.size = 0.4, group.by = c("subtype"))

saveRDS(fseu, file = paste0(inputdir, "Pat_lineage_", lin, ".rds"))







