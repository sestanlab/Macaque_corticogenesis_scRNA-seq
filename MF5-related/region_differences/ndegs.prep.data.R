source("../scripts/nhpf.fun.R")


if (FALSE){
	## Do Umap & identity clusters
	exn <- readRDS(file = paste0("./load_files/intermediate/", "ZALL.ExN.E54-110.harmony.v1.harmony.rds")) 
	sel_cls <- c("IPC EOMES VIM", "IPC_EOMES_NEUROG1", "IPC_EOMES_NHLH1_deep", "IPC_EOMES_NHLH1_up", 
					"ExN_deep_nascent", "ExN_deep_KIF26A", "ExN_deep_NR4A2_GRID2", "ExN_deep_SYT6", "ExN_deep_OPRK1_NR4A2", "ExN_deep_OPRK1_SULF1", 
					"ExN_up_nascent","ExN_up_KCNV1","ExN_up_ADRA2A","ExN_up_ACTN2")
	exn <- subset(exn, subtype %in% sel_cls)


	set.seed(42)
	exn <- RunUMAP(exn, dims = 1:30, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")


	exn <- FindNeighbors(exn, dims = 1:30, reduction = "harmony", k.param = 25) %>%
	                    FindClusters(., resolution = 1.2, n.iter = 20)
	saveRDS(exn, file = paste0("./load_files/intermediate/", "ExN_all_harmony.rds"))
}



if (FALSE){
	## Do Umap & identity clusters
	exn <- readRDS(file = paste0("./load_files/intermediate/", "ExN_all_harmony.rds"))

	#DimFig(exn, group.by = c("subtype"), file_name = "ZALL.ExN.E54-110.harmony.v1", plot.scale = 1.5)
	#DimFig(exn, group.by = c("seurat_clusters"), file_name = "ZALL.ExN.E54-110.harmony.v1", plot.scale = 1.5)


	exn <- exn[, !exn@meta.data$seurat_clusters %in% "17"]
	set.seed(42)
	exn <- RunUMAP(exn, dims = 1:30, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
	DimFig(exn, group.by = c("subtype"), file_name = "ZALL.ExN.E54-110.harmony.v1.fil", plot.scale = 1.5)
	saveRDS(exn, file = paste0("./load_files/intermediate/", "ExN_all_harmony_filtered.rds"))	
}



exn <- readRDS(file = paste0("./load_files/intermediate/", "ExN_all_harmony_filtered.rds"))


## IT ExNs
it_cls <- c("IPC EOMES VIM", "IPC_EOMES_NEUROG1", "IPC_EOMES_NHLH1_up", "ExN_up_nascent", "ExN_up_KCNV1", "ExN_up_ADRA2A", "ExN_up_ACTN2", "ExN_deep_OPRK1_NR4A2", "ExN_deep_OPRK1_SULF1")
exn_it <- exn[, exn@meta.data$subtype %in% it_cls]
set.seed(42)
exn_it <- RunUMAP(exn_it, dims = 1:30, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
DimFig(exn_it, group.by = c("subtype"), file_name = "ZALL.ExN.E54-110.harmony.v1.fil_IT", plot.scale = 1.5)
exn_it <- FindNeighbors(exn_it, dims = 1:30, reduction = "harmony", k.param = 25) %>%
	                    FindClusters(., resolution = 1, n.iter = 20)
DimFig(exn_it, group.by = c("seurat_clusters"), file_name = "ZALL.ExN.E54-110.harmony.v1.fil_IT", plot.scale = 1.5)
saveRDS(exn_it, file = paste0("./load_files/", "ExN_all_harmony_filtered_IT.rds"))	



## non-IT ExNs
nonit_cls <- c("IPC EOMES VIM", "IPC_EOMES_NEUROG1", "IPC_EOMES_NHLH1_deep", "ExN_deep_nascent", "ExN_deep_KIF26A", "ExN_deep_NR4A2_GRID2", "ExN_deep_SYT6")
exn_nit <- exn[, exn@meta.data$subtype %in% nonit_cls]
set.seed(42)
exn_nit <- RunUMAP(exn_nit, dims = 1:30, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
DimFig(exn_nit, group.by = c("subtype"), file_name = "ZALL.ExN.E54-110.harmony.v1.fil_nonIT", plot.scale = 1.5)
exn_nit <- FindNeighbors(exn_nit, dims = 1:30, reduction = "harmony", k.param = 25) %>%
	                    FindClusters(., resolution = 1, n.iter = 20)
DimFig(exn_nit, group.by = c("seurat_clusters"), file_name = "ZALL.ExN.E54-110.harmony.v1.fil_nonIT", plot.scale = 1.5)
saveRDS(exn_nit, file = paste0("./load_files/", "ExN_all_harmony_filtered_nonIT.rds"))	


















