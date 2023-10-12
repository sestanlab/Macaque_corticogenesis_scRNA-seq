library(Seurat)
library(dplyr)


##--------------------------------------------------------------------------
## Load the full IPC-ExN data
seu <- readRDS(file = "../overview/ExN_data_08312022.rds")



##--------------------------------------------------------------------------
## Find the DEX genes (compare one region to each of other regions)
cls_ord <- c("IPC EOMES NEUROG1", "IPC EOMES NHLH1 up", "ExN CUX2 PALMD", "ExN CUX2 ADRA2A", "ExN CUX2 ACTN2")
nseu <- seu[, (seu@meta.data$cbnage %in% c("E93", "E110")) & (seu@meta.data$subtype %in% cls_ord) & seu@meta.data$lobe != "Insula"]


## Perform pair-wise comparison between regions
nseu@meta.data$tmpcls <- paste0(nseu@meta.data$lobe, "|", nseu@meta.data$subtype)


allcls <- levels(as.factor(nseu@meta.data$tmpcls))
complist <- lapply(allcls, function(x) {
	cls <- strsplit(x, "|", fixed = TRUE) %>% 
				sapply(., function(mm) mm[2])
	bg_cls <- grep(paste0("\\|", cls, "$"), allcls, value = TRUE) %>% 
				setdiff(., x)
	paste0(x, "--vs--", bg_cls)
	}) %>% unlist()
Idents(nseu) <- "tmpcls"
DefaultAssay(nseu) <- "RNA"
allres <- lapply(complist, function(xx) {
		cls1 <- strsplit(xx, "--vs--", fixed = TRUE)[[1]][1]
		cls2 <- strsplit(xx, "--vs--", fixed = TRUE)[[1]][2]
		res <- FindMarkers(nseu, ident.1 = cls1, ident.2 = cls2, only.pos = TRUE, max.cells.per.ident = 1000, min.pct = 0.15, logfc.threshold = 0.2) %>%
				tibble::rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(background = cls2) %>%
				mutate(region = strsplit(cls1, "|", fixed = TRUE)[[1]][1]) %>%
				mutate(cluster = strsplit(cls1, "|", fixed = TRUE)[[1]][2])%>%
				mutate(bgregion = strsplit(cls2, "|", fixed = TRUE)[[1]][1])
		res
		}) %>%
			do.call(rbind, .)
saveRDS(allres, file = paste0("./load_files/", "DEGs_RES_bylobe_upper_E93.rds"))







