library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(viridis)
inputdir <- "./load_files/"
outputdir <- "./report/"



## load data: generated in the script: deg.region.calc.v1.R
rgc <- readRDS(file = paste0("../region_differences/load_files/", "RGC_seu_for_DEG_analysis.rds"))
fcells <- readRDS(file = paste0("../region_differences/load_files/", "Region_DEGs_cells.rds"))
subrgc <- rgc[, fcells]

subrgc@meta.data$tmpcls <- paste0(subrgc@meta.data$lobe, "|", subrgc@meta.data$cluster2)


## Do differential gene expression
allcls <- levels(as.factor(subrgc@meta.data$tmpcls))
complist <- lapply(allcls, function(x) {
	cls <- sapply(strsplit(x, "|", fixed = TRUE), "[", 2)
	bg_cls <- grep(paste0("\\|", cls, "$"), allcls, value = TRUE) %>% 
				setdiff(., x)
	paste0(x, "--vs--", bg_cls)
	}) %>% unlist()

Idents(subrgc) <- "tmpcls"
DefaultAssay(subrgc) <- "RNA"
tb <- as.matrix(table(subrgc$cluster2, subrgc$lobe))
allres <- lapply(complist, function(xx) {
		cls1 <- strsplit(xx, "--vs--", fixed = TRUE)[[1]][1]
		cls2 <- strsplit(xx, "--vs--", fixed = TRUE)[[1]][2]
		cls <- sapply(strsplit(cls1, "|", fixed = TRUE), "[", 2)
		size <- tb[cls, ]
		msize <- min(size[size != 0])
		res <- FindMarkers(subrgc, ident.1 = cls1, ident.2 = cls2, only.pos = TRUE, max.cells.per.ident = msize, min.pct = 0.05, logfc.threshold = 0.1) %>%
				rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(background = cls2) %>%
				mutate(region = strsplit(cls1, "|", fixed = TRUE)[[1]][1]) %>%
				mutate(cluster = strsplit(cls1, "|", fixed = TRUE)[[1]][2])%>%
				mutate(bgregion = strsplit(cls2, "|", fixed = TRUE)[[1]][1])
		res
		}) %>%
			do.call(rbind, .)
saveRDS(allres, file = paste0(inputdir, "Region_DEGs_pairwise_rawres.rds"))





