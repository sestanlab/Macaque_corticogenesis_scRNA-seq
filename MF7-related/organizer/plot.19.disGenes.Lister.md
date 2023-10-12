---
title: expression of the 19 PAT-specific disease genes in Lister data
author: Shaojie Ma
date: Dec 6, 2022
---


## Organize PFC data
- Copy data
```bash
cp /home/sm2726/myshare2/PublicData/Lister_2022_Cell/Lister_processed_RNA_seurat_CorrectStage.rds ./load_files/
## Cluster order
```

- Source functions for visualization
vis.fun.R


- Average expression across periods
```R
library(Seurat)
library(dplyr)

seu <- readRDS(file = "./load_files/Lister_processed_RNA_seurat_CorrectStage.rds")
seu$newcls <- paste0(seu$stage, "|", seu$major_clust)


DefaultAssay(seu) <- "DS"
##allcls <- levels(as.factor(seu$newcls))
allctp <- levels(as.factor(seu$major_clust))
res_list <- lapply(allctp, function(ctp) {
	subseu <- subset(seu, major_clust %in% ctp)
	ctm <- subseu$DS@data
	identity <- setNames(subseu$newcls, colnames(subseu))
	subcls <- levels(as.factor(subseu$newcls))
	avg <- lapply(subcls, function(cls) {
				print(cls)
				subident <- identity[identity %in% cls]
				subctm <- ctm[, names(subident), drop = FALSE]
				subavg <- Matrix::rowMeans(expm1(subctm))
				return(subavg)
				}) %>%
				setNames(., subcls) %>%
				as.data.frame(., check.names = FALSE) %>%
				as.matrix()
	ratio <- lapply(subcls, function(cls) {
				print(cls)
				subident <- identity[identity %in% cls]
				subctm <- ctm[, names(subident), drop = FALSE]
				subavg <- Matrix::rowMeans(subctm != 0)
				return(subavg)
				}) %>%
				setNames(., subcls) %>%
				as.data.frame(., check.names = FALSE) %>%
				as.matrix()
	return(list(avg = avg, ratio = ratio))
	}) %>%
	setNames(., allctp)
avgs <- lapply(res_list, function(x) x$avg) %>%
			do.call(cbind, .)
ratios <- lapply(res_list, function(x) x$ratio) %>%
			do.call(cbind, .)
cls_size <- table(seu$newcls)
save(avgs, ratios, cls_size, file = "./load_files/Lister_Expr_avg_ratio_by_period_downsample.Rdata")

```


- Visualization
```R
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
##library(Seurat)
library(cowplot)

load(file = "./load_files/Lister_Expr_avg_ratio_by_period_downsample.Rdata")
## avgs, ratios, cls_size


stage_ord <- c("Fetal", "Neonatal", "Infancy", "Childhood", "Adolescence", "Adult")
cls_ord <- c("L2-3_CUX2", "L4_RORB", "L5-6_THEMIS", "L5-6_TLE4", "PN_dev", "VIP", "ID2", "LAMP5_NOS1", "CGE_dev", "SST", "PV", "PV_SCUBE3", "MGE_dev", "Astro", "OPC", "Oligo", "Micro", "Vas")


genes <- c("LYN", "FGF8", "SLC35D2", "EFHD1", "SPATA33", "SLC7A6", "SLC7A5", "SLC29A4", "SAMD11", "ORAI3", "MYB", "MET", "MECOM", "KCNJ13", "GRM8", "EMB", "CNN2", "CGNL1", "CADPS2", "BCAM", "ADA", "WNT4", "WNT3", "SHOX2", "RNF43", "CA14") %>%
			intersect(., rownames(avgs))
## Gene SHOX2 was not found in the dataset


type_nodes <- grep("Micro|Vas|Poor-Quality", colnames(ratios), invert = TRUE, value = TRUE) %>%
				intersect(., names(cls_size)[cls_size >= 50])
type_max <- apply(ratios[genes, type_nodes], 1, max)
## type1 genes (Not expressed in neural lineage)
type1_genes <- genes[type_max < 0.05]

## type2 genes (lowly expressed in neural lineage)
type2_genes <- genes[type_max >= 0.05 & type_max < 0.15]


source("./vis.fun.R")
pp <- CirclePlot.horizontal.Lister(avg = avgs, ratio = ratios, features = genes, file_name = "test", dot.min = 0.05, dot.scale = 5, scale.by = "radius", shape = 16, cluster.order = cls_ord, split.order = stage_ord, stroke.col = "black", stroke.size = 0.5, col.min = -2.5, col.max = 2.5, scale = TRUE, scale.min = NA, scale.max = NA, mask.vector = setNames(as.numeric(cls_size), names(cls_size)), return.plot = TRUE, width.scale = 1, height.base = 1.5, font.scale = c(1, 1), height.unit = 0.3, min.cells = 50)


## yy colors
yy_colors <- setNames(rep("#000000", length(genes)), rev(genes))
yy_colors[type1_genes] <- "#00AEEF"
yy_colors[type2_genes] <- "#bf812d"

pp <- pp + 
		theme(axis.text.y = element_text(color = yy_colors))



pdf(paste0("./report/Lister_19_PAT_disease_markers_colored.pdf"), width = 15, height = 4.5, useDingbats = FALSE)
print(pp)
dev.off()


```


