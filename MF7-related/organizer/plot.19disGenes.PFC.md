---
title: expression of the 19 PAT-specific disease genes in PFC data
author: Shaojie Ma
date: Dec 6, 2022
---


## Organize PFC data
- Copy data
```bash
## Average expression
cp ~/project/PFC/Shiny_prep/load_files/Expr_avg_ratio_by_species.rds ./load_files/PFC_Expr_avg_ratio_by_species.rds


## Cluster order
cp ~/project/PFC/Shiny_prep/load_files/Cluster_order.rds ./load_files/PFC_Cluster_order.rds
```

- Source functions for visualization
vis.fun.R



- Visualization
```R
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)

res <- readRDS(file = "./load_files/PFC_Expr_avg_ratio_by_species.rds")[['mres']]
cls_ord <- readRDS(file = "./load_files/PFC_Cluster_order.rds")[["shared"]][["mres"]]

genes <- c("LYN", "FGF8", "SLC35D2", "EFHD1", "SPATA33", "SLC7A6", "SLC7A5", "SLC29A4", "SAMD11", "ORAI3", "MYB", "MET", "MECOM", "KCNJ13", "GRM8", "EMB", "CNN2", "CGNL1", "CADPS2", "BCAM", "ADA", "WNT4", "WNT3", "SHOX2", "RNF43", "CA14") %>%
			intersect(., rownames(res$avg))



type_nodes <- grep("Micro|Immune|Endo|RB$|SMC|VLMC", colnames(res$ratio), invert = TRUE, value = TRUE) %>%
				grep("\\|PC$", ., invert = TRUE, value = TRUE) %>%
				grep("Human|Rhesus", ., value = TRUE)
type_max <- apply(res$ratio[genes, type_nodes], 1, max)
## type1 genes (Not expressed in neural lineage)
type1_genes <- genes[type_max < 0.1]

## type2 genes (lowly expressed in neural lineage)
type2_genes <- genes[type_max >= 0.1 & type_max < 0.2]




source("./vis.fun.R")
pp <- CirclePlot.horizontal(avg = res$avg, ratio = res$ratio, features = genes, dot.min = 0.1, file_name = "test", dot.scale = 4, scale.by = "radius", cluster.order = cls_ord, stroke.size = 0, stroke.matrix = NULL, width.scale = 1, height.base = 2, row.cex.sf = 1, col.cex.sf = 1, return.plot = TRUE)

## yy colors
yy_colors <- setNames(rep("#000000", length(genes)), rev(genes))
yy_colors[type1_genes] <- "#00AEEF"
yy_colors[type2_genes] <- "#bf812d"

pp <- pp + 
		theme(axis.text.y = element_text(color = yy_colors))


pdf(paste0("./report/PFC_19_PAT_disease_markers_colored.pdf"), width = 8, height = 6, useDingbats = FALSE)
print(pp)
dev.off()


```

