library(Seurat)
library(dplyr)


inn <- readRDS(file = paste0("./load_files/", "InN_data_09122022.rds"))


cls_list <- list(MGE = c("InN LHX6 CCK", "InN LHX6 DCN", "InN LHX6 GUCY1A2", "InN LHX6 MAF", "InN LHX6 SST NPY", "InN LHX6 SST RELN"), 
				CGE = c("InN NR2F2 CRH", "InN NR2F2 LAMP5", "InN NR2F2 SP8", "InN NR2F2 SP8 KIT", "InN NR2F2 VIP", "InN SP8 VIP", "InN SP8 CRH"))
## InN NKX2-1 CCND2 not included, largely in GE
nseu <- inn[, inn@meta.data$cbnage %in% c("E93", "E110") & 
			inn@meta.data$subtype %in% unname(unlist(cls_list)) & 
			(!inn@meta.data$region %in% c("LGE", "MGE", "CGE"))]
nseu[["group"]] <- ifelse(nseu$subtype %in% cls_list[["MGE"]], "MGE-InN", "CGE-InN")

##aggregate(nFeature_RNA ~ region + cbnage, nseu@meta.data, mean)
##aggregate(nFeature_RNA ~ region + cbnage, nseu@meta.data, median)



##--------------------------------------------------------------------------
## Find the DEX genes (compare one region to each of other regions)
gp <- c("MGE-InN", "CGE-InN")[2]


seu_use <- subset(nseu, group == gp)
## Control the number of cells in each region
##msize <- table(nseu@meta.data$region, nseu@meta.data$group) %>%
##			min()
msize <- 2000
all_regs <- levels(as.factor(seu_use$region))
set.seed(0)
cells <- lapply(all_regs, function(reg) {
	subc <- colnames(seu_use)[seu_use$region == reg]
	if (sum(seu_use$region == reg) > msize){
		subc <- sample(subc, msize)
	}
	subc
	}) %>%
	unlist() %>%
	unique()
seu_use <- seu_use[, cells]
seu_use@meta.data$avg_cls <- seu_use@meta.data$region


## Calculate DEGs & Avgs
Idents(seu_use) <- "region"
mar_res <- FindAllMarkers(seu_use, max.cells.per.ident = 2000, logfc.threshold = 0.2, min.pct = 0.15, only.pos = TRUE) %>%
			mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01))
mars <- mar_res %>%
			filter(ratio_fc >= 1.2 & pct.1 >= 0.15 & p_val_adj <= 0.01) %>%
			.$gene %>% unique()

Idents(seu_use) <- "avg_cls"
avgs <- log(AverageExpression(seu_use)$RNA + 1)
scale_avg <- avgs[mars, ,drop = FALSE] %>%
            as.matrix() %>%
            t() %>% scale() %>% t() %>%
            MinMax(., min = -1.5, max = 2)

res <- list(mar_raw = mar_res, 
			mars = mars,
			avg = avgs, 
			savg = scale_avg)
saveRDS(res, file = paste0("./load_files/", "DEG_res_", gp, "_v1.rds"))



## Random generate 10 pseudobulks for each region (200 cells per pdbulk)
allavg <- NULL
set.seed(42)
for (reg in all_regs){
	subc <- colnames(seu_use)[seu_use$region == reg]
	subc <- replicate(20, sample(subc, 200), simplify = FALSE)

	sub_avg <- lapply(1:length(subc), function(x) {
		subseu <- seu_use[, subc[[x]]]
		Idents(subseu) <- paste0(reg, "|", x)
		avg <- as.matrix(log(AverageExpression(subseu)$RNA + 1))
		return(avg)
		}) %>%
		do.call(cbind, .)
	allavg <- cbind(allavg, sub_avg)
}
saveRDS(allavg, file = paste0("./load_files/", "Avg_res_", gp, "_v1.rds"))







