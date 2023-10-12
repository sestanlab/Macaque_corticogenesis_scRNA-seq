library(Seurat)
library(dplyr)


## Meta data
fmeta <- readRDS(file = paste0("./load_files/intermediate/Reanno_E37-110.org.meta.rds"))
fmeta$subtype <- fmeta$subtype2


seu <- readRDS(file = paste0("./load_files/", "All.MNN.v1.mnn.rds"))
seu$RNA@scale.data <- matrix(0, nrow = 1, ncol = 1)
seu <- seu[, rownames(fmeta)]
seu@meta.data$subtype <- fmeta[colnames(seu), "subtype"]
seu@meta.data$subclass <- fmeta[colnames(seu), "subclass"]
seu@meta.data$cbnage <- fmeta[colnames(seu), "cbnage"]
seu@meta.data$lobe <- fmeta[colnames(seu), "lobe"]



## Subset cells to calculate average expression
set.seed(0)
allcls <- levels(as.factor(fmeta$subtype))
cells <- lapply(allcls, function(x) {
	subc <- rownames(fmeta)[fmeta$subtype == x]
	if (length(subc) >= 1000){
		subc <- sample(subc, 1000)
	}
	return(subc)
	}) %>%
		unlist()



subseu <- seu[, cells]
subseu@meta.data$subtype <- fmeta[colnames(subseu), "subtype"]

ratios <- lapply(allcls, function(x) {
		message(paste0("Calculate expression ratio for cluster: ", x))
		data <- subseu[, as.character(subseu@meta.data$subtype) == x]
		er <- Matrix::rowMeans(data$RNA@counts != 0)
		er
		}) %>% setNames(., allcls) %>%
			as.data.frame(., check.names = FALSE) %>% 
			as.matrix()
colnames(ratios) <- gsub("_", " ", colnames(ratios))



Idents(subseu) <- "subtype"
avgs <- log(AverageExpression(subseu)$RNA + 1)
colnames(avgs) <- gsub("_", " ", colnames(avgs))
###save(fmeta, avgs, ratios, file = paste0("./load_files/", "MF1_expr_avg_ratio.Rdata"))


UpdateName <- function(x) {
	x[x %in% "PC SFRP1"] <- "Ant SFRP1"
	x[x %in% "PC NKX2-1 LMO1"] <- "AntVen NKX2-1 LMO1"
	x[x %in% "PC NKX2-1 NKX6-2"] <- "AntVen NKX2-1 NKX6-2"
	x[x %in% "PC NKX2-1 RAX"] <- "Pos NKX2-1 RAX"
	return(x)
}

colnames(avgs) <- UpdateName(colnames(avgs))
colnames(ratios) <- UpdateName(colnames(ratios))
fmeta <- readRDS(file = paste0("./load_files/Reanno_E37-110.org.meta.10052022.rds"))
save(fmeta, avgs, ratios, file = paste0("./load_files/", "MF1_expr_avg_ratio_v3.Rdata"))









