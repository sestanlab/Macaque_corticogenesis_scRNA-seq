library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
source("./heatmap.fun.R")

inputdir <- "./load_files/"

##--------------------------------------------------------------------------
## Load the full IPC-ExN data
seu <- readRDS(file = paste0("./load_files/ExN.harmony.spread_region.v1.harmony.rds"))



##--------------------------------------------------------------------------
## Find the DEX genes (compare one region to each of other regions)
cls_ord <- c("IPC EOMES NEUROG1", "IPC EOMES NHLH1 up", "ExN up nascent", "ExN up ADRA2A", "ExN up ACTN2")
nseu <- seu[, (seu@meta.data$cbnage %in% c("E93", "E110")) & (seu@meta.data$subtype %in% cls_ord) & seu@meta.data$lobe != "Insula"]


## Perform pair-wise comparison between regions
nseu@meta.data$tmpcls <- paste0(nseu@meta.data$lobe, "|", nseu@meta.data$subtype)


allcls <- levels(as.factor(nseu@meta.data$tmpcls)) 
complist <- lapply(allcls, function(x) {
	cls <- extract_field(x, 2, "|")
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
				rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(background = cls2) %>%
				mutate(region = strsplit(cls1, "|", fixed = TRUE)[[1]][1]) %>%
				mutate(cluster = strsplit(cls1, "|", fixed = TRUE)[[1]][2])%>%
				mutate(bgregion = strsplit(cls2, "|", fixed = TRUE)[[1]][1])
		res
		}) %>%
			do.call(rbind, .)
save(allres, file = paste0(inputdir, "Module_pairwise_genes_up.Rdata"))




##---------------------------------------------------------------------------------
load(file = paste0(inputdir, "Module_pairwise_genes_up.Rdata"))
## allres



## Get expression ratios
allcls <- levels(as.factor(nseu@meta.data$tmpcls))
exp_ratio <- lapply(allcls, function(x) {
			message(paste0("Calculate expression ratio for cluster: ", x))	
			sub_seu <- nseu[, nseu@meta.data$tmpcls == x]
			er <- Matrix::rowMeans(sub_seu$RNA@data != 0)
			er
			}) %>% setNames(., allcls) %>%
				as.data.frame(., check.names = FALSE) %>% as.matrix()



## Organize the DEX results into array (reg * reg * genes) for each cluster
all_genes <- unique(allres$gene)
slim_dex <- allres %>%
				mutate(regionpair = paste0(region, "|", bgregion)) %>%
				subset(pct.1 >= 0.15 & ratio_fc >= 1.2 & avg_logFC >= 0.2) %>%
				group_by(cluster, regionpair) %>% 
				mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
				filter(p_val_adj <= 0.01) %>%
				ungroup()


res_cluster <- split(slim_dex, slim_dex$cluster)
source("./heatmap.fun.R")
dex_sum <- lapply(res_cluster, function(res) {
				dexi <- split(res$gene, res$regionpair) %>%
							lapply(., function(y) setNames(as.numeric(all_genes %in% y), all_genes)) %>%
							as.data.frame(., check.names = FALSE) %>%
							as.matrix() %>%
							spdf2arry(df = .)
				dexi
				})


piei <- lapply(names(dex_sum), function(cls) {
	data <- dex_sum[[cls]]
	ept_mat <- SummariseArray(ary = data) 
	pre_df <- ept_mat %>%
				as.data.frame(., check.names = FALSE) %>%
				rownames_to_column("gene")

	full_regs <- c("FC", "MSC", "TC", "OC")
	mis_regs <- setdiff(full_regs, colnames(pre_df))
	for (ii in mis_regs){
		pre_df[, ii] <- 0
	}
	isshare <- rowSums(as.matrix(pre_df[, full_regs])) == 0
	pre_df <- pre_df[, c("gene", full_regs)]
	pre_df$Shared <- ifelse(isshare, 1, 0)


	## Get the radius information
	reg_cls <- paste0(dimnames(data)[[1]], "|", cls)
	cls_ratio <- exp_ratio[all_genes, reg_cls] %>% as.matrix()
	colnames(cls_ratio) <- dimnames(data)[[1]]

	pre_df$radius[isshare] <- apply(cls_ratio[which(isshare), ], 1, mean)
	pre_df$radius[!isshare] <- sapply(which(!isshare), function(x) {
		y <- cls_ratio[x, ] * ept_mat[x, ]
		y <- y[y!= 0]
		yy <- mean(y)
		return(yy)
		})

	pre_df$cluster <- cls
	return(pre_df)
	}) %>%
		do.call(rbind, .)%>%
		.[, c("cluster", c("FC", "MSC", "TC", "OC"), "Shared", "gene", "radius")]

save(piei, file = paste0(inputdir, "Module_pairwise_pie.dot_up.Rdata"))





















