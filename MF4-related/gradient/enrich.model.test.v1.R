library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(viridis)
inputdir <- "./load_files/"
outputdir <- "./report/"
source("./rgc.fun.R")


##-----------------------------------------------------------------------------------------------------
## Test enrichment models
allres <- readRDS(file = paste0(inputdir, "Region_DEGs_pairwise_rawres.rds"))


## Organize the DEX results into array (reg * reg * genes) for each cluster
## 1. (expr ratio >= 0.1)
slim_dex1 <- allres %>%
				mutate(regionpair = paste0(region, "|", bgregion)) %>%
				filter(pct.1 >= 0.1 & pct.2 <= 0.75 & ratio_fc >= 1.25 & avg_logFC >= 0.2 & p_val_adj <= 0.01)
## 2. (expr ratio < 0.1 &  0.05)
slim_dex2 <- allres %>%
				mutate(regionpair = paste0(region, "|", bgregion)) %>%
				filter(pct.1 >= 0.05 & pct.1 < 0.1 & ratio_fc >= 3 & avg_logFC >= 0.1 & p_val_adj <= 0.01)
## 3. 
slim_dex3 <- allres %>%
				mutate(regionpair = paste0(region, "|", bgregion)) %>%
				filter(pct.2 > 0.75 & ratio_fc > 1 & avg_logFC >= 1 & p_val_adj <= 0.01)
slim_dex <- rbind(slim_dex1, slim_dex2) %>%
				rbind(., slim_dex3)


all_genes <- unique(slim_dex$gene)
dex_sum <- split(slim_dex, slim_dex$cluster) %>%
			lapply(., function(res) {
				dexi <- split(res$gene, res$regionpair) %>%
							lapply(., function(y) setNames(as.numeric(all_genes %in% y), all_genes)) %>%
							as.data.frame(., check.names = FALSE) %>%
							as.matrix() %>%
							spdf2arry(df = .)
				dexi
				})


piei <- lapply(names(dex_sum), function(cls) {
	data <- dex_sum[[cls]]
	pre_df <- SummariseArray_Combined(ary = data, all_sps = c("FC", "MSC", "TC", "OcC")) %>%
				mutate(cluster = cls)
	return(pre_df)
	}) %>%
		do.call(rbind, .)
saveRDS(piei, file = paste0(inputdir, "Region_DEGs_pairwise_test_enrichpatterns.rds"))







