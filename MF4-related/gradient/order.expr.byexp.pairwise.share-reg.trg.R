library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)
library(viridis)
inputdir <- "./load_files/"
outputdir <- "./report/"

source("./ptime.fun.v2.R")



## First compute the fold changes of expression
load(file = paste0("../region_differences/load_files/", "Smooth_by_region_tRG.Rdata"))  ##trg_smt, trg_meta, 

avgs <- as.matrix(trg_smt)
pmeta <- trg_meta %>%
			rownames_to_column("avgcls") %>%
			as.data.frame()
all_regs <- c("FC", "MSC", "TC", "OcC")


## Load region-enriched genes
degs <- readRDS(file = paste0(inputdir, "Region_DEGs_pairwise_test_summarizedpatterns.rds")) %>%
			.$feature %>% unique()


## Then order genes by regional fold changes
tRG_time <- lapply(all_regs, function(reg) {
	meta <- filter(pmeta, region == reg)
	kvec <- c(20, 10)

	time <- FitImpulseModel(genes = degs, meta = meta, pt_col = "pseudotime", knot_col = "avgcls", data = avgs, kvec = kvec, min.slope = 0.2) %>%
				mutate(region = reg, lineage = "NESC-tRG")
	return(time)
	}) %>%
	setNames(., all_regs)
saveRDS(tRG_time, file = paste0(inputdir, "Order_by_EXPR_Pairwise_Share-region_tRG_v1.rds"))









