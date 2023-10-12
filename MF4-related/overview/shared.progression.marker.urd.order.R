args <- commandArgs(trailingOnly = TRUE)
source("~/project/PFC/scripts/pfc.fun.R")
source("./ptime.fun.R")
library(parallel)
library(foreach)


avgres <- readRDS(file = paste0(inputdir, "Pseudobulk_cbn.rds"))
load(file = paste0(inputdir, "Ptime_mars_cbn_step1.Rdata"))
## res_list, mars_early, mars_org, mars_trg


lin <- args[1]

meta_use <- avgres[["raw"]]$meta[avgres[["raw"]]$meta$lineage == lin, ]
data_use <- as.matrix(avgres[["raw"]]$avg[, rownames(meta_use)])


genes <- switch(lin,
				shared = mars_early, 
				oRG = mars_org, 
				tRG = mars_trg)
genes <- unlist(genes) %>% unique()
pcres <- PulsefitProcess(data_use = data_use, genes = genes, ptime = setNames(meta_use$pseudotime, rownames(meta_use)), k = 50, interpolate = 50, pulse.only = TRUE, nCores = 4)
saveRDS(pcres, file = paste0(inputdir, "Ptime_mars_cbn_URDorder_", lin, ".lev1.rds"))









