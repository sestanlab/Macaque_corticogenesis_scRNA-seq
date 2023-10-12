## Generate pseudobulks of the data FOR gene co-expression
library(Seurat) 
library(ggplot2)
library(dplyr)
library(tibble)

## system("cp ../MF_pat_v3/load_files/PAT_inte.organizer.inte.rds ./load_files/")
pat <- readRDS(file = "./load_files/PAT_inte.organizer.inte.rds")


## For each subtype, generate one-dimension UMAP and calculate average expression along the UMAP axis
all_ctps <- table(pat$subclass) %>% names()
avg_list <- lapply(all_ctps, function(ctp) {
		print(ctp)
		subseu <- subset(pat, subclass %in% ctp)
		source("./net.fun.R")
		subavgs <- CalcUMAPbinAvg(object = subseu, dims = 1:30, nbins = 30) 
		return(subavgs)
		}) %>%
		setNames(., all_ctps)

avg_flist <- lapply(avg_list, function(x) x[rownames(avg_list[[1]]), ])


## Update cluster names
names(avg_flist)[names(avg_flist) %in% "GE_NE_NKX2-1"] <- "GE NKX2-1"

saveRDS(avg_flist, file = "./load_files/Pseudobulk_along_UMAP.rds")



