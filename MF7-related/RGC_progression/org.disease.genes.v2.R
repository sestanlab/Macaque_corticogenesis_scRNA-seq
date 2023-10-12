library(dplyr)
library(tibble)


## This rds file comes from Xoel
alllist <- readRDS("./load_files/all_diseases_list.rds")


## Build a gene table
all_genes <- unlist(alllist) %>% unique()
alltb <- matrix(0, nrow = length(all_genes), ncol = length(alllist), dimnames = list(all_genes, names(alllist)))
for (ii in names(alllist)){
	alltb[, ii] <- ifelse(rownames(alltb) %in% alllist[[ii]], 1, 0)
}


save(alltb, alllist, file = paste0("./load_files/", "Disease_genes_v3.Rdata"))


