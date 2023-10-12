library(Seurat)
library(dplyr)
inputdir <- "./load_files/"
outputdir <- "./report/"


## Cell type average expression
exn <- readRDS(file = "../overview/load_files/ExN_data_08312022.rds")
lexn <- subset(exn, cbnage %in% c("E93", "E110"))
rm(exn)


## RMB683 has very limited nnumber of cells (6397 cells)
hvg <- subset(lexn, inte.batch.spread %in% c("E93", "RMB691")) %>%
				SplitObject(., split.by = "inte.batch.spread") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500)) %>% 
				SelectIntegrationFeatures(., nfeatures = 2000)

saveRDS(hvg, file = paste0("./load_files/Hier-region-cluster.hvg.rds"))















