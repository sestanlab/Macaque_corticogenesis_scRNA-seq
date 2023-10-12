library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
inputdir <- "./load_files/"
outputdir <- "./report/"


exn <- readRDS(file = "../overview/load_files/ExN_data_08312022.rds")


##  Highly variable genes
## Deep layer neurons related 
hvg_e56 <- subset(exn, cbnage %in% c("E54", "E62-64")) %>%
                SplitObject(., split.by = "inte.batch.spread") %>%
                lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500)) %>% 
                SelectIntegrationFeatures(., nfeatures = 2000)
## Upper layer neurons related
hvg_e7 <- subset(exn, cbnage %in% c("E77-78")) %>%
                SplitObject(., split.by = "inte.batch.spread") %>%
                lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000)) %>% 
                SelectIntegrationFeatures(., nfeatures = 1500)
## RMB683 has very limited nnumber of cells (6397 cells)
hvg_e910 <- subset(exn, inte.batch.spread %in% c("E93", "RMB691")) %>%
                SplitObject(., split.by = "inte.batch.spread") %>%
                lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500)) %>% 
                SelectIntegrationFeatures(., nfeatures = 2000)
hvg <- union(hvg_e56, hvg_e7) %>% union(., hvg_e910)
saveRDS(hvg, file = paste0("./load_files/", "ExN_HVG_spread_08312022.rds"))


