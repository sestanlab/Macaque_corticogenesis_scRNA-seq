library(Seurat)
library(dplyr)
source('./inte.fun.R')


obj <- readRDS(file = paste0("./load_files/", "InN_data_09122022.rds"))



## Find highly variable genes
obj$samplename[obj$cbnage == "E77-78"] <- gsub("_.*_", "_", obj$cell_origin[obj$cbnage == "E77-78"])
early_hvg <- subset(obj, cbnage %in% c("E37", "E42-43")) %>%
			SplitObject(., split.by = "samplename") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 3000)) %>%
			SelectIntegrationFeatures(., nfeatures = 2500)
late_hvg <- subset(obj, cbnage %in% c("E54", "E62-64", "E77-78", "E93", "E110")) %>%
			SplitObject(., split.by = "samplename") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 3000)) %>%
			SelectIntegrationFeatures(., nfeatures = 2500)
hvg <- union(early_hvg, late_hvg)



file_name <- "InN_inte_AllAge"
seu <- InteAllSp.harmony(object = obj, split.by = "samplename", hvg = hvg, file_name = file_name, input_dir = inputdir, inte.dims = 30, theta = 2, lambda = 0.8, sigma = 0.1)



