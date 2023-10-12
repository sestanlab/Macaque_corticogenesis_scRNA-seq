## source("../scripts/nhpf.fun.R")
library(Seurat)
library(dplyr)
source("./inte.fun.R")


library(future)
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 50*1000*1024^2)

inputdir <- "./load_files/"


## Pollen data
other <- readRDS(file = paste0("./Pollen_data/Pollen_data_full_seuratv3.rds"))
other <- subset(other, age %in% c("E80", "E90", "E100"))
other$cluster <- other$class
other$subtype <- NA
other$samplename <- other$age
other_hvg <- SplitObject(other, split.by = "samplename") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 3000)) %>%
			SelectIntegrationFeatures(., nfeatures = 2500)


## Prepare the fetal dataset
fetal <- readRDS(file = paste0("./load_files/", "InN_data_09122022.rds"))
fetal <- subset(fetal, cbnage %in% c("E93", "E110"))
obj_list <- SplitObject(fetal, split.by = "cell_origin")
obj_list <- obj_list[sapply(obj_list, ncol) > 1500]
fetal_hvg <- lapply(obj_list, function(x) FindVariableFeatures(x, nfeatures = 3000)) %>%
			SelectIntegrationFeatures(., nfeatures = 2500)
rm(obj_list)


hvg <- union(fetal_hvg, other_hvg)


## Use the shared genes
sh_genes <- intersect(rownames(fetal), rownames(other))
hvg <- intersect(hvg, sh_genes)
seu.list <- merge(x = other[sh_genes, ], y = fetal[sh_genes, ]) %>%
			SplitObject(., split.by = "samplename")


file_name <- paste0("IntegrateWithPollen_InN_seurat")


## Do the integration
seu <- Integratelist.seurat(obj.list = seu.list, hvg = hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:25, cluster.dims = 1:30)



