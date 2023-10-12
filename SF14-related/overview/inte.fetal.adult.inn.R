library(Seurat)
library(dplyr)
source("./inte.fun.R")


library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 45*1000*1024^2)



## Prepare the dataset
##fetal <- readRDS(file = paste0(inputdir, "InN_inte_E93-110.harmony.rds"))
fetal <- readRDS(file = paste0("./load_files/", "InN_data_09122022.rds"))
fetal <- subset(fetal, cbnage %in% c("E93", "E110"))
obj_list <- SplitObject(fetal, split.by = "cell_origin")
obj_list <- obj_list[sapply(obj_list, ncol) > 1500]
fetal_hvg <- lapply(obj_list, function(x) FindVariableFeatures(x, nfeatures = 3000)) %>%
			SelectIntegrationFeatures(., nfeatures = 2500)
rm(obj_list)


adult <- readRDS(file = paste0(inputdir, "Rhesus_adult_InN.rds"))
adult@meta.data$hres <- NA
adult_hvg <- SplitObject(adult, split.by = "samplename") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 3000)) %>%
			SelectIntegrationFeatures(., nfeatures = 2500)


hvg <- union(fetal_hvg, adult_hvg)


## Use the shared genes
sh_genes <- intersect(rownames(fetal), rownames(adult))
hvg <- intersect(hvg, sh_genes)
seu.list <- merge(x = adult[sh_genes, ], y = fetal[sh_genes, ]) %>%
			SplitObject(., split.by = "samplename")


file_name <- paste0("AdultFetal_InN")


## Do the integration
seu <- Integratelist.seurat(obj.list = seu.list, hvg = hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:25, cluster.dims = 1:30)

