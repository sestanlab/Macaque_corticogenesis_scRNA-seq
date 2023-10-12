args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(dplyr)

source("./inte.fun.R")


## parallelization of CCA
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 45*1000*1024^2)


lb <- args[1] ## region


## Get ExN from all ages (just dlPFC)
seu <- readRDS(file = paste0("../overview/load_files/", "ExN_data_raw_05242022.rds"))
fetal <- subset(seu, lobe %in% lb & cbnage %in% c("E110", "E54", "E62-64", "E77-78", "E93") & samplename != "RMB683")
fetal$samplename[fetal$cbnage == "E77-78"] <- gsub("_.*_", "_", fetal$cell_origin[fetal$cbnage == "E77-78"])


## Get highly variable genes
obj_list <- SplitObject(fetal, split.by = "samplename")
obj_list <- obj_list[sapply(obj_list, ncol) >= 2000]
obj_list <- obj_list %>%
                lapply(., function(x) FindVariableFeatures(x, nfeatures = 1500))

fetal_hvg <- SelectIntegrationFeatures(obj_list[grep("E54|E62|E64", names(obj_list))], nfeatures = 1200) %>%
        c(., SelectIntegrationFeatures(obj_list[grep("E77|78", names(obj_list))], nfeatures = 1200)) %>%
        c(., SelectIntegrationFeatures(obj_list[c("E93", "RMB691")], nfeatures = 1200)) %>%
        unique()
fetal@meta.data$inte.batch <- gsub("E62|E64", "E62-64", fetal@meta.data$samplename) %>%
            gsub("E77_0501|E77_1222|E78_0501", "E77-78", .)


## Adult samples
rmb196 <- readRDS(file = paste0("./load_files/", "RMB196_ExN.rds"))
rmb196@meta.data$inte.batch <- "AdultFR"
adult_hvg <- FindVariableFeatures(rmb196, nfeatures = 2000) %>%
                VariableFeatures()

hvg <- union(fetal_hvg, adult_hvg)


## Use the shared genes
sh_genes <- intersect(rownames(fetal), rownames(rmb196))
hvg <- intersect(hvg, sh_genes)
seu.list <- merge(x = rmb196[sh_genes, ], y = fetal[sh_genes, ]) %>%
            SplitObject(., split.by = "inte.batch")


file_name <- paste0("AdultFetal_ExN_", lb, "_v1")


## Do the integration
seu <- Integratelist.seurat(obj.list = seu.list, hvg = hvg, file_name = file_name, input_dir = inputdir, inte.dims = 1:30, cluster.dims = 1:30, reference = NULL) ##which(names(seu.list) %in% c("E93", "E62-64", "AdultFR")))



##codes <- paste0("Rscript inte.fetal.adult.seurat.v1.R ", c("FC", "MSC", "TC", "OC"))
##writeLines(codes, con = "inte.fetal.adult.seurat.task.txt")


##dSQ generate the sh file
##module load dSQ
##dsq --job-file inte.fetal.adult.seurat.task.txt --batch-file inte.fetal.adult.seurat.sh -c 4 -p scavenge --mem-per-cpu=45G -J seurat --max-jobs 4 -N 1






