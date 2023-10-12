library(Seurat)

pfc_exn <- readRDS(file = "~/project/PFC/data/Final_PFC_HPRC.ExN.rm1680.07102021.rds")
rhe_exn <- subset(pfc_exn, species == "Rhesus")

## Find the variable genes
obj.list <- SplitObject(rhe_exn, split.by = "samplename") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000))
saveRDS(obj.list, file = paste0("./load_files/", "Rhesus_adult_ExN.rds"))


## Only use one sample to speed up the integration
rmb196 <- obj.list[["RMB196"]]
saveRDS(rmb196, file = paste0("./load_files/", "RMB196_ExN.rds"))


