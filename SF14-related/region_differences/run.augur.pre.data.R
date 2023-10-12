library(Seurat)
library(dplyr)

##--------------------------------------------------------------------------
## Organize InN data for comparisons
inn <- readRDS(file = paste0("../overview/load_files/", "InN_data_09122022.rds"))


cls_list <- list(MGE = c("InN LHX6 CCK", "InN LHX6 DCN", "InN LHX6 GUCY1A2", "InN LHX6 MAF", "InN LHX6 SST NPY", "InN LHX6 SST RELN"), 
				CGE = c("InN NR2F2 CRH", "InN NR2F2 LAMP5", "InN NR2F2 SP8", "InN NR2F2 SP8 KIT", "InN NR2F2 VIP", "InN SP8 VIP", "InN SP8 CRH"))
## InN NKX2-1 CCND2 not included, largely in GE
nseu <- inn[, inn@meta.data$cbnage %in% c("E93", "E110") & 
			inn@meta.data$subtype %in% unname(unlist(cls_list)) & 
			(!inn@meta.data$region %in% c("LGE", "MGE", "CGE"))]
nseu[["group"]] <- ifelse(nseu$subtype %in% cls_list[["MGE"]], "MGE-InN", "CGE-InN")
nseu$avgcls <- paste0(nseu$cbnage, "|", nseu$region, "|", nseu$group)



## Balance region numbers
msize <- 1000
all_cls <- levels(as.factor(nseu$avgcls))
set.seed(0)
cells <- lapply(all_cls, function(cls) {
	subc <- colnames(nseu)[nseu$avgcls == cls]
	if (length(subc) > msize){
		subc <- sample(subc, msize)
	}
	subc
	}) %>%
	unlist() %>%
	unique()
seu_use <- nseu[, cells]
saveRDS(seu_use, file = paste0("./load_files/", "Augur_predata_seu.rds"))



## Select HVGs for Augur analysis
hvg_e93 <- subset(seu_use, cbnage == "E93") %>%
			FindVariableFeatures(., nfeatures = 1000) %>%
			VariableFeatures()
hvg_e110 <- subset(seu_use, cbnage == "E110") %>%
			SplitObject(., split.by = "samplename") %>%
			lapply(., function(x) FindVariableFeatures(x, nfeatures = 1500)) %>%
			SelectIntegrationFeatures(., nfeatures = 1000)
hvg <- union(hvg_e93, hvg_e110)
saveRDS(hvg, file = paste0("./load_files/", "Augur_predata_hvg.rds"))



























