library(Seurat)
library(Matrix)
library(dplyr)

inputdir <- "./load_files/"
outputdir <- "./report/"


## The following dataset has a "DS" assay containing subset data acorss ages
rgc <- readRDS(file = paste0("../overview/load_files/", "RGC_filtered_seu_04112022_addsubset.rds"))
rgc@meta.data$cluster2 <- gsub(" DIRAS3| TEX15", "", rgc@meta.data$cluster)



## Remove certain cells
## 1. cells with smaller cluster size
## 2. tRG E93-110 cells
## 3. frontal-specific subtype
## 4. regions not belong to FC, MSC, TC, OC (we used "region" here, not "lobe")
## 5. some low-quality ependymal cells
rm_cls <- c("Insula|Ependymal", 
		"MSC|Ependymal", "MSC|NEP RSPO3", "MSC|tRG CRYAB FNDC1",
		"TC|Ependymal", "TC|NEP RSPO3", "TC|vRG HMGA2 CCND1", "TC|tRG CRYAB FNDC1")
rm_list <- list(colnames(rgc)[paste0(rgc@meta.data$lobe, "|", rgc@meta.data$cluster2) %in% rm_cls], 
		rownames(rgc@meta.data)[rgc@meta.data$cluster2 == "tRG CRYAB MEST" & rgc@meta.data$cbnage %in% c("E93", "E110")], 
		colnames(rgc)[rgc@meta.data$cluster2 == "RGC FABP7 PMP22"], 
		colnames(rgc)[rgc@meta.data$region %in% c("IPC", "A1C", "PCC", "Insula")],
		colnames(rgc)[rgc@meta.data$cluster2 == "Ependymal" & rgc@meta.data$nFeature_RNA < 1250]
		)
rm_cells <- Reduce("union", rm_list)

rgc <- rgc[, setdiff(colnames(rgc), rm_cells)]
saveRDS(rgc, file = paste0(inputdir, "RGC_seu_for_DEG_analysis.rds"))



## Have a balanced number of cells across regions
tb <- as.matrix(table(rgc$cluster2, rgc$lobe))
set.seed(0)
fcells <- lapply(rownames(tb), function(cls) {
	size <- tb[cls, ]
	msize <- min(size[size != 0])
	cls_cells <- lapply(colnames(tb), function(reg) {
		subc <- colnames(rgc)[rgc$lobe == reg & rgc$cluster2 == cls]
		if (length(subc) > msize){
			subc <- sample(subc, msize)
		}
		return(subc)
		}) %>%
		unlist() %>% unique()
	return(cls_cells)
	}) %>%
	unlist() %>%
	unique()
subrgc <- rgc[, fcells]
table(subrgc$cluster2, subrgc$lobe)
saveRDS(fcells, file = paste0(inputdir, "Region_DEGs_cells.rds"))




## Do differential gene expression
rgc <- readRDS(file = paste0(inputdir, "RGC_seu_for_DEG_analysis.rds"))
fcells <- readRDS(file = paste0(inputdir, "Region_DEGs_cells.rds"))
subrgc <- rgc[, fcells]

subrgc@meta.data$tmpcls <- paste0(subrgc@meta.data$lobe, "|", subrgc@meta.data$cluster2)


complist1 <- levels(as.factor(subrgc@meta.data$tmpcls)) 
complist2 <- lapply(complist1, function(x) {
	cls <- strsplit(x, "|", fixed) %>%
				sapply(., "[", 2) %>% 
				setNames(., NULL)
	bg_cls <- grep(paste0("\\|", cls, "$"), complist1, value = TRUE) %>% 
				setdiff(., x)
	return(bg_cls)
	})
Idents(subrgc) <- "tmpcls"
DefaultAssay(subrgc) <- "RNA"
allres <- lapply(1:length(complist1), function(idx) {
		cls1 <- complist1[idx]
		cls2 <- complist2[[idx]]
		res <- FindMarkers(subrgc, ident.1 = cls1, ident.2 = cls2, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1, max.cells.per.ident = sum(subrgc$tmpcls == cls1)) %>%
				tibble::rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(region = strsplit(cls1, "|", fixed = TRUE)[[1]][1]) %>%
				mutate(cluster = strsplit(cls1, "|", fixed = TRUE)[[1]][2])
		res
		}) %>%
			do.call(rbind, .)
saveRDS(allres, file = paste0(inputdir, "Region_DEGs_rawres.rds"))


