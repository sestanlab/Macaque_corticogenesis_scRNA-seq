## Calculate region-average expression (all the genes) for (RGearly, oRG, ExN-up, ExN-deep)
source("../scripts/nhpf.fun.R")


##------------------------------------------------------------------------
## Load all data
fmeta <- readRDS(file = paste0("../../MF1/overview/load_files/intermediate/Reanno_E37-110.org.meta.rds"))
fmeta$subtype <- fmeta$subtype2


seu <- readRDS(file = paste0("../../MF1/overview/load_files/", "All.MNN.v1.mnn.rds"))
seu$RNA@scale.data <- matrix(0, nrow = 1, ncol = 1)
seu <- seu[, rownames(fmeta)]

## Update meta.data
sel_cls <- c("subtype", "subclass", "cbnage", "lobe", "region", "cell_origin", "samplename", "age")
for (ii in sel_cls){
	seu@meta.data[, ii] <- fmeta[colnames(seu), ii]
}


##------------------------------------------------------------------------
## Specify cells for each lineage
merge_list <- list(NESC = c("NESC RSPO3 DIRAS3", "NESC RSPO3 TEX15"),
				vRG_early = "vRG HMGA2 CCND1",
				vRG_late = "vRG SAT1 STMN2",
				oRG = c("oRG HOPX APOE"),
				`IPC EOMES NEUROG1` = c("IPC EOMES NEUROG1"),
				`IPC EOMES NHLH1` = c("IPC EOMES NHLH1 deep", "IPC EOMES NHLH1 up"), 
				`ExN L6B` = "ExN SOX5 NR4A2 GRID2",
				`ExN L6CT` = "ExN SOX5 SYT6",
				`ExN upper` = "ExN CUX2 ACTN2",
				`MGE-InN` = c("InN LHX6 CCK", "InN LHX6 DCN", "InN LHX6 GUCY1A2", "InN LHX6 MAF", "InN LHX6 SST NPY", "InN LHX6 SST RELN"), 
				`CGE-InN` = c("InN NR2F2 CRH", "InN NR2F2 LAMP5", "InN NR2F2 SP8", "InN NR2F2 SP8 KIT", "InN NR2F2 VIP", "InN SP8 VIP", "InN SP8 CRH"), 
				gIPC = "gIPC EGFR",
				`gIPC cycling` = "gIPC EGFR MKI67",
				aIPC = "aIPC IGFBP2",
				oIPC = "oIPC DLL1",
				Astro = c("Astro EGFR", "Astro GFAP", "Astro MFGE8"),
				OPC = "OPC PDGFRA",
				Microglia = "Microglia PTPRC C1QC",
				Endothelial = c("Endothelial CLDN5 FLT1", "Endothelial CLDN5 RGS5")
				)
merge_age <- list(NESC = c("E37", "E42-43"),
				vRG_early = c("E37", "E42-43"),
				vRG_late = c("E54", "E62-64", "E77-78"),
				oRG = c("E93", "E110"),
				`IPC EOMES NEUROG1` = c("E54", "E62-64", "E77-78"),
				`IPC EOMES NHLH1` = c("E54", "E62-64", "E77-78"), 
				`ExN L6B` = c("E62-64", "E77-78"),
				`ExN L6CT` = c("E62-64", "E77-78"),
				`ExN upper` = c("E93", "E110"),
				`MGE-InN` = c("E93", "E110"), 
				`CGE-InN` = c("E93", "E110"), 
				gIPC = c("E93", "E110"),
				`gIPC cycling` = c("E93", "E110"),
				aIPC = c("E93", "E110"),
				oIPC = c("E93", "E110"),
				Astro = c("E93", "E110"),
				OPC = c("E93", "E110"),
				Microglia = c("E93", "E110"),
				Endothelial = c("E93", "E110")
				)


seu@meta.data$subtype2 <- seu@meta.data$subtype
for (ii in names(merge_list)){
	seu@meta.data$subtype2[seu@meta.data$subtype %in% merge_list[[ii]]] <- ii
}



## Balance cell number across regions
cells <- lapply(names(merge_list), function(gp) {
		## Specify ages for the given cell groups
		ag_use <- merge_age[[gp]]
		subc <- colnames(seu)[seu@meta.data$subtype2 %in% gp & seu@meta.data$cbnage %in% ag_use & seu@meta.data$lobe %in% c("FC", "MSC", "TC", "OC")]


		## Subset regions for certain cell types
		reg_size <- table(seu@meta.data[subc, "lobe"])
		if (gp %in% c("NESC")){
			reg_size <- reg_size[c("FC", "OC")]
		}
		if (gp %in% c("vRG_early")){
			reg_size <- reg_size[c("FC", "MSC", "OC")]
		}
		if (gp %in% c("oIPC")){
			reg_size <- reg_size[c("FC", "MSC", "TC")]
		}

		if (min(reg_size) < 100){
			maxcells <- 100
		} else {
			maxcells <- min(reg_size)
		}


		## Balance cell number across regions
		set.seed(0)
		subc_subset <- lapply(names(reg_size), function(reg) {
			regc <- subc[seu@meta.data[subc, "lobe"] == reg]
			if (length(regc) > maxcells){
				cc <- sample(regc, maxcells)
			} else {
				cc <- regc
			}
			return(cc)
			}) %>%
			unlist()
		print(gp)
		print(table(seu@meta.data[subc_subset, "lobe"]))
		return(subc_subset)
		}) %>%
		unlist()
seu_use <- seu[, cells]
saveRDS(seu_use, file = paste0(inputdir, "Reg-DEG_data_seu.rds"))




## Identify region-specific genes
res_list <- lapply(names(merge_list), function(gp) {
	print(gp)
	subseu <- seu_use[, seu_use@meta.data$subtype2 %in% gp]
	reg_size <- table(subseu$lobe)
	maxcells <- max(reg_size)
	# Generate expressed genes as well (for odds background)
	genes <- lapply(names(reg_size), function(reg) {
		regseu <- subseu[, subseu@meta.data$lobe %in% reg]
		gg <- rownames(regseu$RNA)[Matrix::rowMeans(regseu$RNA@counts != 0) >= 0.1]
		gg
		}) %>%
		unlist() %>% unique()


	Idents(subseu) <- "lobe"
	res <- FindAllMarkers(subseu, features = genes, max.cells.per.ident = maxcells, only.pos = TRUE, logfc.threshold = 0.1, min.pct = 0.1) %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(region = cluster) %>%
				mutate(cluster = gp)
	return(list(res = res, genes = genes))
	}) %>%
	setNames(., names(merge_list))


exp_list <- lapply(res_list, function(x) x$genes)
saveRDS(exp_list, file = paste0(inputdir, "Expgenes_res.rds"))


allres <- lapply(res_list, function(x) x$res) %>% do.call(rbind, .)
saveRDS(allres, file = paste0(inputdir, "DEG_res.rds"))



##----------------------------------------------------------------------
## Get expression ratios
CalcMeanRatio <- function(object, group.by, assay = "RNA") {
	allcls <- levels(as.factor(object@meta.data[, group.by]))
	mat <- lapply(allcls, function(cls){
		print(paste0("Working on cluster: ", cls))
		subobj <- object[[assay]]@data[, object@meta.data[, group.by] == cls, drop = FALSE]
		expr <- Matrix::rowMeans(subobj != 0)
		return(expr)
		}) %>%
		do.call(cbind, .) %>%
		as.matrix()

	rownames(mat) <- rownames(object[[assay]])
	colnames(mat) <- allcls
	return(mat)
}

seu_use@meta.data$avgcls <- paste0(seu_use$subtype2, "|", seu_use$lobe)
exprs <- CalcMeanRatio(object = seu_use, group.by = "avgcls", assay = "RNA")
saveRDS(exprs, file = paste0(inputdir, "ExprRatio_res.rds"))






















