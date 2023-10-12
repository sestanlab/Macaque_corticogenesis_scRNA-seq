## Calculate region-average expression (all the genes) for (RGearly, oRG, ExN-up, ExN-deep)
source("../scripts/nhpf.fun.R")


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





###--------------------------------------------------------------------------------
## IPC-ExNs

seu <- readRDS(file = paste0("../ExNv3/load_files/ExN.harmony.spread_region.v1.harmony.rds"))
ipcs <- colnames(seu)[seu@meta.data$subtype %in% c("IPC EOMES NEUROG1", "IPC EOMES NHLH1 deep", "IPC EOMES NHLH1 up") & seu@meta.data$cbnage %in% c("E54", "E62-64", "E77-78") & seu@meta.data$lobe != "Insula"]
ups <- colnames(seu)[seu@meta.data$subtype %in% c("ExN up ACTN2") & seu@meta.data$cbnage %in% c("E93", "E110") & seu@meta.data$lobe != "Insula"]
deeps <- colnames(seu)[seu@meta.data$subtype %in% c("ExN deep NR4A2 GRID2", "ExN deep SYT6") & seu@meta.data$cbnage %in% c("E62-64", "E77-78") & seu@meta.data$lobe != "Insula"]

exn <- seu[, c(ipcs, ups, deeps)]


## Remove Insula, PCC, IPC, AIC
exn <- exn[, !exn@meta.data$region %in% c("A1C", "IPC", "PCC")]


cls_use <- c("IPC EOMES NEUROG1", "IPC EOMES NHLH1 deep", "IPC EOMES NHLH1 up", "ExN deep NR4A2 GRID2", "ExN deep SYT6", "ExN up ACTN2")
res_list <- lapply(cls_use, function(cls) {
	## Identify cells
	cells <- colnames(exn)[exn$subtype == cls]
	reg_size <- table(exn@meta.data[cells, "lobe"])
	if (min(reg_size) < 100){
		maxcells <- 100
	} else {
		maxcells <- min(reg_size)
	}

	set.seed(0)
	subcells <- lapply(names(reg_size), function(reg) {
		regc <- cells[exn@meta.data[cells, "lobe"] == reg]
		if (length(regc) > maxcells){
			cc <- sample(regc, maxcells)
		} else {
			cc <- regc
		}
		return(cc)
		}) %>%
		unlist()

	subseu <- exn[, subcells]

	# Generate expressed genes as well (for odds background)
	genes <- lapply(names(reg_size), function(reg) {
		regseu <- subset(subseu, lobe == reg)
		gg <- rownames(regseu$RNA)[Matrix::rowMeans(regseu$RNA@counts != 0) >= 0.1]
		gg
		}) %>%
		unlist() %>% unique()


	Idents(subseu) <- "lobe"
	res <- FindAllMarkers(subseu, max.cells.per.ident = maxcells, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1) %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(region = cluster) %>%
				mutate(cluster = cls)
	return(list(res = res, genes = genes))
	}) %>%
	setNames(., cls_use)

exp_list <- lapply(res_list, function(x) x$genes)
saveRDS(exp_list, file = paste0(inputdir, "Expgenes_res_IPC-ExN.rds"))


allres <- lapply(res_list, function(x) x$res) %>% do.call(rbind, .)
saveRDS(allres, file = paste0(inputdir, "DEG_res_IPC-ExN.rds"))


## Get expression ratios
exn[["avgcls"]] <- paste0(exn$subtype, "|", exn$lobe)
exprs <- CalcMeanRatio(object = exn, group.by = "avgcls", assay = "RNA")
saveRDS(exprs, file = paste0(inputdir, "ExprRatio_res_IPC-ExN.rds"))





###--------------------------------------------------------------------------------
## RGCs
rgc <- readRDS(file = paste0("../RGCanalysis/load_files/RGC_filtered_seu_04112022.rds"))
sel_cls <- c("NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "oRG HOPX TNC", "oRG HOPX APOE", "tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal")

rgc@meta.data$subtype <- case_when(
	rgc@meta.data$cluster %in% c("NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15") ~ "NESC",
	rgc@meta.data$cluster %in% c("vRG HMGA2 CCND1") ~ "vRG early",
	rgc@meta.data$cluster %in% c("vRG SAT1 STMN2") ~ "vRG late",
	rgc@meta.data$cluster %in% c("oRG HOPX TNC", "oRG HOPX APOE") ~ "oRG", 
	)
rgc@meta.data$subtype[is.na(rgc@meta.data$subtype)] <- "others"

## Remove Insula, PCC, IPC, AIC
rgc <- rgc[, !rgc@meta.data$region %in% c("A1C", "IPC", "PCC", "Insula")]



cls_use <- c("NESC", "vRG early", "vRG late", "oRG")
res_list <- lapply(cls_use, function(cls) {
	## Identify cells
	cells <- colnames(rgc)[rgc$subtype == cls]
	reg_size <- table(rgc@meta.data[cells, "lobe"])
	if (cls == "vRG early"){
		reg_size <- reg_size[c("FC", "MSC", "OcC")]
	} else if (cls == "NESC") {
		reg_size <- reg_size[c("FC", "OcC")]
	}

	if (min(reg_size) < 100){
		maxcells <- 100
	} else {
		maxcells <- min(reg_size)
	}

	set.seed(0)
	subcells <- lapply(names(reg_size), function(reg) {
		regc <- cells[rgc@meta.data[cells, "lobe"] == reg]
		if (length(regc) > maxcells){
			cc <- sample(regc, maxcells)
		} else {
			cc <- regc
		}
		return(cc)
		}) %>%
		unlist()

	subseu <- rgc[, subcells]
	

	# Generate expressed genes as well (for odds background)
	genes <- lapply(names(reg_size), function(reg) {
		regseu <- subset(subseu, lobe == reg)
		gg <- rownames(regseu$RNA)[Matrix::rowMeans(regseu$RNA@counts != 0) >= 0.1]
		gg
		}) %>%
		unlist() %>% unique()


	Idents(subseu) <- "lobe"
	res <- FindAllMarkers(subseu, max.cells.per.ident = maxcells, only.pos = TRUE, logfc.threshold = 0.2, min.pct = 0.1) %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(region = cluster) %>%
				mutate(cluster = cls)
	return(list(res = res, genes = genes))
	}) %>%
	setNames(., cls_use)


exp_list <- lapply(res_list, function(x) x$genes)
saveRDS(exp_list, file = paste0(inputdir, "Expgenes_res_RGC.rds"))


allres <- lapply(res_list, function(x) x$res) %>% do.call(rbind, .)
allres$region <- gsub("OcC", "OC", allres$region)
saveRDS(allres, file = paste0(inputdir, "DEG_res_RGC.rds"))


## Get expression ratios
rgc[["avgcls"]] <- paste0(rgc$subtype, "|", rgc$lobe)
exprs <- CalcMeanRatio(object = rgc, group.by = "avgcls", assay = "RNA")
exprs <- exprs[, setdiff(colnames(exprs), c("NESC|MSC", "NESC|TC", "vRG early|TC", "others|FC", "others|MSC", "others|OcC", "others|TC"))]
saveRDS(exprs, file = paste0(inputdir, "ExprRatio_res_RGC.rds"))




###--------------------------------------------------------------------------------
## Astrocytes, OPCs and InNs (no GE NSCs)
mon <- readRDS(file = paste0("../inte_all/load_files/", "All_filter_seu_v05222022.rds"))

mon$subtype2 <- ifelse(mon$subtype %in% c("Astro EGFR", "Astro GFAP", "Astro MFGE8"), "Astro", mon$subtype)
mon$subtype2[mon$subtype %in% c("InN LHX6", "InN LHX6 SST", "InN LHX6 PVALB", "InN LHX6 NKX2-1")] <- "InN MGE"
mon$subtype2[mon$subtype %in% c("InN LAMP5", "InN NR2F2 ADARB2", "InN NR2F2 RELN", "InN NR2F2", "InN NR2F1 RND3", "InN PROX1 VIP", "InN PROX1 CRH", "InN PROX1 CALB2", "InN SP8 PAX6")] <- "InN CLGE"


## Remove Insula, PCC, IPC, AIC
mon <- mon[, !mon@meta.data$region %in% c("A1C", "IPC", "PCC", "Insula")]


sel_cls <- c("gIPC", "oIPC", "aIPC", "Astro", "OPC", "InN MGE", "InN CLGE")
int <- mon[, mon$subtype2 %in% sel_cls & mon$lobe %in% c("FC", "MSC", "TC", "OC") & mon$cbnage %in% c("E93", "E110")]


cls_use <- sel_cls
res_list <- lapply(cls_use, function(cls) {
	print(cls)
	## Identify cells
	cells <- colnames(int)[int$subtype2 == cls]
	reg_size <- table(int@meta.data[cells, "lobe"])

	if (min(reg_size) < 100){
		maxcells <- 100
	} else {
		maxcells <- min(reg_size)
	}

	set.seed(0)
	subcells <- lapply(names(reg_size), function(reg) {
		regc <- cells[int@meta.data[cells, "lobe"] == reg]
		if (length(regc) > maxcells){
			cc <- sample(regc, maxcells)
		} else {
			cc <- regc
		}
		return(cc)
		}) %>%
		unlist()

	subseu <- int[, subcells]
	
	# Generate expressed genes as well (for odds background)
	genes <- lapply(names(reg_size), function(reg) {
		regseu <- subset(subseu, lobe == reg)
		gg <- rownames(regseu$RNA)[Matrix::rowMeans(regseu$RNA@counts != 0) >= 0.1]
		gg
		}) %>%
		unlist() %>% unique()

	Idents(subseu) <- "lobe"
	res <- FindAllMarkers(subseu, max.cells.per.ident = maxcells, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1) %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(region = cluster) %>%
				mutate(cluster = cls)
	return(list(res = res, genes = genes))
	}) %>%
	setNames(., cls_use)


exp_list <- lapply(res_list, function(x) x$genes)
saveRDS(exp_list, file = paste0(inputdir, "Expgenes_res_GliaInN.rds"))



allres <- lapply(res_list, function(x) x$res) %>% do.call(rbind, .)
saveRDS(allres, file = paste0(inputdir, "DEG_res_GliaInN.rds"))



## Get expression ratios
int[["avgcls"]] <- paste0(int$subtype2, "|", int$lobe)
expr_int <- CalcMeanRatio(object = int, group.by = "avgcls", assay = "RNA")
saveRDS(expr_int, file = paste0(inputdir, "ExprRatio_res_GliaInN.rds"))


#all_genes <- rownames(mon)
#saveRDS(all_genes, file = paste0(inputdir, "All_monkey_annotation_genes.rds"))


###--------------------------------------------------------------------------------
## Combine

cbnres <- lapply(c("GliaInN", "RGC", "IPC-ExN"), function(type) {
	res <- readRDS(file = paste0(inputdir, "DEG_res_", type, ".rds"))
	res
	}) %>%
		do.call(rbind, .)
saveRDS(cbnres, file = paste0(inputdir, "DEG_res.rds"))














