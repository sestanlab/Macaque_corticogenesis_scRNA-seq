source("~/project/PFC/scripts/pfc.fun.R")


## Calculate Region Fold Changes
rgc <- readRDS(file = paste0("./load_files/", "RGC_seu_for_DEG_analysis.rds"))


## Calculate Average expression in each region
all_regs <- c("FC", "MSC", "TC", "OcC")
subrgc <- subset(rgc, cluster2 %in% c("NEP RSPO3", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "oRG HOPX APOE", "oRG HOPX TNC"))
reg_avgs <- lapply(all_regs, function(reg) {
	print(paste0("Working on region:", reg))
	reg_seu <- subrgc
	reg_seu$newlobe <- ifelse(reg_seu$lobe == reg, reg, "bg")
	reg_seu$avgcls <- paste0(reg_seu$newlobe, "|", reg_seu$cluster2)

	Idents(reg_seu) <- "avgcls"
	avg <- as.matrix(AverageExpression(reg_seu, assay = "RNA")$RNA)
	return(avg)
	}) %>%
	setNames(., all_regs)


## Calculate log Fold changes in each region
ptval <- 0.1
fc_res <- lapply(all_regs, function(reg) {
	avg <- reg_avgs[[reg]]
	sel_cls <- c("NEP RSPO3", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "oRG HOPX APOE", "oRG HOPX TNC")
	que_cols <- paste0(reg, "|", sel_cls) %>%
				intersect(., colnames(avg))
	ref_cols <- gsub(paste0("^", reg, "\\|"), "bg|", que_cols) %>%
				intersect(., colnames(avg))
	fc <- (avg[, que_cols] + ptval)/(avg[, ref_cols] + ptval)
	colnames(fc) <- gsub(paste0("^", reg, "\\|"), "", que_cols)

	## Merge the FCs
	merge_list <- list(`vRG early` = c("NEP RSPO3", "vRG HMGA2 CCND1"),
						`vRG late` = c("vRG SAT1 STMN2"),
						`oRG` = c("oRG HOPX APOE", "oRG HOPX TNC"))
	merge_list <- merge_list[sapply(merge_list, function(x) sum(x %in% colnames(fc)) >= 1)] ## Incase TC doesn't have vRG early

	mergefc <- lapply(names(merge_list), function(cls) {
		nmat <- log(apply(fc[, intersect(merge_list[[cls]], colnames(fc)), drop = FALSE], 1, max))
		return(nmat)
		}) %>%
		setNames(., names(merge_list)) %>%
		as.data.frame(., check.names = FALSE) 
	colnames(mergefc) <- paste0("logFC|", colnames(mergefc))
	mergefc <- mergefc %>%
		rownames_to_column("gene") %>%
		mutate(region = reg)

	return(mergefc)
	}) %>%
	setNames(., all_regs)

save(reg_avgs, fc_res, file = paste0(inputdir, "Shared.region-markers.acrosssubtypes.avg.Rdata"))



##------------------------------------------------------------------------------------------
## First Identify the list of Shared region-specific genes across early and late RGC subtypes. 
allres <- readRDS(file = paste0(inputdir, "Region_DEGs_rawres.rds"))

## Consider two situations
## 1. (expr ratio >= 0.1)
slim_dex1 <- allres %>%
				filter(pct.1 >= 0.1 & pct.2 <= 0.75 & ratio_fc >= 1.4 & avg_logFC >= 0.2 & p_val_adj <= 0.001)
## 2. (expr ratio < 0.1 &  0.05)
slim_dex2 <- allres %>%
				filter(pct.1 >= 0.05 & pct.1 < 0.1 & ratio_fc >= 4 & avg_logFC >= 0.1 & p_val_adj <= 0.001)
## 3. 
slim_dex3 <- allres %>%
				filter(pct.2 > 0.75 & ratio_fc >= 1.1 & avg_logFC >= 1 & p_val_adj <= 0.001)


slim_dex <- rbind(slim_dex1, slim_dex2) %>%
				rbind(., slim_dex3)



load(file = paste0(inputdir, "Shared.region-markers.acrosssubtypes.avg.Rdata"))
## reg_avgs, fc_res, 



## Shared DEGs
all_regs <- c("FC", "MSC", "TC", "OcC")
share_sigs <- lapply(all_regs, function(reg) {
	reg_dex <- slim_dex %>%
			filter(cluster %in% c("NEP RSPO3", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "oRG HOPX APOE", "oRG HOPX TNC") & region == reg) %>%
			mutate(group = case_when(
				cluster %in% c("NEP RSPO3", "vRG HMGA2 CCND1") ~ "vRG early", 
				cluster %in% c("vRG SAT1 STMN2") ~ "vRG late", 
				cluster %in% c("oRG HOPX APOE", "oRG HOPX TNC") ~ "oRG"
				)) %>%
			group_by(group, gene) %>%
			top_n(1, wt = avg_logFC) %>%
			ungroup() %>%
			group_by(gene) %>%
			mutate(nhits = n()) %>%
			ungroup()

	all_genes <- unique(reg_dex$gene)
	enr_mat <- matrix(NA, nrow = length(all_genes), ncol = 3, dimnames = list(all_genes, c("vRG early", "vRG late", "oRG")))
	for (ii in c("vRG early", "vRG late", "oRG")) {
		enr_mat[, ii] <- sapply(all_genes, function(gg) {
			value <- ifelse(ii %in% reg_dex$group[reg_dex$gene == gg], 1, 0)
			value
			})
	}

	reg_fc <- fc_res[[reg]] %>%
				filter(gene %in% rownames(enr_mat))


	corres <- reg_fc %>%
				select(-region) %>%
				column_to_rownames("gene") %>%
				cor(., method = "p")
	print(corres)


	reg_fc <- reg_fc[match(rownames(enr_mat), reg_fc$gene), ]
	rownames(enr_mat) <- NULL
	df <- cbind(reg_fc, enr_mat) %>%
			as.data.frame(., check.names= FALSE) %>%
			mutate(Shared = ifelse(rowSums(enr_mat) == 0, 1, 0)) %>%
			tidyr::gather(., "latetype", "xlogFC", intersect(colnames(reg_fc), paste0("logFC|", c("vRG late", "oRG"))))

	colnames(df) <- gsub("\\|", "..", colnames(df)) %>%
			gsub(" ", ".", .)

	df$radius <- 0.05
	df
	}) %>%
	setNames(., all_regs)


source("~/project/PFC/MF7_contact/pie.fun.R")
source("./ptime.fun.v2.R")


p1 <- PlotWeightedScatterPie(pie.data = share_sigs[["FC"]], x.col = "logFC..vRG.early", y.col = "xlogFC", r.col = "radius", cls_use = c("vRG.early", "vRG.late", "oRG"), rsf = 1, scale.expression = FALSE, fc_limit = c(-1, 2.5))
p2 <- PlotWeightedScatterPie(pie.data = share_sigs[["MSC"]], x.col = "logFC..vRG.early", y.col = "xlogFC", r.col = "radius", cls_use = c("vRG.early", "vRG.late", "oRG"), rsf = 1, scale.expression = FALSE, fc_limit = c(-1, 2.5))
p3 <- PlotWeightedScatterPie(pie.data = share_sigs[["OcC"]], x.col = "logFC..vRG.early", y.col = "xlogFC", r.col = "radius", cls_use = c("vRG.early", "vRG.late", "oRG"), rsf = 1, scale.expression = FALSE, fc_limit = c(-1, 2.5))

pdf(paste0(outputdir, "Shared_region_markers.weightedPie.", "all", ".pdf"), width = 12, height = 8)
##plot <- patchwork::wrap_plots(list(p1, p2, p3), nrow = 1, ncol = 3)
plot <- plot_grid(p1, p2, p3, nrow = 1, ncol = 3)
print(plot)
dev.off()



## FC
#                logFC|vRG early logFC|vRG late  logFC|oRG
#logFC|vRG early      1.00000000     0.06888459 0.06920558
#logFC|vRG late       0.06888459     1.00000000 0.73750692
#logFC|oRG            0.06920558     0.73750692 1.00000000


## MSC
#                logFC|vRG early logFC|vRG late  logFC|oRG
#logFC|vRG early       1.0000000     -0.2286106 -0.3021543
#logFC|vRG late       -0.2286106      1.0000000  0.5248541
#logFC|oRG            -0.3021543      0.5248541  1.0000000


## OC
#                logFC|vRG early logFC|vRG late  logFC|oRG
#logFC|vRG early       1.0000000     -0.2750821 -0.2348342
#logFC|vRG late       -0.2750821      1.0000000  0.2332785
#logFC|oRG            -0.2348342      0.2332785  1.0000000







