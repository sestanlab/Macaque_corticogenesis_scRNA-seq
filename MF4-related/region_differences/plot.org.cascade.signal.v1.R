source("../scripts/nhpf.fun.R")
source("./ptime.fun.v2.R")



## Genearte region order based on expression
load(file = paste0(inputdir, "Smooth_by_region_oRG.Rdata"))  ##org_smt, org_meta, 

## Load the enrichment results
load(file = paste0(inputdir, "Region_DEGs_res_v2.rds"))
# vrg_deg, org_deg, vrg_res, org_res, cbn_res
res_use <- org_res %>%
			mutate(cluster_idx = as.numeric(factor(as.character(cluster), levels = c("oRG HOPX TNC", "oRG HOPX APOE"))))
meta_use <- org_meta
smt_use <- org_smt


timeres <- readRDS(file = paste0(inputdir, "Order_by_EXPR_oRG_v4.rds"))
order_genes <- lapply(names(timeres), function(reg) {
	xx <- timeres[[reg]]
	maxval <- max(meta_use$pseudotime[meta_use$region == reg])

	xx$time.on[xx$time.on < 0] <- -0.001
	xx$time.off[xx$time.off > ceiling(maxval)] <- ceiling(maxval) + 1

	## Add the cell type each gene enriched in
	sub_res <- res_use %>%
				filter(region == reg)
	xx$cluster_idx <- sapply(1:nrow(xx), function(idx) min(sub_res$cluster_idx[sub_res$gene == xx$gene[[idx]]])) 

	xx$start.bin <- as.numeric(cut(xx$time.on, 8))
	xx$end.bin <- as.numeric(cut(xx$time.off, 8))

	yy <- xx$gene[order(xx$start.bin, xx$end.bin, xx$time.on, xx$time.off, decreasing = FALSE)] ##xx$cluster_idx, 
	yy
	}) %>%
	setNames(., names(timeres)) %>%
	.[c("FC", "MSC", "TC", "OcC")]



## Load signaling-related genes & intersect
path_gset <- readRDS("/home/sm2726/project/cortex_development/dorsal_ana/genecluster/load_files/5.go/Pathway_gset.slim.rds")
all_siggenes <- unlist(path_gset) %>% unique()
order_genes <- lapply(order_genes, function(x) intersect(x, all_siggenes))
order_genes <- order_genes[sapply(order_genes, length) > 0]
all_visgenes <- unlist(order_genes) %>% unique()


## Generate scaled average expression & change row names
svg <- smt_use[all_visgenes, ,drop = FALSE] %>%
            as.matrix() %>%
            t() %>% scale() %>% t() %>%
            MinMax(., min = -1.5, max = 2.5)
svg[is.na(svg)] <- -1.5
exp_mat <- lapply(names(order_genes), function(reg) {
	mat <- as.matrix(svg[order_genes[[reg]], ,drop = FALSE])
	rownames(mat) <- paste0(reg, "|", rownames(mat))
	return(mat)
	}) %>%
	do.call(rbind, .)

rsplit <- extract_field(rownames(exp_mat), 1, "|") %>%
			setNames(., NULL) %>%
			factor(., levels = c("FC", "MSC", "TC", "OcC"))


## Plot
plot_heatmap.RGCcascade(mat = exp_mat, meta = meta_use, label_genes = NULL, color_breaks = seq(-1.5, 2.5, 0.5), file_name = paste0("Expr_byExpr_oRG_Signaling"), pdf_height = 6, row_split = rsplit, highlight_genes = NULL, fontsize = 6)



##------------------------------------------------------------------------------------------
## Plot annotation
cls_cols <- c("#e25a9a", "#1ac9a9", "#a3c10d", 
	"#efbc88", "#dbb75c", "#6700e5", "#7133aa", "#003563"
	) %>%
			setNames(., c("NEP RSPO3", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", 
				"oRG HOPX TNC", "oRG HOPX APOE", 
				"tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))
sel_cls <- c("oRG HOPX TNC", "oRG HOPX APOE")
reggene_ord <- lapply(names(order_genes), function(reg) paste0(reg, "|", order_genes[[reg]])) %>%
				unlist()
ctp_annot <- res_use %>%
				mutate(cluster = factor(as.character(cluster), levels = sel_cls)) %>%
				mutate(reggene = paste0(region, "|", gene)) %>%
				filter(reggene %in% reggene_ord) %>%
				mutate(reggene = factor(reggene, levels = rev(reggene_ord))) %>%
				mutate(region = factor(region, levels = c("FC", "MSC", "TC", "OcC"))) %>%
				mutate(fillcol = cls_cols[as.character(cluster)])
p1 <- ggplot(ctp_annot, aes_string(x = "cluster", y = "reggene", fill = "fillcol")) + 
			geom_tile(color = NA, size = 1) +
			##scale_fill_manual(values = c(`0` = "lightgrey", `1` = "black")) +
			scale_fill_identity() +
			theme_classic() + 
			scale_y_discrete(labels = extract_field(rev(reggene_ord), 2, "|")) +
			RotatedAxis() + 
			##facet_grid(rows = vars(region), space = "free_y", scales = "free_y") +
			theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 8, face = "italic"), axis.text.x = element_text(size = 10))
pdf(paste0(outputdir, "Expr_byExpr_oRG_Signaling_v2.annot-region.pdf"), width = 4, height = 6)
print(p1)
dev.off()



## Plot pathway annotation
anno_mat <- lapply(names(path_gset), function(sig) {
	genes <- setNames(extract_field(reggene_ord, 2, "|"), NULL)
	value <- setNames(ifelse(genes %in% path_gset[[sig]], 1, 0), reggene_ord)
	value
	}) %>%
	setNames(., names(path_gset)) %>%
	as.data.frame(., check.names = FALSE) %>%
	rownames_to_column("reggene") %>%
	tidyr::gather(., "pathway", "value", names(path_gset)) %>%
	mutate(pathway = factor(pathway, levels = c("BMP", "EPH", "FGF", "NOTCH", "RA", "WNT"))) %>%
	mutate(reggene = factor(reggene, levels = rev(reggene_ord))) %>%
	mutate(fillcol = ifelse(value == 1, "black", "lightgrey"))

p2 <- ggplot(anno_mat, aes_string(x = "pathway", y = "reggene", fill = "fillcol")) + 
			geom_tile(color = "white", size = 1) +
			scale_fill_identity() +
			theme_classic() + 
			scale_y_discrete(labels = extract_field(rev(reggene_ord), 2, "|")) +
			RotatedAxis() + 
			##facet_grid(rows = vars(region), space = "free_y", scales = "free_y") +
			theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 8, face = "italic"), axis.text.x = element_text(size = 10))
pdf(paste0(outputdir, "Expr_byExpr_oRG_Signaling_v2.annot-region.pdf"), width = 6, height = 6)
plot_grid(p1, p2, align = "h") %>% print()
dev.off()








