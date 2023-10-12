source("../scripts/nhpf.fun.R")
source("./ptime.fun.v2.R")



## Genearte region order based on expression
load(file = paste0(inputdir, "Smooth_by_region_oRG.Rdata"))  ##org_smt, org_meta, 
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



## Label genes
tfs <- read.csv(file = "~/project/public_data/MEME_db/custom_motif/DatabaseExtract_v_1.01.csv", stringsAsFactors = FALSE) %>%
			subset(Is.TF. == "Yes") %>%
			.$HGNC.symbol %>% 
			unique()
cus_genes <- readLines(paste0(inputdir, "genes2label.oRG.txt"))
reg_genes <- lapply(c("FC", "MSC", "TC", "OcC"), function(reg) {
	interest <- union(tfs, cus_genes)
	if (sum(interest %in% order_genes[[reg]]) > 0){
		return(paste0(reg, "|", intersect(order_genes[[reg]], interest)))
	}
	}) %>%
	setNames(., c("FC", "MSC", "TC", "OcC")) %>%
	lapply(., function(x) extract_field(x, 2, "|"))
sel_genes <- list(left = c(reg_genes[["FC"]], reg_genes[["TC"]]),
			right = c(reg_genes[["MSC"]], reg_genes[["OcC"]]))




## Highlight TFs as red text
high_genes <- lapply(c("FC", "MSC", "TC", "OcC"), function(reg) paste0(reg, "|", tfs)) %>%
	unlist()
rsplit <- extract_field(rownames(exp_mat), 1, "|") %>%
			setNames(., NULL) %>%
			factor(., levels = c("FC", "MSC", "TC", "OcC"))



## Plot
plot_heatmap.RGCcascade.regionshare(mat = exp_mat, meta = meta_use, label_genes = sel_genes, color_breaks = seq(-1.5, 2.5, 0.5), file_name = paste0("Expr_byExpr_oRG"), pdf_height = 6, row_split = rsplit, highlight_genes = high_genes, fontsize = 6)

##------------------------------------------------------------------------------------------
## Plot annotation
#"#fc19cf", 
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
				mutate(reggene = factor(reggene, levels = rev(reggene_ord))) %>%
				mutate(region = factor(region, levels = c("FC", "MSC", "TC", "OcC"))) %>%
				mutate(fillcol = cls_cols[as.character(cluster)])
p1 <- ggplot(ctp_annot, aes_string(x = "cluster", y = "reggene", fill = "fillcol")) + 
			geom_tile(color = NA, size = 1) +
			##scale_fill_manual(values = c(`0` = "lightgrey", `1` = "black")) +
			scale_fill_identity() +
			theme_classic() + 
			RotatedAxis() + 
			facet_grid(rows = vars(region), space = "free_y", scales = "free_y") +
			theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size = 8, face = "italic"), axis.text.x = element_text(size = 10))
pdf(paste0(outputdir, "Expr_byExpr_oRG_v2.annot-region.pdf"), width = 4, height = 6)
print(p1)
dev.off()


##------------------------------------------------------------------------------------------
## Save the results as table
allcls <- c("oRG HOPX TNC", "oRG HOPX APOE")
enr_mat <- matrix(0, nrow = nrow(exp_mat), ncol = length(allcls), dimnames = list(rownames(exp_mat), allcls))
for (ii in rownames(enr_mat)){
	cls <- ctp_annot$cluster[ctp_annot$reggene %in% ii]
	if (length(cls) >= 1){
		for (jj in cls){
			enr_mat[ii, jj] <- 1
		} 
	}
}



df <- data.frame(reg_gene = rownames(exp_mat),
			stringsAsFactors = FALSE) %>%
			mutate(region = extract_field(reg_gene, 1, "|"), gene = extract_field(reg_gene, 2, "|")) %>%
			mutate(isTF = ifelse(gene %in% tfs, "Y", "N")) %>%
			cbind(., enr_mat)
df$time.on <- NA
df$time.off <- NA
for (ii in 1:nrow(df)){
	reg <- df$region[ii]
	gene <- df$gene[ii]
	timedata <- timeres[[reg]]

	df$time.on[ii] <- round(timedata$time.on[timedata$gene %in% gene], digits = 5)
	df$time.off[ii] <- round(timedata$time.off[timedata$gene %in% gene], digits = 5)
	if (df$time.on[ii] > df$time.off[ii]){
		df$time.off[ii] <- 100
	}
}
df$time.on <- MinMax(df$time.on, min = -100, max = 100)
df$time.off <- MinMax(df$time.off, min = 0, max = 100)


df <- df[, c("region", "gene", "isTF", colnames(enr_mat), "time.on", "time.off")]
df$region <- gsub("^FC$", "FR", df$region) %>%
				gsub("MSC", "MS", .) %>%
				gsub("TC", "Tem", .) %>%
				gsub("OcC", "OC", .)

colnames(df) <- c("Enriched region", "Gene", "is TF?", paste0("Enriched in ", colnames(enr_mat)), "Pseudotime on", "Pseudotime off")
write.table(df, file = "./report/table_Sxx_NSC_progression_oRG.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)





