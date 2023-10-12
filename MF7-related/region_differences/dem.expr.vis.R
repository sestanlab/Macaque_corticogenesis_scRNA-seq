library(cowplot)
library(dplyr)
library(ggplot2)
library(tibble)
library(Seurat)


###--------------------------------------------------------------------------------
## Combine the DEG res
load("../RGC_progression/load_files/Disease_genes_v3.Rdata") ## alltb, alllist
load(file = "./load_files/Reg-DEG_expr_avgs.Rdata")
## avgs, ratios, 
alltb <- alltb[intersect(rownames(alltb), rownames(avgs)), ]

degres <- readRDS(file = paste0("./load_files/", "DEG_res.rds")) %>%
			filter(gene %in% rownames(alltb)) %>%
			filter(p_val_adj <= 0.01 & avg_logFC >= 0.25 & ratio_fc >= 1.25 & pct.1 >= 0.2)
marres <- readRDS(file = paste0("./load_files/", "Marker_res.rds"))%>%
			filter(gene %in% rownames(alltb)) %>%
			filter(p_val_adj <= 0.01 & avg_logFC >= 0.25 & ratio_fc >= 1.25 & pct.1 >= 0.2)



## Generate dem martrices
allgps <- c("NESC", "vRG_early", "vRG_late", "oRG", "IPC EOMES NEUROG1", "IPC EOMES NHLH1", "ExN L6B", "ExN L6CT", "ExN upper", "MGE-InN", "CGE-InN", "gIPC", "aIPC", "oIPC", "Astro", "OPC")
dem_mats <- lapply(c("FC", "MSC", "OC", "TC"), function(lb) {
	dem_mat <- lapply(allgps, function(gp) {
		mars <- marres %>%
				filter(cluster %in% gp & region %in% lb) %>%
				filter(!grepl("^LOC", gene))
		degs <- degres %>%
				filter(cluster %in% gp & region %in% lb) %>%
				filter(!grepl("^LOC", gene))
		dem <- intersect(mars$gene, degs$gene)
		dem_value <- ifelse(rownames(alltb) %in% dem, 1, 0) %>%
				setNames(., rownames(alltb))
		return(dem_value)
		}) %>%
		setNames(., allgps) %>%
		as.data.frame(., check.names = FALSE) %>%
		as.matrix()
	colnames(dem_mat) <- paste0(lb, "|", colnames(dem_mat))
	return(dem_mat)
	}) %>%
	do.call(cbind, .)


## Top dems
demtops <- lapply(c("FC", "MSC", "OC", "TC"), function(lb) {
	lb_genes <- lapply(allgps, function(gp) {
		mars <- marres %>%
				filter(cluster %in% gp & region %in% lb) %>%
				filter(!grepl("^LOC", gene)) %>%
				top_n(25, wt = ratio_fc)
		degs <- degres %>%
				filter(cluster %in% gp & region %in% lb) %>%
				filter(!grepl("^LOC", gene)) %>%
				top_n(25, wt = ratio_fc)
		dem <- intersect(mars$gene, degs$gene)
		return(dem)
		}) %>%
		unlist() %>% unique()
	print(length(lb_genes))
	return(lb_genes)
	}) %>%
	setNames(., c("FC", "MSC", "TC", "OC"))
print(length(unique(unlist(demtops))))


## Plot gene expression
genes <- unlist(demtops) %>% unique()

source("./reg-deg.fun.R")
CirclePlot.horizontal(avg = avgs[rownames(alltb), ], ratio = ratios[rownames(alltb), ], features = genes, file_name = "Region_DEGs_disease_genes.pdf", dot.min = 0, dot.scale = 4, scale.by = "radius", shape = 16, cluster.order = allgps, stroke.col = "black", stroke.size = 0.5, col.min = -2.5, col.max = 2.5, scale = TRUE, scale.min = NA, scale.max = NA, mask.matrix = dem_mats[rownames(alltb), ], return.plot = FALSE, width.scale = 0.5, height.base = 1.5, font.scale = c(0.9, 0.6), height.unit = 0.065)


dis_df <- read.csv(file = paste0("./load_files/", "DiseaseOrderFinal.csv"), stringsAsFactors = FALSE)
dis_ord <- dis_df$Disease.listname
top_anno <- alltb[genes, dis_ord]
top_anno <- top_anno[, colSums(top_anno) > 0]
colnames(top_anno) <- extract_field(colnames(top_anno), -1, "::")

## Disease gene annotation
pdf(paste0("./report/", "Region_DEGs_disease_genes_annot.pdf"), width = 4, height = 7)
pheatmap::pheatmap(top_anno, cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(c("#f0f0f0", "#000000"))(30), border_color = "white", show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 8, fontsize_row = 8, angle_col = 90)
dev.off()


###--------------------------------------------------------------------------

### Plot all DEM, not just top

###--------------------------------------------------------------------------
dems <- lapply(c("FC", "MSC", "OC", "TC"), function(lb) {
	lb_genes <- lapply(allgps, function(gp) {
		mars <- marres %>%
				filter(cluster %in% gp & region %in% lb) %>%
				filter(!grepl("^LOC", gene))# %>%
				#top_n(25, wt = ratio_fc)
		degs <- degres %>%
				filter(cluster %in% gp & region %in% lb) %>%
				filter(!grepl("^LOC", gene)) #%>%
				#top_n(25, wt = ratio_fc)
		dem <- intersect(mars$gene, degs$gene)
		return(dem)
		}) %>%
		unlist() %>% unique()
	print(length(lb_genes))
	return(lb_genes)
	}) %>%
	setNames(., c("FC", "MSC", "OC", "TC"))
genes2 <- unlist(dems) %>% unique()



CirclePlot.horizontal(avg = avgs[rownames(alltb), ], ratio = ratios[rownames(alltb), ], features = genes2, file_name = "Region_DEGs_disease_genes_all.pdf", dot.min = 0, dot.scale = 4, scale.by = "radius", shape = 16, cluster.order = allgps, stroke.col = "black", stroke.size = 0.5, col.min = -2.5, col.max = 2.5, scale = TRUE, scale.min = NA, scale.max = NA, mask.matrix = dem_mats[rownames(alltb), ], return.plot = FALSE, width.scale = 0.5, height.base = 1.5, font.scale = c(1, 0.6), height.unit = 0.08)


all_anno <- alltb[genes2, dis_ord]
all_anno <- all_anno[, colSums(all_anno) > 0]
colnames(all_anno) <- extract_field(colnames(all_anno), -1, "::")

## Disease gene annotation
pdf(paste0("./report/", "Region_DEGs_disease_genes_all_annot.pdf"), width = 4, height = 12)
pheatmap::pheatmap(all_anno, cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(c("#f0f0f0", "#000000"))(30), border_color = "white", show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 8, fontsize_row = 8, angle_col = 90)
dev.off()








