source("../scripts/nhpf.fun.R")
library(CellChat) 
library(patchwork) 
library(ggpubr)
library(circlize)
options(stringsAsFactors = FALSE)
library(tidyr) ## gather
library(ggsankey) 



res <- readRDS(file = paste0(inputdir, "/old/Cellchat_res_custom.rds"))
##res <- readRDS(file = paste0(inputdir, "/Cellchat_res_custom.rds"))


interactions <- read.csv(file = paste0("./load_files/", "cellchat_custom/", "interaction_input_CellChatDB.csv"), row.names = 1)
interactions["WNT5A_ROR2", "interaction_name"] <- "WNT5A_ROR2"
inter_pre <- rownames(res@LR$LRsig) %>%
				intersect(., interactions$interaction_name)



## Set interaction pairs
pc_clusters <- c("PC FGF17", "PC NKX2-1", "PC RSPO3", "PC TTR")
rgc_clusters <- c("FC NERG-early", "GE NERG-early", "OcC NERG-early")
cls_pairs <- expand.grid(rgc_clusters, pc_clusters) %>%
                mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
                subset(Var1 != Var2) %>%
        mutate(pair = paste0(Var2, "|", Var1)) %>%
        .$pair



## Get the interaction df
probs <- res@net$prob
pvals <- res@net$pval

mats <- lapply(inter_pre, function(ii) {
	mat <- pvals[,,ii]; probmat <- probs[,,ii]
	mat[mat > 0.05] <- 1
	mat[probmat < 1e-6] <- 1
	diag(mat) <- NA

	vec <- reshape2::melt(mat, value.name = "prob") %>%
				setNames(., c("rowcls", "colcls", "prob")) %>%
				mutate(rowcls = as.character(rowcls), colcls = as.character(colcls)) %>%
				mutate(pair = paste0(rowcls, "|", colcls)) %>%
				filter(!is.na(prob)) %>%
				column_to_rownames("pair") %>%
				.[cls_pairs, "prob"]
	vec
	}) %>%
		setNames(., inter_pre) %>%
		as.data.frame(., check.names = FALSE) %>%
		t() %>%
		as.matrix()
colnames(mats) <- cls_pairs
mats <- mats[rowSums(mats <= 0.05) > 0, ]


inter_final <- rownames(mats)
sub_anno <- interactions %>%
				subset(interaction_name %in% inter_final) %>%
				filter(annotation != "ECM-Receptor") %>%
				select(interaction_name, pathway_name, pathway_name, ligand, receptor, evidence)
mats <- mats[rownames(sub_anno), ]


## Load cellphone DB results
load(file = paste0(inputdir, "Cellphone_filtered_res.rds"))
##save(pval_cpb, mean_cpb, meta_cpb, 


## use setdiff(meta_cpb$interacting_pair, sub_anno$interaction_name) to get the cellphoneDB exclusive pairs
sigmeta <- read.table(file = paste0(inputdir, "cpb.exclusive.pairs.txt"), header = TRUE)
rownames(sigmeta) <- sigmeta$interaction_name
sigmeta$ligand <- extract_field(sigmeta$interaction_name, 1, "_")
sigmeta$receptor <- extract_field(sigmeta$interaction_name, "rm_start", "_")
sigmeta$annotation <- ifelse(meta_cpb[rownames(sigmeta), "secreted"] == "True", "Secreted Signaling", "Cell-Cell Contact")
sigmeta$evidence <- "SigCellphoneDB"


pval_cpb[pval_cpb > 0.05] <- 1
newmats <- rbind(mats, pval_cpb[sigmeta$interaction_name, colnames(mats), drop = FALSE])
newmeta <- rbind(sub_anno, sigmeta[, colnames(sub_anno)])
newmeta <- newmeta[rownames(newmats), ]



## Do PCA on the matrix
smat <- newmats %>%
			t() %>% scale()
pcres <- RunPCA(smat, npcs = 10, approx = FALSE)
new_path <- newmeta[rownames(newmats), "pathway_name"]
rare_path <- table(new_path) %>% .[. <= 3] %>% names() ##%>%
				##setdiff(., "BMP")
new_path[new_path %in% rare_path] <- "rare"
cc <- RSKC::RSKC(newmats, ncl = 10, alpha = 0)

pdata <- data.frame(pairs = rownames(newmats),
					pathway = new_path,
					stringsAsFactors = FALSE) %>%
				cbind(., pcres@cell.embeddings[, 1:5])
pdata$cluster <- as.character(cc$labels)
plist <- lapply(paste0("PC_", 2:5), function(xx) {
	p <- ggplot(pdata, aes_string(x = "PC_1", y = xx, color = "cluster")) +
				geom_point(size = 3) +
				theme_bw() +
				theme(axis.line = element_line(size = 0.2), axis.ticks = element_line(0.2))
	p
	})
pdf(paste0(outputdir, "LR_PCs.pdf"), width = 16, height = 18)
patchwork::wrap_plots(plist, nrow = 2, ncol = 2, guides = "collect") & theme(legend.position = "bottom")
dev.off()


tsneres <- RunTSNE(t(smat), check_duplicates= FALSE)
tsnedata <- data.frame(pairs = colnames(smat),
					rawpath = newmeta[rownames(newmats), "pathway_name"],
					pathway = new_path,
					stringsAsFactors = FALSE) %>%
				cbind(., tsneres@cell.embeddings)
tsnedata$cluster <- as.character(cc$labels)

if (FALSE){
## overlapped points
all_dups <- lapply(1:nrow(tsnedata), function(idx) {
	other_idx <- setdiff(1:nrow(tsnedata), idx)
	dups <- other_idx[sapply(other_idx, function(idx2) {
		xdif <- abs(tsnedata$tSNE_1[idx] - tsnedata$tSNE_1[idx2])
		ydif <- abs(tsnedata$tSNE_2[idx] - tsnedata$tSNE_2[idx2])
		test <- ifelse(xdif <= 0.1 & ydif <= 0.1, TRUE, FALSE)
		return(test)
		})]
	if (length(dups) > 0){
		return(union(idx, dups))
	} else {
		return(NULL)
	}
	}) %>%
	unlist() %>% unique()
}

## Add randomization to the coordinates
set.seed(42)
tsnedata$tSNE_1 <- tsnedata$tSNE_1 + rnorm(nrow(tsnedata), mean = 0, sd = 0.2)
tsnedata$tSNE_2 <- tsnedata$tSNE_2 + rnorm(nrow(tsnedata), mean = 0, sd = 0.2)


p2 <- ggplot(tsnedata, aes_string(x = "tSNE_1", y = "tSNE_2", color = "pathway")) +
				geom_jitter(size = 3.5, shape = 16) +
				theme_bw() +
				theme(axis.line = element_line(size = 0.2), axis.ticks = element_line(0.2), legend.position = "bottom")
p3 <- ggplot(tsnedata, aes_string(x = "tSNE_1", y = "tSNE_2", color = "cluster")) +
				geom_jitter(size = 3.5, shape = 16) +
				theme_bw() +
				theme(axis.line = element_line(size = 0.2), axis.ticks = element_line(0.2), legend.position = "bottom")
pdf(paste0(outputdir, "LR_TSNE.pdf"), width = 8, height = 5, useDingbats = FALSE)
patchwork::wrap_plots(list(p2, p3), nrow = 1, ncol = 2, guides = "collect") & theme(legend.position = "bottom")
dev.off()



## Plot all pairs
cls_ord <- tsnedata$cluster %>%
				unique() %>%
				as.numeric() %>% sort() %>% as.character()
cls_cols <- viridis(length(cls_ord)) %>% setNames(., cls_ord)
p_ord <- split(tsnedata$pairs, tsnedata$cluster) %>%
			.[cls_ord] %>% unlist()
dotdata <- -log10(newmats) %>%
			as.matrix() %>%
			reshape2::melt() %>% 
			setNames(., c("ID", "cluster_pair", "mlogp")) %>%
			mutate(mlogp = MinMax(mlogp, min = 0, max = 4)) %>%
			mutate(ID = factor(as.character(ID), levels = p_ord)) %>%
			mutate(cluster_pair = factor(as.character(cluster_pair), levels = rev(cls_pairs)))
dotdata$mlogp[dotdata$mlogp == 0] <- NA
anno <- data.frame(cluster = tsnedata$cluster[match(p_ord, tsnedata$pairs)],
						x_left = 0:(length(p_ord)-1),
						x_right = 1:length(p_ord),
						y_min = 13,
						y_max = 13.5,
						stringsAsFactors = FALSE) %>%
						mutate(fill = cls_cols[cluster]) %>%
						mutate(x_left = x_left + 0.5, x_right = x_right + 0.5)
pdot <- ggplot(dotdata, aes_string(y = "cluster_pair", x = "ID", size = "mlogp", color = "mlogp"))+
			geom_point(shape = 16)+
			annotate("rect", xmin = anno$x_left, xmax = anno$x_right, ymin = anno$y_min, ymax = anno$y_max, fill = anno$fill) +
			scale_size(range = c(0, 3)) +
			scale_color_gradient(low = "lightgrey", high = "red") +
			scale_fill_identity() +
			theme_bw()+
			##RotatedAxis() + 
			theme(panel.grid = element_line(size = 0.2), axis.line = element_line(size = 0.2), axis.ticks = element_line(0.2), legend.position = "bottom", axis.text.x = element_text(size = rel(0.6), angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank(), axis.text.y = element_text(size = rel(0.5)))
pdf(paste0(outputdir, "LR_all_interactions_new.pdf"), width = 10, height = 4, useDingbats= FALSE)
print(pdot)
dev.off()


## Plot interaction pair annotation
mmdata <- tsnedata[match(p_ord, tsnedata$pairs), ] %>%
			mutate(ID = factor(as.character(pairs), levels = p_ord))

annomat <- matrix(0, nrow = length(unique(tsnedata$rawpath)), ncol = length(p_ord), dimnames = list(levels(as.factor(tsnedata$rawpath)), p_ord))
for (i in 1:nrow(tsnedata)){
	annomat[tsnedata$rawpath[i], tsnedata$pairs[i]] <- 1
}
pdf(paste0(outputdir, "LR_all_interactions_annotation.pdf"), width = 10, height = 4)
pheatmap::pheatmap(annomat, cluster_rows =FALSE, cluster_cols = FALSE, color = colorRampPalette(c("#f0f0f0", "#000000"))(30), border_color = "white", show_rownames = TRUE, show_colnames = TRUE, fontsize_col = 6, fontsize_row = 6)
dev.off()



pp <- ggplot(mmdata, aes_string(y = "rawpath", x = "ID"))+
			geom_tile(fill = "black", color = "lightgrey")+
			theme_classic()+
			theme(panel.grid = element_line(size = 0.2), axis.line = element_line(size = 0.2), axis.ticks = element_line(0.2), legend.position = "bottom", axis.text.x = element_text(size = rel(0.6), angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank(), axis.text.y = element_text(size = rel(0.5)))
pdf(paste0(outputdir, "LR_all_interactions_annotation.pdf"), width = 10, height = 4)
print(pp)
dev.off()



## Plot eigengenes (average p values)
newmlp <- MinMax(-log10(newmats), min = 0, max = 4)
avgdata <- by(newmlp, paste0("c", tsnedata$cluster), colMeans) %>%
		do.call(rbind, .) %>%
		as.matrix() %>%
		reshape2::melt() %>% 
		setNames(., c("ID", "cluster_pair", "mlogp")) %>%
		mutate(cluster_pair = factor(as.character(cluster_pair), levels = cls_pairs)) %>%
		mutate(ID = factor(as.character(ID), levels = paste0("c", sort(as.numeric(unique(tsnedata$cluster)), decreasing = TRUE))))
pdot2 <- ggplot(avgdata, aes_string(x = "cluster_pair", y = "ID", size = "mlogp", color = "mlogp"))+
			geom_point(shape = 16)+
			scale_size(range = c(0, 4)) +
			scale_color_gradient(low = "lightgrey", high = "red") +
			scale_fill_identity() +
			theme_bw()+
			##RotatedAxis() + 
			theme(panel.grid = element_line(size = 0.2), axis.line = element_line(size = 0.2), axis.ticks = element_line(0.2), legend.position = "right", axis.text.x = element_text(size = rel(0.6), angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank(), axis.text.y = element_text(size = rel(0.8)))
pdf(paste0(outputdir, "LR_all_interactions_mes.pdf"), width = 3, height = 4, useDingbats= FALSE)
print(pdot2)
dev.off()







## Sankey plots showing the relationship between pathways and modules
source("./chat.fun.R")
cbn_meta <- tsnedata[, c("cluster", "pathway", "rawpath")]
rownames(cbn_meta) <- tsnedata$pairs
cbn_meta$cluster <- paste0("c", cbn_meta$cluster)

all_paths <- table(cbn_meta$pathway)
all_cls <- gsub("c", "", cbn_meta$cluster) %>%
				unique() %>%
				as.numeric() %>% sort() %>% paste0("c", .)

path_list <- list(p1 = c("COLLAGEN", "EPHA", "EPHB", "FGF", "JAM"), 
					p2 = c("LAMININ", "ncWNT", "NOTCH", "NRXN", "PTN", "WNT"), 
					p3 = c("rare"))
s1 <- plot_ggsankey(meta = cbn_meta[cbn_meta$pathway %in% path_list$p1, ], x_col = "pathway", x_ord = path_list$p1, next_col = "cluster", next_ord = intersect(all_cls, unique(cbn_meta$cluster[cbn_meta$pathway %in% path_list$p1])) %>% rev(), type = "sankey", space= 0.01,  width = 0.1)
s2 <- plot_ggsankey(meta = cbn_meta[cbn_meta$pathway %in% path_list$p2, ], x_col = "pathway", x_ord = path_list$p2, next_col = "cluster", next_ord = intersect(all_cls, unique(cbn_meta$cluster[cbn_meta$pathway %in% path_list$p2])) %>% rev(), type = "sankey", space= 0.01,  width = 0.1)
s3 <- plot_ggsankey(meta = cbn_meta[cbn_meta$pathway %in% path_list$p3, ], x_col = "rawpath", x_ord = NULL, next_col = "cluster", next_ord = intersect(all_cls, unique(cbn_meta$cluster[cbn_meta$pathway %in% path_list$p3])) %>% rev(), type = "sankey", space= 0.01,  width = 0.1)
pdf(paste0(outputdir, "LR.sankey.S1.pdf"), width = 10, height = 10)
print(s1)
print(s2)
print(s3)
dev.off() 




## Plot Sankey for just some selected pathways
source("./chat.fun.R")

sel_paths <- c("FGF", "EPHA", "EPHB", "NOTCH", "WNT", "ncWNT", "RSPO")

cbn_meta <- tsnedata[, c("cluster", "pathway", "rawpath")]
rownames(cbn_meta) <- tsnedata$pairs
cbn_meta$cluster <- paste0("c", cbn_meta$cluster)
cbn_meta <- cbn_meta[cbn_meta$rawpath %in% sel_paths, ]
all_cls <- gsub("c", "", cbn_meta$cluster) %>%
				unique() %>%
				as.numeric() %>% sort() %>% paste0("c", .)

ss <- plot_ggsankey(meta = cbn_meta, x_col = "rawpath", x_ord = sel_paths, next_col = "cluster", next_ord = rev(all_cls), type = "sankey", space= 0.01,  width = 0.1)
pdf(paste0(outputdir, "LR.sankey.MF1.pdf"), width = 4, height = 7)
print(ss)
dev.off()










## Individual circular plot
allcls <- c("OcC NERG-early", "PC RSPO3", "PC TTR", "GE NERG-early", "PC NKX2-1", "PC FGF17", "FC NERG-early")
cls.cols <- c("#89DA59", "#0c6e8c", "#6ba4f4", "#4000d4", "#dc59e5", "#d40063", "#FF420E") %>% 
				setNames(., allcls)


source("./chat.fun.R")
mat_vis <- newmats
pairs <- rownames(mat_vis)
for (pair in pairs){
	pdata <- PreCircleData(mat = mat_vis, LR = pair, log10transfer = TRUE, maxp = 3, minp = NULL)
	PlotInterCircle_individual(pdata = pdata, output_dir = paste0(outputdir, "sepcircle/"), file_name = paste0("LR-circle-sep-", pair), cluster.cols = cls.cols, cls.ord = allcls, title_lab = NULL)
}



















































