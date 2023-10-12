source("../scripts/nhpf.fun.R")


seu <- readRDS(file = paste0(inputdir, "ExN_all_harmony_filtered_nonIT.rds"))	


## Remove cycling cells to avoid bias
load(file = paste0("./load_files/intermediate/", "IE_curve_harmony_IT_ptime.Rdata")) ##it_ptime, it_cycle
load(file = paste0("./load_files/intermediate/", "IE_curve_harmony_nonIT_ptime.Rdata")) ##nit_ptime, nit_cycle
cyc_cells <- union(it_cycle, nit_cycle)

seu <- seu[, setdiff(colnames(seu), cyc_cells)]
seu$pseudotime <- nit_ptime[colnames(seu), 1]



## Assign cells to bins
min_pt <- 0
max_pt <- quantile(seu$pseudotime, 0.999, na.rm = TRUE)
nPoints <- 50
half_inter <- 0.25##(max_pt - min_pt)/(2 * (nPoints - 1))
cut_bks <- c(seq(min_pt, max_pt, half_inter * 2) - half_inter, Inf)
seu$bin <- cut(seu$pseudotime, breaks = cut_bks) %>% as.numeric()



## There are too many ExN_deep_nascent cells
set.seed(0)
rm_cells1 <- sample(colnames(seu)[seu$subtype %in% c("ExN_deep_nascent") & seu$bin %in% 21:31], 5000)
rm_cells2 <- sample(colnames(seu)[seu$subtype %in% c("ExN_deep_SYT6") & seu$bin %in% 48:53], 2300)
rm_cells3 <- sample(colnames(seu)[seu$subtype %in% c("ExN_deep_NR4A2_GRID2") & seu$bin %in% 45:47], 200)
nseu <- seu[, setdiff(colnames(seu), unique(c(rm_cells1, rm_cells2, rm_cells3)))]
nseu <- subset(nseu, lobe != "Insula")

saveRDS(nseu, file = paste0(inputdir, "ExN_sliding-window_nonIT.rds"))



##------------------------------------------------------------------------
## Get bin information
nseu <- readRDS(file = paste0(inputdir, "ExN_sliding-window_nonIT.rds"))


df <- aggregate(pseudotime ~ bin, nseu@meta.data, mean)
df$cluster <- sapply(1:nrow(df), function(xx) {
	time <- df$pseudotime[[xx]]
	lin_meta <- nseu@meta.data 

	high_cells <- lin_meta[lin_meta$pseudotime > time, ] %>%
					arrange(pseudotime) %>%
					.$subtype %>% .[1:25]
	low_cells <- lin_meta[lin_meta$pseudotime < time, ] %>%
					arrange(desc(pseudotime)) %>%
					.$subtype %>% .[1:25]
	cls <- c(high_cells, low_cells) %>%
				table() %>% sort() %>% 
				rev() %>% .[1] %>% names()
	cls
	})
cls_cols <- setNames(c("#f1a340", "#bf812d", 
                    "#b35806", "#e89bc4", "#de77ae", "#c51b7d", "#D3D3D3", "#D3D3D3",
                    "#b35806", "#91ebe2", "#48a1e0", "#0868ac", "#c6dbef", "#c6dbef", "#4eb3d3"),
                    c("IPC EOMES VIM", "IPC_EOMES_NEUROG1", 
                    "IPC_EOMES_NHLH1_up", "ExN_up_nascent", "ExN_up_ADRA2A", "ExN_up_ACTN2", "ExN_PCC_R4A3", "ExN_up_KCNV1",
                    "IPC_EOMES_NHLH1_deep", "ExN_deep_nascent", "ExN_deep_KIF26A", "ExN_deep_SYT6", "ExN_deep_OPRK1_SULF1","ExN_deep_OPRK1_NR4A2", "ExN_deep_NR4A2_GRID2"))
anno_df <- df %>%
			mutate(x_left = c(0, df$pseudotime[1:(nrow(df) - 1)]),
					x_right = df$pseudotime,
					y_min = 0,
					y_max = 1,
					yy = 1.5) %>%
			mutate(fill = cls_cols[cluster])
p <- ggplot(anno_df, aes_string(x = "pseudotime", y = "yy", color = "cluster"))+
			geom_point(size = 1, shape = 16) +
			annotate("rect", xmin = anno_df$x_left, xmax = anno_df$x_right, ymin = anno_df$y_min, ymax = anno_df$y_max, fill = anno_df$fill) +
			scale_color_manual(values = cls_cols) +
			scale_fill_identity() +
			theme_classic()
pdf(paste0(outputdir, "ExN_sliding-window_nonIT_xaxis_lable.pdf"), width = 6, height = 4)
print(p)
dev.off()




CalcDEGRes <- function(object, group.by, max.cells = 100, assay = "RNA") {
	all_regs <- levels(as.factor(object@meta.data[, group.by]))
	complist <- lapply(all_regs, function(x) {
		bg_cls <- setdiff(all_regs, x)
		paste0(x, "--vs--", bg_cls)
		}) %>% unlist()
	Idents(object) <- group.by
	DefaultAssay(object) <- assay
	allres <- lapply(complist, function(xx) {
		cls1 <- strsplit(xx, "--vs--", fixed = TRUE)[[1]][1]
		cls2 <- strsplit(xx, "--vs--", fixed = TRUE)[[1]][2]
		res <- FindMarkers(object, ident.1 = cls1, ident.2 = cls2, only.pos = TRUE, max.cells.per.ident = max.cells, min.pct = 0.1, logfc.threshold = 0.2) %>%
				rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(region1 = cls1, region2 = cls2)
		res
		}) %>%
			do.call(rbind, .)
	return(allres)
}



allbins <- 1:53
allregs <- c("FC", "MSC", "TC", "OcC")
allres <- lapply(allbins, function(bb) {
	msize <- min(table(nseu$bin, nseu$lobe)[bb, ])
	if (msize < 100){
		if (bb == 1){
			cc <- c(bb, bb + 1)
		} else if (bb == 47){
			cc <- c(bb - 1, bb)
		} else {
			cc <- c(bb - 1, bb)
		}
	} else {
		cc <- bb
	}
	print(paste0("Bin-", bb, " has min-size as : ", msize))
	sub_cbn <- nseu[, nseu@meta.data$bin %in% cc]


	subcells <-	lapply(allregs, function(reg) {
			cells <- colnames(sub_cbn)[sub_cbn@meta.data$lobe == reg]
			if (length(cells) > 100){
				cells <- sample(cells, 100)
			} 
			return(cells)
		}) %>%
			unlist()

	res <- CalcDEGRes(object = sub_cbn[, subcells], group.by = "lobe", max.cells = 100, assay = "RNA") %>%
			mutate(bin = bb)
	return(res)
	}) %>%
		do.call(rbind, .)

saveRDS(allres, file = paste0(inputdir, "ExN_sliding-window_nonIT_DEGres.rds"))



##------------------------------------------------------------------------
## Visualize the results in 2D plots
allbins <- 1:53
slim_res <- allres %>%
				mutate(regionpair = paste0(region1, "|", region2)) %>%
				subset(pct.1 >= 0.1 & pct.2 <= 0.75 & ratio_fc >= 1.2 & avg_logFC >= 0.2) %>%
				mutate(binchar = as.character(bin)) %>%
				group_by(binchar, regionpair) %>% 
				mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
				filter(p_val_adj <= 0.01) %>%
				ungroup()
reg_div <- lapply(c("FC", "MSC", "TC", "OcC"), function(reg) {
	subdata <- slim_res %>%
				filter(region1 == reg | region2 == reg) %>%
				group_by(binchar, regionpair) %>%
				summarize(ndex = n(), bin = unique(bin))
	subdata$uniq_pair <- sapply(1:nrow(subdata), function(x) {
		paste(sort(strsplit(subdata$regionpair[x], "|", fixed = TRUE)[[1]]), collapse = "|")
		})
	subdata <- subdata %>%
				ungroup() %>% group_by(binchar, uniq_pair) %>%
				summarize(ndex = sum(ndex)) %>%
				ungroup() %>% group_by(binchar) %>%
				summarize(mdex = mean(ndex), bin = as.integer(unique(binchar)), region = reg)
	return(subdata)
	}) %>%
		do.call(rbind, .)
p <- ggplot(reg_div, aes(x = bin, y = mdex, color = region, fill = region)) +
		geom_smooth(size = 1, method = "loess", span = 0.4) +
		geom_point(size = 1)+
		##geom_errorbar(aes(ymin=mdiv-sediv, ymax=mdiv+sediv), width=.2) +
		theme_classic() +
		theme(legend.position = "bottom")
pdf(paste0(outputdir, "ExN_sliding-window_nonIT_nDEGs_byregion.pdf"), width = 10, height = 10)
print(p)
dev.off()





###--------------------------------------------------------------------------------------
### Write out the region-specific divergence scores
### for 3D visualization in matlab
times <- aggregate(pseudotime ~ bin, nseu@meta.data, mean)
pt_min <- round(min(times$pseudotime), digits = 6)
if (pt_min < min(times$pseudotime)){
	pt_min <- pt_min + 5e-6
} 
pt_max <- round(max(times$pseudotime), digits = 6)
if (pt_max > max(times$pseudotime)){
	pt_max <- pt_max - 5e-6
}

ntime <- seq(pt_min, pt_max, length.out = 100) %>%
				sapply(., function(x) round(x, digits = 6))
div_df <- data.frame(row.names = as.character(ntime))
for (reg in c("FC", "MSC", "TC", "OcC")){
	subdata <- reg_div %>% 
				filter(region == reg) %>%
				mutate(pseudotime = times$pseudotime[match(bin, times$bin)]) %>%
				mutate(logndex = log2(mdex + 1))

	## The predicted df rownames will be the column names of the predicted data
	time_use <- ntime[ntime >= min(subdata$pseudotime) & ntime <= max(subdata$pseudotime)]
	predict_df <- data.frame(row.names = paste0(reg, "|", time_use), 
			pseudotime = time_use)	

	fitval <- loess(logndex ~ pseudotime, subdata, span = 0.6) %>%
				predict(., predict_df, se = FALSE) %>%
				MinMax(., min = 0, max = 2000)
	div_df[as.character(time_use), reg] <- fitval
}

rownames(div_df) <- as.character(1:nrow(div_df))
out_df <- t(div_df)
write.table(out_df, file = paste0(inputdir, "ExN_sliding-window_nonIT_nDEGs_log.txt"), row.names = TRUE, col.names = TRUE, quote = TRUE) 














