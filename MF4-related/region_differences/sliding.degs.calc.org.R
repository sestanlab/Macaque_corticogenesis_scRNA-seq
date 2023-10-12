library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)

inputdir <- "./load_files/"
outputdir <- "./report/"

source("./ptime.fun.v2.R")


## Load data
load(file = paste0(inputdir, "RGC_combine_oRG_seu.Rdata"))
## rgc_cbn



##------------------------------------------------------------------------
## Get bin information
df <- aggregate(pseudotime ~ bin, rgc_cbn@meta.data, mean)
df$cluster <- sapply(1:nrow(df), function(xx) {
	time <- df$pseudotime[[xx]]
	lin_meta <- rgc_cbn@meta.data 

	high_cells <- lin_meta[lin_meta$pseudotime > time, ] %>%
					arrange(pseudotime) %>%
					.$cluster %>% .[1:25]
	low_cells <- lin_meta[lin_meta$pseudotime < time, ] %>%
					arrange(desc(pseudotime)) %>%
					.$cluster %>% .[1:25]
	cls <- c(high_cells, low_cells) %>%
				table() %>% sort() %>% 
				rev() %>% .[1] %>% names()
	cls
	})
cls_cols <- c("#990939", "#e25a9a", "#fc19cf", "#1ac9a9", "#a3c10d", 
	"#efbc88", "#dbb75c", "#6700e5", "#7133aa", "#003563"
	) %>%
			setNames(., c("RGC FABP7 PMP22", 
				"NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", 
				"oRG HOPX TNC", "oRG HOPX APOE", 
				"tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))
anno_df <- df %>%
			mutate(x_left = c(9, df$pseudotime[1:(nrow(df) - 1)]),
					x_right = df$pseudotime,
					y_min = 0,
					y_max = 1) %>%
			mutate(fill = cls_cols[cluster])
p <- ggplot(anno_df, aes(x = "pseudotime", y = 1))+
			annotate("rect", xmin = anno_df$x_left, xmax = anno_df$x_right, ymin = anno_df$y_min, ymax = anno_df$y_max, fill = anno_df$fill) +
			theme_classic() +
			scale_x_continuous(breaks = seq(9, 19, 2))+
			scale_fill_identity()
pdf(paste0(outputdir, "RGC_slidingDEGs_oRG_xaxis_lable.pdf"), width = 6, height = 4)
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



allbins <- df$bin
allregs <- levels(as.factor(rgc_cbn$lobe))
allres <- lapply(allbins, function(bb) {
	print(paste0("Working on bin: ", bb))
	msize <- min(table(rgc_cbn$bin, rgc_cbn$lobe)[as.character(bb), ])
	cc <- c(bb, bb + 1)
	sub_cbn <- rgc_cbn[, rgc_cbn@meta.data$bin %in% cc]


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

saveRDS(allres, file = paste0(inputdir, "RGC_slidingDEGs_oRG.rds"))



##------------------------------------------------------------------------
allbins <- df$bin
allres <- readRDS(file = paste0(inputdir, "RGC_slidingDEGs_oRG.rds"))
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
pdf(paste0(outputdir, "RGC_diffs_nDEGs_oRG_byregion.pdf"), width = 10, height = 10)
print(p)
dev.off()



 


###--------------------------------------------------------------------------------------
### Write out the region-specific divergence scores
times <- aggregate(pseudotime ~ bin, rgc_cbn@meta.data, mean)
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
write.table(out_df, file = paste0(outputdir, "RGC_diffs_nDEGs_oRG_log_v2.txt"), row.names = TRUE, col.names = TRUE, quote = TRUE) 





























