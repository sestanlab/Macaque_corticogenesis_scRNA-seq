library(dplyr)
library(tibble)

##-----------------------------------------------------------------------------------------------------
## Summarize enrichment models
piei <- readRDS(file = paste0("./load_files/", "Region_DEGs_pairwise_test_enrichpatterns.rds")) %>%
			mutate(model_idx = as.numeric(factor(model, levels = c("m1", "m2", "m3", "m1-n3", "m3-n3", "m1-n2")))) %>%
			group_by(feature, cluster) %>%
			filter(model_idx == max(model_idx))

sel_cls <- c("NEP RSPO3", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "tRG CRYAB MEST")
slimp <- piei %>%
			ungroup() %>%
			mutate(regsum = FC + MSC + TC + OcC) %>%
			filter(cluster %in% sel_cls & regsum != 0 & regsum != 1) %>%
			select(-c(model_idx, regsum)) %>%
			mutate(pattern = paste0(FC, "-", MSC, "-", TC, "-", OcC))# %>%
			#mutate(clusteridx = as.numeric(factor(cluster, levels = sel_cls))) %>%
			##group_by(pattern) %>%
			##filter(clusteridx == min(clusteridx)) %>%
			##select(-clusteridx)
##all_pats <- table(slimp$pattern) %>% .[. >= 1] %>% names()
saveRDS(slimp, file = paste0("./load_files/", "Region_DEGs_pairwise_test_summarizedpatterns.rds"))



##-----------------------------------------------------------------------------------------------------
## Further remove genes with low associationns
slimp <- readRDS(file = paste0("./load_files/", "Region_DEGs_pairwise_test_summarizedpatterns.rds"))
load(file = paste0("../region_differences/load_files/", "Smooth_by_region_tRG.Rdata"))  ##trg_smt, trg_meta, 


all_pats <- levels(as.factor(slimp$pattern))
all_regs <- c("FC", "MSC", "TC", "OcC")
filter_genes <- lapply(all_pats, function(pat) {
	genes <- slimp$feature[slimp$pattern == pat]
	enr_regs <- all_regs[strsplit(pat, "-", fixed = TRUE)[[1]] == "1"]

	sub_meta <- trg_meta %>%
				mutate(pseudotime = as.character(pseudotime)) %>%
				filter(region %in% enr_regs)
	knots <- table(sub_meta$pseudotime) %>%
				.[. == length(enr_regs)] %>%
				names()

	screen_g <- sapply(genes, function(gg) {
		expr <- lapply(enr_regs, function(reg) trg_smt[gg, paste0(reg, "|", knots)]
			) %>%
			setNames(., enr_regs) %>%
			as.data.frame() %>%
			as.matrix()

		minmax_fc <- apply(expr, 2, function(x) (max(x) + 0.05)/(min(x) + 0.05))
		if (max(minmax_fc) < 1.5){
			keepg <- TRUE
		} else {
			test_cor <- cor(expr, method = "s")
			test_cor[upper.tri(test_cor, diag = TRUE)] <- NA
			mcor <- mean(test_cor, na.rm = TRUE)
			keepg <- ifelse(mcor >= 0.4, TRUE, FALSE)
		}
		return(keepg)
		})
	return(genes[screen_g])
	}) %>%
	setNames(., all_pats)


fil_pie <- lapply(all_pats, function(pat) {
	sub_pie <- slimp %>%
			filter(feature %in% filter_genes[[pat]] & pattern == pat)
	return(sub_pie)
	}) %>%
	do.call(rbind, .)
saveRDS(fil_pie, file = paste0("./load_files/", "Region_DEGs_pairwise_res.rds"))








