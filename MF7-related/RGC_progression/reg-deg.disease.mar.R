source("../scripts/nhpf.fun.R")


seu_use <- readRDS(file = paste0("./load_files/", "Reg-DEG_data_seu.rds"))
exp_list <- readRDS(file = paste0("./load_files/", "Expgenes_res.rds"))



## Identify cell=type markers
allgps <- levels(as.factor(seu_use@meta.data$subtype2))
alllbs <- levels(as.factor(seu_use@meta.data$lobe))
allres <- lapply(alllbs, function(lb) {
	lb_seu <- subset(seu_use, lobe %in% lb)
	gp_use <- levels(as.factor(as.character(lb_seu@meta.data$subtype2))) %>%
				intersect(allgps, .)


	lb_res <- lapply(gp_use, function(gp) {
		Idents(lb_seu) <- "subtype2"
		res <- FindMarkers(lb_seu, features = exp_list[[gp]], ident.1 = gp, max.cells.per.ident = 1000, only.pos = TRUE, logfc.threshold = 0.1, min.pct = 0.1) %>%
				rownames_to_column("gene") %>%
				mutate(cluster = gp) %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(region = lb)
		return(res)
		}) %>%
		do.call(rbind, .)

	return(lb_res)
	}) %>%
	do.call(rbind, .)
saveRDS(allres, file = paste0(inputdir, "Marker_res.rds"))















