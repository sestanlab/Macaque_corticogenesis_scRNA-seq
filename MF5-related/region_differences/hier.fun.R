DoHierCluster <- function(data, method) {
	if (method %in% "cor"){
		hc <- as.dist(1 - cor(data, method = "p")) %>%
				hclust(., method = "ward.D2")
	} else if (method %in% "dist") {
		hc <- dist(scale(t(data))) %>%
				hclust(., method = "ward.D2")
	} else {
		stop("Please enter the right method: 'cor' or 'dist'")
	}
	return(hc)
}



plot_circlize_dendrogram <- function(data, dist.method = c("cor", "dist")) {
	colnames(data) <- strsplit(colnames(data), "|", fixed = TRUE) %>%
					sapply(., function(x) x[2])

	hc <- DoHierCluster(data = data, method = dist.method)
	dend <- hc %>% as.dendrogram
	dend <- dend %>%
		set("labels_colors", value = reg_cols[labels(dend)]) %>%
		branches_attr_by_labels(., labels = c("OFC", "DFC", "MFC", "VFC"), TF_value = c("#FF420E", Inf), type = "all") %>%
		branches_attr_by_labels(., labels = c("M1C", "S1C", "IPC", "PCC"), TF_value = c("#FFBB00", Inf), type = "all") %>%
		branches_attr_by_labels(., labels = c("A1C", "ITC", "STC"), TF_value = c("#4CB5F5", Inf), type = "all") %>%
		branches_attr_by_labels(., labels = c("V1C"), TF_value = c("#89DA59", Inf), type = "all")

	circlize_dendrogram(dend,
	                    labels_track_height = 0.2,
	                    dend_track_height = 0.6)
}



percent_match <- function(x, n = 100) {
  return(Re(x) / (n - Im(x)))
}


Bootstrap_cocluster <- function(data, hvg.list, dist.method = c("cor", "dist")[1], reg_ord) {
	colnames(data) <- strsplit(colnames(data), "|", fixed = TRUE) %>%
					sapply(., function(x) x[2])

	empt_id <- rep(NA, length(reg_ord)) %>%
				setNames(., reg_ord)

	mtch_mat <- lapply(hvg.list, function(hvguse) {
		data_sub <- data[hvguse, ,drop = FALSE]
		id <- DoHierCluster(data = data_sub, method = dist.method) %>%
				cutree(., k = 3)
		id_use <- empt_id
		id_use[names(id)] <- id

		mtch <- outer(id_use, id_use, "==")
		mtch[is.na(mtch)] <- 1i ## non records drop as imaginary, matches as 1, non-matches as 0
		return(mtch)
		}) %>%
		Reduce("+", .)
	
	p_match <- percent_match(x = mtch_mat, n = length(hvg.list))
	return(p_match)
}


Plot_Cocluster_Heatmap <- function(mat, reg_ord, rm_na = TRUE, title = "") {
	if (rm_na){
	 	na_cols <- rownames(mat)[apply(mat, 1, function(x) all(is.na(x)))]
	 	reg_ord <- setdiff(reg_ord, na_cols)
	 	mat <- mat[reg_ord, reg_ord]
	}


	col_bks <- seq(0, 1, 0.1)
	col_fun <- colorRamp2(col_bks, viridis(length(col_bks)))
	meta <- data.frame(cell = reg_ord, 
							stringsAsFactors = FALSE) %>%
						mutate(lobe = reg2lb[cell]) %>%
						column_to_rownames("cell")
	column_ha <- HeatmapAnnotation(df = meta, col = list(lobe = lb_cols), name = "region", 
									annotation_height = unit(0.15, "in"))
	right_ha <- HeatmapAnnotation(df = meta, col = list(lobe = lb_cols), which = "row", name = "region",
									annotation_width = unit(0.15, "in"))
	htlist <- Heatmap(mat, column_title = title, name = "Expression",  
	        cluster_columns = FALSE,
	        na_col = "lightgrey",
	        top_annotation = column_ha, 
	        show_column_dend = FALSE,
	        cluster_column_slices = FALSE,
	        cluster_row_slices = FALSE,
	        column_gap = unit(0.5, "mm"),
	        cluster_rows = FALSE,
	        show_row_dend = FALSE,
	        col = col_fun,
	        row_names_gp = gpar(fontsize = 1),
	        show_column_names = TRUE,
	        show_row_names = TRUE,
	        use_raster = FALSE,
	        width = unit(2, "in"),
	        height = unit(2, "in")) +
			right_ha
	return(htlist)
}



Permutate_cocluster <- function(data.list, dist.method = c("cor", "dist")[1], reg_ord) {
	data.list <- lapply(data.list, function(data) {
		colnames(data) <- strsplit(colnames(data), "|", fixed = TRUE) %>%
					sapply(., function(x) x[2])
		return(data)
		})

	empt_id <- rep(NA, length(reg_ord)) %>%
				setNames(., reg_ord)

	mtch_mat <- lapply(data.list, function(data) {
		id <- DoHierCluster(data = data, method = dist.method) %>%
				cutree(., k = 3)
		id_use <- empt_id
		id_use[names(id)] <- id

		mtch <- outer(id_use, id_use, "==")
		mtch[is.na(mtch)] <- 1i ## non records drop as imaginary, matches as 1, non-matches as 0
		return(mtch)
		}) %>%
		Reduce("+", .)
	
	p_match <- percent_match(x = mtch_mat, n = length(data.list))
	return(p_match)
}










