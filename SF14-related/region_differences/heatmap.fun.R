plot_heatmap.end <- function(mat, label_genes, color_breaks = seq(-1.5, 2.5, 0.5), module_labs = NULL, pdf_height = 5, row_dend = FALSE) {
	## Set color breaks
	col_fun1 = colorRamp2(color_breaks, viridis(length(color_breaks)))


	lobe_cols <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4") %>% setNames(., c("FC", "MSC", "TC", "OC", "Insula", "GE"))
	reg_cols <- rep(c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787"), c(4, 4, 3, 1, 1)) %>% 
				setNames(., c("DFC", "MFC", "OFC", "VFC", 
					"M1C", "S1C", "IPC", "PCC",
					"A1C", "ITC", "STC", 
					"V1C", 
					"Insula"))
	cls_cols <- setNames(c("#d8b365", "#5ab4ac"),
                    c("MGE-InN", "CGE-InN"))


	## The columns of the meta data has to be in "factor" format
	meta <- data.frame(cell = colnames(mat), 
						stringsAsFactors = FALSE) %>%
					mutate(region = sapply(strsplit(cell, "|", fixed = TRUE), function(mm) mm[1])) %>%
					mutate(cluster = sapply(strsplit(cell, "|", fixed = TRUE), function(mm) mm[2])) %>%
					column_to_rownames("cell")
	if ("FC" %in% meta$region){
		reg_col_use <- lobe_cols
	} else {
		reg_col_use <- reg_cols
	}

	reg_ord <- intersect(names(reg_col_use), meta$region)
	cls_ord <- intersect(names(cls_cols), meta$cluster)
	col_list <- list(region = reg_col_use[reg_ord], 
					cluster = cls_cols[cls_ord])
	meta$region <- factor(meta$region, levels = reg_ord)
	meta$cluster <- factor(meta$cluster, levels = cls_ord)


	## Plot the column annotation
	#column_ha <- HeatmapAnnotation(df = meta, col = col_list, annotation_height = unit(rep(0.01, ncol(meta)), "in"))
	column_ha <- HeatmapAnnotation(link = anno_mark(at = seq(from = 10, to = nrow(meta), by = 20), 
						labels = unique(meta$region),
						which = "column", 
						side = "top",
						labels_gp = gpar(fontsize = 12, col = "black"), 
						padding = unit(1, "mm")),
								df = meta, col = col_list,
								name = "ii")
	##column_ha <- HeatmapAnnotation(df = meta, col = col_list, name = "region")

	
	## Prepare row annotation
	gord <- intersect(rownames(mat), label_genes)
	if (!is.null(module_labs)){
		gmeta <- data.frame(row.names = rownames(mat), module = module_labs, stringsAsFactors = FALSE)
		module_labs <- setNames(module_labs, rownames(mat))
		morder <- unique(module_labs)
		mcols <- setNames(gg_color_hue(length(morder)), morder)
		lb_mod <- module_labs[rownames(mat)[which(rownames(mat) %in% gord)]]
		text_cols <- mcols[lb_mod]


		right_ha2 <- HeatmapAnnotation(df = gmeta,
								col = list(module = mcols),
								which = "row",
								name = "module",
								annotation_width = unit(0.1, "in"))
	} else {
		text_cols <- "black"
	}


	right_ha1 <- HeatmapAnnotation(link = anno_mark(at = which(rownames(mat) %in% gord), 
				                    labels = rownames(mat)[which(rownames(mat) %in% gord)], 
				                    side = "right",
				                    labels_gp = gpar(fontsize = 12, col = text_cols), padding = unit(1, "mm")),
								which = "row",
								name = "ii")
	htlist <- Heatmap(mat, name = "Expression",  
        cluster_columns = FALSE,
        top_annotation = column_ha, 
        show_column_dend = FALSE,
        cluster_column_slices = FALSE,
        cluster_row_slices = FALSE,
        column_title_gp = gpar(fontsize = 0),
        row_title_gp = gpar(fontsize = 0),
        column_gap = unit(0.5, "mm"),
        show_row_dend = TRUE,
        col = col_fun1,
        row_names_gp = gpar(fontsize = 1),
        column_title_rot = 90,
        show_column_names = FALSE,
        show_row_names = TRUE,
        use_raster = TRUE,
        raster_quality = 4, 
        cluster_rows = row_dend,
        width = unit(4, "in"))
		
	if (!is.null(label_genes)){
		htlist <- htlist + 
					right_ha1
	}

	return(htlist)
	##pdf(paste0(outputdir, file_name, "_heatmap.pdf"), width = 8, height = pdf_height)
	##draw(htlist)
	##dev.off() 
}


