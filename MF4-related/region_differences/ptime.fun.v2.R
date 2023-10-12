##library(URD)
library(circlize)
library(ComplexHeatmap)
library(scatterpie)
library(ggforce)
library(tidyr)

FitImpulseModel <- function(genes, meta, pt_col = "pseudotime", knot_col = "avgcls", data, kvec = c(20, 15), min.slope = 0.1) {
	ptime <- meta[, pt_col]
	knots <- meta[, knot_col]
	expr <- as.matrix(data)[genes, knots]


	#sub_fit <- lapply(genes, function(gg) {
	#	print(paste0("Working on gene: ", which(genes == gg), " / ", length(genes)))
	#	ff <- tryCatch(expr = {
	#				impulseFit(x = ptime, y = expr[gg, ], k = kvec[1], interpolate = 50)
	#			}, 
	#			error = function(cond) {
	#				impulseFit(x = ptime, y = expr[gg, ], k = kvec[2], interpolate = 50)
	#			})
	#	return(ff)
	#}) %>%
	#	setNames(., genes)


	kchoics <- kvec[1]:kvec[2]
	sub_fit <- lapply(genes, function(gg) {
		print(paste0("Working on gene: ", which(genes == gg), " / ", length(genes)))
		modsuccess <- FALSE
		idx <- 0L
		while(modsuccess == FALSE)
			{
				tryCatch(expr = {
						idx <- idx + 1L
						if (gg %in% c("ACTN1", "MAGEF1")){
							print(paste0("Current IDX is ", idx))
						}
			    		ff <- impulseFit(x = ptime, y = expr[gg, ], k = kchoics[idx], interpolate = 100, min.slope = min.slope)
			    		modsuccess <- TRUE
			  		}, error=function(e) {
			  		},
			  		finally={})
			}
		return(ff)
		})


	# Get out onset/offset times  
    timing <- data.frame(
                time.on=unlist(lapply(sub_fit, function(x) if(is.list(x)){return(min(x[['time.on']]))}else{return(x['time.on'])})),
                time.off=unlist(lapply(sub_fit, function(x) if(is.list(x)){return(max(x[['time.off']]))}else{return(x['time.off'])})),
                row.names = genes, 
                gene = genes, 
                stringsAsFactors = FALSE)
    timing[intersect(which(is.na(timing$time.on)), which(is.infinite(timing$time.off))), "time.on"] <- Inf
    return(timing)
}



plot_heatmap.RGCcascade <- function(mat, meta, label_genes, color_breaks = seq(-1.5, 2.5, 0.5), file_name, pdf_height = 5, row_split = NULL, highlight_genes = "HHHHH", fontsize) {
	## Set color breaks
	col_fun1 = colorRamp2(color_breaks, viridis(length(color_breaks)))


	lobe_cols <- c("#B037C4", "#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5") %>% setNames(., c("GE", "FC", "MSC", "OcC", "Insula", "TC"))
	reg_cols <- rep(c("#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5"), c(4, 4, 1, 1, 3)) %>% 
				setNames(., c("MFC", "OFC", "DFC", "VFC", 
					"M1C", "S1C", "IPC", "PCC",
					"V1C",
					"Insula",
					"A1C", "ITC", "STC"))
	cls_cols <- c("#990939", "#e25a9a", "#fc19cf", "#1ac9a9", "#a3c10d", 
	"#efbc88", "#dbb75c", "#6700e5", "#7133aa", "#003563"
	) %>%
			setNames(., c("RGC FABP7 PMP22", 
				"NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", 
				"oRG HOPX TNC", "oRG HOPX APOE", 
				"tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))



	## The columns of the meta data has to be in "factor" format
	meta <- meta[, c("region", "cluster", "pseudotime")]
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
	meta <- meta[order(meta$region, meta$pseudotime), c("region", "cluster")]


	## Plot the column annotation
	column_ha <- HeatmapAnnotation(df = meta, col = col_list, annotation_height = unit(rep(0.01, ncol(meta)), "in"))

	## Prepare row annotation
	htlist <- Heatmap(mat, name = "Expression",  
        cluster_columns = FALSE,
        column_split = meta$region,
        row_split = row_split,
        top_annotation = column_ha, 
        show_column_dend = FALSE,
        cluster_column_slices = FALSE,
        cluster_row_slices = FALSE,
        column_title_gp = gpar(fontsize = 12),
        row_title_gp = gpar(fontsize = 0),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun1,
        row_names_gp = gpar(fontsize = 4),
        column_title_rot = 90,
        show_column_names = FALSE,
        show_row_names = TRUE,
        use_raster = TRUE,
        raster_quality = 4, 
        width = unit(4, "in"))

	match_order <- function(all_genes, selected_genes) {
			selected_genes[match(all_genes[all_genes %in% names(selected_genes)], names(selected_genes))]
	}
	if (!is.null(label_genes)){
		new_list <- list(c(label_genes[["FC"]], label_genes[["TC"]]), 
						c(label_genes[["MSC"]], label_genes[["OcC"]]))
		htlist <- rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% names(new_list[[1]])), 
                        labels = match_order(rownames(mat), new_list[[1]]), 
                        side = "left", 
                        labels_gp = gpar(fontsize = fontsize, col = ifelse(rownames(mat)[rownames(mat) %in% names(new_list[[1]])] %in% highlight_genes, "red", "black")), padding = unit(1, "mm"))) +
				htlist +
				rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% names(new_list[[2]])), 
                        labels = match_order(rownames(mat), new_list[[2]]), 
                        side = "right", 
                        labels_gp = gpar(fontsize = fontsize, col = ifelse(rownames(mat)[rownames(mat) %in% names(new_list[[2]])] %in% highlight_genes, "red", "black")), padding = unit(1, "mm")))
	}


	return(htlist)
	##pdf(paste0(outputdir, file_name, "_heatmap.pdf"), width = 10, height = pdf_height)
	##draw(htlist)
	##dev.off() 
}



plot_heatmap.RGCcascade.regionshare <- function(mat, meta, label_genes, color_breaks = seq(-1.5, 2.5, 0.5), file_name, pdf_height = 5, row_split = NULL, highlight_genes = "HHHHH", fontsize) {
	## Set color breaks
	col_fun1 = colorRamp2(color_breaks, viridis(length(color_breaks)))


	lobe_cols <- c("#B037C4", "#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5") %>% setNames(., c("GE", "FC", "MSC", "OcC", "Insula", "TC"))
	reg_cols <- rep(c("#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5"), c(4, 4, 1, 1, 3)) %>% 
				setNames(., c("MFC", "OFC", "DFC", "VFC", 
					"M1C", "S1C", "IPC", "PCC",
					"V1C",
					"Insula",
					"A1C", "ITC", "STC"))
	cls_cols <- c("#990939", "#e25a9a", "#fc19cf", "#1ac9a9", "#a3c10d", 
	"#efbc88", "#dbb75c", "#6700e5", "#7133aa", "#003563"
	) %>%
			setNames(., c("RGC FABP7 PMP22", 
				"NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", 
				"oRG HOPX TNC", "oRG HOPX APOE", 
				"tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))



	## The columns of the meta data has to be in "factor" format
	meta <- meta[, c("region", "cluster", "pseudotime")]
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
	meta <- meta[order(meta$region, meta$pseudotime), c("region", "cluster")]


	## Plot the column annotation
	column_ha <- HeatmapAnnotation(df = meta, col = col_list, annotation_height = unit(rep(0.01, ncol(meta)), "in"))

	## Prepare row annotation
	htlist <- Heatmap(mat, name = "Expression",  
        cluster_columns = FALSE,
        column_split = meta$region,
        row_split = row_split,
        top_annotation = column_ha, 
        show_column_dend = FALSE,
        cluster_column_slices = FALSE,
        cluster_row_slices = FALSE,
        column_title_gp = gpar(fontsize = 12),
        row_title_gp = gpar(fontsize = 0),
        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun1,
        row_names_gp = gpar(fontsize = 1),
        column_title_rot = 90,
        show_column_names = FALSE,
        show_row_names = TRUE,
        use_raster = TRUE,
        raster_quality = 4, 
        width = unit(4, "in"))

	match_order <- function(all_genes, selected_genes) {
			selected_genes[match(all_genes[all_genes %in% names(selected_genes)], names(selected_genes))]
	}
	if (!is.null(label_genes)){
		new_list <- label_genes
		htlist <- rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% names(new_list[[1]])), 
                        labels = match_order(rownames(mat), new_list[[1]]), 
                        side = "left", 
                        labels_gp = gpar(fontsize = fontsize, col = ifelse(rownames(mat)[rownames(mat) %in% names(new_list[[1]])] %in% highlight_genes, "red", "black")), padding = unit(1, "mm"))) +
				htlist +
				rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% names(new_list[[2]])), 
                        labels = match_order(rownames(mat), new_list[[2]]), 
                        side = "right", 
                        labels_gp = gpar(fontsize = fontsize, col = ifelse(rownames(mat)[rownames(mat) %in% names(new_list[[2]])] %in% highlight_genes, "red", "black")), padding = unit(1, "mm")))
	}


	return(htlist)
	##pdf(paste0(outputdir, file_name, "_heatmap.pdf"), width = 10, height = pdf_height)
	##draw(htlist)
	##dev.off() 
}



plot_region_line <- function(meta, data, features, reg_col = "lobe", pseudotime_col = "pseudotime", file_name, output_dir = outputdir, reg_order = c("FC", "MSC", "TC", "OcC"), return.plot = FALSE, plot.scale = 1, ncol = NULL){
	reg_cols <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4") %>% 
				setNames(., c("FC", "MSC", "TC", "OcC", "Insula", "GE"))
	cls_cols <- c("#990939", "#e25a9a", "#fc19cf", "#1ac9a9", "#a3c10d", 
	"#efbc88", "#dbb75c", "#6700e5", "#7133aa", "#003563"
	) %>%
			setNames(., c("RGC FABP7 PMP22", 
				"NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", 
				"oRG HOPX TNC", "oRG HOPX APOE", 
				"tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))


	plot_data <- data[features, ,drop = FALSE] %>% 
					t() %>%
					cbind(meta[colnames(data), c(pseudotime_col, reg_col, "cluster")], .) %>%
					tidyr::gather(., "gene", "exp", features) %>%
					filter(!!sym(reg_col) %in% reg_order) %>%
					mutate(!!sym(reg_col) := factor(!!sym(reg_col), levels = reg_order))
	cur_regs <- unique(as.character(meta[, reg_col]))


	## Add 
    p <- ggplot(plot_data, aes_string(x = pseudotime_col, y = "exp", color = reg_col)) +
                geom_smooth(size = 1.5, span = 0.75, method = "loess", se = TRUE) +
                geom_tile(mapping = aes_string(fill = "cluster", y = max(plot_data$exp)), interpolate = TRUE, color = NA) +
                scale_color_manual(values = reg_cols) +
                scale_fill_manual(values = cls_cols) +
                theme_classic() +
                facet_grid(cols = vars(gene)) + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                		axis.text = element_text(size = rel(0.8)), 
                		axis.title.x = element_blank(), 
                		axis.title.y = element_text(size = rel(0.9)), 
                		plot.title = element_blank(), 
                		axis.line = element_line(size = 0.2),
                		axis.ticks = element_line(size = 0.2),
                		legend.position = "bottom",
                		strip.background = element_rect(size = 0.2), strip.text = element_text(face = "bold"))

    if (return.plot){
        return(p)
    } else {
    	#if (is.null(ncol)){
    	#	ncol <- ifelse(length(features) >= 2, 2, 1)
    	#} 
    	ncol = length(features)
        pdf_width <- plot.scale * 3 * ncol
		pdf_height <- plot.scale * 4 + ceiling(length(features)/ncol)

		pdf(paste0(output_dir, file_name, "_Line.pdf"), width = pdf_width + 1, height = 4.5)
		print(p)
		dev.off()
    }
}


geom_scatterpie_new <- function(mapping = NULL, data, cols, legend_name = "species", scale.expression = TRUE, ...) {
      names(mapping)[match(c("x", "y"), names(mapping))] <- c("x0", "y0")
      yvar <- get_aes_var(mapping, "y0")
      xvar <- get_aes_var(mapping, "x0")
      
      
      col_use <- names(cols)
      df <- data[rowSums(data[, col_use]) > 0, ] %>%
               gather(., "type", "value", !!enquo(col_use)) %>% 
               group_by(!!sym(xvar), !!sym(yvar)) %>%
               mutate(scaleexp = value/mean(value)) %>%
               mutate(scaleexp = MinMax(scaleexp, min = 0.5, max = 2)) %>% 
               ungroup()


      nbks <- 1000
      df$border_color <- cols[as.character(df$type)]      
      ## set the border & fill color:
      if (scale.expression){
         df$fill_color <- mapply(FUN = function(color, value) {
                           return(colorRampPalette(colors = c("grey95", color))(nbks)[value])
                  }, color = df$border_color, value = as.numeric(cut(df$scaleexp, nbks)))
         df$border_color[df$scaleexp < 1.5] <- NA
      } else{
         df$fill_color <- df$border_color
         df$border_color <- NA
      }

      
      ## Factorize the type column
      df$type <- factor(df$type, levels = col_use)
      names(df)[which(names(df) == "type")] = legend_name
      df$r0 <- 0
      return(geom_arc_bar(mapping, data = df, stat = "pie", inherit.aes = FALSE, ...))
}
environment(geom_scatterpie_new) <- environment(geom_scatterpie)



PlotWeightedScatterPie <- function(pie.data, x.col = "cluster", y.col = "gene", r.col = "radius", cls_use = c("vRG.early", "vRG.late", "oRG"), rsf = 2, scale.expression = TRUE, fc_limit = c(-1, 2.5)){ 
	## rsf:radius-scale-factor
	## set the colors 
	cls_cols <- c("#fc19cf", "#a3c10d", "#dbb75c") %>% setNames(., c("vRG.early", "vRG.late", "oRG"))
	pie.plot <- pie.data
	cls_use <- intersect(names(cls_cols), cls_use)
	cls_cols <- cls_cols[cls_use]


	## legend data
	##lg.ra <- pie.plot[, r.col]
	##lg.ra[lg.ra == 0] <- 0.00001

	pie.plot[, x.col] <- MinMax(pie.plot[, x.col], min = fc_limit[1], max = fc_limit[2])
	pie.plot[, y.col] <- MinMax(pie.plot[, y.col], min = fc_limit[1], max = fc_limit[2])
	pie.plot$latetype <- factor(as.character(pie.plot$latetype), levels = c("logFC|oRG", "logFC|vRG late"))
	##pie.plot$label <- ifelse(pie.plot[, x.col] > 1 & pie.plot[, y.col] > 1 | , pie.plot$gene, NA)
	pie.plot$label <- ifelse(rowSums(pie.plot[, names(cls_cols)]) == length(cls_cols), pie.plot$gene, NA)

	## Plots
	p <- ggplot(data = pie.plot) + 
			geom_scatterpie_new(mapping = aes_string(x = x.col, y = y.col, r0 = "r0", r = r.col, amount = "value", fill = "fill_color", color = "border_color"), data = pie.plot, cols = cls_cols, scale.expression = scale.expression) + 
			ggrepel::geom_text_repel(mapping = aes_string(x = x.col, y = y.col, label = "label"), color = "black", seed = 42, size = 4) +
			theme_cowplot() + 
			scale_fill_identity() + 
			scale_color_identity() + 
			##RotatedAxis() + 
			facet_grid(rows = vars(latetype)) +
			coord_fixed(ratio = 1, xlim = fc_limit, ylim = fc_limit) +
			theme(panel.grid.major = element_line(colour="lightgrey", size = 0.2), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text = element_text(size = 11)) ##+
      ##geom_scatterpie_legend(lg.ra, x= +1, y=ceiling(length(feature.order)/4), n = 4, labeller = function(x) rsf * x)
   return(p)
}


