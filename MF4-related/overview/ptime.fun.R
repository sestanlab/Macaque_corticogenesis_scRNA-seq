library(URD)
library(WGCNA)
library(circlize)
library(ComplexHeatmap)
library(GetoptLong)
library(ggpubr)
PulsefitProcess <- function(data_use, genes, ptime, k = 50, interpolate = 50, pulse.only = TRUE, nCores = 8) {
	# Do impulse model fitting for all genes
	bin.size <- 10
	bin.ind <- ceiling(c(1:length(genes))/bin.size)
	max.bin <- max(bin.ind)
	

	cl = makeCluster(nCores, outfile="")
	doParallel::registerDoParallel(cl);
	tem_idx <- NULL
	res_list <- foreach(tem_idx = 1:max.bin, .combine = c, .packages = c("URD"), .export =c("PulsefitProcess")) %dopar% {
		curgenes <- genes[bin.ind == tem_idx]
		print(paste0(Sys.time(), ":    ", tem_idx))
		flist <- lapply(curgenes, function(g) {
				md <- impulseFit(x = ptime[colnames(data_use)], y = as.numeric(data_use[g, ]), k = k, interpolate = interpolate, pulse.only = pulse.only, min.slope = 0.1)
				return(md)
				})
		names(flist) <- curgenes
		return(flist)
	}
	resg_list <- res_list[genes]
	stopCluster(cl)


	#impulse.fits <- lapply(genes, function(g) {
	#	if (verbose.genes) {
	#		message(paste0(Sys.time(), ":    ", g))
	#	}
	#	md <- impulseFit(x = ptime[colnames(data_use)], 
	#		y = as.numeric(data_use[g, ]), k = k, interpolate = interpolate, pulse.only = pulse.only, min.slope = 0.1#)
	#	return(md)
	#	}) %>%
	#		setNames(., genes)



	# Get out onset/offset times
	timing <- data.frame(
		time.on=unlist(lapply(resg_list, function(x) if(is.list(x)){return(min(x[['time.on']]))}else{return(x['time.on'])})),
		time.off=unlist(lapply(resg_list, function(x) if(is.list(x)){return(max(x[['time.off']]))}else{return(x['time.off'])})),
		row.names=genes, 
		stringsAsFactors=FALSE
	)


	return(list(
		timing = timing,
		impulse.fits = resg_list))
}



###--------------------------------------------------------------------------------------

## Bin average

###--------------------------------------------------------------------------------------
SmoothBinAvg_combine <- function(data, ptime, predict_df, loess_span = 0.6){
	## Examine the column names
	if (!"pseudotime" %in% colnames(predict_df)){
		stop(paste0("Make sure the predict_df contains a pseudotime column"))
	}


	## Separating genes to bins to boost loess fit speed
    print(Sys.time())
    genebin.size <- 100
    bin.ind <- ceiling(c(1:nrow(data))/genebin.size)
    max.bin <- max(bin.ind)


    cbn_fit <- lapply(1:max.bin, function(idx) {
        print(paste0("Smoothing for gene bin: ", idx))
        genes <- rownames(data)[bin.ind == idx]
        newdata <- data.frame(t(data[genes, ]), check.names = TRUE) 
        new_genes <- colnames(newdata)
        newdata$pseudotime <- ptime[rownames(newdata)]
        fit_data <- lapply(new_genes, function(gene) {
            fit_formula <- as.formula(paste0(gene, " ~ ", "pseudotime"))
            new_fit <- loess(fit_formula, newdata, na.action = "na.omit", span = loess_span)
            fitval <- predict(new_fit, predict_df, se = FALSE) %>% setNames(., rownames(predict_df))
            fitval[fitval <0 ] <- 0 ##some fitted values below zero are not meaningful
            return(fitval)
            }) %>%
            setNames(., genes) %>%
            as.data.frame(., check.names = FALSE) %>%
            as.matrix() %>% 
            t()
        return(fit_data)
        }) %>%
        do.call(rbind, .) %>%
        as.matrix()

    cbn_fit <- cbn_fit[rownames(data), , drop = FALSE]
    return(cbn_fit)
}


###--------------------------------------------------------------------------------------

## Find Markers

###--------------------------------------------------------------------------------------


FindCommonMars <- function(obj.list, ident.1, ident.2, regions, ratio_thre = 1.1, p_thre = 0.01) {
	## Calculate markers
	allres <- lapply(regions, function(reg) {
		xx <- obj.list[[reg]]
		res <- FindMarkers(xx, ident.1 = ident.1, ident.2 = ident.2, only.pos = TRUE, max.cells.per.ident = 1000, min.pct = 0.1, logfc.threshold = 0.25) %>%
				rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(region = reg)
		res
		}) %>%
			do.call(rbind, .)

	message(paste0("Finding markers for cluster: ", ident.1))
	res <- allres %>%
				filter(ratio_fc >= ratio_thre & p_val_adj <= p_thre) %>%
				group_by(gene) %>%
				summarize(nhits = n(), mfc = mean(ratio_fc)) %>%
				ungroup() %>%
				filter(nhits >= ceiling(length(regions) * 0.7)) %>%
				mutate(cluster = ident.1) %>%
				arrange(desc(nhits), desc(mfc))
	return(res)
}



FindSpecMars <- function(obj.list, ident, regions, ratio_thre = 1.1, p_thre = 0.01, max.cells.per.ident = 1000) {
	## Calculate markers
	object <- obj.list[[ident]]
	message(paste0("Finding markers for cluster: ", ident))

	allres <- lapply(regions, function(reg) {
		bg.regs <- setdiff(regions, reg)
		reg_res <- lapply(bg.regs, function(bg1) {
			res <- FindMarkers(object, ident.1 = reg, ident.2 = bg1, only.pos = TRUE, max.cells.per.ident = max.cells.per.ident, min.pct = 0.1, logfc.threshold = 0.25) %>%
				rownames_to_column("gene") %>%
				mutate(ratio_fc = (pct.1 + 0.01)/(pct.2 + 0.01)) %>%
				mutate(background = bg1)
			res
			}) %>%
			do.call(rbind, .)

		reg_res <- reg_res %>%
					filter(ratio_fc >= ratio_thre & p_val_adj <= p_thre) %>%
					group_by(gene) %>%
					summarize(nhits = n(), mfc = mean(ratio_fc)) %>%
					ungroup() %>%
					filter(nhits == length(bg.regs)) %>%
					mutate(cluster = ident, region = reg) %>%
					arrange(desc(nhits), desc(mfc))
		return(reg_res)
		}) %>%
			do.call(rbind, .)
	return(allres)
}

###--------------------------------------------------------------------------------------

## Prepare Weighted dot plots

###--------------------------------------------------------------------------------------

spdf2arry <- function(df) {
	## Add empty cols for missing pairs (diagonal Human-Human, Chimp-Chimp, etc)
	all_regs <- c("FC", "MSC", "TC", "OcC")
	cur_regs <- strsplit(colnames(df), "|", fixed = TRUE) %>% unlist() %>% unique() %>% intersect(., all_regs)
	all_pairs <- rep(cur_regs, each = length(cur_regs)) %>%
					paste0(., "|", rep(cur_regs, length(cur_regs)))
	mis_pairs <- setdiff(all_pairs, colnames(df))
	if (length(mis_pairs) > 0){
		mis_df <- matrix(0, nrow = nrow(df), ncol = length(mis_pairs), dimnames = list(rownames(df), mis_pairs))
		df <- cbind(as.matrix(df), mis_df)
	}
	df <- df[, all_pairs] %>%
			as.matrix()

	dexinfo <- lapply(1:nrow(df), function(i) {matrix(df[i, ], nrow = length(cur_regs), ncol = length(cur_regs), byrow = TRUE)}) %>%
			do.call(c, .) %>%
			array(., dim = c(length(cur_regs), length(cur_regs), nrow(df)), dimnames = list(cur_regs, cur_regs, rownames(df)))
	return(dexinfo)
}




SummariseArray <- function(ary) {
	genes <- dimnames(ary)[[3]]
	rSUM <- lapply(genes, function(x) rowSums(ary[, , x])) %>%
			setNames(., genes) %>%
			as.data.frame(., check.names = FALSE) %>%
			as.matrix() %>% t()
	cSUM <- lapply(genes, function(x) colSums(ary[, , x])) %>%
			setNames(., genes) %>%
			as.data.frame(., check.names = FALSE) %>%
			as.matrix() %>% t()


	##------------------------------------------------
	## Transfer rSUM & cSUM to genes * 4sp
	## 3 models: model 1:species exclusively enriched; mod 2: enriched in 2 species; mod3 depleted in one species
	## Incase mod3 & mod1 have overlaps, put mod1 later will put mod1 in high priority
	## model 2 will not overlap with 
	modeltest_n3 <- function(rSUM, cSUM, ept_mat, all_regs) {
			nregs <- length(all_regs)
			ept_mat <- matrix(0, nrow = nrow(rSUM), ncol = length(all_regs), dimnames = list(rownames(rSUM), all_regs))

			if (nregs < 3){
				stop("modeltest_n3 require at least 3 regions")
			}
			thre <- nregs - 1 ## If nregs = 2, thre will be 1
			mod3_idx <- which(apply(cSUM, 1, function(x) sum(x == thre) == 1))
			tem_mat <- cSUM[mod3_idx, ]
			tem_mat[tem_mat != thre] <- 1
			tem_mat[tem_mat == thre] <- 0
			ept_mat[mod3_idx, ] <- tem_mat


			mod1_idx <- which(apply(rSUM, 1, function(x) sum(x == thre) == 1))
			tem_mat <- rSUM[mod1_idx, ]
			tem_mat[tem_mat != thre] <- 0
			tem_mat[tem_mat == thre] <- 1
			ept_mat[mod1_idx, ] <- tem_mat


			thre2 <- nregs - 2
			mod2_idx <- which(apply(cSUM, 1, function(x) sum(x == 2) == thre2) & 
								apply(rSUM, 1, function(x) sum(x == thre2) == 2))
			tem_mat <- rSUM[mod2_idx, ]
			tem_mat[tem_mat != thre2] <- 0
			tem_mat[tem_mat == thre2] <- 1
			ept_mat[mod2_idx, ] <- tem_mat
			return(ept_mat)
	}

	modeltest_n2 <- function(rSUM, cSUM, ept_mat, all_regs) {
			nregs <- length(all_regs)
			ept_mat <- matrix(0, nrow = nrow(rSUM), ncol = length(all_regs), dimnames = list(rownames(rSUM), all_regs))

			if (nregs != 2){
				stop("modeltest_n2 require 2 regions")
			}
			thre <- 1
			mod1_idx <- which(apply(rSUM, 1, function(x) sum(x == thre) == 1))
			tem_mat <- rSUM[mod1_idx, ]
			ept_mat[mod1_idx, ] <- tem_mat
			return(ept_mat)
	}

	nregs <- dim(ary)[[1]]
	if (nregs > 2){
		ept_mat <- modeltest_n3(rSUM = rSUM, cSUM = cSUM, ept_mat = ept_mat, all_regs = dimnames(ary)[[1]])
	} else if (nregs == 2){
		ept_mat <- modeltest_n2(rSUM = rSUM, cSUM = cSUM, ept_mat = ept_mat, all_regs = dimnames(ary)[[1]])
	}
	return(ept_mat)
}


library(scatterpie)
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



PlotScatterPie <- function(pie.data, group.col = "cluster", feature.col = "gene", r.col = "radius", split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), group.order = NULL, feature.order = NULL, rsf = 2, scale.expression = TRUE, bg.col = "lightgrey", x_scale = 1, y_scale = 1) { ## rsf:radius-scale-factor
   ## set the colors 
   sp.colors <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4", "#D3D3D3", "#FFFFFF") %>% 
				setNames(., c("FC", "MSC", "TC", "OcC", "Insula", "GE", "Shared", "empty")) %>%   
            .[split.order]


   ## Set the order of x & y Axis
   if (is.null(group.order)){
      group.order <- levels(as.factor(pie.data[, group.col]))
   }
   if (is.null(feature.order)){
      feature.order <- levels(as.factor(pie.data[, feature.col]))
   }

   pie.plot <- pie.data %>%
               mutate(!!group.col := as.character(!!sym(group.col))) %>%
               mutate(!!feature.col := as.character(!!sym(feature.col))) %>%
               mutate(!!group.col := as.numeric(factor(!!sym(group.col), levels = group.order))) %>%
               mutate(!!feature.col := as.numeric(factor(!!sym(feature.col), levels = rev(feature.order))))


   pie.plot[, group.col] <- pie.plot[, group.col] * x_scale
   pie.plot[, feature.col] <- pie.plot[, feature.col] * y_scale


   ## legend data
   lg.ra <- pie.plot[, r.col]
   lg.ra[lg.ra == 0] <- 0.00001




   ## Plots
   p <- ggplot(data = pie.plot) + 
      geom_scatterpie_new(mapping = aes_string(x = group.col, y = feature.col, r0 = "r0", r = r.col, amount = "value", fill = "fill_color", color = "border_color"), data = pie.plot, cols = sp.colors, scale.expression = scale.expression) + 
      theme_cowplot() + 
      ##coord_fixed(xlim = c(0, length(group.order) + 3), ylim = c(0, length(feature.order) + 1), expand = FALSE) + 
      coord_fixed(xlim = c(-1, length(group.order)*x_scale + 3), ylim = c(0, length(feature.order)*y_scale + 1), expand = FALSE) + 
      scale_x_continuous(breaks = (1:length(group.order)) * x_scale, labels = group.order) +
      ##scale_y_continuous(breaks = 1:length(feature.order), labels = rev(feature.order)) +
      scale_y_continuous(breaks = (1:length(feature.order)) * y_scale, labels = rev(feature.order)) +
      ##scale_fill_manual(values = sp.colors)+
      scale_fill_identity() + 
      scale_color_identity() + 
      RotatedAxis() + 
      theme(panel.grid.major = element_line(colour=bg.col, size = 0.2), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text = element_text(size = 11))+
      geom_scatterpie_legend(lg.ra, x=length(group.order)+1, y=ceiling(length(feature.order)/4), n = 4, labeller = function(x) rsf * x)
   return(p)
}





PlotScatterPieContinuous <- function(pie.data, ptime.col = "cluster", feature.col = "gene", r.col = "radius", split.order = c("FC", "MSC", "TC", "OcC"), feature.order = NULL, rsf = 2, scale.expression = TRUE, bg.col = "lightgrey", x_scale = 1, y_scale = 1) { ## rsf:radius-scale-factor
   ## set the colors 
   sp.colors <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4", "#D3D3D3", "#FFFFFF") %>% 
				setNames(., c("FC", "MSC", "TC", "OcC", "Insula", "GE", "Shared", "empty")) %>%   
            .[split.order]


   ## Set the order of x & y Axis
   if (is.null(feature.order)){
      feature.order <- levels(as.factor(pie.data[, feature.col]))
   }

   pie.plot <- pie.data %>%
               mutate(!!ptime.col := as.numeric(!!sym(ptime.col))) %>%
               mutate(!!feature.col := as.character(!!sym(feature.col))) %>%
               mutate(!!feature.col := as.numeric(factor(!!sym(feature.col), levels = rev(feature.order))))

   pie.plot[, ptime.col] <- pie.plot[, ptime.col] * x_scale
   pie.plot[, feature.col] <- pie.plot[, feature.col] * y_scale

   ## legend data
   lg.ra <- pie.plot[, r.col]
   lg.ra[lg.ra == 0] <- 0.00001



   ## Plots
   p <- ggplot(data = pie.plot) + 
   		##geom_tile(data = pdata, mapping = aes_string(x = "pseudotime", y = "gene", fill = "gradient"), color = NA) +
      geom_scatterpie_new(mapping = aes_string(x = ptime.col, y = feature.col, r0 = "r0", r = r.col, amount = "value", fill = "fill_color", color = "border_color"), data = pie.plot, cols = sp.colors, scale.expression = scale.expression) + 
      theme_cowplot() + 
      coord_fixed(xlim = c(-1, ceiling(max(pie.data[, ptime.col]))*x_scale + 3), ylim = c(0, length(feature.order)*y_scale + 1), expand = FALSE) + 
      ##scale_x_continuous(breaks = seq(0, ceiling(max(pie.data[, ptime.col]))*x_scale + 3, 5), labels = ) +
      scale_y_continuous(breaks = (1:length(feature.order)) * y_scale, labels = rev(feature.order)) +
      scale_fill_identity() + 
      scale_color_identity() + 
      RotatedAxis() + 
      theme(panel.grid.major.y = element_line(colour=bg.col, size = 0.2), 
      		panel.grid.major.x = element_blank(), 
      		axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text = element_text(size = rel(0.6)))+
      geom_scatterpie_legend(lg.ra, x=max(pie.data[, ptime.col])+1, y=ceiling(length(feature.order)/4), n = 4, labeller = function(x) rsf * x)
      print(max(pie.data[, ptime.col])+1)
   return(p)
}








plot_region_line <- function(meta, data, features, reg_col = "region", pseudotime_col = "pseudotime", file_name, output_dir = outputdir, reg_order = c("FC", "MSC", "TC", "OcC"), return.plot = FALSE, plot.scale = 1){
	reg_cols <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4") %>% 
				setNames(., c("FC", "MSC", "TC", "OcC", "Insula", "GE"))
	lin_cols <- c("#1ac9a9", "#7133aa", "#dbb75c") %>%
				setNames(., c("shared", "tRG", "oRG"))


	plot_data <- data[features, ,drop = FALSE] %>% 
					t() %>%
					cbind(meta[colnames(data), c(pseudotime_col, reg_col, "lineage")], .) %>%
					filter(!!sym(reg_col) %in% reg_order) %>%
					mutate(!!sym(reg_col) := factor(!!sym(reg_col), levels = reg_order))
	cur_regs <- unique(as.character(meta[, reg_col]))

    plist <- lapply(features, function(gene){
        p <- ggplot(plot_data, aes_string(x = pseudotime_col, y = gene, color = "lineage")) +
                ##geom_line(size = 1.5) + 
                geom_smooth(size = 1.5, span = 0.25, method = "loess", se = FALSE) +
                scale_color_manual(values = lin_cols) +
                theme_bw() +
                scale_x_continuous(limits = c(0, max(plot_data[, pseudotime_col]))) +
                labs(y = gene) +
                facet_wrap(facets = vars(!!sym(reg_col)), ncol = length(cur_regs)) + 
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                		axis.text = element_text(size = rel(0.8)), 
                		axis.title.x = element_blank(), 
                		axis.title.y = element_text(size = rel(0.9)), 
                		plot.title = element_blank(), 
                		axis.line = element_line(size = 0.2),
                		axis.ticks = element_line(size = 0.2),
                		legend.position = "right", 
                		strip.background = element_blank(), strip.text = element_text(face = "bold"))
        p
        })

    if (return.plot){
        return(plist)
    } else {
        pdf_width <- plot.scale * 4 * 2
		pdf_height <- plot.scale * length(plist) * 2
		plot <- patchwork::wrap_plots(plist, nrow = length(plist), ncol = 1, guides = "collect") & theme(legend.position = "right")
		pdf(paste0(output_dir, file_name, "_Line.pdf"), width = pdf_width + 1, height = pdf_height)
		print(plot)
		dev.off()
    }
}






plot_region_heatmap <- function(meta, data, reg_col = "region", group_col = "cluster", pseudotime_col = "pseudotime", anchor_reg = "FC", file_name, output_dir = outputdir, pdf_width = 12, show_rownames = TRUE, width = 14, height = 10, label_genes = NULL, font_scale = 1, row_meta = NULL, reg_order = c("FC", "MSC", "TC", "OcC"), return.plot = FALSE) {
	
	meta <- meta[meta[, reg_col] %in% reg_order, ]
	## Get the full bins based on the anchor regions
	meta <- lapply(reg_order, function(reg) {
		ref <- meta %>%
				filter(region == anchor_reg)
		ref[, reg_col] <- reg
		rownames(ref) <- paste0(reg, "|", ref[, pseudotime_col])
		ref
		}) %>%
			do.call(rbind, .)

	meta[, reg_col] <- factor(as.character(meta[, reg_col]), levels = reg_order) 
	order_knots <- order(meta$region, meta$pseudotime)
	meta <- meta[order_knots, ]


	## Get the full bin data
	mis_bins <- setdiff(rownames(meta), colnames(data))
	mis_data <- matrix(NA, nrow = nrow(data), ncol = length(mis_bins), dimnames = list(rownames(data), mis_bins))
	data <- cbind(data, mis_data) %>%
				.[, rownames(meta)]


    ## Do the plots
    plot_meta <- meta[, c(group_col, reg_col)]
    regcols <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4") %>% 
				setNames(., c("FC", "MSC", "TC", "OcC", "Insula", "GE"))
	cls_cols <- c("#990939", "#e25a9a", "#fc19cf", "#1ac9a9", "#a3c10d", 
	"#efbc88", "#dbb75c", "#6700e5", "#7133aa", "#003563"
	) %>%
			setNames(., c("RGC FABP7 PMP22", 
				"NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", 
				"oRG HOPX TNC", "oRG HOPX APOE", 
				"tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))

    anno_cols <- list(regcols[plot_meta[, reg_col] %>% as.factor() %>% levels()], cls_cols[plot_meta[, group_col] %>% as.factor() %>% levels()]) %>%
                        setNames(., c(reg_col, group_col))
    

    column_ha <- HeatmapAnnotation(df = plot_meta, col = anno_cols, annotation_height = unit(c(0.01, 0.01), "in"))
    ## Get the range of the heatmap legends
    fun1 <- colorRampPalette(viridis(3))##c("#FF00FF", "black", "#FFFF00")
	color_breaks <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5)
	col_fun = circlize::colorRamp2(color_breaks, fun1(9))

    font_size <- font_scale * 5
    htlist <- Heatmap(data, name = "Expression", 
        col = col_fun, na_col = "white",
        column_split = meta[, reg_col],
        cluster_rows = FALSE, cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 5),
        row_title_rot = 90, column_title = NULL, row_title_gp = gpar(fontsize = 12), 
        show_column_names = FALSE, show_row_names = show_rownames, column_names_gp = gpar(fontsize = 10),column_names_rot = 45, 
        width = unit(8, "in"),
        top_annotation = column_ha, 
        heatmap_legend_param = list(title = "Scaled expr", color_bar = "continuous", at = color_breaks),
        use_raster = TRUE, raster_quality = 3)


    if (!is.null(label_genes)){
    	for (ii in names(label_genes)){
            htlist <- htlist + 
            	rowAnnotation(link = anno_mark(at = which(rownames(data) %in% label_genes[[ii]]), 
                    labels = rownames(data)[which(rownames(data) %in% label_genes[[ii]])], 
                    labels_gp = gpar(fontsize = font_size), padding = unit(1, "mm")))
    	}
    }
        
    if (return.plot) {
        return(htlist)
    } else {
        pdf(paste0(output_dir, file_name, "_heatmap.pdf"), width = width, height = height)
        draw(htlist)
        dev.off() 
    }
}



plot_region_heatmap_left <- function(meta, data, reg_col = "region", group_col = "cluster", pseudotime_col = "pseudotime", anchor_reg = "FC", file_name, output_dir = outputdir, pdf_width = 12, show_rownames = TRUE, width = 14, height = 10, label_genes = NULL, font_scale = 1, row_meta = NULL, reg_order = c("FC", "MSC", "TC", "OcC"), return.plot = FALSE) {
	
	meta <- meta[meta[, reg_col] %in% reg_order, ]
	## Get the full bins based on the anchor regions
	meta <- lapply(reg_order, function(reg) {
		ref <- meta %>%
				filter(region == anchor_reg)
		ref[, reg_col] <- reg
		rownames(ref) <- paste0(reg, "|", ref[, pseudotime_col])
		ref
		}) %>%
			do.call(rbind, .)

	meta[, reg_col] <- factor(as.character(meta[, reg_col]), levels = reg_order) 
	order_knots <- order(meta$region, meta$pseudotime)
	meta <- meta[order_knots, ]


	## Get the full bin data
	mis_bins <- setdiff(rownames(meta), colnames(data))
	mis_data <- matrix(NA, nrow = nrow(data), ncol = length(mis_bins), dimnames = list(rownames(data), mis_bins))
	data <- cbind(data, mis_data) %>%
				.[, rownames(meta)]


    ## Do the plots
    plot_meta <- meta[, c(group_col, reg_col)]
    regcols <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4") %>% 
				setNames(., c("FC", "MSC", "TC", "OcC", "Insula", "GE"))
	cls_cols <- c("#990939", "#e25a9a", "#fc19cf", "#1ac9a9", "#a3c10d", 
	"#efbc88", "#dbb75c", "#6700e5", "#7133aa", "#003563"
	) %>%
			setNames(., c("RGC FABP7 PMP22", 
				"NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", 
				"oRG HOPX TNC", "oRG HOPX APOE", 
				"tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))

    anno_cols <- list(regcols[plot_meta[, reg_col] %>% as.factor() %>% levels()], cls_cols[plot_meta[, group_col] %>% as.factor() %>% levels()]) %>%
                        setNames(., c(reg_col, group_col))
    

    column_ha <- HeatmapAnnotation(df = plot_meta, col = anno_cols, annotation_height = unit(c(0.01, 0.01), "in"))
    ## Get the range of the heatmap legends
    fun1 <- colorRampPalette(viridis(3))##c("#FF00FF", "black", "#FFFF00")
	color_breaks <- c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5)
	col_fun = circlize::colorRamp2(color_breaks, fun1(9))

    font_size <- font_scale * 5

    if (!is.null(label_genes)){
    	htlist <- rowAnnotation(link = anno_mark(at = which(rownames(data) %in% label_genes[[1]]), 
                    labels = rownames(data)[which(rownames(data) %in% label_genes[[1]])], 
                    side = "left",
                    labels_gp = gpar(fontsize = font_size), 
                    padding = unit(1, "mm"))) +
    			Heatmap(data, name = "Expression", 
			        col = col_fun, na_col = "white",
			        column_split = meta[, reg_col],
			        cluster_rows = FALSE, cluster_columns = FALSE, 
			        row_names_gp = gpar(fontsize = 5),
			        row_title_rot = 90, column_title = NULL, row_title_gp = gpar(fontsize = 12), 
			        show_column_names = FALSE, show_row_names = show_rownames, column_names_gp = gpar(fontsize = 10),column_names_rot = 45, 
			        width = unit(8, "in"),
			        top_annotation = column_ha, 
			        heatmap_legend_param = list(title = "Scaled expr", color_bar = "continuous", at = color_breaks),
			        use_raster = TRUE, raster_quality = 3)
    }
        
    if (return.plot) {
        return(htlist)
    } else {
        pdf(paste0(output_dir, file_name, "_heatmap.pdf"), width = width, height = height)
        draw(htlist)
        dev.off() 
    }
}


Anno_gene <- function(gene_label, mm) {
    ## Build the Gene Metadata
    gene_label <- gene_label[rownames(mm)]
    gene_meta <- data.frame(gene = names(gene_label), 
                    module = gene_label, 
                    modulemembership = sapply(names(gene_label), function(g) round(mm[g, gene_label[g]], digits = 4)), 
                    stringsAsFactors = FALSE) %>%
                    mutate(id = paste0(module, "::", gene, "::", as.character(modulemembership)))
    colnames(gene_meta)[colnames(gene_meta) == "modulemembership"] <- "mm"


    ## Get the order of genes within each module
    all_modules <- levels(as.factor(gene_label))
    gene_meta$gorder <- NA
    all_order <- c()
    for (mod in all_modules){
        idx <- which(gene_label == mod)
        cur_order <- gene_meta$mm[idx] %>% setNames(., gene_meta$gene[idx]) %>% sort(decreasing = TRUE) %>% names()
        gene_meta$gorder[idx] <- match(cur_order, gene_meta$gene)
        all_order <- c(all_order, cur_order)
    }
    gene_meta$all_order <- match(all_order, gene_meta$gene)
    rownames(gene_meta) <- gene_meta$gene
    return(gene_meta)
}




get_modules <- function(data, minClusterSize = 5, sensitivity = 4, file_name, output_dir = outputdir){
    ## Store the parameters [in a list]
    tree.method = "ward.D2"
    cor_method = "p"
    para <- list(minClusterSize = minClusterSize, tree.method = tree.method, sensitivity = sensitivity, cor_method = cor_method)


    ##Cluster the genes [Remove the na_cols first]
    data_use <- t(data)
    dissTOM <- 1-cor(data_use, method = cor_method)
    geneTree <- dissTOM %>% as.dist() %>% hclust(method = tree.method)


    ## Get the gene modules
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = sensitivity, pamRespectsDendro = FALSE,
                            minClusterSize = minClusterSize);
    new_mods <- dynamicMods
    names(new_mods)[names(new_mods) == ""] <- "0"
    mes <- moduleEigengenes(data_use, colors = paste0("_", new_mods))$eigengenes
    mes <- mes %>% 
            as.matrix() %>% scale() %>%
            MinMax(., min = -2.5, max = 2.5) %>%
            as.data.frame(., check.names = FALSE)
    ##print(rownames(mes))


    ## Get the gene annotation
    gene_label <- setNames(paste0("ME_", dynamicMods), colnames(data_use))
    gene_label[gene_label == "ME_"] <- "ME_0"
    mm <- as.data.frame(cor(data_use, mes, use = "p"));
    gene_meta <- Anno_gene(gene_label = gene_label, mm = mm)


    ## Plot the tree of hierachical clustering
    hc <- hclust(as.dist(1-cor(mes)), method = tree.method)
    pdf(paste0(output_dir, file_name, "_module_tree.pdf"), width = 5, height = 4)
    p <- ggdendro::ggdendrogram(hc) + geom_hline(yintercept = 0.25, col = "red")
    print(p)
    dev.off()


    out_res <- list(data = data, mes = mes[colnames(data),], mm = mm[rownames(data),], gene_meta = gene_meta[rownames(data),], hc = hc)
    return(out_res)
}



