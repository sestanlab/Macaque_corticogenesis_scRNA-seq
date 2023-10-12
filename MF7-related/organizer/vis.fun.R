## Horizontal dotplot for PFC data (Human & Rhesus macaque)
RotatedAxis <- function (...) {
  rotated.theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1), validate = TRUE, ...)
  return(rotated.theme)
}



extract_field <- function(input_genes, field = 1, split_by = "_") {
  split_strings <- strsplit(input_genes, split_by, fixed=TRUE)
  if (is.numeric(field)){
    if (all(field > 0)){
      if (length(field) == 1){
        idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
        output_strings <- sapply(1:length(split_strings), function(x) split_strings[[x]][idx[x]])
      } else {
        idx <- sapply(split_strings, length) == 1
        output_strings <- sapply(1:length(split_strings), function(x) ifelse(idx[x], input_genes[x], paste(split_strings[[x]][field], collapse = split_by)))
      }
    } else {
      if (length(field) == 1) {
        field = as.integer(abs(field))
        idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
        output_strings <- sapply(1:length(split_strings), function(x) rev(split_strings[[x]])[idx[x]])
      } else {
        stop("currently doesnot support field with length >= 2")
      }
    }
  } else if (is.character(field)) {
    if (field == "rm_start"){
      idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
      output_strings <- sapply(1:length(split_strings), function(x) paste(split_strings[[x]][-1],collapse=split_by))
    } else if (field == "rm_end") {
      output_strings <- sapply(1:length(split_strings), function(x) paste(rev(rev(split_strings[[x]])[-1]),collapse=split_by))
    } else {
      stop("Currently only support rm_start, rm_end for the character field")
    }
  }
  output_strings <- output_strings %>% setNames(., input_genes)
  return(output_strings)
}


MinMax <- function (data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

## Because some clusters are not conserved across species, need to add some NA cols to make them pseudo-conserved
FillNACol <- function(avg, ratio = NULL, subset.ident = NULL) {
	## Metadata 
	plot_meta <- data.frame(spcls = colnames(avg), stringsAsFactors = FALSE) %>%
		mutate(cluster = extract_field(spcls, 2, "|"), species = extract_field(spcls, 1, "|"))
  
	## Subset the data if necessry
	if (!is.null(subset.ident)){
		plot_meta <- plot_meta %>%
				filter(cluster %in% subset.ident) %>%
				column_to_rownames("spcls")

		avg <- avg[, rownames(plot_meta), drop = FALSE]
		if (!is.null(ratio)){
			ratio <- ratio[, rownames(plot_meta), drop = FALSE]
		}

		dif_cls <- table(plot_meta$cluster) %>% .[. < 4] %>% names() %>% 
				intersect(., subset.ident)
	} else {
		plot_meta <- plot_meta %>%
				column_to_rownames("spcls")
		dif_cls <- table(plot_meta$cluster) %>% .[. < 4] %>% names()
	}
  
  
	## Make sure we have NA columns for non-homologous clusters
	sp_use <- unique(plot_meta$species)
  
  
	if (length(dif_cls) > 0) {
	dif_cols <- paste0(rep(sp_use, times = length(dif_cls)), "|", 
					rep(dif_cls, each = 4)) %>%
					setdiff(., colnames(avg))
	dif_mat <- matrix(NA, nrow = nrow(avg), ncol = length(dif_cols), dimnames = list(rownames(avg), dif_cols))
	avg <- cbind(avg, dif_mat)
	if (!is.null(ratio)){
		ratio <- cbind(ratio, dif_mat)
	}
    
	plot_meta <- data.frame(spcls = dif_cols, stringsAsFactors = FALSE) %>%
			mutate(cluster = extract_field(spcls, 2, "|"), species = extract_field(spcls, 1, "|")) %>%
			column_to_rownames("spcls") %>%
			rbind(., plot_meta)
	}
	return(list(avg = avg, ratio = ratio, meta = plot_meta))
}


## Combinatory melt of avg and ratio matrices
CombnMelt <- function(avg, ratio, stroke = NULL) {
	Mat2Df <- function(mat, name) {
		df <- mat %>% 
				as.matrix() %>%
				reshape2::melt() %>%
				setNames(., c("features.plot", "newcls", name)) %>%
				mutate(features.plot = as.character(features.plot), newcls = as.character(newcls))
		return(df)
	}
  
  
	## Combine average and expression ratio
	avg.exp <- Mat2Df(mat = avg, name = "avg.exp.scaled")
	pct_exp <- Mat2Df(mat = ratio, name = "pct.exp")
	pre.plot <- dplyr::full_join(avg.exp, pct_exp, by = c("features.plot", "newcls"))
  
  
	## Add stroke data if present
	if (!is.null(stroke)){
		pre.plot <- Mat2Df(mat = ratio, name = "stroke") %>%
			dplyr::full_join(pre.plot, ., by = c("features.plot", "newcls")) 
	}
  
  
	## Further format the columns
	data.plot <- pre.plot  %>%
			mutate(species = extract_field(newcls, 1, "|"),  cluster = extract_field(newcls, 2, "|"))##%>%
			##select(-newcls)

	return(data.plot)
}



CirclePlot.horizontal <- function(avg, ratio, features, file_name, dot.min = 0, dot.scale = 4, scale.by = "radius", shape = 16, cluster.order = NULL, stroke.col = "black", stroke.size = 0.5, col.min = -2.5, col.max = 2.5, scale = TRUE, scale.min = NA, scale.max = NA, stroke.matrix = NULL, return.plot = TRUE, width.scale = 1, height.base = 1.5, height.unit = 0.3, col.cex.sf = 1, row.cex.sf = 1, split.order = c("Human", "Rhesus")){
  	cols <- c("#FF420E", "#89DA59") %>% 
			setNames(., c("Human", "Rhesus"))
	split.order <- c("Human", "Rhesus")
	ngradient <- 20

  
  	##-------------------------------------------------------------------------
	## Set data
	avg <- as.matrix(avg[features, , drop = FALSE])
	rownames(avg) <- features
	ratio <- as.matrix(ratio[features, , drop = FALSE])
	rownames(ratio) <- features
	if (scale) {
		data.use <- avg %>%
			t() %>% scale() %>% t() %>%
			MinMax(data = ., min = col.min, max = col.max)
	} else {
		data.use <- log(x = avg)
	}
  
  
	## Fill NA cols (for non-conserved clusters)
	res <- FillNACol(avg = data.use, ratio = ratio, subset.ident = cluster.order)
	
	## Matrix to DF
	data.plot <- CombnMelt(avg = res$avg, ratio = res$ratio, stroke = stroke.matrix)
  
	data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = features))
	data.plot$avg.exp.scaled <- as.numeric(x = cut(x = data.plot$avg.exp.scaled, breaks = ngradient))
  

	## Set colors
	if (is.null(names(cols))){
		warning("The cols should be a named vector")
		cols <- setNames(cols, split.order)
	}
	data.plot$colors <- mapply(FUN = function(color, value) {
		return(colorRampPalette(colors = c("lightgrey", color))(ngradient)[value])
	}, color = cols[data.plot$species], value = data.plot$avg.exp.scaled)
	data.plot$colors[is.na(data.plot$colors)] <- "#FFFFFF"


	## Set pct.exp range
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    data.plot$pct.exp[data.plot$pct.exp < (dot.min*100)] <- NA


	## Define the order the species & cluster
	if (!is.null(split.order)){
		data.plot <- filter(data.plot, species %in% split.order)
		data.plot$species <- factor(as.character(data.plot$species), levels = split.order)
	} 
	if (!is.null(cluster.order)){
		data.plot$cluster <- factor(as.character(data.plot$cluster), levels = cluster.order)
	}

	##------------------------------------------------------------------
	## Plot 
	## Set dot scale functions
	scale.func <- switch(EXPR = scale.by, size = scale_size, radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  
	p1 <- ggplot(data = data.plot, mapping = aes_string(x = "cluster", y = "features.plot")) + 
				geom_point(mapping = aes_string(size = "pct.exp", color = "colors"), shape = shape) + 
				scale_color_identity() +
				scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
				##coord_flip() +
				theme_bw() +
				RotatedAxis() +
				coord_fixed(clip = "off") + 
				facet_grid(cols = vars(species)) + ##, nrow = 1, ncol = length(split.order)) +
				theme(panel.grid.major = element_line(size = 0.2),
					panel.grid.minor = element_line(size = 0.2),
					axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_text(size = rel(0.75 * row.cex.sf)), axis.text.x=element_text(size = rel(0.6 * col.cex.sf)), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(0.1, "in"), legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), legend.title = element_text(size = rel(0.7)), legend.text = element_text(size = rel(0.6))) + 
				guides(size = guide_legend(title = "% Expressed"))
	return(p1) 
}



CirclePlot.horizontal.Lister <- function(avg, ratio, features, file_name, dot.min = 0, dot.scale = 4, scale.by = "radius", shape = 16, cluster.order = NULL, split.order = NULL, stroke.col = "black", stroke.size = 0.5, col.min = -2.5, col.max = 2.5, scale = TRUE, scale.min = NA, scale.max = NA, mask.vector = NULL, return.plot = FALSE, width.scale = 1, height.base = 1.5, font.scale = c(1, 1), height.unit = 0.3, min.cells = 50){
    ngradient <- 20

    ##-------------------------------------------------------------------------
	## Set data
	avg <- as.matrix(avg[features, , drop = FALSE])
	rownames(avg) <- features
	ratio <- as.matrix(ratio[features, , drop = FALSE])
	rownames(ratio) <- features


    ##-------------------------------------------------------------------------
    ## Scale data
    if (scale) {
        data.use <- avg[features, , drop = FALSE] %>%
                        t() %>% scale() %>% t() %>%
                        MinMax(data = ., min = col.min, max = col.max)
    } else {
        data.use <- log(x = avg[features, , drop = FALSE])
    }


   ## Fill NA cols (for non-conserved clusters)
	##res <- FillNACol(avg = data.use, ratio = ratio, subset.ident = cluster.order)
	
	## Matrix to DF
	data.plot <- CombnMelt(avg = data.use, ratio = ratio, stroke = NULL)
	print(head(data.plot))
    data.plot$cluster_size <- mask.vector[data.plot$newcls]


    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = features))
    data.plot$avg.exp.scaled <- as.numeric(x = cut(x = data.plot$avg.exp.scaled, breaks = ngradient))


    ## Set colors
    #if (is.null(names(cols))){
    #    warning("The cols should be a named vector")
    #    cols <- setNames(cols, split.order)
    #}
    #data.plot$colorraw <- mapply(FUN = function(color, value) {
    #				return(color)
    #                }, color = cols[data.plot$age], value = data.plot$avg.exp.scaled)
    #data.plot$colors <- ifelse(data.plot$mask == 1, data.plot$colorraw, "lightgrey")

    ## I think you should set the colors using the following codes
    cols_use <- colorRampPalette(colors = c("lightgrey", "red"))(ngradient)
    ##cols_use <- viridis::viridis(ngradient) ##you can pick better color using the above function
    data.plot$colors <- cols_use[data.plot$avg.exp.scaled]


    ## Set pct.exp range
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    data.plot$pct.exp[data.plot$pct.exp < (dot.min*100)] <- NA



    ## Define the order the age & cluster
    if (!is.null(split.order)){
    	data.plot <- data.plot %>%
    					filter(species %in% split.order) %>%
    					mutate(species = factor(as.character(species), levels = split.order))
    } 
    if (!is.null(cluster.order)){
    	data.plot <- data.plot %>%
    					filter(cluster %in% cluster.order) %>%
    					mutate(cluster = factor(as.character(cluster), levels = cluster.order))
    }


    ## Decide the shape of the dot
    data.plot$shape <- ifelse(data.plot$cluster_size > min.cells, 16, 4)
    data.plot$pct.exp[data.plot$shape == 4] <- 25
    data.plot$colors[data.plot$shape == 4] <- cols_use[1]


    ##------------------------------------------------------------------
    ## Plot  
    ## Set dot scale functions
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))

    plot <- ggplot(data = data.plot, mapping = aes_string(x = "cluster", y = "features.plot")) + 
            geom_point(mapping = aes_string(size = "pct.exp", color = "colors", shape = "shape")) + 
            scale_color_identity() +
            scale_shape_identity()
    plot <- plot +
            scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
            ##coord_flip() +
            theme_cowplot() +
            RotatedAxis() +
            facet_grid(cols = vars(species), scales = "free_x") +
            theme(panel.grid.minor = element_line(size = 0.1, color = "grey"),
                panel.grid.major = element_line(size = 0.1, color = "grey"),
                axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_text(size = rel(0.75)*font.scale[2]), axis.text.x=element_text(size = rel(0.6 * font.scale[1])), strip.background = element_blank(), strip.text = element_text(size = rel(0.9)), panel.spacing = unit(0.1, "in"), legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2)) + 
            guides(size = guide_legend(title = "% expressed"))


    #pdf("./report/test_dotplot.pdf", width = 20, height = 4)
    #print(plot)
    #dev.off()

    if (return.plot){
        return(plot)
    } else {
        pdf(paste0(outputdir, file_name, "_DotHori.pdf"), width = 15 * width.scale, height = ceiling(height.unit * length(features)) + height.base, useDingbats = FALSE)
        print(plot)
        dev.off()
    }
    
}


