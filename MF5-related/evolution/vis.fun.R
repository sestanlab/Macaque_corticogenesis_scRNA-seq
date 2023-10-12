###############################################################################################

## Dot Plot (scRNA-seq version)

###############################################################################################
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



#Function in Seurat package, useful
MinMax <- function (data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
}


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

	return(pre.plot)
}


## Extract Gene expression metadata
ExtractDotExpr <- function(avg, ratio, feature, spcols = c(Human = "red", Macaque = "blue")){
	## Some parameters
	col.min = -2.5
	col.max = 2.5
	ngradient = 20
	allctps <- c("ExN deep", "ExN upper")
	allregs <- c("FC", "OC")
	allsps <- c("Human", "Macaque")


	## All cluster order
	kp_cls <- lapply(allsps, function(sp) {
		aa <- rep(allregs, each = length(allctps)) %>%
			paste0(., "|", rep(allctps, by = length(allregs))) %>%
			paste0(sp, "|", .)
		return(aa)
		}) %>%
		unlist()


	## Scale avg data
	avg <- avg[feature, , drop = FALSE] %>%
					t() %>% scale() %>% t() %>%
					MinMax(data = ., min = col.min, max = col.max) %>%
					as.matrix()
	data.plot <- CombnMelt(avg = avg[, kp_cls, drop = FALSE], 
							ratio = ratio[feature, kp_cls, drop = FALSE])
	data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = feature))
	data.plot$avg.exp.scaled <- as.numeric(x = cut(x = data.plot$avg.exp.scaled, breaks = ngradient))
	
	cols_use <- colorRampPalette(colors = c("grey", "red"))(ngradient)
	sp_reg_ord <- rep(allsps, each = length(allregs)) %>%
					paste0(., "|", rep(allregs, length(allsps)))
	##print(sp_reg_ord)
	data.plot <- data.plot%>%
					mutate(species = extract_field(newcls, 1, "|"), 
							region = extract_field(newcls, 2, "|"), 
							cluster = extract_field(newcls, 3, "|"))%>%
					select(-newcls) %>%
					mutate(species = factor(species, levels = allsps)) %>%
					mutate(region = factor(region, levels = allregs)) %>%
					mutate(cluster = factor(cluster, levels = allctps)) %>%
					mutate(species_region = paste0(species, "|", region)) %>%
					mutate(species_region = factor(species_region, levels = sp_reg_ord))

	## Assign colors
	data.plot$colors <- mapply(FUN = function(color, value) {
		return(colorRampPalette(colors = c("lightgrey", color))(ngradient)[value])
	}, color = spcols[data.plot$species], value = data.plot$avg.exp.scaled)
	data.plot$colors[is.na(data.plot$colors)] <- "#FFFFFF"


	## Set pct.exp range
	data.plot$pct.exp <- data.plot$pct.exp * 100
	return(data.plot)
}


PlotDot <- function(data.plot, dot.scale = 5, dot.min = 0, font.scale = c(0.9, 0.8)) {
	data.plot$pct.exp[data.plot$pct.exp < (dot.min * 100)] <- NA

	p <- ggplot(data = data.plot, mapping = aes_string(x = "cluster", y = "features.plot")) + 
			geom_point(mapping = aes_string(size = "pct.exp", color = "colors"), shape = 16) + 
			scale_color_identity() +
			scale_radius(range = c(0, dot.scale)) + 
			theme_bw() +
			RotatedAxis() +
			facet_wrap(vars(species_region), nrow = 1, ncol = 4) +
			theme(panel.grid.minor = element_line(size = 0.1, color = "grey"),
					panel.grid.major = element_line(size = 0.1, color = "grey"),
					axis.title = element_blank(),
					axis.line = element_line(size = 0.2), 
					axis.ticks = element_line(size = 0.2),
					axis.text.x=element_text(size = rel(font.scale[1])), 
					axis.text.y=element_text(size = rel(font.scale[2])), 
					strip.background = element_blank(), strip.text = element_text(size = rel(0.6)), 
					panel.spacing = unit(0.1, "in"), 
					legend.text=element_text(size= rel(0.6)),
					legend.title=element_text(size= rel(0.6)),
					plot.title = element_text(size = rel(0.7)),
					legend.position = "bottom", 
					legend.box = "vertical", 
					legend.margin = margin()) + 
			guides(size = guide_legend(title = "% expressed"))
	return(p)
}



PlotDEGAnno <- function(dglist) {
	features <- unique(unlist(dglist))
	mat <- matrix(0, nrow = length(features), ncol = length(dglist), dimnames = list(features, names(dglist)))
	for (ii in names(dglist)){
		mat[dglist[[ii]], ii] <- 1
	}

	## Transform the data
	pdata <- as.matrix(mat) %>%
				reshape2::melt() %>%
				setNames(., c("gene", "type", "value")) %>%
				mutate(gene = factor(as.character(gene), levels = rev(features))) %>%
				mutate(type = factor(as.character(type), levels = names(dglist))) %>%
				mutate(value = as.character(value))

	p1 <- ggplot(pdata, aes_string(x = "type", y = "gene", fill = "value")) +
			geom_tile(width = 1, height = 1, size = 0.1, color = "white") +
			scale_fill_manual(values = c(`0` = "lightgrey", `1` = "black")) +
			theme_classic() +
			RotatedAxis() + 
			theme(legend.position = "none", 
				axis.text.x = element_text(size = rel(1)), 
				axis.text.y = element_text(size = rel(0.6)), 
				axis.ticks.y = element_line(linewidth = 0.2),
				axis.line = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank())
	return(p1)
}

