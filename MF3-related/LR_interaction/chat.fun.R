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


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}



##--------------------------------------------------------------------------
## Sankey and circle plots
PreCircleData <- function(mat, LR, log10transfer, maxp = 2.5, minp = -log10(0.05)) {
	## Check if the LR pair is present
	if (!LR %in% rownames(mat)){
		stop("Please make sure the input LR is present in the rownames of mat")
	}

	if (log10transfer){
		mat <- -log10(mat)
	}


	if (!is.null(minp)){
		mat[mat < minp] <- 0
	}


	## Prepare the data
	pdata <- mat[LR, , drop = FALSE] %>%
				as.matrix() %>%
				reshape2::melt() %>% 
				setNames(., c("ID", "cluster_pair", "mlogp")) %>%
				mutate(mlogp = MinMax(mlogp, min = 0, max = maxp)) %>%
				mutate(cluster_pair = as.character(cluster_pair)) %>%
				mutate(from = extract_field(cluster_pair, 1, "|")) %>%
				mutate(to = extract_field(cluster_pair, 2, "|")) 
	pdata <- pdata[, c("from", "to", "mlogp")]
	return(pdata)
}



PlotInterCircle_individual <- function(pdata, output_dir = outputdir, file_name, cluster.cols, cls.ord, title_lab = NULL, return.plot = FALSE) {
	max_sec <- 10
	max_sum1 <- tapply(pdata$mlogp, pdata$from, sum) %>% max(., na.rm = TRUE)
	max_sum2 <- tapply(pdata$mlogp, pdata$to, sum) %>% max(., na.rm = TRUE)
	if (max(c(max_sum1, max_sum2)) > 10){
		max_sec <- max(c(max_sum1, max_sum2))
	}
	xmax_values <- rep(max_sec, length(cluster.cols)) %>% setNames(., names(cluster.cols))


	if (!return.plot){
		pdf(paste0(output_dir, file_name, ".pdf"), width = 6, height = 6)
	}
	
		
	circos.par(start.degree = 90, cell.padding = c(0, 0, 0, 0), gap.degree=5.5)
	chordDiagram(pdata, grid.col = cluster.cols, order = cls.ord, xmax = xmax_values, reduce = -1, col = "black", scale = FALSE, annotationTrack = c("grid"), annotationTrackHeight = convert_height(7, "mm"), transparency = 0.6, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", link.arr.length= 0.05, link.arr.width = 0.01)#, diffHeight = -mm_h(2))
	
	for(si in get.all.sector.index()) {
				xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
				ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
				circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1,facing = "bending.inside", niceFacing = TRUE, col = "white", cex = 1.3)
	}

	if (!is.null(title_lab)){
		title(title_lab, cex = 1)
	}
	circos.clear()

	if (!return.plot){
		dev.off()
	}
}



plot_ggsankey <- function(meta, x_col, next_col, x_ord = NULL, next_ord = NULL, type = "sankey", space= 0.01, width= 0.1) {
    let_ord <- paste0(rep(letters, each = 9), rep(1:9, times = length(letters)))
    if (is.null(x_ord)){
        x_ord <- levels(as.factor(as.character(meta[, x_col])))
    }
    x_trans <- paste0(let_ord[1:length(x_ord)], x_ord) %>%
                    setNames(., x_ord)
    meta$newx <- x_trans[as.character(meta[, x_col])]



    if (is.null(next_ord)){
        next_ord <- levels(as.factor(as.character(meta[, next_col])))
    }
    next_trans <- paste0(let_ord[1:length(next_ord)], next_ord) %>%
                    setNames(., next_ord)
    meta$newnext <- next_trans[as.character(meta[, next_col])]


    df <- meta %>%
                make_long(newx, newnext) %>%
                mutate(nlabel = substring(node, 3))

    cols <- c(gg_color_hue(length(x_ord)), rep("#bdbdbd", length(next_trans))) %>%
                setNames(., c(setNames(x_trans, NULL), setNames(next_trans, NULL)))
    p <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = nlabel)) +
            geom_sankey(flow.alpha = .6, node.color = "gray30", type = type, width = width) +
            geom_sankey_text(size = 3, color = "black") +
            scale_fill_manual(values = cols) +
            theme_sankey(base_size = 18) +
            labs(x = NULL) +
            theme(legend.position = "none",  plot.title = element_text(hjust = .5), axis.text = element_blank())
    return(p)
}




