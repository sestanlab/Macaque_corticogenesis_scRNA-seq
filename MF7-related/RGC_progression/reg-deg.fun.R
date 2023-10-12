library(ComplexHeatmap)


###############################################################################################

## Dot Plot (modified version)

###############################################################################################
## Validate if the genes are contained in the dataset
ValidateFeatures <- function(features, all.genes) {
    ## Filter genes
    rm_genes <- setdiff(features, all.genes)

    ## Check features
    if (length(rm_genes) >= 1){
        warning(paste0("The following features are missing: ", paste(rm_genes, collapse = ", ")))
        features <- features[features %in% all.genes]
    }
    if (length(features) == 0){
        stop("No genes were found in the dataset")
    }
    return(features)
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

        avg <- avg[, rownames(plot_meta)]
        if (!is.null(ratio)){
            ratio <- ratio[, rownames(plot_meta)]
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
CombnMelt <- function(avg, ratio, mask = NULL) {
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


    ## Add mask data if present
    if (!is.null(mask)){
        pre.plot <- Mat2Df(mat = mask, name = "mask") %>%
                        dplyr::full_join(pre.plot, ., by = c("features.plot", "newcls")) 
    }


    ## Further format the columns
    data.plot <- pre.plot  %>%
                    mutate(species = extract_field(newcls, 1, "|"),  cluster = extract_field(newcls, 2, "|"))%>%
                    select(-newcls)

    return(data.plot)
}



CirclePlot.horizontal <- function(avg, ratio, features, file_name, dot.min = 0, dot.scale = 4, scale.by = "radius", shape = 16, cluster.order = NULL, stroke.col = "black", stroke.size = 0.5, col.min = -2.5, col.max = 2.5, scale = TRUE, scale.min = NA, scale.max = NA, mask.matrix = NULL, return.plot = FALSE, width.scale = 1, height.base = 1.5, font.scale = c(1, 1), height.unit = 0.3){
    split.order <- c("FC", "MSC", "TC", "OC")
    cols <- setNames(c("#FF420E","#4CB5F5","#89DA59","#FFBB00"), c("FC", "TC", "OC", "MSC"))
    ngradient <- 20

    ## Validate features
    features <- ValidateFeatures(features = features, all.genes = rownames(avg))
    

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
    res <- FillNACol(avg = data.use, ratio = ratio[features, , drop = FALSE], subset.ident = cluster.order)


    ## Matrix to DF
    data.plot <- CombnMelt(avg = res$avg, ratio = res$ratio, mask = mask.matrix[features, , drop = FALSE])


    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = features))
    data.plot$avg.exp.scaled <- as.numeric(x = cut(x = data.plot$avg.exp.scaled, breaks = ngradient))


    ## Set colors
    if (is.null(names(cols))){
        warning("The cols should be a named vector")
        cols <- setNames(cols, split.order)
    }
    data.plot$colorraw <- mapply(FUN = function(color, value) {
    				return(color)
                    #return(colorRampPalette(colors = c("grey", color))(ngradient)[value])
                    }, color = cols[data.plot$species], value = data.plot$avg.exp.scaled)
    data.plot$colors <- ifelse(data.plot$mask == 1, data.plot$colorraw, "lightgrey")

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
        data.plot$species <- factor(as.character(data.plot$species), levels = split.order)
    } 
    if (!is.null(cluster.order)){
        data.plot$cluster <- factor(as.character(data.plot$cluster), levels = cluster.order)
    }


    ##------------------------------------------------------------------
    ## Plot 
    ## Set dot scale functions
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))

    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "cluster")) + 
            geom_point(mapping = aes_string(size = "pct.exp", color = "colors"), shape = shape) + 
            scale_color_identity()
    plot <- plot +
            scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
            coord_flip() +
            theme_cowplot() +
            RotatedAxis() +
            facet_wrap(vars(species), nrow = 1, ncol = 4, scales = "free_x") +
            theme(panel.grid.minor = element_line(size = 0.1, color = "grey"),
                panel.grid.major = element_line(size = 0.1, color = "grey"),
                axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y=element_text(size = rel(0.75)*font.scale[2]), axis.text.x=element_text(size = rel(0.6 * font.scale[1])), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(0.1, "in"), legend.position = "bottom", axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2)) + 
            guides(size = guide_legend(title = "% expressed"))


    if (return.plot){
        return(plot)
    } else {
        pdf(paste0(outputdir, file_name, "_DotHori.pdf"), width = 15 * width.scale, height = ceiling(height.unit * length(features)) + height.base, useDingbats = FALSE)
        print(plot)
        dev.off()
    }
    
}
 
