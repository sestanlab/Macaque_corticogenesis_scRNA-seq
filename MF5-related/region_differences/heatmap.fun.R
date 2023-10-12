library(ComplexHeatmap)
library(circlize)
library(WGCNA)


###############################################################################################

## Some basic functions

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



gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


###############################################################################################

## Cluster genes

###############################################################################################
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



get_modules <- function(data, meta, gene_use = NULL, minClusterSize = 5, tree.method = "ward.D2", sensitivity = 4, cor_method = "p", file_name, plot_regions = c("FR", "MS", "TP", "OX")){
    ## Store the parameters [in a list]
    para <- list(minClusterSize = minClusterSize, tree.method = tree.method, sensitivity = sensitivity, cor_method = cor_method)


    ## Cluster the genes [Remove the na_cols first]
    if (!is.null(gene_use)){
        data_use <- t(data[gene_use, ,drop = FALSE])
    } else {
        data_use <- t(data)
    }


    ## Remove NA samples
    ##na_names <- rownames(data_use)[is.na(data_use[, 1])]
    ##data_use <- data_use[setdiff(rownames(data_use), na_names), ]
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



    gene_label <- setNames(paste0("ME_", dynamicMods), colnames(data_use))
    gene_label[gene_label == "ME_"] <- "ME_0"
    mm <- as.data.frame(cor(data_use, mes, use = "p"));
    gene_meta <- Anno_gene(gene_label = gene_label, mm = mm)


    ## Plot the tree of hierachical clustering
    hc <- hclust(as.dist(1-cor(mes)), method = tree.method)
    ##pdf(paste0(output_dir, file_name, "_module_tree.pdf"), width = 5, height = 4)
    p <- ggdendro::ggdendrogram(hc) + geom_hline(yintercept = 0.25, col = "red")
    print(p)
    ##dev.off()


    out_res <- list(data = data, mes = t(as.matrix(mes[colnames(data), ])), mm = as.matrix(mm[rownames(data),]), gene_meta = gene_meta[rownames(data), ], para = para, hc = hc)
    return(out_res)
}




###############################################################################################

## Heatmap visualization

###############################################################################################
plot_heatmap.end <- function(mat, label_genes, color_breaks = seq(-1.5, 2.5, 0.5), file_name, module_labs = NULL, pdf_height = 5) {
	## Set color breaks
	col_fun1 = colorRamp2(color_breaks, viridis(length(color_breaks)))


	lobe_cols <- c("#B037C4", "#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5") %>% setNames(., c("GE", "FC", "MSC", "OC", "Insula", "TC"))
	reg_cols <- rep(c("#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5"), c(4, 4, 1, 1, 3)) %>% 
				setNames(., c("MFC", "OFC", "DFC", "VFC", 
					"M1C", "S1C", "IPC", "PCC",
					"V1C",
					"Insula",
					"A1C", "ITC", "STC"))
	cls_cols <- setNames(c("#bf812d", "#bf812d", 
                    "#f56122", "#e89bc4", "#de77ae", "#c51b7d", "#D3D3D3", "#D3D3D3",
                    "#b35806", "#91ebe2", "#48a1e0", "#0868ac", "#c6dbef", "#c6dbef", "#4eb3d3"),
                    c("IPC EOMES VIM", "IPC EOMES NEUROG1", 
                    "IPC EOMES NHLH1 up", "ExN up nascent", "ExN up ADRA2A", "ExN up ACTN2", "ExN PCC NR4A3", "ExN up KCNV1",
                    "IPC EOMES NHLH1 deep", "ExN deep nascent", "ExN deep KIF26A", "ExN deep SYT6", "ExN deep OPRK1 SULF1","ExN deep OPRK1 NR4A2", "ExN deep NR4A2 GRID2"))


	## The columns of the meta data has to be in "factor" format
	meta <- data.frame(cell = colnames(mat), 
						stringsAsFactors = FALSE) %>%
					mutate(region = extract_field(cell, 1, "|")) %>%
					mutate(cluster = extract_field(cell, 2, "|")) %>%
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
	column_ha <- HeatmapAnnotation(link = anno_text(x = as.character(meta$region), which = "column", location = 0.5, just = "center", gp = gpar(fontsize = 10, fill = reg_col_use[as.character(meta$region)], col = "white", border = "black")),
								df = meta, col = col_list,
								name = "ii")

	
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
		
	if (!is.null(module_labs)){
		htlist <- htlist + 
					right_ha2 +
					right_ha1
	}

	return(htlist)
	##pdf(paste0(outputdir, file_name, "_heatmap.pdf"), width = 8, height = pdf_height)
	##draw(htlist)
	##dev.off() 
}




##---------------------------------------------------------------
## Prepare metadata for the pseudobulk data
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



## Prepare meta data for pseudobulk data (PB) for further smoothing
PrepPBmeta <- function(object, avg_col, ptime_col, bin_col, group.by) {
	pmeta <- data.frame(avgcls = unique(object@meta.data[, avg_col]), 
					stringsAsFactors = FALSE) %>%
					mutate(lobe = extract_field(avgcls, 1, "|")) %>%
					mutate(bin = as.numeric(extract_field(avgcls, 2, "|")))

	get_most <- function(x) {
	    table(x) %>% sort() %>% rev() %>% .[1] %>% names()
	}
	labels <- aggregate(as.formula(paste0(group.by, " ~ ", bin_col)), object@meta.data, get_most) %>%
				mutate(bin = as.character(!!sym(bin_col))) %>%
				mutate(cluster = !!sym(group.by))
	times <- aggregate(as.formula(paste0(ptime_col, " ~ ", bin_col)), object@meta.data, mean) %>%
				mutate(bin = as.character(!!sym(bin_col))) %>%
				mutate(pseudotime = !!sym(ptime_col))


	pmeta$pseudotime <- times$pseudotime[match(as.character(pmeta$bin), times$bin)]
	pmeta$cluster <- times$cluster[match(as.character(pmeta$bin), times$bin)]
	rownames(pmeta) <- pmeta$avgcls	
	return(pmeta)
}



## Get smooth gene expression along the trajectory
SmoothExprAcrossRegions <- function(avgs, meta, reg_ord = c("FC", "MSC", "TC", "OC"), span = 0.75) {
	reg_col <- "lobe"
	ptime_col <- "pseudotime"
	meta$pseudotime <- meta[, ptime_col]


	## The predicted df rownames will be the column names of the predicted data
	pt_min <- round(min(meta$pseudotime), digits = 6)
	if (pt_min < min(meta$pseudotime)){
		pt_min <- pt_min + 5e-6
	} 
	pt_max <- round(max(meta$pseudotime), digits = 6)
	if (pt_max > max(meta$pseudotime)){
		pt_max <- pt_max - 5e-6
	}


	ntime_base <- seq(pt_min, pt_max, length.out = 100) %>%
					sapply(., function(x) round(x, digits = 6))



	new_avg <- lapply(reg_ord, function(reg) {
		sub_meta <- meta[meta[, reg_col] == reg, ]
		sub_meta <- sub_meta[order(sub_meta$pseudotime), ,drop = FALSE]
		avg_sub <- avgs[, rownames(sub_meta)]


		ntime <- ntime_base[ntime_base > min(sub_meta$pseudotime) & ntime_base < max(sub_meta$pseudotime)]		
		predict_df <- data.frame(row.names = paste0(reg, "|", ntime), 
				pseudotime = ntime)
		ptime <- setNames(sub_meta$pseudotime, rownames(sub_meta))


		sub_smt <- SmoothBinAvg_combine(data = avg_sub, ptime = ptime, predict_df = predict_df, loess_span = span)
	return(sub_smt)
	}) %>%
		do.call(cbind, .)
	return(new_avg)
}




## Prepare the metadata for the smooth gene expresion data
PrepSMmeta <- function(avgs, object, ptime_col, group.by, ncells = 25) {
	object@meta.data$pseudotime <- object@meta.data[, ptime_col]
	meta <- data.frame(avgcls = colnames(avgs), stringsAsFactors = FALSE) %>%
							mutate(region = extract_field(avgcls, 1, "|")) %>%
							mutate(pseudotime = as.numeric(as.character(extract_field(avgcls, 2, "|")))) %>%
							column_to_rownames("avgcls")
	meta$cluster <- sapply(1:nrow(meta), function(xx) {
		high_cells <- object@meta.data[object@meta.data$pseudotime > meta$pseudotime[xx], ] %>%
						arrange(pseudotime) %>%
						.[, group.by] %>% .[1:ncells]
		low_cells <- object@meta.data[object@meta.data$pseudotime < meta$pseudotime[xx], ] %>%
						arrange(desc(pseudotime)) %>%
						.[, group.by] %>% .[1:ncells]
		cls <- c(high_cells, low_cells) %>%
					table() %>% sort() %>% 
					rev() %>% .[1] %>% names()
		cls
		})
	return(meta)
}




## Plot the region heatmap (continuous)
plot_heatmap.continuous <- function(mat, meta, label_genes, color_breaks = seq(-1.5, 2.5, 0.5), file_name, module_labs = NULL, pdf_height = 5, row_split = NULL) {
	## Set color breaks
	col_fun1 = colorRamp2(color_breaks, viridis(length(color_breaks)))


	lobe_cols <- c("#B037C4", "#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5") %>% setNames(., c("GE", "FC", "MSC", "OC", "Insula", "TC"))
	reg_cols <- rep(c("#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5"), c(4, 4, 1, 1, 3)) %>% 
				setNames(., c("MFC", "OFC", "DFC", "VFC", 
					"M1C", "S1C", "IPC", "PCC",
					"V1C",
					"Insula",
					"A1C", "ITC", "STC"))
	cls_cols <- setNames(c("#bf812d", "#bf812d", 
                    "#f56122", "#e89bc4", "#de77ae", "#c51b7d", "#D3D3D3", "#D3D3D3",
                    "#b35806", "#91ebe2", "#48a1e0", "#0868ac", "#c6dbef", "#c6dbef", "#4eb3d3"),
                    c("IPC EOMES VIM", "IPC EOMES NEUROG1", 
                    "IPC EOMES NHLH1 up", "ExN up nascent", "ExN up ADRA2A", "ExN up ACTN2", "ExN PCC NR4A3", "ExN up KCNV1",
                    "IPC EOMES NHLH1 deep", "ExN deep nascent", "ExN deep KIF26A", "ExN deep SYT6", "ExN deep OPRK1 SULF1","ExN deep OPRK1 NR4A2", "ExN deep NR4A2 GRID2"))


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

	if (!is.null(module_labs)){
		htlist <- htlist + 
					right_ha2
	}
	htlist <- htlist +
					right_ha1

	return(htlist)
	##pdf(paste0(outputdir, file_name, "_heatmap.pdf"), width = 10, height = pdf_height)
	##draw(htlist)
	##dev.off() 
}



###------------------------------------------------------------------------

### Impulse fit

###------------------------------------------------------------------------
library(GetoptLong)
library(ggpubr)
library(parallel)
library(foreach)
ParallelImpulseFit <- function(fit.expression, genes, ptime, k = 15, interpolate = NULL, nCores = 6){
    ## # Do impulse model fitting for all genes 
    message(paste0(Sys.time(), ": Fitting impulse model for all genes."))
    message(Sys.time())
    cl = makeCluster(nCores, type = "FORK", outfile="")
    doParallel::registerDoParallel(cl);

    impulse.fits <- foreach(gg = genes, .inorder = TRUE, .packages = c("URD"), .verbose = TRUE) %dopar% {
            impulseFit(x = ptime, y = as.numeric(fit.expression[gg, ]), k = k, interpolate = interpolate)
        }
    names(impulse.fits) <- genes
    print(names(impulse.fits))
    stopCluster(cl)
    message("Finish impulse model fitting analysis")
    message(Sys.time())
    return(impulse.fits)
}

 

GetCascadeFit <- function(meta, data, pseudotime_col = "pseudotime", ref_sps = NULL, scale.data = TRUE, verbose.genes = TRUE, interpolate = NULL, k = 15, nCores = 36) {
    if (scale.data){
        fit.expression <- data %>%
                t() %>% scale() %>% t() %>%
                MinMax(., min = -2.5, max = 2.5)
    } else {
        fit.expression <- data
    }

    ptime <- meta[colnames(fit.expression), pseudotime_col]
    genes <- rownames(data)

    ## # Do impulse model fitting for all genes 
    if (length(genes) > 36){
        g.bin <- ceiling(c(1:length(genes))/36)
        max.bin <- max(g.bin)
        fit.list <- list()
        for (ii in 1:max.bin){
        	print(ii)
            sub.genes <- genes[g.bin == ii]
            print(sub.genes)
            fit.list[[ii]] <- ParallelImpulseFit(fit.expression = fit.expression, genes = sub.genes, ptime = ptime, k = k, interpolate = interpolate, nCores = nCores)
        }
        impulse.fits <- do.call(c, fit.list)

        ## Order the fit by gene names again
        impulse.fits <- impulse.fits[genes]
    } else {
        impulse.fits <- ParallelImpulseFit(fit.expression = fit.expression, genes = genes, ptime = ptime, k = k, interpolate = interpolate, nCores = nCores)
    }



    # Get out onset/offset times  
    timing <- data.frame(
                time.on=unlist(lapply(impulse.fits, function(x) if(is.list(x)){return(min(x[['time.on']]))}else{return(x['time.on'])})),
                time.off=unlist(lapply(impulse.fits, function(x) if(is.list(x)){return(max(x[['time.off']]))}else{return(x['time.off'])})),
                row.names = genes, stringsAsFactors = FALSE)
    timing[intersect(which(is.na(timing$time.on)), which(is.infinite(timing$time.off))), "time.on"] <- Inf
    
    outs <- list(imfit = impulse.fits, timing = timing)
    return(outs)
}



###--------------------------------------------------------------------------------------

## Prepare Weighted dot plots

###--------------------------------------------------------------------------------------
spdf2arry <- function(df) {
	## Add empty cols for missing pairs (diagonal Human-Human, Chimp-Chimp, etc)
	all_regs <- c("FC", "MSC", "TC", "OC")
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



















