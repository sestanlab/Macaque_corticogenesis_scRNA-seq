work_dir <- getwd()
inputdir <- paste0(work_dir, "/load_files/")
outputdir <- paste0(work_dir, "/report/")
dataDir <- "~/project/NHPfetal/data/"


#Load all related packages 

## data organization related
library(tibble)
library(dplyr)
library(Matrix)
library(tidyr)


## plot
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(viridis);


## parallel computing
library(parallel); 


## Seurat
##library(batchelor)
#library(scran); 
if (sum(grepl("Seurat", (.packages()))) < 1){
    library(Seurat)#if source twice, the non-developmental version could be loaded and functions may change.
}


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}



#This scirpt contains basic functions of single cell RNA-seq analysis (as well as many other technologies) 
map_gene <- function(gene_names, input_genes,ignore_case=TRUE){
    input_genes <- unique(input_genes)
  
    if (sum(grepl("\\|",gene_names))==length(gene_names)){
        if (sum(grepl("\\|",input_genes))==length(input_genes)){
              gene_id <- input_genes[input_genes %in% gene_names]
          }else{
            input_genes <- extract_field(input_genes=input_genes, 2, "|")
              gene_id <- unlist(mclapply(input_genes, function(i) ifelse(length(grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
              gene_id <- gene_id[gene_id != "empty"]
          }
    } else if(sum(grepl("\\|",gene_names))==0){
        input_genes <- extract_field(input_genes=input_genes, 2, "|")
          gene_id <- unlist(mclapply(input_genes, function(i) ifelse(length(grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
          gene_id <- gene_id[gene_id != "empty"]
    } else {
        stop("Inconsistent gene name format")
    }
    return(gene_id)
}

#Get the mito genes based on the inquiry genes
get_genes <- function(input_genes, gene_type = c("mito","ribo", "cc")[1], return_list = FALSE, revised = FALSE){
    gene_use <- list()
    if ("mito" %in% gene_type){
        mito.known <- map_gene(gene_names=input_genes, input_genes=c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")) #Refseq annotation 103 (Macaque)
        mito.combine <- grep(pattern = "\\|MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        mito.single <- grep(pattern = "^MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        gene_use[["mito"]] <- unique(c(mito.known, mito.combine, mito.single))
    }

    if ("ribo" %in% gene_type){
        ribo.combine <- c(grep("\\|RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        ribo.single <- c(grep("^RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        gene_use[["ribo"]] <- c(ribo.combine, ribo.single)
    }

    if ("cc" %in% gene_type){
        if (!revised){
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        } else {
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        }
        
        gene_use[["s"]] <- map_gene(gene_names=input_genes, input_genes=regev.s.genes,ignore_case=TRUE)
        gene_use[["g2m"]] <- map_gene(gene_names=input_genes, input_genes=regev.g2m.genes,ignore_case=TRUE)
    }

    if (return_list){
        return(gene_use)
    } else {
        all_genes <- setNames(unlist(gene_use), NULL)
        return(all_genes)
    }
}



seu_prepare <- function(counts, data = NULL, min.cells = 5, normalization.method = c("LogNormalize", "none","scran")[1], nfeatures = 2500, hvg.method = "vst", assay = "RNA") {
    #Check the normalization.method
    if (!is.null(data)){
        message("input-norm is provided and therefore normalization is not required")
        normalization.method <- "none"
    }

    if (normalization.method == "none"){
        message("normalization.method is none and therefore assume the counts has already been normlized")
    } else if (normalization.method == "LogNormalize"){
        message("normalization.method is LogNormalize and therefore will use seurat LogNormalize to transform the dataset")
    } else {
        stop("Unknown normalization.method")
    }

    inseu <- CreateSeuratObject(counts, meta.data = NULL, assay = "RNA", min.cells = min.cells, min.features = 0, names.field = 1, names.delim = "_"); #at least expressed in 5 cells #only RNA assay is supported

    #quality_gene list
    quality_genes <- get_genes(input_genes = rownames(inseu$RNA@data), gene_type = c("mito","ribo", "cc"), return_list = TRUE, revised = FALSE)
    inseu[["percent.mt"]] <- PercentageFeatureSet(inseu, pattern = NULL, features = quality_genes[["mito"]], col.name = NULL)
    inseu[["percent.ribo"]] <- PercentageFeatureSet(inseu, pattern = NULL, features = quality_genes[["ribo"]], col.name = NULL)


    #normalize the dataset if needed
    if (normalization.method == "none"){
        inseu[[assay]] <- CreateAssayObject(data = as(ifelse_check(is.null(data), inseu[["RNA"]]@data, data), "dgCMatrix"), min.cells = 0, min.features = 0)
        DefaultAssay(inseu) <- assay
    } else if (normalization.method == "LogNormalize"){
        inseu <- NormalizeData(inseu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    }
    
    if (!is.null(hvg.method)){
        inseu <- FindVariableFeatures(inseu, selection.method = hvg.method, nfeatures = nfeatures, verbose = FALSE)
    }
    return(inseu)
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


#Function in Seurat package, useful
MinMax <- function (data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
}


remove_duplicates <- function(x, return_index = FALSE){
    dp_index <- !duplicated(x) & !rev(duplicated(rev(x)))
    if(return_index){
        return(dp_index)
    } else {
        return(x[dp_index])
    }
}


ifelse_check <- function(test, yes, no){if (test){return(yes)} else{ return(no) }}


mean.of.logs <- function (x, base = 2) {
    return(log(mean((base^x) - 1) + 1, base = base))
}




library(AUCell)
GetModuleScore <- function (assay.data, features, nbin = 24, ctrl = 100, k = FALSE, seed = 42, method = c("seurat","aucell")[2], input_dir = new_inputdir, file_name, output_dir = outputdir, rethreshold_list = NULL, cellbin.size = 8000) {
     if (is.null(x = features)) {
        stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
        return(intersect(x = x, y = rownames(x = assay.data)))
    })
    cluster.length <- length(x = features) #number of feature list

    if (method == "seurat"){
        set.seed(seed = seed)
        pool <- rownames(x = assay.data)
        data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
        data.avg <- data.avg[order(data.avg)]
        data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
            n = nbin, labels = FALSE, right = FALSE)
        names(x = data.cut) <- names(x = data.avg)
        ctrl.use <- vector(mode = "list", length = cluster.length)
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            for (j in 1:length(x = features.use)) {
                ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                    data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
            }
        }
        ctrl.use <- lapply(X = ctrl.use, FUN = unique)
        ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
            ncol = ncol(x = assay.data))
        for (i in 1:length(ctrl.use)) {
            features.use <- ctrl.use[[i]]
            ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
               ])
        }
        features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
            ncol = ncol(x = assay.data))
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            data.use <- assay.data[features.use, , drop = FALSE]
            features.scores[i, ] <- Matrix::colMeans(x = data.use)
        }
        features.scores.use <- features.scores - ctrl.scores
        features.scores.use <- as.data.frame(x = t(x = features.scores.use))
        rownames(x = features.scores.use) <- colnames(x = assay.data)
        colnames(features.scores.use) <- names(features)

        return(features.scores.use)
    } else if (method == "aucell"){
        library(AUCell)
        if (!file.exists(paste0(input_dir, file_name, "_modulescore_auc_res.rds"))){
            #Split the cells to bins[Sometimes, the matrix is too large for rankings]
            if (ncol(assay.data) < cellbin.size){
                cellbin.size <- ceiling(ncol(assay.data)/2)
            }
            bin.ind <- ceiling(c(1:ncol(assay.data))/cellbin.size)
            max.bin <- max(bin.ind)

            auc_list <- lapply(1:max.bin, function(tembin) {
                tem_matrix <- assay.data[, bin.ind == tembin]
                tem_rankings <- AUCell_buildRankings(tem_matrix, nCores=1, plotStats=FALSE) 
                tem_auc <- AUCell_calcAUC(features, tem_rankings)#, aucMaxRank = 500)
                tem_aucmatrix <- t(as.matrix(getAUC(tem_auc)))
                rm(tem_matrix, tem_rankings)
                return(tem_aucmatrix)
                })

            hauc_matrix <- do.call(rbind, auc_list)

            if (length(auc_list) == 1){ #When the input is not a named list, the colnames will be "Geneset"
                colnames(hauc_matrix) <- names(features)
            }
            saveRDS(hauc_matrix, file = paste0(input_dir, file_name, "_modulescore_auc_res.rds"))
        } else {
            hauc_matrix <- readRDS(paste0(input_dir, file_name, "_modulescore_auc_res.rds"))
        }
        
        set.seed(seed)
        pdf(paste0(output_dir, file_name, "_modulescore_auc_auto_assignment.pdf"), paper="letter")
        par(mfrow=c(2,2)); cells_assignment <- AUCell_exploreThresholds(t(hauc_matrix), plotHist=TRUE, assign=TRUE) 
        dev.off()

        #Generate assignment matrix (rownames as cells)
        default_assign <- hauc_matrix * 0 #build an empty one
        for (gset in colnames(hauc_matrix)){
            default_assign[cells_assignment[[gset]]$assignment, gset] <- 1
        }

        outdata <- list(auc = hauc_matrix, auto = default_assign, custom = default_assign)

        if (!is.null(rethreshold_list)){
            pdf(paste0(output_dir, file_name, "_modulescore_auc_custom_assignment.pdf"), paper="letter")
            for (j in names(rethreshold_list)){
                AUCell_plotHist(t(hauc_matrix)[j,,drop = FALSE], aucThr=rethreshold_list[[j]])
                abline(v=rethreshold_list[[j]])
                cells_assignment[[j]]$assignment <- rownames(hauc_matrix)[hauc_matrix[, j] >= rethreshold_list[[j]]]
            }
            dev.off()

            ##default_assign <- hauc_matrix * 0 #build an empty one
            for (gset in colnames(hauc_matrix)){
                default_assign[, gset] <- 0
                default_assign[cells_assignment[[gset]]$assignment, gset] <- 1
            }
            outdata$custom <- default_assign
        }
        return(outdata)
    }
}





auto_levels <- function(object, split.by) { 
    all_levels <- levels(as.factor(object@meta.data[, split.by]))
    level_list <- list(c("FR", "MS", "TP", "OX"), 
                        c("Frontal", "MotorSensory", "Temporal", "Occipital"), 
                        c("Human", "Chimpanzee", "Macaque", "Marmoset"),
                        c("Human", "Chimpanzee", "Rhesus", "Marmoset"),
                        c("HSB", "PTB", "RMB", "MMB"),
                        c("E37", "E42-43", "E54", "E62-E64", "E77-78"))
    ## get the idx of reference names
    nshare <- sapply(level_list, function(x) sum(all_levels %in% x))
    idx <- which(idx == max(idx))[1]

    if (max(nshare) == 9){
        return(all_levels)
    } else {
        new_levels <- level_list[[idx]][level_list[[idx]] %in% all_levels]
        return(new_levels)
    }
}




FeatureFig <- function(object, features, cols = c("#f5f5dc", "#31a354","#253494"), split.by = NULL, plot.scale = 0.9, ncol = NULL, file_name, output_dir = outputdir, pt.size = "auto", split.order = "auto", output.ggplot = FALSE, ngradient = 10, ...) {
    
    extra_arguments <- list(...)


    if (pt.size == "auto"){
        point.size <- round(10/sqrt(ncol(object)), digits = 1) %>% MinMax(., min = 0.1, max = 1)
    } else {
        point.size <- pt.size
    }

    ##Set the levels of factors in the supervised mode
    if (!is.null(split.by)){
        ##Get the supervised order the split.by
        new_levels <- ifelse_check("auto" %in% split.order, 
                            auto_levels(object = object, split.by = split.by),
                            split.order)
        object@meta.data[, split.by] <- factor(as.character(object@meta.data[, split.by]), levels = new_levels) 
    }


    cols <- colorRampPalette(cols)(ngradient)

    ##Get the plot list
    if (length(features) == 1){
        plot_args <- c(list(object = object, cols = cols, features = features, pt.size = point.size, combine = FALSE, split.by = split.by), extra_arguments[names(extra_arguments) %in% c("order", "reduction", "min.cutoff", "max.cutoff", "shape.by", "slot", "blend","blend.threshold", "label", "label.size", "repel", "combine", "coord.fixed", "by.col", "sort.cell")])
        plist <- do.call(FeaturePlot, plot_args)
        if (!is.null(split.by)){
            names(plist) <- paste0(features, "-", rep(new_levels, length.out = length(plist)))
        } else {
            names(plist) <- features
        }
    } else {
        plot_args <- c(list(object = object, cols = cols, features = features, pt.size = point.size, combine = FALSE, split.by = split.by), extra_arguments[names(extra_arguments) %in% c("order", "reduction", "min.cutoff", "max.cutoff", "shape.by", "slot", "blend","blend.threshold", "label", "label.size", "repel", "combine", "coord.fixed", "by.col", "sort.cell")])
        plist <- do.call(FeaturePlot, plot_args)


        ##Set the names of the plot list 
        if (!is.null(split.by)){
            new_plist <- lapply(1:(length(plist)/length(new_levels)), function(x) {
                gidx <- seq(x, length(plist), length(plist)/length(new_levels))
                plist[gidx]
                }) %>% do.call(c, .)
            plist <- new_plist; rm(new_plist)

            gene_label <- sapply(plist, function(p) ggplot_build(p)$plot$plot_env$plot$layers[[1]]$mapping$colour[2] %>% as.character())
            gene_label <- paste0(gene_label, "-", rep(new_levels, length.out = length(plist)))
        } else {
            gene_label <- sapply(plist, function(p) ggplot_build(p)$plot$plot_env$plot$layers[[1]]$mapping$colour[2] %>% as.character())
        }

        names(plist) <- gene_label
    }


    ##Polish the figures by changing the theme
    plist <- lapply(names(plist), function(x) {
        p <- plist[[x]] + 
            coord_equal(ratio = 1) + 
                theme_classic() + 
                scale_color_gradientn(colors = cols, limits = c(1,ngradient)) + 
                theme(legend.position = "bottom",
                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
                    plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))+ 
                labs(title =  text_wrapper(x, width=20))
        p
        })


    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist)
    } else {
        nlevel <- ifelse(is.null(split.by), 1, length(new_levels))
        ncol <- ifelse_check(is.null(ncol), nlevel, ncol)
        nrow <- ceiling(length(plist)/ncol)
        
        jpeg(paste0(output_dir, file_name, "_feature.jpeg"), width = 10 * plot.scale * ncol, height = 11 * plot.scale * nrow, units = "in", res = 300)
        plot_grid(plotlist = plist, nrow = nrow, ncol = ncol) %>% print() %>% print()
        dev.off()
    }
}



DimFig <- function(object, file_name, group.by = NULL, split.by = NULL, cols = NULL, plot.scale = 0.9, pt.size = "auto", label.size = 4, split.order = "auto", legend.position = "right", output_dir = outputdir, output.ggplot = FALSE, ...){
    p1 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = legend.position, label = TRUE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)
    p2 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = "none", label = TRUE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)
    p3 <- DimFig.default(object = object, file_name = file_name, group.by = group.by, split.by = split.by, cols = cols, plot.scale = plot.scale, legend.position = "none", label = FALSE, pt.size = pt.size, label.size = label.size, split.order = split.order, output_dir = output_dir, output.ggplot = TRUE, ...)


    plist <- lapply(1:length(group.by), function(x) list(p1[[x]], p2[[x]], p3[[x]]))

    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist %>% do.call(c, .))
    } else {

        nlevel <- ifelse_check(is.null(split.by), 1, ifelse_check("auto" %in% split.order, auto_levels(object = object, split.by = split.by),split.order) %>% length())
        nrow <- 3
        

        max_cls <- sapply(group.by, function(x) unique(object@meta.data[, x]) %>% length()) %>% max()
        pdf_width <- ifelse_check(legend.position == "none", 10 * plot.scale * nlevel, 10 * plot.scale * nlevel + ceiling(0.15 * length(max_cls)))

        for (xx in 1:length(group.by)){
            jpeg(paste0(output_dir, file_name, paste0("_", group.by[[xx]]), ifelse(is.null(split.by), "", paste0("_", split.by)), ".jpeg"), width = pdf_width, height = 10 * plot.scale * 3, units = "in", res = 300)
            plot_grid(plotlist = plist[[xx]], nrow = 3, ncol = 1) %>% print() %>% print()
            dev.off()
        }
    }
}



DimFig.default <- function(object, file_name, group.by = NULL, split.by = NULL, cols = NULL, plot.scale = 0.9, legend.position = "right", label = TRUE, pt.size = "auto", label.size = 4, split.order = "auto", output_dir = outputdir, output.ggplot = FALSE, ...){

    extra_arguments <- list(...)


    ##Set the point size
    if (pt.size == "auto"){
        point.size <- round(10/sqrt(ncol(object)), digits = 1) %>% MinMax(., min = 0.1, max = 1)
    } else {
        point.size <- pt.size
    }


    ##Set the levels of factors in the supervised mode
    if (!is.null(split.by)){
        ##Get the supervised order the split.by
        new_levels <- ifelse_check("auto" %in% split.order, 
                            auto_levels(object = object, split.by = split.by),
                            split.order)
        object@meta.data[, split.by] <- factor(as.character(object@meta.data[, split.by]), levels = new_levels) 
    }


    plot_args <- c(list(object = object, group.by = group.by, cols = cols, pt.size = point.size, combine = FALSE, split.by = split.by, label = label, label.size = label.size), extra_arguments[names(extra_arguments) %in% c("shape.by", "reduction", "cells", "repel", "cells.highlight" , "cols.highlight", "sizes.highlight", "na.value")])
    plist <- do.call(DimPlot, plot_args)


    ##Polish the figures by changing the theme
    plist <- lapply(plist, function(p) {
        p <- p + 
            #coord_equal(ratio = 1) + 
            #    theme_classic() + 
                theme(legend.position = legend.position)#,
            #        line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
            #        axis.text.x=element_blank(),axis.text.y=element_blank(), 
            #        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))+ 
            #    labs(title =  text_wrapper(x, width=20))
        p
        })


    ##Output the plot (either a plot object or print it)
    if (output.ggplot){
        return(plist)
    } else {
        nlevel <- ifelse(is.null(split.by), 1, length(new_levels))
        #ncol <- ifelse_check(is.null(ncol), nlevel, ncol)
        nrow <- length(plist)
        

        max_cls <- sapply(group.by, function(x) unique(object@meta.data[, x]) %>% length()) %>% max()
        pdf_width <- ifelse_check(legend.position == "none", 10 * plot.scale * nlevel, 10 * plot.scale * nlevel + ceiling(0.15 * length(max_cls)))

        jpeg(paste0(output_dir, file_name, ifelse(is.null(group.by), "", paste0("_", paste(group.by, collapse = "-"))), ifelse(is.null(split.by), "", paste0("_", split.by)), ifelse(label, "_labeled", ""), ".jpeg"), width = pdf_width, height = 10 * plot.scale * nrow, units = "in", res = 300)
        plot_grid(plotlist = plist, nrow = nrow, ncol = 1) %>% print() %>% print()
        dev.off()
    }
}


DotFig <- function(object, assay = "RNA", features, cols = c("white", "red"), dot.scale = 3, dot.min = 0, group.by = "hres", file_name, output_dir = outputdir) {
    p <- DotPlot(object = object, assay = assay, features = features, cols = cols, dot.scale = dot.scale,dot.min = dot.min, group.by = group.by) + 
        RotatedAxis() + 
        theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))


    ncls <- unique(object@meta.data[, group.by]) %>% length()
    plot_heights <- ceiling(0.2 * ncls + 2) %>% MinMax(., min = 5, max = 20)
    plot_width <- ceiling(0.6 * length(features) + 2) %>% MinMax(., min = 5, max = 20)
    jpeg(paste0(output_dir, file_name, "_markers.jpeg"), width = 4.3, height = 6, units = "in", res = 300)
    print(p)
    dev.off()
}



get_ident <- function(input_ident, ident_col = NULL, file_path, sample_name, label_names = NULL){
    if (grepl("Seurat", class(input_ident), ignore.case = TRUE)){
        tem_ident <- setNames(as.character(input_ident@meta.data[, ident_col]), colnames(input_ident))
        cells <- colnames(input_ident)
    } else {
        tem_ident <- as.character(input_ident)
        cells <- names(input_ident)
    }


    tem_file <- read.table(file_path, sep="\t", stringsAsFactors=FALSE, header=TRUE)

    tem_file$cluster <- as.character(tem_file$cluster) #change the cluster numbers to characters
    sub_file <- tem_file[tem_file$sample == sample_name, ]

    if (is.null(label_names)){
        ident_names <- sub_file$label; names(ident_names) <- sub_file$cluster
        output_ident <- ident_names[tem_ident]; names(output_ident) <- cells
    } else{
        ident_names <- setNames(sub_file[, label_names], sub_file$cluster)
        output_ident <- setNames(ident_names[tem_ident], cells)
    }
    return(output_ident)
}


LinearRegress <- function (input_data, latent_data, vars_regress, genes_regress = NULL,display_progress = TRUE, num.cores = NULL){
    input_data <- as.matrix(input_data) #not implemented for the "Matrix" object, need to update
    genes_regress <- ifelse_check(is.null(genes_regress), rownames(input_data), genes_regress)
    genes_regress <- intersect(x = genes_regress, y = rownames(input_data))
    latent_data <- latent_data[colnames(input_data),vars_regress, drop=FALSE]
    bin.size <- 100
    bin.ind <- ceiling(c(1:length(x = genes_regress))/bin.size)
    max.bin <- max(bin.ind)
    if (display_progress) {
        message(paste("Regressing out:", paste(vars_regress,collapse = ", ")))
        pb <- txtProgressBar(min = 0, max = max.bin, style = 3,file = stderr())
    }
    data.resid <- c()
    
    print("prepare for parallel calculation")
    num.cores <- ifelse_check(is.null(num.cores), detectCores()/2, num.cores)
    cl <- parallel::makeCluster(num.cores)
    doSNOW::registerDoSNOW(cl)
    opts <- list()
    if (display_progress) {
        progress <- function(n) setTxtProgressBar(pb, n)
        opts <- list(progress = progress)
        time_elapsed <- Sys.time()
    }
    
    #print("Set regression variables")
    reg.mat.colnames <- c(colnames(x = latent_data), "GENE")
    fmla_str = paste0("GENE ", " ~ ", paste(vars_regress,collapse = "+"))  #linear regression formula
    regression.mat <- cbind(latent_data, input_data[1, ])
    colnames(regression.mat) <- reg.mat.colnames
    qr = lm(as.formula(fmla_str), data = regression.mat,qr = TRUE)$qr
    rm(regression.mat)

    print("Do parallel regression")
    data.resid <- foreach::foreach(i = 1:max.bin, .combine="c",.options.snow = opts) %dopar% {
            genes.bin.regress <- rownames(x = input_data)[bin.ind ==i]
            gene.expr <- as.matrix(x = input_data[genes.bin.regress, , drop = FALSE])
            empty_char <- "" # character(length = dim(gene.expr)[1])
            new.data <- sapply(X = genes.bin.regress, FUN = function(x) {
                resid <- qr.resid(qr, gene.expr[x, ])
                if (!is.list(resid)) {
                  resid <- list(resid = resid, mode = empty_char)
                }
                return(resid)
            })          
            #print(length(x = new.data))
            new.data.resid <- new.data[seq.int(from = 1, to = length(x = new.data),by = 2)]
            #print(head(new.data.resid, 5))
            new.data.resid <- matrix(unlist(new.data.resid), nrow = length(new.data.resid[[1]]))
            colnames(x = new.data.resid) <- genes.bin.regress
            new.data.mode <- unlist(x = new.data[seq.int(from = 2,to = length(x = new.data), by = 2)])
            names(x = new.data.mode) <- genes.bin.regress
            new.data <- list(resid = new.data.resid, mode = new.data.mode)
            return(new.data)
    }
    print(length(data.resid))
    if (display_progress) {
        time_elapsed <- Sys.time() - time_elapsed
        cat(paste("\nTime Elapsed: ", time_elapsed, units(time_elapsed)))
        close(pb)
    }
    stopCluster(cl)
    modes <- unlist(x = data.resid[seq.int(from = 2, to = length(x = data.resid),by = 2)])
    modes <- modes[modes == "scale"]
    names(x = modes) <- gsub(pattern = "mode.", replacement = "", x = names(x = modes), fixed = TRUE)
    out_matrix <- data.resid[seq.int(from = 1, to = length(x = data.resid),by = 2)]
    out_matrix2 <- as.matrix(x = as.data.frame(x = out_matrix))
    out_matrix2 <- t(x = out_matrix2)
    print(dim(out_matrix2))
    if (length(x = modes)) {
        message("The following genes failed with glm.nb, and fell back to scale(log(y+1))\n\t",
            paste(names(x = modes), collapse = ", "))
    }
    rownames(out_matrix2) <- genes_regress
    colnames(out_matrix2) <- colnames(input_data)
    suppressWarnings(expr = gc(verbose = FALSE))
    return(out_matrix2)
}



#This scirpt wrap the long text to several lines, designed for the title of ggplot
text_wrapper <- function(x, width=15) {
  y <- sapply(x, function(x) paste(strsplit(x, "_", fixed=TRUE)[[1]], collapse=" "))
  return(paste(strwrap(y, width=width), collapse = "\n"))
}




FastInte <- function(object, split.by, nfeatures = 2000, inte.dims = 20, cls.dims = 25) {
    ## Do the integration for the two species.
    seu_list <- SplitObject(object, split.by = split.by) %>%
                    lapply(., function(x) FindVariableFeatures(x, nfeatures = nfeatures))
    hvg_use <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = nfeatures)

    ## Do the CCA integration
    data.anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:inte.dims, assay = NULL, anchor.features = hvg_use)
    data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:inte.dims)
    DefaultAssay(data.integrated) <- "integrated"
    data.integrated <- ScaleData(data.integrated, verbose = FALSE) %>%
                        RunPCA(., npcs = 50, verbose = FALSE) %>%
                        RunUMAP(., dims = 1:cls.dims) %>%
                        FindNeighbors(., dims = 1:cls.dims) %>%
                        FindClusters(., resolution = 1.2)
    return(data.integrated)
}



label_transfer <- function(object, reduction = "pca", dims.use = 40, transfer_cols = "label", k = 20){
    meta_use <- object@meta.data
    data.use <- object[[reduction]]@cell.embeddings[, 1:dims.use]

    
    ##A function to get the neighoring information
    get_most <- function(x){
        return(names(sort(table(x), decreasing = TRUE))[1])
    }


    ##Assume the NA values are the inquiry cell
    anno_list <- list()

    for (anno in transfer_cols){
        ref_cells <- rownames(meta_use)[!is.na(meta_use[, anno])]
        inquiry_cells <- setdiff(rownames(meta_use), ref_cells)


        ## Get the KNN cells for all the inquiry cells
        knn_cells <- FNN::get.knnx(data = data.use[ref_cells, ,drop = FALSE], query = data.use[inquiry_cells, ,drop = FALSE], k = k) %>%
                    .$nn.index


        ref_anno <- meta_use[ref_cells, anno] %>% setNames(., ref_cells)

        ## Do the annotation transfer
        if (class(ref_anno) %in% "character"){
            new_label <- apply(knn_cells, 1, function(x) get_most(ref_anno[x])) %>% 
                            setNames(., inquiry_cells) %>%
                            c(., ref_anno)
        } else if (class(ref_anno) %in% "numeric") {
            new_label <- apply(knn_cells, 1, function(x) median(ref_anno[x])) %>% 
                            setNames(., inquiry_cells) %>%
                            c(., ref_anno)
        } else {
            stop(paste0("column ", column, " has unsupported object type"))
        }

        anno_list[[anno]] <- new_label[rownames(meta_use)]
    }
    

    new_meta <- anno_list %>%  
                    as.data.frame(., stringsAsFactors = FALSE) %>% 
                    setNames(., paste0(names(anno_list), "_new"))
    object@meta.data <- cbind(object@meta.data[, setdiff(colnames(object@meta.data), colnames(new_meta))], new_meta)
    return(object)
}











###############################################################################################

## Dot Plot (modified version)

###############################################################################################
PointPlot <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA, shape = 16, cluster_order = NULL, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset")){
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    if (!is.null(x = group.by)) {
        Idents(object) <- group.by
    }

    data.features$id <- Idents(object = object)
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot ==
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min,
                  max = col.max)
            }
            else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = features))
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = as.character(x = data.plot$id),
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((",
                paste(sort(x = levels(x = object), decreasing = TRUE),
                  collapse = "|"), ")_)"), replacement = "",
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }

    ## Add a column for split
    data.plot$species <- extract_field(as.character(data.plot$id), -1, "_")
    data.plot$cluster <- extract_field(as.character(data.plot$id), "rm_end", "_")


    ## Define the order the species & cluster
    if (!is.null(species_order)){
        data.plot$species <- factor(data.plot$species, levels = species_order)
    } 
    if (!is.null(cluster_order)){
        data.plot$cluster <- factor(data.plot$cluster, levels = cluster_order)
    }
    

    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
        y = "cluster")) + geom_point(mapping = aes_string(size = "pct.exp",
        color = color.by), shape = shape) + scale.func(range = c(0, dot.scale),
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
        labs(x = "Features", y = ifelse(test = is.null(x = split.by),
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    } else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}
environment(PointPlot) <- environment(DotPlot)


## PointPlot2 is designed for stacked plots
##Change 1: data.plot$features.plot <- factor(x = data.plot$features.plot, levels = rev(x = features)), remove rev
##Change 2: ggplot, y axis is "species"
PointPlot2 <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA, shape = 16, cluster_order = NULL, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset")){
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    if (!is.null(x = group.by)) {
        Idents(object) <- group.by
    }

    data.features$id <- Idents(object = object)
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        names(x = cols) <- unique(x = splits)
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot ==
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min,
                  max = col.max)
            }
            else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features)
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = as.character(x = data.plot$id),
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((",
                paste(sort(x = levels(x = object), decreasing = TRUE),
                  collapse = "|"), ")_)"), replacement = "",
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }

    ## Add a column for split
    data.plot$species <- extract_field(as.character(data.plot$id), -1, "_")
    data.plot$cluster <- extract_field(as.character(data.plot$id), "rm_end", "_")


    ## Define the order the species & cluster
    if (!is.null(species_order)){
        data.plot$species <- factor(data.plot$species, levels = species_order)
    } 
    if (!is.null(cluster_order)){
        data.plot$cluster <- factor(data.plot$cluster, levels = cluster_order)
    }
    

    plot <- ggplot(data = data.plot, mapping = aes_string(x = "species",
        y = "cluster")) + geom_point(mapping = aes_string(size = "pct.exp",
        color = color.by), shape = shape) + scale.func(range = c(0, dot.scale),
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
        labs(x = "Features", y = ifelse(test = is.null(x = split.by),
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    } else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}
environment(PointPlot2) <- environment(DotPlot)



SpeciesDot <- function(object, features, cols = color_list$sp, dot.scale = 5, dot.min = 0.05, group.by = "hres", split.by = "species", shape = 16, layout_type = c("v", "h", "hflip")[1], panel.space = 0, cluster_order = NULL, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), scale.by = "radius") {

    if (!is.null(cluster_order) && layout_type %in% c("v", "h")){
        cluster_order <- rev(cluster_order)
    }

    if (layout_type %in% c("h")){
        features <- rev(features)
    }

    p <- PointPlot(object, assay = "RNA", features = features, cols = cols, dot.scale = dot.scale, dot.min = dot.min, group.by = group.by, split.by = split.by, shape = shape, cluster_order = cluster_order, species_order = species_order, scale.by = scale.by)

    if (layout_type == "v"){
        p <- p + facet_wrap(vars(species), nrow = length(species_order), ncol = 1, scales = "free_y") + 
                    scale_x_discrete(position = "top") + 
                    theme(axis.text.x=element_text(size = 8, hjust = 0, angle = 45), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(panel.space, "in"))
    } else if (layout_type == "hflip") {
        p <- p + coord_flip() + 
                    RotatedAxis() + 
                    facet_wrap(vars(species), nrow = 1, ncol = length(species_order)) + 
                    theme(axis.text.x=element_text(size = 8), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(panel.space, "in"))
    } else if (layout_type == "h") {
        p <- p + RotatedAxis() + 
                    facet_wrap(vars(species), nrow = length(species_order), ncol = 1) + 
                    theme(axis.text.x=element_text(size = 8), strip.background = element_blank(), strip.text = element_blank(), panel.spacing = unit(panel.space, "in"))

    } else {
        stop("layout_type should be v, h or hflip")
    }
    return(p)
}
##, scales = "free_x"


SplitGenePointPlot <- function (object, assay = NULL, features, cols = c("lightgrey", "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, group.by = NULL, split.by = NULL, scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA, shape = 16, cluster_order = NULL, species_order = c("Human", "Chimpanzee", "Rhesus", "Marmoset")){
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    scale.func <- switch(EXPR = scale.by, size = scale_size,
                        radius = scale_radius, 
                        stop("'scale.by' must be either 'size' or 'radius'"))
    data.features <- FetchData(object = object, vars = features)
    if (!is.null(x = group.by)) {
        Idents(object) <- group.by
    }

    data.features$id <- Idents(object = object)
    if (!is.factor(x = data.features$id)) {
        data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    if (!is.null(x = split.by)) {
        splits <- object[[split.by, drop = TRUE]]
        if (length(x = unique(x = splits)) > length(x = cols)) {
            stop("Not enought colors for the number of groups")
        }
        cols <- cols[1:length(x = unique(x = splits))]
        if (is.null(names(cols))){
            names(x = cols) <- unique(x = splits)
        } 
        data.features$id <- paste(data.features$id, splits, sep = "_")
        unique.splits <- unique(x = splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
            "_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
            1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
            return(mean(x = expm1(x = x)))
        })
        pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove,
            threshold = 0)
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
        data.use <- as.data.frame(x = data.plot[[x]])
        data.use$features.plot <- rownames(x = data.use)
        data.use$id <- x
        return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
        data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
        scale <- FALSE
        warning("Only one identity present, the expression values will be not scaled.")
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
        FUN = function(x) {
            data.use <- data.plot[data.plot$features.plot ==
                x, "avg.exp"]
            if (scale) {
                data.use <- scale(x = data.use)
                data.use <- MinMax(data = data.use, min = col.min,
                  max = col.max)
            }
            else {
                data.use <- log(x = data.use)
            }
            return(data.use)
        })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    if (!is.null(x = split.by)) {
        avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled,
            breaks = 20))
    }
    data.plot$avg.exp.scaled <- avg.exp.scaled

    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    if (!is.null(x = split.by)) {
        splits.use <- vapply(X = as.character(x = data.plot$id),
            FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((",
                paste(sort(x = levels(x = object), decreasing = TRUE),
                  collapse = "|"), ")_)"), replacement = "",
            USE.NAMES = FALSE)
        data.plot$colors <- mapply(FUN = function(color, value) {
            return(colorRampPalette(colors = c("grey", color))(20)[value])
        }, color = cols[splits.use], value = avg.exp.scaled)
    }
    color.by <- ifelse(test = is.null(x = split.by), yes = "avg.exp.scaled",
        no = "colors")
    if (!is.na(x = scale.min)) {
        data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }

    ## Add a column for split
    data.plot$species <- extract_field(as.character(data.plot$id), -1, "_")
    data.plot$cluster <- extract_field(as.character(data.plot$id), "rm_end", "_")


    ## Set the new features
    all_sps <- levels(as.factor(object@meta.data[, split.by]))
    all_sps <- species_order[species_order %in% all_sps]
    new_features <- paste0(rep(all_sps, lengt.out = length(features) * length(all_sps)), "_", rep(features, each = length(all_sps)))
    print(new_features)
    print(levels(as.factor(data.plot$features.plot)))

    data.plot$features.plot <- factor(x = paste0(data.plot$species, "_", as.character(data.plot$features.plot)),
        levels = rev(x = new_features))

    ## Define the order the cluster
    if (!is.null(cluster_order)){
        data.plot$cluster <- factor(data.plot$cluster, levels = cluster_order)
    }
    

    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot",
        y = "cluster")) + geom_point(mapping = aes_string(size = "pct.exp",
        color = color.by), shape = shape) + scale.func(range = c(0, dot.scale),
        limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) +
        labs(x = "Features", y = ifelse(test = is.null(x = split.by),
            yes = "Identity", no = "Split Identity")) + theme_cowplot()
    if (!is.null(x = split.by)) {
        plot <- plot + scale_color_identity()
    } else if (length(x = cols) == 1) {
        plot <- plot + scale_color_distiller(palette = cols)
    } else {
        plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
    }
    if (is.null(x = split.by)) {
        plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    }
    return(plot)
}
environment(SplitGenePointPlot) <- environment(DotPlot)




