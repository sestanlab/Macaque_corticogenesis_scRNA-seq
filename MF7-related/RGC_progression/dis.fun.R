library(foreach)
library(BiocParallel)
library(AUCell)
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
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




SubsampleData.count.equalUMI <- function(counts, numi, seed.use = 42, nCores = 12){
    if (TRUE) {
    	count_data <- counts

        ## Do the subsampling
        cellbin.size <- 1000
        bin.ind <- ceiling(c(1:ncol(count_data))/cellbin.size)
        max.bin <- max(bin.ind)


        ## prepare the parallel computing
        print(Sys.time())
        cl = makeCluster(nCores, outfile="")
        doParallel::registerDoParallel(cl);
        tem_idx <- NULL


        ## Do the parallel subset
        set.seed(seed.use)
        new_count <- foreach(tem_idx = 1:max.bin, .combine = cbind, .packages = c("Matrix"), .export = c("ifelse_check")) %dopar% {
                    ctm <- count_data[, bin.ind == tem_idx]
                    rownames(ctm) <- NULL


                    ctm_list <- lapply(1:ncol(ctm), function(idx) {
                        x <- ctm[, idx]
                        i <- which(x != 0)
                        cc <- rep(i, times = x[i])
                        dd <- table(cc[sort(sample(1:length(cc), numi))])
                        new_i <- as.numeric(names(dd))
                        new_xx <- setNames(dd, NULL)
                        list(new_i, new_xx)
                        })


                    i_list <- lapply(ctm_list, function(x) x[[1]])
                    x_list <- lapply(ctm_list, function(x) x[[2]])
                    i_counts <- unlist(lapply(i_list, length), use.names = F)

                    ## Get the new dataset
                    new_ctm <- sparseMatrix(
                        i = unlist(i_list, use.names = FALSE),
                        j = rep(1:ncol(ctm), i_counts),
                        x = unlist(x_list, use.names = FALSE),
                        dims = dim(ctm),
                        dimnames = list(rownames(count_data), colnames(ctm))
                        )
                    return(new_ctm)
        }

        stopCluster(cl)
        print(paste0("Finish Subsampling"))
        print(Sys.time())

        #object[["RNA"]] <- CreateAssayObject(counts = new_count, min.cells = 0, min.features = 0)
        #DefaultAssay(object) <- "RNA"
        #object <- NormalizeData(object, normalization.method = "LogNormalize")
        new_count <- new_count[rownames(counts), ]
        return(new_count)
    }
}



ParallelAVG <- function(object, seed.use = 42, nreps = 100, ncells = 50, nCores) {
	set.seed(seed.use)
	cell_list <- replicate(nreps, sample(colnames(object), ncells), simplify = FALSE)

    print(Sys.time())
    cl = makeCluster(nCores, outfile="")
    doParallel::registerDoParallel(cl);

    data <- object$RNA@data

    ## Do the parallel subset
    set.seed(seed.use)
    pbavgs <- foreach(idx = 1:length(cell_list), .combine = cbind, .packages = c("Matrix"), .export = c("ifelse_check")) %dopar% {
    	subdata  <- data[, cell_list[[idx]]]
    	avg <- rowMeans(x = expm1(subdata))
    	return(avg)
    }

    colnames(pbavgs) <- paste0("n", 1:length(cell_list))

    stopCluster(cl)
    print(paste0("Finish Subsampling"))
    print(Sys.time())

    return(pbavgs)
}


ParallelAVG_Allcls <- function(object, group.by, nreps = 4, ncells = 50, nCores = 2, seed.use = 42) {
	all_cls <- levels(as.factor(object@meta.data[, group.by]))

	avgs <- lapply(all_cls, function(cls){
		obj <- object[, object@meta.data[, group.by] == cls]
		aa <- ParallelAVG(object = obj, seed.use = seed.use, nreps = nreps, ncells = ncells, nCores = nCores) %>%
				as.matrix()
		colnames(aa) <- paste0(cls, "_", 1:ncol(aa))
		return(aa)
		}) %>%
		do.call(cbind, .)
	return(avgs)

}


GetModuleScore_NEW <- function (assay.data, features, nbin = 24, ctrl = 100, k = FALSE, seed = 42, method = c("seurat","aucell")[2], input_dir = inputdir, file_name, output_dir = outputdir, rethreshold_list = NULL, cellbin.size = 8000, nCores = 4) {
     if (is.null(x = features)) {
        stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
        return(intersect(x = x, y = rownames(x = assay.data)))
    })
    cluster.length <- length(x = features) #number of feature list

    if (method == "aucell"){
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
                tem_rankings <- AUCell_buildRankings(tem_matrix, nCores=nCores, plotStats=FALSE) 
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



