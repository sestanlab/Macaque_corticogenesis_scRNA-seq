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


