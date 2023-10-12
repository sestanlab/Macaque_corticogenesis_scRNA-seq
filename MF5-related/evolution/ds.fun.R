## This script will subset all cells to the same number of UMIs
library(foreach)
library(parallel)
library(BiocParallel)
ifelse_check <- function(test, yes, no){if (test){return(yes)} else{ return(no) }}


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








