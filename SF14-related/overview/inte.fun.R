Integratelist.seurat <- function(obj.list, hvg, file_name, input_dir = inputdir, inte.dims = 1:30, cluster.dims = 1:30, reference = NULL, do.cluster = TRUE) {
    if (length(obj.list) == 2){
        newseu <- merge(x = obj.list[[1]], y = obj.list[[2]])
    } else {
        newseu <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)])
    }


    inte.slim.file <- paste0(input_dir, file_name, ".slim.rds")
    if (!file.exists(inte.slim.file)){
        dg.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = inte.dims, assay = NULL, anchor.features = hvg, reference = reference)
        seuinte <- IntegrateData(anchorset = dg.anchors, dims = inte.dims)
        DefaultAssay(seuinte) <- "integrated"
        seuinte <- ScaleData(seuinte, verbose = FALSE) %>%
                                RunPCA(., npcs = 50, verbose = FALSE)
        
        newseu[["pca"]] <- CreateDimReducObject(embeddings = seuinte$pca@cell.embeddings[colnames(newseu), ], loadings = seuinte$pca@feature.loadings, stdev = seuinte$pca@stdev, key = "PC_", assay = "RNA")
        rm(seuinte)
        newseu <- RunUMAP(newseu, dims = cluster.dims, umap.method = "umap-learn", metric = "correlation")
        if (do.cluster){
            newseu <- FindNeighbors(newseu, dims = cluster.dims,  k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20)
        } else {
            newseu@meta.data$seurat_clusters <- "empty"
        }
        saveRDS(newseu, file = inte.slim.file)
    }else {
        newseu <- readRDS(file = inte.slim.file)
    }
    return(newseu)
}


## Select integration features
SelectHVG <- function(hvg.list, nfeatures = 2000) {
    var.features <- unname(obj = unlist(x = hvg.list)) %>%
                        table() %>%
                        sort(., decreasing = TRUE)

    tie.val <- var.features[min(nfeatures, length(x = var.features))]
    features <- names(x = var.features[which(x = var.features > tie.val)])
    if (length(x = features) > 0) {
        feature.ranks <- sapply(X = features, FUN = function(x) {
            ranks <- sapply(X = hvg.list, FUN = function(vf) {
                if (x %in% vf) {
                  return(which(x = x == vf))
                }
                return(NULL)
            })
            median(x = unlist(x = ranks))
        })
        features <- names(x = sort(x = feature.ranks))
    }
    features.tie <- var.features[which(x = var.features == tie.val)]
    tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
        ranks <- sapply(X = hvg.list, FUN = function(vf) {
            if (x %in% vf) {
                return(which(x = x == vf))
            }
            return(NULL)
        })
        median(x = unlist(x = ranks))
    })
    features <- c(features, names(x = head(x = sort(x = tie.ranks), 
        nfeatures - length(x = features))))
    return(features)
}



##---------------------------------------------------------------------------------------------
## Function to perform the harmony integration
library(harmony)
InteAllSp.harmony <- function(object, split.by, hvg, file_name, input_dir = inputdir, inte.dims = 30, theta = 2, lambda = 1, sigma = 0.1, do_cluster = FALSE) {
    inte.slim.file <- paste0(input_dir, file_name, ".harmony.rds")
    if (!file.exists(inte.slim.file)){
        ## Do the Harmony integration
        object <- ScaleData(object, split.by = split.by, do.center = FALSE, features = hvg)%>%
                    RunPCA(., features = hvg, verbose = FALSE) %>%
                    RunHarmony(., group.by.vars = split.by, lambda = lambda, theta = theta, dims.use = 1:inte.dims, sigma = sigma)
        object <- RunUMAP(object, dims = 1:inte.dims, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
        if (do_cluster){
        	object <- FindNeighbors(object, dims = 1:inte.dims, reduction = "harmony", k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20)
        }
        saveRDS(object, file = inte.slim.file)
    }else {
        object <- readRDS(file = inte.slim.file)
    }
    return(object)
}










