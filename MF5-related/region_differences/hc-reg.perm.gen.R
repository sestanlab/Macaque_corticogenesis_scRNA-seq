library(Seurat)
library(dplyr)
library(parallel)
library(Matrix)
inputdir <- "./load_files/"
outputdir <- "./report/"


slimexn <- readRDS(file = paste0("./load_files/Hier-region.slimExN.rds"))
hvg <- readRDS(file = paste0("./load_files/Hier-region-cluster.hvg.rds"))
slimexn <- slimexn[hvg, ]


CalcAvgFromData <- function(data, ident) {
	ident <- ident[colnames(data)]
	all.levs <- levels(as.factor(ident))
	avg.all <- list()
	for (ii in all.levs){
		tmp.cells <- which(ident == ii)
		data.tmp <- data[, tmp.cells]
		if (length(tmp.cells) == 1){
			avg.tmp <- expm1(data.tmp)
		} else if (length(tmp.cells) > 1){
			avg.tmp <- rowMeans(expm1(data.tmp))
		}

		avg.all[[ii]] <- avg.tmp
	}
	avg.final <- log(do.call(cbind, avg.all) + 1)
	return(avg.final)
}

raw_data <- slimexn$RNA@data
raw_ident <- setNames(paste0(slimexn@meta.data$subtype2, "|", slimexn@meta.data$region), colnames(slimexn))

npermutations <- 1000
perm_res <- mclapply(1:npermutations, function(idx) {
		print(idx)
		set.seed(idx)
		perm_data <- apply(as.matrix(raw_data), 1, sample) %>% t()
		colnames(perm_data) <- colnames(raw_data)
		perm_avg <- CalcAvgFromData(data = perm_data, ident = raw_ident)
		return(perm_avg)
		}, mc.cores = 12)
saveRDS(perm_res, file = paste0("./load_files/Hier-region.perm.avgs.rds"))





