library(Seurat)
library(dplyr)
inputdir <- "./load_files/"
outputdir <- "./report/"


## Cell type average expression
exn <- readRDS(file = "../overview/load_files/ExN_data_08312022.rds")
lexn <- subset(exn, cbnage %in% c("E93", "E110"))
rm(exn)


## Get raw regions
lexn@meta.data$region <- strsplit(lexn@meta.data$cell_origin, "_", fixed = TRUE) %>%
						sapply(., function(x) x[2])
lexn@meta.data$subtype2 <- ifelse(lexn@meta.data$subtype %in% "ExN SOX5 ID2", "ExN SOX5 PALMD", lexn@meta.data$subtype)



## Mainly use E93 cells, calculating average gene expression
subexn <- subset(lexn, cbnage %in% "E93")


## Furhter subset the data to have equal # cells per subtype/region pair
ncell <- table(subexn@meta.data$subtype2, subexn@meta.data$region)


set.seed(42)
cells <- list()
for (cls in rownames(ncell)){
	if (sum(ncell[cls, ] >= 40) >= 4){
		msize <- 40
		if (sum(ncell[cls, ] >= 100) >= 4) {
			msize <- 100
		}
	} else {
		msize <- 100000
	}

	reg_kp <- colnames(ncell)[ncell[cls, ] >= msize]
	if (length(reg_kp) >= 3){
		for (reg in reg_kp){
			subc <- colnames(subexn)[subexn@meta.data$subtype2 %in% cls & subexn@meta.data$region %in% reg]
			if (length(subc) > msize){
				subc <- sample(subc, msize)
				cells <- c(cells, subc)
			}
		}
	}
}


slimexn <- subexn[, cells]
ncell_slim <- table(slimexn@meta.data$subtype2, slimexn@meta.data$region)

saveRDS(slimexn, file = paste0("./load_files/Hier-region.slimExN.rds"))

## Calculate Average expression
slimexn@meta.data$tmpcls <- paste0(slimexn@meta.data$subtype2, "|", slimexn@meta.data$region)
Idents(slimexn) <- "tmpcls"
avg <- log(AverageExpression(slimexn, assay = "RNA")$RNA + 1)
save(avg, ncell_slim, file = paste0("./load_files/Hier-region.subset.avgs.Rdata"))












