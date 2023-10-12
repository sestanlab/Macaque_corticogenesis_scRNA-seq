library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)

inputdir <- "./load_files/"
outputdir <- "./report/"



load(file = paste0(inputdir, "Pseudobulk_by_region_tRG_oRG.Rdata")) ## avgs, pmeta, rgc_cbn, 



## Get smooth gene expression along the trajectory
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
							mutate(region = sapply(strsplit(avgcls, "|", fixed = TRUE), "[", 1)) %>%
							mutate(pseudotime = as.numeric(as.character(sapply(strsplit(avgcls, "|", fixed = TRUE), "[", 2)))) %>%
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




trg_pmeta <- pmeta[pmeta$lineage %in% c("shared", "tRG"), ]
trg_smt <- SmoothExprAcrossRegions(avgs = as.matrix(avgs), meta = trg_pmeta, reg_ord = c("FC", "MSC", "TC", "OcC"), span = 0.75)
trg_seu <- subset(rgc_cbn, lineage %in% c("shared", "tRG"))
trg_meta <- PrepSMmeta(avgs = trg_smt, object = trg_seu, ptime_col = "pseudotime", group.by = "cluster", ncells = 25)
save(trg_smt, trg_meta, file = paste0(inputdir, "Smooth_by_region_tRG.Rdata"))



org_pmeta <- pmeta[pmeta$lineage %in% c("oRG"), ]
org_smt <- SmoothExprAcrossRegions(avgs = as.matrix(avgs), meta = org_pmeta, reg_ord = c("FC", "MSC", "TC", "OcC"), span = 0.75)
org_seu <- subset(rgc_cbn, lineage %in% c("oRG"))
org_meta <- PrepSMmeta(avgs = org_smt, object = org_seu, ptime_col = "pseudotime", group.by = "cluster", ncells = 25)
save(org_smt, org_meta, file = paste0(inputdir, "Smooth_by_region_oRG.Rdata"))




