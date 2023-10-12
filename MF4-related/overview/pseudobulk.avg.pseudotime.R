source("../scripts/nhpf.fun.R")


# Load data
rgc <- readRDS(file = paste0(inputdir, "RGC_filtered_seu_04112022.rds"))
rgc <- rgc[, !rgc@meta.data$cluster %in% "RGC FABP7 PMP22"]
Idents(rgc) <- factor(rgc@meta.data$cluster, levels = c("NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "oRG HOPX TNC", "oRG HOPX APOE", "tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))



##--------------------------------------------------------------------------------
## Manually remove certain small clusters for more accurate estimation
rm_cls <- c("Insula|Ependymal", "MSC|Ependymal", "MSC|NEP RSPO3 DIRAS3", "MSC|NEP RSPO3 TEX15", "TC|Ependymal", "TC|NEP RSPO3 DIRAS3", "TC|NEP RSPO3 TEX15", "TC|vRG HMGA2 CCND1", "TC|tRG CRYAB FNDC1")
rgc@meta.data$tmpcls <- paste0(rgc@meta.data$lobe, "|", rgc@meta.data$cluster)
rgc <- rgc[, !rgc@meta.data$tmpcls %in% rm_cls]



##--------------------------------------------------------------------------------
## Assign cells to bins
min_pt <- 0
max_pt <- quantile(rgc@meta.data$pseudotime, 0.999, na.rm = TRUE)
nPoints <- 40
half_inter <- (max_pt - min_pt)/(2 * (nPoints - 1))
cut_bks <- c(seq(min_pt, max_pt, length.out = nPoints) - half_inter, Inf)
rgc@meta.data$bin <- cut(rgc@meta.data$pseudotime, breaks = cut_bks) %>% 
                    as.numeric()
rgc@meta.data$bin[rgc@meta.data$bin == 1] <- 2
rgc@meta.data$bin <- rgc@meta.data$bin - 1


## Split to three lineages
rgc@meta.data$lineage <- "shared"
rgc@meta.data$lineage[rgc@meta.data$cluster %in% c("oRG HOPX TNC", "oRG HOPX APOE")] <- "oRG"
rgc@meta.data$lineage[rgc@meta.data$cluster %in% c("tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal")] <- "tRG"

all_regs <- c("FC","MSC","TC","OcC")



rgc_sh <- subset(rgc, lineage == "shared")
##cls_size <- table(rgc_sub@meta.data$bin, rgc_sub@meta.data$lobe) %>% as.matrix()
## 100 cells minimum
rgc_sh@meta.data$bin <- as.integer(rgc_sh@meta.data$bin)
set.seed(42)
rm_cells <- colnames(rgc_sh)[rgc_sh@meta.data$bin %in% c(1:4) & rgc_sh@meta.data$cluster == "vRG HMGA2 CCND1" & rgc_sh@meta.data$lobe %in% c("FC", "OcC")]
rm_cells <- sample(rm_cells, length(rm_cells) - 200)
rgc_sh <- rgc_sh[, setdiff(colnames(rgc_sh), rm_cells)]

sh_list <- list(FC = 1:20, 
				MSC = 3:20,
				OcC = 1:20,
				TC = 12:20)
sh_cls <- rep(names(sh_list), sapply(sh_list, length)) %>%
			paste0(., "|", unlist(sh_list)) %>%
			paste0("shared", "|", .)
rgc_sh@meta.data$avgcls <- paste0(rgc_sh@meta.data$lineage, "|", rgc_sh@meta.data$lobe, "|", rgc_sh@meta.data$bin)
rgc_sh <- subset(rgc_sh, avgcls %in% sh_cls)



rgc_org <- subset(rgc, lineage == "oRG")
rgc_org@meta.data$bin <- as.character(rgc_org@meta.data$bin)
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("28", "29")] <- "28"
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("30", "31", "32")] <- "30"
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("38", "39")] <- "37"
rgc_org@meta.data$bin <- as.integer(rgc_org@meta.data$bin)
rgc_org@meta.data$avgcls <- paste0(rgc_org@meta.data$lineage, "|", rgc_org@meta.data$lobe, "|", rgc_org@meta.data$bin)
org_list <- list(FC = c(18:28, 30, 33:37), 
				MSC = c(18:28, 30, 33:37),
				OcC = c(18:28, 30, 33:37),
				TC = c(18:28, 30, 33:37))
org_cls <- rep(names(org_list), sapply(org_list, length)) %>%
			paste0(., "|", unlist(org_list)) %>%
			paste0("oRG", "|", .)
rgc_org <- subset(rgc_org, avgcls %in% org_cls)



rgc_trg <- subset(rgc, lineage == "tRG")
rgc_trg@meta.data$bin <- as.character(rgc_trg@meta.data$bin)
rgc_trg@meta.data$bin[rgc_trg@meta.data$bin %in% c("16", "17", "18")] <- "18"
rgc_trg@meta.data$bin[rgc_trg@meta.data$bin %in% c("28", "29", "30", "31")] <- "28"
rgc_trg@meta.data$bin[rgc_trg@meta.data$bin %in% c("32", "33", "34", "35")] <- "29"
rgc_trg@meta.data$bin[rgc_trg@meta.data$bin %in% c("36", "37", "38", "39")] <- "30"
rgc_trg@meta.data$bin <- as.integer(rgc_trg@meta.data$bin)
trg_list <- list(FC = 18:30, 
				MSC = 18:30,
				OcC = 18:30,
				TC = 18:30)
trg_cls <- rep(names(trg_list), sapply(trg_list, length)) %>%
			paste0(., "|", unlist(trg_list)) %>%
			paste0("tRG", "|", .)
rgc_trg@meta.data$avgcls <- paste0(rgc_trg@meta.data$lineage, "|", rgc_trg@meta.data$lobe, "|", rgc_trg@meta.data$bin)
rgc_trg <- subset(rgc_trg, avgcls %in% trg_cls)



rgc_cbn <- Reduce(merge, list(rgc_sh, rgc_trg, rgc_org))
Idents(rgc_cbn) <- "avgcls"
avgs <- log(AverageExpression(rgc_cbn)$RNA + 1)


## Get the meta data for bins
rgc_cbn@meta.data$lin_bin <- paste0(rgc_cbn@meta.data$lineage, "|", extract_field(rgc_cbn@meta.data$avgcls, 3, "|"))
times <- aggregate(pseudotime ~ lin_bin, rgc_cbn@meta.data, mean) %>%
			column_to_rownames("lin_bin")

get_most <- function(x) {
    table(x) %>% sort() %>% rev() %>% .[1] %>% names()
}
labels <- aggregate(cluster ~ lin_bin, rgc_cbn@meta.data, get_most) %>%
			column_to_rownames("lin_bin")


pmeta <- data.frame(avgcls = unique(rgc_cbn@meta.data$avgcls), 
				stringsAsFactors = FALSE) %>%
				mutate(lineage = extract_field(avgcls, 1, "|")) %>%
				mutate(lobe = extract_field(avgcls, 2, "|")) %>%
				mutate(bin = as.numeric(extract_field(avgcls, 3, "|"))) %>%
				mutate(lin_bin = paste0(lineage, "|", bin))
pmeta$pseudotime <- times[pmeta$lin_bin, "pseudotime"]
pmeta$cluster <- labels[pmeta$lin_bin, "cluster"]
save(avgs, pmeta, rgc_cbn, file = paste0(inputdir, "Pseudobulk_by_region.Rdata"))





##--------------------------------------------------------------------------------
## Combine regions
source("./ptime.fun.R")
load(file = paste0(inputdir, "Pseudobulk_by_region.Rdata"))
# avgs, pmeta, rgc_cbn
rm_cells <- rownames(rgc_cbn@meta.data)[rgc_cbn@meta.data$cluster == "tRG CRYAB MEST" & rgc_cbn@meta.data$cbnage %in% c("E93", "E110")]
rgc_cbn <- rgc_cbn[, setdiff(colnames(rgc_cbn), rm_cells)]


Idents(rgc_cbn) <- "lin_bin"
avgs <- log(AverageExpression(rgc_cbn)$RNA + 1)


cbn_meta <- pmeta %>%
			group_by(lin_bin) %>%
			summarize(pseudotime = mean(pseudotime), cluster = unique(cluster), lineage = unique(lineage)) %>%
			column_to_rownames("lin_bin")


## Get the fitted values
avg_smt <- lapply(c("shared", "oRG", "tRG"), function(lin) {
	avg_sub <- avgs[, grep(lin, colnames(avgs), value = TRUE)]
	sub_meta <- cbn_meta[colnames(avg_sub), ]
	sub_meta <- sub_meta[order(sub_meta$pseudotime), ,drop = FALSE]


	## The predicted df rownames will be the column names of the predicted data
	pt_min <- round(min(sub_meta$pseudotime), digits = 6)
	if (pt_min < min(sub_meta$pseudotime)){
		pt_min <- pt_min + 5e-6
	} 
	pt_max <- round(max(sub_meta$pseudotime), digits = 6)
	if (pt_max > max(sub_meta$pseudotime)){
		pt_max <- pt_max - 5e-6
	} 


	ntime <- seq(pt_min, pt_max, length.out = switch(lin, shared = 100, oRG = 70, tRG = 70)) %>%
				sapply(., function(x) round(x, digits = 6))
	predict_df <- data.frame(row.names = paste0(lin, "|", ntime), 
			pseudotime = ntime)
	ptime <- setNames(sub_meta$pseudotime, rownames(sub_meta))



	sub_smt <- SmoothBinAvg_combine(data = avg_sub, ptime = ptime, predict_df = predict_df, loess_span = 0.6)
	return(sub_smt)
	}) %>%
		do.call(cbind, .)

## Get the cluster for the specified ptime
meta_smt <- data.frame(avgcls = colnames(avg_smt), 
						stringsAsFactors = FALSE) %>%
						mutate(lineage = extract_field(avgcls, 1, "|")) %>%
						mutate(pseudotime = as.numeric(as.character(extract_field(avgcls, 2, "|")))) %>%
						column_to_rownames("avgcls")
meta_smt$cluster <- sapply(colnames(avg_smt), function(xx) {
	lin <- strsplit(xx, "|", fixed = TRUE)[[1]][1]
	time <- as.numeric(strsplit(xx, "|", fixed = TRUE)[[1]][2])

	lin_meta <- rgc_cbn@meta.data %>% 
					subset(lineage == lin)

	high_cells <- lin_meta[lin_meta$pseudotime > time, ] %>%
					arrange(pseudotime) %>%
					.$cluster %>% .[1:25]
	low_cells <- lin_meta[lin_meta$pseudotime < time, ] %>%
					arrange(desc(pseudotime)) %>%
					.$cluster %>% .[1:25]
	cls <- c(high_cells, low_cells) %>%
				table() %>% sort() %>% 
				rev() %>% .[1] %>% names()
	cls
	})



res <- list(raw = list(avg = avgs, meta = cbn_meta),
				fit = list(avg = avg_smt, meta = meta_smt)
				)
saveRDS(res, file = paste0(inputdir, "Pseudobulk_cbn2.rds"))


