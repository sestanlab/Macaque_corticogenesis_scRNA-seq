## Combine Shared & oRG lineage to get smoother gene expression patterns
library(Seurat)
library(Matrix)
library(dplyr)
library(tibble)

inputdir <- "./load_files/"
outputdir <- "./report/"
source("./ptime.fun.v2.R")


# Load data
rgc <- readRDS(file = paste0(inputdir, "RGC_seu_for_DEG_analysis.rds"))



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



##------------------------------------------------------------
## shared lineage
rgc_sh <- subset(rgc, lineage == "shared")
## 100 cells minimum
rgc_sh@meta.data$bin <- as.integer(rgc_sh@meta.data$bin)
set.seed(42)
rm_cells <- colnames(rgc_sh)[rgc_sh@meta.data$bin %in% c(1:4) & rgc_sh@meta.data$cluster == "vRG HMGA2 CCND1" & rgc_sh@meta.data$lobe %in% c("FC", "OcC")]
rm_cells <- sample(rm_cells, length(rm_cells) - 200)
rgc_sh <- rgc_sh[, setdiff(colnames(rgc_sh), rm_cells)]

sh_list <- list(FC = 1:18, 
				MSC = 3:18,
				OcC = 1:18,
				TC = 12:18)
sh_cls <- rep(names(sh_list), sapply(sh_list, length)) %>%
			paste0(., "|", unlist(sh_list)) %>%
			paste0("shared", "|", .)
rgc_sh@meta.data$avgcls <- paste0(rgc_sh@meta.data$lineage, "|", rgc_sh@meta.data$lobe, "|", rgc_sh@meta.data$bin)
rgc_sh <- subset(rgc_sh, avgcls %in% sh_cls)



##------------------------------------------------------------
## tRG lineage
rgc_trg <- subset(rgc, lineage == "tRG")
rm_cells1 <- rownames(rgc_trg@meta.data)[rgc_trg@meta.data$cluster == "tRG CRYAB MEST" & rgc_trg@meta.data$cbnage %in% c("E93", "E110")]
rm_cells2 <- rownames(rgc_trg@meta.data)[rgc_trg@meta.data$cluster %in% c("tRG CRYAB FNDC1", "Ependymal")]
rgc_trg <- rgc_trg[, setdiff(colnames(rgc_trg), c(rm_cells1, rm_cells2))]
rgc_trg@meta.data$bin <- as.character(rgc_trg@meta.data$bin)
rgc_trg@meta.data$bin[rgc_trg@meta.data$bin %in% as.character(c(23,24))] <- "23"
rgc_trg@meta.data$bin[rgc_trg@meta.data$bin %in% as.character(c(25,26))] <- "24"
rgc_trg@meta.data$bin <- as.integer(rgc_trg@meta.data$bin)
trg_list <- list(FC = 19:24, 
				MSC = 19:24,
				OcC = 19:24,
				TC = 19:24)
trg_cls <- rep(names(trg_list), sapply(trg_list, length)) %>%
			paste0(., "|", unlist(trg_list)) %>%
			paste0("tRG", "|", .)
rgc_trg@meta.data$avgcls <- paste0(rgc_trg@meta.data$lineage, "|", rgc_trg@meta.data$lobe, "|", rgc_trg@meta.data$bin)
rgc_trg <- subset(rgc_trg, avgcls %in% trg_cls)



##------------------------------------------------------------
## oRG linage
rgc_org <- subset(rgc, lineage == "oRG")
rgc_org@meta.data$bin <- as.character(rgc_org@meta.data$bin)
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("28", "29")] <- "28"
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("30", "31", "32")] <- "30"
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("38", "39")] <- "37"
rgc_org@meta.data$bin <- as.integer(rgc_org@meta.data$bin)
org_list <- list(FC = c(19:28, 30, 33:37), 
				MSC = c(19:28, 30, 33:37),
				OcC = c(19:28, 30, 33:37),
				TC = c(19:28, 30, 33:37))
org_cls <- rep(names(org_list), sapply(org_list, length)) %>%
			paste0(., "|", unlist(org_list)) %>%
			paste0("oRG", "|", .)
rgc_org@meta.data$avgcls <- paste0(rgc_org@meta.data$lineage, "|", rgc_org@meta.data$lobe, "|", rgc_org@meta.data$bin)
rgc_org <- subset(rgc_org, avgcls %in% org_cls)


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
rownames(pmeta) <- pmeta$avgcls
save(avgs, pmeta, rgc_cbn, file = paste0(inputdir, "Pseudobulk_by_region_tRG_oRG.Rdata"))



##--------------------------------------------------------------------------------
## Also prepare background average expression
all_regs <- c("FC","MSC","TC","OcC")
bg_avg <- lapply(all_regs, function(reg) {
	new_cbn <- rgc_cbn
	new_cbn$lobe <- ifelse(new_cbn$lobe == reg, reg, paste0(reg, ".bg"))
	new_cbn$avgcls <- paste0(new_cbn$lineage, "|", new_cbn$lobe, "|", new_cbn$bin)
	Idents(new_cbn) <- "avgcls"
	avg <- log(AverageExpression(new_cbn)$RNA + 1)
	avg
	}) %>%
	do.call(cbind, .)
saveRDS(bg_avg, file = paste0(inputdir, "Pseudobulk_bgexpr_tRG_oRG.Rdata"))
















