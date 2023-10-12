library(Seurat)
library(dplyr)
library(tibble)
inputdir <- "./load_files/"
outputdir <- "./report/"

source("./ptime.fun.v2.R")


# Load data
rgc <- readRDS(file = paste0("../overview/load_files/", "RGC_filtered_seu_04112022_addsubset.rds"))
rgc <- rgc[, !rgc@meta.data$cluster %in% "RGC FABP7 PMP22"]
Idents(rgc) <- factor(rgc@meta.data$cluster, levels = c("NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "oRG HOPX TNC", "oRG HOPX APOE", "tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"))


##--------------------------------------------------------------------------------
## Manually remove certain clusters for more accurate estimation
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

sh_list <- list(FC = 1:18, 
				MSC = 3:18,
				OcC = 1:18,
				TC = 12:18)
sh_cls <- rep(names(sh_list), sapply(sh_list, length)) %>%
			paste0(., "|", unlist(sh_list)) %>%
			paste0("shared", "|", .)
rgc_sh@meta.data$avgcls <- paste0(rgc_sh@meta.data$lineage, "|", rgc_sh@meta.data$lobe, "|", rgc_sh@meta.data$bin)
rgc_sh <- subset(rgc_sh, avgcls %in% sh_cls)



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



rgc_cbn <- Reduce(merge, list(rgc_sh, rgc_trg))
save(rgc_cbn, file = paste0(inputdir, "RGC_combine_NEP-tRG_seu.Rdata"))
















