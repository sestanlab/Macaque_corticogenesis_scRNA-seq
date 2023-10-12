library(Seurat)
library(dplyr)
library(tibble)
inputdir <- "./load_files/"
outputdir <- "./report/"

source("./ptime.fun.v2.R")


# Load data
rgc <- readRDS(file = paste0(inputdir, "RGC_filtered_seu_04112022_addsubset.rds"))
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



rgc_org <- subset(rgc, lineage == "oRG")
rgc_org@meta.data$bin <- as.character(rgc_org@meta.data$bin)
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("28", "29")] <- "28"
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("30", "31", "32")] <- "30"
rgc_org@meta.data$bin[rgc_org@meta.data$bin %in% c("38", "39")] <- "37"
rgc_org@meta.data$bin <- as.integer(rgc_org@meta.data$bin)
rgc_org@meta.data$avgcls <- paste0(rgc_org@meta.data$lineage, "|", rgc_org@meta.data$lobe, "|", rgc_org@meta.data$bin)
org_list <- list(FC = c(19:28, 30, 33:37), 
				MSC = c(19:28, 30, 33:37),
				OcC = c(19:28, 30, 33:37),
				TC = c(19:28, 30, 33:37))
org_cls <- rep(names(org_list), sapply(org_list, length)) %>%
			paste0(., "|", unlist(org_list)) %>%
			paste0("oRG", "|", .)
rgc_cbn <- subset(rgc_org, avgcls %in% org_cls)


save(rgc_cbn, file = paste0(inputdir, "RGC_combine_oRG_seu.Rdata"))








