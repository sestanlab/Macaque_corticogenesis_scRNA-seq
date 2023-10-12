library(Seurat)
library(Matrix)
library(dplyr)

inputdir <- "./load_files/"
outputdir <- "./report/"



##-----------------------------------------------------------------------
## Extract DEGs
allres <- readRDS(file = paste0(inputdir, "Region_DEGs_rawres.rds"))

## Consider two situations
## 1. (expr ratio >= 0.1)
slim_dex1 <- allres %>%
				filter(pct.1 >= 0.1 & pct.2 <= 0.75 & ratio_fc >= 1.4 & avg_logFC >= 0.2 & p_val_adj <= 0.001)
## 2. (expr ratio < 0.1 &  0.05)
slim_dex2 <- allres %>%
				filter(pct.1 >= 0.05 & pct.1 < 0.1 & ratio_fc >= 4 & avg_logFC >= 0.1 & p_val_adj <= 0.001)
## 3. 
slim_dex3 <- allres %>%
				filter(pct.2 > 0.75 & ratio_fc >= 1.1 & avg_logFC >= 1 & p_val_adj <= 0.001)


slim_dex <- rbind(slim_dex1, slim_dex2) %>%
				rbind(., slim_dex3)



## Organize DEGs
vrg_res <- slim_dex %>%
			filter(cluster %in% c("NEP RSPO3", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "tRG CRYAB MEST"))
vrg_deg <- split(vrg_res$gene, vrg_res$region) %>%
			lapply(., unique)


org_res <- slim_dex %>%
			filter(cluster %in% c("oRG HOPX APOE", "oRG HOPX TNC"))
org_deg <- split(org_res$gene, org_res$region) %>%
			lapply(., unique)

cbn_res <- rbind(vrg_res, org_res)

save(vrg_deg, org_deg, vrg_res, org_res, cbn_res, file = paste0(inputdir, "Region_DEGs_res_v2.rds"))






