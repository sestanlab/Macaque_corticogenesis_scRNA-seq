source("~/project/PFC/scripts/pfc.fun.R")
source("./ptime.fun.R")
##source("./rgc.fun.R")

rgc <- readRDS(file = paste0("./load_files/RGC_filtered_seu_04112022.rds"))



file_name <- "RGC_shmars"
sel_cls <- c("NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1", "vRG SAT1 STMN2", "oRG HOPX TNC", "oRG HOPX APOE", "tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal")
rgc <- rgc[, rgc@meta.data$cluster %in% sel_cls]
rgc@meta.data$cluster <- gsub(" DIRAS3| TEX15", "", rgc@meta.data$cluster)
reg_list <- SplitObject(rgc, split.by = "lobe")
for (ii in names(reg_list)){
	Idents(reg_list[[ii]]) <- "cluster"
}



##--------------------------------------------------------------------------------
## First identity lineage-specific genes (shared, oRG, tRG)
abbr <- list(shared = c("NEP RSPO3", "vRG HMGA2 CCND1", "vRG SAT1 STMN2"), 
			tRG = c("tRG CRYAB MEST", "tRG CRYAB FNDC1", "Ependymal"),
			oRG = c("oRG HOPX TNC", "oRG HOPX APOE"))
comp1 <- list(s1 = "NEP RSPO3", 
			s2 = "vRG HMGA2 CCND1", 
			s3 = "vRG SAT1 STMN2", 
			s4 = "NEP RSPO3", 
			s5 = "vRG HMGA2 CCND1", 
			s6 = "vRG SAT1 STMN2", 
			o1 = "oRG HOPX TNC", 
			o2 = "oRG HOPX APOE",
			o3 = "oRG HOPX TNC", 
			o4 = "oRG HOPX APOE", 
			t1 = "tRG CRYAB MEST", 
			t2 = "tRG CRYAB FNDC1", 
			t3 = "Ependymal", 
			t4 = "tRG CRYAB MEST", 
			t5 = "tRG CRYAB FNDC1", 
			t6 = "Ependymal")
comp2 <- list(s1 = abbr[["tRG"]],
			s2 = abbr[["tRG"]], 
			s3 = abbr[["tRG"]], 
			s4 = abbr[["oRG"]],
			s5 = abbr[["oRG"]], 
			s6 = abbr[["oRG"]], 
			o1 = abbr[["shared"]], 
			o2 = abbr[["shared"]], 
			o3 = abbr[["tRG"]], 
			o4 = abbr[["tRG"]], 
			t1 = abbr[["shared"]], 
			t2 = abbr[["shared"]], 
			t3 = abbr[["shared"]], 
			t4 = abbr[["oRG"]], 
			t5 = abbr[["oRG"]], 
			t6 = abbr[["oRG"]])
comp_regs <- list(s1 = c("FC", "OcC"),
			s2 = c("FC", "OcC"), 
			s3 = c("FC", "MSC", "OcC", "TC"), 
			s4 = c("FC", "OcC"),
			s5 = c("FC", "OcC"), 
			s6 = c("FC", "MSC", "OcC", "TC"), 
			o1 = c("FC", "MSC", "OcC", "TC"), 
			o2 = c("FC", "MSC", "OcC", "TC"), 
			o3 = c("FC", "MSC", "OcC", "TC"), 
			o4 = c("FC", "MSC", "OcC", "TC"), 
			t1 = c("FC", "MSC", "OcC", "TC"), 
			t2 = c("FC", "OcC"), 
			t3 = c("FC", "OcC"), 
			t4 = c("FC", "MSC", "OcC", "TC"), 
			t5 = c("FC", "OcC"), 
			t6 = c("FC", "OcC"))



res_list <- lapply(names(comp1), function(tt) {
	res <- FindCommonMars(obj.list = reg_list, ident.1 = comp1[[tt]], ident.2 = comp2[[tt]], regions = comp_regs[[tt]])
	res$type <- tt
	return(res)
	}) %>%
		setNames(., names(comp1))
##saveRDS(res_list, file = paste0(inputdir, "Ptime_mars_cbn_step1.rds"))


IntersectGenes <- function(res_list, id1, id2, max.genes = 500) {
	gene1 <- res_list[[id1]]$gene
	gene1 <- gene1[1:min(max.genes, length(gene1))]

	gene2<- res_list[[id2]]$gene
	gene2 <- gene2[1:min(max.genes, length(gene2))]

	return(intersect(gene1, gene2))
}


mars_early <- list(`NEP RSPO3` = IntersectGenes(res_list, id1 = "s1", id2 = "s4"),
					`vRG HMGA2 CCND1` = IntersectGenes(res_list, id1 = "s2", id2 = "s5"), 
					`vRG SAT1 STMN2` = IntersectGenes(res_list, id1 = "s3", id2 = "s6"))

mars_org <- list(`oRG HOPX TNC` = IntersectGenes(res_list, id1 = "o1", id2 = "o3"),
					`oRG HOPX APOE` = IntersectGenes(res_list, id1 = "o2", id2 = "o4"))

mars_trg <- list(`tRG CRYAB MEST` = IntersectGenes(res_list, id1 = "t1", id2 = "t4"),
					`tRG CRYAB FNDC1` = IntersectGenes(res_list, id1 = "t2", id2 = "t5"),
					`Ependymal` = IntersectGenes(res_list, id1 = "t3", id2 = "t6"))
save(res_list, mars_early, mars_org, mars_trg, file = paste0("./load_files/intermediate/", "Ptime_mars_cbn_step1.Rdata"))



