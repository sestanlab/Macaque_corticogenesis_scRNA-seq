## Two versions:
source("../scripts/nhpf.fun.R")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
future::plan("multiprocess", workers = 4)



##------------------------------------------------------
## Functions to build the cellchat object
RunCellChat <- function(data, meta, group.by) {
	cellchat <- createCellChat(object = data, meta = meta, group.by = group.by)


	## Load database
	dir_use <- paste0("./load_files/", "cellchat_custom/")

	interaction_input <- read.csv(file = paste0(dir_use, 'interaction_input_CellChatDB.csv'), row.names = 1)
	complex_input <- read.csv(file = paste0(dir_use, 'complex_input_CellChatDB.csv'), row.names = 1)
	cofactor_input <- read.csv(file = paste0(dir_use, 'cofactor_input_CellChatDB.csv'), row.names = 1)
	geneInfo <- read.csv(file = paste0(dir_use, 'geneInfo_input_CellChatDB.csv'), row.names = 1)
	CellChatDB <- list()
	CellChatDB$interaction <- interaction_input
	CellChatDB$complex <- complex_input
	CellChatDB$cofactor <- cofactor_input
	CellChatDB$geneInfo <- geneInfo


	cellchat@DB <- CellChatDB
	message("Working on DE analysis")
	cellchat <- subsetData(cellchat) %>%
					identifyOverExpressedGenes(object = ., thresh.pc = 0.05) %>%
					identifyOverExpressedInteractions(object = .)
	message("Finish DE analysis")


	message("Working on community detection")
	cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "truncatedMean", trim = 0.05, do.fast = FALSE) %>%
					filterCommunication(object = ., min.cells = 10) %>%
					computeCommunProbPathway(object = .)
	message("Finish community detection")
	return(cellchat)
}


dv_slim <- readRDS(file = paste0(inputdir, "PAT_RGC_interactionn_seu.rds"))

ccres <- RunCellChat(data = dv_slim$RNA@data, meta = dv_slim@meta.data, group.by = "chatcls")
saveRDS(ccres, file = paste0(inputdir, "Cellchat_res_custom.rds"))








