library(Seurat)
library(dplyr)

exn <- readRDS(file = "../overview/load_files/ExN_data_08312022.rds")


##---------------------------------------------------------------
## Prepare upper layer data (E77-78)

## remove cycling cells in the calculation
load(file = paste0("./load_files/intermediate/", "IE_curve_harmony_IT_ptime.Rdata")) ##it_ptime, it_cycle


up_cls <- c("IPC EOMES NEUROG1", "IPC EOMES NHLH1 up", "ExN CUX2 PALMD", "ExN CUX2 ADRA2A", "ExN CUX2 ACTN2")
up_meta <- exn@meta.data[exn@meta.data$cbnage %in% c("E77-78", "E93", "E110") & exn@meta.data$subtype %in% up_cls & exn@meta.data$lobe != "Insula", ]
up_meta <- up_meta[setdiff(rownames(up_meta), it_cycle), ]
up_size <- table(up_meta$subtype, up_meta$lobe, up_meta$cbnage)

up_cells <- lapply(c("E77-78", "E93", "E110"), function(ag) {
	sub_meta <- up_meta[up_meta$cbnage %in% ag, ]
	sub_size <- table(sub_meta$subtype, sub_meta$lobe)
	cls_kp <- rownames(sub_size)[apply(sub_size, 1, min) >= 100]
	sub_cells <- rownames(sub_meta)[sub_meta$subtype %in% cls_kp]
	return(sub_cells)
	}) %>%
	unlist()

up_fmeta <- up_meta[up_cells, ]
up_seu <- exn[, rownames(up_fmeta)]

saveRDS(up_seu, file = paste0("./load_files/", "Augur_seu_upper.rds"))




##---------------------------------------------------------------
## Prepare deep layer data (E54-78)

## remove cycling cells in the calculation
load(file = paste0("./load_files/intermediate/", "IE_curve_harmony_nonIT_ptime.Rdata")) ##nit_ptime, nit_cycle


exn@meta.data$subtype <- gsub("ExN SOX5 ID2", "ExN SOX5 PALMD", exn@meta.data$subtype)

## "ExN SOX5 ID2", 
deep_cls <- c("IPC EOMES NEUROG1", "IPC EOMES NHLH1 deep", "ExN SOX5 PALMD", "ExN SOX5 NR4A2 GRID2", "ExN SOX5 SYT6")
deep_meta <- exn@meta.data[exn@meta.data$cbnage %in% c("E54", "E62-64") & exn@meta.data$subtype %in% deep_cls & exn@meta.data$lobe != "Insula", ]
deep_meta <- deep_meta[setdiff(rownames(deep_meta), nit_cycle), ]
deep_size <- table(deep_meta$subtype, deep_meta$lobe, deep_meta$cbnage)


deep_cells <- lapply(c("E54", "E62-64"), function(ag) {
	sub_meta <- deep_meta[deep_meta$cbnage %in% ag, ]
	sub_size <- table(sub_meta$subtype, sub_meta$lobe)
	cls_kp <- rownames(sub_size)[apply(sub_size, 1, min) >= 100]
	sub_cells <- rownames(sub_meta)[sub_meta$subtype %in% cls_kp]
	return(sub_cells)
	}) %>%
	unlist()
deep_fmeta <- deep_meta[deep_cells, ]
table(deep_fmeta$subtype, deep_fmeta$lobe, deep_fmeta$cbnage)
deep_seu <- exn[, rownames(deep_fmeta)]
saveRDS(deep_seu, file = paste0("./load_files/", "Augur_seu_deep.rds"))
















