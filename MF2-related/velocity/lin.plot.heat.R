args <- commandArgs(trailingOnly = TRUE)
source("../scripts/nhpf.fun.R")
library(ComplexHeatmap)


data_name <- args[1]##c("Hem", "FGF17", "NKX21", "Antihem")[1]
data <- read.table(file = paste0(inputdir, "Heatdata.", data_name, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
rownames(data) <- data[, 1]
data <- data[, -1] %>% as.matrix() %>%
                        t()
meta <- read.table(file = paste0(inputdir, "Latentime.", data_name, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
                        column_to_rownames("X") %>%
                        .[colnames(data), ]

mat <- data %>% 
        t() %>% scale() %>% t() %>%
        MinMax(., min = -2.5, max = 2.5) #%>%
        ##.[, order(seu@meta.data$latent_time, decreasing = FALSE)]


set.seed(0)
gene2show <- switch(data_name, 
		AV = c("LHX6", "NKX6-2", "FABP7", "OTX2", "FEZF2", "FOS", "FOXG1", "FOXO1", "CALD1", "NKX2-1", "ASPM", "TOP2A", "SPRY1", "DLL3", "DLX1", "DLL1", "SOX2", "ZIC4", "PROX1", "PAX6", "ARX", "EPHA4", "LHX6", "ZIC1", "FGF19", "MEIS2", "DCX", "FGF14"),
        Hem = c("RSPO2", "RSPO1", "EPHA4", "BAMBI", "RSPO3", "WNT8B", "NR2F1", "EPHA7", "ARX", "GLI3", "WLS", "FEZP2", "MKI67", "TOP2A", "NEUROD1", "TP73", "LHX9", "EFNA5","LHX2", "TBR1", "BCL11A", "CALB2", "FGF14", "RELN", "DCX", "NRXN1"),
        Antihem = c("SFRP2", "LHX2", "CCND1", "ID4", "NR2F1", "NR2F2", "DLX1", "DLX2", "MEIS2", "ISL1", "SLIT1", "NNAT", "SP8", "ARX", "SP9", "PROX1", "ZIC1", "CHL1", "STMN2", "ZIC4", "ID4", "FOS", "COL2A1", sample(rownames(data), 3)),
        FGF17 = c("PAX6", "SP8", "FOS", "BMP7", "FGF17", "FGF8", "DLL1", "ST18", "POU2F2", "ONECUT2", "ONECUT1", "TAGLN3", "DLK1", "LHX9", "TBR1", "EBF1", "ISL1", "NTNG1", "POU6F2", "FGF14", "LHX9", sample(rownames(data), 3)))
col_vec <- read.table(file = paste0(inputdir, "cls.col.txt"), header = FALSE, row.names = 1, stringsAsFactors = FALSE, sep = "\t", comment.char = "") %>% 
			as.matrix() %>%
			.[, 1]
lin_list <- list(AV = c("AntVen NKX2-1 LMO1", "inIPC ASCL1 DLX1", "InN LHX8 ZIC1", "InN GNRH1", "InN HMX1"),
				Hem = c("PC RSPO3", "IPC RSPO3 NEUROG1", "IPC RSPO3 NHLH1", "CR TP73"),
				FGF17 = c("PC FGF17", "IPC FGF17", "Neu TAGLN3 ONECUT2"),
				Antihem = c("PC SFRP2", "IPC SFRP2 ASCL1", "InN MEIS2")
				)
cls_col <- col_vec[lin_list[[data_name]]]
## Initial clustering of the shared markers
col_fun = circlize::colorRamp2(c(-2, -1, 0, 1, 2), viridis(5))
all_cls <- levels(as.factor(meta$cluster))
ht <- Heatmap(mat, name = "Expression",  
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = FALSE,
        column_title_gp = gpar(fontsize = 12),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 4),
        column_title_rot = 90,
        show_column_names = FALSE,
        use_raster = TRUE,
        raster_quality = 4, top_annotation = HeatmapAnnotation(bar = meta$cluster, col = list(bar = cls_col)))
ht <- ht + 
                rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% gene2show), 
                                    labels = rownames(mat)[which(rownames(mat) %in% gene2show)], 
                                    side = "right",
                                    labels_gp = gpar(fontsize = 20), padding = unit(1, "mm")),
                                name = "Custom")
pdf(paste0(outputdir, "Lineage_cascade_custom_", data_name, ".pdf"), width = 10, height = 10)
draw(ht)
dev.off()








