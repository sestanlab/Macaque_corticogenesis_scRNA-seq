## In this analysis, all the cells in human and macaque were downsampled to the same depth (1500UMIs)
## Then average expression/ratios were calculated for dotplot visualization
## Prepare the data for visualization (DotPlots)
library(Seurat)
library(dplyr)


mac <- readRDS(file = paste0("./load_files/Macaque_E77-78_slim.rds"))
human <- readRDS(file = paste0("./load_files/Human_GW18-20_slim.rds"))


mac <- subset(mac, subtype %in% c("ExN CUX2 ACTN2", "ExN SOX5 NR4A2 GRID2", "ExN SOX5 SYT6"))
mac$group <- ifelse(mac$subtype %in% "ExN CUX2 ACTN2", "ExN upper", "ExN deep")
mac$species <- "Macaque"
human <- subset(human, subclass %in% c("L2-5 IT", "L6B ST18", "L6CT SYT6"))
human$group <- ifelse(human$subclass %in% "L2-5 IT", "ExN upper", "ExN deep")
human$region <- gsub("PFC", "FC", human$regioncode) %>%
				gsub("V1C", "OC", .)
human$species <- "Human"


## Merge data from the two species
sh_genes <- intersect(rownames(human), rownames(mac))
hm <- merge(x = human[sh_genes, ], y = mac[sh_genes, ])
hm <- subset(hm, region %in% c("FC", "OC"))
hm$newcls <- paste0(hm$species, "|", hm$region, "|", hm$group)
DefaultAssay(hm) <- "RNA" ## Activate the nCount_RNA column
hm <- NormalizeData(hm, normalization.method = "LogNormalize")
hm <- subset(hm, nCount_RNA >= 1500)



## Downsampling to 1500 UMIs per cell
source("./ds.fun.R")
ctm <- SubsampleData.count.equalUMI(counts = hm$RNA@counts, numi = 1500, seed.use = 42, nCores = 12)
hm[["DS"]] <- CreateAssayObject(counts = ctm[, colnames(hm)], min.cells = 0, min.features = 0)
DefaultAssay(hm) <- "DS"
hm <- NormalizeData(hm, normalization.method = "LogNormalize")
saveRDS(hm, file = paste0("./load_files/Human_Macaque_Neuron_slim_addDS_E77.rds"))



##---------------------------------------------------------------------
## Calculate avgs and ratios for visualization
hm <- readRDS(file = paste0("./load_files/Human_Macaque_Neuron_slim_addDS_E77.rds"))
hm$newcls <- paste0(hm$species, "|", hm$region, "|", hm$group)
Idents(hm) <- "newcls"
avgs <- log(AverageExpression(hm)$DS + 1)


allcls <- levels(as.factor(hm$newcls))
ratios <- lapply(allcls, function(cls) {
    subseu <- subset(hm, newcls %in% cls)
    ratio <- Matrix::rowMeans(subseu$DS@data != 0)
    return(ratio)
    }) %>%
    setNames(., allcls) %>%
    as.data.frame(., check.names = FALSE) %>%
    as.matrix()
save(avgs, ratios, file = paste0("./load_files/Human_Macaque_Neuron_AVG_DS_E77.Rdata"))

