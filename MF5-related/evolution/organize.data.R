library(Seurat)
library(dplyr)
library(SeuratObject)


## copy & paste the previously processed human GW18-19 data from Kriegstein lab
system("cp /home/sm2726/project2/SpatialBrain/singlecell/data/Batch_QC_harmony_GW18-20.harmony.addanno.slim.rds ./load_files/")
system("cp /home/sm2726/project2/SpatialBrain/msc/ds.fun.R ./")


## Load the data
seu <- readRDS(file = paste0("./load_files/Batch_QC_harmony_GW18-20.harmony.addanno.slim.rds"))


## Update region information
seu$regioncode <- gsub("Motor", "MSC", seu$region) %>%
				gsub("Somatosensory", "MSC", .) %>%
				gsub("Parietal", "PC", .) %>%
				gsub("Prefrontal", "PFC", .) %>%
				gsub("Temporal", "TC", .) %>%
				gsub("Visual", "V1C", .)
subseu <- subset(seu, regioncode %in% c("MSC", "PFC", "TC", "V1C"))
saveRDS(subseu, file = paste0("./load_files/Human_GW18-20_slim.rds"))


## Macaque data
load("../inte_all/load_files/Macaque.developing.seurat.Rdata") ##mac, subtype_order, subclass_order
e77 <- subset(mac, cbnage %in% c("E77-78") & region %in% c("FC", "MSC", "TC", "OC"))

all_cls <- table(e77$subtype) %>% .[. >= 20] %>% names()
e77 <- subset(e77, subtype %in% all_cls)
saveRDS(e77, file = paste0("./load_files/Macaque_E77-78_slim.rds"))

