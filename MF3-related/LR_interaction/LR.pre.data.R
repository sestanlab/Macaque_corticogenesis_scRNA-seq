###############################################################################################

### Ligand-receptor interaction pairs between Patterning centers and radial glia cells   

###############################################################################################
source("../scripts/nhpf.fun.R") 


## PAT clusters (GE will be extracted separted in below, with additional vRG included)
pat <- readRDS(file = paste0("../MF_pat_v3/load_files/", "PAT_inte.organizer.inte.rds"))
pat@meta.data$chatcls <- pat@meta.data$subclass
pat <- subset(pat, chatcls %in% c("PC FGF17", "PC NKX2-1", "PC RSPO3", "PC TTR")) ## further add GE RGCs


## GE RGCs
ven <- readRDS(file = paste0("~/project/cortex_development/recover_shh/load_files/", "ventral_seu_addSHH_orgident.rds")) %>%
            subset(cell_origin %in% c("E37_GE_0711", "E42_GE_1018", "E43_GE_0926"))
set.seed(42)
slim_cells <- lapply(c("GE_RG_NKX2-1_DLK1", "GE_RG_NKX2-1_OLIG1", "GE_RG_NKX6-2", "GE_vRG"), function(x) {
	cells <- colnames(ven)[ven@meta.data$hres == x]
	if (length(cells) > 1000){
		cells <- sample(cells, 1000)
	}
	return(cells)
	}) %>%
		unlist() %>%
		unique()
ven <- ven[, slim_cells]
ven@meta.data$chatcls <- "GE NERG-early"




## Test gene expression
p <- DotPlot(object = ven, features = c("ASCL1", "GADD45G", "PPP1R17", "VIM", "NES", "ID3", "CCND1", "HMGA2", "PAX6", "EMX1", "LGR4", "EPHA3", "FZD5", "FGFR1", "FGFR3"), cols = c("lightgrey", "darkred"), dot.scale = 5,dot.min = 0.025, group.by = "hres", scale.by = "size") + 
		coord_flip() + 
        RotatedAxis() + 
        theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))
pdf(paste0(outputdir, "Testst2.pdf"), width = 8, height = 8)
print(p)
dev.off()



rgc <- readRDS(file = paste0("../RGCanalysis/load_files/", "RGC_filtered_seu_04112022.rds"))
rgc <- subset(rgc, cluster %in% c("NEP RSPO3 DIRAS3", "NEP RSPO3 TEX15", "vRG HMGA2 CCND1") & lobe %in% c("FC", "OcC"))
rgc@meta.data$chatcls <- paste0(rgc@meta.data$lobe, " ", "NERG-early")


patrgc <- merge(x = pat, y = list(ven, rgc))


## Subset cluster size
all_cls <- levels(as.factor(patrgc@meta.data$chatcls))
set.seed(0)
slim_cells <- lapply(all_cls, function(x) {
	cells <- colnames(patrgc)[patrgc@meta.data$chatcls == x]
	if (length(cells) > 2000){
		cells <- sample(cells, 2000)
	}
	return(cells)
	}) %>%
		unlist() %>%
		unique()
dv_slim <- patrgc[, slim_cells]
##table(dv_slim@meta.data$chatcls, dv_slim$RNA@data["EPHA3",] >0) %>% apply(., 1, function(x) x/sum(x)) %>%t()
saveRDS(dv_slim, file = paste0(inputdir, "PAT_RGC_interactionn_seu.rds"))



## Write out the meta data
meta_data <- dv_slim@meta.data %>% 
			rownames_to_column("Cell") %>%
			mutate(cell_type = chatcls) %>%
			select(Cell, cell_type)
write.table(meta_data, file = paste0(inputdir, "PAT_RGC_meta.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



## Write out the count matrix
exp_genes <- rownames(patrgc)[Matrix::rowSums(patrgc$RNA@counts != 0) > 20]
count_data <- dv_slim$RNA@counts[exp_genes, ] %>% as.data.frame() %>% rownames_to_column("Gene")
write.table(count_data, file = paste0(inputdir, "PAT_RGC_count.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



## Test gene expression
p <- DotPlot(object = dv_slim, features = c("FGF18", "FGFR1", "FGFR3", "WNT5A", "FZD5", "PTPRK", "EFNA2", "EPHA3", "RSPO2", "LGR4"), cols = c("lightgrey", "darkred"), dot.scale = 5,dot.min = 0.05, group.by = "chatcls", scale.by = "size") + 
		coord_flip() + 
        RotatedAxis() + 
        theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))
pdf(paste0(outputdir, "Testst.pdf"), width = 8, height = 8)
print(p)
dev.off()












