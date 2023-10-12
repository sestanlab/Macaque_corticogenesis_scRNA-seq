library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)


file_name <- paste0("AdultFetal_InN")
seu <- readRDS(file = paste0("./load_files/", file_name, ".slim.rds"))


## Update cluster names
seu$subtype2 <- ifelse(is.na(seu$subtype), "unknown", seu$subtype)
seu$mres2 <- ifelse(is.na(seu$mres), "unknown", seu$mres)
seu$mres2[seu$mres2 %in% c("TH")] <- "SST"
tb <- read.table("./load_files/inn.merge.ident.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE, comment.char = "")
cls_cols <- setNames(tb$color, tb$cluster)


## color schemes
mres_cols <- c("#0d4ea3", "#2a70e8", "#0dabd6", "#a1441d", "#a1441d", "#bdb82a", "#f01fcd", "#b01e34", "#e0e0e0") %>%
		setNames(., c("VIP", "ADARB2", "LAMP5_RELN", "PVALB", "PVALB_CHC", "LAMP5_LHX6", "SST", "SST_NPY", "unknown"))
##dataset_cols <- c("#bf812d", "#35978f") %>% setNames(., c("This study", "Ma et al., 2022"))
dataset_cols <- c("#6cf58c", "#434544") %>% setNames(., c("This study", "Ma et al., 2022"))
col_list <- list(subtype2 = cls_cols, 
				mres2 = mres_cols,
				dataset = dataset_cols) 



## Visualization
pdata <- cbind(seu@meta.data[, c("mres2", "subtype2", "age")], seu$umap@cell.embeddings) %>%
			mutate(dataset = ifelse(!is.na(age), "This study", "Ma et al., 2022"))
set.seed(42)
pdata <- pdata[sample(1:nrow(pdata)), ]


plist <- lapply(c("dataset", "subtype2", "mres2"), function(tp) {
	if (tp %in% c("subtype2", "mres2")){
		subdata <- pdata[pdata[, tp] != "unknown", ]
	} else {
		subdata <- pdata
	}
	tp_size <- ifelse(tp %in% c("subtype2", "mres2"), 0.2, 0.1)
	p <- ggplot(subdata, aes_string(x = "UMAP_1", y = "UMAP_2", color = tp)) +
				ggrastr::rasterise(geom_point(size = tp_size, shape = 16), dpi = 300, scale = 1) +
	            theme_classic() + 
	            scale_color_manual(values = col_list[[tp]]) + 
	            xlim(floor(min(pdata$UMAP_1)),ceiling(max(pdata$UMAP_1)))  + 
	            ylim(floor(min(pdata$UMAP_2)),ceiling(max(pdata$UMAP_2))) + 
	            theme(legend.position = "none",
	                    line = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
	                    axis.text.x=element_blank(),axis.text.y=element_blank(), 
	                    plot.title = element_blank())
	return(p)
	})

pdf(paste0("./report/", "AdultFetal_InN_idents.pdf"), width = 24, height = 8)
plot <- patchwork::wrap_plots(plotlist = plist, nrow = 1, ncol = 3)
print(plot)
dev.off()




































