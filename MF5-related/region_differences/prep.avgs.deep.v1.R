library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tibble)
source("./heatmap.fun.R")

inputdir <- "./load_files/"


##--------------------------------------------------------------------------
## Load the full IPC-ExN data
seu <- readRDS(file = paste0("./load_files/ExN.harmony.spread_region.v1.harmony.rds"))


##-----------------------------------------------------------------------
## Get continuous expression for both early and late stage

## Add pseudotime (also remove some cells at this step)
load(file = paste0("./load_files/intermediate/", "IE_curve_harmony_nonIT_ptime.Rdata")) ##nit_ptime, nit_cycles
seu@meta.data$pseudotime <- NA
sh_deep <- intersect(rownames(nit_ptime), rownames(seu@meta.data))
seu@meta.data[sh_deep, "pseudotime"] <- nit_ptime[sh_deep, 1]



#########################################
## E62-64 to E77-78
# These NA cells are likely to be low-qual (not included in the final analysis)
sel_cls <- c("IPC EOMES NEUROG1", "IPC EOMES NHLH1 deep", "ExN deep nascent", "ExN deep KIF26A", "ExN deep NR4A2 GRID2", "ExN deep SYT6")

e67 <- seu[, (!is.na(seu@meta.data$pseudotime)) & seu@meta.data$subtype %in% sel_cls & seu@meta.data$cbnage %in% c("E62-64", "E77-78")]
min_pt <- 0
max_pt <- quantile(e67@meta.data$pseudotime, 0.999, na.rm = TRUE)
nPoints <- 25
half_inter <- (max_pt - min_pt)/(2 * (nPoints - 1))
cut_bks <- c(seq(min_pt, max_pt, length.out = nPoints) - half_inter, Inf)
e67@meta.data$bin <- cut(e67@meta.data$pseudotime, breaks = cut_bks) %>% 
                    as.numeric()
e67 <- e67[, e67$bin != 1]
e67@meta.data$avg_cls <- paste0(e67@meta.data$lobe, "|", e67@meta.data$bin)


Idents(e67) <- "avg_cls"
e67_avgs <- log(AverageExpression(e67, assay = "RNA")$RNA + 1)

e67_pmeta <- PrepPBmeta(object = e67, avg_col = "avg_cls", ptime_col = "pseudotime", bin_col = "bin", group.by = "subtype")


e67_avg_smt <- SmoothExprAcrossRegions(avgs = e67_avgs, meta = e67_pmeta, reg_ord = c("FC", "MSC", "TC", "OC"), span = 0.75) 
e67_meta_smt <- PrepSMmeta(avgs = e67_avg_smt, object = e67, ptime_col = "pseudotime", group.by = "subtype", ncells = 25)
e67_svg_smt <- e67_avg_smt %>%
            as.matrix() %>%
            t() %>% scale() %>% t() %>%
            MinMax(., min = -1.5, max = 2)
e67_svg_smt[is.na(e67_svg_smt)] <- -1.5

save(e67_avgs, e67_pmeta, e67_avg_smt, e67_svg_smt, e67_meta_smt, file = paste0(inputdir, "Module_continuous_expr_deep_E60-70.Rdata"))





