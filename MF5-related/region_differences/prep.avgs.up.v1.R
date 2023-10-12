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
load(file = paste0("./load_files/intermediate/", "IE_curve_harmony_IT_ptime.Rdata")) ##it_ptime, it_cycles
seu@meta.data$pseudotime <- NA
sh_up <- intersect(rownames(it_ptime), rownames(seu@meta.data))
seu@meta.data[sh_up, "pseudotime"] <- it_ptime[sh_up, 1]



#########################################
## E62-64
# These NA cells are likely to be low-qual (not included in the final analysis)
sel_cls <- c("IPC EOMES VIM", "IPC EOMES NEUROG1", "IPC EOMES NHLH1 up", "ExN up nascent", "ExN up ADRA2A", "ExN up ACTN2")
e67 <- seu[, (!is.na(seu@meta.data$pseudotime)) & seu@meta.data$subtype %in% sel_cls & seu@meta.data$cbnage %in% c("E62-64", "E77-78")]
min_pt <- 0
max_pt <- quantile(e67@meta.data$pseudotime, 0.999, na.rm = TRUE)
nPoints <- 25
half_inter <- (max_pt - min_pt)/(2 * (nPoints - 1))
cut_bks <- c(seq(min_pt, max_pt, length.out = nPoints) - half_inter, Inf)
e67@meta.data$bin <- cut(e67@meta.data$pseudotime, breaks = cut_bks) %>% 
                    as.numeric()
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

save(e67_avgs, e67_pmeta, e67_avg_smt, e67_svg_smt, e67_meta_smt, file = paste0(inputdir, "Module_continuous_expr_up_E60-70.Rdata"))



#########################################
## E93-110
sel_cls <- c("IPC EOMES VIM", "IPC EOMES NEUROG1", "IPC EOMES NHLH1 up", "ExN up nascent", "ExN up ADRA2A", "ExN up ACTN2")
e910 <- seu[, (!is.na(seu@meta.data$pseudotime)) & seu@meta.data$subtype %in% sel_cls & seu@meta.data$cbnage %in% c("E93", "E110")]
min_pt <- 0
max_pt <- quantile(e910@meta.data$pseudotime, 0.999, na.rm = TRUE)
nPoints <- 20
half_inter <- (max_pt - min_pt)/(2 * (nPoints - 1))
cut_bks <- c(seq(min_pt, max_pt, length.out = nPoints) - half_inter, Inf)
e910@meta.data$bin <- cut(e910@meta.data$pseudotime, breaks = cut_bks) %>% 
                    as.numeric()
e910@meta.data$avg_cls <- paste0(e910@meta.data$lobe, "|", e910@meta.data$bin)


Idents(e910) <- "avg_cls"
e910_avgs <- log(AverageExpression(e910, assay = "RNA")$RNA + 1)

e910_pmeta <- PrepPBmeta(object = e910, avg_col = "avg_cls", ptime_col = "pseudotime", bin_col = "bin", group.by = "subtype")

## Remove Insula
e910_pmeta <- e910_pmeta[e910_pmeta$lobe != "Insula", ]
e910_avgs <- e910_avgs[, rownames(e910_pmeta)]


e910_avg_smt <- SmoothExprAcrossRegions(avgs = e910_avgs, meta = e910_pmeta, reg_ord = c("FC", "MSC", "TC", "OC"), span = 0.75) 
e910_meta_smt <- PrepSMmeta(avgs = e910_avg_smt, object = e910, ptime_col = "pseudotime", group.by = "subtype", ncells = 25)
e910_svg_smt <- e910_avg_smt %>%
            as.matrix() %>%
            t() %>% scale() %>% t() %>%
            MinMax(., min = -1.5, max = 2)
e910_svg_smt[is.na(e910_svg_smt)] <- -1.5
save(e910_avgs, e910_pmeta, e910_avg_smt, e910_svg_smt, e910_meta_smt, file = paste0(inputdir, "Module_continuous_expr_up_E93-110.Rdata"))



