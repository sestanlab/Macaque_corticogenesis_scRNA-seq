source("../scripts/nhpf.fun.R")



## Get the means and average for each pair
outDir <- paste0(inputdir, "CPDBraw/PATRGCv1/")
means <- read.table(paste0(outDir, "means.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1, check.names = FALSE)
pvals <- read.table(paste0(outDir, "pvalues.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1, check.names = FALSE) 


## Set interaction pairs
pc_clusters <- c("PC FGF17", "PC NKX2-1", "PC RSPO3", "PC TTR")
rgc_clusters <- c("FC NERG-early", "GE NERG-early", "OcC NERG-early")
cls_pairs <- expand.grid(rgc_clusters, pc_clusters) %>%
                mutate(Var1 = as.character(Var1), Var2 = as.character(Var2)) %>%
                subset(Var1 != Var2) %>%
        mutate(pair = paste0(Var2, "|", Var1)) %>%
        .$pair

means <- means[, cls_pairs]
pvals <- pvals[rownames(means), cls_pairs]


meta <- read.table(paste0(outDir, "/means.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t", row.names = 1)[, 1:10] %>%
	rownames_to_column("pair_ID") %>%
    filter(receptor_a != "True"& receptor_b == "True") %>%
    column_to_rownames("pair_ID")


pos_pairs <- rownames(meta)[apply(pvals[rownames(meta), ], 1, min) <= 0.05]
pval_cpb <- pvals[pos_pairs, ]
mean_cpb <- means[pos_pairs, ]
meta_cpb <- meta[pos_pairs, ]

rownames(pval_cpb) <- meta_cpb[rownames(pval_cpb), "interacting_pair"]
rownames(mean_cpb) <- meta_cpb[rownames(mean_cpb), "interacting_pair"]
rownames(meta_cpb) <- meta_cpb$interacting_pair

save(pval_cpb, mean_cpb, meta_cpb, file = paste0(inputdir, "Cellphone_filtered_res.rds"))



