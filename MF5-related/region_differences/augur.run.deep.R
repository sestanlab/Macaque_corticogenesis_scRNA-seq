args <- commandArgs(trailingOnly = TRUE)
source("../scripts/nhpf.fun.R")
library(Augur)



ag <- args[1] ## E93, E110
seu <- readRDS(file = paste0("./load_files/", "Augur_seu_deep.rds"))
hvg <- readRDS(file = paste0("./load_files/", "ExN_HVG_spread_08312022.rds"))


reg_ord <- c("FC","MSC","TC","OC")
all_pairs <- combn(reg_ord, 2)
reg1 <- all_pairs[1, ][as.numeric(args[2])]
reg2 <- all_pairs[2, ][as.numeric(args[2])]

seu <- subset(seu, cbnage == ag & lobe %in% c(reg1, reg2))


expr <- as.matrix(seu$RNA@data[hvg, ])
meta <- seu@meta.data[, c("lobe", "subtype")]


augres <- calculate_auc(expr, meta, cell_type_col = "subtype", label_col = "lobe", augur_mode = "velocity", n_threads = 6)
saveRDS(augres, file = paste0(inputdir, "Augur_res/", "Augur_rawRes_deep_", ag, ".", reg1, ".", reg2, ".rds"))


##codes <- paste0("Rscript augur.run.deep.R ", rep(c("E54", "E62-64"), each = 6), " ", rep(1:6, by = 2))
##writeLines(codes, con = "augur.run.deep.task.txt")

##dSQ generate the sh file
##module load dSQ
##dsq --job-file augur.run.deep.task.txt --batch-file augur.run.deep.sh -c 6 -p scavenge --mem-per-cpu=15G -J augur --max-jobs 12 -N 1





