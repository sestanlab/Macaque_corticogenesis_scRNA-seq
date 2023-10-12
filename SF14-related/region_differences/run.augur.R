args <- commandArgs(trailingOnly = TRUE)
library(Augur)
library(Seurat)
inputdir <- "./load_files/"


ag <- args[1] ## E93, E110
seu <- readRDS(file = paste0(inputdir, "Augur_predata_seu.rds"))
hvg <- readRDS(file = paste0(inputdir, "Augur_predata_hvg.rds"))


reg_ord <- c("VFC","DFC","MFC","OFC","PCC","M1C","S1C","IPC","A1C","STC","ITC","V1C","Insula")
all_pairs <- combn(reg_ord, 2)
reg1 <- all_pairs[1, ][as.numeric(args[2])]
reg2 <- all_pairs[2, ][as.numeric(args[2])]

seu <- subset(seu, cbnage == ag & region %in% c(reg1, reg2))


expr <- as.matrix(seu$RNA@data[hvg, ])
meta <- seu@meta.data[, c("region", "group")]


augres <- calculate_auc(expr, meta, cell_type_col = "group", label_col = "region", augur_mode = "velocity", n_threads = 6)
saveRDS(augres, file = paste0(inputdir, "Augur_res/", "Augur_rawRes_", ag, ".", reg1, ".", reg2, ".rds"))


##codes <- paste0("Rscript run.augur.R ", rep(c("E93", "E110"), each = 78), " ", rep(1:78, by = 2))
##writeLines(codes, con = "run.augur.task.txt")

##dSQ generate the sh file
##module load dSQ
##dsq --job-file run.augur.task.txt --batch-file run.augur.sh -c 6 -p scavenge --mem-per-cpu=15G -J augur --max-jobs 24 -N 1


