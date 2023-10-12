source("../scripts/nhpf.fun.R")


pfc_inn <- readRDS(file = "~/project/PFC/data/Final_PFC_HPRC.InN.rm1680.07102021.rds")
rhe_inn <- subset(pfc_inn, species == "Rhesus")
saveRDS(rhe_inn, file = paste0(inputdir, "Rhesus_adult_InN.rds"))



##rmb196 <- obj.list[["RMB196"]]
##saveRDS(rmb196, file = paste0(inputdir, "RMB196_ExN.rds"))



## P 




