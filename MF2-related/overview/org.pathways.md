---
title: Organize signaling pathways
author: Shaojie Ma
date: Dec 19, 2022

---


## Read the human gene Ontology data
```R
library(dplyr)
go <- read.csv("./load_files/goa_human.gaf", skip = 41, header = FALSE, stringsAsFactors = FALSE, sep = "\t")

slim_go <- go[, c(3, 4, 5, 6, 10, 12)]
names(slim_go) <- c("gene", "type", "GO_ID", "GO_REF", "fullname", "gene_type")


pathway_ids <- list(BMP = c("GO:0030509"),
                    FGF = c("GO:0008543"),
                    NOTCH = c("GO:0007219"),
                    EPH = c("GO:0048013", "GO:1901187", "GO:1901189", "GO:1901188")[1],
                    WNT = c("GO:0016055", "GO:0060070", "GO:0035567"),
                    SHH = c("GO:0007224"), 
                    RA = c("GO:0048384"))
pathways <- lapply(names(pathway_ids), function(id) {
                    gg <- slim_go$gene[slim_go$GO_ID %in% pathway_ids[[id]]] %>% unique()
                    return(gg)
                    }) %>%
                    setNames(., names(pathway_ids))

## system("cp /home/sm2726/project/cortex_development/dorsal_ana/genecluster/load_files/5.go/Pathway_gset.slim.rds ./load_files/")
orig.pathways <- readRDS("./load_files/Pathway_gset.slim.rds")


new.list.names <- union(names(orig.pathways), names(pathways))
new.list <- lapply(new.list.names, function(x) {
                union(orig.pathways[[x]], pathways[[x]])
                }) %>%
                setNames(., new.list.names)
saveRDS(new.list, file = "./data/Pathway_updated_20221219.rds")

```





