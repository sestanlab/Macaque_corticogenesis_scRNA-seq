library(topGO)
library(org.Hs.eg.db)
printgene <- function (object, whichTerms, simplify = TRUE, geneCutOff = 50) {
        term.genes <- genesInTerm(object, whichTerms)

        retList <- vector("list", length(whichTerms))
        names(retList) <- whichTerms
        for (gt in whichTerms) {
            pval <- sort(geneScore(object, term.genes[[gt]], use.names = TRUE))
            affID <- names(pval)

            length(affID) <- min(length(affID), geneCutOff)
            if (length(affID) == 0) {
                message("\n\t No genes over the cutOff for: ", gt)
                next
            }

            retList[[gt]] <- data.frame(row.names = names(pval), gene_symbol = names(pval), raw_p = format.pval(pval, digits = 3, eps = 1e-30) %>% as.numeric(), stringsAsFactors = FALSE)
        }
        if (simplify && length(whichTerms) == 1) 
            return(retList[[1]])
        return(retList)
}



Plot_Enrich_Bar00 <- function(input_meta, id_name, text_x = 0.1, y_col = "mlogp"){
    input_meta$text_loc <- text_x
    label_name <- ifelse(id_name == "GO.ID", "Term", id_name)
    p <- ggplot(data = input_meta) + 
            geom_bar(aes_string(x = id_name, y = y_col), color = "orange", fill = "orange", stat = "identity", width = 0.6) +
            xlim(limits = input_meta[, id_name] %>% rev()) + 
            geom_text(aes_string(x = id_name, y = "text_loc", label = label_name), hjust = 0, size = 2.5) + 
            coord_flip() + 
            theme_bw() + 
            RotatedAxis()
    p
}


Plot_Enrich_Bar <- function(input_meta, id_name, text_x = 0.1){
    input_meta$text_loc <- text_x
    label_name <- ifelse(id_name == "GO.ID", "Term", id_name)
    input_meta$type <- sapply(input_meta[, id_name], function(x) {
        if (x < 2) {
            return("t3")
        } else if (x >= 2 & x < 4) {
            return("t2")
        } else {
            return("t1")
        }
        })
    color_use <- c("#fff7ec", "#fdd49e", "#fdbb84") %>% setNames(., c("t3", "t2", "t1"))
    ##color_use <- c("red", "green", "blue") %>% setNames(., c("t1", "t2", "t3"))
    color_use <- color_use[levels(as.factor(input_meta$type))]
    p <- ggplot(data = input_meta) + 
            geom_bar(aes_string(x = id_name, y = "mlogp", color = "type", fill = "type"), stat = "identity", width = 0.6) +
            xlim(limits = input_meta[, id_name] %>% rev()) + 
            scale_colour_manual(values = color_use) + 
            scale_fill_manual(values = color_use) + 
            geom_text(aes_string(x = id_name, y = "text_loc", label = label_name), hjust = 0, size = 2.5) + 
            coord_flip() + 
            theme_bw() + 
            RotatedAxis()
    p
}




DoGO <- function(all_genes, interesting_genes, pvalthre = 0.01, file_name, output_dir = outputdir){
    gene_list <- rep(0, length(all_genes)) %>% setNames(., all_genes)
    gene_list[names(gene_list) %in% interesting_genes] <- 1
    gene_list <- as.factor(gene_list)


    GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_list, annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol", nodeSize = 10)
    GO_classic <- runTest(object = GOdata, algorithm = "classic", statistic = "fisher")
    allRes <- GenTable(GOdata, classic = GO_classic, orderBy = "classic", ranksOf = "classic", topNodes = 50)
    rawres <- GenTable(GOdata, classic = GO_classic, orderBy = "classic", ranksOf = "classic", topNodes = length(GO_classic@score))
    rawres$classic[rawres$classic == "< 1e-30"] <- "1e-30"
    rawres$classic <- as.numeric(rawres$classic)


    allRes$classic[allRes$classic == "< 1e-30"] <- "1e-30"
	allRes$classic <- as.numeric(allRes$classic)
    allRes$mlogp <- -log10(as.numeric(allRes$classic))

    sign_terms <- allRes$GO.ID[as.numeric(allRes$classic) <= pvalthre]
    if (length(sign_terms) >= 1){
        gt <- printgene(GOdata, whichTerms = allRes$GO.ID, geneCutOff = 10000)
    }

    p <- Plot_Enrich_Bar(input_meta = allRes, id_name = "GO.ID") + 
            theme(axis.text.y = element_blank()) + 
            labs(y = "-log(p)", x = "GO Terms")

    pdf(paste0(output_dir, file_name, "_GO_enrichment.pdf"), width = 5, height = 3)
    ##par(c(2,2,2,2))
    print(p)
    dev.off()

    rawres$fdr <- p.adjust(rawres$classic, "fdr")
    rawres$bonferroni <- p.adjust(rawres$classic, "bonferroni")    
    return(list(res = allRes, gt = gt, raw = rawres))
}


