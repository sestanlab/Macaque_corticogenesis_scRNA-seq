#This scirpt contains basic functions of single cell RNA-seq analysis (as well as many other technologies) 
library(parallel)
map_gene <- function(gene_names, input_genes,ignore_case=TRUE){
    input_genes <- unique(input_genes)
  
    if (sum(grepl("\\|",gene_names))==length(gene_names)){
        if (sum(grepl("\\|",input_genes))==length(input_genes)){
              gene_id <- input_genes[input_genes %in% gene_names]
          }else{
            input_genes <- extract_field(input_genes=input_genes, 2, "|")
              gene_id <- unlist(mclapply(input_genes, function(i) ifelse(length(grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
              gene_id <- gene_id[gene_id != "empty"]
          }
    } else if(sum(grepl("\\|",gene_names))==0){
        input_genes <- extract_field(input_genes=input_genes, 2, "|")
          gene_id <- unlist(mclapply(input_genes, function(i) ifelse(length(grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
          gene_id <- gene_id[gene_id != "empty"]
    } else {
        stop("Inconsistent gene name format")
    }
    return(gene_id)
}



#Get the mito genes based on the inquiry genes
get_genes <- function(input_genes, gene_type = c("mito","ribo", "cc")[1], return_list = FALSE, revised = FALSE){
    gene_use <- list()
    if ("mito" %in% gene_type){
        mito.known <- map_gene(gene_names=input_genes, input_genes=c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")) #Refseq annotation 103 (Macaque)
        mito.combine <- grep(pattern = "\\|MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        mito.single <- grep(pattern = "^MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        gene_use[["mito"]] <- unique(c(mito.known, mito.combine, mito.single))
    }

    if ("ribo" %in% gene_type){
        ribo.combine <- c(grep("\\|RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        ribo.single <- c(grep("^RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        gene_use[["ribo"]] <- c(ribo.combine, ribo.single)
    }

    if ("cc" %in% gene_type){
        if (!revised){
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        } else {
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        }
        
        gene_use[["s"]] <- map_gene(gene_names=input_genes, input_genes=regev.s.genes,ignore_case=TRUE)
        gene_use[["g2m"]] <- map_gene(gene_names=input_genes, input_genes=regev.g2m.genes,ignore_case=TRUE)
    }

    if (return_list){
        return(gene_use)
    } else {
        all_genes <- setNames(unlist(gene_use), NULL)
        return(all_genes)
    }
}



extract_field <- function(input_genes, field = 1, split_by = "_") {
    split_strings <- strsplit(input_genes, split_by, fixed=TRUE)
    if (is.numeric(field)){
        if (all(field > 0)){
            if (length(field) == 1){
                idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
                output_strings <- sapply(1:length(split_strings), function(x) split_strings[[x]][idx[x]])
            } else {
                idx <- sapply(split_strings, length) == 1
                output_strings <- sapply(1:length(split_strings), function(x) ifelse(idx[x], input_genes[x], paste(split_strings[[x]][field], collapse = split_by)))
            }
        } else {
            if (length(field) == 1) {
                field = as.integer(abs(field))
                idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
                output_strings <- sapply(1:length(split_strings), function(x) rev(split_strings[[x]])[idx[x]])
            } else {
                stop("currently doesnot support field with length >= 2")
            }
        }
    } else if (is.character(field)) {
        if (field == "rm_start"){
            idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
            output_strings <- sapply(1:length(split_strings), function(x) paste(split_strings[[x]][-1],collapse=split_by))
        } else if (field == "rm_end") {
            output_strings <- sapply(1:length(split_strings), function(x) paste(rev(rev(split_strings[[x]])[-1]),collapse=split_by))
        } else {
            stop("Currently only support rm_start, rm_end for the character field")
        }
    }
    output_strings <- output_strings %>% setNames(., input_genes)
    return(output_strings)
}
