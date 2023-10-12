library(filesstrings)
library(PWMEnrich)
library(glmnet)
library(BiocParallel)
library(foreach)


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



###borrowed from dev_cluster_wrapp.R
convert_genes <- function(input_genes) {
    if (FALSE) {
        loc_geneids <- readRDS(paste0(new_inputdir, "E42_FR_1018_raw.rds")) %>% rownames() %>% grep("^LOC", ., value = TRUE)

        convert_table <- read.table("/home/sm2726/project/cortex_development/preprocessing/Ensembl2Refseq_Mmul10.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% .[!is.na(.[, "NCBI.gene.ID"]), ] %>% 
                    mutate(NCBI.gene.ID = as.character(NCBI.gene.ID)) %>% 
                    mutate(NCBI.gene.ID = paste0("LOC", NCBI.gene.ID)) %>%
                    filter(NCBI.gene.ID %in% loc_geneids) %>%
                    filter(Gene.name != "")

        missing_loc <- setdiff(loc_geneids, convert_table$NCBI.gene.ID)
        

        ##mannually check if one NCBI.ID corresponds to multiple gene names in ENSEMBLE
        ids_match <- aggregate(Gene.name ~ NCBI.gene.ID, convert_table, function(x) unique(c(x)), simplify = TRUE) %>%
                    mutate(Gene.name = as.character(Gene.name))
        ids_man <- data.frame(NCBI.gene.ID = c("LOC114675439", "LOC114675553", "LOC711652"), Gene.name = c("DEFB4", "WDR45L", "ZNF222"), stringsAsFactors= FALSE)
        ids_retain <- rbind(ids_match[grep(", ", ids_match$Gene.name, invert = TRUE), ], ids_man)
    

        ##mannually check if multiple NCBI.ID corresponds to the same ENSEMBLE ID
        ##For these genes, attach a "-c<index>" after the gene names
        duplicated_genes <- ids_retain[ids_retain$Gene.name %in% ids_retain$Gene.name[duplicated(ids_retain$Gene.name)], ] %>%
                    arrange(Gene.name) %>% mutate(new_name = "")
        for (i in unique(duplicated_genes$Gene.name)) {
            duplicated_genes[duplicated_genes$Gene.name %in% i, "new_name"] <- paste0("c", as.character(1:sum(duplicated_genes$Gene.name %in% i)))
        }
        duplicated_genes <- duplicated_genes %>% mutate(Gene.name = paste0(Gene.name, "-", new_name)) %>% .[, c("NCBI.gene.ID", "Gene.name")]


        ##Get the final list of genes
        final_genes <- rbind(ids_retain[remove_duplicates(ids_retain$Gene.name, return_index = TRUE), ], duplicated_genes)
        convert_vector <- setNames(final_genes$Gene.name, final_genes$NCBI.gene.ID)
        saveRDS(convert_vector, file = "/home/sm2726/project/cortex_development/preprocessing/convert_vector.rds")
    }

    convert_vector <- readRDS("/home/sm2726/project/cortex_development/preprocessing/convert_vector.rds")

    ##Get the ENSEMBLE gene names with multiple hits in NCBI
    loc_idx <- (input_genes %>% grepl("^LOC", .)) & (input_genes %in% names(convert_vector))
    input_genes[loc_idx] <- convert_vector[input_genes[loc_idx]]

    return(input_genes)
}


#Get the coordinate of the select the regions
get_promoter <- function(input_genes, upstream = 800, downstream = 100, out_file = "aa.bed"){
	## Read the GTF file
	gtf_slim <- readRDS("/gpfs/ycga/project/sestan/sm2726/reference/macaque/refseq_assembly/ncbi.genomes.2019.06.10/Mmul10_slim_gtf.rds")
	gtf_slim$gene_id <- gsub("_", "-", gtf_slim$gene_id) ## Consistent with Seurat analysis
	gtf_slim$gene_id <- convert_genes(input_genes = gtf_slim$gene_id) ## may have duplicated genes
	##rownames(gtf_slim) <- gtf_slim$gene_id

	## Get the mapped genes
	mapped_genes <- intersect(input_genes, gtf_slim$gene_id) %>% unique()
	print(paste0(length(mapped_genes), " genes were mapped when extracting promoters"))

	short_gtf <- gtf_slim[match(mapped_genes, gtf_slim$gene_id), ,drop = FALSE]
	out_df <- data.frame(row.names = mapped_genes, 
			chr = short_gtf$seqnames,
			start = ifelse(short_gtf$strand == "+", short_gtf$start - upstream, short_gtf$end - downstream), 
			end = ifelse(short_gtf$strand == "+", short_gtf$start + downstream, short_gtf$end + upstream), 
			name = extract_field(mapped_genes, 2, "|"),
			score = rep(1, dim(short_gtf)[1]),
			strand = short_gtf$strand,
			stringsAsFactors = FALSE)

	##Gtf to bed conversion
	out_df$start <- out_df$start - 1
	out_df$start[out_df$start < 0]  <- 0

	write.table(out_df, file = out_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
}


MotifEnrich <- function(input_genes, input_dir = new_inputdir, upstream = 800, downstream = 100, nCores = 12) {
	##-----------------------------------------------------------------------------------------------
	## Get the promoter sequences
	##First get the promoter sequences of input genes
	file_name <- sample(letters, 6) %>% paste(., collapse = "")
	seq_file <- paste0(input_dir, file_name, "_seq.fa")
	
	##Get the promoter regions of selected genes
	message(paste0("Extracting promoters from the genome and write out"))
	get_promoter(input_genes = input_genes, upstream = upstream, downstream = downstream, out_file = paste0(input_dir, file_name, ".bed"))
	system_codes <- paste0("bedtools getfasta -fi ", "/home/sm2726/project/reference/macaque/refseq_assembly/ncbi.genomes.2019.06.10/GCF_003339765.1_Mmul_10_genomic.rm_scaffold.fna -bed ", input_dir, file_name, ".bed", " -s -name > ", seq_file)
	system(system_codes)


	##Get the seq to use
	message(paste0("Reading promoters"))
	seq_use <- readDNAStringSet(filepath = seq_file, format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE)
	system(paste0("rm ", input_dir, file_name, ".bed ", seq_file))


	##-----------------------------------------------------------------------------------------------
	##601 pwms, may contain duplicated TFs [Covered in e.build.bg.R script]
	##pwms <- readRDS("/home/sm2726/project/cortex_development/dorsal_ana/module/load_files/5.motif/expressed_motifs.rds")


	## Do the motif enrichment analysis
	message(paste0("Performing motif enrichment"))
	##bg.logn <- readRDS("/home/sm2726/project/cortex_development/dorsal_ana/module/load_files/5.motif/Mmul10_bg.logn.human.algorithm.rds") ## Use the human algorithm
	bg.logn <- readRDS("/home/sm2726/project/NHPfetal/PAT_network/load_files/Mmul10_bg.logn.human.algorithm.20230102.rds")

	## Use parallel computing if there are more than 100 genes
	if (length(seq_use) >= 100){
		message(paste0("Start motif enrichment using ", nCores, " cores"))
		message(Sys.time())
		bin.size <- 50
	    bin.ind <- ceiling(c(1:length(seq_use))/bin.size)
	    max.bin <- max(bin.ind)

		
		cl = makeCluster(nCores, outfile="")
	    doParallel::registerDoParallel(cl);
	    idx <- NULL
	    data_list <- foreach(idx = 1:max.bin, .packages = c("PWMEnrich")) %dopar% {
	    	subres <- motifEnrichment(seq_use[bin.ind == idx], bg.logn, verbose = FALSE, group.only = FALSE, B = 10000)
		    subdata <- list(data = subres$sequence.bg, score = subres$sequence.norm)
			print(paste0("Finish ", toString(max(which(bin.ind == idx))), " sequences"))
		    return(subdata)
		}
		col_fix <- colnames(data_list[[1]][["data"]])
		data <- lapply(data_list, function(x) x[["data"]][, col_fix, drop = FALSE]) %>%
					setNames(., NULL) %>%
					do.call(rbind, .)
		score <- lapply(data_list, function(x) x[["score"]][, col_fix, drop = FALSE]) %>%
					setNames(., NULL) %>%
					do.call(rbind, .)
		stopCluster(cl)
		message("Finish motif enrichment analysis")
		message(Sys.time())
	} else {
		res <- motifEnrichment(seq_use, bg.logn, verbose = FALSE, group.only = FALSE, B = 10000)
		data <- res$sequence.bg
		score <- res$sequence.norm
	}
	print(identical(colnames(data), colnames(score)))

	## Get gene * motif, value is the pval
	rownames(score) <- extract_field(rownames(score), 1, "::") %>% setNames(., NULL)
	score[is.na(score)] <- 1
	rownames(data) <- extract_field(rownames(data), 1, "::") %>% setNames(., NULL)
	data[is.na(data)] <- 1


	## Build the meta data
	tf_gene <- extract_field(colnames(data), "rm_end", "-")
	meta <- data.frame(row.names = colnames(data), gene = extract_field(colnames(data), "rm_end", "-"), stringsAsFactors = FALSE)

	## Subset the meta data
	meta <- meta[meta$gene %in% rownames(data), ,drop = FALSE]  ##please note to use "rownames(data)" here, not input_genes. Some input_genes might be missing.
	data <- data[, rownames(meta)]
	score <- score[, rownames(meta)]

	outs <- list(data = data, meta = meta, score = score)
	return(outs)
}



