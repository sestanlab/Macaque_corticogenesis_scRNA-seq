##source("/home/sm2726/project/common_scripts/private_R_pack/scCodes.R")
##new_inputdir <- paste0(inputdir, "5.motif/")
##new_outputdir <- paste0(outputdir, "5.motifs/")
library(filesstrings)
library(PWMEnrich)
library(glmnet)
library(dplyr)


###################################################################################################

## 1. Find the TF-regulons [Region-specific TF-regulons]

###################################################################################################
##---------------------------------------------------------------------------------------
## Extract the motifs of all the expressed TFs [Intersect expressed TF with tF database]
if (FALSE) {
	## system("cp /home/sm2726/project/public_data/MEME_db/custom_motif/pwm_codes.R ./") 
	source("/home/sm2726/project/public_data/MEME_db/custom_motif/pwm_codes.R")
	load("/home/sm2726/project/cortex_development/dorsal_ana/module/load_files/2.module/1_ratioDorsal_ratio_size.rdata")
	all_genes <- exp_genes
	pwms <- Build_Motif_DB(input_genes  = all_genes, db_type = c("meme", "matrix")[2], db_dir = NULL, all_motifs = FALSE) %>% 
			toPWM(motifs = .)
	saveRDS(pwms, file = "/home/sm2726/project/cortex_development/dorsal_ana/module/load_files/5.motif/expressed_motifs.rds")
	system("cp /home/sm2726/project/cortex_development/dorsal_ana/module/load_files/5.motif/expressed_motifs.rds ./load_files/") 
}



##---------------------------------------------------------------------------------------
## Get a slim version of GTF file
if (FALSE){
	system("cp /gpfs/ycga/project/sestan/sm2726/reference/macaque/refseq_assembly/ncbi.genomes.2019.06.10/GCF_003339765.1_Mmul_10_genomic.filtered_rm_pseudogenes.mitoCDStoExon.gtf ./load_files/") 
	gtf <- rtracklayer::import("./load_files/")
	gtf_df <- as.data.frame(gtf) %>% filter(!type %in% c("CDS", "start_codon", "stop_codon"))

	##only keep the gene line
	gtf_slim <- gtf_df[gtf_df$type == "gene", ]; 
	slim_genes <- gtf_slim$gene_id %>% unique;
	nc_genes <- setdiff(gtf_df$gene_id %>% unique(), slim_genes)
	nc_slim <- gtf_df[gtf_df$gene_id %in% nc_genes & gtf_df$exon_number == "1", ]
	nc_slim <- lapply(nc_genes, function(geneid) {
		sub_df <- nc_slim[nc_slim$gene_id == geneid, ,drop = FALSE]
		sub_df <- sub_df[sub_df$width == max(sub_df$width), ,drop = FALSE]
		return(sub_df)
		}) %>% do.call(rbind, .)
	##manully remove one line of duplication
	nc_slim <- nc_slim[!nc_slim$transcript_id %in% "XR_001442989.2", ]

	gtf_slim <- rbind(gtf_slim, nc_slim)
	saveRDS(gtf_slim, file = paste0("/gpfs/ycga/project/sestan/sm2726/reference/macaque/refseq_assembly/ncbi.genomes.2019.06.10/Mmul10_slim_gtf.rds"))
	##system("cp /gpfs/ycga/project/sestan/sm2726/reference/macaque/refseq_assembly/ncbi.genomes.2019.06.10/Mmul10_slim_gtf.rds ./load_files/")
}



##---------------------------------------------------------------------------------------
## Get the promoters of all the genes in the genome
if (FALSE){
	source("./motif.fun.R")
	##Get the genes from the gtf file to extract the promoter sequences of each gene (as a background)
	promoter_file <- paste0("./load_files/", "Mmul10_promoter.fa")
	promoter_bed <- paste0("./load_files/", "Mmul10_promoter.bed")
	all_genes <- readRDS("/gpfs/ycga/project/sestan/sm2726/reference/macaque/refseq_assembly/ncbi.genomes.2019.06.10/Mmul10_slim_gtf.rds") %>% .$gene_id %>% gsub("_", "-", .)
	get_promoter(input_genes = all_genes, upstream = 2000, downstream = 500, out_file = promoter_bed)
	system_codes <- paste0("bedtools getfasta -fi ", "/home/sm2726/project/reference/macaque/refseq_assembly/ncbi.genomes.2019.06.10/GCF_003339765.1_Mmul_10_genomic.rm_scaffold.fna -bed ", promoter_bed, " -s -name > ", promoter_file)
	system(system_codes)	
}



### Further remove the genes with promoters >= 50% Ns.
if (FALSE) {
	promoter_file <- paste0("./load_files/", "Mmul10_promoter.fa")
	promoters <- readLines(con = promoter_file)

	rm_idx <- c()
	nseqs <- length(promoters)/2
	for (i in 1:nseqs){
		seq_idx <- i * 2
		if (str_count(promoters[seq_idx], regex("N", ignore_case = TRUE)) >= 300){
			rm_idx <- c(rm_idx, c(i * 2 - 1, seq_idx))
		}
	}
	promoters <- promoters[setdiff(1:length(promoters), rm_idx)]
	writeLines(promoters, con = paste0("./load_files/", "Mmul10_filtered_promoter.fa"), sep = "\n", useBytes = FALSE)	
}



## Build the background for each motif 
seq_use <- readDNAStringSet(filepath = "./load_files/Mmul10_filtered_promoter.fa", format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE)
pwms <- readRDS("./load_files/expressed_motifs.rds")


registerCoresPWMEnrich(12)
bg.logn <- makeBackground(motifs = pwms, type = "logn", bg.seq = seq_use, algorithm = "human")
saveRDS(bg.logn, file = paste0("./load_files/", "Mmul10_bg.logn.human.algorithm.20230102.rds"))
registerCoresPWMEnrich(NULL)




