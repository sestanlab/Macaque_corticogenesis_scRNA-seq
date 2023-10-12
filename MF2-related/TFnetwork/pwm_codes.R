Build_Motif_DB <- function(input_genes, db_type = c("meme", "matrix")[2], db_dir = NULL, all_motifs = FALSE){
	if (all_motifs){
		tf_motif <- readRDS(file = paste0("/home/sm2726/project/public_data/MEME_db/custom_motif/load_files/", "Human_TF_anno_slim.rds"))
		pfm_list <- readRDS(file = paste0("/home/sm2726/project/public_data/MEME_db/custom_motif/load_files/", "pfm_list.rds"))	

		return(pfm_list)
	}else {
		tf_motif <- readRDS(file = paste0("/home/sm2726/project/public_data/MEME_db/custom_motif/load_files/", "Human_TF_anno_slim.rds"))
		pfm_list <- readRDS(file = paste0("/home/sm2726/project/public_data/MEME_db/custom_motif/load_files/", "pfm_list.rds"))
		shared_genes <- input_genes %>% intersect(., tf_motif$HGNC.symbol) %>% sort()
		tf_motif <- tf_motif[tf_motif$HGNC.symbol %in% shared_genes, ,drop = FALSE]
		if (dim(tf_motif)[1] == 0){
			stop("None of the transcription factors were matched.")
		}

		if (db_type == "matrix"){
			pfm_list <- pfm_list[intersect(names(pfm_list), tf_motif$id)]
		
			if (length(pfm_list) == 0){
			stop("None of the transcription factors were matched.")
			}
			return(pfm_list)
		}
	}

}



