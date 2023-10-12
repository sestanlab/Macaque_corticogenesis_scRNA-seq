spdf2arry <- function(df, all_regs = c("FC", "MSC", "TC", "OcC")) {
	cur_regs <- strsplit(colnames(df), "|", fixed = TRUE) %>% 
				unlist() %>% unique() %>% 
				intersect(., all_regs)
	all_pairs <- rep(cur_regs, each = length(cur_regs)) %>%
					paste0(., "|", rep(cur_regs, length(cur_regs)))
	mis_pairs <- setdiff(all_pairs, colnames(df))
	if (length(mis_pairs) > 0){
		mis_df <- matrix(0, nrow = nrow(df), ncol = length(mis_pairs), dimnames = list(rownames(df), mis_pairs))
		df <- cbind(as.matrix(df), mis_df)
	}

	df <- df[, all_pairs] %>%
			as.matrix()

	dexinfo <- lapply(1:nrow(df), function(i) {matrix(df[i, ], nrow = length(cur_regs), ncol = length(cur_regs), byrow = TRUE)}) %>%
			do.call(c, .) %>%
			array(., dim = c(length(cur_regs), length(cur_regs), nrow(df)), dimnames = list(cur_regs, cur_regs, rownames(df)))
	return(dexinfo)
}



SummariseArray_n4 <- function(ary, all_sps) {
	features <- dimnames(ary)[[3]]
	##------------------------------------------------
	## Transfer rSUM & cSUM to features * 4sp
	## 3 models: model 1:species exclusively enriched; mod 2: enriched in 2 species; mod3 depleted in one species
	## Incase mod3 & mod1 have overlaps, put mod1 later will put mod1 in high priority
	## model 2 will not overlap with 

	## For each feature, evaluate the species-erichment pattern
	test_mod3 <- function(mat) {
		csum <- colSums(mat)
		emp_vec <- setNames(rep(0, length(all_sps)), all_sps)
		if (sum(csum == 3) == 1){
			emp_vec[colnames(mat)] <- 1
			emp_vec[colnames(mat)[csum == 3]] <- 0
		}
		return(emp_vec)
	}


	test_mod1 <- function(mat) {
		rsum <- rowSums(mat)
		emp_vec <- setNames(rep(0, length(all_sps)), all_sps)
		if (sum(rsum == 3) == 1){
			emp_vec[rownames(mat)[rsum == 3]] <- 1
		}
		return(emp_vec)
	}


	test_mod2 <- function(mat) {
		rsum <- rowSums(mat)
		ridx <- which(rsum == 2)
		bgidx <- setdiff(1:length(rsum), ridx)
		csum <- colSums(mat)
		emp_vec <- setNames(rep(0, length(all_sps)), all_sps)
		if (sum(rsum == 2) == 2 & sum(csum[bgidx] >= 2) == 2){
			emp_vec[rownames(mat)[rsum == 2]] <- 1
		}
		return(emp_vec)
	}

	all_df <- lapply(features, function(x) {
		mat <- ary[, , x]
		pat_df <- list(data.frame(t(test_mod1(mat)), model = "m1", feature = x, stringsAsFactors = FALSE), 
						data.frame(t(test_mod2(mat)), model = "m2", feature = x, stringsAsFactors = FALSE), 
						data.frame(t(test_mod3(mat)), model = "m3", feature = x, stringsAsFactors = FALSE)) %>%
					do.call(rbind, .)
		return(pat_df)
		}) %>%
			do.call(rbind, .)
	return(all_df)
}

SummariseArray_n3 <- function(ary, all_sps) {
	features <- dimnames(ary)[[3]]
	##------------------------------------------------
	## Transfer rSUM & cSUM to features * 4sp
	## 3 models: model 1:species exclusively enriched; mod 2: enriched in 2 species; mod3 depleted in one species
	## Incase mod3 & mod1 have overlaps, put mod1 later will put mod1 in high priority

	## For each feature, evaluate the species-erichment pattern
	test_mod3 <- function(mat) {
		csum <- colSums(mat)
		emp_vec <- setNames(rep(0, length(all_sps)), all_sps)
		if (sum(csum == 2) == 1){
			emp_vec[colnames(mat)] <- 1
			emp_vec[colnames(mat)[csum == 2]] <- 0
		}
		return(emp_vec)
	}


	test_mod1 <- function(mat) {
		rsum <- rowSums(mat)
		emp_vec <- setNames(rep(0, length(all_sps)), all_sps)
		if (sum(rsum == 2) == 1){
			emp_vec[rownames(mat)[rsum == 2]] <- 1
		}
		return(emp_vec)
	}

	all_df <- lapply(features, function(x) {
		mat <- ary[, , x]
		pat_df <- list(data.frame(t(test_mod1(mat)), model = "m1-n3", feature = x, stringsAsFactors = FALSE), 
						data.frame(t(test_mod3(mat)), model = "m3-n3", feature = x, stringsAsFactors = FALSE)) %>%
					do.call(rbind, .)
		return(pat_df)
		}) %>%
			do.call(rbind, .)
	return(all_df)
}



SummariseArray_n2 <- function(ary, all_sps) {
	features <- dimnames(ary)[[3]]
	all_df <- lapply(features, function(x) {
		mat <- ary[, , x]
		emp_vec <- setNames(rep(0, length(all_sps)), all_sps)
		enr_reg <- rownames(mat)[rowSums(mat) == 1]
		if (length(enr_reg) >= 1){
			emp_vec[enr_reg] <- 1
		}
		pat_df <- data.frame(t(emp_vec), model = "m1-n3", feature = x, stringsAsFactors = FALSE)
		return(pat_df)
		}) %>%
			do.call(rbind, .)
	return(all_df)
}


SummariseArray_Combined <- function(ary, all_sps) {
	if (dim(ary)[[1]] == 2){
		df <- SummariseArray_n2(ary = ary, all_sps = all_sps)
	} else if (dim(ary)[[1]] == 3){
		df <- SummariseArray_n3(ary = ary, all_sps = all_sps)
	} else if (dim(ary)[[1]] == 4){
		df <- SummariseArray_n4(ary = ary, all_sps = all_sps)
	}
	return(df)
}



