library(igraph) ## network visualization
library(parallel) ## for function mclapply


CalcUMAPbinAvg <- function(object, dims = 1:30, nbins = 30) {
	object <- RunUMAP(object, dims = dims, n.neighbors = 25, n.components = 1L, verbose = FALSE)
	
	qs <- quantile(object$umap@cell.embeddings[, 1], c(0.025, 0.975))
	pdata <- data.frame(UMAP = object$umap@cell.embeddings[, 1], 
					cell = colnames(object),
					stringsAsFactors = FALSE) %>%
				filter(UMAP <= qs[2] & UMAP >= qs[1]) %>%
				mutate(UMAPbin = as.character(as.numeric(cut(UMAP, nbins))))

	## Remove bins with less than 5 cells to avoid variability
	kp_bins <- table(pdata$UMAPbin) %>% .[. >= 15] %>% names()
	pdata <- filter(pdata, UMAPbin %in% kp_bins)


	subobj <- object[, pdata$cell]
	DefaultAssay(subobj) <- "RNA"
	subobj$UMAPbin <- pdata$UMAPbin
	Idents(subobj) <- "UMAPbin"
	avgs <- log(AverageExpression(subobj, assays = "RNA", verbose = FALSE)$RNA + 1) %>%
			as.matrix()
	return(avgs)
}


## Permutate gene expression across cells (disrupting gene-gene correlation) & recalculate average expression
CalcUMAPbin_PermAvg <- function(object, dims = 1:30, nbins = 30, npermutations = 1000, ncores = 12) {
	CalcAvgFromData <- function(data, ident) {
		ident <- ident[colnames(data)]
		all.levs <- levels(as.factor(ident))
		avg.all <- list()
		for (ii in all.levs){
			tmp.cells <- which(ident == ii)
			data.tmp <- data[, tmp.cells]
			if (length(tmp.cells) == 1){
				avg.tmp <- expm1(data.tmp)
			} else if (length(tmp.cells) > 1){
				avg.tmp <- rowMeans(expm1(data.tmp))
			}

			avg.all[[ii]] <- avg.tmp
		}
		avg.final <- log(do.call(cbind, avg.all) + 1)
		return(avg.final)
	}


	## Run UMAP
	object <- RunUMAP(object, dims = dims, n.neighbors = 25, n.components = 1L, verbose = FALSE)
	
	## One UMAP axis bins
	qs <- quantile(object$umap@cell.embeddings[, 1], c(0.025, 0.975))
	pdata <- data.frame(UMAP = object$umap@cell.embeddings[, 1], 
					cell = colnames(object),
					stringsAsFactors = FALSE) %>%
				filter(UMAP <= qs[2] & UMAP >= qs[1]) %>%
				mutate(UMAPbin = as.character(as.numeric(cut(UMAP, nbins))))

	## Remove bins with less than 5 cells to avoid variability
	kp_bins <- table(pdata$UMAPbin) %>% .[. >= 15] %>% names()
	pdata <- filter(pdata, UMAPbin %in% kp_bins)


	## Permutation on the subset data
	subobj <- object[, pdata$cell]
	subobj$UMAPbin <- pdata[colnames(subobj), "UMAPbin"]

	raw_data <- subobj$RNA@data
	raw_ident <- setNames(subobj@meta.data$UMAPbin, colnames(subobj))

	perm_res <- mclapply(1:npermutations, function(idx) {
		set.seed(idx)
		perm_data <- apply(as.matrix(raw_data), 1, sample) %>% t()
		colnames(perm_data) <- colnames(raw_data)
		perm_avg <- CalcAvgFromData(data = perm_data, ident = raw_ident)
		return(perm_avg)
		}, mc.cores = ncores)
	return(perm_res)
}

##-------------------------------------------------------------------------------------------

## Identify regulons

##-------------------------------------------------------------------------------------------
GetInitialRegulons <- function(mot_res, pval_thre = 0.01, score_thre = 1.3) { 
	motifs <- rownames(mot_res$meta)
	tfs <- unique(mot_res$meta$gene)
	pval_data <- as.matrix(mot_res$data)
	logp_data <- -log10(pval_data)
	score_data <- as.matrix(mot_res$score)
	sub_regu <- lapply(tfs, function(tf_name) {
			motifs <- rownames(mot_res$meta)[mot_res$meta$gene %in% tf_name]
			sign_idx <- lapply(motifs, function(motif) {
				tst <- (pval_data[, motif] <= pval_thre) | (score_data[, motif] >= score_thre)
				return(tst)
				}) %>%
						Reduce("|", .)
			TF <- setNames(1e-5, tf_name)

			if (sum(sign_idx) >= 1){
				targets <- score_data[sign_idx, motifs, drop = FALSE] %>% 
							apply(., 1, max) %>%
							setNames(., rownames(score_data)[sign_idx]) %>% 
							sort(decreasing = TRUE)
				return(c(TF, targets))
			}else {
				return(TF)
			}
			}) %>% setNames(., tfs)
	return(sub_regu)
}



## Filter regulons based on 1. significant correlation (pval based on permutation). 2. correlation > 0.1
FilterRegulon_Correlation <- function(regulon, avg, perm_avg, cor_p_thre = 1.3, cor_thre = 0.25) {
	regulon <- lapply(regulon, function(x) {
		TF <- names(x)[x < 0.01] ## TF has a value of 1e-5
		targets <- names(x)[x > 0.01]
		avg <- t(as.matrix(avg))

		## Correlation in the real data
		cors <- cor(avg[, TF, drop = FALSE], avg[, targets, drop = FALSE], method = "p")
		cors <- cors[1, ]

		## Distribution of correlations in the permuted data
		perm_cors <- lapply(1:length(perm_avg), function(idx) {
			pavg <- t(as.matrix(perm_avg[[idx]]))
			pcors <- cor(pavg[, TF, drop = FALSE], pavg[, targets, drop = FALSE], method = "p")
			return(pcors)
			}) %>%
			do.call(rbind, .) %>%
			as.matrix()

		## Calculate p values
		## distribution(permuted & actual), (sum(permuted > actual) + 1)/(npermutations + 1)
		mlogp <- sapply(1:ncol(perm_cors), function(idx) {
			pp <- (sum(perm_cors[, idx] > cors[idx]) + 1)/(nrow(perm_cors) + 1)
			pp <- -log10(pp)
			return(pp)
			}) %>%
			setNames(., names(cors))

		mlogp <- mlogp[mlogp >= cor_p_thre & cors >= cor_thre] ## significance of 0.05
		final <- setNames(1e-5, names(x)[1]) %>%
				c(., mlogp)
		return(final)
	})
	return(regulon)
}












IdentifyRegulon <- function(mar_list, upstream = 800, downstream = 100, pval_thre = 0.1, min_size = 3, file_name, input_dir = new_inputdir, force.recalc = FALSE, nCores = 12) {
    all_genes <- unlist(mar_list) %>% unlist() %>% unique()
    step1_file <- paste0(input_dir, file_name, "_step1motres.rds")
    if ((!file.exists(step1_file)) || force.recalc) {
        mot_res <- MotifEnrich(input_genes = all_genes, input_dir = input_dir, upstream = upstream, downstream = downstream, nCores = nCores)
        saveRDS(mot_res, file = step1_file)
    } else {
        mot_res <- readRDS(file = step1_file)
    }
    

    mot_res <- GetRegulon.marker(mot_res = mot_res, pval_thre = pval_thre, mar_list = mar_list, min_size = min_size)
    return(mot_res)
}



GetRegulon.marker <- function(mot_res, mar_list, pval_thre = 0.01, min_size = 2) {
	## First, extract the data objects
	motifs <- rownames(mot_res$meta)
	pval_data <- mot_res$data ##%>%
					##apply(., 2, function(x) p.adjust(x, method = "fdr"))
	logp_data <- -log10(pval_data)



	## Get the regulons by cluster (names of mar_list) 
	regulons <- lapply(names(mar_list), function(cls) {
		sub_regu <- lapply(motifs, function(motif) {
			tf_name <- mot_res$meta$gene[which(rownames(mot_res$meta) == motif)]
			sign_idx <- (pval_data[, motif] <= pval_thre) & c(rownames(pval_data) %in% mar_list[[cls]])
			TF <- setNames(1e-5, tf_name)

			## Check if the TF is in the marker list
			TFinMar <- tf_name %in% mar_list[[cls]]


			if (sum(sign_idx) >= 1 && TFinMar){
				targets <- logp_data[sign_idx, motif] %>% setNames(., rownames(logp_data)[sign_idx]) %>% 
								sort(decreasing = TRUE)
				return(c(TF, targets))
			}else {
				return(TF)
			}


			}) %>% setNames(., motifs)
		sub_regu
		}) %>%
			setNames(., names(mar_list))


	## Merge regulons across clusters
	mregu <- lapply(1:length(regulons[[1]]), function(x) {
		tf_name <- regulons[[1]][[x]][1] %>% names()
		newdf <- lapply(regulons, function(y) y[[x]]) %>%
					lapply(., function(y) {
							if (length(y) == 1){
								df <- data.frame(gene = character(), 
											logp = numeric(), 
											stringsAsFactors = FALSE)
							} else {
								df <- data.frame(gene = names(y), 
											logp = y, 
											stringsAsFactors = FALSE)
							}
							return(df)
							}) %>%
					do.call(rbind, .)

		mdf <- newdf %>%
				group_by(gene) %>%
				summarize(logp = max(logp)) ## If TF is also present as a target, it will be discarded
		regu <- setNames(mdf$logp, mdf$gene)
		regu <- c(setNames(1e5, tf_name), regu[setdiff(names(regu), tf_name)])
		return(regu)
		}) %>%
			setNames(., names(regulons[[1]]))


	FilterDupMotifs <- function(x) {
		tfs <- extract_field(names(x), 1, "-")
		non_dupidx <- which(remove_duplicates(tfs, return_index = TRUE))
		dup_tfs <- tfs[duplicated(tfs)] %>% unique()
		keep_idx <- c()
		for (g in dup_tfs){
			idx <- which(tfs == g)
			cur_len <- sapply(idx, function(xxx) length(x[[xxx]]))
			keep_idx <- c(keep_idx, idx[cur_len == max(cur_len)][1]) ## In case there are TFs with same regulon length
		}
		return(x[c(non_dupidx, keep_idx)])
	}

	## Remove only regulons with only TFs
	sign_idx <- sapply(mregu, function(x) length(x) >= min_size)
	if (sum(sign_idx) >= 1){
		## Get the regulons & Remove the duplicated TFs.
		mregu <- mregu[sign_idx] %>% 
						FilterDupMotifs()
		outs <- list(data = pval_data[, names(mregu),drop = FALSE], 
					meta = mot_res$meta[names(mregu), ,drop = FALSE], 
					logp_data = logp_data[, names(mregu),drop = FALSE],
					regulons = mregu)
	} else {
		outs <- NULL
	}
	return(outs)
}





## Build the TF-targets network for each regulon
Regulon2Net <- function(regulons, subset_tfs = NULL) {
	## Subset the regulons based on the provided genes(TF)
	if (!is.null(subset_tfs)){
		regulons <- regulons[grep(paste(subset_tfs, collapse = "|"), names(regulons), ignore.case = TRUE)]
	}


	## Edge data
	link_data <- lapply(regulons, function(x) {
		interactions <- x[-1]
		motif_edges <- data.frame(to = names(interactions), weight = interactions, stringsAsFactors = FALSE) %>%
				mutate(from = names(x[1])) 
		motif_edges[, c("from", "to", "weight")]
		}) %>% do.call(rbind, .)


	## Get the nodes meta data
	all_genes <- lapply(regulons, function(x) names(x)) %>% unlist() %>% unique()
	TFs <- lapply(regulons, function(x) names(x)[1]) %>% unlist() %>% unique()
	nodes_meta <- data.frame(row.names = all_genes, gene = all_genes, gtype = rep("target", length(all_genes)), stringsAsFactors = FALSE)
	nodes_meta$gtype[rownames(nodes_meta) %in% TFs] <- "TF"


	netpre_res <- list(links = link_data, meta = nodes_meta)
	return(netpre_res)
}


##-------------------------------------------------------------------------------------------

## Regulon visualization

##-------------------------------------------------------------------------------------------
PrepareNet <- function(regulon) {
	raw_links <- lapply(names(regulon), function(ctp) {
		ctp_regu <- regulon[[ctp]]
		ctp_link <- lapply(ctp_regu, function(xx) {
			regu <- xx[-1]
			edges <- data.frame(from = names(xx[1]), to = names(regu), weight = regu, stringsAsFactors = FALSE) %>%
					mutate(subtype = ctp)
			return(edges)
			}) %>%
			do.call(rbind, .)
		return(ctp_link)
		}) %>%
		do.call(rbind, .)

	## Update markers (for pie colors)
	nmars <- lapply(names(regulon), function(ctp) {
		union(raw_links$from[raw_links$subtype %in% ctp], raw_links$to[raw_links$subtype %in% ctp])
		}) %>%
		setNames(., names(regulon))
	
	## Merged links
	mg_links <- raw_links %>%
			group_by(from, to) %>%
			summarize(weight = mean(weight)) %>%
			ungroup()

	## remove self loops
	nonloop_idx <- which(mg_links$from != mg_links$to)
	mg_links <- mg_links[nonloop_idx, , drop = FALSE]


	all_tfs <- lapply(regulon, function(x) names(x)) %>% unlist() %>% unique()
	all_genes <- unlist(nmars) %>% unique()
	
	## Node meta
	node_meta <- data.frame(row.names = all_genes, 
				gene = all_genes, gtype = rep("target", length(all_genes)), 
				stringsAsFactors = FALSE)
	node_meta$gtype[rownames(node_meta) %in% all_tfs] <- "TF"

	netpre_res <- list(links = mg_links, meta = node_meta, markers = nmars)
	return(netpre_res)
}


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


library(igraph) ## network visualization
library(parallel) ## for function mclapply



List2Table <- function(gset, input_genes) {
	## convert gene list to a matrix (0/1 value), rows as genes, column as list names
	gtable <- lapply(names(gset), function(x) {
		df <- data.frame(gene = gset[[x]], xx = 1, stringsAsFactors = FALSE) %>% 
				setNames(., c("gene", x))
		return(df)
		}) %>%
			purrr::reduce(full_join, by = "gene") %>%
			replace(is.na(.), 0) %>%
			tibble::column_to_rownames("gene")


	## some input genes might not be annotated, so list them as 0
	mis_genes <- setdiff(input_genes, rownames(gtable))
	if (length(mis_genes) > 0){
			gtable_new <- matrix(0, nrow = length(mis_genes), ncol = ncol(gtable), 
								dimnames = list(mis_genes, colnames(gtable))) %>%
					as.data.frame() %>%
					rbind(., gtable) %>%
					.[input_genes, ,drop = FALSE]
	} else {
			gtable_new <- gtable[input_genes, ] 
	}
	return(gtable_new)
}



plot_tfnet_signal_final <- function(netpre_object, 
	vertice_size = c(0.5, 2.5, 5.5),
	label.scale = c(1.5, 1.5), 
	plot.scale = 1, 
	vertice_anno_gset, label_genes, label_TF = TRUE, 
	vertice_highlight, 
	vertice_cate_cols = NULL, 
	edge_hightlight,
	edge_width = c(0.1, 0.4),
	remove_small_vertices = FALSE, remove_TFs = FALSE, 
	redo_layout = FALSE, 
	file_name, output_dir = "./report/") {
	
	## Some parameters
	empty_color <- "#D3D3D3"
	set.seed(42)
	
	## Build the net object
	net <- graph_from_data_frame(d=netpre_object$links, vertices=netpre_object$meta, directed=TRUE)


	## Vertice annotation (pie charts)
	vt_anno_mat <- List2Table(gset = vertice_anno_gset, input_genes = rownames(netpre_object$meta))
	vt_anno_names <- setdiff(colnames(vt_anno_mat), "EMPTY")
	vt_anno_mat$EMPTY <- ifelse(rowSums(as.matrix(vt_anno_mat)) == 0, 1, 0)
	vt_anno_list <- vt_anno_mat %>%
					t() %>% 
					as.data.frame(check.names = FALSE) %>% 
					as.list()


	##------------------------------------------
	## Vertice attributes

	## Shape
	vt.shape <- ifelse(rowSums(vt_anno_mat) >= 2, "pie", "circle")
	V(net)$shape <- vt.shape


	## Color schmes for vertices
	gg_color_hue <- function(n) {
	    hues = seq(15, 375, length = n + 1)
	    hcl(h = hues, l = 65, c = 100)[1:n]
	}
	if (is.null(vertice_cate_cols)){
		vertice_cate_cols <- gg_color_hue(length(vt_anno_names)) %>% 
					setNames(., vt_anno_names)
	} else {
		vertice_cate_cols <- vertice_cate_cols[vt_anno_names]
	}
	vertice_cate_cols <- c(vertice_cate_cols, c(EMPTY = empty_color))
	

	V(net)$pie.color <- list(vertice_cate_cols)
	circle.cols <- ifelse(vt.shape == "circle", 
		sapply(1:length(vt.shape), function(x) vertice_cate_cols[which(vt_anno_list[[x]] == 1)[1]]), 
		"black")
	circle.cols[is.na(circle.cols)] <- empty_color
	V(net)$color <- circle.cols


	## Size
	V(net)$size <- ifelse(V(net)$gtype == "target", vertice_size[1], vertice_size[3])
	V(net)$size[netpre_object$meta$gene %in% vertice_highlight & 
				V(net)$gtype == "target"] <- vertice_size[2]

	## labels
	if (label_TF){
		label_genes <- union(label_genes,
					netpre_object$meta$gene[netpre_object$meta$gtype == "TF"])
	}
	gg_labels <- ifelse(netpre_object$meta$gene %in% label_genes, 
			netpre_object$meta$gene, 
			"")
	V(net)$label <- gg_labels
	V(net)$label.cex <- ifelse(V(net)$gtype == "target", label.scale[1], label.scale[2])

	text_cols <- ifelse(V(net)$gtype == "target", "#808184", "#000000")
	V(net)$label.color <- text_cols

	##------------------------------------------
	# Edge attributes
	ed_hl_idx <- lapply(1:nrow(edge_hightlight), function(x) {
		idx <- which(netpre_object$links$from == edge_hightlight[x, 1] & 
					netpre_object$links$to == edge_hightlight[x, 2])
		return(idx)
		})
	ed_hl_idx <- ed_hl_idx[sapply(ed_hl_idx, length) > 0] %>% unlist() %>% unique()
	print(ed_hl_idx)
	E(net)$color <- ifelse(1:length(E(net)) %in% ed_hl_idx, "#000000", "#808184")
	E(net)$arrow.size <- ifelse(1:length(E(net)) %in% ed_hl_idx, 1, 0.5)
	E(net)$width <- ifelse(1:length(E(net)) %in% ed_hl_idx, edge_width[2], edge_width[1])


	##--------------------------------------------------------------------------
	## remove the vertices with very small sizes
	if (remove_small_vertices){
		rm_idx <- which(V(net)$size == vertice_size[1])
	} else {
		rm_idx <- nrow(netpre_object$meta) + 100000
	}

	## remove some TF regulons for simplicity
	if (!is.null(remove_TFs)){
		targets_of_rmTFs <- netpre_object$links$to[netpre_object$links$from %in% remove_TFs] %>% unique()
		## among the targets of the removed TFs, identify the targets that are only regulated by the removed TFs, but not other TFs
		orphan_targets <- netpre_object$links[netpre_object$links$to %in% targets_of_rmTFs, , drop = FALSE] %>%
						group_by(to) %>%
						summarize(nrmTFs = sum(from %in% remove_TFs), nTFs = n()) %>%
						ungroup() %>%
						filter(nTFs == nrmTFs) %>%
						.$to
		
		## Further make sure, the removed targets are not TFs, or those small nodes
		rm_targets <- intersect(orphan_targets, names(V(net))[V(net)$size == vertice_size[2]])
		print(rm_targets)
		rm_idx <- union(rm_idx, which(names(V(net)) %in% union(remove_TFs, rm_targets)))
	} else {
		rm_idx <- union(rm_idx, nrow(netpre_object$meta) + 100000)
	}


	##--------------------------------------------------------------------------
	## Visualization
	
	## Layout for the network
	if (redo_layout){
			net2 <- induced_subgraph(net, V(net)[-rm_idx])

			e <- get.edgelist(net2, names=FALSE)
			l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net2), 
							area=5*(vcount(net2)^2), repulse.rad=(vcount(net2)^3)) 
			lmeta <- as.data.frame(l)
			lmeta$keep <- ifelse(1:nrow(lmeta) %in% rm_idx, 0, 1)
			saveRDS(lmeta, file = paste0("./load_files/", file_name, "_layout.rds"))
	} else {
			e <- get.edgelist(net, names=FALSE)
			l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net), 
							area=15*(vcount(net)^2), repulse.rad=(vcount(net)^3)) 
			lmeta <- as.data.frame(l)
			lmeta$keep <- ifelse(1:nrow(lmeta) %in% rm_idx, 0, 1)
			saveRDS(lmeta, file = paste0("./load_files/", file_name, "_layout.rds"))
			
			net2 <- induced_subgraph(net, V(net)[-rm_idx])
			l <- l[-rm_idx, ]
	}


	## Output the network
	pdf(paste0(output_dir, file_name, ".pdf"), width = 6 * plot.scale, height = 6 * plot.scale, useDingbats = FALSE)
	par(mar=c(0,0,0,0))
	plot(net2, edge.curved=0, 
			layout=l, 
			main=NULL, 
			vertex.frame.width = 0, 
			vertex.frame.color=NA, 
			vertex.pie = vt_anno_list[-rm_idx],
			vertex.label.dist = 0, 
			rescale = TRUE)
	legend("bottomleft", inset=.02, title="Categories", c(vt_anno_names, "EMPTY"), fill=vertice_cate_cols, cex=1.5, ncol = 4)
	dev.off()
} 


