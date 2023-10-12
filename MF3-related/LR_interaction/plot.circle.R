






plot_data <- reshape2::melt(as.matrix(pval_use)) %>% 
				setNames(., c("ID", "cluster_pair", "pval")) %>%
				mutate(mlogp = -log10(pval)) %>%
				mutate(mlogp = MinMax(mlogp, min = 0, max = 2.5)) %>%
				mutate(cluster_pair = as.character(cluster_pair)) %>%
				mutate(from = extract_field(cluster_pair, 1, "|")) %>%
				mutate(to = extract_field(cluster_pair, 2, "|"))
plot_data$mlogp[plot_data$mlogp <= -log10(0.05)] <- 0

for (pair in pair_order){
		subdata <- plot_data %>% 
				subset(ID == pair) %>%
				select(from, to, mlogp)
		max_sec <- 10
		max_sum1 <- tapply(subdata$mlogp, subdata$from, sum) %>% max(., na.rm = TRUE)
		max_sum2 <- tapply(subdata$mlogp, subdata$to, sum) %>% max(., na.rm = TRUE)
		if (max(c(max_sum1, max_sum2)) > 10){
			max_sec <- max(c(max_sum1, max_sum2))
		}
		print(max_sec)
		allcls <- c("PC FGF17", "PC NKX2-1", "PC RSPO3", "PC TTR", "FC NERG-early", "FC NERG-early", "OcC NERG-early")
		cls.cols <- c("#4d9221", "#80cdc1", "#74add1", "#fdb863", "#f1b6da", "#dfc27d", "#c51b7d") %>% 
				setNames(., allcls)
		xmax_values <- rep(max_sec, length(allcls)) %>% setNames(., allcls)

		randomcoloR::randomColor(count = 2, hue = "red", luminosity = "bright")

















		pdf(paste0(outputdir, "sep_circle/PAT-RGC_circle_", pair, ".pdf"), width = 6, height = 6)
		circos.par(start.degree = 90, cell.padding = c(0, 0, 0, 0), gap.degree=5.5)
		chordDiagram(subdata, grid.col = cls.cols[allcls], order = allcls, xmax = xmax_values, reduce = -1, col = "black", scale = FALSE, annotationTrack = c("grid"), annotationTrackHeight = convert_height(7, "mm"), transparency = 0.6, directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", link.arr.length= 0.05, link.arr.width = 0.01)#, diffHeight = -mm_h(2))
		for(si in get.all.sector.index()) {
				xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
				ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
				circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1,facing = "bending.inside", niceFacing = TRUE, col = "white", cex = 1.6)
		}
		title(pair, cex = 1)
		circos.clear()
		dev.off()
}
		



