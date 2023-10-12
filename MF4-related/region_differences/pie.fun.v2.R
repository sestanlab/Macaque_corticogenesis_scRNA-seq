library(scatterpie)
library(ggforce)
library(tidyr)
geom_scatterpie_new <- function(mapping = NULL, data, cols, legend_name = "species", scale.expression = TRUE, ...) {
      names(mapping)[match(c("x", "y"), names(mapping))] <- c("x0", "y0")
      yvar <- get_aes_var(mapping, "y0")
      xvar <- get_aes_var(mapping, "x0")
      
      
      col_use <- names(cols)
      df <- data[rowSums(data[, col_use]) > 0, ] %>%
               gather(., "type", "value", !!enquo(col_use)) %>% 
               group_by(!!sym(xvar), !!sym(yvar)) %>%
               mutate(scaleexp = value/mean(value)) %>%
               mutate(scaleexp = MinMax(scaleexp, min = 0.5, max = 2)) %>% 
               ungroup()


      nbks <- 1000
      df$border_color <- cols[as.character(df$type)]      
      ## set the border & fill color:
      if (scale.expression){
         df$fill_color <- mapply(FUN = function(color, value) {
                           return(colorRampPalette(colors = c("grey95", color))(nbks)[value])
                  }, color = df$border_color, value = as.numeric(cut(df$scaleexp, nbks)))
         df$border_color[df$scaleexp < 1.5] <- NA
      } else{
         df$fill_color <- df$border_color
         df$border_color <- NA
      }

      
      ## Factorize the type column
      df$type <- factor(df$type, levels = col_use)
      names(df)[which(names(df) == "type")] = legend_name
      df$r0 <- 0
      return(geom_arc_bar(mapping, data = df, stat = "pie", inherit.aes = FALSE, ...))
}
environment(geom_scatterpie_new) <- environment(geom_scatterpie)



PlotWeightedScatterPie <- function(pie.data, x.col = "cluster", y.col = "gene", r.col = "radius", cls_use = c("vRG.early", "vRG.late", "oRG"), rsf = 2, scale.expression = TRUE, fc_limit = c(-1, 2.5)){ 
	## rsf:radius-scale-factor
	## set the colors 
	cls_cols <- c("#fc19cf", "#a3c10d", "#dbb75c") %>% setNames(., c("vRG.early", "vRG.late", "oRG"))
	pie.plot <- pie.data
	cls_use <- intersect(names(cls_cols), cls_use)
	cls_cols <- cls_cols[cls_use]


	## legend data
	##lg.ra <- pie.plot[, r.col]
	##lg.ra[lg.ra == 0] <- 0.00001

	pie.plot[, x.col] <- MinMax(pie.plot[, x.col], min = fc_limit[1], max = fc_limit[2])
	pie.plot[, y.col] <- MinMax(pie.plot[, y.col], min = fc_limit[1], max = fc_limit[2])
	pie.plot$latetype <- factor(as.character(pie.plot$latetype), levels = c("logFC|oRG", "logFC|vRG late"))
	##pie.plot$label <- ifelse(pie.plot[, x.col] > 1 & pie.plot[, y.col] > 1 | , pie.plot$gene, NA)
	pie.plot$label <- ifelse(rowSums(pie.plot[, names(cls_cols)]) == length(cls_cols), pie.plot$gene, NA)

	## Plots
	p <- ggplot(data = pie.plot) + 
			geom_scatterpie_new(mapping = aes_string(x = x.col, y = y.col, r0 = "r0", r = r.col, amount = "value", fill = "fill_color", color = "border_color"), data = pie.plot, cols = cls_cols, scale.expression = scale.expression) + 
			ggrepel::geom_text_repel(mapping = aes_string(x = x.col, y = y.col, label = "label"), color = "black", seed = 42, size = 4) +
			theme_cowplot() + 
			scale_fill_identity() + 
			scale_color_identity() + 
			##RotatedAxis() + 
			facet_grid(rows = vars(latetype)) +
			coord_fixed(ratio = 1, xlim = fc_limit, ylim = fc_limit) +
			theme(panel.grid.major = element_line(colour="lightgrey", size = 0.2), axis.line = element_line(size = 0.2), axis.ticks = element_line(size = 0.2), axis.title = element_blank(), axis.text = element_text(size = 11)) ##+
      ##geom_scatterpie_legend(lg.ra, x= +1, y=ceiling(length(feature.order)/4), n = 4, labeller = function(x) rsf * x)
   return(p)
}
