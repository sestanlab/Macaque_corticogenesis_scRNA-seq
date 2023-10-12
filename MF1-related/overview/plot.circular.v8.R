library(circlize)
library(dplyr)
library(tibble)

load(file = "./load_files/MF1_expr_avg_ratio_v2.Rdata")
## avgs, ratios, fmeta
fmeta$subtype <- gsub("InN HMX1", "InN HMX2", fmeta$subtype)
colnames(avgs) <- gsub("InN HMX1", "InN HMX2", colnames(avgs))
colnames(ratios) <- gsub("InN HMX1", "InN HMX2", colnames(ratios))


###reg_cols <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4") %>% setNames(., c("Frontal", "Parietal", "Temporal", "Occipital", "InsulaPCC", "GE"))
###age_cols <- c("#8c510a", "#bf812d", "#dfc27d", "#c7eae5", "#80cdc1", "#35978f", "#01665e") %>% setNames(., c("E37", "E42-43", "E54", "E62-64", "E77-78", "E93", "E110"))

ord_tb <- read.table(paste0("./load_files/", "all.nhp.cbn.v7.txt"), sep="\t", stringsAsFactors=FALSE, header=TRUE)
ord_tb$cluster <- gsub("InN HMX1", "InN HMX2", ord_tb$cluster)
allcls <- ord_tb$cluster

seg_idx <- c(1, 5, 15, 19, 26, 30, 35, 44, 52, 61, 67, 73, 78, 82, 85, 89, 91, 93, 95, 100, 105)
seq_lwd <- rep(0.35, length(seg_idx))


all_gps <- unique(ord_tb$label)
group_colors <- c(`Patterning centers` = "#b03c64", 
					`dorsal NSC` = "#f584e4", #821f44
					enIPC = "#7ca4f9",
					`Excitatory neurons` = "#2166ac",
					`CR` = "#bccf42",
					`GE NSC` = "#f1b6da",
					inIPC = "#7fe63e",
					`Inhibitory neurons` = "#0e9c23",
					gIPC = "#ffc277",
					Astro = "#e08214",
					`OPC-Oligo` = "#ad630a",
					Mes = "#6aada3",
					Immune = "#7a7878",
					`RB & Vas` = "#525759",
					`Early subtypes` = "#440a63")#fa3980
gp_size <- c(`Patterning centers` = 1, 
					`dorsal NSC` = 1.5,
					enIPC = 1.5,
					`Excitatory neurons` = 1.5,
					`CR` = 1.5,
					`GE NSC` = 1.5,
					inIPC = 1.5,
					`Inhibitory neurons` = 1.5,
					gIPC = 1.5,
					Astro = 1.5,
					`OPC-Oligo` = 1.2,
					Mes = 1.5,
					Immune = 0.73,
					`RB & Vas` = 1.5,
					`Early subtypes` = 1.5)
lb_xxx <- c(`Patterning centers` = 2.5, 
					`dorsal NSC` = 10,
					enIPC = 17,
					`Excitatory neurons` = 24.5,
					`CR` = 32.5,
					`GE NSC` = 39.5,
					inIPC = 48.2,
					`Inhibitory neurons` = 65.5,
					gIPC = 80.2,
					Astro = 83.5,
					`OPC-Oligo` = 87.05,
					Mes = 90.05,
					Immune = 92.05,
					`RB & Vas` = 96.5,
					`Early subtypes` = 105.7)
cls_cols <- setNames(group_colors[ord_tb$label], allcls)


hwb <- 0.25

#90f49c
##-------------------------------------------------------------------
## Build the plots
pdf(paste0("./report/", "MF1_circle_v13.pdf"), width = 20, height = 20, useDingbats = FALSE)
circos.par(start.degree = 270, cell.padding = c(0, 0, 0, 0), gap.degree=15, circle.margin = 0.35)
circos.initialize(xlim = c(0, length(allcls)), factors = "fasdf")##factors = "a", 


##------------------------------------------------
## Cluster labels
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,panel.fun = function(x, y){
    text.cex <- rep(1.15, length(allcls))
    circos.text(x = 0:(length(allcls)-1)+0.5, 
        y = rep(0,length(allcls)), font = 2, 
        labels = paste0(1:length(allcls), " ", allcls), niceFacing = TRUE, facing="reverse.clockwise", 
        col= cls_cols, adj=c(1, 0.5), cex=text.cex)
    ll <- rep(3.5,length(seg_idx))
    ##circos.segments(seg_idx, rep(0,length(seg_idx)), seg_idx, ll, lty=2, lwd=seq_lwd)
    circos.text(x = -2.7, y = 0.5, font = 2, cex = 1.5, 
        labels = "Subtype", col="black", facing="bending.inside", niceFacing = TRUE) 
})



##---------------------------------------------------------
## Rectangle: label major cluster groups
gp_csum <- table(ord_tb$label) %>% .[all_gps] %>% cumsum()
gp_start <- c(0, gp_csum[1:(length(gp_csum)-1)]) %>% setNames(., NULL)
gp_end <- gp_csum[1:length(gp_csum)] %>% setNames(., NULL)
print(gp_start)
print(gp_end)
print(ceiling((gp_start + gp_end + 1)/2)-0.5)

TxData <- data.frame(xl = gp_start + 0.5 - hwb, 
                    yb = rep(0, length(all_gps)),
                    xr = gp_end - 0.5 + hwb, 
                    yt = rep(1, length(all_gps)), 
                    color = group_colors, 
                    stringsAsFactors = FALSE)
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05,panel.fun = function(x, y){
    circos.rect(xleft = TxData$xl, ybottom = TxData$yb, 
        xright = TxData$xr, ytop = TxData$yt, border=NA, col=TxData$color)  
    for (ij in 1:length(all_gps)){
    	#circos.text(ceiling((gp_start + gp_end + 1)/2)[ij]-0.5,
    	circos.text(lb_xxx[all_gps[ij]],
                0.5,
        all_gps[ij], cex = gp_size[all_gps[ij]], 
        col="white", facing="bending.inside", niceFacing=TRUE)
    }
    #circos.text(ceiling((gp_start+gp_end)/2)-0.15,
    #            rep(0.5, length(all_gps)),
    #    all_gps, cex = gp_size, 
    #    col="white", facing="bending.inside", niceFacing=TRUE)

    circos.text(x = -2.7, y = 0.5, font = 2, cex = 1.5, 
        labels = "Class", col="black", facing="bending.inside", niceFacing = TRUE) 
})


##---------------------------------------------------------
## Stacked Barplot: age composition
RatioRes <- fmeta %>%
                group_by(cbnage, subtype) %>%
                summarize(ncells = n()) %>%
                ungroup() %>%
                group_by(subtype) %>%
                mutate(clsratio = ncells/sum(ncells)) %>%
                ungroup() 
hwb <- 0.35 ## half width of the bar chart
allages <- c("E37", "E42-43", "E54", "E62-64", "E77-78", "E93", "E110")
age_cols <- c("#8c510a", "#bf812d", "#dfc27d", "#c7eae5", "#80cdc1", "#35978f", "#01665e") %>% 
			setNames(., allages)
RatioData <- lapply(1:length(allcls), function(x){
    subratio <- rep(0, length(allages)) %>% setNames(., rev(allages))
    new_ratio <- RatioRes %>% subset(subtype == allcls[x]) %>% column_to_rownames("cbnage")
    subratio[rownames(new_ratio)] <- new_ratio$clsratio

    df <- data.frame(xl = rep(x - 0.5 - hwb, length(allages)), 
                    yb = c(0, cumsum(subratio[1:(length(allages)-1)])), 
                    xr = rep(x - 0.5 + hwb, length(allages)),
                    yt = cumsum(subratio), 
                    color = rev(age_cols), 
                    cluster = allcls[x], stringsAsFactors = FALSE)
    return(df)
    }) %>% do.call(rbind, .)


circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05, panel.fun = function(x, y){
    circos.rect(xleft = RatioData$xl, ybottom = RatioData$yb, 
        xright = RatioData$xr, ytop = RatioData$yt, border=NA, col=RatioData$color)
    ##circos.yaxis(at= seq(0, 1, 0.25), labels= c("0", "", "0.5", "", "1"))
    circos.segments(seg_idx, rep(0,length(seg_idx)), seg_idx, rep(1,length(seg_idx)), lty=2, lwd=seq_lwd)

    circos.text(x = -2.7, y = 0.5, font = 2, cex = 1.5, 
        labels = "Age %", col="black", facing="bending.inside", niceFacing = TRUE)  
})




##---------------------------------------------------------
## Stacked Barplot: lobe composition
RatioRes <- fmeta %>%
                group_by(lobe, subtype) %>%
                summarize(ncells = n()) %>%
                ungroup() %>%
                group_by(subtype) %>%
                mutate(clsratio = ncells/sum(ncells)) %>%
                ungroup() 
alllbs <- c("GE", "FC", "MSC", "OC", "Insula", "TC")
reg_cols <- c("#B037C4", "#FF420E", "#FFBB00", "#89DA59", "#878787", "#4CB5F5") %>% 
				setNames(., alllbs)
RatioData <- lapply(1:length(allcls), function(x){
    subratio <- rep(0, length(alllbs)) %>% setNames(., rev(alllbs))
    new_ratio <- RatioRes %>% subset(subtype == allcls[x]) %>% column_to_rownames("lobe")
    subratio[rownames(new_ratio)] <- new_ratio$clsratio

    df <- data.frame(xl = rep(x - 0.5 - hwb, length(alllbs)), 
                    yb = c(0, cumsum(subratio[1:(length(alllbs)-1)])), 
                    xr = rep(x - 0.5 + hwb, length(alllbs)),
                    yt = cumsum(subratio), 
                    color = rev(reg_cols), 
                    cluster = allcls[x], stringsAsFactors = FALSE)
    return(df)
    }) %>% do.call(rbind, .)


circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.05, panel.fun = function(x, y){
    circos.rect(xleft = RatioData$xl, ybottom = RatioData$yb, 
        xright = RatioData$xr, ytop = RatioData$yt, border=NA, col=RatioData$color)
    ##circos.yaxis(at= seq(0, 1, 0.25), labels= c("0", "", "0.5", "", "1"))
    circos.segments(seg_idx, rep(0,length(seg_idx)), seg_idx, rep(1,length(seg_idx)), lty=2, lwd=seq_lwd)

    circos.text(x = -2.7, y = 0.5, font = 2, cex = 1.5, 
        labels = "Region %", col="black", facing="bending.inside", niceFacing = TRUE)  
})


##------------------------------------------------
#DotPlot
mars <- list(SOX2 = c(1, 16), FGF17 = c(1, 1), RSPO3 = c(2, 5), PAX6 = c(6, 15), 
			EOMES = c(16, 19), SOX5 = c(20, 26), CUX2 = c(27, 30),
			RELN = c(32, 34), CDO1 = c(36, 44), ASCL1 = c(45, 52), 
			LHX6 = c(53, 61), NR2F2 = c(62, 69), SP8 = c(68, 74), MEIS2 = c(70, 78),
			EGFR = c(79, 82), AQP4 = c(83, 85), PDGFRA = c(86, 89), 
			LUM = c(90, 91), PTPRC = c(92, 93), HBA1 = c(94, 95), 
			FN1 = c(96, 100), DLK1 = c(101, 105))

MinMax <- function (data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
}

markers <- names(mars)
DotAvg_scale <- avgs[markers, allcls] %>%
                as.matrix() %>% 
                t() %>% scale() %>% t() %>% 
                MinMax(., min = -2.5, max = 2.5)
DotAvg <- DotAvg_scale %>%
                reshape2::melt() %>%
                setNames(., c("gene", "cluster", "ExpAvg")) %>%
                mutate(gene = as.character(gene), cluster = as.character(cluster))
DotRatio <- ratios[markers, allcls] %>%
                as.matrix() %>%
                reshape2::melt() %>%
                setNames(., c("gene", "cluster", "ExpRatio")) %>%
                mutate(gene = as.character(gene), cluster = as.character(cluster))
expcols <- colorRampPalette(c("lightgray", "#FF420E"))(100)
DotData <- full_join(DotAvg, DotRatio, by = c("gene", "cluster")) %>%
                group_by(gene) %>%
                mutate(color = expcols[as.numeric(cut(ExpAvg, breaks = 100))]) %>%
                ungroup() %>%
                mutate(cluster = factor(cluster, levels = allcls)) %>%
                mutate(xaxis = as.numeric(cluster) - 0.5) %>%
                mutate(gene = factor(gene, levels = markers)) %>%
                mutate(yaxis = length(markers) - as.numeric(gene) + 0.5)
SegData <- data.frame(xleft = sapply(mars, function(xx) xx[1]), 
					xright = sapply(mars, function(xx) xx[2]),
					gene = factor(names(mars), levels = names(mars)),
					stringsAsFactors = FALSE) %>%
					mutate(yaxis = length(mars) - as.numeric(gene) + 0.5) 


##"#bcbddc"
circos.track(ylim = c(0, length(markers)), bg.border = "black", track.height = 0.2, panel.fun = function(x, y){

    ## Plot the dots
    hg <- length(markers)
    dotSF <- 1.8
    circos.segments(SegData$xleft - 1, SegData$yaxis, 
    				SegData$xright, SegData$yaxis, lty=1, lwd=0.15)
    circos.points(DotData$xaxis, DotData$yaxis, col = DotData$color, pch=16, cex = DotData$ExpRatio * dotSF)
    start_idx <- c(0, seg_idx)
    end_idx <- c(seg_idx, length(allcls))
    ng <- length(markers)
    text_h <- (ng - (1:ng)^(1/1.1))-4
    text_h[(1:length(text_h)) > 15] <- text_h[(1:length(text_h)) > 15] + 2.5
    text_h[1:5] <- text_h[1:5] + 2
    text_h[6:10] <- text_h[6:10] + 1
    text_h[1] <- 20

    circos.text(x = sapply(mars, function(xx) (xx[1]- 1 + xx[2])/2),
        y = text_h, font = 3, cex = 0.9, 
        labels = markers, col="black", facing="bending.inside", niceFacing = TRUE)  

    circos.segments(seg_idx, rep(0,length(seg_idx)), seg_idx, rep(length(markers),length(seg_idx)), lty=2, lwd=seq_lwd)
    ##circos.yaxis(at=seq(1, length(markers), 2)-0.5, labels = rev(markers)[seq(1, length(markers), 2)], labels.font = 3, labels.cex=0.65, side = "left")
    circos.yaxis(at=1:length(markers)-0.5, labels = rev(markers), labels.font = 3, labels.cex=0.42)
    ##circos.yaxis(at=seq(2, length(markers), 2)-0.5, labels = rev(markers)[seq(2, length(markers), 2)], labels.font = 3, labels.cex=0.65, side = "right")


    ## Legend Dots
    hdw <- 1.5
    xbase <- 4
    circos.text(x = xbase+5, y = 0.17 * ng, labels = "% expressed", col="black", facing="bending.inside", niceFacing=TRUE)
    x_loc <- c(xbase - 2*hdw, xbase - 1*hdw, xbase, xbase + 1*hdw)
    circos.points(x = x_loc,
            y = rep(0.25 * ng, 4),
            pch=16, col="black", cex= c(0.25, 0.5, 0.75, 1) * dotSF)
    circos.text(x = x_loc,
            y = rep(0.1 * ng, 4), 
            labels = as.character(c(25, 50, 75, 100)), col="black", facing="bending.inside", niceFacing=TRUE)

    ## Legend Heatmap
    if (TRUE){
	    hh <- 0.035 ##half height of the heatmap legend
	    circos.text(x = 27, y = 0.1 * ng, labels = "Scaled mean expr", col="black", facing="bending.inside", niceFacing=TRUE)
	    x_loc <- seq(16, 20, by = 0.04)
	    circos.rect(xleft = x_loc[1:length(expcols)],
	                ybottom = rep((0.25 - hh) * ng, length(expcols)),
	                xright = x_loc[c(1:length(expcols)) + 1], 
	                ytop = rep((0.25 + hh) * ng, length(expcols)), 
	                border = NA, col = expcols)
	    circos.text(x = 18 + c(0 - 1.5*hdw, 0 + 1.5*hdw),
	            y = rep(0.1 * ng, 2), 
	            labels = as.character(c(min(DotAvg_scale), max(DotAvg_scale)) %>% round(., digits = 2)), col="black", facing="bending.inside", niceFacing=TRUE)	
    }
    
})


##------------------------------------------------
# Cluster size

hwb <- 0.4
CratioData <- data.frame(xl = 1:length(allcls) - 0.5 - hwb, 
                    yb = rep(0, length(allcls)),
                    xr = 1:length(allcls) - 0.5 + hwb, 
                    yt = as.numeric(sqrt(table(fmeta$subtype)[allcls]/nrow(fmeta))), 
                    color = "lightgrey", 
                    stringsAsFactors = FALSE)
circos.track(ylim = c(0, max(CratioData$yt)), bg.border = NA, track.height = 0.06,panel.fun = function(x, y){
    circos.rect(xleft = CratioData$xl, ybottom = CratioData$yb, 
        xright = CratioData$xr, ytop = CratioData$yt, border=NA, col=CratioData$color)
    yy<-seq(0,max(CratioData$yt),length.out=5)
    yy_lab <- paste0(sapply(yy, function(mm) round((mm^2) * 100, digits = 1)), "%")
    yy_lab[1] <- 0
    circos.yaxis(at=yy, labels=yy_lab, labels.cex=0.56)
    circos.segments(seg_idx, rep(0,length(seg_idx)), seg_idx, rep(max(CratioData$yt),length(seg_idx)), lty=2, lwd=seq_lwd)
})





circos.clear()

dev.off()














