library(dplyr)
library(ggplot2)


meta_use <- readRDS(file = paste0("./load_files/Reanno_E37-110.org.meta.10052022.rds"))

meta_use$lobe[meta_use$region %in% "M1C"] <- "FC"
meta_use$lobe[meta_use$lobe %in% "MSC"] <- "PC"



reg_cols <- c("#FF420E", "#FFBB00", "#4CB5F5", "#89DA59", "#878787", "#B037C4") %>% setNames(., c("FC", "PC", "TC", "OC", "Insula", "GE"))
age_cols <- c("#8c510a", "#bf812d", "#dfc27d", "#c7eae5", "#80cdc1", "#35978f", "#01665e") %>% setNames(., c("E37", "E42-43", "E54", "E62-64", "E77-78", "E93", "E110"))


size_data <- meta_use %>%
                group_by(cell_origin) %>%
                summarize(ncells = n(), lobe = unique(lobe), cbnage = unique(cbnage), region = unique(region)) %>%
                ungroup() %>%
                mutate(region = gsub("^FR$", "FC", region)) %>%
                mutate(region = gsub("^MS$", "MSC", region)) %>%
                mutate(region = gsub("^TP$", "TC", region)) %>%
                mutate(region = gsub("^OX$", "OC", region))

## Update E37-43 labels
size_data$region[size_data$region == "FC" & size_data$cbnage %in% c("E37", "E42-43")] <- "Anterior"
size_data$region[size_data$region == "MSC" & size_data$cbnage %in% c("E37", "E42-43")] <- "Dorso-Lateral"
size_data$region[size_data$region == "OC" & size_data$cbnage %in% c("E37", "E42-43")] <- "Posterior"
size_data$region[size_data$region == "GE" & size_data$cbnage %in% c("E37", "E42-43")] <- " GE"


size_data$region <- factor(size_data$region, levels = 
	c(" GE", "Anterior", "Dorso-Lateral", "Posterior", "GE", "FC", "MSC", "OC", "TC", "MGE", "LGE", "CGE", "MFC", "OFC", "DFC", "VFC", "M1C", "S1C", "IPC", "PCC", "V1C", "Insula", "A1C", "STC", "ITC")) 
size_data$lobe <- factor(size_data$lobe, levels = names(reg_cols)) 
size_data$cbnage <- factor(size_data$cbnage, levels = names(age_cols))


p2 <- ggplot(size_data, aes_string(x = "region", y = "ncells")) +
        geom_bar(aes_string(fill = "lobe"), size = 0.2, stat = "identity", color = "black", position = "stack") + 
        theme_bw() + 
        scale_fill_manual(values = reg_cols) +
        facet_grid(cols = vars(cbnage), scales = "free_x", space = "free_x")+
        theme(panel.border=element_blank(), strip.background = element_rect(size = 0.2), 
        	panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        	axis.line=element_line(size = 0.25), axis.ticks = element_line(size = 0.25),
        	axis.text.x = element_text(size = 7.5, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), axis.title.x = element_blank(), axis.title.y = element_text(size = rel(1)), legend.position = "bottom", plot.title = element_blank())

pdf(paste0("./report/", "QC_sample_size.pdf"), width = 8, height = 3.5)
print(p2)
dev.off()




type_ncells <- size_data %>% 
                group_by(cbnage) %>%
                summarize(ncells = sum(ncells)) %>%
                ungroup() %>%
                mutate(y_loc = 40000) %>% 
				mutate(cbnage = factor(cbnage, levels = names(age_cols)))
print(type_ncells)


