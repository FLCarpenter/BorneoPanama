##########################
## FAMILY TREE PLOTTING ##
##########################


##Set up##
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)
library(ggplot2)
library(ggtree)
library(cowplot)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_4184 <- read.tree("Final_4k_Pars9NTBIN_ultraLSD2_root.tree")

#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))

#Get tip IDs corresponding to 3 most species rich families
staphylinidae_tips <- (phy_metadata%>%filter(family=="Staphylinidae"))$db_id
curculionidae_tips <- (phy_metadata%>%filter(family=="Curculionidae"))$db_id
chrysomelidae_tips <- (phy_metadata%>%filter(family=="Chrysomelidae"))$db_id

#Prune full tree to create 3 separate trees for species rich families 
phy_staphylinidae <- keep.tip(phy_4184, staphylinidae_tips)
phy_curculionidae <- keep.tip(phy_4184, curculionidae_tips)
phy_chrysomelidae <- keep.tip(phy_4184, chrysomelidae_tips)

#Add column to full metadata to show if a sequence is from panama/borneo site (or non-site100)
phy_metadata$Site100 <- ifelse(phy_metadata$db_id%in%site_borneo_metadata$db_id, "Malaysia", ifelse(phy_metadata$db_id%in%site_panama_metadata$db_id, "Panama", "NA"))

#Remove rows to show site membership for heatmap
staphylinidae_heatmap_data <- phy_metadata%>%filter(family=="Staphylinidae")%>%select(db_id, Site100)%>%column_to_rownames(var = "db_id")
curculionidae_heatmap_data <- phy_metadata%>%filter(family=="Curculionidae")%>%select(db_id, Site100)%>%column_to_rownames(var = "db_id")
chrysomelidae_heatmap_data <- phy_metadata%>%filter(family=="Chrysomelidae")%>%select(db_id, Site100)%>%column_to_rownames(var = "db_id")


##Plot##

#Define heatmap colours
site_colours <- c("Malaysia" = "#BB3AD8", "Panama" = "#57D83A", "NA" = "grey90")

#Plot families individually adding heatmap to show which sequences belong to site-samples
plot_staphylinidae <- ggtree(phy_staphylinidae, branch.length='none', layout = "circular", size = 0.25)
heatmap_staphylinidae <- gheatmap(plot_staphylinidae, staphylinidae_heatmap_data, offset = 0.05, width = 0.08, color = NULL, colnames = F) +   
  scale_fill_manual(values = site_colours, breaks = c("Malaysia", "Panama", "NA"), na.value = NA, name="Site") +
  theme(legend.position = "right") + 
  guides(fill=guide_legend(title.position="top", ncol=1)) + 
  theme(plot.margin = unit(c(-0.5, 0, -1, 0), "cm"))


plot_curculionidae <- ggtree(phy_curculionidae, branch.length='none', layout = "circular", size = 0.25)
heatmap_curculionidae <- gheatmap(plot_curculionidae, curculionidae_heatmap_data, offset = 0.05, width = 0.08, color = NULL, colnames = F) +   
  scale_fill_manual(values = site_colours, breaks = c("Malaysia", "Panama", "NA"), na.value = NA, name="Site") +
  theme(legend.position = "right") + 
  guides(fill=guide_legend(title.position="top", ncol=1)) + 
  theme(plot.margin = unit(c(-0.5, 0, -1, 0), "cm"),
        legend.title = element_text(size=13),
        legend.text = element_text(size=11))


plot_chrysomelidae <- ggtree(phy_chrysomelidae, branch.length='none', layout = "circular", size = 0.25)
heatmap_chrysomelidae <- gheatmap(plot_chrysomelidae, chrysomelidae_heatmap_data, offset = 0.05, width = 0.08, color = NULL, colnames = F) +   
  scale_fill_manual(values = site_colours, breaks = c("Malaysia", "Panama", "NA"), na.value = NA, name="Site") +
  theme(legend.position = "right") + 
  guides(fill=guide_legend(title.position="top", ncol=1)) + 
  theme(plot.margin = unit(c(-0.5, 0, -1, 0), "cm"))

#remove individual legends
plot_heatmap_staphylinidae <- plot_grid(heatmap_staphylinidae + theme(legend.position="none"))
plot_heatmap_curculionidae <- plot_grid(heatmap_curculionidae + theme(legend.position="none"))
plot_heatmap_chrysomelidae <- plot_grid(heatmap_chrysomelidae + theme(legend.position="none"))

#extract legend
plot_legend <- get_legend(heatmap_curculionidae + theme(legend.box.margin = margin(-950, 30, 0, -10)))

#plot three trees, adding legend back in
heatmap_SR_plot <- plot_grid(plot_heatmap_staphylinidae, plot_heatmap_curculionidae, plot_heatmap_chrysomelidae, labels = c("(a)", "(b)", "(c)"), vjust = 4, ncol = 1)
heatmap_SR_plot_legend <- plot_grid(heatmap_SR_plot, plot_legend, nrow = 1, rel_widths = c(8,1))

ggsave("Family_Tree_plot.svg", plot = heatmap_SR_plot_legend, width=6, height = 15, limitsize = FALSE)

dev.off()
