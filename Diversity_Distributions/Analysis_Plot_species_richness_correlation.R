##################################
## SPECIES RICHNESS CORRELATION ##
##################################


##Set up##
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)
library(ggplot2)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_4184 <- read.tree("Final_4k_Pars9NTBIN_ultraLSD2_root.tree")

#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))

#count number of taxa per family in each site
borneo_family_count <- site_borneo_metadata%>%count(family)%>%rename("ntaxa_borneo"="n")
panama_family_count <- site_panama_metadata%>%count(family)%>%rename("ntaxa_panama"="n")
#combine both sites to get counts for only the families found in both sites
family_count <- merge(borneo_family_count, panama_family_count, by="family")
family_count[family_count==""] <- NA
family_count <- family_count%>%filter(!is.na(family))

#identify 3 most species rich families for exclusion from plot - they will still be included in the correlation calculation
family_count$exclusion <- ifelse(family_count$family%in%c("Staphylinidae", "Chrysomelidae", "Curculionidae"), "Exclude", "Plot_Include")


##Analysis##

#pearsons correlation between paired sites
species_richness_correlation <- cor.test(family_count$ntaxa_panama, family_count$ntaxa_borneo, method=c("pearson"))

#get parameters
round(species_richness_correlation$estimate, 3)
species_richness_correlation$p.value


##Plot##

#remove 3 species rich families from plotting
family_count <- family_count%>%filter(exclusion=="Plot_Include")

plot_species_richness_correlation <- ggplot(family_count, aes(x=ntaxa_panama, y=ntaxa_borneo))+
  geom_point(size=1.2) +
  xlab("Family Species Richness Panama")+
  ylab("Family Species Richness Malaysia")+
  scale_x_continuous(limits = c(0, 52.5),  breaks = c(0, 10, 20, 30, 40, 50)) +
  scale_y_continuous(limits = c(0, 43),   breaks = c(0, 10, 20, 30, 40)) +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(), axis.title=element_text(size=10)) + 
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color = "red", linewidth = 0.4) +
  annotate("text", x = 52.5, y = 0, size = 3, label = expression(paste("r = 0.733,  ", italic(p), " < 0.001")), parse = TRUE, hjust=1)

ggsave("Species_Richness_plot.svg", plot_species_richness_correlation, height=4, width=5.5, limitsize = F)

dev.off()
