########################################
## PHYLOGENETIC STRUCTURE CORRELATION ##
########################################


##Set up##
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)
library(picante)
library(ggplot2)
library(cowplot)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_4184 <- read.tree("Final_4k_Pars9NTBIN_ultraLSD2_root.tree")

#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))

#Set up metadata for picante table, adding column to show site membership 
picante_metadata_lineage <- phy_metadata
picante_metadata_lineage$site100 <- ifelse(picante_metadata_lineage$db_id%in%site_borneo_metadata$db_id, "Borneo", ifelse(picante_metadata_lineage$db_id%in%site_panama_metadata$db_id, "Panama", NA))

#List which families are present in each site
family_presence <- picante_metadata_lineage%>%group_by(family, site100)%>%summarise(n = n())%>%filter(!is.na(site100))%>%filter(!is.na(family))
#List only the families which are present in both sites
family_presence <- family_presence[duplicated(family_presence[,1]) | duplicated(family_presence[,1], fromLast=TRUE),]
family_presence[family_presence==""] <- NA
family_presence <- family_presence%>%filter(!is.na(family))
family_presence <- unique(family_presence$family)


##Load required functions
#1. convert metadata into a picante table showing whether a sequence is present/absent in a site for a lineage
picante_family_lineage <- function(metadata, family_lineage){
  picante_metadata <- metadata%>%filter(family==family_lineage)
  picante_metadata <- picante_metadata%>%select(db_id, site100)%>%count(db_id, site100)
  picante_t <- pivot_wider(picante_metadata, names_from = db_id, values_from=n)
  picante_t$site100 <- as.character(picante_t$site100)
  picante_t <- picante_t%>%filter(!is.na(site100))
  #picante_t$site100[is.na(picante_t$site100) == TRUE] <- "N/A"
  picante_t[is.na(picante_t)] <- 0
  #subset if need to here
  picante_t <- column_to_rownames(picante_t, var = "site100")
  picante_t[picante_t > 0] <- 1
  return(picante_t)
}
#2. calculate standard effect sizes for: Faith's PD; Mean Pairwise Distance; Mean Nearest Taxon Distance
lineage_diversity_metrics <- function(picante_table, phy){
  tree_lineage <- prune.sample(picante_table, phy)
  phylo_diversity <- ses.pd(picante_table, tree_lineage, runs = 999, null.model = "taxa.labels") 
  phylompd <- ses.mpd(picante_table, cophenetic(tree_lineage), runs = 999, null.model = "taxa.labels")
  phylomntd <- ses.mntd(picante_table, cophenetic(tree_lineage), runs = 999, null.model = "taxa.labels")
  diversity_metrics <- data.frame(phylo_diversity, phylompd, phylomntd)
  return(diversity_metrics)
}
#3. calculate phylodiversity metrics for all lineages according to taxonomic rank
lineage_resample <- function(metadata, taxonomic_rank, lineage, phy){
  if (taxonomic_rank=="Superfamily") {
    picante_superfamily <- picante_superfamily_lineage(metadata, lineage)
    metrics_superfamily <- lineage_diversity_metrics(picante_superfamily, phy)
    return(metrics_superfamily)
  } else if(taxonomic_rank=="Family"){
    picante_family <- picante_family_lineage(metadata, lineage)
    metrics_family <- lineage_diversity_metrics(picante_family, phy)
    return(metrics_family)
  }
}


##Analysis##

##Calculating metrics

#create empty list
diversity_list <- vector("list", length(family_presence))
#calculate phylogenetic structure of each site within each family
for (i in family_presence) {
  diversity_list[[i]] <- lineage_resample(picante_metadata_lineage, "Family", i, phy_4184)%>%rownames_to_column("Site")%>%mutate(Family = i)
}
diversity_family <- bind_rows(diversity_list, .id = NULL)

#Calculate Net Relatedness Index and Nearest Taxon Index
diversity_family <- diversity_family%>%mutate(NRI = (-1*(mpd.obs - mpd.rand.mean)/mpd.rand.sd))%>%mutate(NTI = (-1*(mntd.obs - mntd.rand.mean)/mntd.rand.sd))%>%select(Family, Site, ntaxa, pd.obs.z, pd.obs.p, mpd.obs.z, mpd.obs.p, mntd.obs.z, mntd.obs.p, NRI, NTI)

#Label samples by their phylogenetic structure
diversity_family$PDses_result <- ifelse(diversity_family$pd.obs.z > 0 & diversity_family$pd.obs.p >= 0.95, "Overdispersed", ifelse(diversity_family$pd.obs.z < 0 & diversity_family$pd.obs.p <= 0.05, "Clustered", "Random"))
diversity_family$NRI_result <- ifelse(diversity_family$NRI < 0 & diversity_family$mpd.obs.p >= 0.95, "Overdispersed", ifelse(diversity_family$NRI > 0 & diversity_family$mpd.obs.p <= 0.05, "Clustered", "Random"))
diversity_family$NTI_result <- ifelse(diversity_family$NTI < 0 & diversity_family$mntd.obs.p >= 0.95, "Overdispersed", ifelse(diversity_family$NTI > 0 & diversity_family$mntd.obs.p <= 0.05, "Clustered", "Random"))


#Separate each metric, and calculate the inverse of PDses 
diversity_family_PDses <- diversity_family%>%group_by(Family)%>%select(Family, Site, pd.obs.z, PDses_result)%>%pivot_wider(names_from = Site, values_from = pd.obs.z:PDses_result)%>%mutate(Inverse_PDses_Borneo = -1*pd.obs.z_Borneo)%>%mutate(Inverse_PDses_Panama = -1*pd.obs.z_Panama)
diversity_family_NRI <- diversity_family%>%group_by(Family)%>%select(Family, Site, NRI, NRI_result)%>%pivot_wider(names_from = Site, values_from = NRI:NRI_result)
diversity_family_NTI <- diversity_family%>%group_by(Family)%>%select(Family, Site, NTI, NTI_result)%>%pivot_wider(names_from = Site, values_from = NTI:NTI_result)

#Note whether each site has significant phylogenetic structure 
diversity_family_PDses$paired_significance <- ifelse((diversity_family_PDses$PDses_result_Borneo=="Random" & diversity_family_PDses$PDses_result_Panama=="Random"), "Neither", ifelse(((diversity_family_PDses$PDses_result_Borneo=="Random" & diversity_family_PDses$PDses_result_Panama=="Clustered") | (diversity_family_PDses$PDses_result_Borneo=="Random" & diversity_family_PDses$PDses_result_Panama=="Overdispersed") | (diversity_family_PDses$PDses_result_Borneo=="Clustered" & diversity_family_PDses$PDses_result_Panama=="Random") | (diversity_family_PDses$PDses_result_Borneo=="Overdispersed" & diversity_family_PDses$PDses_result_Panama=="Random")), "One", ifelse((diversity_family_PDses$PDses_result_Borneo=="Clustered" & diversity_family_PDses$PDses_result_Panama=="Clustered") | (diversity_family_PDses$PDses_result_Borneo=="Clustered" & diversity_family_PDses$PDses_result_Panama=="Overdispersed") | (diversity_family_PDses$PDses_result_Borneo=="Overdispersed" & diversity_family_PDses$PDses_result_Panama=="Clustered") | (diversity_family_PDses$PDses_result_Borneo=="Overdispersed" & diversity_family_PDses$PDses_result_Panama=="Overdispersed"), "Both", NA)))
diversity_family_NRI$paired_significance <- ifelse((diversity_family_NRI$NRI_result_Borneo=="Random" & diversity_family_NRI$NRI_result_Panama=="Random"), "Neither", ifelse(((diversity_family_NRI$NRI_result_Borneo=="Random" & diversity_family_NRI$NRI_result_Panama=="Clustered") | (diversity_family_NRI$NRI_result_Borneo=="Random" & diversity_family_NRI$NRI_result_Panama=="Overdispersed") | (diversity_family_NRI$NRI_result_Borneo=="Clustered" & diversity_family_NRI$NRI_result_Panama=="Random") | (diversity_family_NRI$NRI_result_Borneo=="Overdispersed" & diversity_family_NRI$NRI_result_Panama=="Random")), "One", ifelse((diversity_family_NRI$NRI_result_Borneo=="Clustered" & diversity_family_NRI$NRI_result_Panama=="Clustered") | (diversity_family_NRI$NRI_result_Borneo=="Clustered" & diversity_family_NRI$NRI_result_Panama=="Overdispersed") | (diversity_family_NRI$NRI_result_Borneo=="Overdispersed" & diversity_family_NRI$NRI_result_Panama=="Clustered") | (diversity_family_NRI$NRI_result_Borneo=="Overdispersed" & diversity_family_NRI$NRI_result_Panama=="Overdispersed"), "Both", NA)))
diversity_family_NTI$paired_significance <- ifelse((diversity_family_NTI$NTI_result_Borneo=="Random" & diversity_family_NTI$NTI_result_Panama=="Random"), "Neither", ifelse(((diversity_family_NTI$NTI_result_Borneo=="Random" & diversity_family_NTI$NTI_result_Panama=="Clustered") | (diversity_family_NTI$NTI_result_Borneo=="Random" & diversity_family_NTI$NTI_result_Panama=="Overdispersed") | (diversity_family_NTI$NTI_result_Borneo=="Clustered" & diversity_family_NTI$NTI_result_Panama=="Random") | (diversity_family_NTI$NTI_result_Borneo=="Overdispersed" & diversity_family_NTI$NTI_result_Panama=="Random")), "One", ifelse((diversity_family_NTI$NTI_result_Borneo=="Clustered" & diversity_family_NTI$NTI_result_Panama=="Clustered") | (diversity_family_NTI$NTI_result_Borneo=="Clustered" & diversity_family_NTI$NTI_result_Panama=="Overdispersed") | (diversity_family_NTI$NTI_result_Borneo=="Overdispersed" & diversity_family_NTI$NTI_result_Panama=="Clustered") | (diversity_family_NTI$NTI_result_Borneo=="Overdispersed" & diversity_family_NTI$NTI_result_Panama=="Overdispersed"), "Both", NA)))
#assign levels to each option of paired significance
diversity_family_PDses$paired_significance <- factor(diversity_family_PDses$paired_significance, levels = c("Both", "One", "Neither"))
diversity_family_NRI$paired_significance <- factor(diversity_family_NRI$paired_significance, levels = c("Both", "One", "Neither"))
diversity_family_NTI$paired_significance <- factor(diversity_family_NTI$paired_significance, levels = c("Both", "One", "Neither"))

#remove NA values, stemming from a singular site representative
diversity_family_PDses <- diversity_family_PDses%>%filter(!is.na(Inverse_PDses_Borneo) & !is.na(Inverse_PDses_Panama))
diversity_family_NRI <- diversity_family_NRI%>%filter(!is.na(NRI_Borneo) & !is.na(NRI_Panama))
diversity_family_NTI <- diversity_family_NTI%>%filter(!is.na(NTI_Borneo) & !is.na(NTI_Panama))

#Add a new column to indicate the 3 most species rich families
diversity_family_PDses$Species_Rich <- ifelse(diversity_family_PDses$Family%in%c("Staphylinidae", "Curculionidae", "Chrysomelidae"), diversity_family_PDses$Family, "NA")
diversity_family_NRI$Species_Rich <- ifelse(diversity_family_NRI$Family%in%c("Staphylinidae", "Curculionidae", "Chrysomelidae"), diversity_family_NRI$Family, "NA")
diversity_family_NTI$Species_Rich <- ifelse(diversity_family_NTI$Family%in%c("Staphylinidae", "Curculionidae", "Chrysomelidae"), diversity_family_NTI$Family, "NA")

##Correlation Analysis

#pearsons correlation between paired sites
PDses_correlate <- cor.test(diversity_family_PDses$Inverse_PDses_Panama, diversity_family_PDses$Inverse_PDses_Borneo, method=c("pearson"))
#fit linear model 
linear_PDses <- lm(Inverse_PDses_Borneo ~ Inverse_PDses_Panama, data=diversity_family_PDses)
#get parameters from the correlation analysis
r_value_PDses <- PDses_correlate$estimate
p_value_PDses <- PDses_correlate$p.value
#get parameters from the regression analysis
summary(linear_PDses)
estimates_PDses <- coefficients(linear_PDses)
coeff_PDses <- as.data.frame(estimates_PDses)
coeff_PDses <- coeff_PDses%>%rownames_to_column("Coefficients")
coeff_PDses[1,1] <- "Intercept"
coeff_PDses[2,1] <- "Slope"
coeff_PDses <- coeff_PDses%>%pivot_wider(names_from = Coefficients, values_from = estimates_PDses)


#pearsons correlation between paired sites
NRI_correlate <- cor.test(diversity_family_NRI$NRI_Panama, diversity_family_NRI$NRI_Borneo, method=c("pearson"))
#fit linear model 
linear_NRI <- lm(NRI_Borneo ~ NRI_Panama, data=diversity_family_NRI)
#get parameters from the correlation analysis
r_value_NRI <- NRI_correlate$estimate
p_value_NRI <- NRI_correlate$p.value
#get parameters from the regression analysis
summary(linear_NRI)
estimates_NRI <- coefficients(linear_NRI)
coeff_NRI <- as.data.frame(estimates_NRI)
coeff_NRI <- coeff_NRI%>%rownames_to_column("Coefficients")
coeff_NRI[1,1] <- "Intercept"
coeff_NRI[2,1] <- "Slope"
coeff_NRI <- coeff_NRI%>%pivot_wider(names_from = Coefficients, values_from = estimates_NRI)


#pearsons correlation between paired sites
NTI_correlate <- cor.test(diversity_family_NTI$NTI_Panama, diversity_family_NTI$NTI_Borneo, method=c("pearson"))
#fit linear model 
linear_NTI <- lm(NTI_Borneo~NTI_Panama, data=diversity_family_NTI)
#get parameters from the correlation analysis
r_value_NTI <- NTI_correlate$estimate
p_value_NTI <- NTI_correlate$p.value
#get parameters from the regression analysis
summary(linear_NTI)
estimates_NTI <- coefficients(linear_NTI)
coeff_NTI <- as.data.frame(estimates_NTI)
coeff_NTI <- coeff_NTI%>%rownames_to_column("Coefficients")
coeff_NTI[1,1] <- "Intercept"
coeff_NTI[2,1] <- "Slope"
coeff_NTI <- coeff_NTI%>%pivot_wider(names_from = Coefficients, values_from = estimates_NTI)


##Plot##

#Define colours to represent phylogenetic structure & shapes to represent the 3 most species-rich families 
Colours_3 <- c("Neither" = "white", "One" = "grey70", "Both" = "black")
Shapes <- c("Chrysomelidae" = 22, "Curculionidae" = 23, "Staphylinidae"=24, "NA" = 21)

Plot_PDses_Correlate <- ggplot(diversity_family_PDses, aes(x=Inverse_PDses_Panama, y=Inverse_PDses_Borneo, fill=paired_significance, shape=Species_Rich))+
  geom_point()+
  scale_fill_manual(values = Colours_3, na.value = NA, name="Significance") +
  scale_shape_manual(values = Shapes, breaks = c('Staphylinidae', 'Curculionidae', 'Chrysomelidae'), name="Families") +
  xlab(expression("-"*PD[SES]~"Panama"))+  
  ylab(expression("-"*PD[SES]~"Malaysia"))+
  scale_x_continuous(limits = c(-3, 9.9), breaks = c(-2, 0, 2, 4, 6, 8))+
  scale_y_continuous(limits = c(-1.9, 8.4), breaks = c(-2, 0, 2, 4, 6, 8)) +
  theme_bw() +
  theme(legend.position = "right", 
        legend.justification = "top",
        panel.grid = element_blank(),
        legend.title = element_text(size=12.5),
        legend.text = element_text(size=10.5),
        axis.title = element_text(size=12.5),
        axis.text = element_text(size=10),
        plot.margin = unit(c(5.5,30,5.5,5.5),"pt"))+ 
  guides(fill = guide_legend(title.position = "top", nrow=3, override.aes = list(shape = 21, size=1.5), order = 2), shape = guide_legend(title.position = "top", title="Species-rich families", nrow=3, order = 1, override.aes = list(size=1.5))) +
  geom_abline(data = coeff_PDses, aes(intercept = Intercept, slope = Slope), linetype="dashed", color = "red", linewidth = 0.4) +
annotate("text", x = 9.9, y = -1.9, size = 3, label = expression(paste("r = 0.558,  ", italic(p), " < 0.001")), parse = TRUE, hjust=1)

Plot_NRI_Correlate <- ggplot(diversity_family_NRI, aes(x=NRI_Panama, y=NRI_Borneo, fill=paired_significance, shape=Species_Rich))+
  geom_point()+
  scale_fill_manual(values = Colours_3, na.value = NA, name="Significance") +
  scale_shape_manual(values = Shapes, breaks = c('Staphylinidae', 'Curculionidae', 'Chrysomelidae'), name="Families") +
  xlab("NRI Panama")+
  ylab("NRI Malaysia")+
  scale_x_continuous(limits = c(-1, 8.6), breaks = c(0, 2, 4, 6, 8))+
  scale_y_continuous(limits = c(-1.9, 7.2), breaks = c(0, 2, 4, 6, 8)) +
  theme_bw() +
  theme(legend.position = "right", 
        legend.justification = "top",
        panel.grid = element_blank(),
        legend.title = element_text(size=12.5),
        legend.text = element_text(size=10.5),
        axis.title = element_text(size=12.5),
        axis.text = element_text(size=10),
        plot.margin = unit(c(5.5,82,5.5,212),"pt"))+ 
  guides(fill = guide_legend(title.position = "top", nrow=3, override.aes = list(shape = 21, size=1.5), order = 2), shape = guide_legend(title.position = "top", title="Species-rich families", nrow=3, order = 1, override.aes = list(size=1.5))) +
  geom_abline(data = coeff_NRI, aes(intercept = Intercept, slope = Slope), linetype="dashed", color = "red", linewidth = 0.4) +
  annotate("text", x = 8.4, y = -1.9, size = 3, label = expression(paste("r = 0.078, ", italic(p), " = 0.622")), parse = TRUE, hjust=1)

Plot_NTI_Correlate <- ggplot(diversity_family_NTI, aes(x=NTI_Panama, y=NTI_Borneo, fill=paired_significance, shape=Species_Rich))+
  geom_point()+
  scale_fill_manual(values = Colours_3, na.value = NA, name="Significance") +
  scale_shape_manual(values = Shapes, breaks = c('Staphylinidae', 'Curculionidae', 'Chrysomelidae'), name="Families") +
  xlab("NTI Panama")+
  ylab("NTI Malaysia")+
  scale_x_continuous(limits = c(-0.5, 8.5), breaks = c(0, 2, 4, 6, 8))+
  scale_y_continuous(limits = c(-2.1, 6.6), breaks = c(-2, 0, 2, 4, 6)) +
  theme_bw() +
  theme(legend.position = "right", 
        legend.justification = "top",
        panel.grid = element_blank(),
        legend.title = element_text(size=12.5),
        legend.text = element_text(size=10.5),
        axis.title = element_text(size=12.5),
        axis.text = element_text(size=10),
        plot.margin = unit(c(5.5,30,5.5,5.5),"pt"))+ 
  guides(fill = guide_legend(title.position = "top", nrow=3, override.aes = list(shape = 21, size=1.5), order = 2), shape = guide_legend(title.position = "top", title="Species-rich families", nrow=3, order = 1, override.aes = list(size=1.5))) +
  geom_abline(data = coeff_NTI, aes(intercept = Intercept, slope = Slope), linetype="dashed", color = "red", linewidth = 0.4) +
  annotate("text", x = 8.2, y = -2.1, size = 3, label = expression(paste("r = 0.582, ", italic(p), " < 0.001")), parse = TRUE, hjust=1)

#combine plots
plot_Lineage_Correlation_row1 <- plot_grid(Plot_NTI_Correlate + theme(legend.position="none"), Plot_PDses_Correlate + theme(legend.position="none"), labels = c("(b)", "(c)"), label_size = 12, hjust = -0.1, vjust=2, align = "v", nrow=1)
plot_Lineage_Correlation <- plot_grid(Plot_NRI_Correlate, plot_Lineage_Correlation_row1, labels = c("(a)", ""), label_size = 12, hjust = -14, vjust=2, align = "v", nrow=2, rel_widths = c(1,5))

ggsave("Lineage_Correlation_plot.svg", plot_Lineage_Correlation, height=9, width=11, limitsize = F)
