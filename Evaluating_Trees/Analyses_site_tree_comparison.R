######################################################
## PRUNED VS RECONSTRUCTED TREE COMPARISON ANALYSES ##
######################################################

##Set up##
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)
library(phangorn)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_4184 <- read.tree("Final_4k_Pars9NTBIN_ultraLSD2_root.tree")
phy_borneo_recon <- read.tree("Borneo_Recon_Pars9_ultraLSD2.tree")
phy_panama_recon <- read.tree("Panama_Recon_Pars9_ultraLSD2.tree")
tSDI_clusters <- read.csv("tSDI_raw.csv")

#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))

#Prune site trees
site_borneo_tips <- site_borneo_metadata$db_id
site_panama_tips <- site_panama_metadata$db_id

drop_tips_borneo <- setdiff(phy_metadata$db_id, site_borneo_tips)
drop_tips_panama <- setdiff(phy_metadata$db_id, site_panama_tips)

phy_borneo_prune <- drop.tip(phy_4184, drop_tips_borneo)
phy_panama_prune <- drop.tip(phy_4184, drop_tips_panama)

write.tree(phy_panama_prune, "Borneo_Prune_Pars9_ultraLSD2.tree")
write.tree(phy_panama_prune, "Panama_Prune_Pars9_ultraLSD2.tree")

##Calculations

#Robinson-Foulds metric
B_nrf <- RF.dist(phy_borneo_prune, phy_borneo_recon, normalize = TRUE, check.labels = TRUE, rooted = TRUE)
P_nrf <- RF.dist(phy_panama_prune, phy_panama_recon, normalize = TRUE, check.labels = TRUE, rooted = TRUE)

B_rf <- RF.dist(phy_borneo_prune, phy_borneo_recon, normalize = FALSE, check.labels = TRUE, rooted = TRUE)
P_rf <- RF.dist(phy_panama_prune, phy_panama_recon, normalize = FALSE, check.labels = TRUE, rooted = TRUE)


#Ensemble tSDI

#load required function to calculate the numerator and denominator components for each tSDI calculation
fractional_components <- function(x){
  x*(x-1)
}
#calculate numerator and denominator for each cluster
tSDI_clusters <- tSDI_clusters%>%mutate(across(c(ntips,cluster1:cluster14), fractional_components))
#sum tSDI across clusters for each taxonomic group
tSDI_clusters <- tSDI_clusters %>%
  mutate_at(vars(cluster1:cluster14), ~ . /ntips)%>%mutate(tSDI = rowSums(across(cluster1:cluster14), na.rm = T))
tSDI_clusters$tSDI <- ifelse(tSDI_clusters$cluster1=="NaN", NA, tSDI_clusters$tSDI)   
tSDI_clusters <- tSDI_clusters%>%select(Tree, rank, tSDI)%>%filter(!is.na(tSDI))                                                       
#calculate ensemble tSDI as mean across taxonomic rank
tSDI_ensemble <- tSDI_clusters%>%group_by(Tree, rank)%>%summarise(ensemble = mean(tSDI))

