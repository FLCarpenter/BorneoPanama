########################################
## PHYLOGENETIC DIVERSITY / STRUCTURE ##
########################################


##Set up##
rm(list = ls())
#options(scipen=1000000)

#Load libraries
library(tidyverse)
library(ape)
library(picante)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_4184 <- read.tree("Final_4k_Pars9NTBIN_ultraLSD2_root.tree")
phy_482 <- read.tree("taxonomy_constraint.tree")

#Remove non-Coleopteran tips from trees
phy_4184_outgroupRM <- drop.tip(phy_4184, c("SRAA00099", "SRAA00097", "GBDL01738", "SRAA00098", "GBDL01734", "SRAA00102", "SRAA00101", "SRAA00103", "GBDL01736"))
phy_482_outgroupRM  <- drop.tip(phy_482, c("SRAA00099", "SRAA00097", "GBDL01738", "SRAA00098", "GBDL01734", "SRAA00102", "SRAA00101", "SRAA00103", "GBDL01736"))


#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))

#Filter metadata for maximum-coverage dataset
global_coverage_metadata <- phy_metadata%>%filter(db_id%in%phy_482_outgroupRM$tip.label)

#Set up metadata for picante table, adding column to show sample membership
picante_metadata_global_4175 <- coleoptera_metadata%>%mutate(site100 = "Global (All Taxa)")
picante_metadata_global_coverage <- coleoptera_metadata
picante_metadata_sites <- coleoptera_metadata
picante_metadata_global_coverage$site100 <- ifelse(picante_metadata_global_coverage$db_id%in%global_coverage_metadata$db_id, "Global (Representative Lineages)", NA)
picante_metadata_sites$site100 <- ifelse(picante_metadata_sites$db_id%in%site_borneo_metadata$db_id, "Malaysia", ifelse(picante_metadata_sites$db_id%in%site_panama_metadata$db_id, "Panama", NA))

picante_metadata_PD <- rbind(picante_metadata_global_4175, picante_metadata_global_coverage, picante_metadata_sites)
picante_metadata_cluster <- rbind(picante_metadata_global_coverage, picante_metadata_sites)

##Load required function
#convert metadata into a picante table showing whether a sequence is present/absent in each sample
picante_table <- function(metadata){
  picante_metadata <- metadata%>%dplyr::select(db_id, site100)%>%dplyr::count(db_id, site100)
  picante_t <- pivot_wider(picante_metadata, names_from = db_id, values_from=n)
  picante_t$site100 <- as.character(picante_t$site100)
  picante_t <- picante_t%>%filter(!is.na(site100))
  picante_t[is.na(picante_t)] <- 0
  picante_t <- column_to_rownames(picante_t, var = "site100")
  picante_t[picante_t > 0] <- 1
  return(picante_t)
}

##Analysis##

##phylodiversity metrics
#create picante table
picante_list <- picante_table(picante_metadata_PD)

#Faith's PD; Mean Pairwise Distance; Mean Nearest Taxon Distance
pd <- pd(picante_list, phy_4184_outgroupRM, include.root = T)
mpd <- mpd(picante_list, cophenetic(phy_4184_outgroupRM))
mntd <- mntd(picante_list, cophenetic(phy_4184_outgroupRM))

phylo_metrics <- cbind(pd, mpd, mntd)

##phylodiversity structure

#create picante table
picante_sites <- picante_table(picante_metadata_cluster)

#Standard effect sizes for: Faith's PD; Mean Pairwise Distance; Mean Nearest Taxon Distance
phylo_diversity <- ses.pd(picante_sites, phy_4184_outgroupRM, runs = 999, null.model = "taxa.labels") 
phylompd <- ses.mpd(picante_sites, cophenetic(phy_4184_outgroupRM), runs = 999,  null.model = "taxa.labels")
phylomntd <- ses.mntd(picante_sites, cophenetic(phy_4184_outgroupRM), runs = 999,  null.model = "taxa.labels")

#Calculate Net Relatedness Index and Nearest Taxon Index 
phylo_nri <- phylompd%>%mutate(NRI = (-1*(mpd.obs - mpd.rand.mean)/mpd.rand.sd))%>%dplyr::select(NRI, mpd.obs.p)
phylo_nti <- phylomntd%>%mutate(NTI = (-1*(mntd.obs - mntd.rand.mean)/mntd.rand.sd))%>%dplyr::select(NTI, mntd.obs.p)

phylo_structure <- cbind(phylo_diversity, phylompd, phylomntd, phylo_nri, phylo_nti)

#Label samples by their phylogenetic structure
phylo_structure$PDses_result <- ifelse(phylo_structure$pd.obs.z > 0 & phylo_structure$pd.obs.p >= 0.95, "Overdispersed", ifelse(phylo_structure$pd.obs.z < 0 & phylo_structure$pd.obs.p <= 0.05, "Clustered", "Random"))
phylo_structure$NRI_result <- ifelse(phylo_structure$NRI < 0 & phylo_structure$mpd.obs.p >= 0.95, "Overdispersed", ifelse(phylo_structure$NRI > 0 & phylo_structure$mpd.obs.p <= 0.05, "Clustered", "Random"))
phylo_structure$NTI_result <- ifelse(phylo_structure$NTI < 0 & phylo_structure$mntd.obs.p >= 0.95, "Overdispersed", ifelse(phylo_structure$NTI > 0 & phylo_structure$mntd.obs.p <= 0.05, "Clustered", "Random"))


