#######################
## LRPD CALCULATIONS ##
#######################


##Set up## 

#Load libraries
suppressMessages(require(picante))
suppressMessages(require(ape))

#Read in data
phy_borneo_prune <- read.tree("Borneo_Prune_Pars9_ultraLSD2.tree")
phy_panama_prune <- read.tree("Panama_Prune_Pars9_ultraLSD2.tree")
LRPD_borneo_picante_files <- list.files(".", pattern="picante_table_borneo")
LRPD_panama_picante_files <- list.files(".", pattern="picante_table_panama")
LRPD_borneo_picante_table_replicates <- lapply(LRPD_borneo_picante_files, read.csv)
LRPD_panama_picante_table_replicates <- lapply(LRPD_panama_picante_files, read.csv)

#set row names
for (i in 1:length(LRPD_borneo_picante_table_replicates)) {
  LRPD_borneo_picante_table_replicates[[i]] <- data.frame(LRPD_borneo_picante_table_replicates[[i]], row.names = 1)
}

for (i in 1:length(LRPD_panama_picante_table_replicates)) {
  LRPD_panama_picante_table_replicates[[i]] <- data.frame(LRPD_panama_picante_table_replicates[[i]], row.names = 1)
}


##Load required function
#calculate phylodiversity metrics at each lineage rank
LRPD_diversity_metrics <- function(picante_table, phy){
  phylo_diversity <- ses.pd(picante_table, phy, runs = 999, null.model = "taxa.labels") 
  phylompd <- ses.mpd(picante_table, cophenetic(phy), runs = 999, null.model = "taxa.labels")
  phylomntd <- ses.mntd(picante_table, cophenetic(phy), runs = 999, null.model = "taxa.labels")
  diversity_metrics <- data.frame(phylo_diversity, phylompd, phylomntd)
  return(diversity_metrics)
}

##Calculate LRPD
#create empty list
LRPD_borneo_replicates <- vector("list", length(LRPD_borneo_picante_table_replicates))
#calculate phylodiversity metrics at each lineage rank, iterate for randomly shuffled replicates
for (i in 1:length(LRPD_borneo_picante_table_replicates)) {
  LRPD_borneo_replicates[[i]] <- LRPD_diversity_metrics(LRPD_borneo_picante_table_replicates[[i]], phy_borneo_prune)
}

#create empty list
LRPD_panama_replicates <- vector("list", length(LRPD_panama_picante_table_replicates))
#calculate phylodiversity metrics at each lineage rank, iterate for randomly shuffled replicates
for (i in 1:length(LRPD_panama_picante_table_replicates)) {
  LRPD_panama_replicates[[i]] <- LRPD_diversity_metrics(LRPD_panama_picante_table_replicates[[i]], phy_panama_prune)
}


write.csv(LRPD_borneo_replicates[[1]], "LRPD_borneo_replicate1.csv")
write.csv(LRPD_borneo_replicates[[2]], "LRPD_borneo_replicate2.csv")
write.csv(LRPD_borneo_replicates[[3]], "LRPD_borneo_replicate3.csv")
write.csv(LRPD_borneo_replicates[[4]], "LRPD_borneo_replicate4.csv")
write.csv(LRPD_borneo_replicates[[5]], "LRPD_borneo_replicate5.csv")
write.csv(LRPD_borneo_replicates[[6]], "LRPD_borneo_replicate6.csv")
write.csv(LRPD_borneo_replicates[[7]], "LRPD_borneo_replicate7.csv")
write.csv(LRPD_borneo_replicates[[8]], "LRPD_borneo_replicate8.csv")
write.csv(LRPD_borneo_replicates[[9]], "LRPD_borneo_replicate9.csv")
write.csv(LRPD_borneo_replicates[[10]], "LRPD_borneo_replicate10.csv")

write.csv(LRPD_panama_replicates[[1]], "LRPD_panama_replicate1.csv")
write.csv(LRPD_panama_replicates[[2]], "LRPD_panama_replicate2.csv")
write.csv(LRPD_panama_replicates[[3]], "LRPD_panama_replicate3.csv")
write.csv(LRPD_panama_replicates[[4]], "LRPD_panama_replicate4.csv")
write.csv(LRPD_panama_replicates[[5]], "LRPD_panama_replicate5.csv")
write.csv(LRPD_panama_replicates[[6]], "LRPD_panama_replicate6.csv")
write.csv(LRPD_panama_replicates[[7]], "LRPD_panama_replicate7.csv")
write.csv(LRPD_panama_replicates[[8]], "LRPD_panama_replicate8.csv")
write.csv(LRPD_panama_replicates[[9]], "LRPD_panama_replicate9.csv")
write.csv(LRPD_panama_replicates[[10]], "LRPD_panama_replicate10.csv")
