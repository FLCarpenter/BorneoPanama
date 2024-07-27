###########################
## LRPD DATA PREPARATION ##
###########################


##Set up##
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_4184 <- read.tree("Final_4k_Pars9NTBIN_ultraLSD2_root.tree")

#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))

##Load required functions
#1. identify which families are of equal species richness in a site
duplicate_list <- function(site_metadata_family){
  dupes <- (site_metadata_family%>%group_by(n)%>%summarise(count=sum(n()))%>%filter(count>1))$n
  return(dupes)
}
#2. select families which have unique species richness in a site
unshuffled <- function(site_metadata_family, duplicate_list){
  site_metadata_family$reshuffle <- ifelse(site_metadata_family$n%in%duplicate_list, NA, site_metadata_family$rank)
  site_metadata_family_no_shuffle <- site_metadata_family%>%filter(!is.na(reshuffle))
  return(site_metadata_family_no_shuffle)
}
#3. reshuffle each set of families with equal species richness within their duplicate set
shuffle_duplicates <- function(duplicate_list, site_metadata_family){
  lm <- list()
  for (i in duplicate_list){
    lm[[as.character(i)]] <- site_metadata_family%>%filter(n==i)%>%mutate(reshuffle=sample(rank))
  }
  shuffle_dupes <- bind_rows(lm, .id = NULL)
  return(shuffle_dupes)
}
#4. combine functions 1:3 to get lineage species-richness rankings
shuffle_all <- function(site_metadata){
  #rank species richness for each family
  site_metadata_family <- site_metadata%>%count(family)%>%mutate_at(c('family'), ~na_if(., ""))%>%filter(!is.na(family))%>%arrange(desc(n))
  site_metadata_family$rank <- 1:length(site_metadata_family$n)
  #shuffle ranks for duplicate richness counts
  dupli <- duplicate_list(site_metadata_family)
  keep_unshuffled <- unshuffled(site_metadata_family, dupli)
  shuffled <- shuffle_duplicates(dupli, site_metadata_family)
  shuffled <- rbind(shuffled, keep_unshuffled)%>%arrange(reshuffle)
  return(shuffled)
}
#5. create list of which families are in each cumulative lineage rank
family_concatenate_list <- function(site_metadata_family_reshuffled){
  family_concatenated <- site_metadata_family_reshuffled%>%mutate(family2=family)%>%transmute(family_concatenated = accumulate(family2, ~ paste(.x, .y, sep="_")))
  family_concatenated_df <- as.data.frame(str_split(family_concatenated$family_concatenated, "_", simplify=T))
  family_concatenated_df <- cbind(family_concatenated, family_concatenated_df)
  family_concatenated_df <- family_concatenated_df%>%mutate_all(~na_if(., ""))
  family_concatenated_df <- family_concatenated_df%>%pivot_longer(!family_concatenated, names_to = "t", values_to = "family")%>%filter(!is.na(family))%>%select(family_concatenated, family)
  return(family_concatenated_df)
}
#6. Filter metadata according to families within each cumulative lineage rank
family_concatenate_metadata <- function(family_concatenated_df, site_metadata){
  list_metadata <- list()
  family_concatenated_unique <- unique(family_concatenated_df$family_concatenated)
  for(j in family_concatenated_unique){
    families_list <- family_concatenated_df%>%filter(family_concatenated==j)
    list_metadata[[j]] <- site_metadata%>%filter(family%in%families_list$family)%>%mutate(family_concatenated=j)
  }
  site_metadata_families <- bind_rows(list_metadata, .id = NULL)
  return(site_metadata_families)
}
#7. convert metadata into a picante table showing whether a sequence is present/absent in each cumulative lineage rank
picante_LRPD <- function(site_family_concatenate_metadata){
  picante_metadata <- site_family_concatenate_metadata%>%select(db_id, family_concatenated)%>%count(db_id, family_concatenated)
  picante_t <- pivot_wider(picante_metadata, names_from = db_id, values_from=n)
  picante_t$family_concatenated <- as.character(picante_t$family_concatenated)
  picante_t <- picante_t%>%filter(!is.na(family_concatenated))
  picante_t[is.na(picante_t)] <- 0
  picante_t <- column_to_rownames(picante_t, var = "family_concatenated")
  picante_t[picante_t > 0] <- 1
  return(picante_t)
}


#Get species-richness rankings for all families, randomly shuffling sets of equal species-richness over 10 iterations
ranked_families_borneo_replicates <- replicate(10, shuffle_all(site_borneo_metadata), simplify = FALSE)
ranked_families_panama_replicates <- replicate(10, shuffle_all(site_panama_metadata), simplify = FALSE)

#Get metadata for each iteration with a column indicating sequence membership within each rank
ranked_family_names_borneo_replicates <- list()
ranked_families_borneo_metadata_replicates <- list()
for (i in 1:length(ranked_families_borneo_replicates)){
  ranked_family_names_borneo_replicates[[i]] <- family_concatenate_list(ranked_families_borneo_replicates[[i]])
  ranked_families_borneo_metadata_replicates[[i]] <- family_concatenate_metadata(ranked_family_names_borneo_replicates[[i]], site_borneo_metadata)
}

ranked_family_names_panama_replicates <- list()
ranked_families_panama_metadata_replicates <- list()
for (i in 1:length(ranked_families_panama_replicates)){
  ranked_family_names_panama_replicates[[i]] <- family_concatenate_list(ranked_families_panama_replicates[[i]])
  ranked_families_panama_metadata_replicates[[i]] <- family_concatenate_metadata(ranked_family_names_panama_replicates[[i]], site_panama_metadata)
}

#create empty list
picante_table_panama_replicates <- vector("list", 10)
picante_table_borneo_replicates <- vector("list", 10)

#create picante presence/absence tables for each replicate
for (i in 1:10){
  picante_table_borneo_replicates[[i]] <- picante_LRPD(ranked_families_borneo_metadata_replicates[[i]])
  picante_table_panama_replicates[[i]] <- picante_LRPD(ranked_families_panama_metadata_replicates[[i]])
}


write.csv(picante_table_borneo_replicates[[1]], "picante_table_borneo_replicate1.csv")
write.csv(picante_table_borneo_replicates[[2]], "picante_table_borneo_replicate2.csv")
write.csv(picante_table_borneo_replicates[[3]], "picante_table_borneo_replicate3.csv")
write.csv(picante_table_borneo_replicates[[4]], "picante_table_borneo_replicate4.csv")
write.csv(picante_table_borneo_replicates[[5]], "picante_table_borneo_replicate5.csv")
write.csv(picante_table_borneo_replicates[[6]], "picante_table_borneo_replicate6.csv")
write.csv(picante_table_borneo_replicates[[7]], "picante_table_borneo_replicate7.csv")
write.csv(picante_table_borneo_replicates[[8]], "picante_table_borneo_replicate8.csv")
write.csv(picante_table_borneo_replicates[[9]], "picante_table_borneo_replicate9.csv")
write.csv(picante_table_borneo_replicates[[10]], "picante_table_borneo_replicate10.csv")

write.csv(picante_table_panama_replicates[[1]], "picante_table_panama_replicate1.csv")
write.csv(picante_table_panama_replicates[[2]], "picante_table_panama_replicate2.csv")
write.csv(picante_table_panama_replicates[[3]], "picante_table_panama_replicate3.csv")
write.csv(picante_table_panama_replicates[[4]], "picante_table_panama_replicate4.csv")
write.csv(picante_table_panama_replicates[[5]], "picante_table_panama_replicate5.csv")
write.csv(picante_table_panama_replicates[[6]], "picante_table_panama_replicate6.csv")
write.csv(picante_table_panama_replicates[[7]], "picante_table_panama_replicate7.csv")
write.csv(picante_table_panama_replicates[[8]], "picante_table_panama_replicate8.csv")
write.csv(picante_table_panama_replicates[[9]], "picante_table_panama_replicate9.csv")
write.csv(picante_table_panama_replicates[[10]], "picante_table_panama_replicate10.csv")
