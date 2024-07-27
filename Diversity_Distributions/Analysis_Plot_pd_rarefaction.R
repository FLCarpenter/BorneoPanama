########################################
## PHYLOGENETIC DIVERSITY RAREFACTION ##
########################################


##Set up##
rm(list = ls())
options(scipen=1000000)

#Load libraries
library(tidyverse)
library(ape)
library(ggplot2)
library(cowplot)
library(matrixStats)
library(stringi)
library(tibble)
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

#Filter metadata for species-rich families 
staphylinidae_metadata <- phy_metadata%>%filter(family=="Staphylinidae")
staphylinidae_borneo_metadata <- site_borneo_metadata%>%filter(family=="Staphylinidae")
staphylinidae_panama_metadata <- site_panama_metadata%>%filter(family=="Staphylinidae")
curculionidae_metadata <- phy_metadata%>%filter(family=="Curculionidae")
curculionidae_borneo_metadata <- site_borneo_metadata%>%filter(family=="Curculionidae")
curculionidae_panama_metadata <- site_panama_metadata%>%filter(family=="Curculionidae")
chrysomelidae_metadata <- phy_metadata%>%filter(family=="Chrysomelidae")
chrysomelidae_borneo_metadata <- site_borneo_metadata%>%filter(family=="Chrysomelidae")
chrysomelidae_panama_metadata <- site_panama_metadata%>%filter(family=="Chrysomelidae")

#Filter metadata for maximum-coverage dataset
global_coverage_metadata <- phy_metadata%>%filter(db_id%in%phy_482_outgroupRM$tip.label)


##Load required functions
#1. function to give subsample sizes according to specified interval
list_steps <- function(metadata, steps){
  N <- nrow(metadata)
  m <- seq(0,N,steps)
  if(m[length(m)]!=N) {
    m <- c(m,N)}
  return(m)
}
#2. randomly resample metadata according to subsample sizes
resample_without_replacement <- function(metadata, list_steps, steps){
  data <- metadata$db_id
  results=vector("list",length(list_steps))
  #fill first row with 0 (first sample is always 0-length)
  results[[1]] <- 0
  #fill second (e.g. size of sample = 100)
  results[[2]] <- sample(data, steps)
  #if final step is < step i.e. the full tree is larger than the last sample size using steps
  if(list_steps[length(list_steps)]-list_steps[length(list_steps)-1]<steps){
    for (i in c(3:(length(list_steps)-1))) {
      samplesavail <- data[!data %in%results[[i-1]]]
      results[[i]] <- c(results[[i-1]], sample(samplesavail, steps))
      penulitmate_row <- results[[i]]
      remaining_samples <- data[!data %in%results[[i]]]
      all_samples <- c(penulitmate_row, remaining_samples)
      results[[length(list_steps)]] <- all_samples
    }
  } else if(list_steps[length(list_steps)]-list_steps[length(list_steps)-1]==steps){
    for (i in c(3:length(list_steps))) {
      samplesavail <- data[!data %in%results[[i-1]]]
      results[[i]] <- c(results[[i-1]], sample(samplesavail, steps))
    }
  }
  l <- stri_list2matrix(results, byrow=TRUE)
  return(l)
}

#3. convert resampled metadata into a picante table showing whether a sequence is present/absent in each subsample
picante_resample <- function(resample_table){
  colnames(resample_table) <- resample_table[nrow(resample_table), ] 
  resample_table[is.na(resample_table)] <- 0
  resample_table[!resample_table==0] <- 1
  resample_table <- as.data.frame(resample_table)
  resample_table <- resample_table%>%mutate_if(is.character, ~as.numeric(.))
  resample_table <- resample_table%>%mutate(sample_size = rowSums(.))
  resample_table <- resample_table%>%column_to_rownames(var="sample_size")
  return(resample_table)
}

#4. replicate functions 1:3 for a specified number of runs
picante_replicate <- function(metadata, steps, runs){
  sample_size <- list_steps(metadata, steps)
  resample_list <- replicate(runs, resample_without_replacement(metadata, sample_size, steps), simplify = FALSE)
  picante_list <- lapply(resample_list, picante_resample)
  return(picante_list)
}
#5. calculate phylogenetic diversity (Faith's PD) for all replicates
pd_replicate <- function(picante_list, phy){
  pd_resample <- vector("list",length(picante_list))
  for (i in 1:length(picante_list)){
    pd_resample[[i]] <- pd(picante_list[[i]], phy, include.root = T)
  }
  pd_resample_2 <- bind_rows(pd_resample, .id = NULL)
  rownames(pd_resample_2) <- NULL
  return(pd_resample_2)
}

##Analysis##

#Create picante tables indicating presence/absence (0/1) of taxa in different samples after 100 randomisations of resampling procedure (resampled in steps of 10)
picante_global_4175_rarefy <- picante_replicate(coleoptera_metadata, 10, 100)
picante_borneo_rarefy <- picante_replicate(site_borneo_metadata, 10, 100)
picante_panama_rarefy <- picante_replicate(site_panama_metadata, 10, 100)
picante_global_coverage_rarefy <- picante_replicate(global_coverage_metadata, 10, 100)

picante_staphylinidae_rarefy <- picante_replicate(staphylinidae_metadata, 10, 100)
picante_staphylinidae_borneo_rarefy <- picante_replicate(staphylinidae_borneo_metadata, 10, 100)
picante_staphylinidae_panama_rarefy <- picante_replicate(staphylinidae_panama_metadata, 10, 100)
picante_curculionidae_rarefy <- picante_replicate(curculionidae_metadata, 10, 100)
picante_curculionidae_borneo_rarefy <- picante_replicate(curculionidae_borneo_metadata, 10, 100)
picante_curculionidae_panama_rarefy <- picante_replicate(curculionidae_panama_metadata, 10, 100)
picante_chrysomelidae_rarefy <- picante_replicate(chrysomelidae_metadata, 10, 100)
picante_chrysomelidae_borneo_rarefy <- picante_replicate(chrysomelidae_borneo_metadata, 10, 100)
picante_chrysomelidae_panama_rarefy <- picante_replicate(chrysomelidae_panama_metadata, 10, 100)


#Calculate the phylogenetic diversity metric (Faith's PD) for each subsample
pd_global_4175_rarefy  <- pd_replicate(picante_global_4175_rarefy, phy_4184_outgroupRM)
pd_borneo_rarefy  <- pd_replicate(picante_borneo_rarefy, phy_4184_outgroupRM)
pd_panama_rarefy  <- pd_replicate(picante_panama_rarefy, phy_4184_outgroupRM)
pd_global_coverage_rarefy  <- pd_replicate(picante_global_coverage_rarefy, phy_4184_outgroupRM)

pd_staphylinidae_rarefy <- pd_replicate(picante_staphylinidae_rarefy, phy_4184_outgroupRM)
pd_staphylinidae_borneo_rarefy <- pd_replicate(picante_staphylinidae_borneo_rarefy, phy_4184_outgroupRM)
pd_staphylinidae_panama_rarefy <- pd_replicate(picante_staphylinidae_panama_rarefy, phy_4184_outgroupRM)
pd_curculionidae_rarefy <- pd_replicate(picante_curculionidae_rarefy, phy_4184_outgroupRM)
pd_curculionidae_borneo_rarefy <- pd_replicate(picante_curculionidae_borneo_rarefy, phy_4184_outgroupRM)
pd_curculionidae_panama_rarefy <- pd_replicate(picante_curculionidae_panama_rarefy, phy_4184_outgroupRM)
pd_chrysomelidae_rarefy <- pd_replicate(picante_chrysomelidae_rarefy, phy_4184_outgroupRM)
pd_chrysomelidae_borneo_rarefy <- pd_replicate(picante_chrysomelidae_borneo_rarefy, phy_4184_outgroupRM)
pd_chrysomelidae_panama_rarefy <- pd_replicate(picante_chrysomelidae_panama_rarefy, phy_4184_outgroupRM)

#Add in a column to identify samples
pd_global_4175_rarefy <- pd_global_4175_rarefy%>%mutate(Site="Global (All Taxa)")
pd_borneo_rarefy <- pd_borneo_rarefy%>%mutate(Site="Malaysia")
pd_panama_rarefy <- pd_panama_rarefy%>%mutate(Site="Panama")
pd_global_coverage_rarefy <- pd_global_coverage_rarefy%>%mutate(Site="Global (Representative Lineages)")

pd_staphylinidae_rarefy <- pd_staphylinidae_rarefy%>%mutate(Site="Staphylinidae")
pd_staphylinidae_borneo_rarefy <- pd_staphylinidae_borneo_rarefy%>%mutate(Site="Staphylinidae: Malaysia")
pd_staphylinidae_panama_rarefy <- pd_staphylinidae_panama_rarefy%>%mutate(Site="Staphylinidae: Panama")
pd_curculionidae_rarefy <- pd_curculionidae_rarefy%>%mutate(Site="Curculionidae")
pd_curculionidae_borneo_rarefy <- pd_curculionidae_borneo_rarefy%>%mutate(Site="Curculionidae: Malaysia")
pd_curculionidae_panama_rarefy <- pd_curculionidae_panama_rarefy%>%mutate(Site="Curculionidae: Panama")
pd_chrysomelidae_rarefy <- pd_chrysomelidae_rarefy%>%mutate(Site="Chrysomelidae")
pd_chrysomelidae_borneo_rarefy <- pd_chrysomelidae_borneo_rarefy%>%mutate(Site="Chrysomelidae: Malaysia")
pd_chrysomelidae_panama_rarefy <- pd_chrysomelidae_panama_rarefy%>%mutate(Site="Chrysomelidae: Panama")

#Combine samples into one dataframe
pd_replication <- rbind(pd_global_4175_rarefy, pd_borneo_rarefy, pd_panama_rarefy, pd_global_coverage_rarefy, pd_staphylinidae_rarefy, pd_staphylinidae_borneo_rarefy, pd_staphylinidae_panama_rarefy, pd_curculionidae_rarefy, pd_curculionidae_borneo_rarefy, pd_curculionidae_panama_rarefy, pd_chrysomelidae_rarefy, pd_chrysomelidae_borneo_rarefy, pd_chrysomelidae_panama_rarefy)
pd_replication <- pd_replication%>%rename(Sample_Size = SR)
#Calculate the mean phylogenetic diversity value for subsamples across randomisations
pd_curve <- pd_replication%>%group_by(Site, Sample_Size)%>%summarise(mean_PD = mean(PD), sd_PD = sd(PD))%>%mutate(CI95_low = mean_PD - 1.96*(sd_PD/sqrt(Sample_Size)), CI95_high = mean_PD + 1.96*(sd_PD/sqrt(Sample_Size)))
#order levels of samples
pd_curve$Site <- factor(pd_curve$Site, levels = c("Global (All Taxa)", "Global (Representative Lineages)", "Malaysia", "Panama", "Staphylinidae", "Staphylinidae: Malaysia", "Staphylinidae: Panama", "Curculionidae", "Curculionidae: Malaysia", "Curculionidae: Panama", "Chrysomelidae", "Chrysomelidae: Malaysia", "Chrysomelidae: Panama"))


##Plot##

#Define colours
colours <- c("Global (All Taxa)" = "#D86C3A", "Global (Representative Lineages)" = "#3AA6D8", "Malaysia" = "#BB3AD8", "Panama" = "#57D83A", "Staphylinidae" = "#FE8A5F", "Staphylinidae: Malaysia" = "#FFB937", "Staphylinidae: Panama" = "#FFD27B", "Chrysomelidae" = "#018A7B", "Chrysomelidae: Malaysia" = "#01BDA5", "Chrysomelidae: Panama" = "#04E9CB", "Curculionidae" =  "#C00ECD", "Curculionidae: Malaysia" = "#E000E7", "Curculionidae: Panama" = "#EF98FF") 

#x-axis label
xlabel <- expression(N[taxa])

pd_curve_lineages <- pd_curve%>%filter(Site%in%c("Staphylinidae", "Staphylinidae: Malaysia", "Staphylinidae: Panama", "Curculionidae", "Curculionidae: Malaysia", "Curculionidae: Panama", "Chrysomelidae", "Chrysomelidae: Malaysia", "Chrysomelidae: Panama"))
rarefaction_lineages <- ggplot(data = pd_curve_lineages, aes(x=Sample_Size, y=mean_PD, color=factor(Site), fill = factor(Site)))+
  geom_line(linewidth=0.4) +
  scale_color_manual(values = colours, breaks = c("Staphylinidae", "Staphylinidae: Malaysia", "Staphylinidae: Panama", "Curculionidae", "Curculionidae: Malaysia", "Curculionidae: Panama", "Chrysomelidae", "Chrysomelidae: Malaysia", "Chrysomelidae: Panama"), na.value = NA, name="Sample") +
  geom_vline(xintercept = 122, linetype="dashed", color = "black", linewidth = 0.4) +
  ylab("Phylogenetic Diversity (Faith's PD)") +
  xlab(xlabel) +
  scale_y_continuous(n.breaks = 5) +
  xlim(0, 170) +
  ylim(0, 22000) +
  theme_bw() +
  theme(axis.title.y = element_blank(), 
        legend.key.size = unit(0.6, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.box.background = element_rect(),
        panel.grid = element_blank(), 
        legend.position=c(0.001,0.998),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        legend.box.margin = margin(0.3,7.5,0.5,0.5)) 

pd_curve_Coleoptera <- pd_curve%>%filter(Site%in%c("Global (All Taxa)", "Global (Representative Lineages)", "Malaysia", "Panama"))
rarefaction_Coleoptera <- ggplot(data = pd_curve_Coleoptera, aes(x=Sample_Size, y=mean_PD, color=factor(Site), fill = factor(Site)))+
  geom_line(linewidth=0.4) +
  scale_color_manual(values = colours, breaks = c("Global (All Taxa)", "Global (Representative Lineages)", "Malaysia", "Panama"), na.value = NA, name="Sample") +
  geom_vline(xintercept = 473, linetype="dashed", color = "black", linewidth = 0.4) +
  ylab("Phylogenetic Diversity (Faith's PD)") +
  xlab(xlabel) +
  scale_y_continuous(n.breaks = 5) +
  xlim(0, 600) +
  ylim(0, 83000) +
  theme_bw() +
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.box.background = element_rect(),
        panel.grid = element_blank(), 
        legend.position=c(0.001,0.998),
        axis.title = element_text(size=15),
        axis.text = element_text(size=12),
        legend.box.margin = margin(0.3,7.5,0.5,0.5)) 

#combine plots
plot_rarefaction <- plot_grid(rarefaction_Coleoptera + theme(legend.justification = c(0,1)), rarefaction_lineages + theme(legend.justification = c(0,1)), scale =0.98, labels = c("(a)", "(b)"), label_size = 15, hjust = -0.1, vjust=2, align = "h", ncol=2, rel_widths = c(1.03,1))

ggsave("Rarefaction_plot.svg", plot_rarefaction, height=6, width=14, limitsize = F)
