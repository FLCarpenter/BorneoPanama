###########################################
## LRPD PLOTTING AND BREAKPOINT ANALYSIS ##
###########################################


##Set up##
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)
library(ggplot2)
library(cowplot)
library(segmented)
library(strucchange)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_4184 <- read.tree("Final_4k_Pars9NTBIN_ultraLSD2_root.tree")
LRPD_borneo_files <- list.files(".", pattern="LRPD_borneo")
LRPD_panama_files <- list.files(".", pattern="LRPD_panama")
LRPD_borneo_replicates <- lapply(LRPD_borneo_files, read.csv)
LRPD_panama_replicates <- lapply(LRPD_panama_files, read.csv)

#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))


#Calculate Net Relatedness Index and Nearest Taxon Index
for (i in 1:10){
  LRPD_borneo_replicates[[i]] <- LRPD_borneo_replicates[[i]]%>%mutate(NRI = (-1*(mpd.obs - mpd.rand.mean)/mpd.rand.sd))%>%mutate(NTI = (-1*(mntd.obs - mntd.rand.mean)/mntd.rand.sd))
  LRPD_panama_replicates[[i]] <- LRPD_panama_replicates[[i]]%>%mutate(NRI = (-1*(mpd.obs - mpd.rand.mean)/mpd.rand.sd))%>%mutate(NTI = (-1*(mntd.obs - mntd.rand.mean)/mntd.rand.sd))
}


#separate by metrics
#create empty lists
LRPD_borneo_replicates_PDses <-  vector("list", length(LRPD_borneo_replicates))
LRPD_borneo_replicates_NRI <-  vector("list", length(LRPD_borneo_replicates))
LRPD_borneo_replicates_NTI <-  vector("list", length(LRPD_borneo_replicates))
LRPD_panama_replicates_PDses <-  vector("list", length(LRPD_panama_replicates))
LRPD_panama_replicates_NRI <-  vector("list", length(LRPD_panama_replicates))
LRPD_panama_replicates_NTI <-  vector("list", length(LRPD_panama_replicates))
#select columns relating to each metric
for (i in 1:10){
  LRPD_borneo_replicates_PDses[[i]] <- LRPD_borneo_replicates[[i]]%>%dplyr::select(ntaxa, pd.obs, pd.obs.z, pd.obs.p)
  LRPD_borneo_replicates_NRI[[i]] <- LRPD_borneo_replicates[[i]]%>%dplyr::select(ntaxa, mpd.obs, mpd.obs.z, NRI, mpd.obs.p)
  LRPD_borneo_replicates_NTI[[i]] <- LRPD_borneo_replicates[[i]]%>%dplyr::select(ntaxa, mntd.obs, mntd.obs.z, NTI, mntd.obs.p)
  
  LRPD_panama_replicates_PDses[[i]] <- LRPD_panama_replicates[[i]]%>%dplyr::select(ntaxa, pd.obs, pd.obs.z, pd.obs.p)
  LRPD_panama_replicates_NRI[[i]] <- LRPD_panama_replicates[[i]]%>%dplyr::select(ntaxa, mpd.obs, mpd.obs.z, NRI, mpd.obs.p)
  LRPD_panama_replicates_NTI[[i]] <- LRPD_panama_replicates[[i]]%>%dplyr::select(ntaxa, mntd.obs, mntd.obs.z, NTI, mntd.obs.p)
}


#Rank lineages by species richness 
for (i in 1:10){
  LRPD_borneo_replicates_PDses[[i]] <- LRPD_borneo_replicates_PDses[[i]]%>%arrange(ntaxa)
  LRPD_borneo_replicates_PDses[[i]]$lineage_rank <- seq(1,69,1)
  LRPD_panama_replicates_PDses[[i]] <- LRPD_panama_replicates_PDses[[i]]%>%arrange(ntaxa)
  LRPD_panama_replicates_PDses[[i]]$lineage_rank <- seq(1,63,1)
  
  LRPD_borneo_replicates_NRI[[i]] <- LRPD_borneo_replicates_NRI[[i]]%>%arrange(ntaxa)
  LRPD_borneo_replicates_NRI[[i]]$lineage_rank <- seq(1,69,1)
  LRPD_panama_replicates_NRI[[i]] <- LRPD_panama_replicates_NRI[[i]]%>%arrange(ntaxa)
  LRPD_panama_replicates_NRI[[i]]$lineage_rank <- seq(1,63,1)
  
  LRPD_borneo_replicates_NTI[[i]] <- LRPD_borneo_replicates_NTI[[i]]%>%arrange(ntaxa)
  LRPD_borneo_replicates_NTI[[i]]$lineage_rank <- seq(1,69,1)
  LRPD_panama_replicates_NTI[[i]] <- LRPD_panama_replicates_NTI[[i]]%>%arrange(ntaxa)
  LRPD_panama_replicates_NTI[[i]]$lineage_rank <- seq(1,63,1)
}

#combine replicates into singular dataframes
LRPD_borneo_PDses <- bind_rows(LRPD_borneo_replicates_PDses, .id = "replicate")
rownames(LRPD_borneo_PDses) <- NULL
LRPD_panama_PDses <- bind_rows(LRPD_panama_replicates_PDses, .id = "replicate")
rownames(LRPD_panama_PDses) <- NULL
LRPD_borneo_NRI <- bind_rows(LRPD_borneo_replicates_NRI, .id = "replicate")
rownames(LRPD_borneo_NRI) <- NULL
LRPD_panama_NRI <- bind_rows(LRPD_panama_replicates_NRI, .id = "replicate")
rownames(LRPD_panama_NRI) <- NULL
LRPD_borneo_NTI <- bind_rows(LRPD_borneo_replicates_NTI, .id = "replicate")
rownames(LRPD_borneo_NTI) <- NULL
LRPD_panama_NTI <- bind_rows(LRPD_panama_replicates_NTI, .id = "replicate")
rownames(LRPD_panama_NTI) <- NULL

#calculate the mean value of each phylodiversity metric at each rank level across replicates
LRPD_borneo_PDses_mean <- LRPD_borneo_PDses%>%group_by(lineage_rank)%>%summarise(Mean_PDses=mean(pd.obs.z))
LRPD_panama_PDses_mean <- LRPD_panama_PDses%>%group_by(lineage_rank)%>%summarise(Mean_PDses=mean(pd.obs.z))
LRPD_borneo_NRI_mean <- LRPD_borneo_NRI%>%group_by(lineage_rank)%>%summarise(Mean_NRI=mean(NRI))
LRPD_panama_NRI_mean <- LRPD_panama_NRI%>%group_by(lineage_rank)%>%summarise(Mean_NRI=mean(NRI))
LRPD_borneo_NTI_mean <- LRPD_borneo_NTI%>%group_by(lineage_rank)%>%summarise(Mean_NTI=mean(NTI))
LRPD_panama_NTI_mean <- LRPD_panama_NTI%>%group_by(lineage_rank)%>%summarise(Mean_NTI=mean(NTI))

#calculate the absolute value of each mean phylodiversity metric at each rank level
LRPD_borneo_PDses_mean <- LRPD_borneo_PDses_mean%>%mutate(Site="Malaysia")%>%mutate(Absolute_PDses = sqrt(Mean_PDses^2))
LRPD_panama_PDses_mean <- LRPD_panama_PDses_mean%>%mutate(Site="Panama")%>%mutate(Absolute_PDses = sqrt(Mean_PDses^2))
LRPD_PDses <- rbind(LRPD_borneo_PDses_mean, LRPD_panama_PDses_mean)

LRPD_borneo_NRI_mean <- LRPD_borneo_NRI_mean%>%mutate(Site="Malaysia")%>%mutate(Absolute_NRI = sqrt(Mean_NRI^2))
LRPD_panama_NRI_mean <- LRPD_panama_NRI_mean%>%mutate(Site="Panama")%>%mutate(Absolute_NRI = sqrt(Mean_NRI^2))
LRPD_NRI <- rbind(LRPD_borneo_NRI_mean, LRPD_panama_NRI_mean)

LRPD_borneo_NTI_mean <- LRPD_borneo_NTI_mean%>%mutate(Site="Malaysia")%>%mutate(Absolute_NTI = sqrt(Mean_NTI^2))
LRPD_panama_NTI_mean <- LRPD_panama_NTI_mean%>%mutate(Site="Panama")%>%mutate(Absolute_NTI = sqrt(Mean_NTI^2))
LRPD_NTI <- rbind(LRPD_borneo_NTI_mean, LRPD_panama_NTI_mean)

#Fit linear model
linear_borneo_PDses <- lm(Absolute_PDses ~ lineage_rank, data=LRPD_borneo_PDses_mean)
linear_panama_PDses <- lm(Absolute_PDses ~ lineage_rank, data=LRPD_panama_PDses_mean)

linear_borneo_NRI <- lm(Absolute_NRI ~ lineage_rank, data=LRPD_borneo_NRI_mean)
linear_panama_NRI <- lm(Absolute_NRI ~ lineage_rank, data=LRPD_panama_NRI_mean)

linear_borneo_NTI <- lm(Absolute_NTI ~ lineage_rank, data=LRPD_borneo_NTI_mean)
linear_panama_NTI <- lm(Absolute_NTI ~ lineage_rank, data=LRPD_panama_NTI_mean)


#Estimate the number of breakpoints in the piecewise regression according to BIC criterion
#Kmax was estimated after visual inspection to be an integer larger than the expected no. breakpoints 
selgmented(linear_borneo_PDses, seg.Z = ~lineage_rank, type = "bic", Kmax = 10)
selgmented(linear_panama_PDses, seg.Z = ~lineage_rank, type = "bic", Kmax = 10)

selgmented(linear_borneo_NRI, seg.Z = ~lineage_rank, type = "bic", Kmax = 10)
selgmented(linear_panama_NRI, seg.Z = ~lineage_rank, type = "bic", Kmax = 10)

selgmented(linear_borneo_NTI, seg.Z = ~lineage_rank, type = "bic", Kmax = 10)
selgmented(linear_panama_NTI, seg.Z = ~lineage_rank, type = "bic", Kmax = 10)

#Fit piecewise regression using breakpoints as detected by BIC criterion
piecewise_linear_borneo_PDses <- segmented(linear_borneo_PDses, seg.Z = ~lineage_rank, npsi = 3)
piecewise_linear_panama_PDses <- segmented(linear_panama_PDses, seg.Z = ~lineage_rank, npsi = 6)

piecewise_linear_borneo_NRI <- segmented(linear_borneo_NRI, seg.Z = ~lineage_rank, npsi = 1)
piecewise_linear_panama_NRI <- segmented(linear_panama_NRI, seg.Z = ~lineage_rank, npsi = 2)

piecewise_linear_borneo_NTI <- segmented(linear_borneo_NTI, seg.Z = ~lineage_rank, npsi = 4)
piecewise_linear_panama_NTI <- segmented(linear_panama_NTI, seg.Z = ~lineage_rank, npsi = 6)

#create dataframes of estimated breakpoints for each fitted model
breakpoints_borneo_PDses <- piecewise_linear_borneo_PDses$psi
breakpoints_borneo_PDses <- as.data.frame(breakpoints_borneo_PDses)
breakpoints_borneo_PDses <- breakpoints_borneo_PDses%>%mutate(Site="Malaysia")
breakpoints_panama_PDses <- piecewise_linear_panama_PDses$psi
breakpoints_panama_PDses <- as.data.frame(breakpoints_panama_PDses)
breakpoints_panama_PDses <- breakpoints_panama_PDses%>%mutate(Site="Panama")
breakpoints_PDses <- rbind(breakpoints_borneo_PDses, breakpoints_panama_PDses)

breakpoints_borneo_NRI <- piecewise_linear_borneo_NRI$psi
breakpoints_borneo_NRI <- as.data.frame(breakpoints_borneo_NRI)
breakpoints_borneo_NRI <- breakpoints_borneo_NRI%>%mutate(Site="Malaysia")
breakpoints_panama_NRI <- piecewise_linear_panama_NRI$psi
breakpoints_panama_NRI <- as.data.frame(breakpoints_panama_NRI)
breakpoints_panama_NRI <- breakpoints_panama_NRI%>%mutate(Site="Panama")
breakpoints_NRI <- rbind(breakpoints_borneo_NRI, breakpoints_panama_NRI)

breakpoints_borneo_NTI <- piecewise_linear_borneo_NTI$psi
breakpoints_borneo_NTI <- as.data.frame(breakpoints_borneo_NTI)
breakpoints_borneo_NTI <- breakpoints_borneo_NTI%>%mutate(Site="Malaysia")
breakpoints_panama_NTI <- piecewise_linear_panama_NTI$psi
breakpoints_panama_NTI <- as.data.frame(breakpoints_panama_NTI)
breakpoints_panama_NTI <- breakpoints_panama_NTI%>%mutate(Site="Panama")
breakpoints_NTI <- rbind(breakpoints_borneo_NTI, breakpoints_panama_NTI)

#extract estimated phylodiversity value at each lineage rank
fitted_borneo_PDses <- fitted(piecewise_linear_borneo_PDses)
values_plot_borneo_PDses <- data.frame(lineage_rank = LRPD_borneo_PDses_mean$lineage_rank, PDses = fitted_borneo_PDses)
fitted_panama_PDses <- fitted(piecewise_linear_panama_PDses)
values_plot_panama_PDses <- data.frame(lineage_rank = LRPD_panama_PDses_mean$lineage_rank, PDses = fitted_panama_PDses)
values_plot_borneo_PDses <- values_plot_borneo_PDses%>%mutate(Site="Malaysia")
values_plot_panama_PDses <- values_plot_panama_PDses%>%mutate(Site="Panama")
values_plot_PDses <- rbind(values_plot_borneo_PDses, values_plot_panama_PDses)

fitted_borneo_NRI <- fitted(piecewise_linear_borneo_NRI)
values_plot_borneo_NRI <- data.frame(lineage_rank = LRPD_borneo_NRI_mean$lineage_rank, NRI = fitted_borneo_NRI)
fitted_panama_NRI <- fitted(piecewise_linear_panama_NRI)
values_plot_panama_NRI <- data.frame(lineage_rank = LRPD_panama_NRI_mean$lineage_rank, NRI = fitted_panama_NRI)
values_plot_borneo_NRI <- values_plot_borneo_NRI%>%mutate(Site="Malaysia")
values_plot_panama_NRI <- values_plot_panama_NRI%>%mutate(Site="Panama")
values_plot_NRI <- rbind(values_plot_borneo_NRI, values_plot_panama_NRI)

fitted_borneo_NTI <- fitted(piecewise_linear_borneo_NTI)
values_plot_borneo_NTI <- data.frame(lineage_rank = LRPD_borneo_NTI_mean$lineage_rank, NTI = fitted_borneo_NTI)
fitted_panama_NTI <- fitted(piecewise_linear_panama_NTI)
values_plot_panama_NTI <- data.frame(lineage_rank = LRPD_panama_NTI_mean$lineage_rank, NTI = fitted_panama_NTI)
values_plot_borneo_NTI <- values_plot_borneo_NTI%>%mutate(Site="Malaysia")
values_plot_panama_NTI <- values_plot_panama_NTI%>%mutate(Site="Panama")
values_plot_NTI <- rbind(values_plot_borneo_NTI, values_plot_panama_NTI)


#Calculating thresholds for species-rich and species-poor definitions: Q1 cut off for SR; Q1:Q2 intermediate; Q2:Q4 for SP;
site_borneo_metadata$family[site_borneo_metadata$family==""] <- NA
site_panama_metadata$family[site_panama_metadata$family==""] <- NA
b_counts <- site_borneo_metadata%>%count(family)%>%filter(!is.na(family))%>%arrange(desc(n))
p_counts <- site_panama_metadata%>%count(family)%>%filter(!is.na(family))%>%arrange(desc(n))
as.numeric(quantile(b_counts$n, prob=c(.5), type=1)) #ntaxa 6 -> rank 34
as.numeric(quantile(p_counts$n, prob=c(.5), type=1)) #ntaxa 4 -> rank 32
as.numeric(quantile(b_counts$n, prob=c(.75), type=1)) #ntaxa 20 -> rank 18
as.numeric(quantile(p_counts$n, prob=c(.75), type=1)) #ntaxa 15 -> rank 16

#Create a dataframe to specify species-rich species-poor thresholds
Site <- c("Malaysia", "Panama", "Malaysia", "Panama")
quartile <- c(34, 32, 18, 16)
richness_thresholds <- data.frame(quartile, Site)

##Plot##

#y-axis label
yl <- expression("|"*PD[SES]*"|")

NRI_plot <- ggplot(LRPD_NRI , aes(x=lineage_rank, y=Absolute_NRI)) + 
  geom_point(size = 1) +
  scale_x_continuous(expand = c(0,2)) +
  ylab("|NRI|") +
  xlab("Family Rank by Species Richness") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=11),
        strip.text = element_text(size=13)) + 
  facet_wrap(~ Site, scales = "free_x") +
  geom_line(data = values_plot_NRI, aes(x = lineage_rank, y = NRI), colour = "red") +
  geom_vline(data = breakpoints_NRI, aes(xintercept = Est.), linetype="dashed", color = "black", linewidth = 0.3) +
  geom_vline(data = richness_thresholds, aes(xintercept = quartile), linetype = "solid", color = "blue", linewidth = 0.3)

NTI_plot <- ggplot(LRPD_NTI , aes(x=lineage_rank, y=Absolute_NTI)) + 
  geom_point(size = 1) +
  scale_x_continuous(expand = c(0,2)) +
  ylab("|NTI|") +
  xlab("Family Rank by Species Richness") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=11),
        strip.text = element_text(size=13)) + 
  facet_wrap(~ Site, scales = "free_x") +
  theme(panel.spacing.x = unit(2, "lines")) + 
  geom_line(data = values_plot_NTI, aes(x = lineage_rank, y = NTI), colour = "red") +
  geom_vline(data = breakpoints_NTI, aes(xintercept = Est.), linetype="dashed", color = "black", linewidth = 0.3) +
  geom_vline(data = richness_thresholds, aes(xintercept = quartile), linetype = "solid", color = "blue", linewidth = 0.3)

PDses_plot <- ggplot(LRPD_PDses , aes(x=lineage_rank, y=Absolute_PDses)) + 
  geom_point(size = 1) +
  scale_x_continuous(expand = c(0,2)) +
  ylab(yl) +
  xlab("Family Rank by Species Richness") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=13),
        axis.text = element_text(size=11),
        strip.text = element_text(size=13)) + 
  facet_wrap(~ Site, scales = "free_x") +
  geom_line(data = values_plot_PDses, aes(x = lineage_rank, y = PDses), colour = "red") +
  geom_vline(data = breakpoints_PDses, aes(xintercept = Est.), linetype="dashed", color = "black", linewidth = 0.3) +
  geom_vline(data = richness_thresholds, aes(xintercept = quartile), linetype = "solid", color = "blue", linewidth = 0.3)


Lineage_rank_plot <- plot_grid(NRI_plot, NTI_plot, PDses_plot, labels = c("(a)", "(b)", "(c)"), label_size = 12, hjust=-0.3, vjust=2, align = "v", nrow=3)
ggsave("LRPD_plot.svg", Lineage_rank_plot, height=12, width=10, limitsize = F)


dev.off()

##Breakpoint Analysis##

#Chow test to check significance of breakpoints - do slopes differ around each breakpoint?
sctest(LRPD_borneo_PDses_mean$Absolute_PDses ~ LRPD_borneo_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_borneo_PDses[1,2])
sctest(LRPD_borneo_PDses_mean$Absolute_PDses ~ LRPD_borneo_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_borneo_PDses[2,2])
sctest(LRPD_borneo_PDses_mean$Absolute_PDses ~ LRPD_borneo_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_borneo_PDses[3,2])

sctest(LRPD_panama_PDses_mean$Absolute_PDses ~ LRPD_panama_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_panama_PDses[1,2])
sctest(LRPD_panama_PDses_mean$Absolute_PDses ~ LRPD_panama_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_panama_PDses[2,2])
sctest(LRPD_panama_PDses_mean$Absolute_PDses ~ LRPD_panama_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_panama_PDses[3,2])
sctest(LRPD_panama_PDses_mean$Absolute_PDses ~ LRPD_panama_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_panama_PDses[4,2])
sctest(LRPD_panama_PDses_mean$Absolute_PDses ~ LRPD_panama_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_panama_PDses[5,2])
sctest(LRPD_panama_PDses_mean$Absolute_PDses ~ LRPD_panama_PDses_mean$lineage_rank, type = "Chow", point = breakpoints_panama_PDses[6,2])

sctest(LRPD_borneo_NRI_mean$Absolute_NRI ~ LRPD_borneo_NRI_mean$lineage_rank, type = "Chow", point = breakpoints_borneo_NRI[1,2])

sctest(LRPD_panama_NRI_mean$Absolute_NRI ~ LRPD_panama_NRI_mean$lineage_rank, type = "Chow", point = breakpoints_panama_NRI[1,2])
sctest(LRPD_panama_NRI_mean$Absolute_NRI ~ LRPD_panama_NRI_mean$lineage_rank, type = "Chow", point = breakpoints_panama_NRI[2,2])

sctest(LRPD_borneo_NTI_mean$Absolute_NTI ~ LRPD_borneo_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_borneo_NTI[1,2])
sctest(LRPD_borneo_NTI_mean$Absolute_NTI ~ LRPD_borneo_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_borneo_NTI[2,2])
sctest(LRPD_borneo_NTI_mean$Absolute_NTI ~ LRPD_borneo_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_borneo_NTI[3,2])
sctest(LRPD_borneo_NTI_mean$Absolute_NTI ~ LRPD_borneo_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_borneo_NTI[4,2])

sctest(LRPD_panama_NTI_mean$Absolute_NTI ~ LRPD_panama_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_panama_NTI[1,2])
sctest(LRPD_panama_NTI_mean$Absolute_NTI ~ LRPD_panama_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_panama_NTI[2,2])
sctest(LRPD_panama_NTI_mean$Absolute_NTI ~ LRPD_panama_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_panama_NTI[3,2])
sctest(LRPD_panama_NTI_mean$Absolute_NTI ~ LRPD_panama_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_panama_NTI[4,2])
sctest(LRPD_panama_NTI_mean$Absolute_NTI ~ LRPD_panama_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_panama_NTI[5,2])
sctest(LRPD_panama_NTI_mean$Absolute_NTI ~ LRPD_panama_NTI_mean$lineage_rank, type = "Chow", point = breakpoints_panama_NTI[6,2])

