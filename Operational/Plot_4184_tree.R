######################
## 4K TREE PLOTTING ##
######################


##Set up##
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)
library(ggtree)
library(ggnewscale)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_4184 <- read.tree("Final_4k_Pars9NTBIN_ultraLSD2_root.tree")
clade_assignments <- read.csv("Taxonomy_Labels.csv")
family_tips <- read.csv("Family_Labels.csv")

##set up metadata
#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))

#identify outgroup sequences
outgroup_tips <- (phy_metadata%>%filter(!order=="Coleoptera"))$db_id

#Add column to full metadata to show if a sequence is from panama/borneo site (or non-site100)
phy_metadata$Site100 <- ifelse(phy_metadata$db_id%in%site_borneo_metadata$db_id, "Malaysia", ifelse(phy_metadata$db_id%in%site_panama_metadata$db_id, "Panama", ifelse(phy_metadata$db_id%in%outgroup_tips, "Outgroups", "NA")))
#phy_metadata$Site100 <- ifelse(phy_metadata$db_id%in%site_borneo_metadata$db_id, "Malaysia", ifelse(phy_metadata$db_id%in%site_panama_metadata$db_id, "Panama", ifelse(Full_Heatmap_metadata$db_id%in%outgroups, NA, "NA"))

#Remove rows to show site membership for heatmap layer 1
site_heatmap_data <- phy_metadata%>%select(db_id, Site100)%>%column_to_rownames(var = "db_id")

#Add column to full metadata to show family membership for Coleoptera sequences only
phy_metadata$family_group <- ifelse(phy_metadata$order=="Coleoptera", phy_metadata$family, "Outgroups")
#Remove rows to show family membership for heatmap layer 2
family_heatmap_data <- phy_metadata%>%column_to_rownames(var = "db_id")%>%select(family_group)%>%rename(Family = family_group)


##Load required function
#colour_match assigns colours to branches according to clade membership 
colour_match <- function(dataset, phylogeny, tip1, tip2, colour){
  node <- getMRCA(phylogeny, c(tip1, tip2))
  offspr <- tidytree::offspring(as_tibble(phylogeny), node)
  offspr <- dataset%>%filter(node%in%offspr$node)%>%filter(is.na(colour))
  dataset$colour <- ifelse(dataset$node%in%offspr$node, colour, dataset$colour)
  return(dataset)
}

##Plot##

#Define heatmap colours
site_colours <- c("Malaysia" = "#BB3AD8", "Panama" = "#57D83A", "NA" = "grey90")
family_colours <- c("Artematopodidae" = "grey80", "Boridae" = "grey80", "Bostrichidae" = "grey80", "Brentidae" = "grey80", "Buprestidae" = "grey80", "Byturidae" = "grey80", "Chrysomelidae" = "grey80", "Cleridae" = "grey80", "Coccinellidae" = "grey80", "Elateridae" = "grey80", "Elmidae" = "grey80", "Gyrinidae" = "grey80", "Heteroceridae" = "grey80", "Histeridae" = "grey80", "Hydraenidae" = "grey80", "Hydroscaphidae" = "grey80", "Hygrobiidae" = "grey80", "Iberobaeniidae" = "grey80", "Ischaliidae" = "grey80", "Latridiidae" = "grey80", "Leiodidae" = "grey80", "Limnichidae" = "grey80", "Megalopodidae" = "grey80", "Meloidae" = "grey80", "Mordellidae" = "grey80", "Nitidulidae" = "grey80", "Noteridae" = "grey80", "Oxypeltidae" = "grey80", "Phloeostichidae" = "grey80", "Platypodidae" = "grey80", "Prionoceridae" = "grey80", "Rhagophthalmidae" = "grey80", "Rhinorhipidae" = "grey80", "Salpingidae" = "grey80", "Scirtidae" = "grey80", "Silvanidae" = "grey80", "Tenebrionidae" = "grey80", "Vesperidae" = "grey80", "Aderidae" = "grey60", "Alexiidae" = "grey60", "Amphizoidae" = "grey60", "Anobiidae" = "grey60", "Cantharidae" = "grey60", "Carabidae" = "grey60", "Chaetosomatidae" = "grey60", "Chelonariidae" = "grey60", "Corylophidae" = "grey60", "Cryptophagidae" = "grey60", "Curculionidae" = "grey60", "Cybocephalidae" = "grey60", "Discolomatidae" = "grey60", "Eucinetidae" = "grey60", "Georissidae" = "grey60","Geotrupidae" = "grey60", "Glaresidae" = "grey60", "Haliplidae" = "grey60", "Kateretidae" = "grey60", "Lepiceridae" = "grey60", "Lymexylidae" = "grey60", "Melandryidae" = "grey60", "Ommatidae" = "grey60", "Orsodacnidae" = "grey60", "Passandridae" = "grey60", "Phalacridae" = "grey60", "Protocucujidae" = "grey60", "Psephenidae" = "grey60", "Rhipiceridae" = "grey60", "Scarabaeidae" = "grey60", "Silphidae" = "grey60", "Spercheidae" = "grey60", "Sphindidae" = "grey60", "Throscidae" = "grey60", "Trictenotomidae" = "grey60", "Trogossitidae" = "grey60", "Zopheridae" = "grey60", "Anthicidae" = "grey45", "Anthribidae" = "grey45", "Archeocrypticidae" = "grey45", "Belidae" = "grey45", "Biphyllidae" = "grey45", "Bolboceratidae" = "grey45", "Byrrhidae" = "grey45", "Cerambycidae" = "grey45", "Cerylonidae" = "grey45", "Clambidae" = "grey45", "Cupedidae" = "grey45", "Dascillidae" = "grey45", "Dermestidae" = "grey45", "Dytiscidae" = "grey45", "Erirhinidae" = "grey45", "Erotylidae" = "grey45", "Eucnemidae" = "grey45", "Eulichadidae" = "grey45", "Glaphyridae" = "grey45", "Hydrophilidae" = "grey45", "Lucanidae" = "grey45", "Meruidae" = "grey45", "Oedemeridae" = "grey45", "Omalisidae" = "grey45", "Omethidae" = "grey45", "Phengodidae" = "grey45", "Phloiophilidae" = "grey45", "Phycosecidae" = "grey45", "Propalticidae" = "grey45", "Prostomidae" = "grey45", "Ptiliidae" = "grey45", "Ripiphoridae" = "grey45", "Sphaeriusidae" = "grey45", "Tetratomidae" = "grey45", "Trachypachidae" = "grey45", "Aspidytidae" = "grey30", "Attelabidae" = "grey30", "Bothrideridae" = "grey30", "Brachyceridae" = "grey30", "Callirhipidae" = "grey30", "Cerophytidae" = "grey30", "Cicindelidae" = "grey30", "Ciidae" = "grey30", "Cucujidae" = "grey30", "Derodontidae" = "grey30", "Disteniidae" = "grey30", "Dryopidae" = "grey30", "Endomychidae" = "grey30", "Helophoridae" = "grey30", "Helotidae" = "grey30", "Hybosoridae" = "grey30", "Hydrochidae" = "grey30", "Laemophloeidae" = "grey30", "Lampyridae" = "grey30", "Lycidae" = "grey30", "Melyridae" = "grey30", "Micromalthidae" = "grey30", "Monotomidae" = "grey30", "Mycetophagidae" = "grey30", "Mycteridae" = "grey30", "Nemonychidae" = "grey30", "Nosodendridae" = "grey30", "Passalidae" = "grey30", "Ptilodactylidae" = "grey30", "Ptinidae" = "grey30", "Pyrochroidae" = "grey30", "Scraptiidae" = "grey30", "Staphylinidae" = "grey30", "Stenotrachelidae" = "grey30", "Torridincolidae" = "grey30", "Trogidae" = "grey30")


#plot full tree, adding heatmap to show which sequences belong to site-samples
plot_phy_4184 <- ggtree(phy_4184, branch.length='none') + geom_tiplab() + coord_cartesian(clip = "off") + vexpand(.0005, 1)


##Colouring branches
#Get data of plotted tree
phy_data <- plot_phy_4184$data
#Get taxonomy assignments for all tips
tip_taxonomy <- phy_metadata%>%filter(db_id%in%phy_4184$tip.label)%>%select(db_id, superfamily, suborder)%>%rename(label = db_id)
#Add taxonomy to plotted tree data for tips
plot_terminalnode_data <- merge(phy_data, tip_taxonomy, by="label")
#Assign colours to branches according to terminal branches according to clade membership
plot_terminalnode_data$colour <- ifelse(plot_terminalnode_data$superfamily == "Caraboidea", "#0081FA", ifelse(plot_terminalnode_data$superfamily == "Dytiscoidea", "#06ACE9", ifelse(plot_terminalnode_data$superfamily == "Gyrinoidea", "#382CF3", ifelse(plot_terminalnode_data$superfamily == "Haliploidea", "#54D3FF", ifelse(plot_terminalnode_data$superfamily == "Bostrichoidea", "#47F947", ifelse(plot_terminalnode_data$superfamily == "Buprestoidea", "#FF3352", ifelse(plot_terminalnode_data$superfamily == "Byrrhoidea", "#F88379", ifelse(plot_terminalnode_data$superfamily == "Chrysomeloidea", "#EB4AF0", ifelse(plot_terminalnode_data$superfamily == "Cleroidea", "#EB8B73", ifelse(plot_terminalnode_data$superfamily == "Coccinelloidea", "#25D8C5", ifelse(plot_terminalnode_data$superfamily == "Cucujoidea", "#19B519", ifelse(plot_terminalnode_data$superfamily == "Curculionoidea", "#DEB46B", ifelse(plot_terminalnode_data$superfamily == "Dascilloidea", "#A6E52E", ifelse(plot_terminalnode_data$superfamily == "Derodontoidea", "#E03EF5", ifelse(plot_terminalnode_data$superfamily == "Dryopoidea", "#00E6D3", ifelse(plot_terminalnode_data$superfamily == "Elateroidea", "#46E160", ifelse(plot_terminalnode_data$superfamily == "Histeroidea", "#FEDC00", ifelse(plot_terminalnode_data$superfamily == "Hydrophiloidea", "#C75DE6", ifelse(plot_terminalnode_data$superfamily == "Lymexyloidea", "#FF359A", ifelse(plot_terminalnode_data$superfamily == "Rhinorhipoidea", "#FFF747", ifelse(plot_terminalnode_data$superfamily == "Scarabaeoidea", "#DC3535", ifelse(plot_terminalnode_data$superfamily == "Scirtoidea", "#FF5DCC", ifelse(plot_terminalnode_data$superfamily == "Staphylinoidea", "#FF89A7", ifelse(plot_terminalnode_data$superfamily == "Tenebrionoidea", "#FFBF00", NA))))))))))))))))))))))))
plot_terminalnode_data$colour <- ifelse(plot_terminalnode_data$suborder == "Myxophaga", "#EC8510", ifelse(plot_terminalnode_data$suborder == "Archostemata", "#A60CE5", plot_terminalnode_data$colour))
plot_terminalnode_data$colour <- ifelse(plot_terminalnode_data$label%in%outgroup_tips, "grey40", plot_terminalnode_data$colour)
plot_terminalnode_data$colour[is.na(plot_terminalnode_data$colour)] <- "black"

#Filter plotted tree data for internal nodes
plot_internalnode_data <- phy_data%>%filter(isTip==FALSE)
plot_internalnode_data$superfamily <- NA
plot_internalnode_data$suborder <- NA
plot_internalnode_data$colour <- NA

#Remove rows from node assignment data to exclude taxonomy labels/include only colour assignments
colour_assignments <- clade_assignments%>%filter(ntaxa!=1)%>%select(taxa_1, taxa_2, colour)%>%filter(!is.na(taxa_1))%>%filter(!is.na(taxa_2))%>%filter(!is.na(colour))
#Loop through each row to assign colours to descendants according to clade membership
for (i in 1:nrow(colour_assignments)) {
  plot_internalnode_data <- colour_match(plot_internalnode_data, phy_4184, colour_assignments$taxa_1[i], colour_assignments$taxa_2[i], colour_assignments$colour[i])
}

#Combine plotted tree data for internal and terminal nodes
plot_data <- rbind(plot_terminalnode_data, plot_internalnode_data)

#Manually assign colours to remaining internal branches
plot_data[4186,12] <- "black"
plot_data[4187,12] <- "black"
plot_data[4188,12] <- "black"
plot_data[4189,12] <- "black"
plot_data[4190,12] <- "black"
plot_data[4191,12] <- "black"
plot_data[4192,12] <- "black"
plot_data[4193,12] <- "black"
plot_data[4194,12] <- "black"
plot_data[4195,12] <- "#46E160"
plot_data[4446,12] <- "black"
plot_data[4447,12] <- "#FF3352"
plot_data[4480,12] <- "#F88379"
plot_data[4481,12] <- "#00E6D3"
plot_data[4553,12] <- "#A6E52E"
plot_data[4555,12] <- "black"
plot_data[4556,12] <- "black"
plot_data[4557,12] <- "black"
plot_data[4558,12] <- "#FEDC00"
plot_data[4636,12] <- "#C75DE6"
plot_data[4689,12] <- "black"
plot_data[4690,12] <- "#DC3535"
plot_data[4891,12] <- "#FF89A7"
plot_data[5665,12] <- "black"
plot_data[5666,12] <- "black"
plot_data[5667,12] <- "black"
plot_data[5668,12] <- "black"
plot_data[5669,12] <- "black"
plot_data[5670,12] <- "black"
plot_data[5671,12] <- "black"
plot_data[5672,12] <- "#DEB46B"
plot_data[6427,12] <- "#EB4AF0"
plot_data[7099,12] <- "#19B519"
plot_data[7380,12] <- "#19B519"
plot_data[7384,12] <- "black"
plot_data[7385,12] <- "#FFBF00"
plot_data[7759,12] <- "#FF359A"
plot_data[7763,12] <- "#25D8C5"
plot_data[7961,12] <- "black"
plot_data[7962,12] <- "#EB8B73"
plot_data[8071,12] <- "#47F947"
plot_data[8140,12] <- "#E03EF5"
plot_data[8142,12] <- "black"
plot_data[8143,12] <- "#FF5DCC"
plot_data[8149,12] <- "#FF5DCC"
plot_data[8165,12] <- "black"
plot_data[8166,12] <- "black"
plot_data[8167,12] <- "black"
plot_data[8168,12] <- "#0081FA"
plot_data[8288,12] <- "black"
plot_data[8289,12] <- "#06ACE9"
plot_data[8339,12] <- "#54D3FF"
plot_data[8340,12] <- "#382CF3"
plot_data[8345,12] <- "black"
plot_data[8346,12] <- "#A60CE5"
plot_data[8350,12] <- "#EC8510"
plot_data[8355,12] <- "grey40"


##Plot full tree with coloured branches
plot <- plot_phy_4184 %<+% plot_data + aes(color=I(colour))

##Add heatmap layer 1 showing site membership
plot_heatmap_site <- gheatmap(plot, site_heatmap_data, offset=2.9, width = 0.0125, color = NULL, colnames = T, colnames_position ='top', custom_column_labels = "Site", font.size = 4, hjust = 0, colnames_angle = 35, colnames_offset_x = -0.4) +   
  scale_fill_manual(values = site_colours, breaks = c("Malaysia", "Panama", "NA"), na.value = NA, name="Site") +
  guides(fill=guide_legend(title.position="top", ncol=1)) + 
  theme(legend.position="right", 
        legend.justification = "top", 
        legend.box.margin = margin(45, 0, 0, -35))

##Add heatmap layer 2 showing family membership
plot_heatmap_site <- plot_heatmap_site + new_scale_fill()
plot_heatmap_site_family <-  gheatmap(plot_heatmap_site, family_heatmap_data, offset=3.8, width = 0.0125, color = NULL, colnames = T, colnames_position ='top', font.size = 4, hjust = 0, colnames_angle = 35, colnames_offset_x = -0.45) +   
  scale_fill_manual(values = family_colours, na.value = NA) + 
  guides(fill="none")


##Add clade labels to ancestral node

#identify node of the only monotypic superfamily
rhino_node <- (plot_phy_4184$data%>%filter(label=="SRAA00082"))$node

#Remove rows from node assignment data to exclude colour assignments/include only taxonomic labels
labels <- clade_assignments%>%filter(ntaxa!=1)%>%select(taxonomic_rank, group, taxa_1, taxa_2)%>%filter(taxonomic_rank%in%c("Outgroups", "Suborder", "Infraorder", "Superfamily"))
labels <- labels[!is.na(labels$taxa_2),]
#find ancestor node of clade
labels <- labels%>%filter(!is.na(taxa_2))%>%rowwise()%>%mutate(node = getMRCA(phy_4184, c(taxa_1, taxa_2)))
#unite labels if higher-level classification of clade is equal to lower-level classification
labels <- labels%>%pivot_wider(names_from = "taxonomic_rank", values_from = "group")%>% unite("Label", Suborder:Outgroups, na.rm = TRUE, sep = ": ")
labels <- as.data.frame(labels)

plot_heatmap_site_family_label <- plot_heatmap_site_family + geom_text2(aes(subset=(node == labels[1,3])), cex=5, label=labels[1,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[2,3])), cex=5, label=labels[2,4], color ="black", nudge_x = -0.2, nudge_y = -0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[3,3])), cex=5, label=labels[3,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[4,3])), cex=5, label=labels[4,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) + 
  geom_text2(aes(subset=(node == labels[5,3])), cex=4, label=labels[5,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[6,3])), cex=4, label=labels[6,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[7,3])), cex=4, label=labels[7,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[8,3])), cex=4, label=labels[8,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[9,3])), cex=4, label=labels[9,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[10,3])), cex=4, label=labels[10,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[11,3])), cex=4, label=labels[11,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[12,3])), cex=4, label=labels[12,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[13,3])), cex=4, label=labels[13,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[14,3])), cex=4, label=labels[14,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[15,3])), cex=4, label=labels[15,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[16,3])), cex=4, label=labels[16,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[17,3])), cex=4, label=labels[17,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[18,3])), cex=4, label=labels[18,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[19,3])), cex=4, label=labels[19,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[20,3])), cex=4, label=labels[20,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[21,3])), cex=4, label=labels[21,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[22,3])), cex=4, label=labels[22,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[23,3])), cex=4, label=labels[23,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[24,3])), cex=4, label=labels[24,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[25,3])), cex=4, label=labels[25,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[26,3])), cex=4, label=labels[26,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[27,3])), cex=4, label=labels[27,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[28,3])), cex=4, label=labels[28,4], color ="black", nudge_x = -0.2, nudge_y = -0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[29,3])), cex=4, label=labels[29,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[30,3])), cex=4, label=labels[30,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[31,3])), cex=4, label=labels[31,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == labels[32,3])), cex=4, label=labels[32,4], color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1) +
  geom_text2(aes(subset=(node == rhino_node)), cex=4, label="Rhinorhipoidea", color ="black", nudge_x = -0.2, nudge_y = 0.5, hjust = 1)

##Add family labels to correspond to heatmap layer 2

family_labels <- plot_data%>%filter(label%in%family_tips$label)
family_labels <- merge(family_labels, family_tips, by="label")
family_labels <- as.data.frame(family_labels%>%select(node, family, label))

plot_heatmap_site_family_label_2 <- plot_heatmap_site_family_label + geom_text2(aes(subset=(node == family_labels[1,1])), cex=3, label=family_labels[1,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[2,1])), cex=3, label=family_labels[2,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) + 
  geom_text2(aes(subset=(node == family_labels[3,1])), cex=3, label=family_labels[3,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) + 
  geom_text2(aes(subset=(node == family_labels[4,1])), cex=3, label=family_labels[4,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) + 
  geom_text2(aes(subset=(node == family_labels[5,1])), cex=3, label=family_labels[5,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) + 
  geom_text2(aes(subset=(node == family_labels[6,1])), cex=3, label=family_labels[6,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) + 
  geom_text2(aes(subset=(node == family_labels[7,1])), cex=3, label=family_labels[7,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) + 
  geom_text2(aes(subset=(node == family_labels[8,1])), cex=3, label=family_labels[8,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[9,1])), cex=3, label=family_labels[9,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[10,1])), cex=3, label=family_labels[10,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[11,1])), cex=3, label=family_labels[11,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[12,1])), cex=3, label=family_labels[12,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[13,1])), cex=3, label=family_labels[13,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[14,1])), cex=3, label=family_labels[14,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[15,1])), cex=3, label=family_labels[15,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[16,1])), cex=3, label=family_labels[16,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[17,1])), cex=3, label=family_labels[17,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[18,1])), cex=3, label=family_labels[18,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[19,1])), cex=3, label=family_labels[19,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[20,1])), cex=3, label=family_labels[20,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[21,1])), cex=3, label=family_labels[21,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[22,1])), cex=3, label=family_labels[22,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[23,1])), cex=3, label=family_labels[23,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[24,1])), cex=3, label=family_labels[24,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[25,1])), cex=3, label=family_labels[25,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[26,1])), cex=3, label=family_labels[26,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[27,1])), cex=3, label=family_labels[27,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[28,1])), cex=3, label=family_labels[28,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[29,1])), cex=3, label=family_labels[29,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[30,1])), cex=3, label=family_labels[30,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[31,1])), cex=3, label=family_labels[31,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[32,1])), cex=3, label=family_labels[32,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[33,1])), cex=3, label=family_labels[33,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[34,1])), cex=3, label=family_labels[34,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[35,1])), cex=3, label=family_labels[35,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[36,1])), cex=3, label=family_labels[36,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[37,1])), cex=3, label=family_labels[37,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[38,1])), cex=3, label=family_labels[38,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[39,1])), cex=3, label=family_labels[39,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[40,1])), cex=3, label=family_labels[40,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[41,1])), cex=3, label=family_labels[41,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[42,1])), cex=3, label=family_labels[42,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[43,1])), cex=3, label=family_labels[43,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[44,1])), cex=3, label=family_labels[44,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[45,1])), cex=3, label=family_labels[45,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[46,1])), cex=3, label=family_labels[46,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[47,1])), cex=3, label=family_labels[47,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[48,1])), cex=3, label=family_labels[48,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[49,1])), cex=3, label=family_labels[49,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[50,1])), cex=3, label=family_labels[50,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[51,1])), cex=3, label=family_labels[51,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[52,1])), cex=3, label=family_labels[52,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[53,1])), cex=3, label=family_labels[53,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[54,1])), cex=3, label=family_labels[54,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[55,1])), cex=3, label=family_labels[55,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[56,1])), cex=3, label=family_labels[56,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[57,1])), cex=3, label=family_labels[57,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[58,1])), cex=3, label=family_labels[58,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[59,1])), cex=3, label=family_labels[59,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[60,1])), cex=3, label=family_labels[60,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[61,1])), cex=3, label=family_labels[61,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[62,1])), cex=3, label=family_labels[62,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[63,1])), cex=3, label=family_labels[63,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[64,1])), cex=3, label=family_labels[64,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[65,1])), cex=3, label=family_labels[65,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[66,1])), cex=3, label=family_labels[66,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[67,1])), cex=3, label=family_labels[67,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[68,1])), cex=3, label=family_labels[68,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[69,1])), cex=3, label=family_labels[69,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[70,1])), cex=3, label=family_labels[70,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[71,1])), cex=3, label=family_labels[71,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[72,1])), cex=3, label=family_labels[72,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[73,1])), cex=3, label=family_labels[73,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[74,1])), cex=3, label=family_labels[74,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[75,1])), cex=3, label=family_labels[75,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[76,1])), cex=3, label=family_labels[76,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[77,1])), cex=3, label=family_labels[77,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[78,1])), cex=3, label=family_labels[78,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[79,1])), cex=3, label=family_labels[79,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[80,1])), cex=3, label=family_labels[80,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[81,1])), cex=3, label=family_labels[81,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[82,1])), cex=3, label=family_labels[82,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[83,1])), cex=3, label=family_labels[83,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[84,1])), cex=3, label=family_labels[84,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[85,1])), cex=3, label=family_labels[85,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[86,1])), cex=3, label=family_labels[86,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[87,1])), cex=3, label=family_labels[87,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[88,1])), cex=3, label=family_labels[88,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[89,1])), cex=3, label=family_labels[89,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[90,1])), cex=3, label=family_labels[90,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[91,1])), cex=3, label=family_labels[91,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[92,1])), cex=3, label=family_labels[92,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[93,1])), cex=3, label=family_labels[93,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[94,1])), cex=3, label=family_labels[94,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[95,1])), cex=3, label=family_labels[95,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[96,1])), cex=3, label=family_labels[96,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[97,1])), cex=3, label=family_labels[97,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[98,1])), cex=3, label=family_labels[98,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[99,1])), cex=3, label=family_labels[99,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[100,1])), cex=3, label=family_labels[100,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[101,1])), cex=3, label=family_labels[101,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[102,1])), cex=3, label=family_labels[102,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[103,1])), cex=3, label=family_labels[103,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[104,1])), cex=3, label=family_labels[104,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[105,1])), cex=3, label=family_labels[105,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[106,1])), cex=3, label=family_labels[106,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[107,1])), cex=3, label=family_labels[107,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[108,1])), cex=3, label=family_labels[108,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[109,1])), cex=3, label=family_labels[109,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[110,1])), cex=3, label=family_labels[110,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[111,1])), cex=3, label=family_labels[111,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[112,1])), cex=3, label=family_labels[112,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[113,1])), cex=3, label=family_labels[113,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[114,1])), cex=3, label=family_labels[114,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[115,1])), cex=3, label=family_labels[115,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[116,1])), cex=3, label=family_labels[116,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[117,1])), cex=3, label=family_labels[117,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[118,1])), cex=3, label=family_labels[118,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[119,1])), cex=3, label=family_labels[119,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[120,1])), cex=3, label=family_labels[120,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[121,1])), cex=3, label=family_labels[121,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[122,1])), cex=3, label=family_labels[122,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[123,1])), cex=3, label=family_labels[123,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[124,1])), cex=3, label=family_labels[124,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[125,1])), cex=3, label=family_labels[125,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[126,1])), cex=3, label=family_labels[126,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[127,1])), cex=3, label=family_labels[127,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[128,1])), cex=3, label=family_labels[128,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[129,1])), cex=3, label=family_labels[129,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[130,1])), cex=3, label=family_labels[130,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[131,1])), cex=3, label=family_labels[131,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[132,1])), cex=3, label=family_labels[132,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[133,1])), cex=3, label=family_labels[133,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[134,1])), cex=3, label=family_labels[134,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[135,1])), cex=3, label=family_labels[135,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[136,1])), cex=3, label=family_labels[136,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[137,1])), cex=3, label=family_labels[137,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[138,1])), cex=3, label=family_labels[138,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[139,1])), cex=3, label=family_labels[139,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[140,1])), cex=3, label=family_labels[140,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[141,1])), cex=3, label=family_labels[141,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[142,1])), cex=3, label=family_labels[142,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[143,1])), cex=3, label=family_labels[143,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[144,1])), cex=3, label=family_labels[144,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[145,1])), cex=3, label=family_labels[145,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[146,1])), cex=3, label=family_labels[146,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[147,1])), cex=3, label=family_labels[147,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[148,1])), cex=3, label=family_labels[148,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[149,1])), cex=3, label=family_labels[149,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[150,1])), cex=3, label=family_labels[150,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[151,1])), cex=3, label=family_labels[151,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[152,1])), cex=3, label=family_labels[152,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[153,1])), cex=3, label=family_labels[153,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[154,1])), cex=3, label=family_labels[154,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[155,1])), cex=3, label=family_labels[155,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[156,1])), cex=3, label=family_labels[156,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[157,1])), cex=3, label=family_labels[157,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[158,1])), cex=3, label=family_labels[158,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[159,1])), cex=3, label=family_labels[159,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[160,1])), cex=3, label=family_labels[160,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[161,1])), cex=3, label=family_labels[161,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[162,1])), cex=3, label=family_labels[162,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[163,1])), cex=3, label=family_labels[163,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[164,1])), cex=3, label=family_labels[164,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[165,1])), cex=3, label=family_labels[165,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[166,1])), cex=3, label=family_labels[166,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[167,1])), cex=3, label=family_labels[167,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[168,1])), cex=3, label=family_labels[168,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[169,1])), cex=3, label=family_labels[169,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[170,1])), cex=3, label=family_labels[170,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[171,1])), cex=3, label=family_labels[171,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[172,1])), cex=3, label=family_labels[172,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[173,1])), cex=3, label=family_labels[173,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[174,1])), cex=3, label=family_labels[174,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[175,1])), cex=3, label=family_labels[175,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[176,1])), cex=3, label=family_labels[176,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[177,1])), cex=3, label=family_labels[177,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[178,1])), cex=3, label=family_labels[178,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[179,1])), cex=3, label=family_labels[179,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[180,1])), cex=3, label=family_labels[180,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[181,1])), cex=3, label=family_labels[181,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[182,1])), cex=3, label=family_labels[182,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[183,1])), cex=3, label=family_labels[183,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[184,1])), cex=3, label=family_labels[184,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[185,1])), cex=3, label=family_labels[185,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[186,1])), cex=3, label=family_labels[186,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[187,1])), cex=3, label=family_labels[187,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[188,1])), cex=3, label=family_labels[188,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[189,1])), cex=3, label=family_labels[189,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[190,1])), cex=3, label=family_labels[190,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[191,1])), cex=3, label=family_labels[191,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[192,1])), cex=3, label=family_labels[192,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[193,1])), cex=3, label=family_labels[193,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[194,1])), cex=3, label=family_labels[194,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[195,1])), cex=3, label=family_labels[195,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[196,1])), cex=3, label=family_labels[196,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[197,1])), cex=3, label=family_labels[197,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[198,1])), cex=3, label=family_labels[198,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[199,1])), cex=3, label=family_labels[199,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[200,1])), cex=3, label=family_labels[200,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0) +
  geom_text2(aes(subset=(node == family_labels[201,1])), cex=3, label=family_labels[201,2], color ="black", nudge_x = 5.025, nudge_y = 0, hjust = 0)

##Add label for outgroup taxa alongside heatmap
plot_heatmap_site_family_label_3 <- plot_heatmap_site_family_label_2 + geom_cladelabel(node=labels[33,3], label=labels[33,4], color="grey40", offset=2.115, align=TRUE, fontsize = 5)

ggsave("Heatmap_Tree_plot.pdf", plot = plot_heatmap_site_family_label_3, width=25, height = 750, limitsize = FALSE)

dev.off()