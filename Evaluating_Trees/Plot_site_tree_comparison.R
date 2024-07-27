##################################################
## PRUNED VS RECONSTRUCTED TANGLEGRAM PLOTTING ##
##################################################


##Set up##
rm(list = ls())

#Load libraries
library(tidyverse)
library(ape)
library(ggtree)
library(cowplot)

#Read in data
phy_metadata <- read.csv("metadata_4184.csv")
phy_borneo_prune <- read.tree("Borneo_Prune_Pars9_ultraLSD2.tree")
phy_borneo_recon <- read.tree("Borneo_Recon_Pars9_ultraLSD2.tree")
phy_panama_prune <- read.tree("Panama_Prune_Pars9_ultraLSD2.tree")
phy_panama_recon <- read.tree("Panama_Recon_Pars9_ultraLSD2.tree")
borneo_clade_assignments <- read.csv("Taxonomy_Labels_Borneo.csv")
panama_clade_assignments <- read.csv("Taxonomy_Labels_Panama.csv")


#Filter metadata to sites
coleoptera_metadata <- phy_metadata%>%filter(order=="Coleoptera")
site_borneo_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Malaysia")%>%filter(locality%in%c("Poring", "Danum Valley"))
site_panama_metadata <- coleoptera_metadata%>%filter(grepl("BIOD", db_id))%>%filter(country=="Panama")%>%filter(locality%in%c("Cerro Hoya", "Santa Fe"))


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

##Borneo

#Plot pruned tree
plot_b_p <- ggtree(phy_borneo_prune, branch.length = 'none', linewidth = 0.4) + 
  theme(plot.margin = unit(c(0.5,0,0.5,0),"cm"))

#Get data of plotted tree
data_borneo_prune <- plot_b_p$data
#Get taxonomy assignments for all tips
tip_taxonomy_borneo_prune <- site_borneo_metadata%>%filter(db_id%in%phy_borneo_prune$tip.label)%>%select(db_id, superfamily, family)%>%rename(label = db_id)
#Add taxonomy to plotted tree data for tips
data_borneo_prune_terminalnode <- merge(data_borneo_prune, tip_taxonomy_borneo_prune, by="label")
#Assign colours to branches according to terminal branches according to clade membership
data_borneo_prune_terminalnode$colour <- ifelse(data_borneo_prune_terminalnode$superfamily == "Caraboidea", "#0081FA", ifelse(data_borneo_prune_terminalnode$superfamily == "Dytiscoidea", "#06ACE9", ifelse(data_borneo_prune_terminalnode$superfamily == "Gyrinoidea", "#382CF3", ifelse(data_borneo_prune_terminalnode$superfamily == "Haliploidea", "#54D3FF", ifelse(data_borneo_prune_terminalnode$superfamily == "Bostrichoidea", "#47F947", ifelse(data_borneo_prune_terminalnode$superfamily == "Buprestoidea", "#FF3352", ifelse(data_borneo_prune_terminalnode$superfamily == "Byrrhoidea", "#F88379", ifelse(data_borneo_prune_terminalnode$superfamily == "Chrysomeloidea", "#EB4AF0", ifelse(data_borneo_prune_terminalnode$superfamily == "Cleroidea", "#EB8B73", ifelse(data_borneo_prune_terminalnode$superfamily == "Coccinelloidea", "#25D8C5", ifelse(data_borneo_prune_terminalnode$superfamily == "Cucujoidea", "#19B519", ifelse(data_borneo_prune_terminalnode$superfamily == "Curculionoidea", "#DEB46B", ifelse(data_borneo_prune_terminalnode$superfamily == "Dascilloidea", "#A6E52E", ifelse(data_borneo_prune_terminalnode$superfamily == "Derodontoidea", "#E03EF5", ifelse(data_borneo_prune_terminalnode$superfamily == "Dryopoidea", "#00E6D3", ifelse(data_borneo_prune_terminalnode$superfamily == "Elateroidea", "#46E160", ifelse(data_borneo_prune_terminalnode$superfamily == "Histeroidea", "#FEDC00", ifelse(data_borneo_prune_terminalnode$superfamily == "Hydrophiloidea", "#C75DE6", ifelse(data_borneo_prune_terminalnode$superfamily == "Lymexyloidea", "#FF359A", ifelse(data_borneo_prune_terminalnode$superfamily == "Rhinorhipoidea", "#FFF747", ifelse(data_borneo_prune_terminalnode$superfamily == "Scarabaeoidea", "#DC3535", ifelse(data_borneo_prune_terminalnode$superfamily == "Scirtoidea", "#FF5DCC", ifelse(data_borneo_prune_terminalnode$superfamily == "Staphylinoidea", "#FF89A7", ifelse(data_borneo_prune_terminalnode$superfamily == "Tenebrionoidea", "#FFBF00", NA))))))))))))))))))))))))
data_borneo_prune_terminalnode$colour[is.na(data_borneo_prune_terminalnode$colour)] <- "black"

#Filter plotted tree data for internal nodes
data_borneo_prune_internalnode <- data_borneo_prune%>%filter(isTip==FALSE)
data_borneo_prune_internalnode$superfamily <- NA
data_borneo_prune_internalnode$family <- NA
data_borneo_prune_internalnode$colour <- NA


#Remove rows from node assignment data to exclude taxonomy labels/include only colour assignments
borneo_prune_colour_assignments <- borneo_clade_assignments%>%filter(Tree=="Borneo_Prune")%>%filter(ntaxa!=1)%>%select(taxa_1, taxa_2, colour)%>%filter(!is.na(taxa_1))%>%filter(!is.na(taxa_2))%>%filter(!is.na(colour))
#Loop through each row to assign colours to descendants according to clade membership
for (i in 1:nrow(borneo_prune_colour_assignments)) {
  data_borneo_prune_internalnode <- colour_match(data_borneo_prune_internalnode, phy_borneo_prune, borneo_prune_colour_assignments$taxa_1[i], borneo_prune_colour_assignments$taxa_2[i], borneo_prune_colour_assignments$colour[i])
}

#Combine plotted tree data for internal and terminal nodes
plot_borneo_prune_data <- rbind(data_borneo_prune_terminalnode, data_borneo_prune_internalnode)

#Manually assign colours to remaining internal branches
plot_borneo_prune_data[1311,12] <- "black"
plot_borneo_prune_data[1312,12] <- "black"
plot_borneo_prune_data[1313,12] <- "black"
plot_borneo_prune_data[1314,12] <- "black"
plot_borneo_prune_data[1315,12] <- "black"
plot_borneo_prune_data[1316,12] <- "#46E160"
plot_borneo_prune_data[1388,12] <- "#FF3352"
plot_borneo_prune_data[1393,12] <- "#00E6D3"
plot_borneo_prune_data[1433,12] <- "black"
plot_borneo_prune_data[1434,12] <- "black"
plot_borneo_prune_data[1435,12] <- "black"
plot_borneo_prune_data[1436,12] <- "#FEDC00"
plot_borneo_prune_data[1473,12] <- "#C75DE6"
plot_borneo_prune_data[1487,12] <- "black"
plot_borneo_prune_data[1488,12] <- "#DC3535"
plot_borneo_prune_data[1536,12] <- "#FF89A7"
plot_borneo_prune_data[1939,12] <- "black"
plot_borneo_prune_data[1940,12] <- "black"
plot_borneo_prune_data[1941,12] <- "black"
plot_borneo_prune_data[1942,12] <- "black"
plot_borneo_prune_data[1943,12] <- "black"
plot_borneo_prune_data[1944,12] <- "black"
plot_borneo_prune_data[1945,12] <- "black"
plot_borneo_prune_data[1946,12] <- "#DEB46B"
plot_borneo_prune_data[2104,12] <- "#EB4AF0"
plot_borneo_prune_data[2252,12] <- "#19B519"
plot_borneo_prune_data[2353,12] <- "#19B519"
plot_borneo_prune_data[2355,12] <- "black"
plot_borneo_prune_data[2356,12] <- "#FFBF00"
plot_borneo_prune_data[2459,12] <- "#25D8C5"
plot_borneo_prune_data[2534,12] <- "black"
plot_borneo_prune_data[2535,12] <- "#EB8B73"
plot_borneo_prune_data[2579,12] <- "#47F947"
plot_borneo_prune_data[2586,12] <- "#FF5DCC"
plot_borneo_prune_data[2588,12] <- "#FF5DCC"
plot_borneo_prune_data[2590,12] <- "black"
plot_borneo_prune_data[2591,12] <- "#0081FA"
plot_borneo_prune_data[2616,12] <- "#06ACE9"

##Plot pruned tree with coloured branches
plot_borneo_prune <- plot_b_p %<+% plot_borneo_prune_data + aes(color=I(colour))


#Plot reconstructed tree mirrored in y-axis
#Plot pruned tree
plot_b_r <- ggtree(phy_borneo_recon, branch.length = 'none', linewidth = 0.4) + 
  scale_x_reverse() + 
  theme(plot.margin = unit(c(0.5,0,0.5,0),"cm"))

#Get data of plotted tree
data_borneo_recon <- plot_b_r$data
#Get taxonomy assignments for all tips
tip_taxonomy_borneo_recon <- site_borneo_metadata%>%filter(db_id%in%phy_borneo_recon$tip.label)%>%select(db_id, superfamily, family)%>%rename(label = db_id)
#Add taxonomy to plotted tree data for tips
data_borneo_recon_terminalnode <- merge(data_borneo_recon, tip_taxonomy_borneo_recon, by="label")
#Assign colours to branches according to terminal branches according to clade membership
data_borneo_recon_terminalnode$colour <- ifelse(data_borneo_recon_terminalnode$superfamily == "Caraboidea", "#0081FA", ifelse(data_borneo_recon_terminalnode$superfamily == "Dytiscoidea", "#06ACE9", ifelse(data_borneo_recon_terminalnode$superfamily == "Gyrinoidea", "#382CF3", ifelse(data_borneo_recon_terminalnode$superfamily == "Haliploidea", "#54D3FF", ifelse(data_borneo_recon_terminalnode$superfamily == "Bostrichoidea", "#47F947", ifelse(data_borneo_recon_terminalnode$superfamily == "Buprestoidea", "#FF3352", ifelse(data_borneo_recon_terminalnode$superfamily == "Byrrhoidea", "#F88379", ifelse(data_borneo_recon_terminalnode$superfamily == "Chrysomeloidea", "#EB4AF0", ifelse(data_borneo_recon_terminalnode$superfamily == "Cleroidea", "#EB8B73", ifelse(data_borneo_recon_terminalnode$superfamily == "Coccinelloidea", "#25D8C5", ifelse(data_borneo_recon_terminalnode$superfamily == "Cucujoidea", "#19B519", ifelse(data_borneo_recon_terminalnode$superfamily == "Curculionoidea", "#DEB46B", ifelse(data_borneo_recon_terminalnode$superfamily == "Dascilloidea", "#A6E52E", ifelse(data_borneo_recon_terminalnode$superfamily == "Derodontoidea", "#E03EF5", ifelse(data_borneo_recon_terminalnode$superfamily == "Dryopoidea", "#00E6D3", ifelse(data_borneo_recon_terminalnode$superfamily == "Elateroidea", "#46E160", ifelse(data_borneo_recon_terminalnode$superfamily == "Histeroidea", "#FEDC00", ifelse(data_borneo_recon_terminalnode$superfamily == "Hydrophiloidea", "#C75DE6", ifelse(data_borneo_recon_terminalnode$superfamily == "Lymexyloidea", "#FF359A", ifelse(data_borneo_recon_terminalnode$superfamily == "Rhinorhipoidea", "#FFF747", ifelse(data_borneo_recon_terminalnode$superfamily == "Scarabaeoidea", "#DC3535", ifelse(data_borneo_recon_terminalnode$superfamily == "Scirtoidea", "#FF5DCC", ifelse(data_borneo_recon_terminalnode$superfamily == "Staphylinoidea", "#FF89A7", ifelse(data_borneo_recon_terminalnode$superfamily == "Tenebrionoidea", "#FFBF00", NA))))))))))))))))))))))))
data_borneo_recon_terminalnode$colour[is.na(data_borneo_recon_terminalnode$colour)] <- "black"

#Filter plotted tree data for internal nodes
data_borneo_recon_internalnode <- data_borneo_recon%>%filter(isTip==FALSE)
data_borneo_recon_internalnode$superfamily <- NA
data_borneo_recon_internalnode$family <- NA
data_borneo_recon_internalnode$colour <- NA


#Remove rows from node assignment data to exclude taxonomy labels/include only colour assignments
borneo_recon_colour_assignments <- borneo_clade_assignments%>%filter(Tree=="Borneo_Reconstructed")%>%filter(ntaxa!=1)%>%select(taxa_1, taxa_2, colour)%>%filter(!is.na(taxa_1))%>%filter(!is.na(taxa_2))%>%filter(!is.na(colour))
#Loop through each row to assign colours to descendants according to clade membership
for (i in 1:nrow(borneo_recon_colour_assignments)) {
  data_borneo_recon_internalnode <- colour_match(data_borneo_recon_internalnode, phy_borneo_recon, borneo_recon_colour_assignments$taxa_1[i], borneo_recon_colour_assignments$taxa_2[i], borneo_recon_colour_assignments$colour[i])
}

#Combine plotted tree data for internal and terminal nodes
plot_borneo_recon_data <- rbind(data_borneo_recon_terminalnode, data_borneo_recon_internalnode)

#Manually assign colours to remaining internal branches
plot_borneo_recon_data[1311,12] <- "black"
plot_borneo_recon_data[1312,12] <- "black"
plot_borneo_recon_data[1313,12] <- "black"
plot_borneo_recon_data[1314,12] <- "black"
plot_borneo_recon_data[1315,12] <- "#46E160"
plot_borneo_recon_data[1386,12] <- "#00E6D3"
plot_borneo_recon_data[1426,12] <- "black"
plot_borneo_recon_data[1427,12] <- "black"
plot_borneo_recon_data[1428,12] <- "black" 
plot_borneo_recon_data[1429,12] <- "black" 
plot_borneo_recon_data[1430,12] <- "#FEDC00" 
plot_borneo_recon_data[1467,12] <- "#C75DE6" 
plot_borneo_recon_data[1481,12] <- "black" 
plot_borneo_recon_data[1482,12] <- "#FF89A7" 
plot_borneo_recon_data[1884,12] <- "#DC3535" 
plot_borneo_recon_data[1932,12] <- "black" 
plot_borneo_recon_data[1933,12] <- "black" 
plot_borneo_recon_data[1934,12] <- "black" 
plot_borneo_recon_data[1935,12] <- "black" 
plot_borneo_recon_data[1936,12] <- "black" 
plot_borneo_recon_data[1937,12] <- "black" 
plot_borneo_recon_data[1938,12] <- "black" 
plot_borneo_recon_data[1939,12] <- "black" 
plot_borneo_recon_data[1940,12] <- "black" 
plot_borneo_recon_data[1941,12] <- "#DEB46B" 
plot_borneo_recon_data[2101,12] <- "#EB4AF0" 
plot_borneo_recon_data[2248,12] <- "#19B519"
plot_borneo_recon_data[2312,12] <- "#19B519"
plot_borneo_recon_data[2347,12] <- "#19B519"
plot_borneo_recon_data[2349,12] <- "#FFBF00"
plot_borneo_recon_data[2455,12] <- "#25D8C5" 
plot_borneo_recon_data[2529,12] <- "#EB8B73"
plot_borneo_recon_data[2573,12] <- "#47F947" 
plot_borneo_recon_data[2580,12] <- "#FF3352" 
plot_borneo_recon_data[2585,12] <- "#FF5DCC" 
plot_borneo_recon_data[2587,12] <- "#FF5DCC" 
plot_borneo_recon_data[2589,12] <- "black"
plot_borneo_recon_data[2590,12] <- "black"
plot_borneo_recon_data[2591,12] <- "#0081FA" 
plot_borneo_recon_data[2614,12] <- "black"
plot_borneo_recon_data[2615,12] <- "#0081FA"

##Plot mirrored reconstructed tree with coloured branches
plot_borneo_recon <- plot_b_r %<+% plot_borneo_recon_data + aes(color=I(colour))


##plot cophylo with connecting lines
#extract data from tree pair
data_borneo_tree1 <- plot_borneo_prune$data
data_borneo_tree2 <- plot_borneo_recon$data
data_borneo_tree1$tree <-'plot_borneo_prune'
data_borneo_tree2$tree <-'plot_borneo_recon'

#shift and mirror one of the tree pair plots to create space for connecting lines
data_borneo_tree2$x <- max(data_borneo_tree2$x) - data_borneo_tree2$x + max(data_borneo_tree1$x) +  max(data_borneo_tree1$x)*0.33

#build cophylo plot
plot_borneo_twotrees <- plot_borneo_prune + geom_tree(data=data_borneo_tree2, linewidth = 0.4)

#extract terminal node data from cophylo plot
data_borneo_twotrees <- bind_rows(data_borneo_tree1, data_borneo_tree2)%>%filter(isTip == TRUE)
data_borneo_twotrees <- as.data.frame(data_borneo_twotrees)

#assign colours to terminal node according to family-level classification
data_borneo_twotrees$colour <- ifelse(data_borneo_twotrees$family=="Artematopodidae", "grey80", ifelse(data_borneo_twotrees$family=="Boridae", "grey80", ifelse(data_borneo_twotrees$family=="Bostrichidae", "grey80", ifelse(data_borneo_twotrees$family=="Brentidae", "grey80", ifelse(data_borneo_twotrees$family=="Buprestidae", "grey80", ifelse(data_borneo_twotrees$family=="Byturidae", "grey80", ifelse(data_borneo_twotrees$family=="Chrysomelidae", "grey80", ifelse(data_borneo_twotrees$family=="Cleridae", "grey80", ifelse(data_borneo_twotrees$family=="Coccinellidae", "grey80", ifelse(data_borneo_twotrees$family=="Elateridae", "grey80", ifelse(data_borneo_twotrees$family=="Elmidae", "grey80", ifelse(data_borneo_twotrees$family=="Gyrinidae", "grey80", ifelse(data_borneo_twotrees$family=="Heteroceridae", "grey80", ifelse(data_borneo_twotrees$family=="Histeridae", "grey80", ifelse(data_borneo_twotrees$family=="Hydraenidae", "grey80", ifelse(data_borneo_twotrees$family=="Hydroscaphidae", "grey80", ifelse(data_borneo_twotrees$family=="Hygrobiidae", "grey80", ifelse(data_borneo_twotrees$family=="Iberobaeniidae", "grey80", ifelse(data_borneo_twotrees$family=="Ischaliidae", "grey80", ifelse(data_borneo_twotrees$family=="Latridiidae", "grey80", ifelse(data_borneo_twotrees$family=="Leiodidae", "grey80", ifelse(data_borneo_twotrees$family=="Limnichidae", "grey80", ifelse(data_borneo_twotrees$family=="Megalopodidae", "grey80", ifelse(data_borneo_twotrees$family=="Meloidae", "grey80", ifelse(data_borneo_twotrees$family=="Mordellidae", "grey80", ifelse(data_borneo_twotrees$family=="Nitidulidae", "grey80", ifelse(data_borneo_twotrees$family=="Noteridae", "grey80", ifelse(data_borneo_twotrees$family=="Oxypeltidae", "grey80", ifelse(data_borneo_twotrees$family=="Phloeostichidae", "grey80", ifelse(data_borneo_twotrees$family=="Platypodidae", "grey80", ifelse(data_borneo_twotrees$family=="Prionoceridae", "grey80", ifelse(data_borneo_twotrees$family=="Rhagophthalmidae", "grey80", ifelse(data_borneo_twotrees$family=="Rhinorhipidae", "grey80", ifelse(data_borneo_twotrees$family=="Salpingidae", "grey80", ifelse(data_borneo_twotrees$family=="Scirtidae", "grey80", ifelse(data_borneo_twotrees$family=="Silvanidae", "grey80", ifelse(data_borneo_twotrees$family=="Tenebrionidae", "grey80", ifelse(data_borneo_twotrees$family=="Vesperidae", "grey80", data_borneo_twotrees$colour))))))))))))))))))))))))))))))))))))))
data_borneo_twotrees$colour <- ifelse(data_borneo_twotrees$family=="Aderidae", "grey60", ifelse(data_borneo_twotrees$family=="Alexiidae", "grey60", ifelse(data_borneo_twotrees$family=="Amphizoidae", "grey60", ifelse(data_borneo_twotrees$family=="Anobiidae", "grey60", ifelse(data_borneo_twotrees$family=="Cantharidae", "grey60", ifelse(data_borneo_twotrees$family=="Carabidae", "grey60", ifelse(data_borneo_twotrees$family=="Chaetosomatidae", "grey60", ifelse(data_borneo_twotrees$family=="Chelonariidae", "grey60", ifelse(data_borneo_twotrees$family=="Corylophidae", "grey60", ifelse(data_borneo_twotrees$family=="Cryptophagidae", "grey60", ifelse(data_borneo_twotrees$family=="Curculionidae", "grey60", ifelse(data_borneo_twotrees$family=="Cybocephalidae", "grey60", ifelse(data_borneo_twotrees$family=="Discolomatidae", "grey60", ifelse(data_borneo_twotrees$family=="Eucinetidae", "grey60", ifelse(data_borneo_twotrees$family=="Georissidae", "grey60", ifelse(data_borneo_twotrees$family=="Geotrupidae", "grey60", ifelse(data_borneo_twotrees$family=="Glaresidae", "grey60", ifelse(data_borneo_twotrees$family=="Haliplidae", "grey60", ifelse(data_borneo_twotrees$family=="Kateretidae", "grey60", ifelse(data_borneo_twotrees$family=="Lepiceridae", "grey60", ifelse(data_borneo_twotrees$family=="Lymexylidae", "grey60", ifelse(data_borneo_twotrees$family=="Melandryidae", "grey60", ifelse(data_borneo_twotrees$family=="Ommatidae", "grey60", ifelse(data_borneo_twotrees$family=="Orsodacnidae", "grey60", ifelse(data_borneo_twotrees$family=="Passandridae", "grey60", ifelse(data_borneo_twotrees$family=="Phalacridae", "grey60", ifelse(data_borneo_twotrees$family=="Protocucujidae", "grey60", ifelse(data_borneo_twotrees$family=="Psephenidae", "grey60", ifelse(data_borneo_twotrees$family=="Rhipiceridae", "grey60", ifelse(data_borneo_twotrees$family=="Scarabaeidae", "grey60", ifelse(data_borneo_twotrees$family=="Silphidae", "grey60", ifelse(data_borneo_twotrees$family=="Spercheidae", "grey60", ifelse(data_borneo_twotrees$family=="Sphindidae", "grey60", ifelse(data_borneo_twotrees$family=="Throscidae", "grey60", ifelse(data_borneo_twotrees$family=="Trictenotomidae", "grey60", ifelse(data_borneo_twotrees$family=="Trogossitidae", "grey60", ifelse(data_borneo_twotrees$family=="Zopheridae", "grey60", data_borneo_twotrees$colour)))))))))))))))))))))))))))))))))))))
data_borneo_twotrees$colour <- ifelse(data_borneo_twotrees$family=="Anthicidae", "grey45", ifelse(data_borneo_twotrees$family=="Anthribidae", "grey45", ifelse(data_borneo_twotrees$family=="Archeocrypticidae", "grey45", ifelse(data_borneo_twotrees$family=="Belidae", "grey45", ifelse(data_borneo_twotrees$family=="Biphyllidae", "grey45", ifelse(data_borneo_twotrees$family=="Bolboceratidae", "grey45", ifelse(data_borneo_twotrees$family=="Byrrhidae", "grey45", ifelse(data_borneo_twotrees$family=="Cerambycidae", "grey45", ifelse(data_borneo_twotrees$family=="Cerylonidae", "grey45", ifelse(data_borneo_twotrees$family=="Clambidae", "grey45", ifelse(data_borneo_twotrees$family=="Cupedidae", "grey45", ifelse(data_borneo_twotrees$family=="Dascillidae", "grey45", ifelse(data_borneo_twotrees$family=="Dermestidae", "grey45", ifelse(data_borneo_twotrees$family=="Dytiscidae", "grey45", ifelse(data_borneo_twotrees$family=="Erirhinidae", "grey45", ifelse(data_borneo_twotrees$family=="Erotylidae", "grey45", ifelse(data_borneo_twotrees$family=="Eucnemidae", "grey45", ifelse(data_borneo_twotrees$family=="Eulichadidae", "grey45", ifelse(data_borneo_twotrees$family=="Glaphyridae", "grey45", ifelse(data_borneo_twotrees$family=="Hydrophilidae", "grey45", ifelse(data_borneo_twotrees$family=="Lucanidae", "grey45", ifelse(data_borneo_twotrees$family=="Meruidae", "grey45", ifelse(data_borneo_twotrees$family=="Oedemeridae", "grey45", ifelse(data_borneo_twotrees$family=="Omalisidae", "grey45", ifelse(data_borneo_twotrees$family=="Omethidae", "grey45", ifelse(data_borneo_twotrees$family=="Phengodidae", "grey45", ifelse(data_borneo_twotrees$family=="Phloiophilidae", "grey45", ifelse(data_borneo_twotrees$family=="Phycosecidae", "grey45", ifelse(data_borneo_twotrees$family=="Propalticidae", "grey45", ifelse(data_borneo_twotrees$family=="Prostomidae", "grey45", ifelse(data_borneo_twotrees$family=="Ptiliidae", "grey45", ifelse(data_borneo_twotrees$family=="Ripiphoridae", "grey45", ifelse(data_borneo_twotrees$family=="Sphaeriusidae", "grey45", ifelse(data_borneo_twotrees$family=="Tetratomidae", "grey45", ifelse(data_borneo_twotrees$family=="Trachypachidae", "grey45", data_borneo_twotrees$colour)))))))))))))))))))))))))))))))))))
data_borneo_twotrees$colour <- ifelse(data_borneo_twotrees$family=="Aspidytidae", "grey30", ifelse(data_borneo_twotrees$family=="Attelabidae", "grey30", ifelse(data_borneo_twotrees$family=="Bothrideridae", "grey30", ifelse(data_borneo_twotrees$family=="Brachyceridae", "grey30", ifelse(data_borneo_twotrees$family=="Callirhipidae", "grey30", ifelse(data_borneo_twotrees$family=="Cerophytidae", "grey30", ifelse(data_borneo_twotrees$family=="Cicindelidae", "grey30", ifelse(data_borneo_twotrees$family=="Ciidae", "grey30", ifelse(data_borneo_twotrees$family=="Cucujidae", "grey30", ifelse(data_borneo_twotrees$family=="Derodontidae", "grey30", ifelse(data_borneo_twotrees$family=="Disteniidae", "grey30", ifelse(data_borneo_twotrees$family=="Dryopidae", "grey30", ifelse(data_borneo_twotrees$family=="Endomychidae", "grey30", ifelse(data_borneo_twotrees$family=="Helophoridae", "grey30", ifelse(data_borneo_twotrees$family=="Helotidae", "grey30", ifelse(data_borneo_twotrees$family=="Hybosoridae", "grey30", ifelse(data_borneo_twotrees$family=="Hydrochidae", "grey30", ifelse(data_borneo_twotrees$family=="Laemophloeidae", "grey30", ifelse(data_borneo_twotrees$family=="Lampyridae", "grey30", ifelse(data_borneo_twotrees$family=="Lycidae", "grey30", ifelse(data_borneo_twotrees$family=="Melyridae", "grey30", ifelse(data_borneo_twotrees$family=="Micromalthidae", "grey30", ifelse(data_borneo_twotrees$family=="Monotomidae", "grey30", ifelse(data_borneo_twotrees$family=="Mycetophagidae", "grey30", ifelse(data_borneo_twotrees$family=="Mycteridae", "grey30", ifelse(data_borneo_twotrees$family=="Nemonychidae", "grey30", ifelse(data_borneo_twotrees$family=="Nosodendridae", "grey30", ifelse(data_borneo_twotrees$family=="Passalidae", "grey30", ifelse(data_borneo_twotrees$family=="Ptilodactylidae", "grey30", ifelse(data_borneo_twotrees$family=="Ptinidae", "grey30", ifelse(data_borneo_twotrees$family=="Pyrochroidae", "grey30", ifelse(data_borneo_twotrees$family=="Scraptiidae", "grey30", ifelse(data_borneo_twotrees$family=="Staphylinidae", "grey30", ifelse(data_borneo_twotrees$family=="Stenotrachelidae", "grey30", ifelse(data_borneo_twotrees$family=="Torridincolidae", "grey30", ifelse(data_borneo_twotrees$family=="Trogidae", "grey30", data_borneo_twotrees$colour))))))))))))))))))))))))))))))))))))
data_borneo_twotrees$colour  <- ifelse(data_borneo_twotrees$colour%in%c("grey30", "grey45", "grey60", "grey80"), data_borneo_twotrees$colour, "grey15")

#combine plots
tanglegram_borneo <- plot_borneo_twotrees + geom_line(aes(x, y, group=label), linewidth = 0.4, data=data_borneo_twotrees) 


##Panama

#Plot pruned tree
plot_p_p <- ggtree(phy_panama_prune, branch.length = 'none', linewidth = 0.4) + 
  theme(plot.margin = unit(c(0.5,0,0.5,0),"cm"))

#Get data of plotted tree
data_panama_prune <- plot_p_p$data
#Get taxonomy assignments for all tips
tip_taxonomy_panama_prune <- site_panama_metadata%>%filter(db_id%in%phy_panama_prune$tip.label)%>%select(db_id, superfamily, family)%>%rename(label = db_id)
#Add taxonomy to plotted tree data for tips
data_panama_prune_terminalnode <- merge(data_panama_prune, tip_taxonomy_panama_prune, by="label")
#Assign colours to branches according to terminal branches according to clade membership
data_panama_prune_terminalnode$colour <- ifelse(data_panama_prune_terminalnode$superfamily == "Caraboidea", "#0081FA", ifelse(data_panama_prune_terminalnode$superfamily == "Dytiscoidea", "#06ACE9", ifelse(data_panama_prune_terminalnode$superfamily == "Gyrinoidea", "#382CF3", ifelse(data_panama_prune_terminalnode$superfamily == "Haliploidea", "#54D3FF", ifelse(data_panama_prune_terminalnode$superfamily == "Bostrichoidea", "#47F947", ifelse(data_panama_prune_terminalnode$superfamily == "Buprestoidea", "#FF3352", ifelse(data_panama_prune_terminalnode$superfamily == "Byrrhoidea", "#F88379", ifelse(data_panama_prune_terminalnode$superfamily == "Chrysomeloidea", "#EB4AF0", ifelse(data_panama_prune_terminalnode$superfamily == "Cleroidea", "#EB8B73", ifelse(data_panama_prune_terminalnode$superfamily == "Coccinelloidea", "#25D8C5", ifelse(data_panama_prune_terminalnode$superfamily == "Cucujoidea", "#19B519", ifelse(data_panama_prune_terminalnode$superfamily == "Curculionoidea", "#DEB46B", ifelse(data_panama_prune_terminalnode$superfamily == "Dascilloidea", "#A6E52E", ifelse(data_panama_prune_terminalnode$superfamily == "Derodontoidea", "#E03EF5", ifelse(data_panama_prune_terminalnode$superfamily == "Dryopoidea", "#00E6D3", ifelse(data_panama_prune_terminalnode$superfamily == "Elateroidea", "#46E160", ifelse(data_panama_prune_terminalnode$superfamily == "Histeroidea", "#FEDC00", ifelse(data_panama_prune_terminalnode$superfamily == "Hydrophiloidea", "#C75DE6", ifelse(data_panama_prune_terminalnode$superfamily == "Lymexyloidea", "#FF359A", ifelse(data_panama_prune_terminalnode$superfamily == "Rhinorhipoidea", "#FFF747", ifelse(data_panama_prune_terminalnode$superfamily == "Scarabaeoidea", "#DC3535", ifelse(data_panama_prune_terminalnode$superfamily == "Scirtoidea", "#FF5DCC", ifelse(data_panama_prune_terminalnode$superfamily == "Staphylinoidea", "#FF89A7", ifelse(data_panama_prune_terminalnode$superfamily == "Tenebrionoidea", "#FFBF00", NA))))))))))))))))))))))))
data_panama_prune_terminalnode$colour[is.na(data_panama_prune_terminalnode$colour)] <- "black"

#Filter plotted tree data for internal nodes
data_panama_prune_internalnode <- data_panama_prune%>%filter(isTip==FALSE)
data_panama_prune_internalnode$superfamily <- NA
data_panama_prune_internalnode$family <- NA
data_panama_prune_internalnode$colour <- NA


#Remove rows from node assignment data to exclude taxonomy labels/include only colour assignments
panama_prune_colour_assignments <- panama_clade_assignments%>%filter(Tree=="Panama_Prune")%>%filter(ntaxa!=1)%>%select(taxa_1, taxa_2, colour)%>%filter(!is.na(taxa_1))%>%filter(!is.na(taxa_2))%>%filter(!is.na(colour))
#Loop through each row to assign colours to descendants according to clade membership
for (i in 1:nrow(panama_prune_colour_assignments)) {
  data_panama_prune_internalnode <- colour_match(data_panama_prune_internalnode, phy_panama_prune, panama_prune_colour_assignments$taxa_1[i], panama_prune_colour_assignments$taxa_2[i], panama_prune_colour_assignments$colour[i])
}

#Combine plotted tree data for internal and terminal nodes
plot_panama_prune_data <- rbind(data_panama_prune_terminalnode, data_panama_prune_internalnode)

#Manually assign colours to remaining internal branches
plot_panama_prune_data[1101,12] <- "black"
plot_panama_prune_data[1102,12] <- "black"
plot_panama_prune_data[1103,12] <- "black"
plot_panama_prune_data[1104,12] <- "black"
plot_panama_prune_data[1105,12] <- "black"
plot_panama_prune_data[1106,12] <- "black"
plot_panama_prune_data[1107,12] <- "#46E160"
plot_panama_prune_data[1179,12] <- "#FF3352"
plot_panama_prune_data[1181,12] <- "#00E6D3"
plot_panama_prune_data[1192,12] <- "black"
plot_panama_prune_data[1193,12] <- "black"
plot_panama_prune_data[1194,12] <- "black"
plot_panama_prune_data[1195,12] <- "#FEDC00"
plot_panama_prune_data[1218,12] <- "#C75DE6"
plot_panama_prune_data[1230,12] <- "black"
plot_panama_prune_data[1231,12] <- "#DC3535"
plot_panama_prune_data[1264,12] <- "#FF89A7"
plot_panama_prune_data[1427,12] <- "black"
plot_panama_prune_data[1428,12] <- "black"
plot_panama_prune_data[1429,12] <- "black"
plot_panama_prune_data[1430,12] <- "black"
plot_panama_prune_data[1431,12] <- "black"
plot_panama_prune_data[1432,12] <- "black"
plot_panama_prune_data[1433,12] <- "#DEB46B"
plot_panama_prune_data[1731,12] <- "#EB4AF0"
plot_panama_prune_data[1885,12] <- "#19B519"
plot_panama_prune_data[1978,12] <- "black"
plot_panama_prune_data[1979,12] <- "#FFBF00"
plot_panama_prune_data[2088,12] <- "#25D8C5"
plot_panama_prune_data[2128,12] <- "#EB8B73"
plot_panama_prune_data[2147,12] <- "#47F947"
plot_panama_prune_data[2163,12] <- "#FF5DCC"
plot_panama_prune_data[2168,12] <- "black"
plot_panama_prune_data[2169,12] <- "#0081FA"
##Plot pruned tree with coloured branches
plot_panama_prune <- plot_p_p %<+% plot_panama_prune_data + aes(color=I(colour))


#Plot reconstructed tree mirrored in y-axis
#Plot pruned tree
plot_p_r <- ggtree(phy_panama_recon, branch.length = 'none', linewidth = 0.4) + 
  scale_x_reverse() + 
  theme(plot.margin = unit(c(0.5,0,0.5,0),"cm"))

#Get data of plotted tree
data_panama_recon <- plot_p_r$data
#Get taxonomy assignments for all tips
tip_taxonomy_panama_recon <- site_panama_metadata%>%filter(db_id%in%phy_panama_recon$tip.label)%>%select(db_id, superfamily, family)%>%rename(label = db_id)
#Add taxonomy to plotted tree data for tips
data_panama_recon_terminalnode <- merge(data_panama_recon, tip_taxonomy_panama_recon, by="label")
#Assign colours to branches according to terminal branches according to clade membership
data_panama_recon_terminalnode$colour <- ifelse(data_panama_recon_terminalnode$superfamily == "Caraboidea", "#0081FA", ifelse(data_panama_recon_terminalnode$superfamily == "Dytiscoidea", "#06ACE9", ifelse(data_panama_recon_terminalnode$superfamily == "Gyrinoidea", "#382CF3", ifelse(data_panama_recon_terminalnode$superfamily == "Haliploidea", "#54D3FF", ifelse(data_panama_recon_terminalnode$superfamily == "Bostrichoidea", "#47F947", ifelse(data_panama_recon_terminalnode$superfamily == "Buprestoidea", "#FF3352", ifelse(data_panama_recon_terminalnode$superfamily == "Byrrhoidea", "#F88379", ifelse(data_panama_recon_terminalnode$superfamily == "Chrysomeloidea", "#EB4AF0", ifelse(data_panama_recon_terminalnode$superfamily == "Cleroidea", "#EB8B73", ifelse(data_panama_recon_terminalnode$superfamily == "Coccinelloidea", "#25D8C5", ifelse(data_panama_recon_terminalnode$superfamily == "Cucujoidea", "#19B519", ifelse(data_panama_recon_terminalnode$superfamily == "Curculionoidea", "#DEB46B", ifelse(data_panama_recon_terminalnode$superfamily == "Dascilloidea", "#A6E52E", ifelse(data_panama_recon_terminalnode$superfamily == "Derodontoidea", "#E03EF5", ifelse(data_panama_recon_terminalnode$superfamily == "Dryopoidea", "#00E6D3", ifelse(data_panama_recon_terminalnode$superfamily == "Elateroidea", "#46E160", ifelse(data_panama_recon_terminalnode$superfamily == "Histeroidea", "#FEDC00", ifelse(data_panama_recon_terminalnode$superfamily == "Hydrophiloidea", "#C75DE6", ifelse(data_panama_recon_terminalnode$superfamily == "Lymexyloidea", "#FF359A", ifelse(data_panama_recon_terminalnode$superfamily == "Rhinorhipoidea", "#FFF747", ifelse(data_panama_recon_terminalnode$superfamily == "Scarabaeoidea", "#DC3535", ifelse(data_panama_recon_terminalnode$superfamily == "Scirtoidea", "#FF5DCC", ifelse(data_panama_recon_terminalnode$superfamily == "Staphylinoidea", "#FF89A7", ifelse(data_panama_recon_terminalnode$superfamily == "Tenebrionoidea", "#FFBF00", NA))))))))))))))))))))))))
data_panama_recon_terminalnode$colour[is.na(data_panama_recon_terminalnode$colour)] <- "black"

#Filter plotted tree data for internal nodes
data_panama_recon_internalnode <- data_panama_recon%>%filter(isTip==FALSE)
data_panama_recon_internalnode$superfamily <- NA
data_panama_recon_internalnode$family <- NA
data_panama_recon_internalnode$colour <- NA


#Remove rows from node assignment data to exclude taxonomy labels/include only colour assignments
panama_recon_colour_assignments <- panama_clade_assignments%>%filter(Tree=="Panama_Reconstructed")%>%filter(ntaxa!=1)%>%select(taxa_1, taxa_2, colour)%>%filter(!is.na(taxa_1))%>%filter(!is.na(taxa_2))%>%filter(!is.na(colour))
#Loop through each row to assign colours to descendants according to clade membership
for (i in 1:nrow(panama_recon_colour_assignments)) {
  data_panama_recon_internalnode <- colour_match(data_panama_recon_internalnode, phy_panama_recon, panama_recon_colour_assignments$taxa_1[i], panama_recon_colour_assignments$taxa_2[i], panama_recon_colour_assignments$colour[i])
}

#Combine plotted tree data for internal and terminal nodes
plot_panama_recon_data <- rbind(data_panama_recon_terminalnode, data_panama_recon_internalnode)

#Manually assign colours to remaining internal branches
plot_panama_recon_data[1101,12] <- "black"
plot_panama_recon_data[1102,12] <- "black"
plot_panama_recon_data[1103,12] <- "black"
plot_panama_recon_data[1104,12] <- "black"
plot_panama_recon_data[1105,12] <- "#46E160"
plot_panama_recon_data[1177,12] <- "black"
plot_panama_recon_data[1178,12] <- "black"
plot_panama_recon_data[1179,12] <- "black"
plot_panama_recon_data[1180,12] <- "black"
plot_panama_recon_data[1181,12] <- "black"
plot_panama_recon_data[1182,12] <- "black"
plot_panama_recon_data[1183,12] <- "#C75DE6"
plot_panama_recon_data[1195,12] <- "#FEDC00"
plot_panama_recon_data[1218,12] <- "#FF89A7"
plot_panama_recon_data[1233,12] <- "black"
plot_panama_recon_data[1234,12] <- "#DC3535"
plot_panama_recon_data[1267,12] <- "#FF89A7"
plot_panama_recon_data[1413,12] <- "black"
plot_panama_recon_data[1414,12] <- "black"
plot_panama_recon_data[1415,12] <- "black"
plot_panama_recon_data[1416,12] <- "black"
plot_panama_recon_data[1417,12] <- "black"
plot_panama_recon_data[1418,12] <- "black"
plot_panama_recon_data[1419,12] <- "black"
plot_panama_recon_data[1420,12] <- "black"
plot_panama_recon_data[1421,12] <- "#DEB46B"
plot_panama_recon_data[1719,12] <- "#EB4AF0"
plot_panama_recon_data[1874,12] <- "#19B519"
plot_panama_recon_data[1942,12] <- "#19B519"
plot_panama_recon_data[1964,12] <- "#25D8C5"
plot_panama_recon_data[2005,12] <- "black"
plot_panama_recon_data[2006,12] <- "#FFBF00"
plot_panama_recon_data[2116,12] <- "#EB8B73"
plot_panama_recon_data[2135,12] <- "#47F947"
plot_panama_recon_data[2151,12] <- "#FF3352"
plot_panama_recon_data[2153,12] <- "#00E6D3"
plot_panama_recon_data[2163,12] <- "#FF5DCC"
plot_panama_recon_data[2168,12] <- "black"
plot_panama_recon_data[2169,12] <- "#0081FA"

##Plot mirrored reconstructed tree with coloured branches
plot_panama_recon <- plot_p_r %<+% plot_panama_recon_data + aes(color=I(colour))

##plot cophylo with connecting lines
#extract data from tree pair
data_panama_tree1 <- plot_panama_prune$data
data_panama_tree2 <- plot_panama_recon$data
data_panama_tree1$tree <-'plot_panama_prune'
data_panama_tree2$tree <-'plot_panama_recon'

#shift and mirror one of the tree pair plots to create space for connecting lines
data_panama_tree2$x <- max(data_panama_tree2$x) - data_panama_tree2$x + max(data_panama_tree1$x) +  max(data_panama_tree1$x)*0.33

#build cophylo plot
plot_panama_twotrees <- plot_panama_prune + geom_tree(data=data_panama_tree2, linewidth = 0.4)

#extract terminal node data from cophylo plot
data_panama_twotrees <- bind_rows(data_panama_tree1, data_panama_tree2)%>%filter(isTip == TRUE)
data_panama_twotrees <- as.data.frame(data_panama_twotrees)

#assign colours to terminal node according to family-level classification
data_panama_twotrees$colour <- ifelse(data_panama_twotrees$family=="Artematopodidae", "grey80", ifelse(data_panama_twotrees$family=="Boridae", "grey80", ifelse(data_panama_twotrees$family=="Bostrichidae", "grey80", ifelse(data_panama_twotrees$family=="Brentidae", "grey80", ifelse(data_panama_twotrees$family=="Buprestidae", "grey80", ifelse(data_panama_twotrees$family=="Byturidae", "grey80", ifelse(data_panama_twotrees$family=="Chrysomelidae", "grey80", ifelse(data_panama_twotrees$family=="Cleridae", "grey80", ifelse(data_panama_twotrees$family=="Coccinellidae", "grey80", ifelse(data_panama_twotrees$family=="Elateridae", "grey80", ifelse(data_panama_twotrees$family=="Elmidae", "grey80", ifelse(data_panama_twotrees$family=="Gyrinidae", "grey80", ifelse(data_panama_twotrees$family=="Heteroceridae", "grey80", ifelse(data_panama_twotrees$family=="Histeridae", "grey80", ifelse(data_panama_twotrees$family=="Hydraenidae", "grey80", ifelse(data_panama_twotrees$family=="Hydroscaphidae", "grey80", ifelse(data_panama_twotrees$family=="Hygrobiidae", "grey80", ifelse(data_panama_twotrees$family=="Iberobaeniidae", "grey80", ifelse(data_panama_twotrees$family=="Ischaliidae", "grey80", ifelse(data_panama_twotrees$family=="Latridiidae", "grey80", ifelse(data_panama_twotrees$family=="Leiodidae", "grey80", ifelse(data_panama_twotrees$family=="Limnichidae", "grey80", ifelse(data_panama_twotrees$family=="Megalopodidae", "grey80", ifelse(data_panama_twotrees$family=="Meloidae", "grey80", ifelse(data_panama_twotrees$family=="Mordellidae", "grey80", ifelse(data_panama_twotrees$family=="Nitidulidae", "grey80", ifelse(data_panama_twotrees$family=="Noteridae", "grey80", ifelse(data_panama_twotrees$family=="Oxypeltidae", "grey80", ifelse(data_panama_twotrees$family=="Phloeostichidae", "grey80", ifelse(data_panama_twotrees$family=="Platypodidae", "grey80", ifelse(data_panama_twotrees$family=="Prionoceridae", "grey80", ifelse(data_panama_twotrees$family=="Rhagophthalmidae", "grey80", ifelse(data_panama_twotrees$family=="Rhinorhipidae", "grey80", ifelse(data_panama_twotrees$family=="Salpingidae", "grey80", ifelse(data_panama_twotrees$family=="Scirtidae", "grey80", ifelse(data_panama_twotrees$family=="Silvanidae", "grey80", ifelse(data_panama_twotrees$family=="Tenebrionidae", "grey80", ifelse(data_panama_twotrees$family=="Vesperidae", "grey80", data_panama_twotrees$colour))))))))))))))))))))))))))))))))))))))
data_panama_twotrees$colour <- ifelse(data_panama_twotrees$family=="Aderidae", "grey60", ifelse(data_panama_twotrees$family=="Alexiidae", "grey60", ifelse(data_panama_twotrees$family=="Amphizoidae", "grey60", ifelse(data_panama_twotrees$family=="Anobiidae", "grey60", ifelse(data_panama_twotrees$family=="Cantharidae", "grey60", ifelse(data_panama_twotrees$family=="Carabidae", "grey60", ifelse(data_panama_twotrees$family=="Chaetosomatidae", "grey60", ifelse(data_panama_twotrees$family=="Chelonariidae", "grey60", ifelse(data_panama_twotrees$family=="Corylophidae", "grey60", ifelse(data_panama_twotrees$family=="Cryptophagidae", "grey60", ifelse(data_panama_twotrees$family=="Curculionidae", "grey60", ifelse(data_panama_twotrees$family=="Cybocephalidae", "grey60", ifelse(data_panama_twotrees$family=="Discolomatidae", "grey60", ifelse(data_panama_twotrees$family=="Eucinetidae", "grey60", ifelse(data_panama_twotrees$family=="Georissidae", "grey60", ifelse(data_panama_twotrees$family=="Geotrupidae", "grey60", ifelse(data_panama_twotrees$family=="Glaresidae", "grey60", ifelse(data_panama_twotrees$family=="Haliplidae", "grey60", ifelse(data_panama_twotrees$family=="Kateretidae", "grey60", ifelse(data_panama_twotrees$family=="Lepiceridae", "grey60", ifelse(data_panama_twotrees$family=="Lymexylidae", "grey60", ifelse(data_panama_twotrees$family=="Melandryidae", "grey60", ifelse(data_panama_twotrees$family=="Ommatidae", "grey60", ifelse(data_panama_twotrees$family=="Orsodacnidae", "grey60", ifelse(data_panama_twotrees$family=="Passandridae", "grey60", ifelse(data_panama_twotrees$family=="Phalacridae", "grey60", ifelse(data_panama_twotrees$family=="Protocucujidae", "grey60", ifelse(data_panama_twotrees$family=="Psephenidae", "grey60", ifelse(data_panama_twotrees$family=="Rhipiceridae", "grey60", ifelse(data_panama_twotrees$family=="Scarabaeidae", "grey60", ifelse(data_panama_twotrees$family=="Silphidae", "grey60", ifelse(data_panama_twotrees$family=="Spercheidae", "grey60", ifelse(data_panama_twotrees$family=="Sphindidae", "grey60", ifelse(data_panama_twotrees$family=="Throscidae", "grey60", ifelse(data_panama_twotrees$family=="Trictenotomidae", "grey60", ifelse(data_panama_twotrees$family=="Trogossitidae", "grey60", ifelse(data_panama_twotrees$family=="Zopheridae", "grey60", data_panama_twotrees$colour)))))))))))))))))))))))))))))))))))))
data_panama_twotrees$colour <- ifelse(data_panama_twotrees$family=="Anthicidae", "grey45", ifelse(data_panama_twotrees$family=="Anthribidae", "grey45", ifelse(data_panama_twotrees$family=="Archeocrypticidae", "grey45", ifelse(data_panama_twotrees$family=="Belidae", "grey45", ifelse(data_panama_twotrees$family=="Biphyllidae", "grey45", ifelse(data_panama_twotrees$family=="Bolboceratidae", "grey45", ifelse(data_panama_twotrees$family=="Byrrhidae", "grey45", ifelse(data_panama_twotrees$family=="Cerambycidae", "grey45", ifelse(data_panama_twotrees$family=="Cerylonidae", "grey45", ifelse(data_panama_twotrees$family=="Clambidae", "grey45", ifelse(data_panama_twotrees$family=="Cupedidae", "grey45", ifelse(data_panama_twotrees$family=="Dascillidae", "grey45", ifelse(data_panama_twotrees$family=="Dermestidae", "grey45", ifelse(data_panama_twotrees$family=="Dytiscidae", "grey45", ifelse(data_panama_twotrees$family=="Erirhinidae", "grey45", ifelse(data_panama_twotrees$family=="Erotylidae", "grey45", ifelse(data_panama_twotrees$family=="Eucnemidae", "grey45", ifelse(data_panama_twotrees$family=="Eulichadidae", "grey45", ifelse(data_panama_twotrees$family=="Glaphyridae", "grey45", ifelse(data_panama_twotrees$family=="Hydrophilidae", "grey45", ifelse(data_panama_twotrees$family=="Lucanidae", "grey45", ifelse(data_panama_twotrees$family=="Meruidae", "grey45", ifelse(data_panama_twotrees$family=="Oedemeridae", "grey45", ifelse(data_panama_twotrees$family=="Omalisidae", "grey45", ifelse(data_panama_twotrees$family=="Omethidae", "grey45", ifelse(data_panama_twotrees$family=="Phengodidae", "grey45", ifelse(data_panama_twotrees$family=="Phloiophilidae", "grey45", ifelse(data_panama_twotrees$family=="Phycosecidae", "grey45", ifelse(data_panama_twotrees$family=="Propalticidae", "grey45", ifelse(data_panama_twotrees$family=="Prostomidae", "grey45", ifelse(data_panama_twotrees$family=="Ptiliidae", "grey45", ifelse(data_panama_twotrees$family=="Ripiphoridae", "grey45", ifelse(data_panama_twotrees$family=="Sphaeriusidae", "grey45", ifelse(data_panama_twotrees$family=="Tetratomidae", "grey45", ifelse(data_panama_twotrees$family=="Trachypachidae", "grey45", data_panama_twotrees$colour)))))))))))))))))))))))))))))))))))
data_panama_twotrees$colour <- ifelse(data_panama_twotrees$family=="Aspidytidae", "grey30", ifelse(data_panama_twotrees$family=="Attelabidae", "grey30", ifelse(data_panama_twotrees$family=="Bothrideridae", "grey30", ifelse(data_panama_twotrees$family=="Brachyceridae", "grey30", ifelse(data_panama_twotrees$family=="Callirhipidae", "grey30", ifelse(data_panama_twotrees$family=="Cerophytidae", "grey30", ifelse(data_panama_twotrees$family=="Cicindelidae", "grey30", ifelse(data_panama_twotrees$family=="Ciidae", "grey30", ifelse(data_panama_twotrees$family=="Cucujidae", "grey30", ifelse(data_panama_twotrees$family=="Derodontidae", "grey30", ifelse(data_panama_twotrees$family=="Disteniidae", "grey30", ifelse(data_panama_twotrees$family=="Dryopidae", "grey30", ifelse(data_panama_twotrees$family=="Endomychidae", "grey30", ifelse(data_panama_twotrees$family=="Helophoridae", "grey30", ifelse(data_panama_twotrees$family=="Helotidae", "grey30", ifelse(data_panama_twotrees$family=="Hybosoridae", "grey30", ifelse(data_panama_twotrees$family=="Hydrochidae", "grey30", ifelse(data_panama_twotrees$family=="Laemophloeidae", "grey30", ifelse(data_panama_twotrees$family=="Lampyridae", "grey30", ifelse(data_panama_twotrees$family=="Lycidae", "grey30", ifelse(data_panama_twotrees$family=="Melyridae", "grey30", ifelse(data_panama_twotrees$family=="Micromalthidae", "grey30", ifelse(data_panama_twotrees$family=="Monotomidae", "grey30", ifelse(data_panama_twotrees$family=="Mycetophagidae", "grey30", ifelse(data_panama_twotrees$family=="Mycteridae", "grey30", ifelse(data_panama_twotrees$family=="Nemonychidae", "grey30", ifelse(data_panama_twotrees$family=="Nosodendridae", "grey30", ifelse(data_panama_twotrees$family=="Passalidae", "grey30", ifelse(data_panama_twotrees$family=="Ptilodactylidae", "grey30", ifelse(data_panama_twotrees$family=="Ptinidae", "grey30", ifelse(data_panama_twotrees$family=="Pyrochroidae", "grey30", ifelse(data_panama_twotrees$family=="Scraptiidae", "grey30", ifelse(data_panama_twotrees$family=="Staphylinidae", "grey30", ifelse(data_panama_twotrees$family=="Stenotrachelidae", "grey30", ifelse(data_panama_twotrees$family=="Torridincolidae", "grey30", ifelse(data_panama_twotrees$family=="Trogidae", "grey30", data_panama_twotrees$colour))))))))))))))))))))))))))))))))))))
data_panama_twotrees$colour  <- ifelse(data_panama_twotrees$colour%in%c("grey30", "grey45", "grey60", "grey80"), data_panama_twotrees$colour, "grey15")

#combine plots
tanglegram_panama <- plot_panama_twotrees + geom_line(aes(x, y, group=label), linewidth = 0.4, data=data_panama_twotrees) 


##create legend:
#build an example plot - complete overlap of suborders and superfamilies so either site could be chosen!
legend_plot <- ggtree(phy_panama_prune, branch.length='none') + geom_tiplab() + coord_cartesian(clip = "off") + vexpand(.0005, 1)
legend_metadata <- site_panama_metadata
legend_metadata$superfamily[legend_metadata$superfamily==""] <- "NA"
legend_metadata <- legend_metadata%>%select(db_id, superfamily)%>%column_to_rownames(var = "db_id")
legend_colours <- c("Caraboidea" = "#0081FA", "Dytiscoidea" = "#06ACE9", "Gyrinoidea" = "#382CF3", "Haliploidea" = "#54D3FF", "Bostrichoidea" = "#47F947", "Buprestoidea" = "#FF3352", "Byrrhoidea" = "#F88379", "Chrysomeloidea" = "#EB4AF0", "Cleroidea" = "#EB8B73", "Coccinelloidea" = "#25D8C5", "Cucujoidea" = "#19B519", "Curculionoidea" = "#DEB46B", "Dascilloidea" = "#A6E52E", "Derodontoidea" = "#E03EF5", "Dryopoidea" = "#00E6D3", "Elateroidea" = "#46E160", "Histeroidea" = "#FEDC00", "Hydrophiloidea" = "#C75DE6", "Lymexyloidea" = "#FF359A", "Rhinorhipoidea" = "#FFF747", "Scarabaeoidea" = "#DC3535", "Scirtoidea" = "#FF5DCC", "Staphylinoidea" = "#FF89A7", "Tenebrionoidea" = "#FFBF00", "NA" = "black")
legend_plot2 <- gheatmap(legend_plot, legend_metadata, offset=2.7, width = 0.0125, color = NULL, colnames = F) +   
  scale_fill_manual(values = legend_colours, breaks = c("Bostrichoidea", "Buprestoidea", "Caraboidea", "Chrysomeloidea", "Cleroidea", "Coccinelloidea", "Cucujoidea", "Curculionoidea", "Derodontoidea", "Dryopoidea", "Dytiscoidea", "Elateroidea", "Histeroidea", "Hydrophiloidea", "Lymexyloidea", "Scarabaeoidea", "Scirtoidea", "Staphylinoidea", "Tenebrionoidea", "NA"), na.value = NA, name="Superfamily") +
  guides(fill=guide_legend(title.position="top", ncol=1)) + 
  theme(legend.position="right", 
        legend.justification = "top", 
        legend.title = element_text(size=15),
        legend.text = element_text(size=13))

#extract legend
plot_legend <- get_legend(legend_plot2 + theme(legend.box.margin = margin(30, 0, 0, -50)))

#plot tanglegrams side by side, including legend
plot_tanglegram <- plot_grid(tanglegram_borneo, tanglegram_panama, plot_legend, nrow = 1, align = "v", labels = c("(a)", "(b)", ""), rel_widths = c(5.5,5.5,1), label_size = 15, hjust=-0.3, vjust=2)

ggsave("Tanglegram_plot.svg", plot_tanglegram, height=18, width=13, limitsize = F)


dev.off()
