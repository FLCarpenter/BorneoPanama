# BorneoPanama
Bioinformatics scripts and manuscript figures for *If the tape were played again: Lineage evolution and the causes of phylogenetic similarity in two tropical assemblages of Coleoptera*

Scripts were run locally using R v4.3.2 in RStudio v2023.09.1+494, with the exception of ```Calculation_LRPD.R``` which was run on an HPC server.

## 1. Evaluating_Trees
```Analyses_site_tree_comparison.R```
calculates various measures of taxonomic congruence and topological comparison (tRI analyses carried out using [taxonomicindices.R](https://github.com/tjcreedy/phylostuff/blob/main/taxonomicindices.R); tSDI clusters manually identified)

```Plot_4184_tree.R``` plots the full phylogeny with heatmap

```Plot_family_tree.R``` plots the pruned phylogenies for the 3 most species-rich families with heatmap

```Plot_site_tree_comparison.R``` plots co-phylogenies comparing two tree building methods for each of the Borneo and Panama sites


## 2. Diversity_Distributions
```Analysis_Plot_species_richness_correlation.R``` Compares species richness of families found across both sites

```Analysis_phylodiversity_structure.R``` Calculates phylogenetic diversity indices to detect phylogenetic structuring of each site

```Analysis_Plot_pd_rarefaction.R``` Sample-size-based rarefaction to standardise for unequal richness across samples

```Analysis_Plot_lineage_pd_correlation``` Compares phylogenetic structure and measures of PD at the family-level between sites

```Data_Preparation_LRPD.R``` Set up matrices required for LRPD analysis

```Calculation_LRPD.R``` Calculates cumulative phylogenetic diversity indices as decreasingly species-rich families are included in analysis for each site

```Plot_LRPD.R``` Plots the LRPD results and calculates breakpoint significance
