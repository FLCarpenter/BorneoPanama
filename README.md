# BorneoPanama
Scripts for *If the tape were played again: Lineage evolution and the causes of phylogenetic similarity in two tropical assemblages of Coleoptera*

R scripts were run locally using R v4.3.2 in RStudio v2023.09.1+494, with the exception of ```Calculation_LRPD.R``` which was run on an HPC server.

## 1. Evaluating_Trees
```Analyses_site_tree_comparison.R```
calculates various measures of taxonomic congruence and topological comparison (tRI analyses carried out using [taxonomicindices.R](https://github.com/tjcreedy/phylostuff/blob/main/taxonomicindices.R); tSDI clusters manually identified)

```Plot_4184_tree.R``` plots the full phylogeny with heatmap

```Plot_family_tree.R``` plots the pruned phylogenies for the 3 most species-rich families with heatmap

```Plot_site_tree_comparison.R``` plots a tanglegram/cophylogeny for the Borneo and Panama sites comparing two tree building methods for each site


## 2. Diversity_Distributions
