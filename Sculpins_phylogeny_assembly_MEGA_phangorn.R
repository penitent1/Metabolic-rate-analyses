library(ctv)
library(tidyverse)
library(ape)
library(rotl)
library(geiger)
library(ggthemes)
library(MCMCglmm)
library(kinship2)
library(phytools)
library(phangorn)
library(nlme)
library(visreg)
library(RColorBrewer)
library(caper)
library(Rphylip)



## Estimating a maximum likelihood phylogenetic tree

## Read in aligment data using phangorn

#sculpins <- read.phyDat(file.choose(), format = "nexus")
test <- read.phyDat("renamed_sculpin_practice_CYTB_aligned_new_format.fas", format = "fasta")

dm_test <- dist.ml(test)
treeUPGMA_test <- upgma(dm_test) # rooted tree
treeNJ_test <- NJ(dm_test) # unrooted tree
plot.phylo(treeUPGMA_test)
plot.phylo(treeNJ_test)
