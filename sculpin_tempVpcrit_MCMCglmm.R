library(ctv)
library(tidyverse)
library(ape)
library(rotl)
library(geiger)
library(ggthemes)
library(MCMCglmm)
library(kinship2)
library(phytools)

#####################################################################
## Read sculpin phyologeny into R and prune for species in analysis
#####################################################################

knope <- read.tree("knope_rooted_ultrametric")
#knope.tree.new <- read.nexus("bayes_tree_final_NEW.nex")
#knope2 <- read.nexus("bayes_tree_final_knope.nex")

is.ultrametric(knope)
is.rooted(knope)

plot.phylo(knope, cex = 0.5)

## Species I have Pcrit data for, plus one outgroup to keep tree rooted
keepers <- c("Oligocottus_maculosus",
             "Clinocottus_globiceps",
             "Artedius_harringtoni",
             "Artedius_lateralis",
             "Artedius_fenestralis",
             "Blepsias_cirrhosus",
             "Scorpaenichthys_marmoratus",
             "Enophrys_bison",
             "Hemilepidotus_hemilepidotus")

## Prune phylogeny for species in analysis
knope_pruned <- drop.tip(knope, setdiff(knope$tip.label, keepers))
#knope_pruned_root <- root(knope_pruned, "Hemilepidotus_hemilepidotus", resolve.root = TRUE)
#knope_pruned_root <- multi2di(knope_pruned_root)
#knope_pruned_root <- chronopl(knope_pruned_root, 1)

rr.62<-reroot(tree,62,position=0.5*tree$edge.length[which(tree$edge[,2]==62)])
plotTree(rr.62)

reroot<-function(tree,node.number,position){
  # first, re-root the tree at node.number
  tr<-root(tree,node=node.number,resolve.root=T)
  # now re-allocate branch length to the two edges descending
  # from the new root node
  b<-sum(tr$edge.length[tr$edge==(length(tree$tip)+1)])
  tr$edge.length[tr$edge==(length(tree$tip)+1)]<-c(position,b-position)
  return(tr)
}

tr_knope <- root(knope_pruned, , resolve.root = T)
knope_pruned_root <- reroot(knope_pruned, 10, "Scorpaenichthys_marmoratus")
# Hemilepidotus_hemilepidotus


## Check the structure of the "relatedness" or A var-covar matrix
A_knope <- vcv.phylo(knope_pruned, cor = T)
# A_knope_test <- vcv.phylo(knope_pruned) ## Weird, not scaled to 1,
  ## Need to use "corr = TRUE" argument to make sure relatedness values
  ## are scaled to 1 for MCMCglmm analysis

############################################
## Read sculpin pcrit and meta-data into R
############################################

md <- read.csv(file.choose(), stringsAsFactors = FALSE)
md <- as_tibble(md)

#
md$smr.best[md$smr.best == "na"] <- NA
md$smr.best <- as.double(md$smr.best)
#
md$pcrit.r[md$pcrit.r == "na"] <- NA
md$pcrit.r <- as.double(md$pcrit.r)
#
md$mo2.regress[md$mo2.regress == "na"] <- NA
md$mo2.regress <- as.double(md$mo2.regress)
#
md$pcrit.regress[md$pcrit.regress == "na"] <- NA
md$pcrit.regress <- as.double(md$pcrit.regress)
#
md$pcrit.time.min[md$pcrit.time.min == "na"] <- NA
md$pcrit.time.min <- as.double(md$pcrit.time.min)

plot_pcrit <- ggplot(data = md) +
  geom_point(mapping = aes(x = temp, y = pcrit.r, color = spps))

## Pcrit and temp data, NA's filtered out
md_pcrit_smr <- md %>%
  select(spps, temp, fish.id, mass.g, smr.best, pcrit.r) %>%
  filter(!is.na(smr.best)) %>%
  filter(!is.na(pcrit.r))

## "species" column added for joining up phylogeny and data set
species <- vector()

for(i in 1:length(md_pcrit_smr$spps)){
  if (md_pcrit_smr$spps[i] == "olma"){
    md_pcrit_smr$species[i] <- "Oligocottus_maculosus"
    } else {
  if (md_pcrit_smr$spps[i] == "clgl"){
    md_pcrit_smr$species[i] <- "Clinocottus_globiceps"
    } else {
  if (md_pcrit_smr$spps[i] == "blci"){
    md_pcrit_smr$species[i] <- "Blepsias_cirrhosus"
    } else {
  if (md_pcrit_smr$spps[i] == "arha"){
    md_pcrit_smr$species[i] <- "Artedius_harringtoni"
    } else {
  if (md_pcrit_smr$spps[i] == "arfe"){
    md_pcrit_smr$species[i] <- "Artedius_fenestralis"
    } else {
  if (md_pcrit_smr$spps[i] == "arla"){
    md_pcrit_smr$species[i] <- "Artedius_lateralis"
    } else {
  if (md_pcrit_smr$spps[i] == "hehe"){
    md_pcrit_smr$species[i] <- "Hemilepidotus_hemilepidotus"
    } else {
  if (md_pcrit_smr$spps[i] == "scma"){
    md_pcrit_smr$species[i] <- "Scorpaenichthys_marmoratus"
    } else {
  if (md_pcrit_smr$spps[i] == "enbi"){
    md_pcrit_smr$species[i] <- "Enophrys_bison"
  }    
    }   
    }  
    }  
    }
  }
}
    }
    }
}

## Mean values data
#md_means <- md_pcrit_smr %>%
#  group_by(spps, temp) %>%
#  summarise(mean_pcrit = mean(pcrit.r), mean_smr = mean(smr.best))

#md_means$species <- 

#plot_means <- ggplot(data = md_means) +
#  geom_line(mapping = aes(x = temp, y = mean_pcrit, color = spps)) +
#  labs(
#    x = "Temperature (degrees C)",
#    y = "Pcrit (torr)",
#    colour = "Species"
#  ) +
#  theme_classic()


############################################
### MCMCglmm analysis of temperature effects on Pcrit
############################################

## Need to root tree somehow, maybe by choosing node - try using plot.phylo to pick node

knope_pruned_root <- midpoint.root(knope_pruned)
#knope_pruned_root <- root(knope_pruned, "Hemilepidotus_hemilepidotus")
inv_knope <- inverseA(knope_pruned_root, 
                      nodes = "TIPS",
                      scale = TRUE)

prior <- list(G = list(G1 = list(V = 1, nu = 0.02)),
              R = list(V = 1, nu = 0.02)
              )

## "ginverse" function not working - species missing levels
model <- MCMCglmm(pcrit.r ~ temp + spps + temp*spps,
                  random = ~ species, # Need to define this variable, see Vikram's email
                  ginverse = list(species = inv_knope$Ainv),
                  prior = prior, # Still need to make sure I am specifying my prior appropriately, this is Vikram's prior but looks pretty standard so far from Hadfeld's vignettes
                  nitt = 5000000, # Number of iterations to run
                  burnin = 1000, # sample and store every 1000th iteration
                  thin = 500, #  from help: "thinning interval"
                  data = md_pcrit_smr,
                  verbose = TRUE # Still not sure what this means; from help: "logical: if TRUE MH diagnostics are printed to screen"
                  )

summary(model)
plot(model)
## Extra, unused code, some potentially useful functions

# Make the tree rooted
# root(knope)
# Make the tree ultrametric
# chronopl(knope, )