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

#####################################################################
## Read sculpin phyologeny into R and prune for species in analysis
#####################################################################
ramon_tree <- read.nexus("T5000_FigTree_TimeScaledRoot1.nex")

## NOTE! This does not include Blepsias cirrhosus (not in Ramon tree)
keepers_ramon <- c("Oligocottus_maculosus",
                   "Clinocottus_globiceps",
                   "Artedius_harringtoni",
                   "Artedius_lateralis",
                   "Artedius_fenestaolis", ## Ramon tree has typo here, should be "fenestralis"
                   "Scorpaenichthys_marmoratus",
                   "Enophrys_bison",
                   "Hemilepidotus_hemilepidotus")

## Prune phylogeny for species in analysis
ramon_phy <- drop.tip(ramon_tree, setdiff(ramon_tree$tip.label, keepers_ramon))
## Fix "fenestralis" typo
ramon_phy$tip.label[ramon_phy$tip.label == "Artedius_fenestaolis"] <- "Artedius_fenestralis"
## Compute branch lengths, tree didn't have any
ramon_phy <- compute.brlen(ramon_phy, method = "Grafen")

plot.phylo(ramon_phy)

############################################
## Read sculpin pcrit and meta-data into R
############################################

md <- read.csv("SculpinPcritData_ComparativeAnalysisFormat.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
md <- as_tibble(md)

## Pcrit and temp data, NA's filtered out, only "trial.no == 1" inlcuded
md_pcrit_smr <- md %>%
  #select(species, spps, date, temp, fish.id, mass.g, smr.best, pcrit.r, trial.no) %>%
  filter(!is.na(smr.best)) %>%
  filter(!is.na(pcrit.r)) %>%
  filter(trial.no == 1)

## Averages for Pcrit and SMR per species per temperature, all 3 temps included
md_all_temps <- md_pcrit_smr %>%
  group_by(temp, species) %>%
  summarise(avg_pcrit = mean(pcrit.r), sd_pcrit = sd(pcrit.r), n_pcrit = length(pcrit.r),
            avg_smr = mean(smr.best), sd_smr = sd(smr.best), n_smr = length(smr.best)) %>%
  group_by(temp, species) %>%
  mutate(sem_pcrit = sd_pcrit/n_pcrit,
         sem_smr = sd_smr/n_smr)

## Open ct_max data and add to q10 data
ct_max_df <- read.csv("ct_max_data.csv", stringsAsFactors = FALSE,
                      strip.white = TRUE, na.strings = c("NA","."))
