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
ct_max_df_no_rhri <- ct_max_df[ct_max_df$spps!="rhri",]
ct_max_df_phylo <- ct_max_df_no_rhri[,c(1,3)]

## Reduce sculpin phylo to species in ct_max dataset
ramon_phy_ctmax <- drop.tip(ramon_phy, "Hemilepidotus_hemilepidotus")
plot.phylo(ramon_phy_ctmax)

#############
# # # Data for ct max phenogram
#############
# Using ct_max_df, no q10 values or anything besides ct_max data
rownames(ct_max_df_phylo)<-c() ## For some stupid reason I have to run this
ctmax_phylo_rn <- column_to_rownames(ct_max_df_phylo, var = "species")
ctmax_phylo_matrix <- as.matrix(ctmax_phylo_rn)[,1]
#############
# # # END Data for ct max phenogram
#############

## Phenogram of CTmax on phylogeny
pheno_ctmax <- phenogram(ramon_phy_ctmax, ctmax_phylo_matrix,
                         ylab=expression(paste("CTmax (",degree,C,")")),
                         xlab=NULL,
                         xaxt='n')

### Now I want to get slopes of Pcrit vs Temp and put these in a matrix
### Then check to see how these look on a phenogram
# # # 

## Model: Pcrit ~ temp + spps + temp*spps

arfe <- md_all_temps[md_all_temps$species=="Artedius_fenestralis",]
arha <- md_all_temps[md_all_temps$species=="Artedius_harringtoni",]
arla <- md_all_temps[md_all_temps$species=="Artedius_lateralis",]
blci <- md_all_temps[md_all_temps$species=="Blepsias_cirrhosus",]
clgl <- md_all_temps[md_all_temps$species=="Clinocottus_globiceps",]
enbi <- md_all_temps[md_all_temps$species=="Enophrys_bison",]
hehe <- md_all_temps[md_all_temps$species=="Hemilepidotus_hemilepidotus",]
olma <- md_all_temps[md_all_temps$species=="Oligocottus_maculosus",]
scma <- md_all_temps[md_all_temps$species=="Scorpaenichthys_marmoratus",]

## Linear models of pcrit vs temperature for each species FOR PCRIT
arfe_lm <- lm(arfe$avg_pcrit~arfe$temp)
arha_lm <- lm(arha$avg_pcrit~arha$temp)
arla_lm <- lm(arla$avg_pcrit~arla$temp)
blci_lm <- lm(blci$avg_pcrit~blci$temp)
clgl_lm <- lm(clgl$avg_pcrit~clgl$temp)
enbi_lm <- lm(enbi$avg_pcrit~enbi$temp)
hehe_lm <- lm(hehe$avg_pcrit~hehe$temp)
olma_lm <- lm(olma$avg_pcrit~olma$temp)
scma_lm <- lm(scma$avg_pcrit~scma$temp)

coef(arfe_lm)[2] ## Extracts temp effect slope
coef(arha_lm)[2]
coef(arla_lm)[2]
coef(blci_lm)[2]
coef(clgl_lm)[2]
coef(enbi_lm)[2]
coef(hehe_lm)[2]
coef(olma_lm)[2]
coef(scma_lm)[2]

pcrit_temp_slopes <- c(coef(arfe_lm)[2],
                       coef(arha_lm)[2],
                       coef(arla_lm)[2],
                       coef(blci_lm)[2],
                       coef(clgl_lm)[2],
                       coef(enbi_lm)[2],
                       coef(hehe_lm)[2],
                       coef(olma_lm)[2],
                       coef(scma_lm)[2])

## Linear models of SMR vs temperature for each species
arfe_lm_smr <- lm(arfe$avg_smr~arfe$temp)
arha_lm_smr <- lm(arha$avg_smr~arha$temp)
arla_lm_smr <- lm(arla$avg_smr~arla$temp)
blci_lm_smr <- lm(blci$avg_smr~blci$temp)
clgl_lm_smr <- lm(clgl$avg_smr~clgl$temp)
enbi_lm_smr <- lm(enbi$avg_smr~enbi$temp)
hehe_lm_smr <- lm(hehe$avg_smr~hehe$temp)
olma_lm_smr <- lm(olma$avg_smr~olma$temp)
scma_lm_smr <- lm(scma$avg_smr~scma$temp)

coef(arfe_lm_smr)[2] ## Extracts temp effect slope
coef(arha_lm_smr)[2]
coef(arla_lm_smr)[2]
coef(blci_lm_smr)[2]
coef(clgl_lm_smr)[2]
coef(enbi_lm_smr)[2]
coef(hehe_lm_smr)[2]
coef(olma_lm_smr)[2]
coef(scma_lm_smr)[2]

smr_temp_slopes <- c(coef(arfe_lm_smr)[2],
                     coef(arha_lm_smr)[2],
                     coef(arla_lm_smr)[2],
                     coef(blci_lm_smr)[2],
                     coef(clgl_lm_smr)[2],
                     coef(enbi_lm_smr)[2],
                     coef(hehe_lm_smr)[2],
                     coef(olma_lm_smr)[2],
                     coef(scma_lm_smr)[2])

spps <- c("arfe",
          "arha",
          "arla",
          "blci",
          "clgl",
          "enbi",
          "hehe",
          "olma",
          "scma")

species <- c("Artedius_fenestralis",
             "Artedius_harringtoni",
             "Artedius_lateralis",
             "Blepsias_cirrhosus",
             "Clinocottus_globiceps",
             "Enophrys_bison",
             "Hemilepidotus_hemilepidotus",
             "Oligocottus_maculosus",
             "Scorpaenichthys_marmoratus")

md_phen_df <- data.frame(species, spps, pcrit_temp_slopes, smr_temp_slopes)
md_phen_df <- md_phen_df[md_phen_df$spps!= "blci",]

#############
# # # Data for pcrit-smr phenograms
#############
# Using ct_max_df, no q10 values or anything besides ct_max data
md_phen_df_pcrit <- md_phen_df[,c(1,3)]
rownames(md_phen_df_pcrit)<-c() ## For some stupid reason I have to run this
pcrit_slopes_phen_rn <- column_to_rownames(md_phen_df_pcrit, var = "species")
pcrit_slopes_phen_rn <- pcrit_slopes_phen_rn[match(ramon_phy$tip.label,rownames(pcrit_slopes_phen_rn)),]
pcrit_slopes_phen_matrix <- as.matrix(pcrit_slopes_phen_rn)[,1]

md_phen_df_smr <- md_phen_df[,c(1,4)]
rownames(md_phen_df_smr)<-c() ## For some stupid reason I have to run this
smr_slopes_phen_rn <- column_to_rownames(md_phen_df_smr, var = "species")
smr_slopes_phen_matrix <- as.matrix(smr_slopes_phen_rn)[,1]
#############
# # # END Data for Pcrit/SMR phenograms
#############

## Phenogram of Pcrit-temp slopes on phylogeny
pheno_pcrit <- phenogram(ramon_phy, pcrit_slopes_phen_matrix,
                         ylab=expression(paste("Pcrit-temp slope (torr per ",degree,C,")")),
                         xlab=NULL,
                         xaxt='n')

## Phenogram of SMR-temp slopes on phylogeny
pheno_smr <- phenogram(ramon_phy, smr_slopes_phen_matrix,
                         ylab=expression(paste("SMR-temp slope (",delta,"MO2 per ",degree,C,")")),
                         xlab=NULL,
                         xaxt='n')


###########################################
# Colored by species ID: Pcrit vs temp plot

species_names <- c(
  "Artedius fenestralis",
  "Artedius harringtoni",
  "Artedius lateralis",
  "Blepsias cirrhosus",
  "Clinocottus globiceps",
  "Enophrys bison",
  "Hemilepidotus hemilepidotus",
  "Oligocottus maculosus",
  "Scorpaenichthys marmoratus"
)

md_all_temps$spps_names <- species_names

pd <- position_dodge(0.2)

pcrit_vs_temp_plot <- ggplot(md_all_temps, aes(x=temp, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  geom_line(position=pd, size=1.25)  +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Species"))
pcrit_vs_temp_plot
