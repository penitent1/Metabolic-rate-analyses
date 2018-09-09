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

#####################################################################
## Read sculpin phyologeny into R and prune for species in analysis
#####################################################################
mandic_2013_tree <- read.nexus(file.choose())
plot.phylo(mandic_2013_tree)
is.ultrametric(mandic_2013_tree) ## It is not :(
is.rooted(mandic_2013_tree) ## It's rooted!!

mandic_2013_tree_no_out <- drop.tip(mandic_2013_tree, "Satyrichthys_amiscus")

mandic_um <- chronopl(mandic_2013_tree_no_out, 1)
## Fix "fenestralis" typo
mandic_um$tip.label[mandic_um$tip.label == "EF521369.1_Hemilepidotus_hemilep"] <- "Hemilepidotus_hemilepidotus"
mandic_um$tip.label[mandic_um$tip.label == "Fluffy"] <- "Oligocottus_snyderi"
mandic_um$tip.label[mandic_um$tip.label == "Great_Sculpin_AB114909"] <- "Myoxocephalus_polyacanthocephalus"
mandic_um$tip.label[mandic_um$tip.label == "Mosshead"] <- "Clinocottus_globiceps"
mandic_um$tip.label[mandic_um$tip.label == "Pacific_Staghorn"] <- "Leptocottus_armatus"
mandic_um$tip.label[mandic_um$tip.label == "Padded"] <- "Artedius_fenestralis"
mandic_um$tip.label[mandic_um$tip.label == "Prickly"] <- "Cottus_asper"
mandic_um$tip.label[mandic_um$tip.label == "Tidepool"] <- "Oligocottus_maculosus"
mandic_um$tip.label[mandic_um$tip.label == "Scalyhead"] <- "Artedius_harringtoni"
mandic_um$tip.label[mandic_um$tip.label == "Shorthorn"] <- "Myoxocephalus_scorpius"
mandic_um$tip.label[mandic_um$tip.label == "Silverspotted"] <- "Blepsias_cirrhosus"
mandic_um$tip.label[mandic_um$tip.label == "Smoothhead"] <- "Artedius_lateralis"
mandic_um$tip.label[mandic_um$tip.label == "Cabezon"] <- "Scorpaenichthys_marmoratus"
mandic_um$tip.label[mandic_um$tip.label == "Buffalo"] <- "Enophrys_bison"
mandic_um$tip.label[mandic_um$tip.label == "cottus_bairdii"] <- "Cottus_bairdii"

plot.phylo(mandic_um) ## Ultrametric!

## Drop species not in my study:
keepers_mandic <- c("Oligocottus_maculosus",
                   "Clinocottus_globiceps",
                   "Artedius_harringtoni",
                   "Artedius_lateralis",
                   "Artedius_fenestralis", ## Ramon tree has typo here, should be "fenestralis"
                   "Scorpaenichthys_marmoratus",
                   "Enophrys_bison",
                   "Hemilepidotus_hemilepidotus",
                   "Blepsias_cirrhosus")
mandic_phy <- drop.tip(mandic_um, setdiff(mandic_um$tip.label, keepers_mandic))
plot.phylo(mandic_phy)

## Compute branch lengths, tree didn't have any
##ramon_phy <- compute.brlen(ramon_phy, method = "Grafen")

############################################
## Read sculpin pcrit and meta-data into R
############################################

md <- read.csv("SculpinPcritData_ComparativeAnalysisFormat_withPcritSlopes.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
md <- as_tibble(md)

## Pcrit and temp data, NA's filtered out, only "trial.no == 1" inlcuded
md_pcrit_smr <- md %>%
  #select(species, spps, date, temp, fish.id, mass.g, smr.best, pcrit.r, trial.no) %>%
  filter(!is.na(smr.best.ms)) %>%
  filter(!is.na(pcrit.r)) %>%
  filter(trial.no == 1) %>%
  mutate(smr.best.raw = smr.best.ms*mass.g) # Add NOT mass-specific mo2

## Averages for Pcrit and SMR per species per temperature, all 3 temps included
md_all_temps <- md_pcrit_smr %>%
  group_by(temp, species) %>%
  summarise(avg_pcrit = mean(pcrit.r), sd_pcrit = sd(pcrit.r), n_pcrit = length(pcrit.r),
            avg_smr_ms = mean(smr.best.ms), sd_smr_ms = sd(smr.best.ms), n_smr_ms = length(smr.best.ms),
            avg_smr_raw = mean(smr.best.raw), sd_smr_raw = sd(smr.best.raw), n_smr_raw = length(smr.best.raw)) %>%
  group_by(temp, species) %>%
  mutate(sem_pcrit = sd_pcrit/n_pcrit,
         sem_smr_ms = sd_smr_ms/n_smr_ms,
         sem_smr_raw = sd_smr_raw/n_smr_raw)

### Blocked ANOVA(?)
### ANCOVA
# Can also do a phylogenetic ANCOVA

m_ancova <- lm(pcrit.r ~ temp + spps, data = md_pcrit_smr)
m_ancova_int <- lm(pcrit.r ~ temp * spps, data = md_pcrit_smr) # interaction included

library(car)
Anova(m_ancova_int, type = "III")
summary(m_ancova_int)

library(broom) ## Great package that takes model fit coefficients and puts them in a data frame!

fitted_models <- md_pcrit_smr %>% 
  group_by(species) %>%
  do(model = lm(pcrit.r ~ temp, data = .))

model_df <- fitted_models %>% tidy(model)

library(tidyr)

model_df_spread <- model_df %>%
  dplyr::select(term, estimate) %>%
  spread(term, estimate)

cor(model_df_spread$`(Intercept)`, model_df_spread$temp) # Testing correlation between temp slope and intercept
## I'm not sure if the intercept of the lm: pcrit ~ temp is a fair representation of my hypothesis...

plot(model_df_spread$`(Intercept)`, model_df_spread$temp)
text(model_df_spread$`(Intercept)`, model_df_spread$temp, labels=model_df_spread$species, cex=0.7) #, cex= 0.7

## Visualize without interaction
summary(m_ancova)
visreg(m_ancova, xvar = "temp", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))
visreg(m_ancova, xvar = "temp", by = "spps", 
       points.par = list(cex = 1.1, col = "red"))
visreg(m_ancova, xvar = "temp", by = "spps", whitespace = 0.4, overlay = TRUE, 
       band = FALSE, points.par = list(cex = 1.1))

## Visualize WITH interaction
summary(m_ancova_int)
visreg(m_ancova_int, xvar = "temp", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))
visreg(m_ancova_int, xvar = "temp", by = "spps", 
       points.par = list(cex = 1.1, col = "red"))
visreg(m_ancova_int, xvar = "temp", by = "spps", whitespace = 0.4, overlay = TRUE, 
       band = FALSE, points.par = list(cex = 1.1))




## Open ct_max data and correct raw loe temperatures
ct_max_df <- read.csv("Sculpins_ctmax_data_10spps.csv", stringsAsFactors = FALSE,
                      strip.white = TRUE, na.strings = c("NA","."))
ct_max_df <- ct_max_df[ct_max_df$spps!="rhri",] ## Remove Grunt sculpins
## This is based on the test I did for the temp probe labelled "2"
ct_max_df$loe_temp_corrected <- (ct_max_df$loe_temp - 0.2909)/0.9857
#ct_max_df_phylo <- ct_max_df_no_rhri[,c(1,3)]

ct_max_df <- as_tibble(ct_max_df)
ct_max_avg_df <- ct_max_df %>%
  dplyr::select(loe_temp_corrected,spps_phylo) %>%
  group_by(spps_phylo) %>%
  summarize(avg_ctmax = mean(loe_temp_corrected),
            n_ctmax = length(spps_phylo),
            sd_ctmax = sd(loe_temp_corrected),
            sem_ctmax = sd_ctmax/sqrt(n_ctmax))
ct_max_avg_df

ctmax_plot <- ggplot(ct_max_avg_df, aes(x=spps_phylo, y=avg_ctmax)) +
  geom_errorbar(mapping = aes(x=spps_phylo, ymin=(avg_ctmax-sem_ctmax), ymax=(avg_ctmax+sem_ctmax)), width=0.75, size=1.25) +
  geom_line(size=1.25)  +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.title.y = (element_text(size = 20, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 09),
        axis.title.x = element_text(size = 20, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression(paste("CT"["max"]," (",degree,"C)")))
ctmax_plot

## Reduce sculpin phylo to species in ct_max dataset
mandic_phy_ctmax <- drop.tip(mandic_phy, "Blepsias_cirrhosus")
plot.phylo(mandic_phy_ctmax)

#############
# # # Data for ct max phenogram
#############
# Using ct_max_df, no q10 values or anything besides ct_max data
#rownames(ct_max_df_phylo)<-c() ## For some stupid reason I have to run this
ctmax_phylo_rn <- column_to_rownames(ct_max_avg_df, var = "spps_phylo")
ctmax_phylo_matrix <- as.matrix(ctmax_phylo_rn)[,1]
#############
# # # END Data for ct max phenogram
#############

## Phenogram of CTmax on phylogeny
pheno_ctmax <- phenogram(mandic_phy_ctmax, ctmax_phylo_matrix,
                         ylab=expression(paste("CTmax (",degree,C,")")),
                         #xlab=NULL,
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
