library(tidyverse)
library(ggthemes)
library(mclust)
library(shape)
library(StreamMetabolism)
library(fishMO2)
#library(cowplot) # Decided to manually adjust plots - too confusing with cowplot
library(gridExtra)
library(ape)
library(caper)
library(geiger)
library(phytools)
library(MCMCglmm)
library(visreg)
library(nlme)
library(broom)
library(car)

## *******************************
##
##    Importing data
##
## *******************************

## ******************
## Pcrit and SMR data
## ******************
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses")
pcrit_smr_data_summary <- read.csv("SculpinPcritData_ComparativeAnalysisFormat_withPcritSlopes.csv", stringsAsFactors = FALSE,
                                   strip.white = TRUE, na.strings = c("NA","."))
pcrit_smr_data_summary <- as_tibble(pcrit_smr_data_summary)
pcrit_smr_data_summary <- pcrit_smr_data_summary %>%
  mutate(smr.raw = smr.best.ms*mass.g) %>%
  filter(!is.na(pcrit.r), !is.na(smr.best.ms), trial.no==1)

## *****************************
## Open ct_max data 
## Correct raw loe temperatures
##
## *****************************

#NOTE! I don't have CTmax data for BLCI
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/CT max data")
ct_max_df <- read_csv("CT_max_sculpins_scma_enbi_BOTH_trials_1_2_all_data_included.csv") %>% # Use data file that includes trials 1 and 2 for scma and enbi
  filter(species != "rhri") %>% ## Remove Grunt sculpins
  mutate(loe_temp_corrected = (loe_temp - 0.2909)/0.9857) ## Correcting for calibration; based on test of temp probe labelled "2"

ct_max_data_spps_summary <- ct_max_df %>%
  group_by(species) %>%
  dplyr::summarise(ct_max_avg = mean(loe_temp_corrected))

## *****************************************************************
##
## Read sculpin phyologeny into R and prune for species in analysis
##
## *****************************************************************
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Phylogeny files")
mandic_2013_tree <- read.nexus("Oct62009_Likelihood w branchlengths.tre")
plot.phylo(mandic_2013_tree)
is.ultrametric(mandic_2013_tree) ## It is not :(
is.rooted(mandic_2013_tree) ## It's rooted!!

mandic_2013_tree_no_out <- drop.tip(mandic_2013_tree, "Satyrichthys_amiscus")

mandic_um <- chronopl(mandic_2013_tree_no_out, 1)
#mandic_um <- chronos(mandic_2013_tree_no_out, 1)
## Fix tip names so they actually give species names
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
mandic_ctmax_phy <- drop.tip(mandic_phy, "Blepsias_cirrhosus") # No CTmax for BLCI
plot.phylo(mandic_ctmax_phy)

plot(mandic_phy, cex = 2, no.margin = TRUE, x.lim = 2, label.offset = 0.4, edge.width = 3.5,
     use.edge.length = TRUE)

# 
# printed plot size: W = 1132 H = 932

#######################
#######################                        

pcrit_smr_data_summary
ct_max_df
ct_max_data_spps_summary
plot.phylo(mandic_phy)
plot.phylo(mandic_ctmax_phy)

## *******************************************************************
##
##  Body mass VS smr/Pcrit:  Linear models and AOVs
##
## *******************************************************************

lm_scaling_smr_pcrit <- pcrit_smr_data_summary %>%
  group_by(species) %>%
  nest() %>%
  mutate(spps_scaling_mod_smr = data %>% purrr::map(~ lm(log(smr.raw)~log(mass.g), data=.)),
         spps_scaling_mod_pcrit = data %>% purrr::map(~ lm(log(pcrit.r)~log(mass.g), data=.)),
         aov_scaling_mod_smr = spps_scaling_mod_smr %>% purrr::map(~ Anova(., type = "III")),
         aov_scaling_mod_pcrit = spps_scaling_mod_pcrit %>% purrr::map(~ Anova(., type = "III")),
         median_mass = data %>% purrr::map_dbl(~ median(.$mass.g)),
         mean_mass = data %>% purrr::map_dbl(~ mean(.$mass.g)),
         range_mass = data %>% purrr::map_dbl(~ (max(.$mass.g) - min(.$mass.g))),
         range_fold_mass = data %>% purrr::map_dbl(~ (max(.$mass.g) / min(.$mass.g))),
         tidy_smr_aov = aov_scaling_mod_smr %>% purrr::map(broom::tidy),
         tidy_smr_lm = spps_scaling_mod_smr %>% purrr::map(broom::tidy),
         tidy_pcrit_aov = aov_scaling_mod_pcrit %>% purrr::map(broom::tidy),
         tidy_pcrit_lm = spps_scaling_mod_pcrit %>% purrr::map(broom::tidy),
         f_val_smr_aov = tidy_smr_aov %>% purrr::map_dbl(c(4,2)),
         f_df_num_smr_aov = tidy_smr_aov %>% purrr::map_dbl(c(3,2)),
         f_df_den_smr_aov = tidy_smr_aov %>% purrr::map_dbl(c(3,3)),
         f_val_pcrit_aov = tidy_pcrit_aov %>% purrr::map_dbl(c(4,2)),
         f_df_num_pcrit_aov = tidy_pcrit_aov %>% purrr::map_dbl(c(3,2)),
         f_df_den_pcrit_aov = tidy_pcrit_aov %>% purrr::map_dbl(c(3,3)),
         p_val_smr_aov = tidy_smr_aov %>% purrr::map_dbl(c(5,2)),
         p_val_pcrit_aov = tidy_pcrit_aov %>% purrr::map_dbl(c(5,2)),
         slope_smr = tidy_smr_lm %>% purrr::map_dbl(c(2,2)),
         slope_pcrit = tidy_pcrit_lm %>% purrr::map_dbl(c(2,2)))

## Adjust smr for body mass effects on SMR:

mass_corr_smr_pcrit_data <- 
  lm_scaling_smr_pcrit %>%
  dplyr::select(species, data, median_mass, mean_mass, range_mass, range_fold_mass, 
                p_val_smr_aov, slope_smr, slope_pcrit) %>%
  unnest() %>%
  mutate(smr.mass.corr = if_else(p_val_smr_aov < 0.05,
                                 smr.raw*exp(slope_smr*log(mean_mass/mass.g)),
                                 smr.raw),
         smr.mass.corr.ms = if_else(p_val_smr_aov < 0.05,
                                    smr.mass.corr/mean_mass,
                                    smr.raw/mass.g))

## In Deutsch et al., body size correction of "HYPOXIA TOLERANCE" eg Pcrit
## was only done if body mass range was greater than 3 fold

lm_scaling_smr_pcrit %>% 
  mutate(range_bigger_3fold = if_else(range_fold_mass > 3, "yes", "no")) %>%
  dplyr::select(species, p_val_smr_aov, p_val_pcrit_aov, range_mass, range_fold_mass)

lm_scaling_smr_pcrit[,c("species","p_val_smr_aov","p_val_pcrit_aov")] %>%
  mutate(sig_smr = if_else(p_val_smr_aov<0.05, "yes", "no"),
         sig_pcrit = if_else(p_val_pcrit_aov<0.05, "yes", "no"))

## Artedius harrringtoni is the only species with a significant effect
## of body mass on Pcrit, but the fold-range of body masses is 2.39.

summary(lm_scaling_smr_pcrit$spps_scaling_mod_pcrit[[3]])

lm_scaling_smr_pcrit[3,]

## Looks like there's only 3 data points for ARHA at 16 degrees - need to double check this
lm_scaling_smr_pcrit[3,2][[1]]

# Plot of Pcrit vs mass for ARHA - consider log-log plot?
pcrit_smr_data_summary %>%
  filter(species == "Artedius_harringtoni") %>%
  ggplot(aes(x = mass.g, y = pcrit.r)) +
  geom_point(size = 5) +
  stat_smooth(method = "lm")

# No real difference with log-log plot so probably stick with linear plot
pcrit_smr_data_summary %>%
  filter(species == "Artedius_harringtoni") %>%
  ggplot(aes(x = log(mass.g), y = log(pcrit.r))) +
  geom_point(size = 5) +
  stat_smooth(method = "lm") +
  scale_x_continuous(limits = c(0,3.25),
                     breaks = seq(0,3.25,0.25)) +
  scale_y_continuous(limits = c(0,4.5),
                     breaks = seq(0,4.5,0.5))

## ******************************************************************
##
## Figure 1:
## Raw Pcrit vs temperature, every species in it's own plot
##
## ******************************************************************

mass_corr_smr_pcrit_data %>%
  group_by(species, temp) %>%
  mutate(mean_pcrit = mean(pcrit.r),
         mean_smr = mean(smr.mass.corr.ms),
         species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  ggplot(aes(x = temp, y = pcrit.r)) +
  geom_jitter(width = 0.2, size = 3) +
  geom_line(aes(x = temp, y = mean_pcrit), size = 1.5) +
  scale_x_continuous(name = expression(paste("Temperature (",degree,C,")")),
                     limits = c(9,23)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,100)) +
  facet_wrap("species_plotting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text.x = element_text(size = rel(2.25), colour = "black"),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        strip.text = element_text(face = "bold", size = rel(1.25)),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)))

## ******************************************************************
##
## Figure 1.5:
## Mass-adjusted MO2 vs temperature, every species in it's own plot
##
## ******************************************************************

mass_corr_smr_pcrit_data %>%
  group_by(species, temp) %>%
  mutate(mean_pcrit = mean(pcrit.r),
         mean_smr = mean(smr.mass.corr.ms),
         species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  ggplot(aes(x = temp, y = smr.mass.corr.ms)) +
  geom_jitter(width = 0.2, size = 3) +
  geom_line(aes(x = temp, y = mean_smr), size = 1.5) +
  scale_x_continuous(name = expression(paste("Temperature (",degree,C,")")),
                     limits = c(9,23)) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2][",standard*"]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,6),
                     breaks = seq(0,6,2)) +
  facet_wrap("species_plotting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text.x = element_text(size = rel(2.25), colour = "black"),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        strip.text = element_text(face = "bold", size = rel(1.25)),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)))

## Arrhenius plot: MO2 vs temperature
mass_corr_smr_pcrit_data %>%
  group_by(species, temp) %>%
  mutate(mean_pcrit = mean(pcrit.r),
         mean_smr = mean(smr.mass.corr.ms),
         species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus"),
         temp_k = temp + 273.15) %>%
  ggplot(aes(x = 1/temp_k, y = log(smr.mass.corr.ms))) +
  geom_point() +
  #geom_jitter(width = 0.2, size = 3) +
  #geom_line(aes(x = temp, y = mean_smr), size = 1.5) +
  stat_smooth(method = "lm") +
  #scale_x_continuous(name = expression(paste("Temperature (",degree,C,")")),
  #                   limits = c(9,23)) +
  #scale_y_continuous(name = expression(paste(dot(M),"o"[2][",standard*"]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
  #                   limits = c(0,6),
  #                   breaks = seq(0,6,2)) +
  facet_wrap("species_plotting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text.x = element_text(size = rel(2.25), colour = "black"),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        strip.text = element_text(face = "bold", size = rel(1.25)),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)))

## Not facetted
mass_corr_smr_pcrit_data %>%
  group_by(species, temp) %>%
  mutate(mean_pcrit = mean(pcrit.r),
         mean_smr = mean(smr.mass.corr.ms),
         species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus"),
         tpo = case_when(species == "Oligocottus_maculosus" ~ "Present",
                         species == "Clinocottus_globiceps" ~ "Present",
                         species == "Artedius_harringtoni" ~ "Present",
                         species == "Artedius_lateralis" ~ "Present",
                         species == "Artedius_fenestralis" ~ "Present",
                         species == "Blepsias_cirrhosus" ~ "Absent",
                         species == "Enophrys_bison" ~ "Absent",
                         species == "Hemilepidotus_hemilepidotus" ~ "Absent",
                         species == "Scorpaenichthys_marmoratus" ~ "Absent")) %>%
  ggplot() +
  geom_jitter(aes(x = temp, y = smr.mass.corr.ms, colour = tpo, shape = species_plotting), width = 0.2, size = 3, alpha = 0.75) +
  geom_line(aes(x = temp, y = mean_smr, group = species_plotting, colour = tpo), size = 1.5) +
  scale_x_continuous(name = expression(paste("Temperature (",degree,C,")")),
                     limits = c(9,23)) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2][",standard*"]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,6),
                     breaks = seq(0,6,2)) +
  #facet_wrap("species_plotting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text.x = element_text(size = rel(2.25), colour = "black"),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        #strip.text = element_text(face = "bold", size = rel(1.25)),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)))

## Checking sample size per spps per temp
smr_pcrit_temp_data_check <- 
  mass_corr_smr_pcrit_data %>%
  group_by(species, temp) %>%
  nest()

#################################################################################
##
## Figure(# ?) Comparing Q10 values between species based on tidepool occupancy
##
#################################################################################

## First plot out Arrhenius plot of MEAN mo2 at each temp to see if Q10 appropriate
mass_corr_smr_pcrit_data %>%
  group_by(species,temp) %>%
  dplyr::select(species, temp, smr.mass.corr.ms) %>%
  summarise(mo2_mean = mean(smr.mass.corr.ms),
            mo2_sd = sd(smr.mass.corr.ms)) %>%
  ggplot(aes(x = 1/((temp+273.15)*1000), y = log(mo2_mean))) +
  geom_point() +
  stat_smooth(method = "lm") +
  facet_wrap("species")
  
mass_corr_smr_pcrit_data %>%
  group_by(species,temp) %>%
  dplyr::select(species, temp, smr.mass.corr.ms) %>%
  summarise(mo2_mean = mean(smr.mass.corr.ms)) %>%
  spread(temp, mo2_mean) %>%
  mutate(q10_mo2 = (`20`/`12`)^(10/(20-12)),
         tpo = case_when(species == "Artedius_fenestralis"~"yes",
                         species == "Artedius_harringtoni"~"yes",
                         species == "Artedius_lateralis"~"yes",
                         species == "Blepsias_cirrhosus"~"no",
                         species == "Clinocottus_globiceps"~"yes",
                         species == "Enophrys_bison"~"no",
                         species == "Hemilepidotus_hemilepidotus"~"no",
                         species == "Oligocottus_maculosus"~"yes",
                         species == "Scorpaenichthys_marmoratus"~"no")) %>%
  ggplot(aes(x = tpo, y = q10_mo2)) +
  geom_jitter(width = 0.05, size = 4) +
  scale_x_discrete(name = expression(paste("Tidepool occupancy"))) +
  scale_y_continuous(name = expression(paste("Q"[10]," of ",dot(M)["o"][2][",routine"])),
                     limits = c(0,4),
                     breaks = seq(0,4,0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text.x = element_text(size = rel(2.25), colour = "black"),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)))


## ******************************************************************
##
## Figure 2:
## Beta-Pcrit-low_temps vs Pcrit at 12 C
##
## ******************************************************************

## Get linear models of Pcrit vs temperature:
## 1) lm for between 12-16 degrees
## 2) lm for between 16-20 degrees

lm_pcrit_smr_temp <- mass_corr_smr_pcrit_data %>%
  group_by(species) %>%
  nest() %>%
  mutate(data_low_temps = data %>% purrr::map(~ filter(., temp != 20)),
         data_high_temps = data %>% purrr::map(~ filter(., temp !=12)),
         data_mean_pcrit_mo2_12 = data %>% purrr::map(~ filter(., temp == 12)),
         lm_pcrit_low_temps = data_low_temps %>% purrr::map(~ lm(pcrit.r ~ temp, data=.)),
         lm_pcrit_high_temps = data_high_temps %>% purrr::map(~ lm(pcrit.r ~ temp, data=.)),
         tidy_pcrit_low_temps = lm_pcrit_low_temps %>% purrr::map(broom::tidy),
         tidy_pcrit_high_temps = lm_pcrit_high_temps %>% purrr::map(broom::tidy),
         slope_pcrit_low_temps = tidy_pcrit_low_temps %>% purrr::map_dbl(c(2,2)),
         slope_pcrit_high_temps = tidy_pcrit_high_temps %>% purrr::map_dbl(c(2,2)),
         mean_pcrit_12 = data_mean_pcrit_mo2_12 %>% purrr::map_dbl(~ mean(.$pcrit.r)),
         lm_mo2_low_temps = data_low_temps %>% purrr::map(~ lm(smr.mass.corr.ms ~ temp, data=.)),
         lm_mo2_high_temps = data_high_temps %>% purrr::map(~ lm(smr.mass.corr.ms ~ temp, data=.)),
         tidy_mo2_low_temps = lm_mo2_low_temps %>% purrr::map(broom::tidy),
         tidy_mo2_high_temps = lm_mo2_high_temps %>% purrr::map(broom::tidy),
         slope_mo2_low_temps = tidy_mo2_low_temps %>% purrr::map_dbl(c(2,2)),
         slope_mo2_high_temps = tidy_mo2_high_temps %>% purrr::map_dbl(c(2,2)),
         mean_mo2_12 = data_mean_pcrit_mo2_12 %>% purrr::map_dbl(~ mean(.$smr.mass.corr.ms)))

## Beta pcrit: 12-16 degrees ~ Pcrit at 12 degrees
beta_pcrit_low_temp_plot <- 
  lm_pcrit_smr_temp[,c(1,10,11,12)] %>%
  ggplot(aes(x = mean_pcrit_12, y = slope_pcrit_low_temps)) +
  geom_point(size = 5) +
  #stat_smooth(method = "lm", size = 2) +
  scale_x_continuous(name = expression(paste("P"["crit"]," at 12",degree,C," (Torr)")),
                     limits = c(10, 60),
                     breaks = seq(10,60,10)) +
  scale_y_continuous(name = expression(paste(beta["P"]["crit"][" 12-16"][degree][C]," (Torr  ",degree,C^-1,")")),
                     limits = c(-2,12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)),
        plot.margin = unit(c(0.1,0.2,0.2,0.2), "in"))


## Beta pcrit: 16-20 degrees ~ Pcrit at 12 degrees
beta_pcrit_high_temp_plot <- 
  lm_pcrit_smr_temp[,c(1,10,11,12)] %>%
  ggplot(aes(x = mean_pcrit_12, y = slope_pcrit_high_temps)) +
  geom_point(size = 5) +
  #stat_smooth(method = "lm", size = 2) +
  scale_x_continuous(name = expression(paste("P"["crit"]," at 12",degree,C," (Torr)")),
                     limits = c(10, 60),
                     breaks = seq(10,60,10)) +
  scale_y_continuous(name = expression(paste(beta["P"]["crit "]["16-20"][degree][C]," (Torr  ",degree,C^-1,")")),
                     limits = (c(-2, 12))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)),
        plot.margin = unit(c(0.1,0.2,0.2,0.2), "in"))


grid.arrange(beta_pcrit_low_temp_plot, beta_pcrit_high_temp_plot, ncol=2)
# PNG exported: W = 1132 H = 654

## *******************************************
##
##  Plotting Pcrit against SMR
##
## *******************************************

# Version 1: all species in one plot with species coded by colour
mass_corr_smr_pcrit_data %>%
  mutate(species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  ggplot() +
  #stat_smooth(aes(x = smr.mass.corr.ms, y = pcrit.r), method = "lm", color = "black") +
  stat_smooth(aes(x = smr.mass.corr.ms, y = pcrit.r, colour = species_plotting), 
              method = "lm", se = FALSE, size = rel(2)) + #, colour = "black"
  scale_x_continuous(name = expression(paste(dot(M),"o"[2][",standard*"]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,6),
                     breaks = seq(0,6,1)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,100)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5))) +
  labs(colour = "Species")

###################################

# Version 1.5: Po2 on X, Mo2 on Y

###################################

## Colored by Tidepool occupancy
mass_corr_smr_pcrit_data %>%
  mutate(species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  group_by(species_plotting,temp) %>%
  dplyr::select(species_plotting,temp,pcrit.r,smr.mass.corr.ms) %>%
  summarise(mean_mo2 = mean(smr.mass.corr.ms),
            sd_mo2 = sd(smr.mass.corr.ms),
            n_mo2 = length(smr.mass.corr.ms),
            sem_mo2 = sd_mo2/sqrt(n_mo2),
            mean_pcrit = mean(pcrit.r),
            sd_pcrit = sd(pcrit.r),
            n_pcrit = length(pcrit.r),
            sem_pcrit = sd_pcrit/sqrt(n_pcrit)) %>%
  mutate(tpo = case_when(species_plotting == "Oligocottus maculosus" ~ "yes",
                         species_plotting == "Clinocottus globiceps" ~ "yes",
                         species_plotting == "Artedius harringtoni" ~ "yes",
                         species_plotting == "Artedius lateralis" ~ "yes",
                         species_plotting == "Artedius fenestralis" ~ "yes",
                         species_plotting == "Blepsias cirrhosus" ~ "no",
                         species_plotting == "Enophrys bison" ~ "no",
                         species_plotting == "Hemilepidotus hemilepidotus" ~ "no",
                         species_plotting == "Scorpaenichthys marmoratus" ~ "no")) %>%
  ggplot() +
  geom_point(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting, colour = tpo), size = 5) +
  geom_errorbar(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting, colour = tpo,
                    ymin = mean_mo2-sem_mo2, ymax = mean_mo2+sem_mo2), size = 1.5) +
  geom_errorbarh(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting, colour = tpo,
                     xmin = mean_pcrit-sem_pcrit, xmax = mean_pcrit+sem_pcrit), size = 1.5) +
  geom_line(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting, colour = tpo), size = 1.25) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2][",routine"]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,6),
                     breaks = seq(0,6,1)) +
  scale_x_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,100)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5))) +
  labs(colour = "Tidepool occupancy")

## Colored by species
mass_corr_smr_pcrit_data %>%
  mutate(species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  group_by(species_plotting,temp) %>%
  dplyr::select(species_plotting,temp,pcrit.r,smr.mass.corr.ms) %>%
  summarise(mean_mo2 = mean(smr.mass.corr.ms),
            sd_mo2 = sd(smr.mass.corr.ms),
            n_mo2 = length(smr.mass.corr.ms),
            sem_mo2 = sd_mo2/sqrt(n_mo2),
            mean_pcrit = mean(pcrit.r),
            sd_pcrit = sd(pcrit.r),
            n_pcrit = length(pcrit.r),
            sem_pcrit = sd_pcrit/sqrt(n_pcrit)) %>%
  ggplot() +
  geom_point(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting, colour = species_plotting), size = 5) +
  geom_errorbar(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting, colour = species_plotting,
                    ymin = mean_mo2-sem_mo2, ymax = mean_mo2+sem_mo2), size = 1.5) +
  geom_errorbarh(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting, colour = species_plotting,
                     xmin = mean_pcrit-sem_pcrit, xmax = mean_pcrit+sem_pcrit), size = 1.5) +
  geom_line(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting, colour = species_plotting), size = 1.25) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2][",routine"]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,5),
                     breaks = seq(0,5,1)) +
  scale_x_continuous(name = expression(paste("Po"["2"]," (Torr)")),
                     limits = c(0,100)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5))) +
  labs(colour = "Species")

## Facetted on species
## MO2 ~ Pcrit AS PO2! has color coded data for EACH INDIVIDUAL
lol_curve_mean_data <- mass_corr_smr_pcrit_data %>%
  mutate(species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  group_by(species_plotting,temp) %>%
  dplyr::select(species_plotting,temp,pcrit.r,smr.mass.corr.ms) %>%
  summarise(mean_mo2 = mean(smr.mass.corr.ms),
            sd_mo2 = sd(smr.mass.corr.ms),
            n_mo2 = length(smr.mass.corr.ms),
            sem_mo2 = sd_mo2/sqrt(n_mo2),
            mean_pcrit = mean(pcrit.r),
            sd_pcrit = sd(pcrit.r),
            n_pcrit = length(pcrit.r),
            sem_pcrit = sd_pcrit/sqrt(n_pcrit))

lol_curve_all_data <- 
  mass_corr_smr_pcrit_data %>%
  mutate(species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  dplyr::select(species_plotting, temp, pcrit.r, smr.mass.corr.ms)


  ggplot(lol_curve_mean_data) +
  facet_wrap("species_plotting") +
  geom_point(data = lol_curve_all_data, aes(x = pcrit.r, 
                                              y = smr.mass.corr.ms, 
                                              group = species_plotting, 
                                              colour = temp),
               size = 2) +
  scale_colour_gradient(low = "blue", high = "red", breaks = c(12,16,20)) +
  labs(colour = expression(paste("Temperature (",degree,C,")"))) +
  geom_point(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting), 
             size = 2, stroke = 1.5) +
  geom_errorbar(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting,
                    ymin = mean_mo2-sem_mo2, ymax = mean_mo2+sem_mo2),
                size = 1.5) +
  geom_errorbarh(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting,
                     xmin = mean_pcrit-sem_pcrit, xmax = mean_pcrit+sem_pcrit),
                 size = 1.5) +
  geom_line(aes(x = mean_pcrit, y = mean_mo2, group = species_plotting),
            size = 1.25) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2][",routine"]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,5),
                     breaks = seq(0,5,1)) +
  scale_x_continuous(name = expression(paste("Po"["2"]," (Torr)")),
                     limits = c(0,100)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"),
        strip.text = element_text(face = "bold", size = rel(1.25)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)))


###########################################
# Version 2: each species in their own plot
mass_corr_smr_pcrit_data %>%
  mutate(species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  ggplot(aes(x = smr.mass.corr.ms, y = pcrit.r)) +
  stat_smooth(method = "lm", size = rel(2), colour = "black") +
  geom_point(size = 3) +
  scale_x_continuous(name = expression(paste(dot(M),"o"[2][",standard*"]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,6)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,100)) +
  facet_wrap("species_plotting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text.x = element_text(size = rel(2.25), colour = "black"),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        strip.text = element_text(face = "bold", size = rel(1.25)),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)))


## Linear model: Pcrit ~ SMR + species + SMR*species
pcrit_smr_smr_spps_ancova <- 
  mass_corr_smr_pcrit_data %>%
  lm(pcrit.r ~ smr.mass.corr.ms + species + smr.mass.corr.ms*species, data = .) %>%
  Anova(., type = "III")

# Significant effects of SMR, species, and SMR*species on pcrit
# THIS MEANS!
# SMR is a significant predictor of SMR, but the sensitivity of pcrit
# to variation in SMR varies between species.

pcrit_vs_smr_aov_per_species_data <- 
  mass_corr_smr_pcrit_data %>%
  group_by(species) %>%
  nest() %>%
  mutate(pcrit_vs_smr_lm = data %>% purrr::map(~ lm(pcrit.r ~ smr.mass.corr.ms, 
                                                    data = .)),
         aov_pcrit_vs_smr_lm = pcrit_vs_smr_lm %>% purrr::map(~ Anova(., type = "III")),
         tidy_pcrit_vs_smr_lm = pcrit_vs_smr_lm %>% purrr::map(broom::tidy),
         tidy_pcrit_vs_smr_aov = aov_pcrit_vs_smr_lm %>% purrr::map(broom::tidy),
         slope_pcrit_smr = tidy_pcrit_vs_smr_lm %>% purrr::map_dbl(c(2,2)),
         p_val_pcrit_smr = tidy_pcrit_vs_smr_aov %>% purrr::map_dbl(c(5,2)),
         f_val = tidy_pcrit_vs_smr_aov %>% purrr::map_dbl(c(4,2)),
         f_df_num = tidy_pcrit_vs_smr_aov %>% purrr::map_dbl(c(3,2)),
         f_df_den = tidy_pcrit_vs_smr_aov %>% purrr::map_dbl(c(3,3)))

pcrit_vs_smr_aov_per_species_data[,c("species","slope_pcrit_smr","f_val","f_df_num","f_df_den","p_val_pcrit_smr")]
  
## Artedius harringtoni is the only species 
## which does not have a significant relationship: Pcrit ~ Mo2 (all temp included):
pcrit_vs_smr_aov_per_species_data$tidy_pcrit_vs_smr_lm[[3]]



## *******************************************
##
##  Plotting Delta(SMR:20-12) against species
##
## *******************************************

delta_smr_temp_data <- 
  mass_corr_smr_pcrit_data %>%
  mutate(species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Blepsias_cirrhosus" ~ "Blepsias cirrhosus",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus")) %>%
  group_by(species, temp) %>%
  nest() %>%
  mutate(avg_smr = data %>% purrr::map_dbl(~ mean(.$smr.mass.corr.ms)),
         avg_pcrit = data %>% purrr::map_dbl(~ mean(.$pcrit.r))) %>%
  filter(temp != 16) %>%
  dplyr::select(-data) %>%
  group_by(species) %>%
  nest() %>%
  mutate(delta_12_data = data %>% purrr::map(~ filter(., temp == 12)),
         delta_20_data = data %>% purrr::map(~ filter(., temp == 20)),
         delta_smr_12_20 = purrr::map2_dbl(delta_20_data, 
                                           delta_12_data, 
                                           ~ .x$avg_smr -.y$avg_smr),
         delta_pcrit_12_20 = purrr::map2_dbl(delta_20_data,
                                             delta_12_data,
                                             ~ .x$avg_pcrit - .y$avg_pcrit))

lm_pcrit_vs_smr_data <- 
  mass_corr_smr_pcrit_data %>%
  group_by(species) %>%
  nest() %>%
  mutate(lm_pcrit_smr = data %>% 
           purrr::map(~ lm(pcrit.r ~ smr.mass.corr.ms, data = .)),
         tidy_lm_pcrit_smr = lm_pcrit_smr %>% purrr::map(broom::tidy),
         slope_pcrit_smr = tidy_lm_pcrit_smr %>% purrr::map_dbl(c(2,2)))

lm_pcrit_vs_smr_delta_smr_data <- tibble(species = lm_pcrit_vs_smr_data$species,
                                         delta_smr_12_20 = delta_smr_temp_data$delta_smr_12_20,
                                         delta_pcrit_12_20 = delta_smr_temp_data$delta_pcrit_12_20,
                                         slope_pcrit_smr = lm_pcrit_vs_smr_data$slope_pcrit_smr)

## Plot relationship between how much SMR increases and the SLOPE of Pcrit ~ SMR
ggplot(lm_pcrit_vs_smr_delta_smr_data, aes(x = delta_smr_12_20, y = slope_pcrit_smr)) +
  geom_point(size = 5) +
  stat_smooth(method = "lm")

## Plot relationship between how much SMR increases VS how much Pcrit increases
ggplot(lm_pcrit_vs_smr_delta_smr_data, aes(x = delta_smr_12_20, y = delta_pcrit_12_20)) +
  geom_point(size = 5) +
  #stat_smooth(method = "lm", size = 2) +
  scale_x_continuous(name = expression(paste(Delta,dot(M),"o"["2,standard* "]["12-20"][degree]["C"]," (",mu,"mol O"[2], " g"^-1," h"^-1,")")),
                     limits = c(0,4)) +
  scale_y_continuous(name = expression(paste(Delta,"P"["crit 12-20"][degree]["C"]," (Torr)")),
                     limits = c(0,50)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        plot.margin = unit(c(0.1,0.2,0.2,0.2), "in"))


# Anova on Delta Pcrit ~ Delta SMR across the tested temperature range (12-20 degrees)
#lm_pcrit_vs_smr_delta_smr_data %>%
#  lm(delta_pcrit_12_20 ~ delta_smr_12_20, data = .) %>%
#  Anova(., type = "III")

## >>>>> SEE analysis section at end of script for PGLS (caper) statistical analysis

## *******************************************
##
##  CTmax vs Beta Pcrit
##
## *******************************************

## Body mass effects on CTmax

ct_max_df_scaling <- ct_max_df %>%
  filter(species != "olma",
         mass != ".") %>%
  mutate(mass_dbl = as.double(mass)) %>%
  group_by(species) %>%
  nest() %>%
  mutate(spps_scaling_mod_ctmax = data %>% purrr::map(~ lm(loe_temp_corrected~mass_dbl, data=.)),
         aov_scaling_mod_ctmax = spps_scaling_mod_ctmax %>% purrr::map(~ Anova(., type = "III")),
         median_mass = data %>% purrr::map_dbl(~ median(.$mass_dbl)),
         mean_mass = data %>% purrr::map_dbl(~ mean(.$mass_dbl)),
         range_mass = data %>% purrr::map_dbl(~ (max(.$mass_dbl) - min(.$mass_dbl))),
         range_fold_mass = data %>% purrr::map_dbl(~ (max(.$mass_dbl) / min(.$mass_dbl))),
         tidy_ctmax_aov = aov_scaling_mod_ctmax %>% purrr::map(broom::tidy),
         tidy_ctmax_lm = spps_scaling_mod_ctmax %>% purrr::map(broom::tidy),
         f_val_ctmax_aov = tidy_ctmax_aov %>% purrr::map_dbl(c(4,2)),
         f_df_num_ctmax_aov = tidy_ctmax_aov %>% purrr::map_dbl(c(3,2)),
         f_df_den_ctmax_aov = tidy_ctmax_aov %>% purrr::map_dbl(c(3,3)),
         p_val_ctmax_aov = tidy_ctmax_aov %>% purrr::map_dbl(c(5,2)),
         slope_ctmax = tidy_ctmax_lm %>% purrr::map_dbl(c(2,2)))

# P values and relevant info for ctmax ~ body mass
ct_max_df_scaling[,c(1,5,6,7,8,11,12,13,14)]

## Plot: CTmax for each species, all raw data WITH ENBI outlier
ct_max_df %>%
  mutate(species_plotting_abb = case_when(species=="arfe"~"A. fenestralis",
                                          species=="arha"~"A. harringtoni",
                                          species=="arla"~"A. lateralis",
                                          species=="clgl"~"C. globiceps",
                                          species=="enbi"~"E. bison",
                                          species=="hehe"~"H. hemilepidotus",
                                          species=="olma"~"O. maculosus",
                                          species=="scma"~"S. marmoratus")) %>%
ggplot(aes(x = species_plotting_abb, y = loe_temp_corrected)) + 
  geom_jitter(width = 0.15, size = 3) +
  scale_x_discrete(name = "Species") + 
  scale_y_continuous(name = expression(paste("CT"["max"]," (",degree,C,")"))) + 
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.title.x = element_text(size = rel(2.5),
                                    colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        axis.text.x = element_text(size = rel(1.25), 
                                   angle = 25, 
                                   vjust = 1, hjust = 1,
                                   face = "bold",
                                   colour = "black"),
        axis.title.y = element_text(size = rel(2.5),
                                    hjust = 0.35,
                                    colour = "black"),
        axis.text.y = element_text(size = rel(2.5),
                                   colour = "black"))

## Plot: CTmax for each species, all raw data WITHOUT ENBI outlier
ct_max_df_no_enbi_out <- ct_max_df[-37,]
View(ct_max_df)
View(ct_max_df_no_enbi_out)

ct_max_df_no_enbi_out %>%
  mutate(species_plotting_abb = case_when(species=="arfe"~"Artedius\nfenestralis",
                                          species=="arha"~"Artedius\nharringtoni",
                                          species=="arla"~"Artedius\nlateralis",
                                          species=="clgl"~"Clinocottus\nglobiceps",
                                          species=="enbi"~"Enophrys\nbison",
                                          species=="hehe"~"Hemilepidotus\nhemilepidotus",
                                          species=="olma"~"Oligocottus\nmaculosus",
                                          species=="scma"~"Scorpaenichthys\nmarmoratus")) %>%
  ggplot(aes(x = species_plotting_abb, y = loe_temp_corrected)) + 
  geom_jitter(width = 0.15, size = 4) + 
  scale_x_discrete(name = "Species") + 
  scale_y_continuous(name = expression(paste("CT"["max"]," (",degree,C,")"))) + 
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.title.x = element_text(size = rel(2.5),
                                    colour = "black",
                                    margin = margin(t = 30)),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        axis.text.x = element_text(size = rel(1.5), 
                                   #angle = 20, 
                                   #vjust = 1, hjust = 1,
                                   #face = "bold",
                                   colour = "black"),
        axis.title.y = element_text(size = rel(2.5),
                                    margin = margin(r = 20),
                                    colour = "black"),
        axis.text.y = element_text(size = rel(2.25),
                                   colour = "black"))

## Plot Beta_Pcrit_low/high ~ CTmax

## Beta pcrit: 12-16 degrees ~ CTmax

beta_pcrit_ct_max_data <- 
  ct_max_df_no_enbi_out %>%
  group_by(species) %>%
  dplyr::select(species, mass, species_plotting, distribution, max_depth, loe_temp_corrected) %>%
  summarise(avg_ct_max = mean(loe_temp_corrected)) %>%
  rename(spps = species) %>%
  mutate(species = case_when(spps == "olma"~"Oligocottus_maculosus",
                             spps == "clgl"~"Clinocottus_globiceps",
                             spps == "arha"~"Artedius_harringtoni",
                             spps == "arla"~"Artedius_lateralis",
                             spps == "arfe"~"Artedius_fenestralis",
                             spps == "hehe"~"Hemilepidotus_hemilepidotus",
                             spps == "scma"~"Scorpaenichthys_marmoratus",
                             spps == "enbi"~"Enophrys_bison")) %>%
  full_join(., lm_pcrit_smr_temp[lm_pcrit_smr_temp$species!="Blepsias_cirrhosus",c(1,10,11,12)], 
            by = "species") %>%
  mutate(species_plotting = case_when(species == "Oligocottus_maculosus" ~ "Oligocottus maculosus",
                                      species == "Clinocottus_globiceps" ~ "Clinocottus globiceps",
                                      species == "Artedius_harringtoni" ~ "Artedius harringtoni",
                                      species == "Artedius_lateralis" ~ "Artedius lateralis",
                                      species == "Artedius_fenestralis" ~ "Artedius fenestralis",
                                      species == "Enophrys_bison" ~ "Enophrys bison",
                                      species == "Hemilepidotus_hemilepidotus" ~ "Hemilepidotus hemilepidotus",
                                      species == "Scorpaenichthys_marmoratus" ~ "Scorpaenichthys marmoratus"))
  
## AOV on Beta_pcrit_low/high ~ CTmax
beta_pcrit_ct_max_data %>%
  lm(slope_pcrit_low_temps ~ avg_ct_max, data = .) %>%
  Anova(., type = "III")

beta_pcrit_ct_max_data %>%
  lm(slope_pcrit_high_temps ~ avg_ct_max, data = .) %>%
  Anova(., type = "III")


beta_pcrit_low_temp_ctmax_plot <- 
  beta_pcrit_ct_max_data %>%
  ggplot(aes(x = avg_ct_max, y = slope_pcrit_low_temps)) +
  geom_point(size = 5) +
  #stat_smooth(method = "lm") +
  scale_x_continuous(name = expression(paste("CT"["max"]," (",degree,C,")")),
                     limits = c(20, 30)) +
  scale_y_continuous(name = expression(paste(beta["P"]["crit"][" 12-16"][degree][C]," (Torr  ",degree,C^-1,")")),
                     limits = c(-5,15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)),
        plot.margin = unit(c(0.1,0.2,0.2,0.2), "in"))

beta_pcrit_high_temp_ctmax_plot <- 
  beta_pcrit_ct_max_data %>%
  ggplot(aes(x = avg_ct_max, y = slope_pcrit_high_temps)) +
  stat_smooth(method = "lm", size = 2, colour = "black", se = FALSE) +
  geom_point(size = 5) +
  scale_x_continuous(name = expression(paste("CT"["max"]," (",degree,C,")")),
                     limits = c(20, 30)) +
  scale_y_continuous(name = expression(paste(beta["P"]["crit"][" 16-20"][degree][C]," (Torr  ",degree,C^-1,")")),
                     limits = c(-5,15)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.5)),
        axis.text = element_text(size = rel(2.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(3), colour = "black"),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20)),
        plot.margin = unit(c(0.1,0.2,0.2,0.2), "in"))

grid.arrange(beta_pcrit_low_temp_ctmax_plot, beta_pcrit_high_temp_ctmax_plot, ncol=2)
# PNG exported: W = 1132 H = 654

## **********************************************************
##
## PGLS models: Beta pcrit, CTmax, Pcrit_12, Delta_Pcrit/SMR
##
## **********************************************************

## NOTE: Consider using MCMCglmm instead, should let me feed in RAW data
## instead of 1 mean per species - See Vikram's phylo methods link

## ***********************************
## Beta Pcrit - low temps vs Pcrit_12
## ***********************************

## Using CAPER pacakge
data_beta_pcrit12 <- 
  lm_pcrit_smr_temp %>%
  dplyr::select(species, slope_pcrit_low_temps, slope_pcrit_high_temps, mean_pcrit_12) %>%
  as.data.frame()

comp_data_beta_pcrit12_caper <- 
  comparative.data(mandic_phy, data_beta_pcrit12,
                   names.col = "species",
                   vcv.dim = 2,warn.dropped = TRUE)

# *************************************
# OLS model: beta LOW temps ~ Pcrit_12
ols_beta_low_temps_pcrit12_caper <- pgls(slope_pcrit_low_temps ~ mean_pcrit_12,
                                         data = comp_data_beta_pcrit12_caper)
summary(ols_beta_low_temps_pcrit12_caper)
# PGLS model: with lambda estimated using ML
pgls_beta_low_temps_pcrit12_caper <- pgls(slope_pcrit_low_temps ~ mean_pcrit_12,
                                         data = comp_data_beta_pcrit12_caper,
                                         lambda = "ML")
summary(pgls_beta_low_temps_pcrit12_caper)
lk_pgls_beta_low_temps_pcrit_12_caper <- pgls.profile(pgls_beta_low_temps_pcrit12_caper,
                                                      which = "lambda")
plot(lk_pgls_beta_low_temps_pcrit_12_caper) ## No peak but rises continuously to 0
# PIC model: using caper's CRUNCH function
pic_beta_low_temps_pcrit12_caper <- crunch(slope_pcrit_low_temps ~ mean_pcrit_12,
                                           data = comp_data_beta_pcrit12_caper)
summary(pic_beta_low_temps_pcrit12_caper)
#
# _____________________________________

# *************************************
# OLS model: beta HIGH temps ~ Pcrit_12
ols_beta_high_temps_pcrit12_caper <- pgls(slope_pcrit_high_temps ~ mean_pcrit_12,
                                         data = comp_data_beta_pcrit12_caper)
summary(ols_beta_high_temps_pcrit12_caper)
# PGLS model: with lambda estimated using ML
pgls_beta_high_temps_pcrit12_caper <- pgls(slope_pcrit_high_temps ~ mean_pcrit_12,
                                          data = comp_data_beta_pcrit12_caper,
                                          lambda = "ML")
summary(pgls_beta_high_temps_pcrit12_caper)
lk_pgls_beta_high_temps_pcrit_12_caper <- pgls.profile(pgls_beta_high_temps_pcrit12_caper,
                                                      which = "lambda")
plot(lk_pgls_beta_high_temps_pcrit_12_caper) ## No peak but rises continuously to 0
# PIC model: using caper's CRUNCH function
pic_beta_high_temps_pcrit12_caper <- crunch(slope_pcrit_high_temps ~ mean_pcrit_12,
                                           data = comp_data_beta_pcrit12_caper)
summary(pic_beta_high_temps_pcrit12_caper)
#
# _____________________________________

# *************************************
# Delta Pcrit_12-20 ~ Delta SMR_12-20
# *************************************

data_delta_pcrit_smr <- lm_pcrit_vs_smr_delta_smr_data %>% as.data.frame()

comp_data_delta_pcrit_smr_caper <- 
  comparative.data(mandic_phy, data_delta_pcrit_smr,
                   names.col = "species",
                   vcv.dim = 2,warn.dropped = TRUE)

# *************************************
# OLS model: delta_pcrit ~ delta_smr
ols_delta_pcrit_delta_smr_caper <- pgls(delta_pcrit_12_20 ~ delta_smr_12_20,
                                         data = comp_data_delta_pcrit_smr_caper)
summary(ols_delta_pcrit_delta_smr_caper)
# PGLS model: with lambda estimated using ML
pgls_delta_pcrit_delta_smr_caper <- pgls(delta_pcrit_12_20 ~ delta_smr_12_20,
                                          data = comp_data_delta_pcrit_smr_caper,
                                          lambda = "ML")
summary(pgls_delta_pcrit_delta_smr_caper)
lk_pgls_delta_pcrit_delta_smr_caper <- pgls.profile(pgls_delta_pcrit_delta_smr_caper,
                                                      which = "lambda")
plot(lk_pgls_delta_pcrit_delta_smr_caper) ## No peak but rises continuously to 0
# PIC model: using caper's CRUNCH function
pic_delta_pcrit_delta_smr_caper <- crunch(delta_pcrit_12_20 ~ delta_smr_12_20,
                                           data = comp_data_delta_pcrit_smr_caper)
summary(pic_delta_pcrit_delta_smr_caper)
#
# _____________________________________

# **********************************************
# Comparative dataframe: beta LOW temps ~ CTmax

data_beta_ctmax <- 
  lm_pcrit_smr_temp %>%
  filter(species != "Blepsias_cirrhosus") %>%
  dplyr::select(species, slope_pcrit_low_temps, slope_pcrit_high_temps, mean_pcrit_12) %>%
  mutate(mean_ct_max = case_when(species=="Oligocottus_maculosus"~29.0,
                                 species=="Clinocottus_globiceps"~27.7,
                                 species=="Artedius_harringtoni"~24.7,
                                 species=="Artedius_lateralis"~26.8,
                                 species=="Artedius_fenestralis"~26.1,
                                 species=="Hemilepidotus_hemilepidotus"~27.3,
                                 species=="Scorpaenichthys_marmoratus"~27.1,
                                 species=="Enophrys_bison"~21.6)) %>%
  as.data.frame()

comp_data_beta_ctmax_caper <- 
  comparative.data(mandic_ctmax_phy, data_beta_ctmax,
                   names.col = "species",
                   vcv.dim = 2,warn.dropped = TRUE)

# *************************************
# OLS model: beta LOW temps ~ CTmax
ols_beta_low_temps_ctmax_caper <- pgls(slope_pcrit_low_temps ~ mean_ct_max,
                                          data = comp_data_beta_ctmax_caper)
summary(ols_beta_low_temps_ctmax_caper)
# PGLS model: with lambda estimated using ML
pgls_beta_low_temps_ct_max_caper <- pgls(slope_pcrit_low_temps ~ mean_ct_max,
                                           data = comp_data_beta_ctmax_caper,
                                           lambda = "ML")
summary(pgls_beta_low_temps_ct_max_caper)
lk_pgls_beta_low_temps_ct_max_caper <- pgls.profile(pgls_beta_low_temps_ct_max_caper,
                                                       which = "lambda")
plot(lk_pgls_beta_low_temps_ct_max_caper) ## No peak but rises continuously to 0
# PIC model: using caper's CRUNCH function
pic_beta_low_temps_ct_max_caper <- crunch(slope_pcrit_low_temps ~ mean_ct_max,
                                            data = comp_data_beta_ctmax_caper)
summary(pic_beta_low_temps_ct_max_caper)
#
# _____________________________________

# *************************************
# BM model: beta HIGH  temps ~ CTmax
pgls_bm_beta_high_temps_ctmax_caper <- pgls(slope_pcrit_high_temps ~ mean_ct_max,
                                       data = comp_data_beta_ctmax_caper)
summary(pgls_bm_beta_high_temps_ctmax_caper)
# PGLS model: with lambda estimated using ML
pgls_pagel_beta_high_temps_ct_max_caper <- pgls(slope_pcrit_high_temps ~ mean_ct_max,
                                         data = comp_data_beta_ctmax_caper,
                                         lambda = "ML")
summary(pgls_pagel_beta_high_temps_ct_max_caper)
lk_pgls_pagel_beta_high_temps_ct_max_caper <- pgls.profile(pgls_pagel_beta_high_temps_ct_max_caper,
                                                     "lambda")
plot(lk_pgls_pagel_beta_high_temps_ct_max_caper) ## No peak but rises continuously to 0
# PIC model: using caper's CRUNCH function
pic_beta_high_temps_ct_max_caper <- crunch(slope_pcrit_high_temps ~ mean_ct_max,
                                          data = comp_data_beta_ctmax_caper)
summary(pic_beta_high_temps_ct_max_caper)

## Compare models using AIC (lower is beter)
AIC(pgls_bm_beta_high_temps_ctmax_caper)
AIC(pgls_pagel_beta_high_temps_ct_max_caper)

## Get the AIC weights
## Create a named vector
aic_all <- c(AIC(pgls_bm_beta_high_temps_ctmax_caper),
             AIC(pgls_pagel_beta_high_temps_ct_max_caper))
names(aic_all) <- c("bm", "lambda")
aicw(aic_all)

anova.pgls(pgls_bm_beta_high_temps_ctmax_caper,
           pgls_pagel_beta_high_temps_ct_max_caper)

lrtest <- function(model1, model2){
  lik1 <- logLik(model1)
  lik2 <- logLik(model2)
  LR <- -2*(lik1-lik2)
  degf <-  attr(lik2,"df") - attr(lik1,"df")
  P <- pchisq(LR, df = degf, lower.tail = FALSE)
  cat(paste("Likelihood ratio = ",
            signif(LR, 5),"(df=",degf,") P = ",
            signif(P,4),"\n",sep = ""))
  invisible(list(likelihood.ratio=LR,p=P))
}

lrtest(pgls_bm_beta_high_temps_ctmax_caper,
       pgls_pagel_beta_high_temps_ct_max_caper)

#
# _____________________________________


## *******************
##
## Using APE package
##
## *******************
data_beta_pcrit12_gls <- lm_pcrit_smr_temp %>%
  dplyr::select(species, slope_pcrit_low_temps, slope_pcrit_high_temps, mean_pcrit_12) %>%
  column_to_rownames(var = "species")

## Regular OLS
ols_beta_low_temps_pcrit12 <- gls(slope_pcrit_low_temps ~ mean_pcrit_12,
      data = data_beta_pcrit12_gls,
      method = "ML")

## PGLS assuming a Brownian correlation
pgls_bm_beta_low_temps_pcrit12 <- gls(slope_pcrit_low_temps ~ mean_pcrit_12,
      data=data_beta_pcrit12_gls,
      correlation=corBrownian(phy=mandic_phy),
      method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_lambda_beta_low_temps_pcrit12 <- gls(slope_pcrit_low_temps ~ mean_pcrit_12,
      data=data_beta_pcrit12_gls,
      correlation=corPagel(value=0,phy=mandic_phy,fixed=TRUE),
      method="ML")

lambda <- seq(0,1, length.out = 500)
lik <- sapply(lambda, function(lambda) logLik(gls(slope_pcrit_low_temps ~ mean_pcrit_12,
                                                  correlation = 
                                                    corPagel(value = lambda, 
                                                             phy = mandic_phy, 
                                                             fixed = TRUE),
                                                  data = data_beta_pcrit12_gls)))
plot(lik ~ lambda, type = "l", 
     main = expression(paste("Likelihood Plot for ",lambda)), 
     ylab = "Log Likelihood", 
     xlab = expression(lambda))
tibble(lambda = lambda, likelihood = lik) %>% filter(likelihood == max(likelihood))
#abline(v = pgls_lambda_beta_high_temps_ctmax_gls$modelStruct$corStruct, col = "red")

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(ols_beta_low_temps_pcrit12, pgls_lambda_beta_low_temps_pcrit12)

plot(ols_beta_low_temps_pcrit12)
summary(ols_beta_low_temps_pcrit12)
anova(ols_beta_low_temps_pcrit12)

## ************************************
## Beta Pcrit - high temps vs Pcrit_12
## ************************************

## Regular OLS
ols_beta_high_temps_pcrit12 <- gls(slope_pcrit_high_temps ~ mean_pcrit_12,
      data = data_beta_pcrit12_gls,
      method = "ML")

## PGLS assuming a Brownian correlation
pgls_bm_beta_high_temps_pcrit12 <- gls(slope_pcrit_high_temps ~ mean_pcrit_12,
      data=data_beta_pcrit12_gls,
      correlation=corBrownian(phy=mandic_phy),
      method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_lambda_beta_high_temps_pcrit12 <- gls(slope_pcrit_high_temps ~ mean_pcrit_12,
      data=data_beta_pcrit12_gls,
      correlation=corPagel(value=1,phy=mandic_phy,fixed=TRUE),
      method="ML")

lambda <- seq(0,1, length.out = 500)
lik <- sapply(lambda, function(lambda) logLik(gls(slope_pcrit_high_temps ~ mean_pcrit_12,
                                                  correlation = 
                                                    corPagel(value = lambda, 
                                                             phy = mandic_phy, 
                                                             fixed = TRUE),
                                                  data = data_beta_pcrit12_gls)))
plot(lik ~ lambda, type = "l", 
     main = expression(paste("Likelihood Plot for ",lambda)), 
     ylab = "Log Likelihood", 
     xlab = expression(lambda))
lk_plot_data <- tibble(lambda = lambda, likelihood = lik) 
# filter(likelihood == max(likelihood))
ggplot(lk_plot_data, aes(x = lambda, y = likelihood)) + 
  geom_point(size = rel(0.5))
#abline(v = pgls_lambda_beta_high_temps_ctmax_gls$modelStruct$corStruct, col = "red")

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(ols_beta_high_temps_pcrit12, pgls_lambda_beta_high_temps_pcrit12)

summary(pgls_lambda_beta_high_temps_pcrit12)
plot(pgls_lambda_beta_high_temps_pcrit12)

plot(ols_beta_high_temps_pcrit12)
summary(ols_beta_high_temps_pcrit12)
anova(ols_beta_high_temps_pcrit12)

## ************************************
## Delta Pcrit ~ Delta SMR
## ************************************

lm_pcrit_vs_smr_delta_smr_data_spps_rn <- lm_pcrit_vs_smr_delta_smr_data %>%
  column_to_rownames(var = "species")

## Regular OLS
ols_delta_pcrit_delta_smr_gls <- 
  gls(delta_pcrit_12_20 ~ delta_smr_12_20,
      data = lm_pcrit_vs_smr_delta_smr_data_spps_rn,
      method = "ML")

## PGLS assuming a Brownian correlation
pgls_bm_delta_pcrit_delta_smr_gls <- 
  gls(delta_pcrit_12_20 ~ delta_smr_12_20,
      data=lm_pcrit_vs_smr_delta_smr_data_spps_rn,
      correlation=corBrownian(phy=mandic_phy),
      method="ML")

## PGLS assuming a Pagel's lambda correlation
pgls_delta_pcrit_delta_smr_gls <- 
  gls(delta_pcrit_12_20 ~ delta_smr_12_20,
      data=lm_pcrit_vs_smr_delta_smr_data_spps_rn,
      correlation=corPagel(value=0,phy=mandic_phy,fixed=TRUE),
      method="ML")

lambda <- seq(0,1, length.out = 500)
lik <- sapply(lambda, function(lambda) logLik(gls(delta_pcrit_12_20 ~ delta_smr_12_20,
                                                  correlation = 
                                                    corPagel(value = lambda, 
                                                             phy = mandic_phy, 
                                                             fixed = TRUE),
                                                  data = lm_pcrit_vs_smr_delta_smr_data_spps_rn)))
plot(lik ~ lambda, type = "l", 
     main = expression(paste("Likelihood Plot for ",lambda)), 
     ylab = "Log Likelihood", 
     xlab = expression(lambda))
lk_plot_data <- tibble(lambda = lambda, likelihood = lik) 
lk_plot_data %>% filter(likelihood == max(likelihood))

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(ols_delta_pcrit_delta_smr_gls, pgls_delta_pcrit_delta_smr_gls)

plot(ols_delta_pcrit_delta_smr_gls)
summary(ols_delta_pcrit_delta_smr_gls)
anova(ols_delta_pcrit_delta_smr_gls)

## ************************************
## Beta Pcrit - low temps vs CT_max
## ************************************

data_beta_ctmax_gls <- 
  lm_pcrit_smr_temp %>%
  filter(species != "Blepsias_cirrhosus") %>%
  dplyr::select(species, slope_pcrit_low_temps, slope_pcrit_high_temps, mean_pcrit_12) %>%
  mutate(mean_ct_max = case_when(species=="Oligocottus_maculosus"~29.0,
                                 species=="Clinocottus_globiceps"~27.7,
                                 species=="Artedius_harringtoni"~24.7,
                                 species=="Artedius_lateralis"~26.8,
                                 species=="Artedius_fenestralis"~26.1,
                                 species=="Hemilepidotus_hemilepidotus"~27.3,
                                 species=="Scorpaenichthys_marmoratus"~27.1,
                                 species=="Enophrys_bison"~21.6)) %>%
  column_to_rownames(var = "species")

## Regular OLS
ols_beta_low_temps_ctmax_gls <- gls(slope_pcrit_low_temps ~ mean_ct_max,
                                  data = data_beta_ctmax_gls,
                                  method = "ML")

## PGLS assuming a Brownian correlation
pgls_bm_beta_low_temps_ctmax_gls <- gls(slope_pcrit_low_temps ~ mean_ct_max,
                                      data=data_beta_ctmax_gls,
                                      correlation=corBrownian(phy=mandic_ctmax_phy),
                                      method="ML")

## PGLS assuming a Pagel's lambda correlation
## Model does not converge with `fixed=FALSE` - had to fix `value`
lambda <- seq(0,1, length.out = 500)
lik <- sapply(lambda, function(lambda) logLik(gls(slope_pcrit_low_temps ~ mean_ct_max,
                                                  correlation = 
                                                    corPagel(value = lambda, 
                                                             phy = mandic_ctmax_phy, 
                                                             fixed = TRUE),
                                                  data = data_beta_ctmax_gls)))
plot(lik ~ lambda, type = "l", 
     main = expression(paste("Likelihood Plot for ",lambda)), 
     ylab = "Log Likelihood", 
     xlab = expression(lambda))
lk_plot_data <- tibble(lambda = lambda, likelihood = lik) 
lk_plot_data %>% filter(likelihood == max(likelihood))

pgls_lambda_beta_low_temps_ctmax_gls <- gls(slope_pcrit_low_temps ~ mean_ct_max,
      data=data_beta_ctmax_gls,
      correlation=corPagel(value=0.343,phy=mandic_ctmax_phy,fixed=TRUE),
      method="ML")


## Using a likelihood ratio test to evaluate whether pgls is justified
anova(ols_beta_low_temps_ctmax_gls, pgls_lambda_beta_low_temps_ctmax_gls)

plot(ols_beta_low_temps_ctmax_gls)
summary(ols_beta_low_temps_ctmax_gls)
anova(ols_beta_low_temps_ctmax_gls)
logLik(ols_beta_low_temps_ctmax_gls)
summary(pgls_lambda_beta_low_temps_ctmax_gls)
anova(pgls_lambda_beta_low_temps_ctmax_gls)
logLik(pgls_lambda_beta_low_temps_ctmax_gls)


## ************************************
## Beta Pcrit - high temps vs CT_max
## ************************************

## Regular OLS
ols_beta_high_temps_ctmax_gls <- gls(slope_pcrit_high_temps ~ mean_ct_max,
                                    data = data_beta_ctmax_gls,
                                    method = "ML")

## PGLS assuming a Brownian correlation
pgls_bm_beta_high_temps_ctmax_gls <- gls(slope_pcrit_high_temps ~ mean_ct_max,
                                        data=data_beta_ctmax_gls,
                                        correlation=corBrownian(phy=mandic_ctmax_phy),
                                        method="ML")

## PGLS assuming a Pagel's lambda correlation
## Model does not converge with `fixed=FALSE` - had to fix `value`
lambda <- seq(0,1, length.out = 500)
lik <- sapply(lambda, function(lambda) logLik(gls(slope_pcrit_high_temps ~ mean_ct_max,
                                                  correlation = 
                                                    corPagel(value = lambda, 
                                                             phy = mandic_ctmax_phy, 
                                                             fixed = TRUE),
                                                  data = data_beta_ctmax_gls)))
plot(lik ~ lambda, type = "l", 
     main = expression(paste("Likelihood Plot for ",lambda)), 
     ylab = "Log Likelihood", 
     xlab = expression(lambda))
lk_plot_data <- tibble(lambda = lambda, likelihood = lik) 
lk_plot_data %>% filter(likelihood == max(likelihood))

pgls_lambda_beta_high_temps_ctmax_gls <- gls(slope_pcrit_high_temps ~ mean_ct_max,
                                            data=data_beta_ctmax_gls,
                                            correlation=corPagel(value=0.864,
                                                                 phy=mandic_ctmax_phy,
                                                                 fixed=FALSE),
                                            method="ML")

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(ols_beta_high_temps_ctmax_gls, pgls_lambda_beta_high_temps_ctmax_gls)

plot(ols_beta_high_temps_ctmax_gls)
summary(ols_beta_high_temps_ctmax_gls)
anova(ols_beta_high_temps_ctmax_gls)
logLik(ols_beta_high_temps_ctmax_gls)
plot(pgls_lambda_beta_high_temps_ctmax_gls)
summary(pgls_lambda_beta_high_temps_ctmax_gls)
anova(pgls_lambda_beta_high_temps_ctmax_gls)
logLik(pgls_lambda_beta_high_temps_ctmax_gls)

## ************************************************
##
## Phylogenetic ANOVA: Y ~ Tidepool occupancy
##
## ************************************************

# phylogenetic anova: Pcrit ~ tidepool occupancy (effectively a phylogenetic t-test...???)
sculpins_phy_data <- list(mandic_phy, lm_pcrit_smr_temp[,c("species",
                                                           "slope_pcrit_low_temps",
                                                           "slope_pcrit_high_temps",
                                                           "mean_pcrit_12",
                                                           "slope_mo2_low_temps",
                                                           "slope_mo2_high_temps",
                                                           "mean_mo2_12")])
sculpins_phy_data[[2]]$tp_occ <- c("Present",
                                   "Present",
                                   "Present",
                                   "Present",
                                   "Present",
                                   "Absent",
                                   "Absent",
                                   "Absent",
                                   "Absent")

sculpins_phy_data[[2]] <- column_to_rownames(sculpins_phy_data[[2]], var = "species")
sculpins_phy_data[[2]] <- as.data.frame(sculpins_phy_data[[2]])
## Dependent variables
pcrit12 <- sculpins_phy_data[[2]]$mean_pcrit_12
names(pcrit12) <- rownames(sculpins_phy_data[[2]])
beta_pcrit_12_16 <- sculpins_phy_data[[2]]$slope_pcrit_low_temps
names(beta_pcrit_12_16) <- rownames(sculpins_phy_data[[2]])
beta_pcrit_16_20 <- sculpins_phy_data[[2]]$slope_pcrit_high_temps
names(beta_pcrit_16_20) <- rownames(sculpins_phy_data[[2]])
mo2_12 <- sculpins_phy_data[[2]]$mean_mo2_12
names(mo2_12) <- rownames(sculpins_phy_data[[2]])
beta_mo2_12_16 <- sculpins_phy_data[[2]]$slope_mo2_low_temps
names(beta_mo2_12_16) <- rownames(sculpins_phy_data[[2]])
beta_mo2_16_20 <- sculpins_phy_data[[2]]$slope_mo2_high_temps
names(beta_mo2_16_20) <- rownames(sculpins_phy_data[[2]])
## Independent variable
tpo <- as.factor(sculpins_phy_data[[2]]$tp_occ)
names(tpo) <- rownames(sculpins_phy_data[[2]])
## Anova: Pcrit at 12C ~ Tidepool occupancy
aov_phylo_pcrit_12_tp_occ <- aov.phylo(pcrit12 ~ tpo, sculpins_phy_data[[1]], nsim = 1000)

## Anova: Beta_Pcrit_12-16C ~ Tidepool occupancy
aov_phylo_beta_pcrit12_16_tp_occ <- aov.phylo(beta_pcrit_12_16 ~ tpo, sculpins_phy_data[[1]], nsim = 1000)

## Anova: Beta_Pcrit_12-16C ~ Tidepool occupancy
aov_phylo_beta_pcrit16_20_tp_occ <- aov.phylo(beta_pcrit_16_20 ~ tpo, sculpins_phy_data[[1]], nsim = 1000)

## Anova: MO2 at 12C ~ Tidepool occupancy
aov_phylo_mo2_12_tp_occ <- aov.phylo(mo2_12 ~ tpo, sculpins_phy_data[[1]], nsim = 1000)

## Anova: Beta_mo2_12_16C ~ Tidepool occupancy
aov_phylo_beta_mo2_12_16_tp_occ <- aov.phylo(beta_mo2_12_16 ~ tpo, sculpins_phy_data[[1]], nsim = 1000)

## Anova: Beta_mo2_12-16C ~ Tidepool occupancy
aov_phylo_beta_mo2_16_20_tp_occ <- aov.phylo(beta_mo2_16_20 ~ tpo, sculpins_phy_data[[1]], nsim = 1000)


## ***************************
##     CT MAX ~ TPO
## ***************************

# phylogenetic anova: Pcrit ~ tidepool occupancy (effectively a phylogenetic t-test...???)
data_beta_ctmax_gls <- as.data.frame(data_beta_ctmax_gls)
sculpins_phy_ctmax_data <- list(mandic_ctmax_phy, data_beta_ctmax_gls)

sculpins_phy_ctmax_data[[2]]$tp_occ <- c("Present",
                                   "Present",
                                   "Present",
                                   "Present",
                                   "Present",
                                   "Absent",
                                   "Absent",
                                   "Absent")

## Dependent
ctmax <- sculpins_phy_ctmax_data[[2]]$mean_ct_max
names(ctmax) <- rownames(sculpins_phy_ctmax_data[[2]])

## Independent variable
tpo <- as.factor(sculpins_phy_ctmax_data[[2]]$tp_occ)
names(tpo) <- rownames(sculpins_phy_ctmax_data[[2]])
## Anova: Pcrit at 12C ~ Tidepool occupancy
aov_phylo_ctmax_tp_occ <- aov.phylo(ctmax ~ tpo, sculpins_phy_ctmax_data[[1]], nsim = 1000)

## ************************************************
##
## Phylogenetic signal in Traits and phenograms
##
## ************************************************

# phylogenetic anova: Pcrit ~ tidepool occupancy (effectively a phylogenetic t-test...???)
sculpins_phy_data <- list(mandic_phy, lm_pcrit_smr_temp[,c("species","slope_pcrit_low_temps",
                                                           "slope_pcrit_high_temps",
                                                           "mean_pcrit_12")])
sculpins_phy_data[[2]]$tp_occ <- c("Present",
                                   "Present",
                                   "Present",
                                   "Present",
                                   "Present",
                                   "Absent",
                                   "Absent",
                                   "Absent",
                                   "Absent")

sculpins_phy_data[[2]] <- column_to_rownames(sculpins_phy_data[[2]], var = "species")
sculpins_phy_data[[2]] <- as.data.frame(sculpins_phy_data[[2]])
## Dependent variables
pcrit12 <- sculpins_phy_data[[2]]$mean_pcrit_12
names(pcrit12) <- rownames(sculpins_phy_data[[2]])
beta_pcrit_12_16 <- sculpins_phy_data[[2]]$slope_pcrit_low_temps
names(beta_pcrit_12_16) <- rownames(sculpins_phy_data[[2]])
beta_pcrit_16_20 <- sculpins_phy_data[[2]]$slope_pcrit_high_temps
names(beta_pcrit_16_20) <- rownames(sculpins_phy_data[[2]])
## Independent variable
tpo <- as.factor(sculpins_phy_data[[2]]$tp_occ)
names(tpo) <- rownames(sculpins_phy_data[[2]])

## ************************************************
## Using lambda as indicator of phylogenetic signal
## in traits; 
## but see Kamilar and Cooper 2013 for caveats!
## ************************************************

## Test of phylogenetic signal in Pcrit
fitContinuous(sculpins_phy_data[[1]], pcrit12, model = "lambda")
# Phenogram for Pcrit at 12 degrees for all sculpins in this study
phenogram(sculpins_phy_data[[1]], pcrit12) ## Zero evidence for phylo signal in Pcrit at 12 C
## BUT low sample size (ie < 30) means very low power

## Test of phylogenetic signal in Beta_Pcrit_low_temps
fitContinuous(sculpins_phy_data[[1]], beta_pcrit_12_16, model = "lambda")
# Phenogram for Pcrit at 12 degrees for all sculpins in this study
phenogram(sculpins_phy_data[[1]], beta_pcrit_12_16)

## Test of phylogenetic signal in Beta_Pcrit_high_temps
fitContinuous(sculpins_phy_data[[1]], beta_pcrit_16_20, model = "lambda")
# Phenogram for Pcrit at 12 degrees for all sculpins in this study
phenogram(sculpins_phy_data[[1]], beta_pcrit_16_20)

## *******************************
##     CT MAX phylogenetic signal
## *******************************

# phylogenetic anova: Pcrit ~ tidepool occupancy (effectively a phylogenetic t-test...???)
data_beta_ctmax_gls <- as.data.frame(data_beta_ctmax_gls)
sculpins_phy_ctmax_data <- list(mandic_ctmax_phy, data_beta_ctmax_gls)

sculpins_phy_ctmax_data[[2]]$tp_occ <- c("Present",
                                         "Present",
                                         "Present",
                                         "Present",
                                         "Present",
                                         "Absent",
                                         "Absent",
                                         "Absent")

## Dependent
ctmax <- sculpins_phy_ctmax_data[[2]]$mean_ct_max
names(ctmax) <- rownames(sculpins_phy_ctmax_data[[2]])

## Test of phylogenetic signal in ctmax
fitContinuous(sculpins_phy_ctmax_data[[1]], ctmax, model = "lambda")
# Phenogram for Pcrit at 12 degrees for all sculpins in this study
phenogram(sculpins_phy_ctmax_data[[1]], ctmax)

################################################################
##
## Potential new figure?
##
################################################################
mass_corr_smr_pcrit_data %>% 
  filter(species %in% c("Enophrys_bison","Oligocottus_maculosus")) %>%
  dplyr::select(species,pcrit.r,smr.mass.corr.ms, temp) %>%
  group_by(species, temp) %>%
  mutate(mean_smr = mean(smr.mass.corr.ms),
         sd_smr = sd(smr.mass.corr.ms),
         n_smr = length(smr.mass.corr.ms),
         sem_smr = sd_smr/sqrt(n_smr),
         mean_pcrit = mean(pcrit.r),
         sd_pcrit = sd(pcrit.r),
         n_pcrit = length(pcrit.r),
         sem_pcrit = sd_pcrit/sqrt(n_pcrit)) %>%
  ggplot() + 
  geom_point(aes(x = pcrit.r, y = smr.mass.corr.ms, shape = species, colour = temp), size = 2) +
  geom_line(aes(x = mean_pcrit, y = mean_smr, colour = temp, group = species), size = 2) +
  scale_x_continuous(name = expression(paste("Po"[2]," (torr)")),
                     limits = c(0,160),
                     breaks = seq(0,160,20)) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,6),
                     breaks = seq(0,6,1))
