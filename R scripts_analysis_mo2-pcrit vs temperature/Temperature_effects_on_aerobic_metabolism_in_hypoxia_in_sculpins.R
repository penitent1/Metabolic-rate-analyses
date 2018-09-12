library(tidyverse)
library(ggthemes)
library(mclust)
library(shape)
library(StreamMetabolism)
library(fishMO2)
#library(cowplot)
library(ape)
library(caper)
library(geiger)
library(MCMCglmm)
library(visreg)
library(nlme)
library(broom)
library(car)

setwd("C:/Users/derek/Documents/Metabolic-rate-analyses")
pcrit_smr_data_summary <- read.csv("SculpinPcritData_ComparativeAnalysisFormat_withPcritSlopes.csv", stringsAsFactors = FALSE,
                                   strip.white = TRUE, na.strings = c("NA","."))
pcrit_smr_data_summary <- as_tibble(pcrit_smr_data_summary)
pcrit_smr_data_summary <- pcrit_smr_data_summary %>%
  mutate(smr.raw = smr.best.ms*mass.g) %>%
  filter(!is.na(pcrit.r), !is.na(smr.best.ms), trial.no==1)

## Open ct_max data and correct raw loe temperatures
# NOTE! I don't have CTmax data for BLCI
#setwd("C:/Users/derek/Documents/Metabolic-rate-analyses")
ct_max_df <- read_csv(file.choose()) %>% # Use data file that includes trials 1 and 2 for scma and enbi
  filter(species != "rhri") %>% ## Remove Grunt sculpins
  mutate(loe_temp_corrected = (loe_temp - 0.2909)/0.9857) ## Correcting for calibration; based on test of temp probe labelled "2"

ct_max_data_spps_summary <- ct_max_df %>%
  group_by(species) %>%
  dplyr::summarise(ct_max_avg = mean(loe_temp_corrected))

#######################
#######################                        

pcrit_smr_data_summary
ct_max_data_spps_summary


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
         p_val_smr_aov = tidy_smr_aov %>% purrr::map_dbl(c(5,2)),
         p_val_pcrit_aov = tidy_pcrit_aov %>% purrr::map_dbl(c(5,2)),
         slope_smr = tidy_smr_lm %>% purrr::map_dbl(c(2,2)))

## Adjust smr for body mass effects on SMR:

mass_corr_smr_pcrit_data <- 
  lm_scaling_smr_pcrit %>%
  dplyr::select(species, data, median_mass, mean_mass, range_mass, range_fold_mass, 
                p_val_smr_aov, slope_smr) %>%
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
  geom_jitter(width = 0.15) +
  geom_line(aes(x = temp, y = mean_pcrit), size = 1.25) +
  scale_x_continuous(name = expression(paste("Temperature (",degree,C,")")),
                     limits = c(9,23)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,100)) +
  facet_wrap("species_plotting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(1), colour = "black"),
        strip.text = element_text(face = "bold", size = rel(1.5)))

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
         data_mean_pcrit_12 = data %>% purrr::map(~ filter(., temp == 12)),
         lm_pcrit_low_temps = data_low_temps %>% purrr::map(~ lm(pcrit.r ~ temp, data=.)),
         lm_pcrit_high_temps = data_high_temps %>% purrr::map(~ lm(pcrit.r ~ temp, data=.)),
         tidy_pcrit_low_temps = lm_pcrit_low_temps %>% purrr::map(broom::tidy),
         tidy_pcrit_high_temps = lm_pcrit_high_temps %>% purrr::map(broom::tidy),
         slope_pcrit_low_temps = tidy_pcrit_low_temps %>% purrr::map_dbl(c(2,2)),
         slope_pcrit_high_temps = tidy_pcrit_high_temps %>% purrr::map_dbl(c(2,2)),
         mean_pcrit_12 = data_mean_pcrit_12 %>% purrr::map_dbl(~ mean(.$pcrit.r)))

## Beta pcrit: 12-16 degrees ~ Pcrit at 12 degrees
lm_pcrit_smr_temp[,c(1,10,11,12)] %>%
  ggplot(aes(x = mean_pcrit_12, y = slope_pcrit_low_temps)) +
  geom_point(size = 5) +
  stat_smooth(method = "lm", size = 2) +
  scale_x_continuous(name = expression(paste("P"["crit"]," at 12",degree,C," (Torr)")),
                     limits = c(20, 50)) +
  scale_y_continuous(name = expression(paste(beta["P"]["crit"]," (Torr  ",degree,C^-1,")")),
                     limits = c(-2,12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(1), colour = "black"))

## Beta pcrit: 16-20 degrees ~ Pcrit at 12 degrees
lm_pcrit_smr_temp[,c(1,10,11,12)] %>%
  ggplot(aes(x = mean_pcrit_12, y = slope_pcrit_high_temps)) +
  geom_point(size = 5) +
  stat_smooth(method = "lm", size = 2) +
  scale_x_continuous(name = expression(paste("P"["crit"]," at 12",degree,C," (Torr)")),
                     limits = c(20, 50)) +
  scale_y_continuous(name = expression(paste(beta["P"]["crit"]," (Torr  ",degree,C^-1,")")),
                     limits = (c(-2, 12))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(1), colour = "black"))

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
  ggplot(aes(x = smr.mass.corr.ms, y = pcrit.r, colour = species_plotting)) +
  #geom_point() +
  stat_smooth(method = "lm", se = FALSE, size = rel(2)) +
  scale_x_continuous(name = expression(paste("Standard ",dot(M)[O][2]," (",mu,"mol",O[2]," g"^-1," hr"^-1)),
                     limits = c(0,6)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,100)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(1), colour = "black"),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5))) +
  labs(colour = "Species")

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
  geom_point() +
  stat_smooth(method = "lm") +
  scale_x_continuous(name = expression(paste("Standard ",dot(M)[O][2]," (",mu,"mol",O[2]," g"^-1," hr"^-1)),
                     limits = c(0,6)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,100)) +
  facet_wrap("species_plotting") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(1.5), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(1), colour = "black"),
        strip.text = element_text(face = "bold", size = rel(1.5)))

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
  mutate(delta_smr_12_data = data %>% purrr::map(~ filter(., temp == 12)),
         delta_smr_20_data = data %>% purrr::map(~ filter(., temp == 20)),
         delta_smr_12_20 = purrr::map2_dbl(delta_smr_20_data, delta_smr_12_data, ~ .x$avg_smr -.y$avg_smr))

delta_smr_temp_data  
