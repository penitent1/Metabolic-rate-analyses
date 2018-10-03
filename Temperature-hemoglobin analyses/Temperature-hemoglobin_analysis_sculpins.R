## **********************************************************
## Temperature and CO2 effects on hemoglobin-oxygen binding
## **********************************************************

## Relevant packages
library(tidyverse)
library(ape)
library(caper)
library(geiger)
library(phytools)
library(MCMCglmm)
library(visreg)
library(nlme)
library(broom)
library(car)

## Import data
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Temperature-hemoglobin analyses")
temp_hb_p50 <- read_csv("Temp-CO2_Hb-O2_binding_Collab_PhilMorrison_PreliminaryData.csv")

## Plot P50 by temperature, using CO2 = 5%
pd <- position_dodge(width = 1.5)

temp_hb_p50 %>%
  filter(percent_co2 == 0.5) %>%
  ggplot(aes(x = temperature, y = p50_mean_torr, colour = species)) +
  geom_point(position = pd)

## Summarise data by temperature, CO2, and species
temp_hb_p50_mean_co2_5 <- temp_hb_p50 %>%
  filter(percent_co2 == 0.5) %>%
  group_by(species, temperature, percent_co2) %>%
  summarise(avg_p50 = mean(p50_mean_torr),
            sd_p50 = sd(p50_mean_torr),
            n_p50 = length(p50_mean_torr),
            sem_p50 = sd_p50/sqrt(n_p50),
            avg_n50 = mean(n50_hill_mean),
            sd_n50 = sd(n50_hill_mean),
            n_n50 = length(n50_hill_mean),
            sem_n50 = sd_n50/sqrt(n_n50))

ggplot(temp_hb_p50_mean_co2_5, aes(x = temperature, y = avg_p50, colour = species)) +
geom_point(size = 5, position = pd) +
geom_errorbar(aes(ymin = avg_p50 - sd_p50, ymax = avg_p50 + sd_p50), position = pd)

## Import and summarise MO2, pcrit data for sculpins from temperature - Pcrit study

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

## Summarise by species

mass_corr_smr_pcrit_data_spps_mean <- mass_corr_smr_pcrit_data %>%
  group_by(species, temp) %>%
  dplyr::select(species, temp, pcrit.r, smr.mass.corr.ms) %>%
  summarise(avg_pcrit = mean(pcrit.r),
            avg_rmr = mean(smr.mass.corr.ms))
