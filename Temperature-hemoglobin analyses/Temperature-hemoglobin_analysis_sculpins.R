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

mass_corr_smr_pcrit_data_spps_mean_arfe_blci_olma <- mass_corr_smr_pcrit_data_spps_mean %>%
  filter(species %in% c("Artedius_fenestralis",
                        "Blepsias_cirrhosus",
                        "Oligocottus_maculosus"),
         temp %in% c(12,20)) %>%
  ungroup()

temp_hb_p50 <- temp_hb_p50_mean_co2_5 %>% dplyr::select(species, 
                                                        temperature, 
                                                        avg_p50)

temp_hb_p50_12c_20c <- temp_hb_p50 %>% 
  filter(temperature %in% c(18,23)) %>% 
  spread(temperature, avg_p50) %>%
  mutate(`20` = mean(c(`18`,`23`))) %>%
  dplyr::select(species, `20`) %>%
  rename(avg_p50 = `20`) %>%
  mutate(temperature = 20) %>%
  bind_rows(., temp_hb_p50 %>% filter(temperature == 12)) %>%
  rename(temp = temperature) %>%
  ungroup()

pcrit_rmr_p50 <- left_join(mass_corr_smr_pcrit_data_spps_mean_arfe_blci_olma,
                           temp_hb_p50_12c_20c,
                           by = c("species", "temp"))

ggplot(pcrit_rmr_p50, aes(x = avg_p50, y = avg_pcrit, colour = species)) +
  geom_point(size = 5) +
  geom_line(size = 2) +
  geom_abline(slope = 1, intercept = 0, size = 1.2) +
  scale_x_continuous(name = expression(paste("Hemoglobin-O"[2], " P"[50], " (Torr)")),
                     limits = c(0,60)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,60)) +
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


## Carrying capacity of sculpin blood at 12 degrees
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Temperature-hemoglobin analyses")
blood_o2_capacity <- read_csv("Mandic_et_al_2009_blood_oxygen_capacity_sculpins.csv")

blood_o2_capacity <- 
 blood_o2_capacity %>%
  mutate(max_ml_o2_100ml_blood = hct * mchc * 4 * 0.022391)

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


blood_o2_capacity_pcrit <- lm_pcrit_smr_temp[,c(1,10,11,12)] %>%
  filter(species %in% c("Artedius_fenestralis",
                        "Blepsias_cirrhosus",
                        "Oligocottus_maculosus",
                        "Clinocottus_globiceps",
                        "Enophrys_bison",
                        "Artedius_lateralis")) %>%
  mutate(blood_ml_o2_capacity = case_when(species == "Artedius_fenestralis" ~ 8.46,
                                          species == "Oligocottus_maculosus" ~ 11.6,
                                          species == "Blepsias_cirrhosus" ~ 7.01,
                                          species == "Clinocottus_globiceps" ~ 10.7,
                                          species == "Artedius_lateralis" ~ 7.52,
                                          species == "Enophrys_bison" ~ 4.93))

ggplot(blood_o2_capacity_pcrit, aes(x = blood_ml_o2_capacity, 
                                    y = slope_pcrit_high_temps)) +
  geom_point(size = 7) +
  scale_x_continuous(name = expression(paste("Oxygen capacity (mL O"[2]," 100 mL blood"^-1,")")),
                     limits = c(0,20),
                     breaks = seq(0,20,5)) +
  scale_y_continuous(name = expression(paste(beta["P"]["crit,16-20"][degree]["C"]," (Torr)")),
                     limits = c(-2,14),
                     breaks = seq(-2,14,2)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title = element_text(size = rel(2.75)),
        axis.text.x = element_text(size = rel(3), colour = "black"),
        axis.text.y = element_text(size = rel(3), colour = "black"),
        axis.line = element_line(size = rel(2), colour = "black"),
        axis.ticks = element_line(size = rel(3.5), colour = "black"),
        axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 20))) +
  geom_smooth(method = "lm", formula = (y ~ exp(1/x)), size = rel(2))

log_blood_o2_capacity_pcrit <- blood_o2_capacity_pcrit %>%
  mutate(log10_beta_pcrit = log10(slope_pcrit_high_temps),
         log10_blood_o2_capacity = log10(blood_ml_o2_capacity))

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
                    "Artedius_lateralis",
                    "Artedius_fenestralis", ## Ramon tree has typo here, should be "fenestralis"
                    "Enophrys_bison",
                    "Blepsias_cirrhosus")
mandic_phy <- drop.tip(mandic_um, setdiff(mandic_um$tip.label, keepers_mandic))
plot.phylo(mandic_phy)
#mandic_ctmax_phy <- drop.tip(mandic_phy, "Blepsias_cirrhosus") # No CTmax for BLCI
#plot.phylo(mandic_ctmax_phy)

## *********************************************
## MCMCglmm model: Beta Pcrit vs blood capacity
## *********************************************

mcmcglmm_blood_o2_capacity_df <- log_blood_o2_capacity_pcrit %>% as.data.frame()

inv_phylo_mandic <- inverseA(mandic_phy, nodes = "ALL")

prior <- list(G=list(G1=list(V=1, nu=0.02)), R=list(V=1, nu=0.02))
mcmcglmm_beta_high_temps_blood_O2_capacity <- MCMCglmm(log10_beta_pcrit ~ log10_blood_o2_capacity,
                                            random = ~ species,
                                            family = "gaussian",
                                            ginverse = list(species=inv_phylo_mandic$Ainv),
                                            prior = prior,
                                            data = mcmcglmm_blood_o2_capacity_df,
                                            nitt = 5000000, burnin = 1000, thin = 500)

summary(mcmcglmm_beta_high_temps_blood_O2_capacity)

mcmcglmm_beta_high_temps_blood_O2_capacity_lambda <- (mcmcglmm_beta_high_temps_blood_O2_capacity$VCV[,1])/(mcmcglmm_beta_high_temps_blood_O2_capacity$VCV[,1]+
                                                                                                             mcmcglmm_beta_high_temps_blood_O2_capacity$VCV[,2])
plot(mcmcglmm_beta_high_temps_blood_O2_capacity_lambda);mean(mcmcglmm_beta_high_temps_blood_O2_capacity_lambda)


blood_o2_capacity_pcrit %>%
  mutate(log10_beta_pcrit = log10(slope_pcrit_high_temps),
         log10_blood_o2_capacity = log10(blood_ml_o2_capacity)) %>%
  ggplot(aes(x = log10_blood_o2_capacity, y = log10_beta_pcrit)) +
  geom_point() +
  geom_abline(intercept = 2.2108, slope = -1.8495) ## From MCMCglmm model

blood_o2_capacity_pcrit %>%
  ggplot(aes(x = blood_ml_o2_capacity, y = slope_pcrit_high_temps)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ 162.48 + x^-1.8495)
