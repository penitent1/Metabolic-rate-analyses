library(tidyverse)
library(ggthemes)
library(mclust)
library(shape)
library(StreamMetabolism)
library(fishMO2)
library(cowplot)
library(ape)
library(caper)
library(geiger)
library(MCMCglmm)
library(visreg)
library(nlme)
library(broom)

setwd("C:/Users/derek/Documents/Metabolic-rate-analyses")
pcrit_smr_data_summary <- read.csv("SculpinPcritData_ComparativeAnalysisFormat_withPcritSlopes.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
pcrit_smr_data_summary <- as_tibble(pcrit_smr_data_summary)
pcrit_smr_data_summary <- pcrit_smr_data_summary %>%
  mutate(smr.raw = smr.best.ms*mass.g) %>%
  filter(!is.na(pcrit.r), !is.na(smr.best.ms), trial.no==1)

setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Raw MO2 slope data")
mo2_data <- read.csv("Pcrit-vs-temp_Sculpin_study_IndividualMO2data.csv", stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c("NA","."))
mo2_data <- as_tibble(mo2_data)
mo2_data

#mo2_data <- transform(mo2_data, date = paste(date.day, date.month, date.year))
#mo2_data <- as_tibble(mo2_data)
mo2_data <- left_join(mo2_data, pcrit_smr_data_summary, by = c("probe" = "probe",
                                                      "spps" = "spps",
                                                      "date" = "date"))
mo2_data <- mo2_data %>%
  filter(!is.na(pcrit.r), !is.na(smr.best.ms)) %>%
  mutate(mo2.raw = mo2.ms*mass.g,
         species = if_else(spps=="olma","Oligocottus_maculosus",
                           if_else(spps=="scma", "Scorpaenichthys_marmoratus",
                                   if_else(spps=="blci", "Blepsias_cirrhosus",
                                           if_else(spps=="enbi", "Enophrys_bison",
                                                   if_else(spps=="hehe", "Hemilepidotus_hemilepidotus",
                                                           if_else(spps=="arla", "Artedius_lateralis",
                                                                   if_else(spps=="arha", "Artedius_harringtoni",
                                                                           if_else(spps=="arfe", "Artedius_fenestralis",
                                                                                   "Clinocottus_globiceps")))))))))


pcrit_folder <- "C:/Users/derek/Documents/Metabolic-rate-analyses/R fish MO2 O2crit csv files"
pcrit_file_list <- list.files(path = pcrit_folder)
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/R fish MO2 O2crit csv files")
## This works to read in all the pcrit data files
## From following website:
## https://stackoverflow.com/questions/26273555/inserting-file-names-as-column-values-in-a-data-frame?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
data2 <- lapply(pcrit_file_list, read.table, header = TRUE, sep = ",")
for (i in 1:length(data2)) {
  data2[[i]] <- cbind(data2[[i]], pcrit_file_list[i])
}
pcrit_data <- do.call("rbind", data2)
colnames(pcrit_data)[c(1,2,3)] <- c("do_percent", "mo2", "file")

pcrit_data <- as_tibble(pcrit_data)
## Alternative (didn't work for me) from the following website:
## http://www.reed.edu/data-at-reed/resources/R/reading_and_writing.html
#pcrit_data <- do.call("rbind",
#                      lapply(pcrit_file_list,
#                             function(x)
#                               read.csv(paste(pcrit_folder, x, sep = ''),
#                               stringsAsFactors = FALSE)))

pcrit_data

# Extract the date from filename
# See r4ds online e-book, Section 14.3: Matching patterns with regular expressions
#str_extract(pcrit_data$file, "(\\d|\\d\\d)[A-Z][a-z][a-z]\\d\\d\\d\\d")[1:10] # extracts date from filename
#str_extract(pcrit_data$file, "(\\d|\\d{2})[A-Z][a-z]{2}\\d{4}")[1:10] # same thing but shorter code

pcrit_data$trial_date <- str_extract(pcrit_data$file, "(\\d|\\d{2})[A-Z][a-z]{2}\\d{4}")

pcrit_data$trial_date <- case_when(
  pcrit_data$trial_date == str_extract(pcrit_data$trial_date, "\\d{2}[A-Z][a-z]{2}\\d{4}") ~
    str_c(
      str_sub(pcrit_data$trial_date, start = 1, end = 2),
      "-",
      str_sub(pcrit_data$trial_date, start = 3, end = 5),
      "-",
      str_sub(pcrit_data$trial_date, start = 6, end = 9)
    ),
  pcrit_data$trial_date == str_extract(pcrit_data$trial_date, "\\d[A-Z][a-z]{2}\\d{4}") ~
    str_c(
      str_sub(pcrit_data$trial_date, start = 1, end = 1),
      "-",
      str_sub(pcrit_data$trial_date, start = 2, end = 4),
      "-",
      str_sub(pcrit_data$trial_date, start = 5, end = 8)
    )
)

ex2 <- "20Nov2020"
ex2_dash <- str_c(
  str_sub(ex2, start = 1, end = 2),
  "-",
  str_sub(ex2, start = 3, end = 5),
  "-",
  str_sub(ex2, start = 6, end = 9)
)


# Extract the probe from filename
str_extract(pcrit_data$file, "NFB00[0-9]{2}")[1:10]

pcrit_data$probe <- str_extract(pcrit_data$file, "NFB00[0-9]{2}")
pcrit_data

# Extract the spps from filename
str_extract(pcrit_data$file, "[A-Z]{4}")[1:10]

pcrit_data$spps <- str_extract(pcrit_data$file, "[A-Z]{4}")
pcrit_data$spps <- tolower(pcrit_data$spps) ## Convert spps names "to lower" case letters

# Extract the temperature from filename
str_extract(pcrit_data$file, "\\s[0-9]{2}")[1:10]

pcrit_data$temp_c <- str_extract(pcrit_data$file, "\\s[0-9]{2}")
pcrit_data$temp_c <- as.numeric(str_extract(pcrit_data$temp_c, "[0-9]{2}"))

pcrit_data <- pcrit_data %>%
  mutate(species = if_else(spps=="olma","Oligocottus_maculosus",
                           if_else(spps=="scma", "Scorpaenichthys_marmoratus",
                                   if_else(spps=="blci", "Blepsias_cirrhosus",
                                           if_else(spps=="enbi", "Enophrys_bison",
                                                   if_else(spps=="hehe", "Hemilepidotus_hemilepidotus",
                                                           if_else(spps=="arla", "Artedius_lateralis",
                                                                   if_else(spps=="arha", "Artedius_harringtoni",
                                                                           if_else(spps=="arfe", "Artedius_fenestralis",
                                                                                   "Clinocottus_globiceps")))))))))

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

mo2_data
pcrit_data
pcrit_smr_data_summary
ct_max_data

## SMR allometry, with species-specific linear models including at all temperatures
# Colours = species
ggplot(pcrit_smr_data_summary, aes(x = mass.g, y = smr.raw, colour = spps)) +
  geom_point()

ggplot(pcrit_smr_data_summary, aes(x = mass.g, y = log(smr.raw), colour = spps)) +
  geom_point()

ggplot(pcrit_smr_data_summary, aes(x = log(mass.g), y = log(smr.raw), colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(pcrit_smr_data_summary, aes(x = log(mass.g), y = log(smr.raw), colour = spps, group = spps)) +
  stat_smooth(method = "lm")

## Get linear models of ln(smr) ~ ln(mass), as well as tests of significance
## for model terms
smr_allo_exponents <- pcrit_smr_data_summary %>%
  group_by(species) %>%
  do(spps_scaling_mod = lm(log(smr.raw)~log(mass.g), data=.)) %>%
  tidy(spps_scaling_mod) %>%
  dplyr::select(term, estimate, std.error, p.value) %>%
  filter(term == "log(mass.g)") %>%
  rename(slope.b = estimate) %>%
  mutate(sig_allo_exp = if_else(p.value < 0.05, "yes", "no"))

ggplot(smr_allo_exponents, aes(x = species, 
                               y = slope.b, 
                               colour = sig_allo_exp)) +
  geom_point(size = 5) +
  geom_errorbar(aes(x = species, 
                    ymin = slope.b-std.error, 
                    ymax = slope.b+std.error)) +
  labs(x = expression(paste("Species")),
       y = expression(paste("Scaling coefficient")),
       colour = expression(paste("Significant \nexponent?")))

## ***************************************************************
##
##    Linear models and AOVs on scaling of body mass and smr/pcrit
##
## ***************************************************************

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
         tidy_pcrit_lm = spps_scaling_mod_pcrit %>% purrr::map(broom::tidy))

p_values_lm_smr_aov <- lm_scaling_smr_pcrit %>% dplyr::select(species, tidy_smr_aov) %>% unnest() %>%
  spread(term, p.value) %>% dplyr::select(species, `log(mass.g)`) %>% filter(!is.na(`log(mass.g)`)) %>%
  mutate(p_val_smr_sig = if_else(`log(mass.g)` < 0.05, "yes", "no")) %>%
  rename(scaling_slope_smr_p_value = `log(mass.g)`)

p_values_lm_pcrit_aov <- lm_scaling_smr_pcrit %>% dplyr::select(species, tidy_pcrit_aov) %>% unnest() %>%
  spread(term, p.value) %>% dplyr::select(species, `log(mass.g)`) %>% filter(!is.na(`log(mass.g)`)) %>%
  mutate(p_val_pcrit_sig = if_else(`log(mass.g)` < 0.05, "yes", "no")) %>%
  rename(scaling_slope_pcrit_p_value = `log(mass.g)`)

p_values_scaling_aov_smr_pcrit <- bind_cols(p_values_lm_smr_aov, p_values_lm_pcrit_aov) %>%
  dplyr::select(-species1)

## In Deutsch et al., body size correction of "HYPOXIA TOLERANCE" was only done if
## body mass range was greater than 3 fold
lm_scaling_smr_pcrit %>% 
  dplyr::select(species, range_fold_mass) %>% 
  mutate(range_bigger_3fold = if_else(range_fold_mass > 3, "yes", "no"))

## Hemilepidotus hemilepidotus is the only species with a significant effect
## of body mass on Pcrit, but the fold-range of body masses is 1.8.

summary(lm_scaling_smr_pcrit$spps_scaling_mod_pcrit[[6]])

# Checking the 95% CI for the slope: Estimate +/- (Std. Error * t value)
#> 0.5514 - (0.2537*2.173)
#[1] 0.0001099
#> 0.5514 + (0.2537*2.173)
#[1] 1.10269

## The 95% CI range effectively includes 0 so no body mass correction for Pcrit!


## Raw Pcrit vs temperature data

pcrit_smr_data_summary %>%
  group_by(species, temp) %>%
  mutate(mean_pcrit = mean(pcrit.r),
         mean_smr = mean(smr.raw),
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
                     limits = c(10,22)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," (Torr)")),
                     limits = c(0,100)) +
  facet_wrap("species_plotting")

lm_pcrit_temp <- pcrit_smr_data_summary %>%
  group_by(species) %>%
  nest() %>%
  mutate(data_low_temps = data %>% purrr::map(~ filter(., temp != 20)),
         data_high_temps = data %>% purrr::map(~ filter(., temp !=12)),
         data_mean_pcrit_12 = data %>% purrr::map(~ filter(., temp == 12)),
         lm_low_temps = data_low_temps %>% purrr::map(~ lm(pcrit.r ~ temp, data=.)),
         lm_high_temps = data_high_temps %>% purrr::map(~ lm(pcrit.r ~ temp, data=.)),
         tidy_low_temps = lm_low_temps %>% purrr::map(broom::tidy),
         tidy_high_temps = lm_high_temps %>% purrr::map(broom::tidy),
         slope_low_temps = tidy_low_temps %>% purrr::map_dbl(c(2,2)),
         slope_high_temps = tidy_high_temps %>% purrr::map_dbl(c(2,2)),
         mean_pcrit_12 = data_mean_pcrit_12 %>% purrr::map_dbl(~ mean(.$pcrit.r)))
  
lm_pcrit_temp[,c(1,10,11,12)] %>%
  ggplot(aes(x = mean_pcrit_12, y = slope_low_temps)) +
  geom_point(size = 3) +
  stat_smooth(method = "lm") +
  scale_x_continuous(name = expression(paste("P"["crit"]," at 12",degree,C," (Torr)")),
                     limits = c(20, 50)) +
  scale_y_continuous(name = expression(paste("P"["crit"],"-Temperature sensitivity (Torr ",degree,C^-1,")")))

lm_pcrit_temp[,c(1,10,11,12)] %>%
  ggplot(aes(x = mean_pcrit_12, y = slope_high_temps)) +
  geom_point(size = 3) +
  stat_smooth(method = "lm") +
  scale_x_continuous(name = expression(paste("P"["crit"]," at 12",degree,C," (Torr)")),
                     limits = c(20, 50)) +
  scale_y_continuous(name = expression(paste("P"["crit"],"-Temperature sensitivity (Torr ",degree,C^-1,")")))

## Versus CTmax

lm_pcrit_temp[,c(1,10,11,12)] %>%
  ggplot(aes(x = mean_pcrit_12, y = slope_low_temps)) +
  geom_point(size = 3) +
  stat_smooth(method = "lm") +
  scale_x_continuous(name = expression(paste("P"["crit"]," at 12",degree,C," (Torr)")),
                     limits = c(20, 50)) +
  scale_y_continuous(name = expression(paste("P"["crit"],"-Temperature sensitivity (Torr ",degree,C^-1,")")))

lm_pcrit_temp[,c(1,10,11,12)] %>%
  ggplot(aes(x = mean_pcrit_12, y = slope_high_temps)) +
  geom_point(size = 3) +
  stat_smooth(method = "lm") +
  scale_x_continuous(name = expression(paste("P"["crit"]," at 12",degree,C," (Torr)")),
                     limits = c(20, 50)) +
  scale_y_continuous(name = expression(paste("P"["crit"],"-Temperature sensitivity (Torr ",degree,C^-1,")")))







ggplot() +
stat_smooth(data = pcrit_smr_data_summary[pcrit_smr_data_summary$temp!=20,], 
              aes(x = temp, y = pcrit.r),
            method = "lm") +
stat_smooth(data = pcrit_smr_data_summary[pcrit_smr_data_summary$temp!=12,],
              aes(x = temp, y = pcrit.r),
            method = "lm") +
facet_wrap("species")


## Arrhenius plot: ln(pcrit.r) ~ 1/T

pcrit_smr_data_summary %>%
  mutate(temp_k = temp + 273) %>%
  ggplot(aes(x = (1000/(temp_k)), y = log(pcrit.r))) + 
  geom_point() +
  stat_smooth(method = "lm") +
  facet_wrap("species")



pcrit_smr_data_summary %>%
  group_by(species, temp) %>%
  mutate(mean_pcrit = mean(pcrit.r),
         mean_smr = mean(smr.best.ms)) %>%
  ggplot(aes(x = temp, y = smr.best.ms)) +
  geom_jitter(width = 0.15) +
  geom_line(aes(x = temp, y = mean_smr), size = 1.25) +
  facet_wrap("species")


## Raw Pcrit vs temperature data >>> CHECK HIGH ARLA, same individual???

## Bootstrap pcrit's, replot and see how slopes look! - VIKRAM

pcrit_smr_data_summary %>%
  group_by(species, temp) %>%
  mutate(mean_pcrit = mean(pcrit.r),
         mean_smr = mean(smr.raw)) %>%
  ggplot() +
  geom_line(aes(x = temp, y = pcrit.r, group = finclip.id)) +
  facet_wrap("species")

## ln(Pcrit) vs temperature

pcrit_smr_data_summary %>%
  ggplot(aes(x = temp, y = log(pcrit.r))) +
  geom_point() +
  facet_wrap("species")

## ln(Pcrit) vs ln(temperature)

pcrit_smr_data_summary %>%
  ggplot(aes(x = log(temp), y = log(pcrit.r))) +
  geom_point() +
  stat_smooth(method = "lm") +
  facet_wrap("species")

##

##  calculate q10's for Pcrit

q10_data_pcrit <- 
  pcrit_smr_data_summary %>%
  dplyr::select(species, temp, pcrit.r) %>%
  group_by(species, temp) %>%
  summarise(avg_pcrit = mean(pcrit.r)) %>%
  spread(temp, avg_pcrit) %>%
  mutate(q10_low = (`16`/`12`)^(10/(16-12)),
         q10_high = (`20`/`16`)^(10/(20-16)),
         q10_all = (`20`/`12`)^(10/(20-12))) %>%
  ggplot(aes(x = `12`, y = q10_all)) +
  geom_point(size = 3) +
  stat_smooth(method = "lm")

anova(lm(q10_low ~ `12`, data = q10_data_pcrit))
anova(lm(q10_high ~ `12`, data = q10_data_pcrit))
anova(lm(q10_all ~ `12`, data = q10_data_pcrit))

ct_max_q10_pcrit <- 
ct_max_data %>%
  dplyr::select(species, loe_temp_corrected) %>%
  group_by(species) %>%
  summarise(avg_ct_cmax = mean(loe_temp_corrected)) %>%
  bind_cols(., q10_data_pcrit[q10_data_pcrit$species!="Blepsias_cirrhosus",]) %>%
  dplyr::select(-species1) %>%
  ggplot(aes(x = avg_ct_cmax, y = q10_high)) +
  geom_point(size = 3) +
  stat_smooth(method = "lm")

anova(lm(q10_low ~ avg_ct_cmax, data = ct_max_q10_pcrit))
anova(lm(q10_high ~ avg_ct_cmax, data = ct_max_q10_pcrit))
anova(lm(q10_all ~ avg_ct_cmax, data = ct_max_q10_pcrit))

##

lm_scaling_smr_pcrit$tidy_smr_aov[[1]]$p.value[2]

lm_scaling_smr_pcrit$tidy_smr_aov[[1]] %>% magrittr::use_series(p.value) %>% magrittr::extract(2)

lm_scaling_smr_pcrit$tidy_smr_aov %$% p.value
  purrr::map(magrittr::use_series(., p.value)) %>% 
  purrr::map_dbl(maggritr::extract(2))

 
##
##
##
##
##
##
##
##
##
## UNDER CONSTRUCTION!
##
##
##
##
##
##
##
##
##







## Pcrit vs body mass
ggplot(pcrit_smr_data_summary, aes(x = mass.g, y = pcrit.r, colour = spps)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(pcrit_smr_data_summary, aes(x = log(mass.g), y = pcrit.r, colour = spps)) +
  geom_point() +
  stat_smooth(method = "lm")

ggplot(pcrit_smr_data_summary, aes(x = log(mass.g), y = log(pcrit.r), colour = spps)) +
  geom_point() +
  stat_smooth(method = "lm")

## Get linear models of ln(pcrit) ~ ln(mass), as well as tests of significance
## for model terms
pcrit_allo_exponents <- pcrit_smr_data_summary %>%
  group_by(species) %>%
  do(spps_scaling_mod = lm(log(pcrit.r)~log(mass.g), data=.)) %>%
  tidy(spps_scaling_mod) %>%
  dplyr::select(term, estimate, std.error, p.value) %>%
  filter(term == "log(mass.g)") %>%
  rename(slope.b = estimate) %>%
  mutate(sig_allo_exp = if_else(p.value < 0.05, "yes", "no"))

ggplot(pcrit_allo_exponents, aes(x = species, 
                               y = slope.b, 
                               colour = sig_allo_exp)) +
  geom_point(size = 5) +
  geom_errorbar(aes(x = species, 
                    ymin = slope.b-std.error, 
                    ymax = slope.b+std.error)) +
  labs(x = expression(paste("Species")),
       y = expression(paste("Scaling coefficient")),
       colour = expression(paste("Significant \nexponent?")))





## Plot of ln(smr) ~ ln(mass), all temperatures included, and a single all-species regression
ggplot(pcrit_smr_data_summary, aes(x = log(mass.g), y = log(smr.raw))) +
  geom_point() +
  stat_smooth(method = "lm")

scaling_lm <- lm(log(smr.raw) ~ log(mass.g), data = pcrit_smr_data_summary)
coef_scaling_lm <- coef(scaling_lm)
residuals(scaling_lm)

fish_35_logmass <- data.frame(mass.g = 35)

predict.lm(lm(log(smr.raw) ~ log(mass.g), data = pcrit_smr_data_summary), fish_35_logmass)
coef_scaling_lm[1] + coef_scaling_lm[2]*(log(35))


pcrit_smr_data_summary <- pcrit_smr_data_summary %>%
  mutate(smr_35g = smr.raw*exp(coef_scaling_lm[1]+coef_scaling_lm[2]*log(35/mass.g)),
         smr_35g_noint = smr.raw*exp(coef_scaling_lm[2]*log(35/mass.g)))
pcrit_smr_data_summary$smr_35g_add_resid <- (predict.lm(
  lm(log(smr.raw) ~ log(mass.g), data = pcrit_smr_data_summary), fish_35_logmass)) +
  residuals(scaling_lm)

fun_scaling <- function(x){
  y = coef_scaling_lm[1] + coef_scaling_lm[2]*x
}

ggplot(pcrit_smr_data_summary) +
  geom_point(aes(x = log(mass.g), y = log(smr.raw))) +
  geom_point(aes(x = log(mass.g), y = log(smr_35g_noint)), colour = "red") +
  geom_hline(aes(yintercept = 4.355586), colour = "red") +
  geom_vline(aes(xintercept = 3.555348)) +
  stat_function(fun = fun_scaling, lwd = 1.15) +
  geom_point(aes(x = log(mass.g), y = smr_35g_add_resid), colour = "purple")

# 35-g mass adjusted smr vs temperature
ggplot(pcrit_smr_data_summary) +
  geom_jitter(aes(x = temp, y = smr_35g_noint, colour = spps), width = 0.2)



## SMR allometry, All species at 12 degrees, colored by species
smr_allo_12c_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==12,], aes(x = mass.g, y = smr.raw, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE)
# at 16 degrees, colored by species
smr_allo_16c_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==16,], aes(x = mass.g, y = smr.raw, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE)
# at 20 degrees, colored by species
smr_allo_20c_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==20,], aes(x = mass.g, y = smr.raw, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = TRUE)

plot_grid(smr_allo_12c_plot, smr_allo_16c_plot, smr_allo_20c_plot, labels = c(12,16,20))

## Zooming in on smaller species
smr_allo_12c_small_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==12&pcrit_smr_data_summary$mass.g<200,], aes(x = mass.g, y = smr.raw, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)
smr_allo_16c_small_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==16&pcrit_smr_data_summary$mass.g<200,], aes(x = mass.g, y = smr.raw, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)
smr_allo_20c_small_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==20&pcrit_smr_data_summary$mass.g<200,], aes(x = mass.g, y = smr.raw, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)

plot_grid(smr_allo_12c_small_plot, smr_allo_16c_small_plot, smr_allo_20c_small_plot)

## CTmax allometry
## Checking whether body mass affects CTmax - doesn't look like it at all!
ggplot(ct_max_data[!is.na(ct_max_data$mass_g),], aes(x = mass_g, y = loe_temp_corrected, group = spps, colour = spps)) + geom_point()

ggplot(ct_max_data[(!is.na(ct_max_data$mass_g)&ct_max_data$mass_g<200),], aes(x = mass_g, y = loe_temp_corrected, group = spps, colour = spps)) + 
  geom_point() +
  stat_smooth(method = "lm")

## Pcrit allometry, All species at 12 degrees
pcrit_allo_12c_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==12,], aes(x = mass.g, y = pcrit.r, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)
# at 16 degrees
pcrit_allo_16c_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==16,], aes(x = mass.g, y = pcrit.r, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)
# at 20 degrees
pcrit_allo_20c_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==20,], aes(x = mass.g, y = pcrit.r, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)

plot_grid(pcrit_allo_12c_plot, pcrit_allo_16c_plot, pcrit_allo_20c_plot, labels = c(12,16,20))

## Zooming in on smaller species
pcrit_allo_12c_small_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==12&pcrit_smr_data_summary$mass.g<200,], aes(x = mass.g, y = pcrit.r, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)
pcrit_allo_16c_small_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==16&pcrit_smr_data_summary$mass.g<200,], aes(x = mass.g, y = pcrit.r, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)
pcrit_allo_20c_small_plot <- ggplot(pcrit_smr_data_summary[pcrit_smr_data_summary$temp==20&pcrit_smr_data_summary$mass.g<200,], aes(x = mass.g, y = pcrit.r, colour = spps, group = spps)) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE)

plot_grid(pcrit_allo_12c_small_plot, pcrit_allo_16c_small_plot, pcrit_allo_20c_small_plot)

############################################################################################

#         Log-linear allometric model: smr vs mass

smr_allo_model <- lm(log10(smr.raw)~log10(mass.g), data = smr_data[smr_data$temp.x==12,])
smr_allo_model
summary(smr_allo_model)
shapiro.test(smr_allo_model$residuals)
plot(smr_allo_model)

#smr_allo_nls_model <- nls(smr.raw~a*mass.g^b, data = smr_data[smr_data$temp.x==12,], start = list(a=1,b=1))
#smr_allo_nls_model
#summary(smr_allo_nls_model)
#shapiro.test(residuals(smr_allo_nls_model))
#plot(smr_allo_nls_model)

fun_sculpins_12c_loglin <- function(x){2.330269*x^0.894206}
#fun_sculpins_12c_nls <- function(x){4.79797*x^0.76553} Don't use, the error is clearly multiplicative so log-log model is appropriate

# Allometric exponent here is 0.89

smr_allo_12c_all_plot <- ggplot(smr_data[smr_data$temp.x==12,], aes(x = mass.g, y = smr.raw)) +
  geom_point() +
  stat_function(fun = fun_sculpins_12c_loglin, colour = "red") #+
  #stat_function(fun = fun_sculpins_12c_nls, colour = "blue")
smr_allo_12c_all_plot

smr_allo_12c_all_small_plot <- ggplot(smr_data[smr_data$temp.x==12&smr_data$mass.g<200,], aes(x = mass.g, y = smr.raw)) +
  geom_point() +
  stat_function(fun = fun_sculpins_12c_loglin, colour = "red") #+
  #stat_function(fun = fun_sculpins_12c_nls, colour = "blue")
smr_allo_12c_all_small_plot

#############
#############
#############
#############
########        PCRIT VS TEMP: Plots with exponential functions
#############
#############
#############
#############

pcrit_smr_data_summary_grp_spps <- pcrit_smr_data_summary %>% 
  group_by(species)


## Plotted with all species on a single plot
ggplot(pcrit_smr_data_summary_grp_spps, 
       aes(x = smr.best.ms,
           y = pcrit.r,
           colour = species)) +
  geom_point() +
  stat_smooth(method = "lm")

## Plotted with each species on a separate plot
ggplot(pcrit_smr_data_summary_grp_spps, 
       aes(x = smr.best.ms,
           y = pcrit.r)) +
  geom_point() +
  stat_smooth(method = "lm") +
  facet_wrap(~species)


# Now using only 1 data point per species
pcrit_smr_data_spps_summary <- pcrit_smr_data_summary %>%
  group_by(species, temp) %>%
  dplyr::summarise(smr.ms = mean(smr.best.ms),
                   smr.raw = mean(smr.raw),
                   pcrit.avg = mean(pcrit.r),
                   mass.g.avg = mean(mass.g))


#####################################################################
## Read sculpin phyologeny into R and prune for species in analysis
#####################################################################
mandic_2013_tree <- read.nexus(file.choose())
plot.phylo(mandic_2013_tree)
is.ultrametric(mandic_2013_tree) ## It is not :(
is.rooted(mandic_2013_tree) ## It's rooted!!

mandic_2013_tree_no_out <- drop.tip(mandic_2013_tree, "Satyrichthys_amiscus")

mandic_um <- chronopl(mandic_2013_tree_no_out, 1)
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

## Compute branch lengths, tree didn't have any
##ramon_phy <- compute.brlen(ramon_phy, method = "Grafen")


#####################################################

######
######
######
######
######
###
###         LINEAR MODELS: PCRIT VS TEMP, 1 model per species
###
######
######
######
######
######

###   ************************************************************
###   Set up dataframe for Pcrit-Temp slope (b_pt) VS Pcrit at 12C
###   ************************************************************

pcrit_temp_models <- pcrit_smr_data_summary %>%
  mutate(pcrit.r.kpa = pcrit.r*0.133) %>%
  group_by(species) %>%
  do(pcrit_temp_lm = lm(pcrit.r.kpa~temp, data = .))

library(broom) ## Great package that takes model fit coefficients and puts them in a data frame!

pcrit_temp_models_df <- pcrit_temp_models %>%
  tidy(pcrit_temp_lm) %>%
  dplyr::select(term, estimate) %>%
  spread(term, estimate) %>%
  rename(slope = temp)

pcrit_12c_all_spps <- pcrit_smr_data_spps_summary %>%
  filter(temp == 12) %>%
  transmute(pcrit.r.kpa = pcrit.avg*0.133)

pcrit_temp_slope_pcrit_12c_analysis_df <- full_join(pcrit_temp_models_df,
                                                    pcrit_12c_all_spps,
                                                    by = "species")
### pgls goes here!

beta_pcrit_12c_df_rn <- pcrit_temp_slope_pcrit_12c_analysis_df %>%
  ungroup(pcrit_temp_slope_pcrit_12c_analysis_df) %>%
  column_to_rownames(var = "species")

## Regular OLS
#ols <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, method="ML")
ols_beta_pcrit12 <- gls(slope~pcrit.r.kpa, 
                        data=beta_pcrit_12c_df_rn, 
                        method="ML")

## PGLS assuming a Brownian correlation
#pgls_bm <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, correlation = corBrownian(phy=anolis_phy), method="ML")
pgls_bm_beta_pcrit12 <- gls(slope~pcrit.r.kpa, 
                            data=beta_pcrit_12c_df_rn, 
                            correlation=corBrownian(phy=mandic_phy), 
                            method="ML")

## PGLS assuming a Pagel's lambda correlation
#pgls_lambda <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, 
#                   correlation = corPagel(value=1, phy=anolis_phy, fixed=FALSE), method="ML")

pgls_lambda_beta_pcrit12 <- gls(slope~pcrit.r.kpa, 
                                data=beta_pcrit_12c_df_rn,
                                correlation=corPagel(value=1, 
                                                     phy=mandic_phy, 
                                                     fixed=FALSE),
                                method="ML")

## Compare models using AIC (lower is beter)
AIC(ols_beta_pcrit12)
AIC(pgls_bm_beta_pcrit12)
AIC(pgls_lambda_beta_pcrit12)

## Get the AIC weights
## Create a named vector
aic_all <- c(AIC(ols_beta_pcrit12), 
             AIC(pgls_bm_beta_pcrit12), 
             AIC(pgls_lambda_beta_pcrit12))
names(aic_all) <- c("ols", "bm", "lambda")
aicw(aic_all)
## More weight on bm BUT delta AIC < 2 for all, 
## so no real difference bw models

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(ols_beta_pcrit12, pgls_bm_beta_pcrit12)

plot(aov(slope ~ pcrit.r.kpa, data = beta_pcrit_12c_df_rn))
anova(ols_beta_pcrit12)
summary.aov(with(ols_beta_pcrit12))

## Plot of Beta Pcrit vs Pcrit at 12C, with ols linear model and CI

ggplot(pcrit_temp_slope_pcrit_12c_analysis_df,
       aes(x = pcrit.r.kpa,
           y = slope)) +
  geom_point(size = 7) +
  stat_smooth(method = "lm") +
  labs(x = expression(paste("P"["crit"]," at 12",degree,C," (kPa)")),
       y = expression(paste(beta["P"["crit"]]," (kPa ",degree,C^-1,")"))) +
  ggplot2::theme(axis.text.y = element_text(size = 20),
                 axis.title.y = (element_text(size = 26, margin = margin(t = 0, r = 20, b = 0, l = 0))),
                 axis.text.x = element_text(size = 20),
                 axis.title.x = element_text(size = 26, margin = margin(t = 20, r = 0, b = 0, l = 0)))

## Without Enorphrys bison
ggplot(pcrit_temp_slope_pcrit_12c_analysis_df[pcrit_temp_slope_pcrit_12c_analysis_df$species!="Enophrys_bison",],
       aes(x = pcrit.r.kpa,
           y = slope)) +
  geom_point(size = 7) +
  stat_smooth(method = "lm") +
  labs(x = expression(paste("P"["crit"]," at 12",degree,C," (kPa)")),
       y = expression(paste(beta["P"["crit"]]," (kPa ",degree,C^-1,")"))) +
  ggplot2::theme(axis.text.y = element_text(size = 18),
                 axis.title.y = (element_text(size = 26, margin = margin(t = 0, r = 20, b = 0, l = 0))),
                 axis.text.x = element_text(size = 18),
                 axis.title.x = element_text(size = 26, margin = margin(t = 20, r = 0, b = 0, l = 0)))

######
######
######
######
######
###       Beta Pcrit versus CTmax
######
######
######
######
######

pcrit_slope_ctmax_df <- full_join(pcrit_temp_models_df, ct_max_data_spps_summary, by = "species")
pcrit_slope_ctmax_df <- pcrit_slope_ctmax_df %>%
  ungroup() %>%
  filter(species != "Blepsias_cirrhosus")

pcrit_slope_ctmax_df_rn <- column_to_rownames(pcrit_slope_ctmax_df, var = "species")

mandic_phy_ctmax <- drop.tip(mandic_phy, "Blepsias_cirrhosus")

## Regular OLS
#ols <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, method="ML")
ols_beta_ctmax <- gls(slope~ct_max_avg, 
                        data=pcrit_slope_ctmax_df_rn, 
                        method="ML")

## PGLS assuming a Brownian correlation
#pgls_bm <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, correlation = corBrownian(phy=anolis_phy), method="ML")
pgls_bm_beta_ctmax <- gls(slope~ct_max_avg, 
                            data=pcrit_slope_ctmax_df_rn, 
                            correlation=corBrownian(phy=mandic_phy_ctmax), 
                            method="ML")

## PGLS assuming a Pagel's lambda correlation
#pgls_lambda <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, 
#                   correlation = corPagel(value=1, phy=anolis_phy, fixed=FALSE), method="ML")

pgls_lambda_beta_ctmax <- gls(slope~ct_max_avg,
                              data=pcrit_slope_ctmax_df_rn,
                              correlation=corPagel(value=1, 
                                                   phy=mandic_phy_ctmax,
                                                   fixed=FALSE),
                              method="ML")

## PGLS assuming a Blomberg correlation
#pgls_blomberg <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, 
#                   correlation = corPagel(value=1, phy=anolis_phy, fixed=FALSE), method="ML")

pgls_blomberg_beta_ctmax <- gls(slope~ct_max_avg,
                              data=pcrit_slope_ctmax_df_rn,
                              correlation=corBlomberg(value=1,
                                                      phy=mandic_phy_ctmax,
                                                      fixed=FALSE),
                              method="ML")


## Compare models using AIC (lower is beter)
AIC(ols_beta_ctmax)
AIC(pgls_bm_beta_ctmax)
#AIC(pgls_lambda_beta_ctmax) Model did not converge!

## Get the AIC weights
## Create a named vector
aic_all <- c(AIC(ols_beta_ctmax), 
             AIC(pgls_bm_beta_ctmax))
names(aic_all) <- c("ols", "bm")
aicw(aic_all)
## More weight on OLS! 

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(ols_beta_ctmax, pgls_bm_beta_ctmax)

anova(ols_beta_ctmax)

## Plotting Beta Pcrit vs CTmax
ggplot(pcrit_slope_ctmax_df,
       aes(x = ct_max_avg,
           y = slope)) +
  geom_point(size = 7) +
  stat_smooth(method = "lm") +
  labs(x = expression(paste("CT"["max"]," (",degree,C,")")),
       y = expression(paste(beta["P"["crit"]]," (kPa ",degree,C^-1,")"))) +
  ggplot2::theme(axis.text.y = element_text(size = 20),
                 axis.title.y = (element_text(size = 26, margin = margin(t = 0, r = 20, b = 0, l = 0))),
                 axis.text.x = element_text(size = 20),
                 axis.title.x = element_text(size = 26, margin = margin(t = 20, r = 0, b = 0, l = 0)))

## Plotting Beta Pcrit vs CTmax WITHOUT BUFFALO SCULPIN (Enophrys bison)
ggplot(pcrit_slope_ctmax_df,
       aes(x = ct_max_avg,
           y = slope)) +
  geom_point(size = 7) +
  stat_smooth(method = "lm") +
  labs(x = expression(paste("CT"["max"]," (",degree,C,")")),
       y = expression(paste(beta["P"["crit"]]," (kPa ",degree,C^-1,")"))) +
  ggplot2::theme(axis.text.y = element_text(size = 20),
                 axis.title.y = (element_text(size = 26, margin = margin(t = 0, r = 20, b = 0, l = 0))),
                 axis.text.x = element_text(size = 20),
                 axis.title.x = element_text(size = 26, margin = margin(t = 20, r = 0, b = 0, l = 0)))

















# Using code from Melissa Guzman!

# Can also do a phylogenetic ANCOVA

m_ancova_int_mass_pcrit <- lm(pcrit.r ~ temp * species + mass.g, data = pcrit_smr_data_summary)

Anova(m_ancova_int_mass_pcrit)

m_ancova <- lm(pcrit.avg ~ temp + species, data = pcrit_smr_data_spps_summary)
m_ancova_int <- lm(pcrit.avg ~ temp * species, data = pcrit_smr_data_spps_summary) # interaction included
m_ancova_int_raw <- lm(pcrit.r ~ temp * species, data = pcrit_smr_data_summary)

m_ancova_int_raw_smr <- lm(log10(smr.raw) ~ temp * species + log10(mass.g), data = pcrit_smr_data_summary)
Anova(m_ancova_int_raw_smr)

library(car)
Anova(m_ancova_int, type = "III")
summary(m_ancova_int)
summary(m_ancova)

Anova(m_ancova_int_raw, type = "III")


library(broom) ## Great package that takes model fit coefficients and puts them in a data frame!

# For Pcrit models: NO MASS ADJUSTMENT NEEDED!
fitted_pcrit_models <- pcrit_smr_data_summary %>% 
  group_by(species) %>%
  do(model = lm(pcrit.r ~ temp, data = .))

model_pcrit_df <- fitted_pcrit_models %>% tidy(model)

model_pcrit_df_spread <- model_pcrit_df %>%
  dplyr::select(term, estimate) %>%
  spread(term, estimate)
model_pcrit_df_spread ## THIS has the SLOPES of the lm(pcrit ~ temp) FOR EACH SPPS

# For SMR models: MASS ADJUSTMENT WITHIN SPECIES?!?
fitted_smr_massadj_models <- pcrit_smr_data_summary %>%
  group_by(species) %>%
  do(model = lm(log10(smr.raw)~log10(mass.g)+temp, data = .))

model_smr_massadj_df <- fitted_smr_massadj_models %>% tidy(model)

model_smr_massadj_df_spread <- model_smr_massadj_df %>%
  dplyr::select(term, estimate) %>%
  spread(term, estimate)
model_smr_massadj_df_spread ## THIS has the SLOPES of the lm(pcrit ~ temp) FOR EACH SPPS


##########################################################################################
library(phytools)
## datasets for phenograms
pcrit_slopes <- model_pcrit_df_spread %>%  # 9 species for Pcrit vs temp slopes
  dplyr::select(species, temp)
smr_slopes <- model_smr_massadj_df_spread %>%  # 9 species for SMR Vs temp slopes
  dplyr::select(species, temp)
ct_max_data_spps_summary # 8 species for CTmax

#############
# # # Data for Pcrit-temp slope phenogram
#############
# Using ct_max_df, no q10 values or anything besides ct_max data
#rownames(ct_max_df_phylo)<-c() ## For some stupid reason I have to run this
pcrit_slopes_pheno_rn <- column_to_rownames(pcrit_slopes, var = "species")
pcrit_slopes_pheno_matrix <- as.matrix(pcrit_slopes_pheno_rn)[,1]
#############
# # # END Data for Pcrit-temp slope phenogram
#############

#par()
#orig_par <- par()
#par(cex.lab = 2)
#par(orig_par)

## Phenogram of Pcrit ~ temp slopes on phylogeny
pheno_pcrit_slopes <- phenogram(mandic_phy, pcrit_slopes_pheno_matrix,
                         ylab=expression(paste("P"["crit"],"~temperature: slope")),
                         xlab="Relative time",
                         xaxt='n',
                         cex.lab = 2)
contMap(mandic_phy, pcrit_slopes_pheno_matrix)

#############
# # # Data for  SMR ~ temp slopes on phylogeny
#############
# Using ct_max_df, no q10 values or anything besides ct_max data
#rownames(ct_max_df_phylo)<-c() ## For some stupid reason I have to run this
smr_slopes_pheno_rn <- column_to_rownames(smr_slopes, var = "species")
smr_slopes_pheno_matrix <- as.matrix(smr_slopes_pheno_rn)[,1]
#############
# # # END Data for SMR ~ temp slopes phenogram
#############

## Phenogram of SMR ~ temp slopes on phylogeny
pheno_smr_slopes <- phenogram(mandic_phy, smr_slopes_pheno_matrix,
                         ylab=expression(paste("SMR ~ temperature: slope")),
                         xlab="Relative time",
                         xaxt='n')
contMap(mandic_phy, smr_slopes_pheno_matrix)

#############
# # # Data for ct max phenogram
#############
# Using ct_max_df, no q10 values or anything besides ct_max data
#rownames(ct_max_df_phylo)<-c() ## For some stupid reason I have to run this
ctmax_pheno_rn <- column_to_rownames(ct_max_data_spps_summary, var = "species")
ctmax_pheno_matrix <- as.matrix(ctmax_pheno_rn)[,1]
#############
# # # END Data for ct max phenogram
#############
mandic_phy_ctmax <- drop.tip(mandic_phy, "Blepsias_cirrhosus")
## Phenogram of CTmax on phylogeny
pheno_ctmax <- phenogram(mandic_phy_ctmax, ctmax_pheno_matrix,
                         ylab=expression(paste("CTmax (",degree,C,")")),
                         xlab="Relative time",
                         xaxt='n')
contMap(mandic_phy_ctmax, ctmax_pheno_matrix)

#############
## Pcrit at 12 versus Pcrit-temp slope
pcrit12_vs_pcritslope <- full_join(model_pcrit_df_spread, pcrit_smr_data_spps_summary[pcrit_smr_data_spps_summary$temp==12,],by = "species")
  
ggplot(pcrit12_vs_pcritslope, aes(x=pcrit.avg, y=temp.x)) +
  geom_point(size = 4) +
  stat_smooth(method="lm", se=TRUE) +
  labs(x = expression(paste("P"["crit"]," at 12",degree,C," (torr)")),
       y = expression(paste("P"["crit"]," ~ temperature: slope (torr ",degree,C^-1,")"))) +
  ggplot2::theme(axis.text.y = element_text(size = 18),
                 axis.title.y = (element_text(size = 26, margin = margin(t = 0, r = 20, b = 0, l = 0))),
                 axis.text.x = element_text(size = 18),
                 axis.title.x = element_text(size = 26, margin = margin(t = 20, r = 0, b = 0, l = 0)))

##################################################################################################

##      CTmax vs Pcrit~temp slope

ctmax_vs_pcritslope <- full_join(pcrit12_vs_pcritslope, ct_max_data_spps_summary, by = "species") %>%
  filter(!is.na(ct_max_avg))

ggplot(ctmax_vs_pcritslope, aes(x=ct_max_avg, y=temp.x)) +
  geom_point(size = 4) +
  stat_smooth(method="lm", se=TRUE) +
  labs(x = expression(paste("CT"["max"]," (",degree,C,")")),
       y = expression(paste("P"["crit"]," ~ temperature: slope (torr ",degree,C^-1,")"))) +
  ggplot2::theme(axis.text.y = element_text(size = 18),
                 axis.title.y = (element_text(size = 26, margin = margin(t = 0, r = 20, b = 0, l = 0))),
                 axis.text.x = element_text(size = 18),
                 axis.title.x = element_text(size = 26, margin = margin(t = 20, r = 0, b = 0, l = 0)))

#model_pcritSlopes_ctmax <- lm(temp.x~ct_max_avg, data = ctmax_vs_pcritslope)
#Anova(model_pcritSlopes_ctmax, type = "III")

## Load the anolis tree and dataset
#anolis_phy <- read.tree("anolis.tre.R")
mandic_phy_ctmax
## Tip: Setting `row.names=1` assigns the rownames to be equal to the first column (in this case the species names)
#anolis_dat <- read.csv("anolis.csv", row.names=1)
ctmax_vs_pcritslope_rn <- column_to_rownames(ctmax_vs_pcritslope, var="species")
## Using gls to test for an association of body size with hind limb (dumb, i know but bear with me)

## Regular OLS
#ols <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, method="ML")
ols_pslopes_ctmax <- gls(temp.x~ct_max_avg, data=ctmax_vs_pcritslope_rn, method="ML")
  
## PGLS assuming a Brownian correlation
#pgls_bm <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, correlation = corBrownian(phy=anolis_phy), method="ML")
pgls_bm_pslopes_ctmax <- gls(temp.x~ct_max_avg, data=ctmax_vs_pcritslope_rn, correlation=corBrownian(phy=mandic_phy_ctmax), method="ML")

## PGLS assuming a Pagel's lambda correlation
#pgls_lambda <- gls(AVG_hl ~ AVG_SVL, data=anolis_dat, 
#                   correlation = corPagel(value=1, phy=anolis_phy, fixed=FALSE), method="ML")

pgls_lambda_pslopes_ctmax <- gls(temp.x~ct_max_avg, data=ctmax_vs_pcritslope_rn,
                                 correlation=corPagel(value=1, phy=mandic_phy_ctmax, fixed=FALSE),
                                 method="ML")

## Compare models using AIC (lower is beter)
AIC(ols_pslopes_ctmax)
AIC(pgls_bm_pslopes_ctmax)
# AIC(pgls_lambda) DID NOT CONVERGE

## Using a likelihood ratio test to evaluate whether pgls is justified
anova(pgls_bm_pslopes_ctmax, ols_pslopes_ctmax)

## Is lambda significantly better than Brownian motion
anova(ols_pslopes_ctmax, type="marginal")
summary(ols_pslopes_ctmax)

m_linear_pslopes_ctmax <- lm(temp.x~ct_max_avg, data=ctmax_vs_pcritslope_rn)
Anova(m_linear_pslopes_ctmax, type="III")

## Bop









cor(model_df_spread$`(Intercept)`, model_df_spread$temp) # Testing correlation between temp slope and intercept
## I'm not sure if the intercept of the lm: pcrit ~ temp is a fair representation of my hypothesis...

plot(model_df_spread$`(Intercept)`, model_df_spread$temp)
text(model_df_spread$`(Intercept)`, model_df_spread$temp, labels=model_df_spread$species, cex=0.7) #, cex= 0.7


##################################
####        PCRIT VS TEMP

## Visualize WITH interaction using INDIVIDUAL fish data, group by species
Anova(m_ancova_int_raw, type = "III")
visreg(m_ancova_int_raw, xvar = "temp", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))
visreg(m_ancova_int_raw, xvar = "temp", by = "species", 
       points.par = list(cex = 1.1, col = "red"))
visreg(m_ancova_int_raw, xvar = "temp", by = "species", whitespace = 0.4, overlay = TRUE, 
       band = FALSE, points.par = list(cex = 1.1))

##################################
####        SMR VS TEMP

## Visualize WITH interaction using INDIVIDUAL fish data, group by species
Anova(m_ancova_int_raw_smr, type = "III")
visreg(m_ancova_int_raw_smr, xvar = "temp", whitespace = 0.4, 
       points.par = list(cex = 1.1, col = "red"))
visreg(m_ancova_int_raw_smr, xvar = "temp", by = "species", 
       points.par = list(cex = 1.1, col = "red"))
visreg(m_ancova_int_raw_smr, xvar = "temp", by = "species", whitespace = 0.4, overlay = TRUE, 
       band = FALSE, points.par = list(cex = 1.1))





#### PLOTS PLOTS PLOTS

pcrit_temp_plot_global_line <- ggplot(pcrit_smr_data_summary, aes(x = temp, y = pcrit.r)) +
  geom_jitter(size = 3, width = 0.2) +
  stat_smooth(method = "lm", se = TRUE) +
  labs(x = expression(paste("Temperature (",degree,C,")")),
       y = expression(paste("P"["crit"]," (torr)")))


pcrit_temp_plot_spps_line <- ggplot(pcrit_smr_data_summary, aes(x = temp, y = pcrit.r, colour = species, group = species)) +
  geom_jitter(size = 3, width = 0.2) +
  stat_smooth(method = "lm", se = FALSE) +
  labs(x = expression(paste("Temperature (",degree,C,")")),
       y = expression(paste("P"["crit"]," (torr)")))

pcrit_temp_plot <- ggplot(pcrit_smr_data_spps_summary, aes(x = temp, y = pcrit.avg, colour = species)) +
  geom_jitter(size = 3, width = 0.2) +
  stat_smooth(method = "lm", se = FALSE) +
  labs(x = expression(paste("Temperature (",degree,C,")")),
       y = expression(paste("P"["crit"]," (torr)")))














######################
######################
######################
######################
######################            EFFECT OF TEMPERATURE ON HYPOXIA TOLERANCE:
######################
######################                Using 9 species/3 temperatures as (N=27) data points
######################
######################
######################
######################

geo=get(data(geospiza))
dat=geo$dat
d1=dat[,1]
grp<-as.factor(c(rep(0, 7), rep(1, 6)))
names(grp)=rownames(dat)
## MANOVA
x=aov.phylo(dat~grp, geo$phy, nsim=50, test="Wilks")
print(attributes(x)$summary) # summary table
## ANOVA
x1=aov.phylo(d1~grp, geo$phy, nsim=50)















## model:
## Y = mu + Beta*X + a + s + e
# mu = intercept
# Beta = slope
# a = phylogenetic effect << This let's you do a regression on the full set of ALL species data
# s = "species-specific effect" << This "accounts for the variability that has been caused by species-specific effects"
# e = residual error

inv_phylo <- inverseA(mandic_phy, nodes = "TIPS", scale = TRUE)
#inv_phyloAsReml <- sm2asreml(inv_phylo$Ainv, inv_phylo$node.names)

prior <- list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

m_allSpps_simple <- MCMCglmm(pcrit.avg~temp,random = ~species,
                             family="gaussian",ginverse = list(species=inv_phylo$Ainv),
                             prior=prior,data=pcrit_smr_data_spps_summary,
                             nitt = 5000000,burnin = 1000,thin=500)

pcrit_smr_data_spps_summary


## End of work in progress



## Example graph
scaling_olma_mo2_plot <- ggplot(plot_check_olma_df,
                                aes(x = ind_mass_g, y = ind_avg_umol_hr)) +
  geom_point(size = 3) +
  stat_function(fun = fn_olma, colour = "red", lwd = 1) +
  stat_function(fun = fn_olma_mean, colour = "blue", lwd = 1) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        axis.text.y = element_text(size = 24),
        axis.title.y = (element_text(size = 28, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 28, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Body mass (g)")),
       y = expression(paste(dot(M),"o"[2]," (",mu,"mol O"[2]," hr"^-1,")"))) +
  annotate ("segment", x = 4.95, xend = 5.2, y = 40.176, yend = 40.176, size = 1.25, colour = "black", arrow = arrow ())

