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
library(car)

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
ct_max_data_spps_summary

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
         tidy_pcrit_lm = spps_scaling_mod_pcrit %>% purrr::map(broom::tidy),
         p_val_smr_aov = tidy_smr_aov %>% purrr::map_dbl(c(5,2)),
         p_val_pcrit_aov = tidy_pcrit_aov %>% purrr::map_dbl(c(5,2)))

