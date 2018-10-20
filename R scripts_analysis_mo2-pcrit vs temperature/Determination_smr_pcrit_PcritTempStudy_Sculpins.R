##################################
# Load required packages
##################################
library(tidyverse)
library(broom)
library(mclust)
library(shape)
library(StreamMetabolism)
library(fishMO2)

##################################
# Import MO2 data: SMR estimation
##################################
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Raw MO2 slope data")
smr_data <- read_csv("Spring2017_Sculpins_SMR_3day_ForR.csv")

##################################
# Organize data for SMR analysis
# and calculate SMR
##################################
smr_data_post5hr <- smr_data %>%
  group_by(probe,spps,date.day,date.month,date.year,trial.no,pcrit.type,temp) %>%
  filter(time.hrs > 5) %>%
  nest() %>%
  mutate(smr_ms_list = data %>% purrr::map(~ fishMO2::calcSMR(.$mo2.ms)),
         smr_ms_cv_mlnd = smr_ms_list %>% purrr::map_dbl("CVmlnd"),
         smr_ms_quantiles = smr_ms_list %>% purrr::map("quant"),
         smr_ms_smr_q2 = smr_ms_list %>% purrr::map_dbl(4),
         smr_ms_smr_mlnd = smr_ms_list %>% purrr::map_dbl("mlnd"),
         smr_ms = if_else(smr_ms_cv_mlnd > 5.4, smr_ms_smr_q2, smr_ms_smr_mlnd))

##################################
# Zoom in on ARLA, 12C
##################################
smr_data_post5hr %>%
  filter(spps == "arla",
         temp == 12) %>%
  dplyr::select(probe,spps,date.day,date.month,date.year,trial.no,data,smr_ms) %>%
  unnest() %>%
  filter(!is.na(smr_ms),
         trial.no == 1) %>%
  ggplot(aes(x = time.hrs, y = mo2.ms)) +
  geom_point() +
  facet_wrap(c("probe","date.day","date.month","date.year","trial.no"))

##################################
##################################
#
# Import and organize Pcrit data
#
##################################
##################################
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

# Extract the probe from filename
pcrit_data$probe <- str_extract(pcrit_data$file, "NFB00[0-9]{2}")
pcrit_data

# Extract the spps from filename
pcrit_data$spps <- str_extract(pcrit_data$file, "[A-Z]{4}")
pcrit_data$spps <- tolower(pcrit_data$spps) ## Convert spps names "to lower" case letters

# Extract the temperature from filename
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

pcrit_data

pcrit_data %>%
  group_by(trial_date,probe,species,temp_c) %>%
  filter(is.na(temp_c)) %>%
  nest()
