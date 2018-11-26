##################################
# Load required packages
##################################
library(tidyverse)
library(broom)
library(lubridate)
library(mclust)
library(shape)
#library(StreamMetabolism)
#library(fishMO2)

##################################
##################################
#
# Import and organize Pcrit data
#
##################################
##################################
mo2_folder <- "C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/KenChu MMR hypoxia scaling Tidepool Sculpin study/Raw foxy csv files"
mo2_file_list <- list.files(path = mo2_folder)
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/KenChu MMR hypoxia scaling Tidepool Sculpin study/Raw foxy csv files")
## This works to read in all the data files
## From following website:
## https://stackoverflow.com/questions/26273555/inserting-file-names-as-column-values-in-a-data-frame?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
mo2_data_files <- lapply(mo2_file_list, read.table, header = TRUE, sep = ",")
for (i in 1:length(mo2_data_files)) {
  mo2_data_files[[i]] <- cbind(mo2_data_files[[i]], mo2_file_list[i])
}
mo2_data <- do.call("rbind", mo2_data_files)
colnames(mo2_data)[c(1,2,3,4,5,6)] <- c("foxy-date-time", "do", "tau",
                                  "foxy_temp", "foxy_air_pressure", "file")

#pcrit_data <- as_tibble(pcrit_data)
#pcrit_data

# Extract the date from filename
# See r4ds online e-book, Section 14.3: Matching patterns with regular expressions
#str_extract(pcrit_data$file, "(\\d|\\d\\d)[A-Z][a-z][a-z]\\d\\d\\d\\d")[1:10] # extracts date from filename
#str_extract(pcrit_data$file, "(\\d|\\d{2})[A-Z][a-z]{2}\\d{4}")[1:10] # same thing but shorter code

mo2_data$trial_date <- str_extract(mo2_data$file, "(\\d|\\d{2})[A-Z][a-z]{2}\\d{4}")

mo2_data$trial_date <- case_when(
  mo2_data$trial_date == str_extract(mo2_data$trial_date, "\\d{2}[A-Z][a-z]{2}\\d{4}") ~
    str_c(
      str_sub(mo2_data$trial_date, start = 1, end = 2),
      "-",
      str_sub(mo2_data$trial_date, start = 3, end = 5),
      "-",
      str_sub(mo2_data$trial_date, start = 6, end = 9)
    ),
  mo2_data$trial_date == str_extract(mo2_data$trial_date, "\\d[A-Z][a-z]{2}\\d{4}") ~
    str_c(
      str_sub(mo2_data$trial_date, start = 1, end = 1),
      "-",
      str_sub(mo2_data$trial_date, start = 2, end = 4),
      "-",
      str_sub(mo2_data$trial_date, start = 5, end = 8)
    )
)

# Extract the probe from filename
mo2_data$probe <- str_extract(mo2_data$file, "NFB00[0-9]{2}")
head(mo2_data)

# Extract the spps from filename          #### WHERE I LEFT OFF!!
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
