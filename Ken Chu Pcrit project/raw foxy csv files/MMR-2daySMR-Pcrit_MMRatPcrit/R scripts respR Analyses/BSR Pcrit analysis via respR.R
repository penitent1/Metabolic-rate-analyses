library(tidyverse)
library(respR)
library(data.table)
library(lubridate)

## Define filepath to folder containing pcrit data files (.csv from NEOFOX software)
data_path <- "C:/Users/derek/OneDrive/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit/Raw foxy Pcrit files respR analyses"

## Grab file names and add .csv to the end of each to generate a character vector of filenames
files <- dir(data_path, pattern = "*.csv")

## Januar's function for parsing individual datafiles from NEOFOX (.csv files)
parse_neofox <- function(x) {
  raw <- fread(x)
  dtime <- ymd_hms(raw$Time)
  ntime <- as.numeric(difftime(dtime, dtime[1], units = "secs"))
  out <- data.table("Elapsed (s)" = ntime, raw)
  return(out)
}

## Read data from "files" into a tibble
pcrit_data <- tibble(filename = files) %>%
  mutate(file_contents = filename %>% purrr::map(~ parse_neofox(file.path(data_path, .)))) %>%
  unnest() %>%
  dplyr::select(filename, `Elapsed (s)`, Oxygen) %>%
  group_by(filename) %>%
  nest()

## Inspect data for each fish' pcrit trial and 
pcrit_data$data[[1]]
inspect(pcrit_data$data[[1]])
pcrit(pcrit_data$data[[1]])

pcrit_data$data[[2]]
inspect(pcrit_data$data[[2]])
pcrit(pcrit_data$data[[2]])

pcrit_data$data[[3]]
parse_neofox("C:/Users/derek/OneDrive/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit/Raw foxy Pcrit files respR analyses/Ken_Tidepool_InTopInBottomFinClip_Pcrit_12C_NFB0008_2Nov2018.csv")
inspect(pcrit_data$data[[3]])
pcrit(pcrit_data$data[[3]])

pcrit_data$data[[4]]
inspect(pcrit_data$data[[4]])
pcrit(pcrit_data$data[[4]])

pcrit_data$data[[5]]
inspect(pcrit_data$data[[5]])
pcrit(pcrit_data$data[[5]])

pcrit_data$data[[6]]
inspect(pcrit_data$data[[6]])
pcrit(pcrit_data$data[[6]])
