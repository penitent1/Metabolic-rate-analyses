getwd()
data_path <- "C:/Users/derek/OneDrive/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit/Raw foxy Pcrit files respR analyses"
files <- dir(data_path, pattern = "*.csv")
## Januar's function for parsing individual datafiles from NEOFOX (.csv files)
parse_neofox <- function(x) {
  raw <- fread(x)
  dtime <- ymd_hms(raw$Time)
  ntime <- as.numeric(difftime(dtime, dtime[1], units = "secs"))
  out <- data.table("Elapsed (s)" = ntime, raw)
  return(out)
}
pcrit_data <- tibble(filename = files) %>%
  mutate(file_contents = filename %>% purrr::map(~ parse_neofox(file.path(data_path, .))),
         inspect_pcrits = file_contents %>% purrr::map(~ inspect(.)))

pcrit_data$inspect_pcrits[[3]]
  

library(data.table)
library(lubridate)

#parse_neofox <- function(path) {
#  raw <- fread(path)
#  dtime <- ymd_hms(raw$Time)
#  ntime <- as.numeric(difftime(dtime, dtime[1], units = "secs"))
#  out <- data.table("Elapsed (s)" = ntime, raw)
#  return(out)
#}

parse_neofox <- function(x) {
  raw <- fread(x)
  dtime <- ymd_hms(raw$Time)
  ntime <- as.numeric(difftime(dtime, dtime[1], units = "secs"))
  out <- data.table("Elapsed (s)" = ntime, raw)
  return(out)
}

parse_neofox("./Ken_Tidepool_BottomFinClip_Pcrit_12C_NFB0013_05Oct2018.csv")
