##################################
# Load required packages
##################################
library(tidyverse)
library(broom)
library(lubridate)
library(mclust)
library(shape)
library(StreamMetabolism)
library(fishMO2)

##################################
# Import data: Pre Labchart
##################################
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit")
md <- read_csv("Ken_Tidepool_CollatedMO2Data_MMR-SMR-Pcrit.csv") # Time not importing properly from weird foxy format, use seq()
md <- md %>%
  select(species,finclip_id,probe,resp,date_start,sampling_rate_sec,
         expt_po2_type,expt_period,Oxygen,`Air Pressure_GovCan`) %>%
  group_by(finclip_id,probe,date_start,expt_po2_type,expt_period) %>%
  mutate(time_sec = 1:n(),
         time_hr = time_sec/3600) %>%
  filter(row_number() %/% 15 == 0) %>%
  nest()

ggplot(md$data[[1]], aes(x = time, y = Oxygen)) +
  geom_point()

probe10_data <- tibble(pre = c(0,100), post = c(2.2,101.4))
lm(post ~ pre, data = probe10_data)

##################################
# Import data: POST Labchart
##################################
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit")
md_lc <- read_csv("Ken_Tidepool_CollatedMO2Data_MMR-SMR-Pcrit_PostLabChart.csv") # Time not importing properly from weird foxy format, use seq()
md_lc <- md_lc %>%
  mutate(mo2_raw = umol_o2_per_sec*3600*-1) %>%
  group_by(finclip_id) %>%
  filter(finclip_id == "top",
         probe == "NFB0008") # Update per fish

# Add "time order" bc I did not add a time column for labchart :)
md_lc_mmr_smr <- md_lc %>%
  filter(expt_period == "mmr_smr") %>%
  mutate(time = as.double(Time))
  #mutate(time_order = 1:n())

# Visual estimate of when to cut off MMR recovery
ggplot(md_lc_mmr_smr, aes(x = time, y = mo2_raw)) +
  geom_point()

# Max MO2
max(md_lc_mmr_smr$mo2_raw)

# For SMR estimate
md_lc_smr <- md_lc_mmr_smr %>%
  filter(time > 800)

smr <- calcSMR(md_lc_smr$mo2_raw)
smr_best <- ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)

# For Pcrit estimate
md_lc_pcrit_pretibble <- md_lc %>% filter(expt_period != "mmr_at_pcrit")
md_lc_pcrit <- tibble(MO2 = md_lc_pcrit_pretibble$mo2_raw, DO = md_lc_pcrit_pretibble$do_percent_sat)
plotO2crit(calcO2crit(md_lc_pcrit, SMR = 4.32))

ggplot() +
  geom_point(md_lc_pcrit, aes(x = DO, y = MO2)) +
  geom_point(md_lc_mmr_pcrit_summ, aes(x = mean_do, y = mean_mo2))

# For MMR at Pcrit
md_lc_mmr_pcrit <- md_lc %>%
  filter(expt_period == "mmr_at_pcrit")
max(md_lc_mmr_pcrit$mo2_raw)

ggplot(md_lc_mmr_pcrit, aes(x = do_percent_sat, y = mo2_raw)) +
  geom_point() +
  geom_point(aes(x = mean(do_percent_sat), y = mean(mo2_raw)), colour = "red")

md_lc_mmr_pcrit_summ <-  md_lc_mmr_pcrit %>% 
  select(do_percent_sat, mo2_raw) %>% 
  dplyr::summarise(mean_do = mean(do_percent_sat),
                   mean_mo2 = mean(mo2_raw),
                   max_mo2 = max(mo2_raw),
                   max_do = 12.2)
