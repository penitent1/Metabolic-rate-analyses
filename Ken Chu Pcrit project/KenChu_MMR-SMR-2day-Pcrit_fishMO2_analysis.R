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
  mutate(mo2_raw = umol_o2_per_sec*3600*-1,
         time = as.double(time)) %>%
  group_by(date, finclip_id, expt_period) %>%
  nest()

## Take a look at the dataframe
md_lc

ggplot(md_lc$data[[13]], aes(x = time, y = mo2_raw/mass_g))+geom_point()+geom_line()
ggplot(md_lc$data[[14]], aes(x = do_percent_sat, y = mo2_raw/mass_g))+geom_point(size = 3)

smr_data <- md_lc$data[[13]] %>% dplyr::select(do_percent_sat,po2_torr,mo2_raw,mass_g) %>% 
  mutate(MO2 = mo2_raw/mass_g, DO = do_percent_sat) %>% as.data.frame(.)
head(smr_data)

max(smr_data$MO2)

smr_est <- calcSMR(smr_data$MO2) # Top-bottom SMR = 1.13

ifelse(smr_est$CVmlnd > 5.4, smr_est$quant, smr_est$mlnd)

ggplot(md_lc$data[[13]], aes(x = time, y = mo2_raw/mass_g))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = 1.85)

pcrit_data <- md_lc$data[[14]] %>% dplyr::select(do_percent_sat,po2_torr,mo2_raw,mass_g) %>% 
  mutate(MO2 = mo2_raw/mass_g, DO = do_percent_sat) %>% as.data.frame(.)
head(pcrit_data)

#write_excel_csv(tibble(do = pcrit_data$DO, mo2 = pcrit_data$MO2),"C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/olma_pcrit_intop-inbottom.csv")

plotO2crit(calcO2crit(pcrit_data, SMR = 1.85))

smr_pcrit_data <- bind_rows(smr_data, pcrit_data)
head(smr_pcrit_data)

plotO2crit(calcO2crit(smr_pcrit_data, SMR = 1.85))

## MMR at pcrit

mmr_at_pcrit <- md_lc$data[[12]] %>% 
  mutate(mo2_ms = mo2_raw/mass_g) %>% 
  dplyr::select(mass_g, time, mo2_raw, mo2_ms)

ggplot(mmr_at_pcrit, aes(x = time, y = mo2_ms)) +
  geom_point(size = 3)

mmr_at_pcrit

mean(mmr_at_pcrit$mo2_ms)
