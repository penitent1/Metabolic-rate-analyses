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
# Import data: POST Labchart
##################################
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Roughback sculpin prelim analyses/Post-Labchart processed data")
md_lc <- read_csv("Post_labchart_processed_mo2_data.csv")
md_lc_nest <- md_lc %>%
  mutate(mo2_raw = umol_per_sec*3600*-1,
         mo2_ms = mo2_raw/mass_g) %>%
  group_by(date, finclip, expt_period) %>%
  nest()

## Take a look at the dataframe
md_lc_nest

## Plot out mo2 vs time RMR data
md_lc_nest$data[[1]] %>% ggplot() +
  geom_point(aes(x = time_min, y = mo2_ms))

## Calc "RMR" using fishMO2
rmr <- md_lc_nest$data[[1]] %>% as.data.frame(.)
calcSMR(rmr$mo2_ms)

## Plot out Pcrit data
md_lc_nest$data[[2]] %>% ggplot() +
  geom_point(aes(x = `po2 torr`, y = mo2_ms), size = 2) #+
  #geom_line(aes(x = `po2 torr`, y = mo2_ms))

## Calc Pcrit using fishMO2
pcrit <- md_lc_nest$data[[2]] %>% rename(DO =`corr % sat`, MO2 = mo2_ms) %>%
  dplyr::select(MO2,DO) %>% as.data.frame(.)
calcO2crit(pcrit, SMR = 1.64)
plotO2crit(calcO2crit(pcrit, SMR = 1))

## Plot out Pcrit data: 10 min interval
md_lc_nest$data[[3]] %>% ggplot() +
  geom_point(aes(x = `po2 torr`, y = mo2_ms), size = 2) #+
#geom_line(aes(x = `po2 torr`, y = mo2_ms))

## Calc Pcrit using fishMO2: 10 min interval
pcrit_10min <- md_lc_nest$data[[3]] %>% rename(DO =`corr % sat`, MO2 = mo2_ms) %>%
  dplyr::select(MO2,DO) %>% as.data.frame(.)
calcO2crit(pcrit_10min, SMR = 1.64)
plotO2crit(calcO2crit(pcrit_10min, SMR = 1.64))

## Plot out Pcrit data: 2 min interval
md_lc_nest$data[[4]] %>% ggplot() +
  geom_point(aes(x = `po2 torr`, y = mo2_ms), size = 2) +
  geom_line(aes(x = `po2 torr`, y = mo2_ms))

## Calc Pcrit using fishMO2: 2 min interval
pcrit_2min <- md_lc_nest$data[[4]] %>% rename(DO =`corr % sat`, MO2 = mo2_ms) %>%
  dplyr::select(MO2,DO) %>% as.data.frame(.)
calcO2crit(pcrit_2min, SMR = 1.64)
plotO2crit(calcO2crit(pcrit_2min, SMR = 1.64))
