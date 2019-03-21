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
#setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit")
setwd("C:/Users/derek/OneDrive/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit")
md_lc <- read_csv("Ken_Tidepool_CollatedMO2Data_MMR-SMR-Pcrit_PostLabChart.csv") # Time not importing properly from weird foxy format, use seq()
md_lc <- md_lc %>%
  mutate(mo2_raw = umol_o2_per_sec*3600*-1,
         time = as.double(time)) %>%
  group_by(date, finclip_id, expt_period) %>%
  nest()

## Take a look at the dataframe
md_lc

## Visualize MMR-SMR data
md_lc %>% dplyr::select(-c(date)) %>% filter(expt_period == "mmr_smr") %>% unnest() %>%
  ggplot(aes(x = time, y = mo2_raw)) +
  geom_point() +
  geom_line() +
  facet_wrap("finclip_id")

## Remove first "middle_stumpy" finclip mmr_smr trial
notes_check <- md_lc %>% filter(finclip_id != "middle_stumpy") %>% mutate(notes = data %>% purrr::map("notes"))

md_lc %>% filter(!(finclip_id == "middle_stumpy" & date == "24-Oct-18")) %>% filter(expt_period == "mmr_smr") %>%
  unnest() %>%
  ggplot(aes(x = time, y = mo2_raw)) +
  geom_point() +
  geom_line() +
  facet_wrap("finclip_id") ## Looks good!

## Calculate MMR and SMR!
mmr_smr_md_lc <- md_lc %>% filter(!(finclip_id == "middle_stumpy" & date == "24-Oct-18")) %>% filter(expt_period == "mmr_smr") %>%
  mutate(MMR = data %>% purrr::map_dbl(~ max(.$mo2_raw)),
         mass_g = data %>% purrr::map_dbl(~ mean(.$mass_g)),
         smr_object = data %>% purrr::map(~ calcSMR(.$mo2_raw)),
         SMR_mlnd = smr_object %>% purrr::map_dbl("mlnd"))

# Does body size affect MMR?
mmr_smr_md_lc %>% lm(log(.$MMR) ~ log(.$mass_g), data = .) %>% anova() # Yes
mmr_smr_md_lc %>% lm(log(.$MMR) ~ log(.$mass_g), data = .) # Parameter estimates from regression

# Does body size affect SMR?
mmr_smr_md_lc %>% lm(log(.$SMR_mlnd) ~ log(.$mass_g), data = .) %>% anova() # Yes
mmr_smr_md_lc %>% lm(log(.$SMR_mlnd) ~ log(.$mass_g), data = .) # Parameter estimates from the linear regression

# Visualizing relationship: log(SMR and MMR) ~ log(body mass)
mmr_smr_md_lc %>% ggplot() +
  geom_point(aes(x = log(mass_g), y = log(SMR_mlnd)), size = 4, colour = "red") +
  stat_smooth(aes(x = log(mass_g), y = log(SMR_mlnd)), method = "lm", colour = "black") +
  geom_point(aes(x = log(mass_g), y = log(MMR)), size = 4, colour = "blue") +
  stat_smooth(aes(x = log(mass_g), y = log(MMR)), method = "lm", colour = "black") +
  scale_x_continuous(name = expression(paste("ln(mass) (ln(g))")),
                     limits = c(0,3),
                     breaks = seq(0,3,1)) +
  scale_y_continuous(name = expression(paste("ln(",dot(M),"o"[2],") (ln(",mu,"mol h"^-1,"))")),
                     limits = c(0,5),
                     breaks = seq(0,5,1))
