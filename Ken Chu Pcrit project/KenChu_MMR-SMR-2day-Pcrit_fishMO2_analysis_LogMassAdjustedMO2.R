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
  filter(!(finclip_id == "middle_stumpy" & date == "24-Oct-18")) %>%
  filter(!(finclip_id == "top" & date == "3-Oct-18")) %>%
  filter(!(finclip_id == "top" & date == "5-Oct-18")) %>%
  nest()

## Take a look at the dataframe
md_lc

## Visualize MMR-SMR data
md_lc %>% dplyr::select(-c(date)) %>% filter(expt_period == "mmr_smr") %>% unnest() %>%
  ggplot(aes(x = time, y = mo2_raw)) +
  geom_point() +
  geom_line() +
  facet_wrap("finclip_id")

## Calculate MMR and SMR!
mmr_smr_md_lc <- md_lc %>% 
  filter(expt_period == "mmr_smr") %>%
  mutate(MMR = data %>% purrr::map_dbl(~ max(.$mo2_raw)),
         mass_g = data %>% purrr::map_dbl(~ mean(.$mass_g)),
         smr_data = data %>% purrr::map(~ filter(., time > 1000)),
         smr_object = smr_data %>% purrr::map(~ calcSMR(.$mo2_raw)),
         SMR_mlnd = smr_object %>% purrr::map_dbl("mlnd"),
         SMR_quantiles = smr_object %>% purrr::map("quant"),
         SMR_q20 = SMR_quantiles %>% purrr::map_dbl(4))

mmr_smr_md_lc %>% dplyr::select(MMR,mass_g,SMR_mlnd,SMR_q20)

# Does body size affect MMR?
mmr_smr_md_lc %>% lm(log(.$MMR) ~ log(.$mass_g), data = .) %>% anova() # Yes
mmr_smr_md_lc %>% lm(log(.$MMR) ~ log(.$mass_g), data = .) # Parameter estimates from regression

# Does body size affect SMR, estimated using mlnd?
mmr_smr_md_lc %>% lm(log(.$SMR_mlnd) ~ log(.$mass_g), data = .) %>% anova() # Yes
mmr_smr_md_lc %>% lm(log(.$SMR_mlnd) ~ log(.$mass_g), data = .) # Parameter estimates from the linear regression

# Does body size affect SMR, estimated using the 20th percentile estimate?
mmr_smr_md_lc %>% lm(log(.$SMR_q20) ~ log(.$mass_g), data = .) %>% anova() # Yes
mmr_smr_md_lc %>% lm(log(.$SMR_q20) ~ log(.$mass_g), data = .) # Parameter estimates from the linear regression


# Visualizing relationship: log(SMR and MMR) ~ log(body mass)
mmr_smr_md_lc %>% ggplot() +
  geom_point(aes(x = log(mass_g), y = log(SMR_mlnd)), size = 4, colour = "red") +
  stat_smooth(aes(x = log(mass_g), y = log(SMR_mlnd)), method = "lm", colour = "black") +
  geom_point(aes(x = log(mass_g), y = log(MMR)), size = 4, colour = "blue") +
  stat_smooth(aes(x = log(mass_g), y = log(MMR)), method = "lm", colour = "black") +
  geom_point(aes(x = log(mass_g), y = log(SMR_q20)), size = 4, colour = "magenta") +
  stat_smooth(aes(x = log(mass_g), y = log(SMR_q20)), method = "lm", colour = "black") +
  scale_x_continuous(name = expression(paste("ln(mass) (ln(g))")),
                     limits = c(0,3),
                     breaks = seq(0,3,1)) +
  scale_y_continuous(name = expression(paste("ln(",dot(M),"o"[2],") (ln(",mu,"mol h"^-1,"))")),
                     limits = c(0,5),
                     breaks = seq(0,5,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)))

## Visualize Pcrit data
md_lc %>% dplyr::select(-c(date)) %>% filter(expt_period == "pcrit") %>% unnest() %>%
  ggplot(aes(x = po2_torr, y = mo2_raw)) +
  geom_point() +
  geom_line() +
  facet_wrap("finclip_id")

## Caculate Pcrit!        
pcrit_md_lc <- md_lc %>%
  filter(expt_period == "pcrit") %>%
  unnest() %>%
  rename(MO2 = mo2_raw,
         DO = po2_torr) %>%
  group_by(date, finclip_id, expt_period) %>%
  nest() %>%
  left_join((mmr_smr_md_lc %>% dplyr::select(finclip_id,MMR,SMR_mlnd,
                                            SMR_q20)),
            by = "finclip_id") %>%
  mutate(data_o2crit = data %>% purrr::map(~ dplyr::select(.,MO2,DO)),
         pcrit_mlnd_object = purrr::map2(data_o2crit, SMR_mlnd, ~ calcO2crit(.x, .y, .y)),
         pcrit_mlnd = pcrit_mlnd_object %>% purrr::map_dbl("o2crit"),
         pcrit_q20_object = purrr::map2(data_o2crit, SMR_q20, ~ calcO2crit(.x, .y, .y)),
         pcrit_q20 = pcrit_q20_object %>% purrr::map_dbl("o2crit"))

pcrit_md_lc %>% dplyr::select(MMR,SMR_mlnd,SMR_q20,pcrit_mlnd,pcrit_q20)
pcrit_md_lc$pcrit_mlnd_object[[1]]
plotO2crit(pcrit_md_lc$pcrit_q20_object[[1]])
