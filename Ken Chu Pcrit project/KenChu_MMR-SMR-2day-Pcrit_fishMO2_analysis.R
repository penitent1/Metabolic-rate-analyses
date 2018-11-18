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

ggplot(md_lc$data[[19]], aes(x = time, y = mo2_raw/mass_g))+geom_point()+geom_line()
ggplot(md_lc$data[[20]], aes(x = do_percent_sat, y = mo2_raw/mass_g))+geom_point(size = 3)

smr_data <- md_lc$data[[19]] %>% dplyr::select(do_percent_sat,po2_torr,mo2_raw,mass_g) %>% 
  mutate(MO2 = mo2_raw/mass_g, DO = do_percent_sat) %>% as.data.frame(.)
head(smr_data)

max(smr_data$MO2)

smr_est <- calcSMR(smr_data$MO2) # Top-bottom SMR = 1.13

ifelse(smr_est$CVmlnd > 5.4, smr_est$quant, smr_est$mlnd)

ggplot(md_lc$data[[19]], aes(x = time, y = mo2_raw/mass_g))+
  geom_point()+
  geom_line()+
  geom_hline(yintercept = 1.38)

pcrit_data <- md_lc$data[[20]] %>% dplyr::select(do_percent_sat,po2_torr,mo2_raw,mass_g) %>% 
  mutate(MO2 = mo2_raw/mass_g, DO = do_percent_sat) %>% as.data.frame(.)
head(pcrit_data)

#write_excel_csv(tibble(do = pcrit_data$DO, mo2 = pcrit_data$MO2),"C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/olma_pcrit_intop-inbottom.csv")

plotO2crit(calcO2crit(pcrit_data, SMR = 1.38))

smr_pcrit_data <- bind_rows(smr_data, pcrit_data)
head(smr_pcrit_data)

plotO2crit(calcO2crit(smr_pcrit_data, SMR = 1.38))

## MMR at pcrit

mmr_at_pcrit <- md_lc$data[[18]] %>% 
  mutate(mo2_ms = mo2_raw/mass_g) %>% 
  dplyr::select(mass_g, time, mo2_raw, mo2_ms, do_percent_sat, po2_torr)

ggplot(mmr_at_pcrit, aes(x = time, y = mo2_ms)) +
  geom_point(size = 3)

mmr_at_pcrit

mean(mmr_at_pcrit$mo2_ms)

## MMR at Pcrit: Mo2 vs time for all trials

mmr_at_pcrit_vs_time <- md_lc %>%
  filter(expt_period == "mmr_at_pcrit")

# Look at dataset
mmr_at_pcrit_vs_time

# Unnest
mmr_at_pcrit_vs_time <- mmr_at_pcrit_vs_time %>% unnest() %>% 
  mutate(mo2_ms = mo2_raw/mass_g,
         smr_ms = case_when(finclip_id == "bottom" ~ 1.87,
                            finclip_id == "top" ~ 1.42,
                            finclip_id == "top_bottom" ~ 1.13,
                            finclip_id == "intop_inbottom" ~ 1.85,
                            finclip_id == "middle_stumpy" ~ 1.33,
                            finclip_id == "intop_bottom" ~ 1.38),
         mmr_ms = case_when(finclip_id == "bottom" ~ 8.82,
                            finclip_id == "top" ~ 7.75,
                            finclip_id == "top_bottom" ~ 6.83,
                            finclip_id == "intop_inbottom" ~ 8.71,
                            finclip_id == "middle_stumpy" ~ 9.40,
                            finclip_id == "intop_bottom" ~ 6.34),
         fish_name = case_when(finclip_id == "bottom" ~ "Fish 1",
                               finclip_id == "top" ~ "Fish 2",
                               finclip_id == "top_bottom" ~ "Fish 3",
                               finclip_id == "intop_inbottom" ~ "Fish 4",
                               finclip_id == "middle_stumpy" ~ "Fish 5",
                               finclip_id == "intop_bottom" ~ "Fish 6"))

# Plot: MO2 vs time, facet by individual
ggplot(mmr_at_pcrit_vs_time) +
  facet_wrap("fish_name") +
  geom_hline(aes(yintercept = mmr_ms), color = "red", size = 2.5) +
  geom_hline(aes(yintercept = smr_ms), color = "blue", size = 2.5) +
  geom_point(aes(x = time, y = mo2_ms), size = 7) +
  geom_line(aes(x = time, y = mo2_ms), size = 2.5) +
  scale_x_continuous(name = "Time post-chase (min)",
                     limits = c(0,35),
                     breaks = seq(0,30,10)) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol O"[2]," g"^-1," h"^-1,")")),
                     limits = c(0,10),
                     breaks = seq(0,10,2)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = rel(3), face = "bold"),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        axis.text = element_text(size = rel(3), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"))
## Print at W = 2000, H = 1239

## Raw Pcrit data, faceted by individual
raw_po2_vs_mo2 <- md_lc %>%
  filter(expt_period == "pcrit") %>% 
  unnest() %>%
  mutate(mo2_ms = mo2_raw/mass_g,
         smr_ms = case_when(finclip_id == "bottom" ~ 1.87,
                            finclip_id == "top" ~ 1.42,
                            finclip_id == "top_bottom" ~ 1.13,
                            finclip_id == "intop_inbottom" ~ 1.85,
                            finclip_id == "middle_stumpy" ~ 1.33,
                            finclip_id == "intop_bottom" ~ 1.38),
         mmr_ms = case_when(finclip_id == "bottom" ~ 8.82,
                            finclip_id == "top" ~ 7.75,
                            finclip_id == "top_bottom" ~ 6.83,
                            finclip_id == "intop_inbottom" ~ 8.71,
                            finclip_id == "middle_stumpy" ~ 9.40,
                            finclip_id == "intop_bottom" ~ 6.34),
         fish_name = case_when(finclip_id == "bottom" ~ "Fish 1",
                               finclip_id == "top" ~ "Fish 2",
                               finclip_id == "top_bottom" ~ "Fish 3",
                               finclip_id == "intop_inbottom" ~ "Fish 4",
                               finclip_id == "middle_stumpy" ~ "Fish 5",
                               finclip_id == "intop_bottom" ~ "Fish 6"))


ggplot(raw_po2_vs_mo2, aes(x = po2_torr, y = mo2_ms)) +
  facet_wrap("fish_name") +
  geom_point(size = 2) +
  scale_x_continuous(name = expression(paste("Po"[2], " (Torr)")),
                     limits = c(0,165),
                     breaks = seq(0,160,20)) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol O"[2]," g"^-1," h"^-1,")")),
                     limits = c(0,5),
                     breaks = seq(0,5,1)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        strip.text = element_text(size = rel(3), face = "bold"),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        axis.text = element_text(size = rel(3), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"))
## Print at W = 2000, H = 1239
