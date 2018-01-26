library(tidyverse)
library(ggthemes)
library(phytools)
library(phangorn)
library(nlme)
library(visreg)
library(RColorBrewer)
library(caper)

############################################
## Read sculpin pcrit and meta-data into R
############################################

md <- read.csv("SculpinPcritData_ComparativeAnalysisFormat_withPcritSlopes.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
md <- as_tibble(md)

## Add a column with slopes transformed to a "per torr"
## instead of a "per percent saturation" denominator

md <- md %>%
  mutate(slope_torr = pcrit.slope.r.percentsat*100/160)

## Pcrit, temp, pcrit slope data, NA's filtered out, all trials inlcuded
md_all_trials <- md %>%
  filter(!is.na(smr.best)) %>%
  filter(!is.na(pcrit.r)) %>%
  filter(!is.na(pcrit.slope.r.percentsat))

## Pcrit, temp, pcrit slope data, only "trial.no == 1" inlcuded
md_trial_1 <- md_all_trials %>%
  filter(trial.no == 1)

## Averages for Pcrit, SMR, and pcrit slopes per species per temperature, all 3 temps included
md_temps_avg <- md_trial_1 %>%
  group_by(temp, species) %>%
  summarise(avg_pcrit = mean(pcrit.r), sd_pcrit = sd(pcrit.r), n_pcrit = length(pcrit.r),
            avg_smr = mean(smr.best), sd_smr = sd(smr.best), n_smr = length(smr.best),
            avg_slope = mean(pcrit.slope.r.percentsat), sd_slope = sd(pcrit.slope.r.percentsat), n_slope = length(pcrit.slope.r.percentsat),
            avg_slope_torr = mean(slope_torr), sd_slope_torr = sd(slope_torr), n_slope_torr = length(slope_torr)) %>%
  group_by(temp, species) %>%
  mutate(sem_pcrit = sd_pcrit/n_pcrit,
         sem_smr = sd_smr/n_smr,
         sem_slope = sd_slope/n_slope,
         sem_slope_torr = sd_slope_torr/n_slope_torr)

species_names <- c(
  "Artedius fenestralis",
  "Artedius harringtoni",
  "Artedius lateralis",
  "Blepsias cirrhosus",
  "Clinocottus globiceps",
  "Enophrys bison",
  "Hemilepidotus hemilepidotus",
  "Oligocottus maculosus",
  "Scorpaenichthys marmoratus"
)

md_temps_avg$spps_names <- species_names

#################################################################
#################################################################
#################################################################
#################################################################

## Plots

########
########
########
########

## Pcrit slopes vs temperature

plot_slopes_temp <- ggplot(data = md_temps_avg, aes(x=temp, y=avg_slope, colour=spps_names)) +
  geom_line(position = pd, size=1.25) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression("Slope"~"P"["crit"]~("Torr")))
plot_slopes_temp

## Pcrit slopes in torr not % saturation vs temperature

plot_slopes_torr_temp <- ggplot(data = md_temps_avg, aes(x=temp, y=avg_slope_torr, colour=spps_names)) +
  geom_line(position = pd, size=1.25) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression("Slope"~"P"["crit"]~("umol O2/g/hr/Torr")))
plot_slopes_torr_temp

## Pcrit slopes in torr not % saturation vs temperature with error bars

plot_slopes_torr_temp_ebars <- ggplot(data = md_temps_avg, aes(x=temp, y=avg_slope_torr, colour=spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_slope_torr-sem_slope_torr), ymax=(avg_slope_torr+sem_slope_torr)), width=0.75, size=1.25) +
  geom_line(size=1.25) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression("Slope"~"P"["crit"]~("umol O2/g/hr/Torr")))
plot_slopes_torr_temp_ebars

## Pcrit slopes vs Pcrit

plot_slopes_pcrit <- ggplot(data = md_temps_avg, aes(x=avg_pcrit, y=avg_slope, colour=spps_names)) +
  geom_errorbar(mapping = aes(x=avg_pcrit, ymin=(avg_slope-sem_slope), ymax=(avg_slope+sem_slope)), width=0.75, size=1.25) +
  geom_errorbarh(mapping = aes(xmax = avg_pcrit+sem_pcrit, xmin = avg_pcrit-sem_pcrit, height = 0.0025), size=1.5) +
  geom_line(size=1.25) +
  geom_point(size = 4) +
theme(panel.background = element_rect(fill = "white"),
      axis.line = element_line(size = 1, colour = "black"),
      panel.border = element_rect(linetype = "blank", fill = NA),
      axis.text.y = element_text(size = 28),
      axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
      axis.text.x = element_text(size = 28),
      axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("P"["crit"]~("Torr")),
       y = expression("Slope"~"P"["crit"]~("umol O2/g/hr/Torr")))
plot_slopes_pcrit 

## Pcrit slopes in torr not % sat vs Pcrit

plot_slopes_torr_pcrit <- ggplot(data = md_temps_avg, aes(x=avg_pcrit, y=avg_slope_torr, colour=spps_names)) +
  geom_errorbar(mapping = aes(x=avg_pcrit, ymin=(avg_slope_torr-sem_slope_torr), ymax=(avg_slope_torr+sem_slope_torr)), width=0.75, size=1.25) +
  geom_errorbarh(mapping = aes(xmax = avg_pcrit+sem_pcrit, xmin = avg_pcrit-sem_pcrit, height = 0.0025), size=1.5) +
  geom_line(size=1.25) +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("P"["crit"]~("Torr")),
       y = expression("Slope"~"P"["crit"]~("umol O2/g/hr/Torr")))
plot_slopes_torr_pcrit 
