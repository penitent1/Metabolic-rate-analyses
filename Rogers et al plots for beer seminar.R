library(ctv)
library(tidyverse)
library(ape)
library(rotl)
library(geiger)
library(ggthemes)
library(MCMCglmm)
library(kinship2)
library(phytools)
library(phangorn)
library(nlme)
library(visreg)
library(RColorBrewer)
library(caper)

### Import Rogers et al Pcrit data for fishes

df <- read.csv("Rogers et al 2016_Pcrit-Mo2-temp_database.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
df <- as_tibble(df)

### Trim dataset down to acute temp change, seawater fish pcrit and mo2 data

acute_sw <- df %>%
  dplyr::select(species,sw_or_fw,reference,acute_or_acclimated,acclimation_temp_c,
                test_temp_c,pcrit_mean,pcrit_sd,pcrit_units,pcrit_n,mo2_mean,
                mo2_sd,mo2_units,mo2_n,mass_mean_g,mass_sd_g,comments) %>%
  filter(sw_or_fw == "sw", acute_or_acclimated == "acute")

acute_sw_pcrit_mo2 <- acute_sw %>%
  dplyr::select(species,reference,test_temp_c,pcrit_mean,pcrit_sd,pcrit_units,pcrit_n,mo2_mean,
                mo2_sd,mo2_units,mo2_n) %>%
  filter(!is.na(pcrit_mean))

acute_sw_hilton <- acute_sw_pcrit_mo2 %>%
  filter(reference == "Hilton et al. (2008)") %>%
  mutate(pcrit_mean_torr = pcrit_mean*160/6.4,
         pcrit_sd_torr = pcrit_sd*160/6.4,
         mo2_mean_umol = mo2_mean*1000/32,
         mo2_sd_umol = mo2_sd*1000/32)

acute_sw_nilsson <- acute_sw_pcrit_mo2 %>%
  filter(reference == "Nilsson et al. (2010)") %>%
  mutate(pcrit_mean_torr = pcrit_mean*160/100,
         pcrit_sd_torr = pcrit_sd*160/100,
         mo2_mean_umol = mo2_mean/32,
         mo2_sd_umol = mo2_sd/32)

acute_sw_torr_umol <- bind_rows(acute_sw_hilton,acute_sw_nilsson)

# Colored by species ID: Pcrit 
pcrit_vs_temp_plot <- ggplot(acute_sw_torr_umol, aes(x=test_temp_c, y=pcrit_mean_torr, colour = species)) +
  geom_errorbar(mapping = aes(x=test_temp_c, ymin=(pcrit_mean_torr-pcrit_sd_torr), ymax=(pcrit_mean_torr+pcrit_sd_torr)), width=0.75, size=1.25) + #position=pd, 
  geom_line(size=1.25)  + #position=pd, 
  geom_point(size = 4) + #position=pd, 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression("P"["crit"]~("Torr"))) #x = expression(paste("Temperature (",degree~C,")"))

# Colored by species ID: MO2min
mo2_vs_temp_plot <- ggplot(acute_sw_torr_umol, aes(x=test_temp_c, y=mo2_mean_umol, colour = species)) +
  geom_errorbar(mapping = aes(x=test_temp_c, ymin=(mo2_mean_umol-mo2_sd_umol), ymax=(mo2_mean_umol+mo2_sd_umol)), width=0.75, size=1.25) + #position=pd, 
  geom_line(size=1.25)  + #position=pd, 
  geom_point(size = 4) + #position=pd, 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")))

# Colored by species ID: MO2min vs Pcrit
mo2_vs_temp_plot <- ggplot(acute_sw_torr_umol, aes(x=mo2_mean_umol, y=pcrit_mean_torr, colour = species)) +
  geom_errorbar(mapping = aes(x=mo2_mean_umol, ymin=(pcrit_mean_torr-pcrit_sd_torr), ymax=(pcrit_mean_torr+pcrit_sd_torr)), width=0.75, size=1.25) + #position=pd,
  geom_errorbarh(mapping = aes(xmax = mo2_mean_umol+mo2_sd_umol, xmin = mo2_mean_umol-mo2_sd_umol, height = 2.5), size=1.5) + 
  geom_line(size=1.25)  + #position=pd, 
  geom_point(size = 4) + #position=pd, 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       y = expression("P"["crit"]~("Torr")))
