## Load appropriate packages
library(tidyverse)
library(fishMO2)

## Import data - someday I'll learn how to do this taking advantage of project format...
# Alter as appropriate
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/JEB_Pcrit_RI_analysis/Data")
md_o2 <- read_csv("MO2_PO2_data_from_temp_vs_pcrit_study.csv")

## Plot data: MO2 vs PO2 for SMR and Pcrit for each individual

# SMR visualization
md_o2 %>%
  mutate(mo2_raw_per_hr = umol_o2_per_sec*-1*3600,
         mo2_ms_per_hr = mo2_raw_per_hr/mass_g,
         time_hr = time_min/60) %>%
  select(species_plotting,finclip_id,time_hr,mo2_ms_per_hr, expt_period) %>%
  dplyr::filter(expt_period == "smr",
                mo2_ms_per_hr > 0) %>% # Couple of negative values from a mosshead pcrit trial
  ggplot(aes(x = time_hr, y = mo2_ms_per_hr)) +
  facet_wrap(species_plotting ~ finclip_id) +
  geom_point(size = 3) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,11),
                     breaks = seq(0,10,1)) +
  scale_x_continuous(name = expression(paste("Time (hr)")),
                     limits = c(0,26),
                     breaks = seq(0,24,2)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        strip.text = element_text(face = "bold", size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.text.x = element_text(size = rel(1.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"))

# Pcrit visualization by individual fish
md_o2 %>%
  mutate(mo2_raw_per_hr = umol_o2_per_sec*-1*3600,
         mo2_ms_per_hr = mo2_raw_per_hr/mass_g) %>%
  select(species_plotting,finclip_id,po2_torr,mo2_ms_per_hr, expt_period) %>%
  dplyr::filter(expt_period == "pcrit",
                mo2_ms_per_hr > 0) %>% # Couple of negative values from a mosshead pcrit trial
  ggplot(aes(x = po2_torr, y = mo2_ms_per_hr)) +
  facet_wrap(species_plotting ~ finclip_id) +
  geom_point(size = 3) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,8),
                     breaks = seq(0,6,1)) +
  scale_x_continuous(name = expression(paste("Po"[2]," (torr)")),
                     limits = c(0,170),
                     breaks = seq(0,160,10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        strip.text = element_text(face = "bold", size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.text.x = element_text(size = rel(1.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"))

# Pcrit visualization by SPECIES: All individuals on a single panel per species
md_o2 %>%
  mutate(mo2_raw_per_hr = umol_o2_per_sec*-1*3600,
         mo2_ms_per_hr = mo2_raw_per_hr/mass_g) %>%
  select(species_plotting,finclip_id,po2_torr,mo2_ms_per_hr, expt_period) %>%
  dplyr::filter(expt_period == "pcrit",
                mo2_ms_per_hr > 0) %>% # Couple of negative values from a mosshead pcrit trial
  ggplot(aes(x = po2_torr, y = mo2_ms_per_hr)) +
  facet_wrap("species_plotting") +
  geom_point(size = 3) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,7),
                     breaks = seq(0,7,1)) +
  scale_x_continuous(name = expression(paste("Po"[2]," (torr)")),
                     limits = c(0,170),
                     breaks = seq(0,160,10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        strip.text = element_text(face = "bold", size = rel(1.5)),
        axis.text.y = element_text(size = rel(1.5), colour = "black"),
        axis.text.x = element_text(size = rel(1.25), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"))
