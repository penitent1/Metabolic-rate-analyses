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

############################################
## Read sculpin pcrit and meta-data into R
############################################

md <- read.csv("SculpinPcritData_ComparativeAnalysisFormat_withPcritSlopes.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
md <- as_tibble(md)

md_slopesTorr <- md %>%
  mutate(slopes_torr = pcrit.slope.r.percentsat*100/160) ## 100% DO = 160 torr O2

## Pcrit and temp data, NA's filtered out, only "trial.no == 1" inlcuded
md_pcrit_smr <- md_slopesTorr %>%
  #select(species, spps, date, temp, fish.id, mass.g, smr.best, pcrit.r, trial.no) %>%
  filter(!is.na(smr.best)) %>%
  filter(!is.na(pcrit.r)) %>%
  filter(trial.no == 1)

## Averages for Pcrit and SMR per species per temperature, all 3 temps included
md_all_temps <- md_pcrit_smr %>%
  group_by(temp, species) %>%
  summarise(avg_pcrit = mean(pcrit.r), sd_pcrit = sd(pcrit.r), n_pcrit = length(pcrit.r),
            avg_smr = mean(smr.best), sd_smr = sd(smr.best), n_smr = length(smr.best),
            avg_pcrit_slope = mean(slopes_torr), sd_pcrit_slope = sd(slopes_torr), n_pcrit_slope = length(slopes_torr)) %>%
  group_by(temp, species) %>%
  mutate(sem_pcrit = sd_pcrit/n_pcrit,
         sem_smr = sd_smr/n_smr,
         sem_pcrit_slope = sd_pcrit_slope/n_pcrit_slope)

################
################


### Plots


###############

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

md_all_temps$spps_names <- species_names

# Colored by species ID: Pcrit vs Temp with legend
pcrit_vs_temp_plot <- ggplot(md_all_temps, aes(x=temp, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  geom_line(position=pd, size=1.25)  +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Species"))
pcrit_vs_temp_plot


######################################################

# Colored by species ID: ln(Pcrit) vs Temp with legend

kb <- 0.000086173324 ## Boltzmann constant in electron volts (eV)
kelvin <- 273.15 ## for converting Celsius to Kelvin

md_ln_pcrit <- md_all_temps %>%
  mutate(ln_pcrit = log(avg_pcrit),
         inv_temp = 1/((temp+273.15)*kb))

ln_pcrit_vs_temp_plot <- ggplot(md_ln_pcrit, aes(x=inv_temp, y=ln_pcrit, colour = spps_names)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  geom_line(position=pd, size=1.25)  +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Inverse temperature (1/T*kB, eV)"),
       y = expression("log P"["crit"]~("Torr")),
       colour = expression("Species"))
ln_pcrit_vs_temp_plot


# Colored by species ID: Pcrit vs MO2 at 12 C
md_12C <- md_all_temps %>%
  filter(temp == 12)

pcrit_vs_mo2_12c_plot <- ggplot(md_12C, aes(x=avg_smr, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=avg_smr, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.25) +
  geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 1.5), size=1.5) +
  geom_line(size=1.25)  +
  geom_point(size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Species"))
pcrit_vs_mo2_12c_plot


# Colored by species ID: Below-Pcrit slope from MO2 vs PO2 plot, vs temperature
pcrit_slope_vs_temp_plot <- ggplot(md_all_temps, aes(x=temp, y=avg_pcrit_slope, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit_slope-sem_pcrit_slope), ymax=(avg_pcrit_slope+sem_pcrit_slope)), position=pd, width=0.75, size=1.25) +
  geom_line(position=pd, size=1.25)  +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 20),
        axis.title.y = (element_text(size = 22, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 22, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression(paste("Oxyconforming slope"," (Mo"[2]," torr"^-1,")")),
       colour = expression("Species"))
pcrit_slope_vs_temp_plot



##############################

### Example Pcrit figure (MO2 vs PO2) from OLMA 12C trial on 14 Jun 2017 (Probe 14)


olma_pcrit <- read.csv("OLMA_12C_Pcrit_example_trace.csv", stringsAsFactors = FALSE,
                       strip.white = TRUE, na.strings = c("NA","."))
olma_pcrit <- as_tibble(olma_pcrit)
olma_pcrit_torr <- olma_pcrit[,c("po2","mo2_ms")]

olma_pcrit_plot <- ggplot(olma_pcrit_torr, aes(x=po2, y=mo2_ms)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd, size=1.25)  +
  geom_point(size = 3) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 26),
        axis.title.y = (element_text(size = 30, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 26),
        axis.title.x = element_text(size = 30, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Po"[2]~"(torr)"),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")"))) +
  #geom_abline(slope = 0.03985455, intercept = 29.61927908, 
  #            color="black", size=1) +
  #geom_abline(slope = 0.05706547, intercept = 28.81887854,
  #            color="blue", size=1) +
  coord_cartesian(ylim = c(0,2.5)) #+
  geom_segment(aes(x=x_min_1, y=y_min_1, xend=x_max_1, yend=y_max_1),
               color = "red", size=1) +
  geom_segment(aes(x=x_min_2, y=y_min_2, xend=x_max_2, yend=y_max_2),
               color = "red", size=1)

