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
library(Hmisc)

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

## CTmax data

ctmax <- read.csv("ct_max_maxdepth_data.csv", stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
ctmax <- as_tibble(ctmax)


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

# Colored by species ID: MO2min vs Temp with legend
smr_vs_temp_plot <- ggplot(md_all_temps, aes(x=temp, y=avg_smr, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_smr-sem_smr), ymax=(avg_smr+sem_smr)), position=pd, width=1, size=1.5) +
  geom_line(position=pd, size=1.25) +
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
       y = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Species"))
smr_vs_temp_plot

# Colored by species ID: MO2min vs Temp without legend
smr_vs_temp_plot_noLeg <- ggplot(md_all_temps, aes(x=temp, y=avg_smr, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_smr-sem_smr), ymax=(avg_smr+sem_smr)), position=pd, width=1, size=1.5) +
  geom_line(position=pd, size=1.25) +
  geom_point(position=pd, size = 4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression(paste("Mo"[2~"min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")))
smr_vs_temp_plot_noLeg

######################################################

# Colored by species ID: ln(Pcrit) vs Temp with legend

kb <- 0.000086173324 ## Boltzmann constant in electron volts (eV)
kelvin <- 273.15 ## for converting Celsius to Kelvin

md_ln_pcrit <- md_all_temps %>%
  mutate(ln_pcrit = log(avg_pcrit),
         inv_temp = 1/((temp+273.15)*kb))

ln_pcrit_vs_temp_plot_withPD <- ggplot(md_ln_pcrit, aes(x=inv_temp, y=ln_pcrit, colour = spps_names)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd, size=1.25)  +
  geom_smooth(method = lm, se = FALSE, position = pd, size=1.25) +
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
       y = expression("ln"("P"["crit"]~("Torr"))),
       colour = expression("Species"))
ln_pcrit_vs_temp_plot_withPD

ln_pcrit_vs_temp_plot_withoutPD <- ggplot(md_ln_pcrit, aes(x=inv_temp, y=ln_pcrit, colour = spps_names)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd, size=1.25)  +
  geom_smooth(method = lm, se = FALSE, size=1.25) +
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
  labs(x = expression("Inverse temperature (1/T*kB, eV)"),
       y = expression("ln"("P"["crit"]~("Torr"))),
       colour = expression("Species"))
ln_pcrit_vs_temp_plot_withoutPD


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
       colour = expression("Species")) +
  coord_cartesian(ylim = c(20,80),
                  xlim = c(0,5))
pcrit_vs_mo2_12c_plot

# Colored by species ID: Pcrit vs MO2 at 16 C
md_16C <- md_all_temps %>%
  filter(temp == 16)

pcrit_vs_mo2_16c_plot <- ggplot(md_16C, aes(x=avg_smr, y=avg_pcrit, colour = spps_names)) +
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
       colour = expression("Species")) +
  coord_cartesian(ylim = c(20,80),
                  xlim = c(0,5))
pcrit_vs_mo2_16c_plot

# Colored by species ID: Pcrit vs MO2 at 20 C
md_20C <- md_all_temps %>%
  filter(temp == 20)

pcrit_vs_mo2_20c_plot <- ggplot(md_20C, aes(x=avg_smr, y=avg_pcrit, colour = spps_names)) +
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
       colour = expression("Species")) +
  coord_cartesian(ylim = c(20,80),
                  xlim = c(0,5))
pcrit_vs_mo2_20c_plot

####################################################

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

# Colored by Pcrit: Below-Pcrit slope from MO2 vs PO2 plot, vs temperature
pd1 <- position_dodge(0.4)

pcrit_slope_vs_temp_pcritColored_plot <- ggplot(md_all_temps, aes(x=temp, y=avg_pcrit_slope, colour = avg_pcrit, group = spps_names)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit_slope-sem_pcrit_slope), ymax=(avg_pcrit_slope+sem_pcrit_slope)), position=pd1, width=0.75, size=1.25) +
  geom_line(position=pd1, size=1.25)  +
  geom_point(position=pd1, size = 6) +
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
       colour = expression(paste("P"["crit"]~"(torr)"))) +
  scale_colour_gradient(low = "#056BF1", high = "#F70828")
pcrit_slope_vs_temp_pcritColored_plot



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

##############################################################
##############################################################
  
  ### Fitting an OLS and mixed effects model to test association of MO2 and Pcrit
  
model1 <- lm(avg_pcrit~avg_smr, data = md_all_temps)

res1 <- residuals(model1)
plot(model1)
hist(res1)
qqnorm(model1)
qqline(model1)
  
summary(model1)  
anova(model1)
coef(model1)

plot(md_all_temps$avg_pcrit~md_all_temps$avg_smr)
abline(model1)

coef(model1)
x_min <- 0.956 #min(md_all_temps$avg_smr)
x_max <- 4.232 #max(md_all_temps$avg_smr)
y_min <- 24.677122 + 6.955594*x_min #using model1 coefficients
y_max <- 24.677122 + 6.955594*x_max #using model1 coefficients

x_min_lmm <- 0.956
x_max_lmm <- 4.232
y_min_lmm <- 11.28970 + 12.66513*x_min_lmm
y_max_lmm <- 11.28970 + 12.66513*x_max_lmm

### Mixed model with spps as random effect
  ## The idea here is that, while I don't 
model1_lmm <- lme(avg_pcrit ~ avg_smr, random = ~ avg_smr|spps_names, data = md_all_temps)
plot(model1_lmm)
summary(model1_lmm)
intervals(model1_lmm)
fitted(model1_lmm)
VarCorr(model1_lmm)
coef(model1_lmm)

anova.lme(model1_lmm)

# Colored by species ID: MO2min vs Pcrit all temperatures with legend
mo2_vs_pcrit_alltemps_plot <- ggplot(md_all_temps, aes(x=avg_smr, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=avg_smr, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 1.5), size=1.5) +
  geom_line(size=1.25) +
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
mo2_vs_pcrit_alltemps_plot

# Colored by species ID: MO2min vs Temp with legend
mo2_vs_pcrit_alltemps_globalLine_plot <- ggplot(md_all_temps, aes(x=avg_smr, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=avg_smr, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 1.5), size=1.5) +
  geom_line(size=1.25) +
  #geom_abline(intercept = 24.677122, slope = 6.955594, size = 1.25) +
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
  labs(x = expression(paste("Mo"["2min"]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Species")) +
  #coord_cartesian(ylim = c(20,80),
  #                xlim = c(0,5)) +
  geom_segment(aes(x=x_min_lmm, y=y_min_lmm, xend=x_max_lmm, yend=y_max_lmm),
               color = "black", size=1.5)
mo2_vs_pcrit_alltemps_globalLine_plot

# Colored by species ID: MO2min vs Pcrit with Global regression line and geom_smooth for species
mo2_vs_pcrit_alltemps_globalLine_sppsLine_plot <- ggplot(md_all_temps, aes(x=avg_smr, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=avg_smr, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.05, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = avg_smr+sem_smr, xmin = avg_smr-sem_smr, height = 1.5), size=1.5) +
  geom_smooth(method = "lm", se = FALSE, size=1.25) +
  #geom_line(size=1.25) +
  #geom_abline(intercept = 24.677122, slope = 6.955594, size = 1.25) +
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
       colour = expression("Species")) +
  coord_cartesian(ylim = c(20,80),
                  xlim = c(0,5)) +
  geom_segment(aes(x=x_min, y=y_min, xend=x_max, yend=y_max),
               color = "black", size=1.5)
mo2_vs_pcrit_alltemps_globalLine_sppsLine_plot

######################################################################
######################################################################

### CTmax versus Pcrit at 12, 16, and 20 C

md_all_temps
ctmax
md_ctmax <- full_join(md_all_temps, ctmax, by = "species")
md_ctmax <- md_ctmax %>%
  filter(!is.na(ct_max_avg)) %>%
  filter(!is.na(avg_pcrit))

## CTmax vs Pcrit at 12C
md_ctmax_12 <- md_ctmax %>%
  filter(temp == 12)

pcrit_ctmax_12_df <- data.frame(md_ctmax_12$ct_max_avg,md_ctmax_12$avg_pcrit)
pcrit_ctmax_12_matrix <- as.matrix(pcrit_ctmax_12_df)
rcorr(pcrit_ctmax_12_matrix)
#cor(x = md_ctmax_12$ct_max_avg, y = md_ctmax_12$avg_pcrit)

ctmax_vs_pcrit_12C_plot <- ggplot(md_ctmax_12, aes(x=ct_max_avg, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=ct_max_avg, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.1, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = ct_max_avg+ct_max_sem, xmin = ct_max_avg-ct_max_sem, height = 1), size=1.5) +
  #geom_smooth(method = "lm", se = FALSE, size=1.25) +
  #geom_line(size=1.25) +
  #geom_abline(intercept = 24.677122, slope = 6.955594, size = 1.25) +
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
  labs(x = expression(paste("CT"["max"]," (",degree,"C)")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Species")) +
  coord_cartesian(ylim = c(20,80))
ctmax_vs_pcrit_12C_plot

## CTmax vs Pcrit at 16C
md_ctmax_16 <- md_ctmax %>%
  filter(temp == 16)

pcrit_ctmax_16_df <- data.frame(md_ctmax_16$ct_max_avg,md_ctmax_16$avg_pcrit)
pcrit_ctmax_16_matrix <- as.matrix(pcrit_ctmax_16_df)
rcorr(pcrit_ctmax_16_matrix)

ctmax_vs_pcrit_16C_plot <- ggplot(md_ctmax_16, aes(x=ct_max_avg, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=ct_max_avg, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.1, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = ct_max_avg+ct_max_sem, xmin = ct_max_avg-ct_max_sem, height = 1), size=1.5) +
  #geom_smooth(method = "lm", se = FALSE, size=1.25) +
  #geom_line(size=1.25) +
  #geom_abline(intercept = 24.677122, slope = 6.955594, size = 1.25) +
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
  labs(x = expression(paste("CT"["max"]," (",degree,"C)")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Species")) +
  coord_cartesian(ylim = c(20,80))
ctmax_vs_pcrit_16C_plot

## CTmax vs Pcrit at 20C
md_ctmax_20 <- md_ctmax %>%
  filter(temp == 20)

pcrit_ctmax_20_df <- data.frame(md_ctmax_20$ct_max_avg,md_ctmax_20$avg_pcrit)
pcrit_ctmax_20_matrix <- as.matrix(pcrit_ctmax_20_df)
rcorr(pcrit_ctmax_20_matrix)

ctmax_vs_pcrit_20C_plot <- ggplot(md_ctmax_20, aes(x=ct_max_avg, y=avg_pcrit, colour = spps_names)) +
  geom_errorbar(mapping = aes(x=ct_max_avg, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), width=0.1, size=1.5) +
  geom_errorbarh(mapping = aes(xmax = ct_max_avg+ct_max_sem, xmin = ct_max_avg-ct_max_sem, height = 1), size=1.5) +
  #geom_smooth(method = "lm", se = FALSE, size=1.25) +
  #geom_line(size=1.25) +
  #geom_abline(intercept = 24.677122, slope = 6.955594, size = 1.25) +
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
  labs(x = expression(paste("CT"["max"]," (",degree,"C)")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression("Species")) +
  coord_cartesian(ylim = c(20,80))
ctmax_vs_pcrit_20C_plot