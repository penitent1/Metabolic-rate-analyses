library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(mclust)
library(shape)
library(StreamMetabolism)
library(tidyverse)
library(fishMO2)

probeCorrData <- read_csv("Sculpins_MMR-SMR-Pcrit_Probe_correction_Data_2018.csv")
probeCorrData

probeCorr <- lm(probeCorrData$post_percent_sat[probeCorrData$date=="5-Mar-18"]~probeCorrData$pre_percent_sat[probeCorrData$date=="5-Mar-18"])
summary(probeCorr)
coef(probeCorr)
plot(probeCorrData$post_percent_sat~probeCorrData$pre_percent_sat)
abline(probeCorr)
  # post = 1.015*pre
  # pre = (post)/1.015

md <- read_csv("Sculpins_MMR-SMR-Pcrit_Data_2018.csv")

md_mo2 <- md %>%
  mutate(mo2 = umol_O2_per_sec*-1*3600,
         mo2_ms = umol_O2_per_sec*-1*3600/mass_g,
         timepoint_measurement_hrs = timepoint_measurement_minutes/60)

md_mo2_hypoxia <- md %>%
  filter(date == "5-Mar-18") %>%
  mutate(mo2 = umol_O2_per_sec*-1*3600,
         mo2_ms = umol_O2_per_sec*-1*3600/mass_g,
         timepoint_measurement_hrs = timepoint_measurement_minutes/60)

md_olma_mmr <- md %>%
  filter(measurement_period == "mo2_max") %>%
  mutate(mo2_max = umol_O2_per_sec*-1*3600,
         mo2_max_ms = umol_O2_per_sec*-1*3600/mass_g)
md_olma_mmr$mo2_max_ms

md_olma_smr <- md %>%
  filter(measurement_period == "mo2_routine") %>%
  mutate(mo2_routine = umol_O2_per_sec*-1*3600,
         mo2_routine_ms = umol_O2_per_sec*-1*3600/mass_g,
         timepoint_measurement_hrs = timepoint_measurement_minutes/60)

md_olma_smr_hypoxia <- md %>%
  filter(measurement_period == "mo2_routine") %>%
  filter(date == "5-Mar-18") %>%
  mutate(mo2_routine = umol_O2_per_sec*-1*3600,
         mo2_routine_ms = umol_O2_per_sec*-1*3600/mass_g,
         timepoint_measurement_hrs = timepoint_measurement_minutes/60)

md_olma_pcrit <- md %>%
  filter(measurement_period == "pcrit") %>%
  mutate(mo2_pcrit = umol_O2_per_sec*-1*3600,
         MO2 = umol_O2_per_sec*-1*3600/mass_g, ## Mass-specific MO2
         DO = (po2_torr*101.325*100)/(baro_pressure_kPa*0.2095*760),
         timepoint_measurement_hrs = timepoint_measurement_minutes/60)

md_olma_pcrit_hypoxia <- md %>%
  filter(date == "5-Mar-18") %>%
  filter(measurement_period == "pcrit") %>%
  mutate(mo2_pcrit = umol_O2_per_sec*-1*3600,
         MO2 = umol_O2_per_sec*-1*3600/mass_g, ## Mass-specific MO2
         DO = (po2_torr*101.325*100)/(baro_pressure_kPa*0.2095*760),
         timepoint_measurement_hrs = timepoint_measurement_minutes/60)

md_olma_pcrit_Rchabot <- md_olma_pcrit %>%
  dplyr::select(MO2, DO)

md_olma_pcrit_Rchabot_hypoxia <- md_olma_pcrit_hypoxia %>%
  dplyr::select(MO2, DO)

md_olma_pcrit_regress <- md_olma_pcrit %>%
  dplyr::select(po2_torr, MO2)

md_olma_pcrit_regress_hypoxia <- md_olma_pcrit %>%
  dplyr::select(po2_torr, MO2)

#write.table(mtcars, file = "mtcars.txt", sep = "\t",
#            row.names = TRUE, col.names = NA)

write.table(md_olma_pcrit_regress, file = "MMR-SMR-Pcrit_Pcrit_OLMA_12C_finclip-middle-new-wt_2Mar2018.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE)

write.table(md_olma_pcrit_regress_hypoxia, file = "MMR-SMR-Pcrit-PostHypoxia_Pcrit_OLMA_12C_finclip-middle-new-wt_5Mar2018.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE)

### SMR estimate function

five.hr.plus.data <- md_olma_smr[md_olma_smr$timepoint_measurement_hrs > 5, ]

ten.hr.plus.data_hypoxia <- md_olma_smr_hypoxia[md_olma_smr_hypoxia$timepoint_measurement_hrs > 10, ]

smr <- calcSMR(five.hr.plus.data$mo2_routine_ms) # Update as appropriate
smr

smr_hypoxia <- calcSMR(md_olma_smr_hypoxia$mo2_routine_ms)
smr_hypoxia

#smr <- calcSMR(probe10.enbi.19oct$mo2) # Use if I need ALL the data eg ENBI!
#smr

smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best

smr.check.best_hypoxia <- as.numeric(ifelse(smr_hypoxia$CVmlnd > 5.4, smr_hypoxia$quant[4], smr_hypoxia$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best_hypoxia

##############################################################

## Pcrit

pcrit <- calcO2crit(md_olma_pcrit_Rchabot, 1.76, lowestMO2 = 1.76) # Enter value of SMR obtained above here, after "pcrit.data"
pcrit
pcrit$o2crit ## Show only the calculated Pcrit

plotO2crit(calcO2crit(md_olma_pcrit_Rchabot, 1.76))#, lowestMO2 = 2.55))

### In torr
# (O2crit.%sat/100)*P.ATM.KPA*760*0.2095/101.325

## Pcrit in hypoxia

pcrit_hypoxia <- calcO2crit(md_olma_pcrit_Rchabot_hypoxia, 1.76, lowestMO2 = 1.76) # Enter value of SMR obtained above here, after "pcrit.data"
pcrit_hypoxia
pcrit_hypoxia$o2crit ## Show only the calculated Pcrit

plotO2crit(calcO2crit(md_olma_pcrit_Rchabot_hypoxia, 2.08, lowestMO2 = 2.08))#, lowestMO2 = 2.55))

### In torr
# (O2crit.%sat/100)*P.ATM.KPA*760*0.2095/101.325

##########################################################################
# Colored by measurement period: Mo2 vs time colored by measurement period
mo2_vs_time_plot <- ggplot(md_mo2, aes(x=timepoint_measurement_hrs, y=mo2_ms, colour = measurement_period)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1.76, size = 1.25) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Time (hrs)"),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Period"))
mo2_vs_time_plot

# In hypoxia: Colored by measurement period: Mo2 vs time colored by measurement period
mo2_vs_time_hypoxia_plot <- ggplot(md_mo2_hypoxia, aes(x=timepoint_measurement_hrs, y=mo2_ms, colour = measurement_period)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1.76, size = 1, colour = "red") +
  geom_hline(yintercept = 2.08, size = 1) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Time (hrs)"),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Period"))
mo2_vs_time_hypoxia_plot

# In hypoxia through Pcrit trial: Colored by measurement period: Mo2 vs time colored by measurement period
  ## Subset hypoxia data by time point of last Pcrit-period MO2
md_mo2_hypoxia_toPcrit <- md_mo2_hypoxia[md_mo2_hypoxia$timepoint_measurement_minutes < 194,]

mo2_vs_time_hypoxia_throughPcrit_plot <- ggplot(md_mo2_hypoxia_toPcrit, aes(x=timepoint_measurement_hrs, y=mo2_ms, colour = measurement_period)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1.76, size = 1, colour = "red") +
  geom_hline(yintercept = 2.08, size = 1) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Time (hrs)"),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Period")) +
  coord_cartesian(ylim = c(0,3))
mo2_vs_time_hypoxia_throughPcrit_plot

##########################################################################
pcrit_data_regress <- read.table("MMR-SMR-Pcrit_Pcrit_OLMA_12C_finclip-middle-new-wt_2Mar2018_withColNames.txt",
                                 header = TRUE)
pcrit_data_regress <- as_tibble(pcrit_data_regress)

for (i in 1:nrow(pcrit_data_regress)){
  if(pcrit_data_regress$mo2[i] >= 1.49){
    pcrit_data_regress$reg_or_conform[i] <- "regulating"
  } else{
    pcrit_data_regress$reg_or_conform[i] <- "conforming"
  }
}

## From REGRESS

#There are 4 points to the left, 25 points to the right.
#Regression equation for line 1 is: 
#  y = 3.36451804103635E-02 + 7.87511036792425E-02 x.
#Statistics for the first line:
#  Coefficient of determination: 0.98833456490998
#Slope: 7.87511036792425E-02
#Standard error of slope: 2.00509680132866E-02
#t-statistic: 3.92754622255936
#Degrees of freedom: 2

#Regression equation for line 2 is: 
#  y = 1.45409486798753 + 2.18525562152906E-03 x.
#Statistics for second line:
#  Coefficient of determination: 0.463965047557502
#Slope: 2.18525562152906E-03
#Standard error of slope: 8.56750740527285E-04
#t-statistic: 2.55063172771132
#Degrees of freedom: 23

# For below Pcrit line
x_min_1 <- 20.0717
x_max_1 <- 54.4
y_min_1 <- 0.56974 + 0.03401*x_min_1
y_max_1 <- 0.56974 + 0.03401*x_max_1

x_s <- c(20.0717,54.4)
y_s <- c(1.25242, 2.42)
x_y_dat <- data.frame(x_s, y_s)
lm_hehe_pcrit <- lm(y_s~x_s, data = x_y_dat)

# For above Pcrit line
x_min_2 <- 54.4
x_max_2 <- 149.5455
y_min_2 <- 2.42
y_max_2 <- 2.42

# Colored by pcrit relative level: above or below Pcrit
pcrit_plot <- ggplot(pcrit_data_regress, aes(x=po2, y=mo2, colour = reg_or_conform)) +
  geom_point(size = 3) +
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
  labs(x = expression("P"["crit"]~("Torr")),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Regulating or conforming")) #+
  #geom_segment(aes(x=x_min_1, y=y_min_1, xend=x_max_1, yend=y_max_1),
  #           color = "red", size=1) +
  #geom_segment(aes(x=x_min_2, y=y_min_2, xend=x_max_2, yend=y_max_2),
  #             color = "red", size=1)
pcrit_plot

# Colored by pcrit relative level: above or below Pcrit
pcrit_hypoxia_plot <- ggplot(pcrit_data_regress, aes(x=po2, y=mo2, colour = reg_or_conform)) +
  geom_point(size = 3) +
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
  labs(x = expression("P"["crit"]~("Torr")),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       colour = expression("Regulating or conforming")) #+
#geom_segment(aes(x=x_min_1, y=y_min_1, xend=x_max_1, yend=y_max_1),
#           color = "red", size=1) +
#geom_segment(aes(x=x_min_2, y=y_min_2, xend=x_max_2, yend=y_max_2),
#             color = "red", size=1)
pcrit_hypoxia_plot













##############################

# Extra functions

#md$time <- seq(0, (nrow(md)-1)*15, by=15)