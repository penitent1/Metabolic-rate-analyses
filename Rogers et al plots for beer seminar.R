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

df <- read.csv(file.choose(), stringsAsFactors = FALSE,
               strip.white = TRUE, na.strings = c("NA","."))
df <- as_tibble(df)

### Trim dataset down to acute temp change, seawater fish pcrit and mo2 data

acute_sw <- df %>%
  dplyr::select(species,sw_or_fw,reference,acute_or_acclimated,acclimation_temp_c,
                test_temp_c,pcrit_mean,pcrit_sd,pcrit_units,pcrit_n,mo2_mean,
                mo2_sd,mo2_units,mo2_n,mass_mean_g,mass_sd_g,comments) %>%
  filter(sw_or_fw == "sw", acute_or_acclimated == "acute")

acute_sw_pcrit_mo2 <- acute_sw %>%
  dplyr::select(species,reference,acclimation_temp_c,test_temp_c,pcrit_mean,pcrit_sd,pcrit_units,pcrit_n,mo2_mean,
                mo2_sd,mo2_units,mo2_n) %>%
  filter(!is.na(pcrit_mean))

## Create a vector in dataframe for converted pcrit values (to torr)
acute_sw_pcrit_mo2$pcrit_mean_torr <- c(1:18)
acute_sw_pcrit_mo2$pcrit_sd_torr <- c(1:18)

## Note! mg O2/L conversion based on conversion from Loligo oxygen calculator
## For mgO2_per_L: Hilton et al, 15 or 25 C
for (i in 1:nrow(acute_sw_pcrit_mo2)) {
  if(acute_sw_pcrit_mo2$pcrit_units[i]=="mgO2_per_L") {
    if(acute_sw_pcrit_mo2$test_temp_c[i]==15){
    acute_sw_pcrit_mo2$pcrit_mean_torr[i] <- acute_sw_pcrit_mo2$pcrit_mean[i]*(156.46753/8.18666)
    acute_sw_pcrit_mo2$pcrit_sd_torr[i] <- acute_sw_pcrit_mo2$pcrit_sd[i]*(156.46753/8.18666)
    } else {
      acute_sw_pcrit_mo2$pcrit_mean_torr[i] <- acute_sw_pcrit_mo2$pcrit_mean[i]*(154.17018/6.74635)
      acute_sw_pcrit_mo2$pcrit_sd_torr[i] <- acute_sw_pcrit_mo2$pcrit_sd[i]*(154.17018/6.74635)
    }
  } else if(acute_sw_pcrit_mo2$pcrit_units[i]=="percent_sat"){
    acute_sw_pcrit_mo2$pcrit_mean_torr[i] <- (acute_sw_pcrit_mo2$pcrit_mean[i]/100)*101*760*0.2095/101.325 # 101 kPa is a typical air pressure at Lizard Island where Nilsson et al study occurred
    acute_sw_pcrit_mo2$pcrit_sd_torr[i] <- (acute_sw_pcrit_mo2$pcrit_sd[i]/100)*101*760*0.2095/101.325
    }
  else {
    acute_sw_pcrit_mo2$pcrit_mean_torr[i] <- acute_sw_pcrit_mo2$pcrit_mean[i]
    acute_sw_pcrit_mo2$pcrit_sd_torr[i] <- acute_sw_pcrit_mo2$pcrit_sd[i]
    }
}

## Create a vector in dataframe for converted MO2 values (to umol O2/g/hr)
acute_sw_pcrit_mo2$mo2_mean_umol <- c(1:18)
acute_sw_pcrit_mo2$mo2_sd_umol <- c(1:18)

## mgO2_kg_hr: O2 = 32 mg/mmol, (mg/kg/hr)*(mmol/32mg)*(1000umol/mmol)*(1kg/1000g)=(umol O2/g/hr)
for (i in 1:nrow(acute_sw_pcrit_mo2)) {
  if(acute_sw_pcrit_mo2$mo2_units[i]=="mgO2_g_hr") {
    acute_sw_pcrit_mo2$mo2_mean_umol[i] <- acute_sw_pcrit_mo2$mo2_mean[i]*1000/32
    acute_sw_pcrit_mo2$mo2_sd_umol[i] <- acute_sw_pcrit_mo2$mo2_sd[i]*1000/32
  } else {
    acute_sw_pcrit_mo2$mo2_mean_umol[i] <- acute_sw_pcrit_mo2$mo2_mean[i]/32
    acute_sw_pcrit_mo2$mo2_sd_umol[i] <- acute_sw_pcrit_mo2$mo2_sd[i]/32
  }
}

## Dumb way but can't figure out how else to get Hilton et al's 20C acclimation temp data
acute_sw_pcrit_mo2_no20 <- acute_sw_pcrit_mo2[acute_sw_pcrit_mo2$acclimation_temp_c != 20,]

acute_sw_pcrit_mo2_no20 <- acute_sw_pcrit_mo2_no20 %>%
  mutate(pcrit_sem_torr = pcrit_sd_torr/sqrt(pcrit_n),
         mo2_sem_umol = mo2_sd_umol/sqrt(mo2_n))

# Colored by species ID: Pcrit 
pcrit_vs_temp_sw_plot <- ggplot(acute_sw_pcrit_mo2_no20, aes(x=test_temp_c, y=pcrit_mean_torr, colour = species)) +
  geom_errorbar(mapping = aes(x=test_temp_c, ymin=(pcrit_mean_torr-pcrit_sem_torr), ymax=(pcrit_mean_torr+pcrit_sem_torr)), width=0.75, size=1.25) + #position=pd, 
  geom_line(size=1.25)  + #position=pd, 
  geom_point(size = 4) + #position=pd, 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression("P"["crit"]~("Torr"))) #x = expression(paste("Temperature (",degree~C,")"))
pcrit_vs_temp_sw_plot

# Colored by species ID: MO2min
mo2_vs_temp_sw_plot <- ggplot(acute_sw_pcrit_mo2_no20, aes(x=test_temp_c, y=mo2_mean_umol, colour = species)) +
  geom_errorbar(mapping = aes(x=test_temp_c, ymin=(mo2_mean_umol-mo2_sem_umol), ymax=(mo2_mean_umol+mo2_sem_umol)), width=0.75, size=1.25) + #position=pd, 
  geom_line(size=1.25)  + #position=pd, 
  geom_point(size = 4) + #position=pd, 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")))
mo2_vs_temp_sw_plot

# Colored by species ID: MO2min vs Pcrit
mo2_vs_pcrit_sw_plot <- ggplot(acute_sw_pcrit_mo2_no20, aes(x=mo2_mean_umol, y=pcrit_mean_torr, colour = species)) +
  geom_errorbar(mapping = aes(x=mo2_mean_umol, ymin=(pcrit_mean_torr-pcrit_sem_torr), ymax=(pcrit_mean_torr+pcrit_sem_torr)), width=0.75, size=1.25) + #position=pd,
  geom_errorbarh(mapping = aes(xmax = mo2_mean_umol+mo2_sem_umol, xmin = mo2_mean_umol-mo2_sem_umol, height = 2.5), size=1.5) + 
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
mo2_vs_pcrit_sw_plot








##################################################################################
##################################################################################



##      FRESHWATER AND SEAWATER SPECIES

##      ACUTE TEMP EFFECTS ONLY


## Pcrit, MO2 vs temp analysis in FW and SW spps combined
## Still acute temp effects  only




### Trim dataset down to acute temp change, pcrit and mo2 data

acute <- df %>%
  dplyr::select(species,sw_or_fw,reference,acute_or_acclimated,acclimation_temp_c,
                test_temp_c,pcrit_mean,pcrit_sd,pcrit_units,pcrit_n,mo2_mean,
                mo2_sd,mo2_units,mo2_n,mass_mean_g,mass_sd_g,comments) %>%
  filter(acute_or_acclimated == "acute")

acute_pcrit_mo2 <- acute %>%
  dplyr::select(species,sw_or_fw,reference,acclimation_temp_c,test_temp_c,pcrit_mean,pcrit_sd,pcrit_units,pcrit_n,mo2_mean,
                mo2_sd,mo2_units,mo2_n) %>%
  filter(!is.na(pcrit_mean))

## Create a vector in dataframe for converted pcrit values (to torr)
length(acute_pcrit_mo2$species)
acute_pcrit_mo2$pcrit_mean_torr <- c(1:27)
acute_pcrit_mo2$pcrit_sd_torr <- c(1:27)

## Note! mg O2/L conversion based on conversion from Loligo oxygen calculator
## For mgO2_per_L: Hilton et al, 15 or 25 C
for (i in 1:nrow(acute_pcrit_mo2)) {
  if(acute_pcrit_mo2$pcrit_units[i]=="mgO2_per_L") {
    if(acute_pcrit_mo2$test_temp_c[i]==15){
      acute_pcrit_mo2$pcrit_mean_torr[i] <- acute_pcrit_mo2$pcrit_mean[i]*(156.46753/8.18666)
      acute_pcrit_mo2$pcrit_sd_torr[i] <- acute_pcrit_mo2$pcrit_sd[i]*(156.46753/8.18666)
    } else {
      acute_pcrit_mo2$pcrit_mean_torr[i] <- acute_pcrit_mo2$pcrit_mean[i]*(154.17018/6.74635)
      acute_pcrit_mo2$pcrit_sd_torr[i] <- acute_pcrit_mo2$pcrit_sd[i]*(154.17018/6.74635)
    }
  } else if(acute_pcrit_mo2$pcrit_units[i]=="percent_sat"){
    acute_pcrit_mo2$pcrit_mean_torr[i] <- (acute_pcrit_mo2$pcrit_mean[i]/100)*101*760*0.2095/101.325 # 101 kPa is a typical air pressure at Lizard Island where Nilsson et al study occurred
    acute_pcrit_mo2$pcrit_sd_torr[i] <- (acute_pcrit_mo2$pcrit_sd[i]/100)*101*760*0.2095/101.325
  } else if(acute_pcrit_mo2$pcrit_units[i]=="kpa"){
    acute_pcrit_mo2$pcrit_mean_torr[i] <- acute_pcrit_mo2$pcrit_mean[i]*7.5 # 1 kPa = 7.5 torr
    acute_pcrit_mo2$pcrit_sd_torr[i] <- acute_pcrit_mo2$pcrit_sd[i]*7.5 # 1 kPa = 7.5 torr
  } else {
    acute_pcrit_mo2$pcrit_mean_torr[i] <- acute_pcrit_mo2$pcrit_mean[i]
    acute_pcrit_mo2$pcrit_sd_torr[i] <- acute_pcrit_mo2$pcrit_sd[i]
  }
}

## Create a vector in dataframe for converted MO2 values (to umol O2/g/hr)
acute_pcrit_mo2$mo2_mean_umol <- c(1:27)
acute_pcrit_mo2$mo2_sd_umol <- c(1:27)

## mgO2_kg_hr: O2 = 32 mg/mmol, (mg/kg/hr)*(mmol/32mg)*(1000umol/mmol)*(1kg/1000g)=(umol O2/g/hr)
for (i in 1:nrow(acute_pcrit_mo2)) {
  if(acute_pcrit_mo2$mo2_units[i]=="mgO2_g_hr") {
    acute_pcrit_mo2$mo2_mean_umol[i] <- acute_pcrit_mo2$mo2_mean[i]*1000/32
    acute_pcrit_mo2$mo2_sd_umol[i] <- acute_pcrit_mo2$mo2_sd[i]*1000/32
  } else {
    acute_pcrit_mo2$mo2_mean_umol[i] <- acute_pcrit_mo2$mo2_mean[i]/32
    acute_pcrit_mo2$mo2_sd_umol[i] <- acute_pcrit_mo2$mo2_sd[i]/32
  }
}

acute_pcrit_mo2 <- acute_pcrit_mo2 %>%
  mutate(pcrit_sem_torr = pcrit_sd_torr/sqrt(pcrit_n),
         mo2_sem_umol = mo2_sd_umol/sqrt(mo2_n))

# Colored by species ID: Pcrit 
pcrit_vs_temp_plot <- ggplot(acute_pcrit_mo2, aes(x=test_temp_c, y=pcrit_mean_torr, colour = species)) +
  geom_errorbar(mapping = aes(x=test_temp_c, ymin=(pcrit_mean_torr-pcrit_sem_torr), ymax=(pcrit_mean_torr+pcrit_sem_torr)), width=0.75, size=1.25) + #position=pd, 
  geom_line(size=1.25)  + #position=pd, 
  geom_point(size = 4) + #position=pd, 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression("P"["crit"]~("Torr"))) #x = expression(paste("Temperature (",degree~C,")"))
pcrit_vs_temp_plot

# Colored by species ID: MO2min
mo2_vs_temp_plot <- ggplot(acute_pcrit_mo2, aes(x=test_temp_c, y=mo2_mean_umol, colour = species)) +
  geom_errorbar(mapping = aes(x=test_temp_c, ymin=(mo2_mean_umol-mo2_sem_umol), ymax=(mo2_mean_umol+mo2_sem_umol)), width=0.75, size=1.25) + #position=pd, 
  geom_line(size=1.25)  + #position=pd, 
  geom_point(size = 4) + #position=pd, 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Temperature (",degree,"C)")),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")))
mo2_vs_temp_plot

# Colored by species ID: MO2min vs Pcrit
mo2_vs_pcrit_plot <- ggplot(acute_pcrit_mo2, aes(x=mo2_mean_umol, y=pcrit_mean_torr, colour = species)) +
  geom_errorbar(mapping = aes(x=mo2_mean_umol, ymin=(pcrit_mean_torr-pcrit_sem_torr), ymax=(pcrit_mean_torr+pcrit_sem_torr)), width=0.75, size=1.25) + #position=pd,
  geom_errorbarh(mapping = aes(xmax = mo2_mean_umol+mo2_sem_umol, xmin = mo2_mean_umol-mo2_sem_umol, height = 2.5), size=1.5) + 
  geom_line(size=1.25)  + #position=pd, 
  geom_point(size = 4) + #position=pd, 
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 28),
        axis.title.y = (element_text(size = 32, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 28),
        axis.title.x = element_text(size = 32, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")")),
       y = expression("P"["crit"]~("Torr")))
mo2_vs_pcrit_plot
