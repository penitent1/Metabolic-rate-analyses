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
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses")
clgl_md <- read_csv("beer_seminar_31Oct2018_Pcrit_temp_mosshead_12to20_rmr_pcrit_data.csv")
clgl_md <- clgl_md %>%
  mutate(mo2_raw = d_umol_o2*3600*-1,
         MO2 = mo2_raw/mass_g) %>%
  rename(DO = do_percent) %>%
  group_by(date, temp, expt_period)

clgl_md_nested <- clgl_md %>% nest()

# Pcrit at 12
clgl_md %>% filter(expt_period == "pcrit", temp == 12) %>%
ggplot(aes(x = po2_torr*101.325/760, y = mo2_raw/mass_g, colour = temp)) +
  scale_colour_gradient(low = "blue", high = "red", breaks = c(12,20),
                        guide = FALSE) +
  geom_point(size = 6, colour = "blue") +
  labs(colour = expression(paste("Temperature (",degree,C,")"))) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,5),
                     breaks = seq(0,5,1)) +
  scale_x_continuous(name = expression(paste("Po"[2]," (kPa)")),
                     limits = c(0,22),
                     breaks = seq(0,21,3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        axis.text = element_text(size = rel(3), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"))
## Print 12 C figure: W = 1920, H = 1200

# Pcrit at 12 and 16, no guide
clgl_md %>% filter(expt_period == "pcrit", temp %in% c(12,16)) %>%
  ggplot(aes(x = po2_torr*101.325/760, y = mo2_raw/mass_g, colour = temp)) +
  geom_point(size = 6) +
  scale_colour_gradient(low = "blue", high = "maroon3", breaks = c(12,16),
                        guide = FALSE) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,5),
                     breaks = seq(0,5,1)) +
  scale_x_continuous(name = expression(paste("Po"[2]," (kPa)")),
                     limits = c(0,22),
                     breaks = seq(0,21,3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        axis.text = element_text(size = rel(3), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"))
#
#legend.title = element_text(size = rel(2.5)),
#legend.text = element_text(size = rel(2))

# Pcrit at 12, 16, and 20, no legend
clgl_md %>% filter(expt_period == "pcrit") %>%
  ggplot(aes(x = po2_torr*101.325/760, y = mo2_raw/mass_g, colour = temp)) +
  geom_point(size = 6) +
  scale_colour_gradient(low = "blue", high = "red", breaks = c(12,16,20),
                        guide = FALSE) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,5),
                     breaks = seq(0,5,1)) +
  scale_x_continuous(name = expression(paste("Po"[2]," (kPa)")),
                     limits = c(0,22),
                     breaks = seq(0,21,3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        axis.text = element_text(size = rel(3), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"))

# Pcrit at 12, 16, and 20, WITH LEGEND - Cut/paste legend for powerpoint
clgl_md %>% filter(expt_period == "pcrit") %>%
  ggplot(aes(x = po2_torr*101.325/760, y = mo2_raw/mass_g, colour = temp)) +
  geom_point(size = 6) +
  scale_colour_gradient(low = "blue", high = "red", breaks = c(12,16,20)) +
  labs(colour = expression(paste("Temperature (",degree,C,")"))) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol ",O[2]," g"^-1," hr"^-1,")")),
                     limits = c(0,5),
                     breaks = seq(0,5,1)) +
  scale_x_continuous(name = expression(paste("Po"[2]," (kPa)")),
                     limits = c(0,22),
                     breaks = seq(0,21,3)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = rel(5.5), margin = margin(t = 30)),
        axis.title.y = element_text(size = rel(5.5), margin = margin(r = 30)),
        axis.text = element_text(size = rel(3), colour = "black"),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"),
        legend.title = element_text(size = rel(2.5)),
        legend.text = element_text(size = rel(2)))


## Pcrit!
clgl_md_smr_pcrit <- clgl_md_nested %>%
  mutate(data_smr_pcrit = data %>% purrr::map(~ dplyr::select(., MO2, DO)))

smr_12_data <- clgl_md_smr_pcrit$data_smr_pcrit[[1]] %>% as.data.frame(.)
calcSMR(smr_12_data) # 1.43
smr_16_data <- clgl_md_smr_pcrit$data_smr_pcrit[[3]] %>% as.data.frame(.)
calcSMR(smr_16_data) # 2.79
smr_20_data <- clgl_md_smr_pcrit$data_smr_pcrit[[5]] %>% as.data.frame(.)
calcSMR(smr_20_data) # 4.10

pcrit_12_data <- clgl_md_smr_pcrit$data_smr_pcrit[[2]] %>% as.data.frame(.)
plotO2crit(calcO2crit(pcrit_12_data, SMR = 1.43)) # 1.43
smr_16_data <- clgl_md_smr_pcrit$data_smr_pcrit[[4]] %>% as.data.frame(.)
plotO2crit(calcO2crit()) # 2.79
smr_20_data <- clgl_md_smr_pcrit$data_smr_pcrit[[6]] %>% as.data.frame(.)
calcO2crit() # 4.10
