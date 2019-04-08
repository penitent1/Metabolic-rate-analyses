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
library(respR)
##################################
# Import data: POST Labchart
##################################
#setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit")
setwd("C:/Users/derek/OneDrive/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/raw foxy csv files/MMR-2daySMR-Pcrit_MMRatPcrit")
md_lc <- read_csv("Ken_Tidepool_CollatedMO2Data_MMR-SMR-Pcrit_PostLabChart.csv") # Time not importing properly from weird foxy format, use seq()
md_lc <- md_lc %>%
  mutate(mo2_raw = umol_o2_per_sec*3600*-1,
         time = as.double(time),
         mo2_ms = mo2_raw/mass_g) %>%
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

## Visualize ms PCrit data
md_lc %>% dplyr::select(-c(date)) %>% filter(expt_period == "pcrit") %>% unnest() %>%
  ggplot(aes(x = po2_torr, y = mo2_ms)) +
  geom_point() +
  geom_line() +
  facet_wrap("finclip_id")


## Caculate Pcrit! - Both Claireaux-Chabot and BSR methods        
pcrit_md_lc <- md_lc %>%
  filter(expt_period == "pcrit") %>%
  unnest() %>%
  rename(MO2 = mo2_ms,
         DO = po2_torr) %>%
  group_by(date, finclip_id, expt_period) %>%
  nest() %>%
  left_join((mmr_smr_md_lc %>% dplyr::select(finclip_id,MMR,SMR_mlnd,
                                            SMR_q20)),
            by = "finclip_id") %>%
  mutate(data_o2crit = data %>% purrr::map(~ dplyr::select(.,MO2,DO)),
         mass_g = data %>% purrr::map_dbl(~ mean(.$mass_g)),
         pcrit_mlnd_object = purrr::map2(data_o2crit, SMR_mlnd, ~ calcO2crit(.x, .y, .y)),
         pcrit_mlnd = pcrit_mlnd_object %>% purrr::map_dbl("o2crit"),
         pcrit_q20_object = purrr::map2(data_o2crit, SMR_q20, ~ calcO2crit(.x, .y, .y)),
         pcrit_q20 = pcrit_q20_object %>% purrr::map_dbl("o2crit"),
         pcrit_yu_object = data %>% purrr::map(~ dplyr::select(po2_torr, time)) %>%
           purrr::map(~ pcrit(., has.rate = TRUE, plot = FALSE)))

## MMR, SMR_mlnd, SMR_q20, Pcrit_mlnd, and Pcrit_q20 data summarized in a table:
pcrit_md_lc %>% dplyr::select(mass_g,MMR,SMR_mlnd,SMR_q20,pcrit_mlnd,pcrit_q20)

## Visualizing all the Pcrit plots, for both Pcrit calculations
plotO2crit(pcrit_md_lc$pcrit_mlnd_object[[1]])
plotO2crit(pcrit_md_lc$pcrit_q20_object[[1]])

plotO2crit(pcrit_md_lc$pcrit_mlnd_object[[2]])
plotO2crit(pcrit_md_lc$pcrit_q20_object[[2]])

plotO2crit(pcrit_md_lc$pcrit_mlnd_object[[3]])
plotO2crit(pcrit_md_lc$pcrit_q20_object[[3]])

plotO2crit(pcrit_md_lc$pcrit_mlnd_object[[4]])
plotO2crit(pcrit_md_lc$pcrit_q20_object[[4]])

plotO2crit(pcrit_md_lc$pcrit_mlnd_object[[5]])
plotO2crit(pcrit_md_lc$pcrit_q20_object[[5]])

plotO2crit(pcrit_md_lc$pcrit_mlnd_object[[6]])
plotO2crit(pcrit_md_lc$pcrit_q20_object[[6]])

## Broken stick regression - Pcrit estimates
pcrit_bsr_md_lc <- md_lc %>%
  filter(expt_period == "pcrit") %>%
  unnest() %>%
  group_by(date, finclip_id, expt_period)  %>%
  rename(MO2 = mo2_raw,
         DO = po2_torr) %>%
  nest() %>%
  mutate(respr_pcrit_ob = data %>% dplyr::select(MO2, DO)) #%>%
  #left_join((mmr_smr_md_lc %>% dplyr::select(finclip_id,MMR,SMR_mlnd,
  #                                           SMR_q20)),
  #          by = "finclip_id") %>%
  mutate(data_o2crit = data %>% purrr::map(~ dplyr::select(.,MO2,DO)),
         mass_g = data %>% purrr::map_dbl(~ mean(.$mass_g)),
         pcrit_mlnd_object = purrr::map2(data_o2crit, SMR_mlnd, ~ calcO2crit(.x, .y, .y)),
         pcrit_mlnd = pcrit_mlnd_object %>% purrr::map_dbl("o2crit"),
         pcrit_q20_object = purrr::map2(data_o2crit, SMR_q20, ~ calcO2crit(.x, .y, .y)),
         pcrit_q20 = pcrit_q20_object %>% purrr::map_dbl("o2crit"))

# Visualizing relationship: log(Pcrit_mlnd or Pcrit_q20) ~ log(body mass)
pcrit_md_lc %>% ggplot() +
  geom_point(aes(x = log(mass_g), y = log(pcrit_mlnd)), size = 4, colour = "red") +
  stat_smooth(aes(x = log(mass_g), y = log(pcrit_mlnd)), method = "lm", colour = "black") +
  geom_point(aes(x = log(mass_g), y = log(pcrit_q20)), size = 4, colour = "blue") +
  stat_smooth(aes(x = log(mass_g), y = log(pcrit_q20)), method = "lm", colour = "black") +
  scale_x_continuous(name = expression(paste("ln(mass) (ln(g))")),
                     limits = c(0,3),
                     breaks = seq(0,3,1)) +
  scale_y_continuous(name = expression(paste("ln(P"["crit"]," (ln(Torr))")),
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

# Visualizing Pcrit vs SMR, Pcrit vs MMR
pcrit_md_lc %>% ggplot() +
  geom_point(aes(x = MMR, y = pcrit_mlnd), size = 4, colour = "red") +
  geom_point(aes(x = SMR_mlnd, y = pcrit_mlnd), size = 4, colour = "blue") +
  scale_x_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol O"[2]," h"^-1,")")),
                     limits = c(0,50),
                     breaks = seq(0,50,10)) +
  scale_y_continuous(name = expression(paste("P"["crit"]," Torr)")),
                     limits = c(0,40),
                     breaks = seq(0,40,10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)))

# Visualizing MMR, (MO2 @ Pcrit), and MMR at Pcrit on a single plot with MO2 on y axis
mmr_at_pcrit_md_lc <- md_lc %>% filter(expt_period == "mmr_at_pcrit") %>% unnest() %>% group_by(finclip_id) %>%
  summarise(mmr_at_pcrit_mo2_mean = mean(mo2_raw),
            mmr_at_pcrit_mo2_sd = sd(mo2_raw))

pd <- position_dodge(0.1)
mmr_at_pcrit_plotting <- mmr_at_pcrit_md_lc %>% left_join(pcrit_md_lc, by = "finclip_id") %>% 
  gather(key = metabolic_state, value = mo2, MMR, SMR_mlnd, mmr_at_pcrit_mo2_mean) %>%
  mutate(met_state_plot_names = case_when(metabolic_state == "MMR" ~ expression(paste("Mo"[2]["MAX"])),
                                          metabolic_state == "mmr_at_pcrit_mo2_mean" ~ expression(paste("Mo"[2]["MAX"],"at P"["crit"])),
                                          metabolic_state == "SMR_mlnd" ~ expression(paste("Mo"[2]["STANDARD"]))))

mmr_at_pcrit_plotting %>%
  ggplot(aes(x = met_state_plot_names, y = mo2, group = finclip_id)) +
  geom_point(position = pd, size = 5, shape = 1, stroke = 2) +
  geom_line(position = pd, size = 1.5) +
  scale_x_discrete(name = expression(paste("Metabolic state"))) +
  scale_y_continuous(name = expression(paste(dot(M),"o"[2]," (",mu,"mol O"[2]," h"^-1,")")),
                     limits = c(0,50),
                     breaks = seq(0,50,10)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(size = rel(1.5), colour = "black"),
        axis.ticks = element_line(size = rel(5), colour = "black"),
        axis.ticks.length = unit(0.25, "cm"),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2), colour = "black"))

