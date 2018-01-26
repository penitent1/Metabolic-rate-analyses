library(tidyverse)
library(ggthemes)
library(nlme)

yamanaka_mo2 <- read.csv("yamanaka_etal_2013_mo2_vs_temp.csv", stringsAsFactors = FALSE,
                     strip.white = TRUE, na.strings = c("NA","."))

yamanaka_mo2 <- as_tibble(yamanaka_mo2)

yamanaka_pcrit <- read.csv("yamanaka_etal_2013_pcrit_vs_temp.csv", stringsAsFactors = FALSE,
                         strip.white = TRUE, na.strings = c("NA","."))

yamanaka_pcrit <- as_tibble(yamanaka_pcrit)

mo2_mean <- yamanaka_mo2 %>%
  group_by(temp) %>%
  summarise(mass_avg = mean(body_weight_g),
            mass_sd = sd(body_weight_g),
            mo2_avg = mean(mo2_mgO2_perg_perh),
            mo2_sd = sd(mo2_mgO2_perg_perh))

pcrit_mean <- yamanaka_pcrit %>%
  group_by(temp) %>%
  summarise(mass_avg = mean(body_weight_g),
            mass_sd = sd(body_weight_g),
            pcrit_avg = mean(pcrit_mgO2_perL),
            pcrit_sd = sd(pcrit_mgO2_perL))

#### MO2 plots
mo2_vs_mass_plot <- ggplot(yamanaka_mo2, aes(x=body_weight_g, y=mo2_mgO2_perg_perh)) +
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
  labs(x = expression("Mass (g)"),
       y = expression(paste("Mo"[2]," (mg O"[2]," g"^-1," hr"^-1,")"))) #+
  #geom_abline(slope = 0.03985455, intercept = 29.61927908, 
  #            color="black", size=1) +
  #geom_abline(slope = 0.05706547, intercept = 28.81887854,
  #            color="blue", size=1) +
  #coord_cartesian(ylim = c(0,5))

mo2_vs_temp_plot <- ggplot(yamanaka_mo2, aes(x=temp, y=mo2_mgO2_perg_perh)) +
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
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression(paste("Mo"[2]," (mg O"[2]," g"^-1," hr"^-1,")"))) #+
#geom_abline(slope = 0.03985455, intercept = 29.61927908, 
#            color="black", size=1) +
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(0,5))

mo2mean_vs_temp_plot <- ggplot(mo2_mean, aes(x=temp, y=mo2_avg)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(mo2_avg-mo2_sd), ymax=(mo2_avg+mo2_sd)), width=0.75, size=1.25) +
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
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression(paste("Mo"[2]," (mg O"[2]," g"^-1," hr"^-1,")"))) #+
#geom_abline(slope = 0.03985455, intercept = 29.61927908, 
#            color="black", size=1) +
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(0,5))

#### Pcrit plots
pcrit_vs_mass_plot <- ggplot(yamanaka_pcrit, aes(x=body_weight_g, y=pcrit_mgO2_perL)) +
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
  labs(x = expression("Mass (g)"),
       y = expression(paste("P"["crit"]," (mg O"[2]," g"^-1," L"^-1,")"))) #+
#geom_abline(slope = 0.03985455, intercept = 29.61927908, 
#            color="black", size=1) +
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(0,5))

pcrit_vs_temp_plot <- ggplot(yamanaka_pcrit, aes(x=temp, y=pcrit_mgO2_perL)) +
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
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression(paste("P"["crit"]," (mg O"[2]," g"^-1," L"^-1,")"))) #+
#geom_abline(slope = 0.03985455, intercept = 29.61927908, 
#            color="black", size=1) +
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(0,5))

pcritmean_vs_temp_plot <- ggplot(pcrit_mean, aes(x=temp, y=pcrit_avg)) +
  geom_errorbar(mapping = aes(x=temp, ymin=(pcrit_avg-pcrit_sd), ymax=(pcrit_avg+pcrit_sd)), width=0.75, size=1.25) +
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
  labs(x = expression(paste("Temperature (",degree~C,")")),
       y = expression(paste("P"["crit"]," (mg O"[2]," g"^-1," L"^-1,")"))) #+
#geom_abline(slope = 0.03985455, intercept = 29.61927908, 
#            color="black", size=1) +
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
#coord_cartesian(ylim = c(0,5))
