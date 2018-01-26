hehe_raw <- read.csv("FOXY 16Aug2017 NFB0008 HEHE-II RespMGilbert1 20degC_Pcrit_MinMO2.csv", stringsAsFactors = FALSE,
                     strip.white = TRUE, na.strings = c("NA","."))
hehe_raw$time <- seq(0, (nrow(hehe_raw)-1)*15, by=15)
hehe_raw$time.hr <- hehe_raw$time/3600

hehe_mo2 <- read.csv("HEHE_20C_MO2_beauty_trace.csv", stringsAsFactors = FALSE,
                           strip.white = TRUE, na.strings = c("NA","."))
hehe_mo2 <- as_tibble(hehe_mo2)

hehe_pcrit <- read.csv("HEHE_20C_Pcrit_beauty_trace.csv", stringsAsFactors = FALSE,
                       strip.white = TRUE, na.strings = c("NA","."))
hehe_pcrit <- as_tibble(hehe_pcrit)
hehe_pcrit_torr <- hehe_pcrit[,c("po2","mo2_ms")]

hehe_raw_plot <- ggplot(hehe_raw, aes(x=time.hr, y=Oxygen)) +
  #geom_errorbar(mapping = aes(x=temp, ymin=(avg_pcrit-sem_pcrit), ymax=(avg_pcrit+sem_pcrit)), position=pd, width=0.75, size=1.25) +
  #geom_line(position=pd, size=1.25)  +
  geom_point(size = 1) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 0.5, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        axis.text.y = element_text(size = 26),
        axis.title.y = (element_text(size = 30, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 26),
        axis.title.x = element_text(size = 30, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression("Time (hrs)"),
       y = expression("Oxygen (% air saturation)")) #+
  #geom_abline(slope = 0.03985455, intercept = 29.61927908, 
  #            color="black", size=1) +
  #geom_abline(slope = 0.05706547, intercept = 28.81887854,
  #            color="blue", size=1) +
  #coord_cartesian(ylim = c(10,60))

hehe_mo2_plot <- ggplot(hehe_mo2, aes(x=time.hrs, y=mo2_ms)) +
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
  labs(x = expression("Time (hrs)"),
       y = expression(paste("Mo"[2]," (",mu,"mol O"[2]," g"^-1," hr"^-1,")"))) +
#geom_abline(slope = 0.03985455, intercept = 29.61927908, 
#            color="black", size=1) +
#geom_abline(slope = 0.05706547, intercept = 28.81887854,
#            color="blue", size=1) +
  coord_cartesian(ylim = c(0,5))

md[md$date=="16-Aug-17",c("species","date","pcrit.r","smr.best")]

# For below Pcrit line
x_min_1 <- 20.0717
x_max_1 <- 54.4
y_min_1 <- 0.56974 + 0.03401*x_min_1
y_max_1 <- 0.56974 + 0.03401*x_max_1
#y_min_1 <- 0.06747*x_min_1
#y_max_1 <- 0.06747*x_max_1
#y_min_1 <- 0.78214 + 0.02343*x_min_1
#y_max_1 <- 0.78214 + 0.02343*x_max_1 

x_s <- c(20.0717,54.4)
y_s <- c(1.25242, 2.42)
x_y_dat <- data.frame(x_s, y_s)
lm_hehe_pcrit <- lm(y_s~x_s, data = x_y_dat)

# For above Pcrit line
x_min_2 <- 54.4
x_max_2 <- 149.5455
y_min_2 <- 2.42
y_max_2 <- 2.42

hehe_pcrit_plot <- ggplot(hehe_pcrit_torr, aes(x=po2, y=mo2_ms)) +
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
  coord_cartesian(ylim = c(0,4)) +
  geom_segment(aes(x=x_min_1, y=y_min_1, xend=x_max_1, yend=y_max_1),
               color = "red", size=1) +
  geom_segment(aes(x=x_min_2, y=y_min_2, xend=x_max_2, yend=y_max_2),
               color = "red", size=1)
  