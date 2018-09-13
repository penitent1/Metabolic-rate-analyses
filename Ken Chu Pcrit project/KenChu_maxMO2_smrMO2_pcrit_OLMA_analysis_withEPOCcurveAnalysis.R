library(tidyverse)
library(fishMO2)

## Import MO2 data: each file is one trial, one fish. This may need to change later.
md <- read_csv(file.choose()) %>%
  mutate(time_hrs = time_min/60,
         o2_kpa = o2_torr/0.133,
         mo2 = -1*umol_per_sec*3600,
         mo2_ms = mo2/4.5) %>% # 4.4 g = body mass
  group_by(experiment_period)

ggplot(md, aes(x = time_hrs, y = mo2_ms, colour = experiment_period)) +
  geom_point() +
  #geom_hline(aes(yintercept = 2.08)) +
  theme_classic()

md_smr <- md %>%
  filter(experiment_period == "mmr_smr",
         time_hrs > 10) %>%
  select(experiment_period, time_hrs, o2_torr, o2_kpa, mo2, mo2_ms)

ggplot(md_smr, aes(x = time_hrs, y = mo2_ms, colour = experiment_period)) +
  geom_point() +
  geom_hline(aes(yintercept = 2.08)) +
  theme_classic()

calcSMR(md_smr$mo2_ms)

ifelse(calcSMR(md_smr$mo2_ms)$CVmlnd > 5.4, calcSMR(md_smr$mo2_ms)$quant, calcSMR(md_smr$mo2_ms)$mlnd)

# Confidence interval multiplier for standard error
# Calculate t-statistic for confidence interval: 
# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
N_smr <- length(calcSMR(md_smr$mo2_ms)$cl)
conf.interval <- 0.95
ciMult <- qt(conf.interval/2 + .5, N_smr-1)
datac$ci <- datac$se * ciMult

ggplot(md_smr, aes(x = time_hrs, y = mo2_ms, colour = experiment_period)) +
  geom_point() +
  geom_hline(aes(yintercept=1.99), colour = "blue") +
  stat_smooth(method = "lm", colour = "red") +
  theme_classic()
  
md_pcrit <- md %>%
  filter(experiment_period == "pcrit") %>%
  select(experiment_period, time_hrs, o2_sat, mo2_ms) %>%
  rename(DO = o2_sat,
         MO2 = mo2_ms)

plotO2crit(calcO2crit(md_pcrit, SMR = 1.99, lowestMO2 = 1.99)) # 2.33 = Round 1 SMR

ggplot(md_pcrit, aes(x = DO, y = MO2)) +
  geom_point() +
  geom_vline(aes(xintercept = 13.4)) +
  theme_classic()

# Note that SMR was estimated using mo2 values after 10 hrs, not after 5 hrs
## Round 1: MMR=9.90, SMR=2.14, pcrit=35.2 torr (22% sat), AAS=7.76, FAS=4.63
## Round 2: MMR=9.12, SMR=2.22, pcrit=23.36 torr (14.6% sat), AS=6.95, FAS=4.20


#################################################################
###
### Attempt to calculate and plot EPOC 
### 
### for Round 1 of Tidepool sculpin MMR-SMR
###

md_epoc <- md %>%
  filter(experiment_period == "mmr_smr")

#write_csv(md_epoc[1:44,],
#          path = "C:/Users/derek/Documents/Metabolic-rate-analyses/Ken 4Jul2018 OLMA EPOC mo2 vs time hrs.csv")

epoc_nlm <- nls(mo2_ms ~ a*exp(k_1*time_hrs)+b*exp(k_2*time_hrs)+2.08,
                data = md_epoc[md_epoc$time_hrs<10,],
                start = list(a = 175,
                             k_1 = -3,
                             b = 20,
                             k_2 = -0.3))

epoc_func <- function(x) {7.528*exp(-2.074*x)+1.494*exp(-0.173*x)+2.08}

ggplot(md_epoc[md_epoc$time_hrs<10,], aes(x=time_hrs, y=mo2_ms)) +
  geom_point() +
  geom_hline(aes(yintercept=2.08), colour = "blue") +
  geom_hline(aes(yintercept=9.82), colour = "red") +
  theme_classic() +
  stat_function(fun = epoc_func)

## Old code for reference  
#  fn_olma <- function(x) {7.201119 * x ^ 0.7963} ## Using allometric equation

# Allometric scaling of MO2 vs body mass in Tidepool sculpins
#scaling_olma_mo2_plot <- ggplot(plot_check_olma_df,
#                                aes(x = ind_mass_g, y = ind_avg_umol_hr)) +
#  geom_point(size = 3) +
#  stat_function(fun = fn_olma, colour = "red", lwd = 1) +
#  stat_function(fun = fn_olma_mean, colour = "blue", lwd = 1) +
#  theme(panel.background = element_rect(fill = "white"),
#        axis.line = element_line(size = 1, colour = "black"),
#        panel.border = element_rect(linetype = "blank", fill = NA),
#        axis.text.y = element_text(size = 24),
#        axis.title.y = (element_text(size = 28, margin = margin(t = 0, r = 20, b = 0, l = 0))),
#        axis.text.x = element_text(size = 24),
#        axis.title.x = element_text(size = 28, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
#  labs(x = expression(paste("Body mass (g)")),
#       y = expression(paste(dot(M),"o"[2]," (",mu,"mol O"[2]," hr"^-1,")"))) +
#  annotate ("segment", x = 4.95, xend = 5.2, y = 40.176, yend = 40.176, size = 1.25, colour = "black", arrow = arrow ())

scaling_olma_mo2_plot

