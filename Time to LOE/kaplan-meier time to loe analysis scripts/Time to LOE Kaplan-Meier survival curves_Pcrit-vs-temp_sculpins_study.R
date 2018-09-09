library(tidyverse)
library(ggfortify)
library(survival)

loe_data <- read_csv(file.choose())
loe_data

loe_data <- loe_data %>%
  filter(date != "28-Apr-18") %>%
  group_by(temp_c, species)
loe_data

loe_data_olma <- loe_data %>%
  filter(species == "Oligocottus_maculosus")

loe_data_clgl <- loe_data %>%
  filter(species == "Clinocottus_globiceps")

# Kaplan Meier Survival Curve: Mossheads
loe_km_clgl <- with(loe_data_clgl, Surv(time_loe_min, observed))
loe_km_clgl

loe_km_fit_clgl <- survfit(Surv(time_loe_min, observed) ~ 1, data = loe_data_clgl)
#loe_km_fit_clgl <- survfit(Surv(time_loe_min, observed) ~ temp_c, data=loe_data_clgl)
summary(loe_km_fit_clgl)
print(loe_km_fit_clgl) # This gives me the median! Yay!
autoplot(loe_km_fit_clgl)

## Tidepool sculpin data: Pcrit vs Time to LOE median by temperature

pcrit_loe_olma <- data.frame(temp_c = c(12,16,20),
                             pcrit_avg_torr = c(21.7,32.7,40.6),
                             loe_median_min = c(860,162,62))

## For tidepool sculpins: N = 6 per temperature
## Linear model of avg pcrit vs time to loe
lm_pcrit_loe <- lm(pcrit_loe_olma$pcrit_avg_torr~pcrit_loe_olma$loe_median_min)
summary(lm_pcrit_loe) # Adjusted R squared = 0.81

## Power model of avg pcrit vs time to loe
lm_pcrit_loe_power <- lm(log(pcrit_loe_olma$pcrit_avg_torr)~log(pcrit_loe_olma$loe_median_min))
summary(lm_pcrit_loe_power) # Adjusted R squared = 0.999
# y = a*x^b therefore ln(y) = ln(a) + b*ln(x)
# ln(a) = 4.695451 therefore a = exp(4.695451) = 109.4482

## For adding to plot
fn_linear_olma <- function(x) {39.17056 + -0.02077 * x}
fn_power_olma <- function(x) {109.4482 * x ^ -0.239042}

loe_raw_plots_temp <- ggplot(pcrit_loe_olma, aes(x = loe_median_min, y = pcrit_avg_torr, color = temp_c)) +
  geom_point(size=4) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(size = 1, colour = "black"),
        panel.border = element_rect(linetype = "blank", fill = NA),
        #legend.position = "none",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.y = (element_text(size = 20, margin = margin(t = 0, r = 20, b = 0, l = 0))),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 20, margin = margin(t = 20, r = 0, b = 0, l = 0))) +
  labs(x = expression(paste("Time to LOE"[50]," (min)")),
       y = expression("P"["crit"]~("Torr")),
       colour = expression(paste("Test temperature (",degree,C,")"))) +
  stat_function(fun = fn_power_olma, lwd = 0.5)
loe_raw_plots_temp

