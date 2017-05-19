#install.packages("ggthemes", dependencies = TRUE)

library(ggplot2)
library(ggthemes)
library(mclust)
library(shape)
library(StreamMetabolism)

library(fishMO2)
#?fishMO2

### Use the following script, changing names and specs of objects,
  # to get SMR and "O2crit" estimates

smr.data <-read.csv(file.choose())
head(smr.data)

### NFB00--: Specify probe and species for analysis

probe12 <- smr.data[smr.data$probe == 'NFB0012', ] # Update as appropriate
probe12.arha <- probe12[probe12$spps == 'arha', ] # Update as appropriate
probe12.arha.trial1 <- probe12.arha[probe12.arha$trial.no == '1', ] # Update as appropriate
head(probe12.arha.trial1) # Update as appropriate
str(probe12.arha.trial1) # Update as appropriate

### SMR estimate function

smr <- calcSMR(probe12.arha.trial1$mo2) # Update as appropriate
smr
#> smr
#$mlnd
#1 
#1.330412

## $CVmlnd
## 1 
## 7.251364
smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best
## > smr.check.best
## [1] 1.315862

### Estimate Pcrit using the SMR obtained above

pcrit.data <-read.csv(file.choose())
head(pcrit.data)

calcO2crit(pcrit.data, 1.18) # Enter value of SMR obtained above here, after "pcrit.data
#?calcO2crit

plotO2crit(calcO2crit(pcrit.data, 1.18))

### In torr
#> (O2crit.%sat/100)*P.ATM.KPA*760*0.2095/101.325

###################################################
###################################################
###################################################
##
### Plotting data for visualization
##
###################################################
###################################################
###################################################

pcrit.data <-read.csv(file.choose())
head(pcrit.data)
str(pcrit.data)
#pcrit.data$pcrit.bestsmr <- as.numeric(pcrit.data$pcrit.bestsmr)

fish.id.mass <- pcrit.data$fish.id.mass
spps <- pcrit.data$spps
pcrit.type <- pcrit.data$p.crit.type
pcrit.bestsmr <- pcrit.data$pcrit.bestsmr
pcrit.regress <- pcrit.data$p.crit.regress

### Plot data to visualize whether "pcrit type" affects pcrit measurements.

### Stripchart of the Pcrit versus type of Pcrit - thanks Melissa!

pcrit.bestsmr_plot <- ggplot(pcrit.data, 
                                    aes(x=pcrit.type, 
                                        y=pcrit.bestsmr, 
                                        color=fish.id.mass
                                    )) + 
  geom_jitter(position=position_jitter(0.2))+
  labs(title="Pcrit type for each fish tested - olma and clgl",
       x = "Pcrit type",
       y = "Pcrit (torr)")
pcrit.bestsmr_plot
#pcrit.bestsmr_plot + stat_summary(fun.data = mean_sdl, 
#                                         fun.args = list(mult = 1),
#                                         geom = "pointrange",
#                                         color = "black")

### Stripchart of the Ventilation amplitude vs salinity with lines connecting fish

pcrit.bestsmr.lines_plot <- ggplot(pcrit.data, 
                            aes(x=pcrit.type, 
                                y=pcrit.bestsmr, 
                            )) + 
  geom_line(aes(group = fish.id.mass, 
                color = factor(fish.id.mass)),
            size = 2) +   ### Specify SIZE (ie thickness) of lines OUTSIDE of "aes()" function
  labs(title="Pcrit type for each fish tested - Tidepool and Mosshead sculpins",
       x = "Pcrit type",
       y = "Pcrit (torr)")

## Print plot
pcrit.bestsmr.lines_plot +
  scale_shape_discrete(name="Fish\nIdentity",
                          breaks=c("clgl.big", "clgl.small", "olma.big"),
                          labels=c("Mosshead 1", "Mosshead 2", "Tidepool 1")) + 
  theme_base()


##################################################################################

# Specify colour and shape
plot.pcrit.bestsmr <- ggplot(data=pcrit.data, aes(x=pcrit.type, 
                                   y=pcrit.bestsmr, 
                                   group=fish.id.mass, 
                                   shape=fish.id.mass, 
                                   colour=fish.id.mass)) + 
  geom_line() + 
  geom_point() +
  labs(title="Pcrit type for each fish tested - Tidepool and Mosshead sculpins",
       x = "Pcrit type",
       y = "Pcrit (torr)")
plot.pcrit.bestsmr +
  scale_colour_discrete(name = "Fish ID",
                        breaks = c("clgl.big", "clgl.small", "olma.big"),
                        labels = c("Mosshead 1", " Mosshead 2", "Tidepool 1")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("clgl.big", "clgl.small", "olma.big"),
                       labels = c("Mosshead 1", " Mosshead 2", "Tidepool 1")) +
  theme_base()
