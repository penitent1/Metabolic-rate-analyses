#install.packages("ggthemes", dependencies = TRUE)

library(dplyr)
library(tidyr)
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

probe13 <- smr.data[smr.data$probe == 'NFB0013', ] # Update as appropriate
probe13.enbi <- probe13[probe13$spps == 'enbi', ] # Update as appropriate
probe13.enbi.13aug <- probe13.enbi[probe13.enbi$date.day == '13' &
                                     probe13.enbi$date.month == 'aug' , ] # Update as appropriate
#probe14.olma.trial2.16c <- probe14.olma.trial2[probe14.olma.trial2$temp == '16', ] # Update as appropriate
head(probe13.enbi.13aug) # Update as appropriate
str(probe13.enbi.13aug) # Update as appropriate

### SMR estimate function

five.hr.plus.data <- probe13.enbi.13aug[probe13.enbi.13aug$time.hrs > 5, ]

smr <- calcSMR(five.hr.plus.data$mo2) # Update as appropriate
smr

smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best

  #############################

plot.smr <- ggplot(data=probe13.enbi.13aug, aes(x=probe13.enbi.13aug$time.hrs, 
                                                  y=probe13.enbi.13aug$mo2)) 
plot.smr + (geom_point(size = 1.5)) +
  geom_hline(yintercept = 2.07) +
  labs(x = "Time (hrs)",
       y = "MO2 (umol O2/g/hr)") +
  theme_base()

#plot.pcrit.bestsmr +
#  scale_colour_discrete(name = "Fish ID",
#                        breaks = c("clgl.big", "clgl.small", "olma.big"),
#                        labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1")) +
#  scale_shape_discrete(name = "Fish ID",
#                       breaks = c("clgl.big", "clgl.small", "olma.big"),
#                       labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1")) +
#  theme_base()

  ###############################################
  ###############################################
  ###
  ###   Estimate Pcrit using the SMR obtained above
  ###
  ###############################################
  ###############################################

pcrit.data <-read.csv(file.choose())
head(pcrit.data)

calcO2crit(pcrit.data, 2.07) # Enter value of SMR obtained above here, after "pcrit.data
#?calcO2crit

plotO2crit(calcO2crit(pcrit.data, 2.07))

### In torr
# (O2crit.%sat/100)*P.ATM.KPA*760*0.2095/101.325



  #######################
  ## Subset and create datasets for comparing 16 hr vs "16+ hr" SMR estimates
  #######################

smr.data.16hr <- smr.data[smr.data$time.hrs < 16, ]
head(smr.data.16hr)

# smr.data.16hr$time.hrs[smr.data.16hr$spps == "arha"] # Update as appropriate

### NFB00--: Specify probe and species for analysis

probe14 <- smr.data.16hr[smr.data.16hr$probe == 'NFB0014', ] # Update as appropriate
probe14.olma <- probe14[probe14$spps == 'olma', ] # Update as appropriate
probe14.olma.trial2 <- probe14.olma[probe14.olma$trial.no == '2', ] # Update as appropriate
head(probe14.olma.trial2) # Update as appropriate
str(probe14.olma.trial2) # Update as appropriate

### SMR estimate function

smr <- calcSMR(probe14.olma.trial2$mo2) # Update as appropriate
smr

smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best

### Estimate Pcrit using the SMR obtained above

pcrit.data <-read.csv(file.choose())
head(pcrit.data)

calcO2crit(pcrit.data, 2.21) # Enter value of SMR obtained above here, after "pcrit.data
#?calcO2crit

plotO2crit(calcO2crit(pcrit.data, 2.21))

### In torr
# (O2crit.%sat/100)*P.ATM.KPA*760*0.2095/101.325

###################################################
###################################################
###################################################
##
### Plotting data for visualization
##
###################################################
###################################################
###################################################

md <-read.csv(file.choose())
head(md)
str(md)

## Variable objects for use in plotting

#fish.id.mass <- pcrit.data$fish.id.mass
spps <- md$spps
#pcrit.type <- pcrit.data$p.crit.type
pcrit.r <- md$pcrit.r.avg
pcrit.regress <- md$pcrit.regress.avg
#smr.time <- pcrit.data$smr.time
#smr.best <- pcrit.data$best.smr
pcrit.r.sem <- md$pcrit.r.sem
pcrit.regress.sem <- md$pcrit.sem
temp <- md$temp

##################################################################################
### Plot data to visualize whether "pcrit type" affects pcrit measurements.
##################################################################################

  ### Plot: Pcrit via "best smr" fishMO2 o2crit estimate

plot.pcrit.r <- ggplot(data=md, aes(x=temp, 
                                   y=pcrit.r, 
                                   group=spps, 
                                   shape=spps, 
                                   colour=spps)) + 
  geom_errorbar(aes(ymin=pcrit.r-pcrit.r.sem, ymax=pcrit.r+pcrit.r.sem), width=0.25) +
  geom_line(size = 1.5, linetype = 2) + 
  geom_point(size = 3) +
  labs(title="Pcrit vs Temperature",
       x = "Temperature (degrees C)",
       y = "Pcrit (torr)")

plot.pcrit.r +
  scale_colour_discrete(name = "Species",
                        breaks = c("arfe", "arha", "arla", "clgl", "olma"),
                        labels = c("Padded", "Scalyhead", "Smoothhead", "Mosshead", "Tidepool")) +
  scale_shape_discrete(name = "Species",
                       breaks = c("arfe", "arha", "arla", "clgl", "olma"),
                       labels = c("Padded", "Scalyhead", "Smoothhead", "Mosshead", "Tidepool")) +
  theme_base()

  ### Plot: Pcrit via "REGRESS" estimate

plot.pcrit.regress <- ggplot(data=pcrit.data, aes(x=pcrit.type, 
                                                  y=pcrit.regress, 
                                                  group=fish.id.mass, 
                                                  shape=fish.id.mass, 
                                                  colour=fish.id.mass)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) +
  labs(title="Pcrit via the REGRESS program",
       x = "Pcrit type",
       y = "Pcrit (torr)")
plot.pcrit.regress +
  scale_colour_discrete(name = "Fish ID",
                        breaks = c("clgl.big", "clgl.small", "olma.big"),
                        labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("clgl.big", "clgl.small", "olma.big"),
                       labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1")) +
  theme_base()
