#install.packages("ggthemes", dependencies = TRUE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(mclust)
library(shape)
library(StreamMetabolism)

library(fishMO2)
?fishMO2

### Use the following script, changing names and specs of objects,
  # to get SMR and "O2crit" estimates

smr.data <-read.csv(file.choose())
head(smr.data)

### NFB00--: Specify probe and species for analysis

probe14 <- smr.data[smr.data$probe == 'NFB0014', ] # Update as appropriate
probe14.arla <- probe14[probe14$spps == 'arla', ] # Update as appropriate
probe14.arla.22jun <- probe14.arla[probe14.arla$date.day == '22' &
                                     probe14.arla$date.month == 'jun' , ] # Update as appropriate
#probe14.olma.trial2.16c <- probe14.olma.trial2[probe14.olma.trial2$temp == '16', ] # Update as appropriate
head(probe14.arla.22jun) # Update as appropriate
str(probe14.arla.22jun) # Update as appropriate

### SMR estimate function

five.hr.plus.data <- probe14.arla.22jun[probe14.arla.22jun$time.hrs > 5, ]

smr <- calcSMR(five.hr.plus.data$mo2) # Update as appropriate
smr

smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best

  #############################

plot.smr <- ggplot(data=probe14.arla.22jun, aes(x=probe14.arla.22jun$time.hrs, 
                                                  y=probe14.arla.22jun$mo2)) 
plot.smr + (geom_point(size = 1.47)) +
  geom_hline(yintercept = 1.47) +
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

calcO2crit(pcrit.data, 1.47) # Enter value of SMR obtained above here, after "pcrit.data
#?calcO2crit

plotO2crit(calcO2crit(pcrit.data, 1.47))

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

pcrit.data <-read.csv(file.choose())
head(pcrit.data)
str(pcrit.data)

# Variable objects for use in plotting

fish.id.mass <- pcrit.data$fish.id.mass
spps <- pcrit.data$spps
pcrit.type <- pcrit.data$p.crit.type
pcrit.bestsmr <- pcrit.data$pcrit.bestsmr
pcrit.regress <- pcrit.data$p.crit.regress
smr.time <- pcrit.data$smr.time
smr.best <- pcrit.data$best.smr

##################################################################################
### Plot data to visualize whether "pcrit type" affects pcrit measurements.
##################################################################################

  ### Plot: Pcrit via "best smr" fishMO2 o2crit estimate

plot.pcrit.bestsmr <- ggplot(data=pcrit.data, aes(x=pcrit.type, 
                                   y=pcrit.bestsmr, 
                                   group=fish.id.mass, 
                                   shape=fish.id.mass, 
                                   colour=fish.id.mass)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) +
  labs(title="Pcrit via SMR estimate",
       x = "Pcrit type",
       y = "Pcrit (torr)")
plot.pcrit.bestsmr +
  scale_colour_discrete(name = "Fish ID",
                        breaks = c("clgl.big", "clgl.small", "olma.big"),
                        labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("clgl.big", "clgl.small", "olma.big"),
                       labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1")) +
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

  ###
  ### Plot: SMR estimates following 16hr vs 45hr measurement period - DIFFERENT TRIALS!
  ###

plot.smr.time <- ggplot(data=smr.16v45.data, aes(x=smr.time, 
                                                  y=best.smr, 
                                                  group=fish.id.mass[smr.time.cat$date == '6-May-17'], 
                                                  shape=fish.id.mass[smr.time.cat$date == '6-May-17'], 
                                                  colour=fish.id.mass[smr.time.cat$date == '6-May-17'])) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) +
  labs(title="SMR estimated after 16 or 45 hours of measurement",
       x = "Measurement time (hrs)",
       y = "SMR (umol O2/g/hr)")
plot.smr.time +
  scale_colour_discrete(name = "Fish ID",
                        breaks = c("clgl.big", "clgl.small", "olma.big"),
                        labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("clgl.big", "clgl.small", "olma.big"),
                       labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1")) +
  theme_base()

######################################
######################################
######################################
    ###
    ### OLMA trial 2: 18 - April - 2017
    ###
######################################
######################################
######################################

# Y values for plot - OLMA trial 2
smr.16.16plus.olma.trial2 <- smr.16v45.data$best.smr[smr.16v45.data$date == "18-Apr-17"]
# X values for plot - OLMA trial 2
time.category.smr.olma.trial2 <- smr.16v45.data$smr.time.cat[smr.16v45.data$date == "18-Apr-17"]
# Grouping variable - OLMA trial 2
fishid.forplot.olma.trial2 <- smr.16v45.data$fish.id.mass[smr.16v45.data$date == "18-Apr-17"]
######

df.16hr.16plus.olma.t2 <- data.frame(smr.16.16plus.olma.trial2,
                                     time.category.smr.olma.trial2,
                                     fishid.forplot.olma.trial2)

plot.smr.time <- ggplot(data=df.16hr.16plus.olma.t2, aes(x=time.category.smr.olma.trial2, 
                                             y=smr.16.16plus.olma.trial2, 
                                             group=fishid.forplot.olma.trial2, 
                                             shape=fishid.forplot.olma.trial2, 
                                             colour=fishid.forplot.olma.trial2)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) +
  labs(title="SMR estimated after 16\nor 16+ hours of measurement",
       x = "Measurement time (hrs)",
       y = "SMR (umol O2/g/hr)")
plot.smr.time +
  scale_colour_discrete(name = "Fish ID",
                        breaks = c("olma.big", "olma.med", "olma.small"),
                        labels = c("Tidepool 1", "Tidepool 2", "Tidepool 3")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("olma.big", "olma.med", "olma.small"),
                       labels = c("Tidepool 1", "Tidepool 2", "Tidepool 3")) +
  theme_base()

######################################
    ###
    ### END plot for OLMA trial 2: 18-Apr-2017
    ###
######################################

######################################
######################################
######################################
###
### OLMA trial 4/CLGL trial 2: 6 - May - 2017
###
######################################
######################################
######################################

# Y values for plot - OLMA trial 4/CLGL trial 2
smr.16.16plus.6may <- smr.16v45.data$best.smr[smr.16v45.data$date == "6-May-17"]
# X values for plot - OLMA trial 4/CLGL trial 2
time.category.smr.6may <- smr.16v45.data$smr.time.cat[smr.16v45.data$date == "6-May-17"]
# Grouping variable - OLMA trial 4/CLGL trial 2
fishid.forplot.6may <- smr.16v45.data$fish.id.mass[smr.16v45.data$date == "6-May-17"]
######

df.16hr.16plus.6may <- data.frame(smr.16.16plus.6may,
                                     time.category.smr.6may,
                                     fishid.forplot.6may)

plot.smr.time <- ggplot(data=df.16hr.16plus.6may, aes(x=time.category.smr.6may, 
                                                         y=smr.16.16plus.6may, 
                                                         group=fishid.forplot.6may, 
                                                         shape=fishid.forplot.6may, 
                                                         colour=fishid.forplot.6may)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) +
  labs(title="SMR estimated after 16\nor 16+ hours of measurement",
       x = "Measurement time (hrs)",
       y = "SMR (umol O2/g/hr)")
plot.smr.time +
  scale_colour_discrete(name = "Fish ID",
                        breaks = c("clgl.big", "clgl.small", "olma.big", "olma.small"),
                        labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1", "Tidepool 2")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("clgl.big", "clgl.small", "olma.big", "olma.small"),
                       labels = c("Mosshead 1", "Mosshead 2", "Tidepool 1", "Tidepool 2")) +
  theme_base()

######################################
###
### END plot for OLMA trial 4/CLGL trial 2: 6-May-2017
###
######################################

######################################
######################################
######################################
###
### ARHA trial 1: 5 - May - 2017
###
######################################
######################################
######################################

# Y values for plot - ARHA trial 1
smr.16.16plus.5may <- smr.16v45.data$best.smr[smr.16v45.data$date == "5-May-17"]
# X values for plot - ARHA trial 1
time.category.smr.5may <- smr.16v45.data$smr.time.cat[smr.16v45.data$date == "5-May-17"]
# Grouping variable - ARHA trial 1
fishid.forplot.5may <- smr.16v45.data$fish.id.mass[smr.16v45.data$date == "5-May-17"]
######

df.16hr.16plus.5may <- data.frame(smr.16.16plus.5may,
                                  time.category.smr.5may,
                                  fishid.forplot.5may)

plot.smr.time <- ggplot(data=df.16hr.16plus.5may, aes(x=time.category.smr.5may, 
                                                      y=smr.16.16plus.5may, 
                                                      group=fishid.forplot.5may, 
                                                      shape=fishid.forplot.5may, 
                                                      colour=fishid.forplot.5may)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) +
  labs(title="SMR estimated after 16\nor 16+ hours of measurement",
       x = "Measurement time (hrs)",
       y = "SMR (umol O2/g/hr)")
plot.smr.time +
  scale_colour_discrete(name = "Fish ID",
                        breaks = c("arha.1", "arha.2", "arha.3", "arha.4"),
                        labels = c("Scalyhead 1", "Scalyhead 2", "Scalyhead 3", "Scalyhead 4")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("arha.1", "arha.2", "arha.3", "arha.4"),
                       labels = c("Scalyhead 1", "Scalyhead 2", "Scalyhead 3", "Scalyhead 4")) +
  theme_base()

######################################
###
### END plot for ARHA trial 1: 5-May-2017
###
######################################

### Plot: Pcrit via "best smr" fishMO2 o2crit estimate: 16 hrs VS 16+ hrs

pcrit.data <-read.csv(file.choose())
head(pcrit.data)
str(pcrit.data)

# Variable objects for use in plotting

fish.id.mass <- pcrit.data$fish.id.mass
mass <- pcrit.data$mass
spps <- pcrit.data$spps
pcrit.calc <- pcrit.data$pcrit.calc
pcrit <- pcrit.data$pcrit

#pcrit.bestsmr.16hr <- pcrit.bestsmr[pcrit.data$smr.time.cat == "16hr"]
#pcrit.bestsmr.16hrplus <- pcrit.bestsmr[pcrit.data$smr.time.cat == "more16hr"]

plot.pcrit <- ggplot(data=pcrit.data, aes(x=pcrit.calc, 
                                                  y=pcrit, 
                                                  group=fish.id.mass, 
                                                  colour=fish.id.mass)) + 
  geom_line(size = 1.5) + 
  geom_point(size = 3) +
  labs(title="Pcrit via 16hrs SMR, 16+hrs SMR, REGRESS",
       x = "Pcrit calculation method",
       y = "Pcrit (torr)")
plot.pcrit +
  scale_colour_discrete(name = "Fish ID",
                        breaks = c("olma.big.1", "arha.1", "arha.2", "arha.4",
                                   "clgl.big", "clgl.small", "olma.big",
                                   "olma.small"),
                        labels = c("Tidepool 1\nday 1", "Scalyhead 1",
                                   "Scalyhead 2", "Scalyhead 4",
                                   "Mosshead 1", "Mosshead 2", 
                                   "Tidepool 1\n day 2", "Tidepool 2")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("olma.big.1", "arha.1", "arha.2", "arha.4",
                                  "clgl.big", "clgl.small", "olma.big",
                                  "olma.small"),
                       labels = c("Tidepool 1\nday 1", "Scalyhead 1",
                                  "Scalyhead 2", "Scalyhead 4",
                                  "Mosshead 1", "Mosshead 2", 
                                  "Tidepool 1\n day 2", "Tidepool 2")) +
  theme_base()


################################################
################################################
  ###
  ### Testing for differences in Pcrit technique 
  ### and SMR measurement period length : Paired t-test
  ###
################################################
################################################

pcrit.tech.data <-read.csv(file.choose())
pcrit.tech.data

# paired t-test: Pcrit technique, best-smr estimated Pcrit

pcrit.closed.bestsmr <- pcrit.tech.data$pcrit.bestsmr[pcrit.tech.data$p.crit.type == "closed"]
pcrit.semiclosed.bestsmr <- pcrit.tech.data$pcrit.bestsmr[pcrit.tech.data$p.crit.type == "semi-closed"]

t.test(pcrit.closed.bestsmr,pcrit.semiclosed.bestsmr,paired=TRUE) # where y1 & y2 are numeric

# paired t-test: Pcrit technique, REGRESS estimated Pcrit
pcrit.closed.regress <- pcrit.tech.data$p.crit.regress[pcrit.tech.data$p.crit.type == "closed"]
pcrit.semiclosed.regress <- pcrit.tech.data$p.crit.regress[pcrit.tech.data$p.crit.type == "semi-closed"]

t.test(pcrit.closed.regress,pcrit.semiclosed.regress,paired=TRUE) # where y1 & y2 are numeric

  ### paired t-test: 16 vs 16+ hour SMR estimates from same trial

smr.tech.data <-read.csv(file.choose())
smr.tech.data

SMR.16hr.stats <- smr.tech.data$best.smr[smr.tech.data$smr.time.cat == "16hr"]
SMR.16hrplus.stats <- smr.tech.data$best.smr[smr.tech.data$smr.time.cat == "more16hr"]
head(SMR.16hr.stats)
head(SMR.16hrplus.stats)

delta.smr <- SMR.16hr.stats - SMR.16hrplus.stats

qqnorm(delta.smr)
qqline(delta.smr,col="red")
shapiro.test(delta.smr)
boxplot(delta.smr)

log.delta.smr <- log(delta.smr + 1)
qqnorm(log.delta.smr)
qqline(log.delta.smr, col="red")
shapiro.test(log.delta.smr)
boxplot(delta.smr)

sqrt.delta.smr <- sqrt(delta.smr + 1)
qqnorm(sqrt.delta.smr)
qqline(sqrt.delta.smr, col="red")
shapiro.test(sqrt.delta.smr)
boxplot(sqrt.delta.smr)

wilcox.test(delta.smr)
# t.test(SMR.16hr.stats,SMR.16hrplus.stats,paired=TRUE) # where y1 & y2 are numeric

### paired t-test: 16 vs 16+ hour SMR-PCRIT estimates from same trial

pcrit.time.data <-read.csv(file.choose())
pcrit.time.data

pcrit.16hr.stats <- pcrit.time.data$pcrit.bestsmr[pcrit.time.data$smr.time.cat == "16hr"]
pcrit.16hrplus.stats <- pcrit.time.data$pcrit.bestsmr[pcrit.time.data$smr.time.cat == "more16hr"]
head(pcrit.16hr.stats)
head(pcrit.16hrplus.stats)

delta.pcrit.time <- pcrit.16hr.stats - pcrit.16hrplus.stats

qqnorm(delta.pcrit.time)
qqline(delta.pcrit.time,col="red")
shapiro.test(delta.pcrit.time)
boxplot(delta.pcrit.time)

log.delta.pcrit.time <- log(delta.pcrit.time + 1)
qqnorm(log.delta.pcrit.time)
qqline(log.delta.pcrit.time, col="red")
shapiro.test(log.delta.pcrit.time)
boxplot(log.delta.pcrit.time)

#sqrt.delta.smr <- sqrt(delta.pcrit.time + 1)
#qqnorm(sqrt.delta.smr)
#qqline(sqrt.delta.smr, col="red")
#shapiro.test(sqrt.delta.smr)
#boxplot(sqrt.delta.smr)

#wilcox.test(delta.smr)
t.test(log.delta.pcrit.time) # where y1 & y2 are numeric
t.test()

######################################
######################################
######################################
###
### OLMA trial 1: 18 - Apr - 2017 MO2 values over 70 hour SMR measurment period
###
######################################
######################################
######################################

## dataset for "Tidepool 1": NFB0014

head(probe12.olma.trial4) ## Need to update this object at head of script for each fish
str(probe12.olma.trial4)

mo2.all <- probe12.olma.trial4$mo2
time.all <- probe12.olma.trial4$time.hrs
time.16VS16plus <- ifelse(probe12.olma.trial4$time.hrs < 16, "16hr", "16hrPlus")

df.mo2.time.all <- data.frame(mo2.all,
                              time.all,
                              time.16VS16plus)


plot.mo2.all <- ggplot(data=df.mo2.time.all, aes(x=time.all, 
                                                 y=mo2.all, 
                                                 colour=time.16VS16plus)) + 
  geom_point(size = 2) +
  geom_smooth()+
  labs(title = "45 hour measurement period: Tidepool 2",
       x = "Measurement time (hrs)",
       y = "MO2 (umol O2/g/hr)")
plot.mo2.all +
  scale_colour_discrete(name = "Measurement\nperiod",
                        breaks = c("16hr", "16hrPlus"),
                        labels = c("16 hours", "> 16 hours")) +
  scale_shape_discrete(name = "Fish ID",
                       breaks = c("16hr", "16hrPlus"),
                       labels = c("16 hours", "> 16 hours")) +
  theme_base()

######################################
###
### END MO2 vs time plot: "Tidepool 1" NFB0014 18-Apr-2017
###
######################################
