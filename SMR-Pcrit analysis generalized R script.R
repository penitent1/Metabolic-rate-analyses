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

probe09 <- smr.data[smr.data$probe == 'NFB0009', ] # Update as appropriate
probe09.clgl <- probe09[probe09$spps == 'clgl', ] # Update as appropriate
probe09.clgl.trial1 <- probe09.clgl[probe09.clgl$trial.no == '1', ] # Update as appropriate
head(probe09.clgl.trial1) # Update as appropriate
str(probe09.clgl.trial1) # Update as appropriate

### SMR estimate function

smr <- calcSMR(probe09.clgl.trial1$mo2) # Update as appropriate
smr
#> smr
#$mlnd
#1 
#2.273519

## $CVmlnd
## 1 
## 7.251364
smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best
## > smr.check.best
## [1] 2.142857

### Estimate Pcrit using the SMR obtained above

pcrit.data <-read.csv(file.choose())
head(pcrit.data)

calcO2crit(pcrit.data, 2.14) # Enter value of SMR obtained above here, after "pcrit.data
#?calcO2crit

plotO2crit(calcO2crit(pcrit.data, 2.14))

### In torr
#> (O2crit.%sat/100)*P.ATM.KPA*760*0.2095/101.325

