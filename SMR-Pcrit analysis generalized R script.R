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
probe12.olma <- probe12[probe12$spps == 'olma', ] # Update as appropriate
probe12.olma.trial3 <- probe12.olma[probe12.olma$trial.no == '3', ] # Update as appropriate
head(probe12.olma.trial3) # Update as appropriate
str(probe12.olma.trial3) # Update as appropriate

### SMR estimate function

smr <- calcSMR(probe12.olma.trial3$mo2) # Update as appropriate
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

calcO2crit(pcrit.data, 1.32) # Enter value of SMR obtained above here, after "pcrit.data
#?calcO2crit

plotO2crit(calcO2crit(pcrit.data, 1.32))

### In torr
#> (O2crit.%sat/100)*P.ATM.KPA*760*0.2095/101.325

