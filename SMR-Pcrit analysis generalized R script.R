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

probe10 <- smr.data[smr.data$probe == 'NFB0010', ] # Update as appropriate
probe10.olma <- probe10[probe10$spps == 'olma', ] # Update as appropriate
probe10.olma.trial4 <- probe10.olma[probe10.olma$trial.no == '4', ] # Update as appropriate
head(probe10.olma.trial4) # Update as appropriate
str(probe10.olma.trial4) # Update as appropriate

### SMR estimate function

smr <- calcSMR(probe10.olma.trial4$mo2) # Update as appropriate
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

calcO2crit(pcrit.data, 2.13) # Enter value of SMR obtained above here, after "pcrit.data
#?calcO2crit

plotO2crit(calcO2crit(pcrit.data, 2.13))

### In torr
#> (O2crit.%sat/100)*P.ATM.KPA*760*0.2095/101.325

