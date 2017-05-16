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
head(probe09.clgl) # Update as appropriate
str(probe09.clgl) # Update as appropriate

### Leaving off data recorded at less than 15 hrs in respirometry trial

#probe09.15plus.clgl <- probe09.clgl[probe09.clgl$time.hrs > 15, ] # Update as appropriate
#probe09.15plus.clgl # Update as appropriate
#head(probe09.15plus.clgl) # Update as appropriate
#str(probe09.15plus.clgl) # Update as appropriate

### SMR estimate function

smr <- calcSMR(probe09.clgl$mo2) # Update as appropriate
smr
#> smr
#$mlnd
#1 
#2.113142

## $CVmlnd
## 1 
## 9.648905
smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best
## > smr.check.best
## [1] 1.955556

### Estimate Pcrit using the SMR obtained above

pcrit.data <-read.csv(file.choose())
head(pcrit.data)

calcO2crit(pcrit.data, 1.96) # Enter value of SMR obtained above here, after "pcrit.data
#?calcO2crit

plotO2crit(calcO2crit(pcrit.data, 1.96))

#> calcO2crit(pcrit.data, 2.20)
#$o2crit
#[1] 18.2 >> 18.2% Air Sat USING!! ALL data (not post-2sigma analysis)
### In torr
#> (18.2/100)*102.1*760*0.2095/101.325
#[1] 29.19968

#> calcO2crit(pcrit.data, 1.96)
#$o2crit
#[1] 18 >> 18% Air Sat USING!! ALL data (not post-2sigma analysis)
### In torr
#> (18/100)*102.1*760*0.2095/101.325
#[1] 28.87881

#> calcO2crit(pcrit.data, 2.20)
#$o2crit
#[1] 15.1 >> 15.1% Air sat using 2-sigma analysis data - no difference!

### In torr
#> (15.1/100)*102.4*760*0.2095/101.325
#[1] 24.29729
