getwd()

#install.packages("fishMO2.tar.gz")

#install.packages("mclust", dependencies = TRUE)
#install.packages("shape", dependencies = TRUE)
#install.packages("StreamMetabolism", dependencies = TRUE)

library(mclust)
library(shape)
library(StreamMetabolism)

library(fishMO2)

#?fishMO2

smr.data <-read.csv(file.choose())
head(smr.data)

### NFB0009: Tipepool sculpin (OLMA) #1: Fin clip ID = Top clip

probe9 <- smr.data[smr.data$probe == 'NFB0009', ]
str(probe9)

  ### Leaving off data recorded at less than 15 hrs in respirometry trial

probe9.15plus <- probe9[probe9$time.hrs > 15, ]
probe9.15plus
head(probe9.15plus)
str(probe9.15plus)

smr <- calcSMR(probe9.15plus$mo2)
smr
smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best

### NFB0012: Tipepool sculpin (OLMA) #2: Fin clip ID = Low clip

probe12 <- smr.data[smr.data$probe == 'NFB0012', ]
str(probe12)

### Leaving off data recorded at less than 15 hrs in respirometry trial

probe12.15plus <- probe12[probe12$time.hrs > 15, ]
probe12.15plus
head(probe12.15plus)
str(probe12.15plus)

smr.12 <- calcSMR(probe12.15plus$mo2)
smr.12
smr.12.check.best <- as.numeric(ifelse(smr.12$CVmlnd > 5.4, smr.12$quant[4], smr.12$smr.12mlnd)) # as recommended in Chabot et al. 2016
smr.12.check.best

### NFB0014: Tipepool sculpin (OLMA) #3

probe14 <- smr.data[smr.data$probe == 'NFB0014', ]
str(probe14)

### Leaving off data recorded at less than 15 hrs in respirometry trial

probe14.15plus <- probe14[probe14$time.hrs > 15, ]
probe14.15plus
head(probe14.15plus)
str(probe14.15plus)

smr.14 <- calcSMR(probe14.15plus$mo2)
smr.14
smr.14.check.best <- as.numeric(ifelse(smr.14$CVmlnd > 5.4, smr.14$quant[4], smr.14$smr.14mlnd)) # as recommended in Chabot et al. 2016
smr.14.check.best

plotMO2fdis(probe14.15plus$mo2)

  ### Pcrit data: Semi-closed Pcrit for "No-Clip" OLMA

pcrit.data <-read.csv(file.choose())
head(pcrit.data)

calcO2crit(pcrit.data, 1.94)
?calcO2crit

noclip.persat <- pcrit.data$percent.sat
noclip.mo2 <- pcrit.data$mo2

###########################################
### ARHA #1
###########################################

smr.data <-read.csv(file.choose())
head(smr.data)

### NFB0009: Scalyhead sculpin (ARHA) #1: Fin clip ID = Top clip

probe9 <- smr.data[smr.data$probe == 'NFB0009', ]
str(probe9)

### Leaving off data recorded at less than 15 hrs in respirometry trial

probe9.15plus <- probe9[probe9$time.hrs > 15, ]
probe9.15plus
head(probe9.15plus)
str(probe9.15plus)

smr <- calcSMR(probe9.15plus$mo2)
smr
smr.check.best <- as.numeric(ifelse(smr$CVmlnd > 5.4, smr$quant[4], smr$mlnd)) # as recommended in Chabot et al. 2016
smr.check.best

### Pcrit data: Semi-closed Pcrit for "Top-Clip" ARHA

pcrit.data <-read.csv(file.choose())
head(pcrit.data)

calcO2crit(pcrit.data, 1.67)
 ?calcO2crit

do <- pcrit.data$DO
mo2<- pcrit.data$MO2

plot(do, mo2)
