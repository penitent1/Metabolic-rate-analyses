## My sculpin data
md <- read.csv(file.choose())
respVols <- read.csv(file.choose())

### 'md' is the data frame

### Probe drift correction data

pre <- c(0,100)
post <- c(0,100.8) # Update for each trial
probeDrift <- lm(post~pre)
probeDrift ## "pre" is the slope of the line

wsize <- 20

respVol <- respVols$resp.volume_mL[respVols$respirometer.id == "e"] # Update each trial
mass <- 48.2 # Update each trial
alpha <- 1.4842 ## Change as appropriate for different temperatures

md$time <- seq(0, (nrow(md)-1)*15, by=15)
# Update "po2" for each trial
md$po2 <- (((md$Oxygen)/1.008)/100)*101.6*760*0.2095/101.325
# Update "oxygen.umol" for each trial
md$oxygen.umol <- (((md$Oxygen)/1.008)/100)*101.6*760*0.2095*alpha*(respVol/1000-mass/1000)/101.325 ### convert from %sat to umol O2

o2.rmr <- numeric(150000)
o2.pcrit <- numeric(150000)
po2.rmr <- numeric(150000)
po2.pcrit <- numeric(150000)
time.rmr <- numeric(150000)
time.pcrit <- numeric(150000)
j <- 1
k <- 1

i <- 1
while (i <= nrow(md)-wsize*5) # Mult by 10 to extend window to 50 minutes (all Pcrits are at least an hour long)
{
  time <- md$time[i:(i+wsize*5)]
  o2 <- md$oxygen.umol[i:(i+wsize*5)]
  z <- lm(o2~time)
  
  if (summary(z)$adj.r.squared<0.1 & coef(z)[2]>0) 
  {
    o2.rmr[j] <- md$oxygen.umol[md$time==[i]] #coef(z)[2]
    po2.rmr[j] <- md$po2[md$time==[i]]
    j <- j+1
    i <- i+1
  } 
  else 
  {
    o2.pcrit <- md$oxygen.umol[[i]:nrow(md)]
    #k <- k+1
    #i <- i+1
  }
}
output <- output[output!=0]
### could modify output to final MO2 units here if you wanted
output <- output*-1*3600/mass
meano2 <- meano2[meano2!=0] ## NOTE These are umol units - convert to PO2
meano2 <- meano2/alpha/(respVol/1000-respVol/1000) # This converts back to torr
check <- check[check!=0]

### extend below so the last plus value matches 'wsize'
regRegions <- c(check, check+1, check+2, check+3, check+4, check+5, check+6, check+7, check+8, check+9, check+10, check+11, check+12, check+13, check+14, check+15, check+16, check+17, check+18, check+19, check+20)
regRegions <- regRegions[order(regRegions)]

plot(md$time, md$oxygen, pch=20, cex=0.67, xlab="Time (sec)", ylab="Oxygen concentratio\n (umol O2)")
points(md$time[check], md$oxygen[check], pch=1, col="darkorange", cex=1.2)
segments(md$time[check],0,md$time[check],110,col="grey67",lty=2)
