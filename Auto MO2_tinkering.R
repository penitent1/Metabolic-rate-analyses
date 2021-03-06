## My sculpin data
md <- read.csv(file.choose())
respVols <- read.csv(file.choose())

### 'md' is the data frame

### Probe drift correction data
  ## NEW NOTE: I should make a new datafile for the probe drift data
    ## EG 
pre <- c(0,100)
post <- c(0,100.8) # UPDATE for each trial
probeDrift <- lm(post~pre)
probeDrift ## "pre" is the slope of the line -> coef(probeDrift)[2] == "pre"

wsize <- 20 ## 20 rows * 15 sec/row = 300 sec / 60 sec / min = 5 min windows

respVol <- respVols$resp.volume_mL[respVols$respirometer.id == "a"] # Update each trial
mass <- 18.4 # Update each trial
alpha <- 1.4842 ## Change as appropriate for different temperatures
pAtm <- 101.9 # Update each trial: atmospheric pressure

md$time <- seq(0, (nrow(md)-1)*15, by=15)
# Update "po2" for each trial
md$po2 <- (((md$Oxygen)/(coef(probeDrift)[2]))/100)*pAtm*760*0.2095/101.325
# Update "oxygen.umol" for each trial
md$oxygen.umol <- (((md$Oxygen)/(coef(probeDrift)[2]))/100)*pAtm*760*0.2095*alpha*(respVol/1000-mass/1000)/101.325 ### convert from %sat to umol O2


for (i in 1:(nrow(md)-wsize*5)) # Mult by 10 to extend window to 50 minutes (all Pcrits are at least an hour long)
{
  time <- md$time[i:(i+wsize*5)]
  o2 <- md$oxygen.umol[i:(i+wsize*5)]
  z <- lm(o2~time)
  
  if (summary(z)$adj.r.squared > 0.975 & coef(z)[2] < 0) 
  {
    pcritCaliData <- md[i:nrow(md), 6:8]
    rmrData <- md[1:(i-1), c(6,8)]
    break()
  } 
}

for (i in 1:(nrow(pcritCaliData)-10))
{
  if (all(pcritCaliData[i:(i+10),3]==((pcritCaliData[i:(i+10),3])[order(pcritCaliData[i:(i+10),3])])))
  {
    caliData <- pcritCaliData[(i-7):nrow(pcritCaliData),]
    pcritData <- pcritCaliData[1:(i-8),]
    break()
  }
}

# Reset time to start at zero
pcritData$time <- seq(0, (nrow(pcritData)-1)*15, by=15)
# plot Pcrit with all data included
plot(pcritData$oxygen.umol~pcritData$time, cex = 0.25)
# plot Pcrit with first 5 min of data removed (respirometer still mixing - MAY REQUIRE PLAYING)
plot(pcritData$oxygen.umol[pcritData$time>300]~pcritData$time[pcritData$time>300], cex = 0.25)

### Get slopes & Po2s for Pcrit

pcritData <- pcritData[pcritData$time>300,] # Use the data AFTER the mixing reaches equilibrium
  ## Determine how much time needs to be lopped off visually for now,
  ## maybe figure out a consistent criterion later

pcritWsize <- 7 ## CHANGE depending on desired window size ie 7 = 20degC, 15 =< 16degC

pcritZ <- list()
pcritSlopes <- numeric(200)
pcritPo2s <- numeric(200)

k <- 1
j <- 1

while (k <= (nrow(pcritData)-pcritWsize))
{
  pcritO2umol <- pcritData$oxygen.umol[k:(k+pcritWsize)]
  pcritTime <- pcritData$time[k:(k+pcritWsize)]
  pcritPo2 <-  pcritData$po2[k:(k+pcritWsize)]
  pcritZ <- lm(pcritO2umol~pcritTime)
  pcritSlopes[j] <- coef(pcritZ)[2]
  pcritPo2s[j] <- mean(pcritPo2)
  k <- k+pcritWsize
  j <- j+1
}

pcritSlopes <- pcritSlopes[pcritSlopes!=0]
pcritPo2s <- pcritPo2s[pcritPo2s!=0]

pcrit.df <- data.frame(pcritSlopes, pcritPo2s)

# plot Pcrit MO2's vs PO2's
plot((-1*pcrit.df$pcritSlopes*3600/mass)~pcrit.df$pcritPo2s, 
     cex = 1.5,
     xlab = "PO2 (torr)",
     ylab = "MO2 (umol O2/g/hr)")


##########################################
##########################################
##
## END PCRIT DATA DIVIONS FROM MIN MO2 DATA
##
##########################################
##########################################

head(rmrData) ## These are the data for Min MO2 assessment


output <- numeric(150)
check <- numeric(150)
glms <- list()
j <- 1

i <- 1
while (i <= nrow(rmrData)-wsize)
{
  time <- rmrData$time[i:(i+wsize)]
  o2 <- rmrData$oxygen.umol[i:(i+wsize)]
  z <- lm(o2~time)
  
  if (summary(z)$adj.r.squared>=0.97 & coef(z)[2]<0) #  & coef(z)[2]>-0.1 
  {
    output[j] <- coef(z)[2]
    glms[[j]] <- z
    check[j] <- i
    j <- j+1
    i <- i+wsize
  } 
  else 
  {
    i <- i+1
  }
}
output <- output[output!=0]
### could modify output to final MO2 units here if you wanted
output <- output*-1*3600/mass ## units = umol O2/g/hr
#meano2 <- meano2[meano2!=0] ## NOTE These are umol units - convert to PO2
#meano2 <- meano2/alpha/(respVol/1000-respVol/1000) # This converts back to torr
check <- check[check!=0]

### extend below so the last plus value matches 'wsize'
regRegions <- c(check, check+1, check+2, check+3, check+4, check+5, check+6, check+7, check+8, check+9, check+10, check+11, check+12, check+13, check+14, check+15, check+16, check+17, check+18, check+19, check+20)
regRegions <- regRegions[order(regRegions)]

plot(rmrData$time, rmrData$oxygen, pch=20, cex=0.67, xlab="Time (sec)", ylab="Oxygen concentratio\n (umol O2)")
points(rmrData$time[check], rmrData$oxygen[check], pch=1, col="darkorange", cex=1.2)
segments(rmrData$time[check],0,rmrData$time[check],110,col="grey67",lty=2)

## WHERE I LEFT OFF!