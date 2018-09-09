library(tidyverse)
library(lubridate)
#install.packages("devtools")

#devtools::install_github("januarharianto/respR")

library(respR)

getwd()
setwd("C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project")

urchin_int <- inspect_data(intermittent.rd)
head(intermittent.rd)

### Probe drift correction data

pre <- c(0,100)
post <- c(0,99.4) # Update for each trial
probeDrift <- lm(post~pre)
probeDrift ## "pre" is the slope of the line

alpha <- 1.7206 ## At 12 C and 35 ppt salinity
body_mass <- 4.4 # grams >> 0.0044 "L"
respVols <- read.csv(file.choose())

### MO2 data
mo2 <- read_csv(file.choose()) %>%
  dplyr::select(-c(`Tau - Phase Method`,`Sensor Temperature`)) %>%
  mutate(resp_id = "b",
         resp_vol = (respVols$resp.volume_mL[respVols$respirometer.id == "b"])/1000, # Update each trial
         time = seq(0, (nrow(mo2)-1)*15, by=15),
         po2 = (((Oxygen)/probeDrift$coefficients[2])/100)*`Air Pressure`*760*0.2095/101.325,
         umol = po2*alpha*(resp_vol-(body_mass/1000))) %>%
  filter(Time < as_datetime("2018-07-05 11:05:05"))
  

# Plot of raw Time (dttm) vs Oxygen (%)
ggplot(mo2, aes(x = time, y = umol)) +
  geom_point()

## Get MO2's!

wsize <- 20

output <- numeric(150)
check <- numeric(150)
#meano2 <- numeric(150)
glms <- list()
j <- 1

i <- 1
while (i <= nrow(mo2)-wsize)
{
  time <- mo2$time[i:(i+wsize)]
  umol <- mo2$umol[i:(i+wsize)]
  z <- lm(umol~time)
  
  if (summary(z)$adj.r.squared>=0.9 & coef(z)[2]<0) #  & coef(z)[2]>-0.1 
  {
    output[j] <- coef(z)[2]
    glms[[j]] <- z
    check[j] <- i
    #meano2[j] <- mean(o2)
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
output <- output*-1*3600/body_mass
#meano2 <- meano2[meano2!=0] ## NOTE These are umol units - convert to PO2
#meano2 <- meano2/alpha/(respVol/1000-respVol/1000) # This converts back to torr
check <- check[check!=0]

### extend below so the last plus value matches 'wsize'
regRegions <- c(check, check+1, check+2, check+3, check+4, check+5, check+6, check+7, check+8, check+9, check+10, check+11, check+12, check+13, check+14, check+15, check+16, check+17, check+18, check+19, check+20)
regRegions <- regRegions[order(regRegions)]

plot(mo2$time, mo2$umol, pch=20, cex=0.67, xlab="Time (sec)", ylab="Oxygen concentratio\n (umol O2)")
points(mo2$time[check], mo2$umol[check], pch=1, col="darkorange", cex=1.2)
segments(mo2$time[check],0,mo2$time[check],110,col="grey67",lty=2)

plot(mo2$time[regRegions],mo2$umol[regRegions], pch=20, cex=0.67, ylim=c(40,70), xlab="Time", ylab="[oxygen] (umol)")

for (i in 1:length(glms))
{
  lines(mo2$time[(check[i]-10):(check[i]+30)],predict(glms[[i]], newdata=mo2[(check[i]-10):(check[i]+30),]),col="magenta",lty=2)
}

ind <- numeric(length(check)-1)
for (i in 1:length(ind))
{
  ind[i] <- check[i+1]-check[i]
}

ind2 <- logical(length(ind)-4)
for (i in 1:length(ind))
{
  ind2[i] <- (ind[i]==wsize & ind[i+1]==wsize & ind[i+2]==wsize & ind[i+3]==wsize)
}

rmr <- output[1:(grep(TRUE, ind2)[1]-1)]
extraRmrCheck <- ind[1:(grep(TRUE, ind2)[1]-1)]
rmr <- rmr[extraRmrCheck!=wsize]
