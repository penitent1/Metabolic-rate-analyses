library(ggplot2)

mydata <-read.csv(file.choose())
head(mydata)
str(mydata)

time <- mydata$time.min
percent.sat <- mydata$sat
baro.kpa <- mydata$pressure.kpa
baro.torr <- (760/101.325)*baro.kpa
po2.torr <- baro.torr*.2095*(percent.sat/100)
mydata$po2.torr <- po2.torr
head(mydata)

plot(po2.torr ~ time)



