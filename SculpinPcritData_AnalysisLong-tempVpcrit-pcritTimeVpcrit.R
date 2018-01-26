library(tidyverse)

## Code R studio uses to bring in a dataset:
#library(readr)

md <- read_csv("SculpinPcritData_DealingWithData.csv")
View(md)

# Visualize pcrit vs pcrit trial time: Is there an effect?
md.ptimes <- md %>%
  select(spps, temp, fish.id, pcrit.r, trial.no, pcrit.time.min) %>%
  filter(trial.no == 1, pcrit.r != "na", pcrit.time.min != "na") %>%
  group_by(spps,temp)

# Plot of pcrit vs pcrit time by spps
par(mfrow = c(3,3))
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="olma"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="olma"],
     xlab = "time", ylab = "Pcrit", main = "olma")
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="clgl"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="clgl"],
     xlab = "time", ylab = "Pcrit", main = "clgl")
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="arha"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="arha"],
     xlab = "time", ylab = "Pcrit", main = "arha")
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="arla"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="arla"],
     xlab = "time", ylab = "Pcrit", main = "arla")
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="arfe"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="arfe"],
     xlab = "time", ylab = "Pcrit", main = "arfe")
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="blci"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="blci"],
     xlab = "time", ylab = "Pcrit", main = "blci")
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="hehe"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="hehe"],
     xlab = "time", ylab = "Pcrit", main = "hehe")
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="scma"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="scma"],
     xlab = "time", ylab = "Pcrit", main = "scma")
plot(md.ptimes$pcrit.r[md.ptimes$temp==12&md.ptimes$spps=="enbi"]~
       md.ptimes$pcrit.time.min[md.ptimes$temp==12&md.ptimes$spps=="enbi"],
     xlab = "time", ylab = "Pcrit", main = "enbi")

