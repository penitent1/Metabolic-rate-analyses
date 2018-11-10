## Non-linear regression (NLR) a la Marshall et al. 2013 JEB
library(tidyverse)
library(fishMO2)
## Pcrit: 19 Oct 2018, Tidepool sculpin at 12 degrees C

## gnls() not a recognized function?

pd <- read_csv("C:/Users/derek/Documents/Metabolic-rate-analyses/Ken Chu Pcrit project/Ken_Pcrit_OLMA_TopFinClip_Attempt2PostCompFreeze_19Oct2018.csv")

ggplot(pd, aes(x = DO, y = MO2)) +
  geom_point()

pd_trunc <- pd %>% filter(DO < 75)

plotO2crit(calcO2crit(pd, SMR = 4.6, lowestMO2 = 4.6))

nls(MO2 ~ a*(1-exp((-1*(DO/b))^c)), data = pd, start = list(a = 6, b = 5, c = 1))

x <- seq(0,100,2)
a <- 6
b <- 5
c <- 1
y <- a*(1-exp((-1*(x/b))^c))
nlrd <- tibble(x = x, y = y)
ggplot(nlrd, aes(x = x, y = y)) +
  geom_point()
