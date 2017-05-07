library(ggplot2)
library(knitr)

smr.data <-read.csv(file.choose())
head(smr.data)

smr.data$time.bin <- factor(smr.data$time.bin)
str(smr.data)

time.min <- smr.data$time.min
mo2 <- smr.data$mo2
time.bin <- smr.data$time.bin
probe <- smr.data$probe

### NFB0009

probe9 <- smr.data[smr.data$probe == 'NFB0009', ]
str(probe9)

  ## Linear model: continuous, test hypothesis: regression slope is not different
  ## from zero AFTER (including) 15 HR time bin

probe9.15plus <- probe9[probe9$time.hrs > 15, ]
probe9.15plus
head(probe9.15plus)
str(probe9.15plus)

lm.p9plus15 <- lm(probe9.15plus$mo2 ~ probe9.15plus$time.hrs, probe9.15plus)
summary(lm.p9plus15)
anova(lm.p9plus15)

#time.bin.vs.mo2 <- lm(mo2 ~ time.bin, smr.data)
#summary(time.bin.vs.mo2)
#anova(time.bin.vs.mo2)

#a1 <- aov(mo2 ~ time.bin)
#posthoc.probe9 <- TukeyHSD(x=a1, 'time.bin', conf.level=0.95)
#posthoc.probe9
#plot(time.bin.vs.mo2)

### NFB0012

probe12 <- smr.data[smr.data$probe == 'NFB0012', ]
str(probe12)

## Linear model: continuous, test hypothesis: regression slope is not different
## from zero AFTER (including) 15 HR time bin

probe12.15plus <- probe12[probe12$time.hrs > 15, ]
probe12.15plus
head(probe12.15plus)
str(probe12.15plus)

lm.p12plus15 <- lm(probe12.15plus$mo2 ~ probe12.15plus$time.hrs, probe12.15plus)
summary(lm.p12plus15)
anova(lm.p12plus15)


#plot(probe12$time.hrs, probe12$mo2)
#plot(probe12$time.bin, probe12$mo2)
#probe12$time.hrs
#p12.time.v.mo2 <- lm(probe12$mo2 ~ probe12$time.bin, smr.data)
#summary(p12.time.v.mo2)
#anova(p12.time.v.mo2)

#a1.probe9 <- aov(probe12$mo2 ~ probe12$time.bin)
#posthoc <- TukeyHSD(x=a1.probe9, 'time.bin', conf.level=0.95)
#posthoc

### NFB0014

probe14 <- smr.data[smr.data$probe == 'NFB0014', ]
str(probe14)

## Linear model: continuous, test hypothesis: regression slope is not different
## from zero AFTER (including) 15 HR time bin

probe14.15plus <- probe14[probe14$time.hrs > 15, ]
probe14.15plus
head(probe14.15plus)
str(probe14.15plus)

lm.p9plus15 <- lm(probe9.15plus$mo2 ~ probe9.15plus$time.hrs, probe9.15plus)
summary(lm.p9plus15)
anova(lm.p9plus15)

#plot(probe14$time.hrs, probe14$mo2)
#plot(probe14$time.bin, probe14$mo2)
#probe14$time.hrs
#p14.time.v.mo2 <- lm(probe14$mo2 ~ probe14$time.bin, smr.data)
#summary(p14.time.v.mo2)
#anova(p14.time.v.mo2)

#a1.probe14 <- aov(probe14$mo2 ~ probe14$time.bin)
#posthoc.probe14 <- TukeyHSD(x=a1.probe14, 'time.bin', conf.level=0.95)
#posthoc.probe14

