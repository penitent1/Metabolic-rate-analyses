library(mclust)
library(shape)
library(StreamMetabolism)
library(fishMO2)

## Functions requiring fixin
?calcO2crit()

## Fixin this damn pcrit calculation function

pcrit.data <-read.csv(file.choose()) ## Working on 12-Sep 2017 NFB0008 SCMA
pcrit.data
str(pcrit.data)
nrow(pcrit.data) ## 26

calcO2crit(pcrit.data, SMR = 3.90) ## 3.90 
plotO2crit(calcO2crit(pcrit.data, SMR = 3.90))

## Chabot: IF lowest MO2 not directly stated, calculated as
## 5th percentile lowest MO2 from MO2 data calculated at DO above 80%
lowestMO2_check <- quantile(pcrit.data$MO2[pcrit.data$DO >= 80], p = 0.05)
lowestMO2_check ## 3.241269 = 5% quantile
length(pcrit.data$MO2[pcrit.data$DO >= 80]) ## 4

## Define whether an MO2 is >= "lowest MO2" == TRUE
## OR MO2 < "lowest MO2" == FALSE
geqSMR_check <- pcrit.data$MO2 >= lowestMO2_check
geqSMR_check
length(geqSMR_check) ## 26

## Minimum DO where MO2 is >= the 5th percentile lowest MO2 from quantile
## function estimate - so if I have a lowest MO2 estimate from the SMR
## data collection period I SHOULD DEFINITELY USE THE "lowestMO2 = "
## ARGUMENT!
pivotDO_check <- min(pcrit.data$DO[geqSMR_check])
pivotDO_check



calcO2crit_somo <- function (Data, SMR, lowestMO2 = NA, gapLimit = 4, max.nb.MO2.for.reg = 20) 
{
  Data = na.omit(Data)
  method = "LS_reg"
  if (is.na(lowestMO2)) 
    lowestMO2 = quantile(Data$MO2[Data$DO >= 80], p = 0.05)
  geqSMR = Data$MO2 >= lowestMO2
  pivotDO = min(Data$DO[geqSMR])
  lethal = Data$DO < pivotDO
  N_under_SMR = sum(lethal)
  final_N_under_SMR = lethal
  lastMO2reg = Data$MO2[Data$DO == pivotDO]
  if (N_under_SMR > 1) 
    theMod = lm(MO2 ~ DO, data = Data[lethal, ])
  if (N_under_SMR < 3) {
    missing = 3 - sum(lethal)
    not.lethal = Data$DO[geqSMR]
    DOlimit = max(sort(not.lethal)[1:missing])
    addedPoints = Data$DO <= DOlimit
    lethal = lethal | addedPoints
    theMod = lm(MO2 ~ DO, data = Data[lethal, ])
  }
  if (N_under_SMR >= 3) {
    lethalB = Data$DO <= pivotDO
    regA = theMod
    regB = lm(MO2 ~ DO, data = Data[lethalB, ])
    large_slope_drop = (coef(regA)[2]/coef(regB)[2]) > 1.1
    large_DO_gap = (max(Data$DO[lethalB]) - max(Data$DO[lethal])) > 
      gapLimit
    tooSmallMO2 = lastMO2reg < SMR
    if (!large_slope_drop & !large_DO_gap & !tooSmallMO2) {
      lethal = lethalB
      theMod = regB
    }
  }
  if (!is.na(max.nb.MO2.for.reg) & sum(lethal) > max.nb.MO2.for.reg) {
    Ranks = rank(Data$DO)
    lethal = Ranks <= max.nb.MO2.for.reg
    theMod = lm(MO2 ~ DO, data = Data[lethal, ])
    final_N_under_SMR = max.nb.MO2.for.reg
  }
  predMO2 = as.numeric(predict(theMod, data.frame(DO = Data$DO)))
  Data$delta = (Data$MO2 - predMO2)/predMO2 * 100
  Data$delta[Data$DO < pivotDO | lethal] = 0
  tol = 0
  HighValues = Data$delta > tol
  Ranks = rank(-1 * Data$delta)
  HighMO2 = HighValues & Ranks == min(Ranks)
  if (sum(HighValues) > 0) {
    nblethal = sum(lethal)
    Data$W = NA
    Data$W[lethal] = 1/nblethal
    Data$W[HighMO2] = 1
    theMod = lm(MO2 ~ DO, weight = W, data = Data[lethal | 
                                                    HighMO2, ])
    predMO2_2 = as.numeric(predict(theMod, data.frame(DO = Data$DO)))
    Data$delta2 = (Data$MO2 - predMO2_2)/predMO2_2 * 100
    Data$delta2[Data$DO < pivotDO] = 0
    tol = Data$delta2[HighMO2]
    HighValues2 = Data$delta2 > tol
    if (sum(HighValues2) > 0) {
      Ranks2 = rank(-1 * Data$delta2)
      HighMO2_2 = HighValues2 & Ranks2 == 1
      nblethal = sum(lethal)
      Data$W = NA
      Data$W[lethal] = 1/nblethal
      Data$W[HighMO2_2] = 1
      theMod2 = lm(MO2 ~ DO, weight = W, data = Data[lethal | 
                                                       HighMO2_2, ])
      if (theMod2$coef[2] > theMod$coef[2]) {
        theMod = theMod2
        HighMO2 = HighMO2_2
      }
    }
  }
  Coef = coefficients(theMod)
  AboveOrigin = Coef[1] > 0
  if (AboveOrigin) {
    theMod = lm(MO2 ~ DO - 1, data = Data[lethal, ])
    Coef = c(0, coefficients(theMod))
    method = "through_origin"
    HighMO2 = rep(FALSE, nrow(Data))
  }
  po2crit = as.numeric(round((SMR - Coef[1])/Coef[2], 1))
  sum_mod = summary(theMod)
  anov_mod = anova(theMod)
  O2CRIT = list(o2crit = po2crit, SMR = SMR, Nb_MO2_conforming = N_under_SMR, 
                Nb_MO2_conf_used = final_N_under_SMR, High_MO2_required = sum(HighMO2) == 
                  1, origData = Data, Method = method, mod = theMod, 
                r2 = sum_mod$r.squared, P = anov_mod$"Pr(>F)", lethalPoints = which(lethal), 
                AddedPoints = which(HighMO2))
  return(O2CRIT)
}