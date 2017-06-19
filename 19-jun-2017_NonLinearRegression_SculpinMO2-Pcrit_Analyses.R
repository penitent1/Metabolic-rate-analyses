
pcrit.data <- read.table(file.choose())
head(pcrit.data)
pcrit.data$po2 <- pcrit.data$V1
pcrit.data$mo2 <- pcrit.data$V2

plot(mo2 ~ po2,
     data = pcrit.data)

pcritFctWeibel <- function(po2,a,b,c)
                      {a * (1 - exp(-(po2/b)^c))}
pcritFctWeiWithInt <- function(po2,a,b,c,d)
                      {a * (1 - exp(-(po2/b)^c)) + d}
pcritFctMM <-  function(po2,a,b)
                      {(a * po2)/(b + po2)}
pcritFctPareto <- function(po2,a,b)
                      {1 - ((a * po2) ^ b)}
pcritFctPower <- function(po2,a,b)
                      {a * po2 ^ b}

plot(mo2 ~ po2,
     data = pcrit.data,
     ylim = c(0,5),
     xlim = c(0,155))

curve(pcritFct(x, a = 1, b = 1, c = 1),
      add = TRUE)
curve(pcritFct(x, a = 2, b = 1, c = 1),
      add = TRUE)
curve(pcritFct(x, a = 4, b = 1, c = 1),
      add = TRUE)
curve(pcritFct(x, a = 3.75, b = 1, c = 1),
      add = TRUE)

curve(pcritFct(x, a = 3.75, b = 2, c = 1),
      add = TRUE)
curve(pcritFct(x, a = 3.75, b = 10, c = 1),
      add = TRUE)
curve(pcritFct(x, a = 3.75, b = 20, c = 1),
      add = TRUE)
curve(pcritFct(x, a = 3.75, b = 20, c = 1),
      add = TRUE)

curve(pcritFctMM(x, a = 3.75, b = 2),
      add = TRUE)

curve(pcritFctPareto(x, a = -0.01, b = 0.5),
      add = TRUE)

curve(pcritFctPower(x, a = 1, b = 0.3),
      add = TRUE)

curve(pcritFctWeiWithInt(x, a = 3.75, b = 1, c = 1, d = 1),
      add = TRUE)
curve(pcritFctWeiWithInt(x, a = 2, b = 1, c = 1, d = 1),
      add = TRUE)
curve(pcritFctWeiWithInt(x, a = 3.75, b = 1, c = 1, d = 0),
      add = TRUE)
curve(pcritFctWeiWithInt(x, a = 3.75, b = 1, c = 4, d = 0),
      add = TRUE)
curve(pcritFctWeiWithInt(x, a = 3.75, b = 1, c = 0.5, d = 0),
      add = TRUE)
curve(pcritFctWeiWithInt(x, a = 3.75, b = 1, c = 0.01, d = 0),
      add = TRUE)
curve(pcritFctWeiWithInt(x, a = 3.75, b = 0.11, c = 0.1, d = 0.25),
      add = TRUE)

mo2.m1 <- nls(a(1-exp(-(po2/b)^c)),
              data = pcrit.data,
              start = list(a = ,
                           b = ,
                           c = ),
              trace = TRUE)

################################################################
################################################################

  ### NLR info below

################################################################
?nls()
#install.packages("nlrwr", dependencies = TRUE)
plot(num.fish ~ spawn.biomass,
     data = M.merluccius, xlab = "Spawning biomass (1000 tonnes)",
     ylab = "Recruitment (million fish)")


conc <- c(2,5,7,22,27,39,45,203)
rate <- c(14,24,31,72,77,96,96,108)
L.minor <- data.frame(conc, rate)

L.minor.m1 <- nls(rate ~ Vm * conc/(K+conc),
                  data = L.minor, start = list(K = 20, Vm = 120),
                  trace = TRUE)

## Measures of model fit
deviance(L.minor.m1) # estimated residual standard error or variance
logLik(L.minor.m1) # LogLikelihood of estimated mean function at parameter estimates B_hat

## List parameter estimates
coef(L.minor.m1)

## Detailed summary of model fit
summary(L.minor.m1)

## Fitted values from predictor values used to estimate model
fitted(L.minor.m1)

## Example of how to obtain predicted values
concVal <- with(L.minor, # From the dataframe "L.minor"...
                seq(min(conc), max(conc), length.out = 10)) # generate a sequence of 10 values equally spaced out starting and ending with the min and max from the variable "conc"

predict(L.minor.m1, # using our estimated model
        data.frame(conc = concVal)) # using values in a data.frame containing the "concVal" values
### Note! Need to supply data.frame with "x" vector named the same as the "x" vector in the data.frame that was used to estimate the model
### I.e., the L.minor.m1 model only understands "x" values from a vector "conc"

## How to data and model
plot(rate ~ conc, data = L.minor,
     ylim = c(10,130), ylab = "Uptake rate (weight/h)", # ylim argument used to ensure sufficiently large plot window to encompass both original data and entire fitted regression curve
     xlab = Substrate ~ concentration ~ (mmol ~ m^-3)) # NOTE didn't use "" marks because want to be able to use mathematical notation!

concVal <- with(L.minor, seq(min(conc),
                             max(conc),
                             length.out = 100)) # Increase number of predicted values
# used to plot model line - the line is plotted by linear interpolation
# i.e. connecting predicted values with line segments

lines(concVal, predict(L.minor.m1,
                       newdata = data.frame(conc=concVal))) # the function that connects the "predicted value dots"

abline(h = coef(L.minor.m1)[2],
       lty = 2) # Add a dashed line to show the horizontal asymptote
