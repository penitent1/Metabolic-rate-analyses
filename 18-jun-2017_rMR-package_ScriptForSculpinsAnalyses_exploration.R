#install.packages("rMR", dependencies = TRUE)
library(rMR)

data(fishMR)
head(fishMR)
str(fishMR)

## create time variable in POSIXct format ##
fishMR$std.time <- as.POSIXct(fishMR$Date.time,
                              format = "%d/%m/%Y %I:%M:%S %p")

## My sculpin data
sculpinMO2 <- read.csv(file.choose())
sculpinMO2$std.time <- as.POSIXct(sculpinMO2$Time)#,
                             # format = "%d/%m/%Y %I:%M:%S %p")

?as.POSIXct
plot(sculpinMO2$Oxygen ~ sculpinMO2$std.time)



### Example from package help info

## load data ##
data(fishMR)

## create time variable in POSIXct format ##
fishMR$std.time <- as.POSIXct(fishMR$Date.time,
                              format = "%d/%m/%Y %I:%M:%S %p")

## calc background resp rate
bgd.resp <- 
  background.resp(fishMR, "DO.mgL", 
                  start.time = "2015-07-02 16:05:00",
                  end.time = "2015-07-02 16:35:00",
                  ylab = "DO (mg/L)", xlab = "time (min)")

bg.slope.a <- bgd.resp$mat[2]


starts <- c("2015-07-03 01:15:00", "2015-07-03 02:13:00",
            "2015-07-03 03:02:00", "2015-07-03 03:50:00",
            "2015-07-03 04:50:00")

stops <- c("2015-07-03 01:44:00", "2015-07-03 02:35:30",
           "2015-07-03 03:25:00", "2015-07-03 04:16:00",
           "2015-07-03 05:12:00")

metR <- MR.loops(data = fishMR, DO.var.name ="DO.mgL",
                 start.idx = starts, time.units = "hr",
                 stop.idx = stops, time.var.name = "std.time",
                 temp.C = "temp.C", elevation.m = 1180,
                 bar.press = NULL, in.DO.meas = "mg/L", 
                 background.consumption = bg.slope.a,
                 ylim=c(6, 8))

metR$MR.summary


## now lets assume we ran a control loop for background rate
## before and after we ran the MR loops
## let:

bg.slope.b <-bg.slope.a -0.0001 
metRa <- MR.loops(data = fishMR, DO.var.name ="DO.mgL",
                  start.idx = starts, time.units = "hr",
                  stop.idx = stops, time.var.name = "std.time",
                  temp.C = "temp.C", elevation.m = 1180,
                  bar.press = NULL, in.DO.meas = "mg/L", 
                  background.consumption = c(bg.slope.a, bg.slope.b),
                  background.indices = c("2015-07-02 16:20:00",
                                         "2015-07-03 06:00:00"),
                  ylim=c(6, 8))


metRa$MR.summary

# note that the calculated slopes
# diverge as time increases. This is
# because the background respiration
# rate is increasing.

metR$MR.summary-metRa$MR.summary

## This looks great, but you need to check your start and
## stop vectors, otherwise, you could end up with some
## atrocious loops, e.g.:

starts <- c("2015-07-03 01:15:00", "2015-07-03 02:13:00",
            "2015-07-03 03:02:00", "2015-07-03 03:50:00",
            "2015-07-03 04:50:00")

stops <- c("2015-07-03 01:50:00", "2015-07-03 02:35:30",
           "2015-07-03 03:25:00", "2015-07-03 04:16:00",
           "2015-07-03 05:12:00")

metRb <- MR.loops(data = fishMR, DO.var.name ="DO.mgL",
                  start.idx = starts,
                  stop.idx = stops, time.var.name = "std.time",
                  temp.C = "temp.C", elevation.m = 1180,
                  bar.press = NULL, in.DO.meas = "mg/L", 
                  background.consumption = bg.slope.a,
                  ylim=c(6,8))



####################################################
####################################################
####################################################

for (i in -10:10) {
  if (i %% 2){
    next
  }
  print(i)
}


# based on variable values
#newdata <- mydata[ which(mydata$gender=='F' 
#                         & mydata$age > 65), ]

head(sculpinMO2)

orig.time <- sculpinMO2$Time
std.time <- sculpinMO2$std.time

std.time[2] - std.time[1]
# orig.time[2] - orig.time[1] ## Doesn't work!!




## create fake data
tmp <- seq(as.POSIXct('2011-08-01 13:00'), as.POSIXct('2011-08-05 03:00'),
           len=42)
df <- data.frame(tm=tmp, x=seq(42))
## subset examples

## everything in hours  11am through 1pm inclusive
df3 <- subset(df, format(tm,'%H') %in% c('11','12','13'))

## 11 am through 3:59 pm on the 2nd
df4 <- subset(df, tm >= as.POSIXct('2011-08-02 11:00') & tm <=
                as.POSIXct('2011-08-02 15:59'))

#newdata <- subset(mydata, age >= 20 | age < 10, 
#                  select=c(ID, Weight))

sculpinMO2min <- subset(sculpinMO2, 
                        std.time >= as.POSIXct("2017-06-06 16:00:00") 
                        & std.time <= as.POSIXct("2017-06-07 08:00:00"))

#sculpinMO2trunc <- sculpinMO2[sculpinMO2$std.time > '2017-06-06 16:00:00'
#                              & sculpinMO2$std.time < '2017-06-06 17:00:00', ]
#plot(sculpinMO2min$Oxygen ~ sculpinMO2min$std.time)
#head(sculpinMO2trunc)
#tail(sculpinMO2trunc)

startSlope <- subset(sculpinMO2min,
                     std.time = as.POSIXct("2017-06-06 16:00:05"))


## Example for loop

j <- 1
for (i in 1:10) {
  j[i] = 10
}

j <- rep(NA, 10)
for (i in 1:10) {
  j[i] = 10
}
?rep
j
###########

head(fishMR)
str(fishMR)
fishMRtrunc <- fishMR[fishMR$std.time < "2015-07-02 22:00:00", ]
tail(fishMR)

plot(fishMR$DO.mgL ~ fishMR$std.time)





