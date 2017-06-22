## My sculpin data
md <- read.csv(file.choose())

### 'md' is the data frame

#wsize <- 20

md$time <- seq(0, (nrow(md)-1)*15, by=15)
md$oxygen.umol <- (((md$Oxygen-0.05)/0.999)/100)*101.85*760*1.7206*(0.46324-0.0353)/101.325 ### convert from %sat to umol O2

output <- numeric(150)
check <- numeric(150)
glms <- list()
j <- 1

i <- 1
while (i <= nrow(md)-20)
{
	time <- md$time[i:(i+20)]
	o2 <- md$oxygen[i:(i+20)]
	z <- lm(o2~time)
	
	if (summary(z)$adj.r.squared>=0.97 & coef(z)[2]<0) 
	{
		output[j] <- coef(z)[2]
		glms[[j]] <- z
		check[j] <- i
		j <- j+1
		i <- i+20
	} 
	else 
	{
		i <- i+1
	}
}
output <- output[output!=0]
### could modify output to final MO2 units here if you wanted
check <- check[check!=0]

### extend below so the last plus value matches 'wsize'
regRegions <- c(check, check+1, check+2, check+3, check+4, check+5, check+6, check+7, check+8, check+9, check+10, check+11, check+12, check+13, check+14, check+15, check+16, check+17, check+18, check+19, check+20)
regRegions <- regRegions[order(regRegions)]

plot(md$time, md$oxygen, pch=20, cex=0.67, xlab="Time (sec)", ylab="Oxygen concentratio\n (umol O2)")
points(md$time[check], md$oxygen[check], pch=1, col="darkorange", cex=1.2)
segments(md$time[check],0,md$time[check],110,col="grey67",lty=2)

plot(md$time[regRegions],md$oxygen[regRegions], pch=20, cex=0.67, ylim=c(15,110), xlab="Time", ylab="Oxygen saturation")
for (i in 1:length(glms))
{
	lines(md$time[(check[i]-10):(check[i]+30)],predict(glms[[i]], newdata=md[(check[i]-10):(check[i]+30),]),col="magenta",lty=2)
}

ind <- numeric(length(check)-1)
for (i in 1:length(ind))
{
	ind[i] <- check[i+1]-check[i]
}

ind2 <- logical(length(ind)-4)
for (i in 1:length(ind))
{
	ind2[i] <- (ind[i]==20 & ind[i+1]==20 & ind[i+2]==20 & ind[i+3]==20)
}

rmr <- output[1:(grep(TRUE, ind2)[1]-1)]
extraRmrCheck <- ind[1:(grep(TRUE, ind2)[1]-1)]
rmr <- rmr[extraRmrCheck!=20]
pcrit <- output[grep(TRUE, ind2)[1]:(length(output)-1)]
