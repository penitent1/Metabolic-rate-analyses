## My sculpin data
md <- read.csv(file.choose())
respVols <- read.csv(file.choose())

### 'md' is the data frame

### Probe drift correction data

pre <- c(0,100)
post <- c(0.7,101) # Update for each trial
probeDrift <- lm(post~pre)
probeDrift ## "pre" is the slope of the line

wsize <- 20

respVol <- respVols$resp.volume_mL[respVols$respirometer.id == "d"] # Update each trial
mass <- 8.4 # Update each trial
alpha <- 1.7206 ## Change as appropriate for different temperatures

md$time <- seq(0, (nrow(md)-1)*15, by=15)
# Update "oxygen.umol" for each trial
md$oxygen.umol <- (((md$Oxygen-0.7)/1.003)/100)*102.9*760*alpha*(respVol/1000-mass/1000)/101.325 ### convert from %sat to umol O2

output <- numeric(150)
check <- numeric(150)
meano2 <- numeric(150)
time_smr_plot <- numeric(150)
glms <- list()
j <- 1

i <- 1
while (i <= nrow(md)-wsize)
{
	time <- md$time[i:(i+wsize)]
	o2 <- md$oxygen.umol[i:(i+wsize)]
	z <- lm(o2~time)
	
	if (summary(z)$adj.r.squared>=0.97 & coef(z)[2]<0) #  & coef(z)[2]>-0.1 
	{
		output[j] <- coef(z)[2]
		glms[[j]] <- z
		check[j] <- i
		meano2[j] <- mean(o2)
		time_smr_plot[j] <- mean(time)
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
output <- output*-1*3600 #/mass
output_ms <- output/mass
meano2 <- meano2[meano2!=0] ## NOTE These are umol units - convert to PO2
meano2 <- meano2/alpha/(respVol/1000-respVol/1000) # This converts back to torr
check <- check[check!=0]
time_smr_plot <- time_smr_plot[time_smr_plot!=0]

### extend below so the last plus value matches 'wsize'
regRegions <- c(check, check+1, check+2, check+3, check+4, check+5, check+6, check+7, check+8, check+9, check+10, check+11, check+12, check+13, check+14, check+15, check+16, check+17, check+18, check+19, check+20)
regRegions <- regRegions[order(regRegions)]

plot(md$time, md$oxygen, pch=20, cex=0.67, xlab="Time (sec)", ylab="Oxygen concentratio\n (umol O2)")
points(md$time[check], md$oxygen[check], pch=1, col="darkorange", cex=1.2)
segments(md$time[check],0,md$time[check],500,col="grey67",lty=2)

plot(md$time[regRegions],md$oxygen[regRegions], pch=20, cex=0.67, ylim=c(200,350), xlab="Time", ylab="[Dioxygen] (umol)")
for (i in 1:length(glms))
{
	lines(md$time[(check[i]):(check[i]+30)], #-10 , 
	      predict(glms[[i]],newdata=md[(check[i]):(check[i]+30),]), #-10
	      col="magenta",lty=2)
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
pcrit <- output[grep(TRUE, ind2)[1]:(length(output)-1)]
#rmr <- rmr*-1*3600/21.2
#pcrit <- pcrit*-1*3600/21.2

#write.excel <- function(x,row.names=FALSE,col.names=TRUE,...) {
#  write.table(x,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
#}

#write.excel(rmr)
#write.excel(pcrit)

