#install.packages("tidyverse")
library("tidyverse")

data <- read.csv(file.choose(), stringsAsFactors = FALSE,
                  strip.white = TRUE, na.strings = c("NA","") )

str(data)

	### Add trial.no for silverspotted sculpins
#> data$trial.no[data$spps=="blci"] <- 1
#> data$trial.no[data$spps=="blci"]
# [1] 1 1 1 1 1 1 1 1 1 1 1 1
#> data1 <- data[data$trial.no == 1,]
	### Table of pcrits available by spps and temp
#> table(data1$spps, data1$temp) ## Checking how many trials per spps per temp
#      					## NOTE: date today: 27 Sep 2017
#       12 16 20				## DOES NOT INCLUDE:
#  arfe 12  3  3
#  arha  5  3  3	## 6 at 20 - no more fish! all dead except 1 from old batch
#  arla 13  3  3
#  blci  6  3  3
#  clgl 10  3  3	## +3 mossheads at 20C
#  enbi  4  3  3	## +3 buffalo at 20C
#  hehe  5  2  3	## +3 irish lords at 20C
#  olma 10  4  3
#  scma  6  1  3	## +3 cabezon (BIG) at 20C

data1 <- data[data$trial.no == 1,]
head(data1)






