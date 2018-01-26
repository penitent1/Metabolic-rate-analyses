## Data, copy data in excel then run script below
read.excel <- function(header=FALSE,...) {
  read.table("clipboard",sep="\t",header=header,...)
}
data <- read.excel()
x <- data$V1
y <- data$V2

## Mac version

## Data, copy data in excel then run script below
read.excel <- function(header=FALSE,...) {
  read.table(pipe("pbpaste"),sep="\t",header=header,...)
}
data <- read.excel()
x <- data$V1
y <- data$V2