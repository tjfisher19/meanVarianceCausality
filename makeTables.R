#########################################
##
## These functions allow us to make the tables for the paper easily.
##
## It outputs the rejection rates in the order which matches the paper
##
makeTable <- function(file="n50M4and10empiricalSizeNormalData.csv") {
  full.file <- paste("resultCsvFiles/", file, sep="")
  out <- read.csv(full.file)
  row.order <- c(1,3,7,5,9,2,4,8,6,10,
                 11,13,17,15,19,12,14,18,16,20)
  col.order <- c(1,2,3,4,6,7,8)
  out[row.order,col.order]
}

makeTableRobust <- function(file="stableN251M6and9and16null.csv") {
  full.file <- paste("resultCsvFiles/", file, sep="")
  out <- read.csv(full.file)
  row.order <- c(1,3,7,5,9,2,4,8,6,10,
                 11,13,17,15,19,12,14,18,16,20)
  out[row.order,]
}

makeTablePvals <- function(file="causalityInMeanBoot.csv", mean=TRUE) {
  full.file <- paste("empiricalExample/", file, sep="")
  out <- read.csv(full.file)
  if(mean)
    row.order <- c(1,3,7,5,9)
  else row.order <- c(2,4,8,6,10)
  out[row.order,]
}