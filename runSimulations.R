#####################################################
##
## Here is the main code run.
##
## 
##  $> R CMD BATCH runSimulations.R
##
##
##  Below you will see RNG seeds. There is no rhyme or reason behind them (family birthdays)
##  but we wanted to be able to reproduce results
##
##  You will also note some run times for specific tables.
##   These are specific to the machine used for the paper.
##
##   A Dell Precision T3600
##    hex-core Intel Xeon CPU E5-1650
##    24 GB of RAM
##    R version 3.2.2
##    Linux Kernel 3.13.0-66
##
##  18000 seconds is 5-hours...
##
## except for the 'trivial' example, nearly everything takes ~5 hours to run.


source("sourceFiles/testStatFunctions.R")
source("sourceFiles/bootstrapSimFunctions.R")
source("sourceFiles/simulationCodes.R")
source("sourceFiles/nonMatrixTest.R")
source("sourceFiles/matrixBasedTests.R")
source("sourceFiles/trivialSimFunctions.R")
source("sourceFiles/varSimFunctions.R")
source("sourceFiles/newVarimaSim.R")
source("sourceFiles/stableSimulationCodes.R")
source("sourceFiles/stableNullSimFunctions.R")
source("sourceFiles/stableVarSimFunctions.R")
source("sourceFiles/sizeDetermineSimFunctions.R")

num.nodes <- 11;
run.times <- 10000
boot.run <- 1000
error.file <- paste("outErrors.txt")
status.file <- paste("outStatus.txt")

library(snow)
options(width=120);

cl <- makeSOCKcluster(spec=num.nodes, names=rep("localhost", num.nodes) );

checkCluster(cl);
clusterSetupRNG(cl, seed=rep(22681, 6) );

library(fGarch)
# library(rugarch)
library(portes)
library(forecast)

invisible(clusterEvalQ(cl, library(fGarch)));
# invisible(clusterEvalQ(cl, library(rugarch)));
invisible(clusterEvalQ(cl, library(portes)));
invisible(clusterEvalQ(cl, library(forecast)));
clusterExport(cl, ls() );


clusterSetupRNG(cl, seed=rep(41682, 6) );
cat("Run Status\n*************************\n", file=status.file)
cat("Error Status\n************************\n", file=error.file)
counter <- 0
sub.counter <- 0
clusterExport(cl, ls() )


#######################################
## iid Normal data...
#######################################
# counter <- 0
# ptm <- proc.time()
# set.seed(41682)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(41682, 6) );  ## Set cluster RNG
# out <- simEmpiricalTrivialSize(cluster=cl, n=50, M=c(4,10), run.times=run.times, boot.runs=boot.run)
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/n50M4and10empiricalSizeNormalData.csv")

### > proc.time()
### user   system  elapsed 
### 27.574    0.773 2562.031 


############################################################
##
##  Type I error rates for paper...
############################################################
# counter <- 0
# ptm <- proc.time()
# set.seed(41682)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(41682, 6) );  ## Set cluster RNG
# out <- simEmpiricalSize(cluster=cl, n=100, M=c(5,8), run.times=run.times, boot.runs=boot.run, alternative=0)
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/n100M5and8empiricalSize.csv")
# print(proc.time()-ptm)
### > print(proc.time()-ptm)
### user    system   elapsed 
### 275.406     3.824 19027.922


# counter <- 0
# ptm <- proc.time()
# set.seed(22681)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(22681, 6) );  ## Set cluster RNG
# out <- simEmpiricalSize(cluster=cl, n=50, M=c(4,7), run.times=run.times, boot.runs=boot.run, alternative=0)
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/n50M4and7empiricalSize.csv")
# print(proc.time()-ptm)
### > print(proc.time()-ptm)
### user    system   elapsed 
### 371.518     4.129 17903.346 

# counter <- 0
# ptm <- proc.time()
# set.seed(41682)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(41682, 6) );  ## Set cluster RNG
# out <- simEmpiricalSize(cluster=cl, n=250, M=c(9,16), run.times=run.times, boot.runs=boot.run, alternative=0, one.way = "right")
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/n250M9and16empiricalSizeOneWayRight.csv")
# print(proc.time()-ptm)
### > print(proc.time()-ptm)
### user    system   elapsed 
### 272.641     4.776 24037.350


#################################
## Causality in mean...
#################################

# counter <- 0
# ptm <- proc.time()
# set.seed(41682)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(41682, 6) );  ## Set cluster RNG
# out <- simEmpiricalSize(cluster=cl, n=100, M=c(5,8), run.times=run.times, boot.runs=boot.run, alternative=1)
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/n100M5and8empiricalPowerInMeanLag1.csv")
# print(proc.time()-ptm)

### print(proc.time()-ptm)
### user    system   elapsed 
### 671.660     4.492 19460.227 

# counter <- 0
# ptm <- proc.time()
# set.seed(41682)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(41682, 6) );  ## Set cluster RNG
# out <- simEmpiricalVARSize(cluster=cl, n=100, M=c(5,8), run.times=run.times, boot.runs=boot.run)
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/varN100M5and8empiricalSize.csv")
# print(proc.time()-ptm)

###> print(proc.time()-ptm)
###  user    system   elapsed 
###  114.131     3.319 12869.092   
### runs quicker because less parameters fit within garchFit() function

##########################################
## Causality in Variance
##########################################
# counter <- 0
# ptm <- proc.time()
# set.seed(22681)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(22681, 6) );  ## Set cluster RNG
# out <- simEmpiricalSize(cluster=cl, n=50, M=c(4,7), run.times=run.times, boot.runs=boot.run, alternative=3)
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/n50M4and7empiricalPowerAlt3.csv")
# print(proc.time()-ptm)
### > print(proc.time()-ptm)
### user        system   elapsed 
### 1194.409     5.625 18674.051  

# counter <- 0
# ptm <- proc.time()
# set.seed(41682)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(41682, 6) );  ## Set cluster RNG
# out <- simEmpiricalSize(cluster=cl, n=250, M=c(9,16), run.times=run.times, boot.runs=boot.run, alternative=2, one.way = "right")
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/n250M9and16empiricalPowerAlt2OneWayRight.csv")
# print(proc.time()-ptm)
### > print(proc.time()-ptm)
### user    system   elapsed 
### 685.407     5.661 24302.626 

###############################################
## Robustness Study
###############################################
# 
# counter <- 0
# ptm <- proc.time()
# set.seed(22681)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(22681, 6) );  ## Set cluster RNG
# out <- simStableEmpiricalSize(cluster=cl, n=251, M=c(6, 9, 16), run.times=run.times, boot.runs=boot.run, one.way="no")
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/stableN251M6and9and16null.csv")
# print(proc.time()-ptm)
###> proc.time()
###     user    system   elapsed 
###  162.672     2.988 24433.895 

# counter <- 0
# ptm <- proc.time()
# set.seed(22681)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(22681, 6) );  ## Set cluster RNG
# out <- simStableVarEmpiricalSize(cluster=cl, n=251, M=c(6, 9, 16), run.times=run.times, boot.runs=boot.run, one.way="no")
# out <- 100*out
# out
# write.csv(out, file="resultCsvFiles/stableN251M6and9and16VarNull.csv")
# print(proc.time()-ptm)
### > print(proc.time()-ptm)
### user    system   elapsed 
### 140.065     2.736 20991.414 


#######################################
## size determination for iid Normal data
#######################################
# counter <- 0
# ptm <- proc.time()
# set.seed(41682)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(41682, 6) );  ## Set cluster RNG
# out <- sampleSizeDetermination(cluster=cl, n=seq(25,1500,1), M=5, run.times=run.times)
# write.csv(out[[2]], file="resultCsvFiles/sizeDeterminationAlpha05Mfixed5.csv")
# print(proc.time()-ptm)
### > print(proc.time()-ptm)
### user    system   elapsed 
### 41.228     1.556 13308.617 

# ptm <- proc.time()
# set.seed(41682)   ## Set local RNG
# clusterSetupRNG(cl, seed=rep(41682, 6) );  ## Set cluster RNG
# out <- sampleSizeDetermination(cluster=cl, n=seq(25,1500,1), run.times=run.times)
# write.csv(out[[2]], file="resultCsvFiles/sizeDeterminationAlpha05MgrowsLogN.csv")
# print(proc.time()-ptm)
### > print(proc.time()-ptm)
### user    system   elapsed 
### 40.532     1.376 11091.981 

stopCluster(cl);


