################################
##
## Basically a modification to modelFitsPlots.R
##
## Here we know which models we are using and we bootstrap the p-values
## for the various test.  This code is run from the command line.
##

library(fGarch)
library(WeightedPortTest)
library(forecast)
library(rugarch)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)
library(portes)
library(xtable)
library(snow)

source("../sourceFiles/nonMatrixTest.R")
source("../sourceFiles/matrixBasedTests.R")
source("../sourceFiles/stableNullSimFunctions.R")
source("../sourceFiles/testStatFunctions.R")
source("../sourceFiles/stableVarSimFunctions.R")
source("../sourceFiles/newVarimaSim.R")



num.nodes <- 11;
boot.run <- 1000
error.file <- paste("outErrors.txt")
status.file <- paste("outStatus.txt")
cat("Run Status\n*************************\n", file=status.file)
cat("Errors\n***************************\n", file=error.file)

library(snow)
options(width=120);

cluster <- makeSOCKcluster(spec=num.nodes, names=rep("localhost", num.nodes) );

checkCluster(cluster);
set.seed(41682)
clusterSetupRNG(cluster, seed=rep(41682, 6) );
invisible(clusterEvalQ(cluster, library(fGarch)));
invisible(clusterEvalQ(cluster, library(rugarch)));
invisible(clusterEvalQ(cluster, library(portes)));
invisible(clusterEvalQ(cluster, library(forecast)));
counter <- 0
sub.counter <- 0
clusterExport(cluster, ls() )


data.sp500 <- read.csv("YAHOO-INDEX_GSPC.csv")
data.sp500<-data.sp500[,c(1,7)]
data.sp500$Date<- ymd(data.sp500$Date)
colnames(data.sp500)<-c("Date","sp500")

data.EUR.2015 <- read.csv("EUR2015.csv")
data.EUR.2015$Date<- ymd(data.EUR.2015$Date)

data.full<- left_join(data.EUR.2015,data.sp500)

colnames(data.full)<-c("Date","Exchange Rate (EUR/USD)","S&P500 Index")
data.full$Closed <- as.factor(is.na(data.full$'S&P500 Index')*1)


### REMOVE MISSING VALUES
colnames(data.full)<-c("Date","EUR","sp500","Closed")
data.full <- na.omit(data.full)

n <- dim(data.full)[1]

## LOG DIFFERENCING
sp500<- diff(log(data.full$sp500))
EUR<- diff(log(data.full$EUR))

## Explore features about sp500
arma.fit.sp500 <- Arima(sp500, order=c(1,0,2), include.mean=FALSE)
Weighted.Box.test(arma.fit.sp500$residuals, lag=30, fitdf=3)

### GARCH(1,1)???
Weighted.Box.test(arma.fit.sp500$residuals, lag=30, fitdf=3, sqrd.res = TRUE)

sp500.fit <- garchFit(~arma(1,2)+garch(1,1), data=sp500, trace=FALSE, include.mean=FALSE)
sp500.res <- sp500.fit@residuals/sqrt(sp500.fit@h.t) 
#####
## Marginall fit with an ARMA(1,2)+GARCH(1,1)



########################
# EUR

arma.fit.EUR <- Arima(EUR, order=c(0,0,1), include.mean=FALSE)
Weighted.Box.test(arma.fit.EUR$residuals, lag=30, fitdf=1)

Weighted.Box.test(arma.fit.EUR$residuals, lag=30, fitdf=1, sqrd.res = TRUE)

#### No GARCH here...
EUR.res <- arma.fit.EUR$residuals
EUR.fit <- arma.fit.EUR

EUR.res <- scale(EUR.res)[,1]

########## Just an MA(1)...



## set up the model
model1 <- ~arma(1,2)+garch(1,1)
model2 <- c(0,0,1)

#### Fit the data, get the parameters and set up the simulation
x1.fit <- sp500.fit
x2.fit <- EUR.fit
coef1 <- coef(x1.fit)
name1 <- substring(names(coef1),1,2)
coef2 <- coef(x2.fit)
name2 <- substring(names(coef2),1,2)

fit.order.x1<-x1.fit@fit$series$order
fit.order.x2 <- x2.fit$arma[1:2]
spec1 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                    mean.model=list(armaOrder=c(1,2), include.mean=FALSE),
                    fixed.pars=list(omega=coef1[name1=="om"][[1]], 
                                    alpha1=  ifelse(fit.order.x1[[3]],coef1[name1=="al"][[1]],0), 
                                    beta1= ifelse(fit.order.x1[[4]],coef1[name1=="be"][[1]],0),
                                    ar1= ifelse(fit.order.x1[[1]],coef1[name1=="ar"][[1]],0),
                                    ma1= ifelse(fit.order.x1[[2]],coef1[name1=="ma"][[1]],0),
                                    ma2= ifelse(fit.order.x1[[2]]>1,coef1[name1=="ma"][[2]],0)) )
spec2 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(0,0)),
                    mean.model=list(armaOrder=c(0,1), include.mean=FALSE),
                    fixed.pars=list(omega=x2.fit$sigma2, 
                                    ma1= ifelse(fit.order.x2[[2]],coef2[name2=="ma"][[1]],0)) )

full.stat.m6 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=6)
full.stat.m9 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=9)
full.stat.m16 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=16)
full.stats <- c(full.stat.m6[,1], full.stat.m9[,1], full.stat.m16[,1] )

right.stat.m6 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=6, one.way="right")
right.stat.m9 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=9, one.way="right")
right.stat.m16 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=16, one.way="right")
right.stats <- c(right.stat.m6[,1], right.stat.m9[,1], right.stat.m16[,1] )

left.stat.m6 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=6, one.way="left")
left.stat.m9 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=9, one.way="left")
left.stat.m16 <- get.test.stats(r1=sp500.res, r2=EUR.res, M=16, one.way="left")
left.stats <- c(left.stat.m6[,1], left.stat.m9[,1], left.stat.m16[,1] )


M <- c(6,9,16)

if(is.null(cluster)) {   # run in sequence
  full.null.distro <- sapply(rep(n,boot.run), one.iteration.stable.null.boot, spec1=spec1, spec2=spec2, 
                        M=M, model1=model1, model2=model2, one.way="no" )
} else {   # run in parallel
  full.null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.stable.null.boot, spec1=spec1, spec2=spec2, 
                           M=M, model1=model1, model2=model2, one.way="no" )
}

if(is.null(cluster)) {   # run in sequence
  right.null.distro <- sapply(rep(n,boot.run), one.iteration.stable.null.boot, spec1=spec1, spec2=spec2, 
                             M=M, model1=model1, model2=model2, one.way="right" )
} else {   # run in parallel
  right.null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.stable.null.boot, spec1=spec1, spec2=spec2, 
                                M=M, model1=model1, model2=model2, one.way="right" )
}
 
if(is.null(cluster)) {   # run in sequence
  left.null.distro <- sapply(rep(n,boot.run), one.iteration.stable.null.boot, spec1=spec1, spec2=spec2, 
                              M=M, model1=model1, model2=model2, one.way="left" )
} else {   # run in parallel
  left.null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.stable.null.boot, spec1=spec1, spec2=spec2, 
                                 M=M, model1=model1, model2=model2, one.way="left" )
}


out.full <- rbind(sapply(1:10, get.bootstrap.pval, null.distro=full.null.distro, test.stats=full.stats),
              sapply(11:20, get.bootstrap.pval, null.distro=full.null.distro, test.stats=full.stats),
              sapply(21:30, get.bootstrap.pval, null.distro=full.null.distro, test.stats=full.stats))
out.full <- t(out.full)

right.full <- rbind(sapply(1:10, get.bootstrap.pval, null.distro=right.null.distro, test.stats=right.stats),
                  sapply(11:20, get.bootstrap.pval, null.distro=right.null.distro, test.stats=right.stats),
                  sapply(21:30, get.bootstrap.pval, null.distro=right.null.distro, test.stats=right.stats))
right.full <- t(right.full)

left.full <- rbind(sapply(1:10, get.bootstrap.pval, null.distro=left.null.distro, test.stats=left.stats),
                  sapply(11:20, get.bootstrap.pval, null.distro=left.null.distro, test.stats=left.stats),
                  sapply(21:30, get.bootstrap.pval, null.distro=left.null.distro, test.stats=left.stats))
left.full <- t(left.full)

full.pval <- cbind(out.full, left.full, right.full)
rownames(full.pval) <- c("Boot BP Mean", "Boot BP Var",
                         "Boot LB Mean", "Boot LB Var",
                         "Boot WLB Mean", "Boot WLB Var",
                         "Boot Dan Mean", "Boot Dan Var",
                         "Boot Mat Mean", "Boot Mat Var" )
colnames(full.pval) <- c("m=6", "m=9", "m=16","m=6", "m=9", "m=16","m=6", "m=9", "m=16")
full.pval
write.csv(full.pval, file="causalityInMeanBoot.csv")



X <- cbind(sp500, EUR)
var.fit <- ar(X, aic=FALSE, order.max=4)
var.fit
resids <- var.fit$resid[-(1:4),]

var.sp500.garch.fit <- garchFit(~garch(1,1), data=resids[,1], include.mean=FALSE, trace=FALSE, cond.dist="QMLE")
coef1 <- coef(var.sp500.garch.fit)
name1 <- substring(names(coef1),1,2)
var.sp500.garch.res <- var.sp500.garch.fit@residuals/sqrt(var.sp500.garch.fit@h.t)
var(resids[,2])
var.eur.res <- scale(resids[,2])[,1]



full.stat.m6 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=6)
full.stat.m9 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=9)
full.stat.m16 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=16)
full.stats <- c(full.stat.m6[,1], full.stat.m9[,1], full.stat.m16[,1] )

right.stat.m6 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=6, one.way="right")
right.stat.m9 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=9, one.way="right")
right.stat.m16 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=16, one.way="right")
right.stats <- c(right.stat.m6[,1], right.stat.m9[,1], right.stat.m16[,1] )

left.stat.m6 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=6, one.way="left")
left.stat.m9 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=9, one.way="left")
left.stat.m16 <- get.test.stats(r1=var.sp500.garch.res, r2=var.eur.res, M=16, one.way="left")
left.stats <- c(left.stat.m6[,1], left.stat.m9[,1], left.stat.m16[,1] )


spec1 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1) ),
                    mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                    fixed.pars=list(omega=coef1[name1=="om"][[1]], 
                                    alpha1=  ifelse(fit.order.x1[[3]],coef1[name1=="al"][[1]],0), 
                                    alpha2=  ifelse(fit.order.x1[[3]]>1,coef1[name1=="al"][[2]],0),
                                    alpha3=  ifelse(fit.order.x1[[3]]>2,coef1[name1=="al"][[3]],0),
                                    beta1= ifelse(fit.order.x1[[4]],coef1[name1=="be"][[1]],0),
                                    beta2= ifelse(fit.order.x1[[4]]>1,coef1[name1=="be"][[2]],0),
                                    beta3= ifelse(fit.order.x1[[4]]>2,coef1[name1=="be"][[3]],0)) )
spec2 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(0,0)),
                    mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                    fixed.pars=list(omega=var(resids[,2]) ) )  ## A GARCH(0,0)


M <- c(6,9,16)

if(is.null(cluster)) {   # run in sequence
  full.null.distro <- sapply(rep(n,boot.run), one.iteration.stable.var.boot, var.fit=var.fit, spec1=spec1, spec2=spec2, 
                             M=M, model1=model1, model2=model2, one.way="no" )
} else {   # run in parallel
  full.null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.stable.var.boot, var.fit=var.fit, spec1=spec1, spec2=spec2, 
                                M=M, model1=model1, model2=model2, one.way="no" )
}

if(is.null(cluster)) {   # run in sequence
  right.null.distro <- sapply(rep(n,boot.run), one.iteration.stable.var.boot, var.fit=var.fit, spec1=spec1, spec2=spec2, 
                              M=M, model1=model1, model2=model2, one.way="right" )
} else {   # run in parallel
  right.null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.stable.var.boot, var.fit=var.fit, spec1=spec1, spec2=spec2, 
                                 M=M, model1=model1, model2=model2, one.way="right" )
}

if(is.null(cluster)) {   # run in sequence
  left.null.distro <- sapply(rep(n,boot.run), one.iteration.stable.var.boot, var.fit=var.fit, spec1=spec1, spec2=spec2, 
                             M=M, model1=model1, model2=model2, one.way="left" )
} else {   # run in parallel
  left.null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.stable.var.boot, var.fit=var.fit, spec1=spec1, spec2=spec2, 
                                M=M, model1=model1, model2=model2, one.way="left" )
}

out.full <- rbind(sapply(1:10, get.bootstrap.pval, null.distro=full.null.distro, test.stats=full.stats),
                  sapply(11:20, get.bootstrap.pval, null.distro=full.null.distro, test.stats=full.stats),
                  sapply(21:30, get.bootstrap.pval, null.distro=full.null.distro, test.stats=full.stats))
out.full <- t(out.full)

right.full <- rbind(sapply(1:10, get.bootstrap.pval, null.distro=right.null.distro, test.stats=right.stats),
                    sapply(11:20, get.bootstrap.pval, null.distro=right.null.distro, test.stats=right.stats),
                    sapply(21:30, get.bootstrap.pval, null.distro=right.null.distro, test.stats=right.stats))
right.full <- t(right.full)

left.full <- rbind(sapply(1:10, get.bootstrap.pval, null.distro=left.null.distro, test.stats=left.stats),
                   sapply(11:20, get.bootstrap.pval, null.distro=left.null.distro, test.stats=left.stats),
                   sapply(21:30, get.bootstrap.pval, null.distro=left.null.distro, test.stats=left.stats))
left.full <- t(left.full)

full.pval <- cbind(out.full, left.full, right.full)
rownames(full.pval) <- c("Boot BP Mean", "Boot BP Var",
                         "Boot LB Mean", "Boot LB Var",
                         "Boot WLB Mean", "Boot WLB Var",
                         "Boot Dan Mean", "Boot Dan Var",
                         "Boot Mat Mean", "Boot Mat Var" )
colnames(full.pval) <- c("m=6", "m=9", "m=16","m=6", "m=9", "m=16","m=6", "m=9", "m=16")
full.pval
write.csv(full.pval, file="causalityInVarBoot.csv")


