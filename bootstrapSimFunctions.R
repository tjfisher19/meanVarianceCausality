##########################################
## generateData 
##
## This function builds two dataset of length n
## 
## The underlying innovations are generated following
## four structures:
##    alternative==1 - correlated at lag 1 in mean
##    alternative==2 - correlated at lag 1 in variance
##    alternative==3 - correlation at lags 0, 1, and 2 in variance
##    alternative==0 - No correlation
##
## Each time series is generated as an ARMA(1,1)+GARCH(1,1) with parameters
##
##  (omega=0.1, alpha1=0.4, beta1=0.15, ar1=-0.5, ma1=0.3)
generateData <- function(n=100, alternative=0) {
  #########################
  ## Generate underlying error series
  N <- n + 100000
  if(alternative==1) {  ### Correlation at lag 1 in mean
    a1t <- rnorm(N+1)
    a2t <- rep(0,N)
    for(i in 1:(N)) {
      a2t[i] <- 0.3*a1t[i+1] + rnorm(1)
    }
    a1t <- a1t[1:N]
  } else if(alternative==2) { ### Correlation at lag 1 in variance
    a1t <- rnorm(N+1)
    a2t <- rep(0,N)
    for(i in 1:(N)) {
      a2t[i] <- sqrt(1+0.4*a1t[i+1]^2)*rnorm(1)
    }
    a1t <- a1t[1:N]
  } else if(alternative==3) { ### correlation at lags 0, 1, and 2 in variance
    a1t <- rnorm(N+2)
    a2t <- rep(0,N)
    for(i in 1:(N)) {
      a2t[i] <- sqrt(1+0.1*a1t[i]^2)*sqrt(1+0.2*a1t[i+1]^2)*sqrt(1+0.3*a1t[i+2]^2)*rnorm(1)
    }
    a1t <- a1t[1:N]
  }
  else { # alternative=0
    a1t <- rnorm(N)
    a2t <- rnorm(N)
  }
  a1t <- scale(a1t)[,1]
  a2t <- scale(a2t)[,1]
  
  ## Generate data using ugarchspec
  ## 
  ## These parameters are based off the GARCH model fit in Tol (1996)
  ## These are based on meterological series in De Bilt, Holland summer series
  # spec1 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
  #                     mean.model=list(armaOrder=c(1,1), include.mean=FALSE),
  #                     fixed.pars=list(omega=0.1, alpha1=0.4, beta1=0.15, ar1=-0.5, ma1=0.3) )
  # spec2 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
  #                     mean.model=list(armaOrder=c(1,1), include.mean=FALSE),
  #                     fixed.pars=list(omega=0.1, alpha1=0.4, beta1=0.15, ar1=-0.5, ma1=0.3) )
  
  spec1 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                      mean.model=list(armaOrder=c(1,1), include.mean=FALSE),
                      fixed.pars=list(omega=0.1, alpha1=0.1, beta1=0.75, ar1=-0.5, ma1=0.3) )
  spec2 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                      mean.model=list(armaOrder=c(1,1), include.mean=FALSE),
                      fixed.pars=list(omega=0.1, alpha1=0.1, beta1=0.75, ar1=-0.5, ma1=0.3) )

  path.sgarch1 <- ugarchpath(spec1, n.sim=n, n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) )
  path.sgarch2 <- ugarchpath(spec2, n.sim=n, n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) )
  
  x1 <- as.vector(path.sgarch1@path$seriesSim)
  x2 <- as.vector(path.sgarch2@path$seriesSim)
  list(x1=x1, x2=x2)
}


###########################################
## This function performs one replication of the main simulations
## 
## First, we generate two datasets of length n
## Second, fit each series as an ARMA(1,1)+GARCH(1,1)
##    ** We get errors sometimes, especially with small samples... non-stationary estimates
##       We handle this with some exception handling, only interested in stationary series for now
##
## Third, now calculate the test statistics
## Then bootstrap the sampling distribution using the fitted parameters
##
## This code does quite a bit of I/O.  This was initially needed for debugging. We decided
## to keep it here as by reading the files, outStatus.txt and outErrors.txt, you get
## an idea of how the code is running
do.one.simulation <- function(cluster=NULL, n=100, M=5, alternative=FALSE, boot.run=1000, one.way="no") {

  #############
  ## Some counting so we know we are making progress
  counter <<- counter + 1
  if(counter%%10==0) {
    time <- (proc.time()-ptm)[3]
    cat(counter, "  ", time, "\n", file=status.file, append=TRUE)
  }
  
  #### Generate the data
  cat(counter, ": generating initial data", file=status.file, append=TRUE);
  tmp <- generateData(n=n, alternative=alternative)
  x1 <- tmp$x1
  x2 <- tmp$x2
  model1 <- ~arma(1,1) + garch(1,1)
  model2 <- ~arma(1,1) + garch(1,1)
  #####################
  ## Fit data with specified models, ARMA(1,1)+GARCH(1,1)
  ##
  ## With smaller sample sizes, we run into issues of the fitted models not being
  ## stationary and/or having other numerical issues
  ##
  ## We are interested in stationary fits, so use some exception handling
  ## If there is an error with a fit. Generate new data and fit it.
  cat(" -- attempting model fits", file=status.file, append=TRUE);
  fit1.success <- FALSE
  fit2.success <- FALSE

  fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                           fit1.success = TRUE},
            error = function(e) { cat(counter, "Fitting x1", file=error.file, append=TRUE); print(e) } )

  fit2 <- tryCatch( { x2.fit=suppressWarnings(garchFit(formula=model2, data=x2, trace=FALSE, include.mean=FALSE, cont.dist="QMLE") )
                      fit2.success = TRUE},
            error = function(e) { cat(counter, "Fitting x2", file=error.file, append=TRUE); print(e) } )

  ##################
  ## If the fits are weird, generate new data!
  ## Technically, the GARCH(1,1) is stationary as long as the two fitted parameters, alpha & beta
  ## sum to less than 1. BUT we use these fitted parameters to generate more data
  ##
  ## We found that if these fitted parameters were very close to 1, the bootstrapping part below
  ## would run a very long time before finally getting stationary samples. So we cutoff a little
  ## less than 1 here. This improves run time in the code.
  error1 <- FALSE;
  error2 <- FALSE;
  if((!fit1.success) || (!fit2.success)) error1 <- TRUE    ## This happens if the AR part is not-stationary
  else if ((sum(coef(x1.fit)[4:5])>=0.99) || (sum(coef(x2.fit)[4:5])>=0.99 ) ) error2 <- TRUE  ## GARCH stationary error
  else FALSE;
  while( error1 || error2 ) { 
    cat(" -- ERROR:", error1, " ", error2, file=status.file, append=TRUE);
    cat(counter, ": Initial Fits\n", append=TRUE, file=error.file)
    if(!error1) {
      cat("    x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
      cat("    x2.fit: ", coef(x2.fit), "\n", append=TRUE, file=error.file)
    }
    ## Generate a new set of data
    tmp <- generateData(n=n, alternative=alternative)
    x1 <- tmp$x1
    x2 <- tmp$x2    
    #####################
    ## Fit data with specified models
    fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                        fit1.success = TRUE},
                      error = function(e) { cat(counter, "Fitting x1", file=error.file, append=TRUE); print(e) } )
    
    fit2 <- tryCatch( { x2.fit=suppressWarnings(garchFit(formula=model2, data=x2, trace=FALSE, include.mean=FALSE, cont.dist="QMLE") )
                        fit2.success = TRUE},
                      error = function(e) { cat(counter, "Fitting x2", file=error.file, append=TRUE); print(e) } )
    ## More error checking
    if((!fit1.success) || (!fit2.success)) error1 <- TRUE
    else if ((sum(coef(x1.fit)[4:5])>=0.99) || (sum(coef(x2.fit)[4:5])>=0.99 ) ) error2 <- TRUE
    else {
      cat("New Fits\n", append=TRUE, file=error.file)
      cat(" x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
      cat(" x2.fit: ", coef(x2.fit), "\n", append=TRUE, file=error.file)
      error1 <- FALSE
      error2 <- FALSE
    }
  }

  cat("\n", append=TRUE, file=status.file);
  
  ############################
  ## At this point, the fits should be okay.
  
  #####################
  ## get standardized residuals
  x1.res <- x1.fit@residuals/sqrt(x1.fit@h.t)
  x2.res <- x2.fit@residuals/sqrt(x2.fit@h.t)
  
  ##Calculate Test Statistics
  if(length(M)==1) {
    tmp <- get.test.stats(r1=x1.res, r2=x2.res, M=M, one.way=one.way)
    stats <- tmp[,1]
    pvals1 <- tmp[,2]
  } else {
    tmp1 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[1], one.way=one.way)
    tmp2 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[2], one.way=one.way)
    stats <- c(tmp1[,1], tmp2[,1] )
    pvals1 <- tmp1[,2]
    pvals2 <- tmp2[,2]
  }
 
  ### Now bootstrap the distribution!
  ## First I need the parameter values (I do not know them, so use estimated)
  coef1 <- coef(x1.fit)
  name1 <- substring(names(coef1),1,2)
  coef2 <- coef(x2.fit)
  name2 <- substring(names(coef2),1,2)

  fit.order.x1<-x1.fit@fit$series$order
  fit.order.x2<-x2.fit@fit$series$order
  # The code below is setup to handle up to an ARMA(3,3)+GARCH(3,3) but we only use (1,1)+(1,1)
  spec1 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(fit.order.x1[[3]],fit.order.x1[[4]])),
                      mean.model=list(armaOrder=c(fit.order.x1[[1]],fit.order.x1[[2]]), include.mean=FALSE),
                      fixed.pars=list(omega=coef1[name1=="om"][[1]], 
                                      alpha1=  ifelse(fit.order.x1[[3]],coef1[name1=="al"][[1]],0), 
                                      alpha2=  ifelse(fit.order.x1[[3]]>1,coef1[name1=="al"][[2]],0),
                                      alpha3=  ifelse(fit.order.x1[[3]]>2,coef1[name1=="al"][[3]],0),
                                      beta1= ifelse(fit.order.x1[[4]],coef1[name1=="be"][[1]],0),
                                      beta2= ifelse(fit.order.x1[[4]]>1,coef1[name1=="be"][[2]],0),
                                      beta3= ifelse(fit.order.x1[[4]]>2,coef1[name1=="be"][[3]],0),
                                      ar1= ifelse(fit.order.x1[[1]],coef1[name1=="ar"][[1]],0),
                                      ar2= ifelse(fit.order.x1[[1]]>1,coef1[name1=="ar"][[2]],0),
                                      ar3= ifelse(fit.order.x1[[1]]>2,coef1[name1=="ar"][[3]],0),
                                      ma1= ifelse(fit.order.x1[[2]],coef1[name1=="ma"][[1]],0),
                                      ma2= ifelse(fit.order.x1[[2]]>1,coef1[name1=="ma"][[2]],0),
                                      ma3= ifelse(fit.order.x1[[2]]>2,coef1[name1=="ma"][[3]],0)))
  
  spec2 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(fit.order.x2[[3]],fit.order.x2[[4]])),
                      mean.model=list(armaOrder=c(fit.order.x2[[1]],fit.order.x2[[2]]), include.mean=FALSE),
                      fixed.pars=list(omega=coef2[name2=="om"][[1]], 
                                      alpha1=  ifelse(fit.order.x2[[3]],coef2[name2=="al"][[1]],0), 
                                      alpha2=  ifelse(fit.order.x2[[3]]>1,coef2[name2=="al"][[2]],0),
                                      alpha3=  ifelse(fit.order.x2[[3]]>2,coef2[name2=="al"][[3]],0),
                                      beta1= ifelse(fit.order.x2[[4]],coef2[name2=="be"][[1]],0),
                                      beta2= ifelse(fit.order.x2[[4]]>1,coef2[name2=="be"][[2]],0),
                                      beta3= ifelse(fit.order.x2[[4]]>2,coef2[name2=="be"][[3]],0),
                                      ar1= ifelse(fit.order.x2[[1]],coef2[name2=="ar"][[1]],0),
                                      ar2= ifelse(fit.order.x2[[1]]>1,coef2[name2=="ar"][[2]],0),
                                      ar3= ifelse(fit.order.x2[[1]]>2,coef2[name2=="ar"][[3]],0),
                                      ma1= ifelse(fit.order.x2[[2]],coef2[name2=="ma"][[1]],0),
                                      ma2= ifelse(fit.order.x2[[2]]>1,coef2[name2=="ma"][[2]],0),
                                      ma3= ifelse(fit.order.x2[[2]]>2,coef2[name2=="ma"][[3]],0)))

  cat(counter, ": About to start the bootstrap distro code\n", file=status.file, append=TRUE);

  ## Get the bootstrapped sampling distribution
  sub.counter <<- 0;
  if(is.null(cluster)) {   # run in sequential model
    null.distro <- sapply(rep(n,boot.run), one.iteration.boot, spec1=spec1, spec2=spec2, 
                          M=M, model1=model1, model2=model2, one.way=one.way )
  } else {   # run in parallel
    null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.boot, spec1=spec1, spec2=spec2, 
                             M=M, model1=model1, model2=model2, one.way=one.way )
  }
  ## now calculate the p-values
  if(length(M)==1) {
    out <- c(sapply(1:10, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals1)
  } else {
    out <- c(sapply(1:10, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals1,
             sapply(11:20, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals2 )
  }
  out
}


####################################
## Here is the function that does the bootstrapping distro
##
## We are always generating data under the null distro so the
## two series are uncorrelated
## Code is largely a replication of the above but a bit simpler.
one.iteration.boot <- function(spec1, spec2, model1, model2, n, M, one.way="no") {
  N <- n+1000
  sub.counter <<- sub.counter + 1;
  ### Generate data fom the specied models
  path.sgarch1 <- try(ugarchpath(spec1, n.sim=n, n.start=(N-n), m.sim=1))
  path.sgarch2 <- try(ugarchpath(spec2, n.sim=n, n.start=(N-n), m.sim=1))
  ### occassionally we get an error here, so just retry
  while(inherits(path.sgarch1, "try-error") || inherits(path.sgarch2, "try-error")) {
    cat(counter, ":", sub.counter, ": Generation error\n", append=TRUE, file=error.file)
    cat("   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
    cat("   spec2: ", unlist(spec2@model$fixed.pars), "\n", append=TRUE, file=error.file)
    path.sgarch1 <- try(ugarchpath(spec1, n.sim=n, n.start=(N-n), m.sim=1))
    path.sgarch2 <- try(ugarchpath(spec2, n.sim=n, n.start=(N-n), m.sim=1))
  }
  ## Data is generated...
  x1 <- as.vector(path.sgarch1@path$seriesSim)
  x2 <- as.vector(path.sgarch2@path$seriesSim)
  
  
  ## Fit with the correct models. Much like in the above, for smaller n
  ## we run into computational/stationarity issues on occasion.
  x1.fit <- try(garchFit(model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
  x2.fit <- try(garchFit(model2, data=x2, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") ) 
  while(inherits(x1.fit, "try-error") || inherits(x2.fit, "try-error") )
  {
    cat(counter, ":", sub.counter,  ": Boot Fit error\n", append=TRUE, file=error.file)
    cat("   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
    cat("   spec2: ", unlist(spec2@model$fixed.pars), "\n", append=TRUE, file=error.file)
    ## Generate new data
    path.sgarch1 <- try(ugarchpath(spec1, n.sim=n, n.start=(N-n), m.sim=1))
    path.sgarch2 <- try(ugarchpath(spec2, n.sim=n, n.start=(N-n), m.sim=1))
    while(inherits(path.sgarch1, "try-error") || inherits(path.sgarch2, "try-error")) {
      cat(counter, ":", sub.counter, ": generate error after fit\n", append=TRUE, file=error.file)
      cat("   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
      cat("   spec2: ", unlist(spec2@model$fixed.pars), "\n", append=TRUE, file=error.file)
    }
    
    x1 <- as.vector(path.sgarch1@path$seriesSim)
    x2 <- as.vector(path.sgarch2@path$seriesSim)
    ## and fit the new data
    x1.fit <- try(garchFit(model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
    x2.fit <- try(garchFit(model2, data=x2, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") ) 
  }
  ## At this point, we have a good fit.
  ## save the standardized residuals
  x1.res <- x1.fit@residuals/sqrt(x1.fit@h.t)
  x2.res <- x2.fit@residuals/sqrt(x2.fit@h.t)
  
  if(length(M)==1){
    tmp <- get.test.stats(r1=x1.res, r2=x2.res, M=M, one.way=one.way)
    out <- tmp[,1]
  } else {
    tmp1 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[1], one.way=one.way)
    tmp2 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[2], one.way=one.way)
    out <- c(tmp1[,1], tmp2[,1])
  }
  out
}

