####################
## The main call function that performs the VAR robustness studies in the paper
##
## A modification of stableNullSimFunctions.R
##
## Here we modify the stableNullSimsFunction.R code to include a VAR(p) fit.
##
generateStableVarData <- function(n=100) {
  N <- n + 1000 + 100
  ### SP500 Residuals...
  a1t <- rsged(N, mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)
  ### EUR/USD residuals
  a2t <- rsged(N, mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498)
  
  a1t <- scale(a1t)[,1]
  a2t <- scale(a2t)[,1]
  
  ## Generate data using ugarchspec
  ## 
  ## These parameters are based off the GARCH model fit in Tol (1996)
  ## These are based on meterological series in De Bilt, Holland summer series
  spec1 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                      mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                      fixed.pars=list(omega=1.078e-05, alpha1=1.537e-01, beta1=7.174e-01) )
  spec2 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(0,0)),
                      mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                      fixed.pars=list(omega=2.846e-05) )

  ### Generate the two "error" series
  path.sgarch1 <- ugarchpath(spec1, n.sim=(n+100), n.start=(N-n-100), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) )
  path.sgarch2 <- ugarchpath(spec2, n.sim=(n+100), n.start=(N-n-100), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) )
  
  #### From the VAR fit...

  x1 <- as.vector(path.sgarch1@path$seriesSim)
  x2 <- as.vector(path.sgarch2@path$seriesSim)
  Innov <- cbind(x1, x2)  ## Combine
  
  ### Phi is the from fitted VAR(4) model
  Phi <- array(c(0.045,-0.116,0.0211,0.358, -0.133,0.0413,0.054,-0.072, -0.0446,0.0702,-0.1126,-0.039, -0.1720,-0.0101,0.2296,0.1093), dim=c(2,2,4))
  ### Sigma is needed but not used... The new.varima.sim() function allows you to 
  ### specify the innovations, so the Signma paramter is never used.
  ### These are the values from the VAR though
  Sigma <- matrix(c(9.2886e-05, -6.61129e-05, -6.61129e-05,2.9549e-05), nrow=2,ncol=2)
  X <- varima.sim(model=list(ar=Phi), innov=Innov, k=ncol(Innov), n=(n+100), sigma=Sigma )
  x1 <- X[-(1:100),1]
  x2 <- X[-(1:100),2]
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
do.one.var.stable.var.simulation <- function(cluster=NULL, n=100, M=5, boot.run=1000, one.way="no") {

  #############
  ## Some counting so we know we are making progress
  counter <<- counter + 1
  if(counter%%5==0) {
    time <- (proc.time()-ptm)[3]
    cat(counter, "  ", time, "\n", file=status.file, append=TRUE)
  }
  
  #### Generate the data
  cat(counter, ": generating initial data", file=status.file, append=TRUE);
  tmp <- generateStableVarData(n=n)
  x1 <- tmp$x1
  x2 <- tmp$x2
  
  cat(" -- done\n", file=status.file, append=TRUE)
    
  cat(counter, " -- attempting model fits", file=status.file, append=TRUE);
  
  X <- cbind(x1, x2)
  fit <- ar(X, aic=FALSE, order.max=4)  ## Fit a VAR(4)
  p <- fit$order
  resid <- fit$resid
  x1 <- resid[-(1:p),1]   # Remove the NAs...
  x2 <- resid[-(1:p),2]
  
  model1 <- ~garch(1,1)
  model2 <- ~garch(1,1)
  #####################
  ## Fit data with specified models, ARMA(1,1)+GARCH(1,1)
  ##
  ## With smaller sample sizes, we run into issues of the fitted models not being
  ## stationary and/or having other numerical issues
  ##
  ## We are interested in stationary fits, so use some exception handling
  ## If there is an error with a fit. Generate new data and fit it.

  fit1.success <- FALSE
  fit2.success <- TRUE    ## Nothing to fit, so a success!

  fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                           fit1.success = TRUE},
            error = function(e) { cat(counter, "Fitting x1", file=error.file, append=TRUE); print(e) } )

  ##################
  ## If the fits are weird, generate new data!
  error1 <- FALSE;
  error2 <- FALSE;
  if((!fit1.success) || (!fit2.success)) error1 <- TRUE    ## This happens if the AR part is not-stationary
  else if ((sum(coef(x1.fit)[2:3])>=0.99)  ) error2 <- TRUE  ## GARCH stationary error
  else FALSE;
  while( error1 || error2 ) { 
    cat(" -- ERROR:", error1, " ", error2, file=status.file, append=TRUE);
    cat(counter, ": Initial Fits\n", append=TRUE, file=error.file)
    if(!error1) {
      cat("    x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
    }
    ## Generate a new set of data
    tmp <- generateStableVarData(n=n)
    x1 <- tmp$x1
    x2 <- tmp$x2
    
    X <- cbind(x1, x2)
    fit <- ar(X, aic=FALSE, order.max=4)
    p <- fit$order
    resid <- fit$resid
    x1 <- resid[-(1:p),1]   # Remove the NAs...
    x2 <- resid[-(1:p),2]
 
    #####################
    ## Fit data with specified models
    fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                        fit1.success = TRUE},
                      error = function(e) { cat(counter, "Fitting x1", file=error.file, append=TRUE); print(e) } )
    error1 <- FALSE;
    error2 <- FALSE;
    if((!fit1.success) || (!fit2.success)) error1 <- TRUE
    else if ((sum(coef(x1.fit)[2:3])>=0.99)  ) error2 <- TRUE
    else {
      cat("New Fits\n", append=TRUE, file=error.file)
      cat(" x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
      error1 <- FALSE
      error2 <- FALSE
    }
  }

  cat(" -- done\n", append=TRUE, file=status.file);
  
  ############################
  ## At this point, the fits should be okay.
  
  #####################
  ## get standardized residuals
  x1.res <- x1.fit@residuals/sqrt(x1.fit@h.t)
  #x2.res <- x2.fit@residuals/sqrt(x2.fit@h.t)
  x2.res <- scale(x2)[,1]  ### Standardize...

  ##Calculate Test Statistics
  if(length(M)==1) {
    tmp <- get.test.stats(r1=x1.res, r2=x2.res, M=M, one.way=one.way)
    stats <- tmp[,1]
    pvals1 <- tmp[,2]
  } else {
    tmp1 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[1], one.way=one.way)
    tmp2 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[2], one.way=one.way)
    tmp3 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[3], one.way=one.way)
    stats <- c(tmp1[,1], tmp2[,1], tmp3[,1] )
    pvals1 <- tmp1[,2]
    pvals2 <- tmp2[,2]
    pvals3 <- tmp3[,2]
  }
  ### Now bootstrap the distribution!
  ## First I need the parameter values (I do not know them, so use estimated)
  coef1 <- coef(x1.fit)
  name1 <- substring(names(coef1),1,2)
  
  fit.order.x1<-x1.fit@fit$series$order
  
  # The code below is setup to handle an ARMA(3,3)+GARCH(3,3) but we only use (1,1)+(1,1)
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
  spec2 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(0,0)),
                      mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
                      fixed.pars=list(omega=fit$var.pred[2,2] ) )  ## A GARCH(0,0)

  cat(counter, ": About to start the bootstrap distro code", file=status.file, append=TRUE);

  sub.counter <<- 0;
  if(is.null(cluster)) {   # run in sequence
    null.distro <- sapply(rep(n,boot.run), one.iteration.stable.var.boot, var.fit=fit, spec1=spec1, spec2=spec2,
                          M=M, model1=model1, model2=model2, one.way=one.way )
    out <- c(sapply(1:10, get.boostrap.pval, null.distro=null.distro, test.stats=stats),
             pvals1) 
  } else {   # run in parallel
    null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.stable.var.boot, var.fit=fit, spec1=spec1, spec2=spec2,
                             M=M, model1=model1, model2=model2, one.way=one.way )

    out <- c(sapply(1:10, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals1,
             sapply(11:20, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals2,
             sapply(21:30, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals3 )
  }
  cat(" -- done\n", file=status.file, append=TRUE)

  out
}


####################################
## Here is the function that does the bootstrapping distro
##
## We are always generating data under the null distro so the
## two series are uncorrelated
one.iteration.stable.var.boot <- function(var.fit, spec1, spec2, model1, model2, n, M, one.way="no") {
  N <- n+1000+100
  sub.counter <<- sub.counter + 1;
  a1t <- rsged(N, mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)
  a2t <- rsged(N, mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498)
  a1t <- scale(a1t)[,1]
  a2t <- scale(a2t)[,1]
  path.sgarch1 <- try(ugarchpath(spec1, n.sim=(n+100), n.start=(N-n-100), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) ) )
  path.sgarch2 <- try(ugarchpath(spec2, n.sim=(n+100), n.start=(N-n-100), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) ) )
  
  ### occassionally we get an error here, so just retry
  while(inherits(path.sgarch1, "try-error") || inherits(path.sgarch2, "try-error")) {
    cat(counter, ":", sub.counter, ": Generation error\n", append=TRUE, file=error.file)
    cat("   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
    cat("   spec2: ", unlist(spec2@model$fixed.pars), "\n", append=TRUE, file=error.file)
    a1t <- rsged(N, mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)
    a2t <- rsged(N, mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498)
    a1t <- scale(a1t)[,1]
    a2t <- scale(a2t)[,1]
    path.sgarch1 <- try(ugarchpath(spec1, n.sim=(n+100), n.start=(N-n-100), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) ) )
    path.sgarch2 <- try(ugarchpath(spec2, n.sim=(n+100), n.start=(N-n-100), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) ) )
  }
  ## Data is generated...
  x1 <- as.vector(path.sgarch1@path$seriesSim)    ## GARCH part, the "errors" in the VAR...
  x2 <- as.vector(path.sgarch2@path$seriesSim)    ## of length n+100... need extra to generate VAR
  
  Innov <- cbind(x1, x2)
  p <- var.fit$order
  
  if(p > 0) {
    Phi <- aperm(var.fit$ar, c(2,3,1))
    ##  As before, the sigma is not used because we specify the errors
    ## which follow our GARCH processes
    X <- varima.sim(model=list(ar=Phi), innov=Innov, n=(n+100), k=ncol(Innov), sigma=var.fit$var.pred )
    X <- X[-(1:100),]
    fit <- ar(X, aic=FALSE, order.max=p)
    resid <- fit$resid
    x1 <- resid[-(1:p),1]   # Remove the NAs...
    x2 <- resid[-(1:p),2]
  } else {
    resid <- Innov
    x1 <- resid[,1]
    x2 <- resid[,2]
  }
  
  ##
  ## Fit with the correct models. Much like in the above, for smaller n
  ## we run into computational issues on occasion.
  
  fit1.success <- FALSE
  fit2.success <- TRUE    ## Nothing to fit, so a success!
  
  fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                      fit1.success = TRUE},
                    error = function(e) { cat(counter, "Fitting x1", file=error.file, append=TRUE); print(e) } )
  
  ##################
  ## If the fits are weird, generate new data!
  error1 <- FALSE;
  error2 <- FALSE;
  if((!fit1.success) || (!fit2.success)) error1 <- TRUE    ## This happens if the AR part is not-stationary
  else if ((sum(coef(x1.fit)[2:3])>=1)  ) error2 <- TRUE  ## GARCH stationary error
  else FALSE;
  while( error1 || error2 ) { 
    cat(counter, ":", sub.counter,  ": Boot Fit error\n", append=TRUE, file=error.file)
    cat("   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
    a1t <- rsged(N, mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)
    a2t <- rsged(N, mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498)
    a1t <- scale(a1t)[,1]
    a2t <- scale(a2t)[,1]
    path.sgarch1 <- try(ugarchpath(spec1, n.sim=(n+100), n.start=(N-n-100), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) ) )
    path.sgarch2 <- try(ugarchpath(spec2, n.sim=(n+100), n.start=(N-n-100), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) ) )
    while(inherits(path.sgarch1, "try-error") || inherits(path.sgarch2, "try-error")) {
      cat(counter, ":", sub.counter, ": generate error after fit\n", append=TRUE, file=error.file)
      cat("   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
      cat("   spec2: ", unlist(spec2@model$fixed.pars), "\n", append=TRUE, file=error.file)
    }
    
    x1 <- as.vector(path.sgarch1@path$seriesSim)
    x2 <- as.vector(path.sgarch2@path$seriesSim)
    
    Innov <- cbind(x1, x2)
    p <- 4
    Phi <- aperm(var.fit$ar, c(2,3,1))
    X <- varima.sim(model=list(ar=Phi), innov=Innov, n=(n+100), k=ncol(Innov), sigma=var.fit$var.pred )
    X <- X[-(1:100),]
    fit <- ar(X, aic=FALSE, order.max=p)
    resid <- fit$resid
    x1 <- resid[-(1:p),1]   # Remove the NAs...
    x2 <- resid[-(1:p),2]
    
    fit1.success <- FALSE
    fit2.success <- TRUE 
    fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                        fit1.success = TRUE},
                      error = function(e) { cat(counter,":",sub.counter, " Fitting x1 in bootstrap error: ", file=error.file, append=TRUE); cat(as.character.error(e), file=error.file, append=TRUE) } )
    
    if( (!fit1.success) ) error1 <- TRUE
    else if ((sum(coef(x1.fit)[2:3])>=1)  ) error2 <- TRUE
    else {
      cat(counter, ":", sub.counter, " New Good Fits\n", append=TRUE, file=error.file)
      cat(counter, ":", sub.counter, " x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
      error1 <- FALSE
      error2 <- FALSE
    }
    sub.counter <<- sub.counter + 1  
  }
  ## save the standardized residuals
  x1.res <- x1.fit@residuals/sqrt(x1.fit@h.t)
  #x2.res <- x2.fit@residuals/sqrt(x2.fit@h.t)
  x2.res <- scale(x2)[,1]

  if(length(M)==1){
    tmp <- get.test.stats(r1=x1.res, r2=x2.res, M=M, one.way=one.way)
    out <- tmp[,1]
  } else {
    tmp1 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[1], one.way=one.way)
    tmp2 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[2], one.way=one.way)
    tmp3 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[3], one.way=one.way)
    out <- c(tmp1[,1], tmp2[,1], tmp3[,1] )
  }
  out
}

