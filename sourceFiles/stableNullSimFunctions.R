####################
## The main call function that performs the robustness studies in the paper
##
## A modification of bootstrapSimFunctions.R
##
## The word "stable" is used because we originally wrote these functions to work with
## alpha-stable errors but scratched that idea as the paper progressed
## (alpha-stable does not have finite moments generally). 
## 
generateStableNullData <- function(n=100) {
  #########################
  ## Generate underlying error series
  N <- n + 1000
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
                      mean.model=list(armaOrder=c(1,2), include.mean=FALSE),
                      fixed.pars=list(omega=1.1171e-05, alpha1=0.1609, beta1=0.7079, ar1=-0.1795, ma1=0.1453, ma2=-0.1058) )
  spec2 <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(0,0)),
                      mean.model=list(armaOrder=c(0,1), include.mean=FALSE),
                      fixed.pars=list(ma1=0.35, omega=3.1135e-05) )

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
do.one.var.stable.null.simulation <- function(cluster=NULL, n=100, M=5, boot.run=1000, one.way="no") {

  #############
  ## Some counting so we know we are making progress
  counter <<- counter + 1
  sub.counter <<- 0
  if(counter%%5==0) {
    time <- (proc.time()-ptm)[3]
    cat(counter, "  ", time, "\n", file=status.file, append=TRUE)
  }
  
  #### Generate the data
  cat(counter, ": generating initial data", file=status.file, append=TRUE);
  tmp <- generateStableNullData(n=n)
  x1 <- tmp$x1
  x2 <- tmp$x2
  
  cat(" -- done\n", file=status.file, append=TRUE)
    
  cat(counter, " -- attempting model fits", file=status.file, append=TRUE);
    
  model1 <- ~arma(1,2)+garch(1,1)
  #model2 <- ~arma(0,1)+garch(0,0)
  model2 <- c(0,0,1)
  #####################
  ## Fit data with specified models, ARMA(1,1)+GARCH(1,1)
  ##
  ## With smaller sample sizes, we run into issues of the fitted models not being
  ## stationary and/or having other numerical issues
  ##
  ## We are interested in stationary fits, so use some exception handling
  ## If there is an error with a fit. Generate new data and fit it.

  fit1.success <- FALSE
  fit2.success <- FALSE

  fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                           fit1.success = TRUE},
            error = function(e) { cat(counter, "Fitting x1:", file=error.file, append=TRUE); cat(as.character.error(e), file=error.file, append=TRUE) } )

  fit2 <- tryCatch( { x2.fit= suppressWarnings(Arima(x2, order=model2, include.mean=FALSE) )
                      fit2.success = TRUE},
            error = function(e) { cat(counter, "Fitting x2:", file=error.file, append=TRUE); cat(as.character.error(e), file=error.file, append=TRUE) } )

  ##################
  ## If the fits are weird, generate new data!
  error1 <- FALSE;
  error2 <- FALSE;
  if((!fit1.success) || (!fit2.success)) error1 <- TRUE    ## This happens if the AR part is not-stationary
  else if (sum(coef(x1.fit)[5:6])>=0.99)  error2 <- TRUE  ## GARCH stationary error
  else FALSE;
  while( error1 || error2 ) { 
    cat(" -- ERROR:", error1, " ", error2, file=status.file, append=TRUE);
    cat(counter, ": Initial Fits\n", append=TRUE, file=error.file)
    if(!error1) {
      cat("    x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
      cat("    x2.fit: ", coef(x2.fit), "\n", append=TRUE, file=error.file)
    }
    ## Generate a new set of data
    tmp <- generateStableNullData(n=n)
    x1 <- tmp$x1
    x2 <- tmp$x2
    
    #####################
    ## Fit data with specified models
    fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                        fit1.success = TRUE},
                      error = function(e) { cat(counter, "Fitting x1", file=error.file, append=TRUE); cat(as.character.error(e), file=error.file, append=TRUE) } )
    
    fit2 <- tryCatch( { x2.fit= suppressWarnings(Arima(x2, order=model2, include.mean=FALSE) )
                        fit2.success = TRUE},
                      error = function(e) { cat(counter, "Fitting x2", file=error.file, append=TRUE); cat(as.character.error(e), file=error.file, append=TRUE) } )
    
    if((!fit1.success) || (!fit2.success)) error1 <- TRUE
    else if ((sum(coef(x1.fit)[5:6])>=0.99)  ) error2 <- TRUE
    else {
      cat("New Fits\n", append=TRUE, file=error.file)
      cat(" x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
      cat(" x2.fit: ", coef(x2.fit), "\n", append=TRUE, file=error.file)
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
  x2.res <- scale(x2.fit$residuals)[,1]
  
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
  coef2 <- coef(x2.fit)
  name2 <- substring(names(coef2),1,2)

  fit.order.x1<-x1.fit@fit$series$order
  #fit.order.x2<-x2.fit@fit$series$order
  fit.order.x2 <- x2.fit$arma[1:2]
  # The code below is setup to handle an ARMA(3,3)+GARCH(3,3) but we only use (1,1)+(1,1)
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
  cat(counter, ": About to start the bootstrap distro code", file=status.file, append=TRUE);

  sub.counter <<- 0;
  if(is.null(cluster)) {   # run in sequence
    null.distro <- sapply(rep(n,boot.run), one.iteration.stable.null.boot, spec1=spec1, spec2=spec2, 
                          M=M, model1=model1, model2=model2, one.way=one.way )
    out <- c(sapply(1:10, get.boostrap.pval, null.distro=null.distro, test.stats=stats),
             pvals1) 
  } else {   # run in parallel
    null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.iteration.stable.null.boot, spec1=spec1, spec2=spec2, 
                             M=M, model1=model1, model2=model2, one.way=one.way )

    out <- c(sapply(1:10, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals1,
             sapply(11:20, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals2,
             sapply(21:30, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals3)
  }
  cat(" -- done\n", file=status.file, append=TRUE)

  out
}


####################################
## Here is the function that does the bootstrapping distro
##
## We are always generating data under the null distro so the
## two series are uncorrelated
one.iteration.stable.null.boot <- function(spec1, spec2, model1, model2, n, M, one.way="no") {
  N <- n+1000
  sub.counter <<- sub.counter + 1;
  a1t <- rsged(N, mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)
  a2t <- rsged(N, mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498)
  a1t <- scale(a1t)[,1]
  a2t <- scale(a2t)[,1]
  path.sgarch1 <- try(ugarchpath(spec1, n.sim=(n), n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) ) )
  path.sgarch2 <- try(ugarchpath(spec2, n.sim=(n), n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) ) )
  ### occassionally we get an error here, so just retry
  while(inherits(path.sgarch1, "try-error") || inherits(path.sgarch2, "try-error")) {
    cat(counter, ":", sub.counter, ": Generation error\n", append=TRUE, file=error.file)
    cat("   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
    cat("   spec2: ", unlist(spec2@model$fixed.pars), "\n", append=TRUE, file=error.file)
    a1t <- rsged(N, mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)
    a2t <- rsged(N, mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498)
    a1t <- scale(a1t)[,1]
    a2t <- scale(a2t)[,1]
    path.sgarch1 <- try(ugarchpath(spec1, n.sim=(n), n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) ) )
    path.sgarch2 <- try(ugarchpath(spec2, n.sim=(n), n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) ) )
  }
  ## Data is generated...
  x1 <- as.vector(path.sgarch1@path$seriesSim)    ## GARCH part, the "errors" in the VAR...
  x2 <- as.vector(path.sgarch2@path$seriesSim)    ## of length n+100... need extra to generate VAR
  
  ##
  ## Fit with the correct models. Much like in the above, for smaller n
  ## we run into computational issues on occasion.
  fit1.success <- FALSE
  fit2.success <- FALSE

  fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                           fit1.success = TRUE},
            error = function(e) { cat(counter,":", sub.counter, "  Fitting x1 in bootstrap error: ", file=error.file, append=TRUE); cat(as.character.error(e), file=error.file, append=TRUE) } )

  fit2 <- tryCatch( { x2.fit= Arima(x2, order=model2, include.mean=FALSE)
                      fit2.success = TRUE},
            error = function(e) { cat(counter,":",sub.counter, "  Fitting x2 in bootstrap error: ", file=error.file, append=TRUE); cat(as.cahracter.error(e), file=error.file, append=TRUE) } )

  ##################
  ## If the fits are weird, generate new data!
  error1 <- FALSE;
  error2 <- FALSE;
  if((!fit1.success) || (!fit2.success)) error1 <- TRUE    ## This happens if the AR part is not-stationary
  else if (sum(coef(x1.fit)[5:6])>=1)  error2 <- TRUE  ## GARCH stationary error
  else FALSE;
  while( error1 || error2 ) {
    #cat(counter, ":", sub.counter, " -- Boot fit error:", error1, " ", error2, "\n", file=status.file, append=TRUE);
    cat(counter, ":", sub.counter, " -- Boot fit error:", error1, " ", error2, "\n", file=error.file, append=TRUE);
    cat(counter, ":", sub.counter, ": Initial Fits, n=", n, "\n", append=TRUE, file=error.file)
    if(!error1) {
      cat(counter, ":", sub.counter,"    x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
      cat(counter, ":", sub.counter,"    x2.fit: ", coef(x2.fit), "\n", append=TRUE, file=error.file)
    }
    cat(counter,":", sub.counter,"   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
    cat(counter,":", sub.counter,"   spec2: ", unlist(spec2@model$fixed.pars), "\n", append=TRUE, file=error.file)
    #cat(counter,":", sub.counter," x1.mod: ", as.character(model1), "\n", append=TRUE, file=error.file)
    #cat(counter,":", sub.counter," x2.mod: ", model2, "\n", append=TRUE, file=error.file)

    a1t <- rsged(N, mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)
    a2t <- rsged(N, mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498)
    a1t <- scale(a1t)[,1]
    a2t <- scale(a2t)[,1]
    path.sgarch1 <- try(ugarchpath(spec1, n.sim=(n), n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) ) )
    path.sgarch2 <- try(ugarchpath(spec2, n.sim=(n), n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) ) )
    while(inherits(path.sgarch1, "try-error") || inherits(path.sgarch2, "try-error")) {
      cat(counter, ":", sub.counter, ": generate error after fit\n", append=TRUE, file=error.file)
      cat("   spec1: ", unlist(spec1@model$fixed.pars), "\n", append=TRUE, file=error.file)
      cat("   spec2: ", unlist(spec2@model$fixed.pars), "\n", append=TRUE, file=error.file)
      a1t <- rsged(N, mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)
      a2t <- rsged(N, mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498)
      a1t <- scale(a1t)[,1]
      a2t <- scale(a2t)[,1]
      path.sgarch1 <- try(ugarchpath(spec1, n.sim=(n), n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a1t, ncol=1)) ) )
      path.sgarch2 <- try(ugarchpath(spec2, n.sim=(n), n.start=(N-n), m.sim=1, startMethod="sample", custom.dist=list(name="sample", distfit=matrix(a2t, ncol=1)) ) )
      
    }
    
    x1 <- as.vector(path.sgarch1@path$seriesSim)
    x2 <- as.vector(path.sgarch2@path$seriesSim)
   
    fit1.success <- FALSE
    fit2.success <- FALSE 
    fit1 <- tryCatch( { x1.fit=suppressWarnings(garchFit(formula=model1, data=x1, trace=FALSE, include.mean=FALSE, cond.dist="QMLE") )
                               fit1.success = TRUE},
               error = function(e) { cat(counter,":",sub.counter, " Fitting x1 in bootstrap error: ", file=error.file, append=TRUE); cat(as.character.error(e), file=error.file, append=TRUE) } )

    fit2 <- tryCatch( { x2.fit= Arima(x2, order=model2, include.mean=FALSE)
                               fit2.success = TRUE},
                error = function(e) { cat(counter,":",sub.counter, " Fitting x2 in bootstrap error: ", file=error.file, append=TRUE); cat(as.cahracter.error(e), file=error.file, append=TRUE) } )

    if((!fit1.success) || (!fit2.success)) error1 <- TRUE
    else if ((sum(coef(x1.fit)[5:6])>=1)  ) error2 <- TRUE
    else {
      cat(counter, ":", sub.counter, " New Good Fits\n", append=TRUE, file=error.file)
      cat(counter, ":", sub.counter, " x1.fit: ", coef(x1.fit), "\n", append=TRUE, file=error.file)
      cat(counter, ":", sub.counter, " x2.fit: ", coef(x2.fit), "\n", append=TRUE, file=error.file)
      error1 <- FALSE
      error2 <- FALSE
    }
    sub.counter <<- sub.counter + 1
  }
  ## save the standardized residuals
  x1.res <- x1.fit@residuals/sqrt(x1.fit@h.t)
  x2.res <- scale(x2.fit$residuals)[,1]
  
  if(length(M)==1){
    tmp <- get.test.stats(r1=x1.res, r2=x2.res, M=M, one.way=one.way)
    out <- tmp[,1]
  } else {
    tmp1 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[1], one.way=one.way)
    tmp2 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[2], one.way=one.way)
    tmp3 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[3], one.way=one.way)
    out <- c(tmp1[,1], tmp2[,1], tmp3[,1])
  }
  out
}

