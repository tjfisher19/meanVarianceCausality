##########################
### 
### Modification of the bootstrapSimFunctions.R code
###
### It says trivial because everything is simple here...
### 
### Data is iid standard normals, not fitting anything!
### I took the bootstrapSimFunctions.R code and edited...
### This consisted of mostly shortening it up
### This is the simulation in the motivating section of the paper
###  Section 2.2 I believe...

generateNormalData <- function(n=100) {
  #########################
  ## Generate underlying error series
  N <- n + 1000
  a1t <- rnorm(N)
  a2t <- rnorm(N)

  a1t <- a1t/sd(a1t)
  a2t <- a2t/sd(a2t)
  
  list(x1=a1t[1001:N], x2=a2t[1001:N])
}


do.one.trivial.simulation <- function(cluster=NULL, n=100, M=5, boot.run=1000, one.way="no") {
  #############
  ## Some counting so we know we are making progress
  counter <<- counter + 1
  if(counter%%5==0) {
    time <- (proc.time()-ptm)[3]
    cat(counter, "  ", time, "\n", file=status.file, append=TRUE)
  }
  
  #### Generate the data
  cat(counter, ": generating initial data", file=status.file, append=TRUE);
  tmp <- generateNormalData(n=n)
  x1 <- tmp$x1
  x2 <- tmp$x2
  
  #####################
  ## get standardized residuals
  x1.res <- scale(x1)[,1]
  x2.res <- scale(x2)[,1]
  
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

  cat(counter, ": About to start the bootstrap distro code\n", file=status.file, append=TRUE);

  sub.counter <<- 0;
  if(is.null(cluster)) {
    null.distro <- sapply(rep(n,boot.run), one.normal.iteration, M=M, one.way=one.way )
    out <- c(sapply(1:10, get.boostrap.pval, null.distro=null.distro, test.stats=stats),
             pvals1) 
  } else {
    null.distro <- parSapply(cl=cluster, X=rep(n,boot.run), FUN=one.normal.iteration,  
                             M=M, one.way=one.way )

    out <- c(sapply(1:10, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals1,
             sapply(11:20, get.bootstrap.pval, null.distro=null.distro, test.stats=stats),
             pvals2 )
  }

  out
}


one.normal.iteration <- function(n, M, one.way="no") {
  N <- n+1000
  #cat("here", file=status.file, append=TRUE)
 
  x1 <- rnorm(N)[1001:N]
  x2 <- rnorm(N)[1001:N]
  
  x1.res <- scale(x1)[,1]
  x2.res <- scale(x2)[,1]
  
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

