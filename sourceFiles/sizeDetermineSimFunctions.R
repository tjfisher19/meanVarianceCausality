##########################
# 
## A simple code that generates iid Normal data
## and calculates the emprical size of the statistics
## We do this over a bunch of sample sizes for Figure 2
## in the paper.


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


do.one.trivial.noboot.simulation <- function(n=100, M=5, one.way="no") {
  #############
  ## Some counting so we know we are making progress
  
  #### Generate the data
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
    out <- pvals1
  } else {
    tmp1 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[1], one.way=one.way)
    tmp2 <- get.test.stats(r1=x1.res, r2=x2.res, M=M[2], one.way=one.way)
    stats <- c(tmp1[,1], tmp2[,1] )
    pvals1 <- tmp1[,2]
    pvals2 <- tmp2[,2]
    out <- c(pvals1, pvals2)
  }

  out
}



