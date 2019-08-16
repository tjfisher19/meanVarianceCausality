#################
## Nothing fancy here, just some functions we use quite a bit in the simulation
##
## get.test.stats is a wrapper for all the test we consider
## The output is a big matrix following the form
##   BP.stat,    BP.pval
##   BP*.stat,   BP*.pval
##   LB.stat,    LB.pval
##   LB*.stat,   LB*.pval
##   FG.stat,    FG.pval
##   FG*.stat,   FG*.pval
##   Hong.stat,  Hong.pval
##   Hong*.stat, Hong*.pval
##   Det.stat,   Det.pval    (optional, matrix=TRUE)
##   Det*.stat,  Det*.pval   (optional, matrix=TRUE)
## 
## where the * represents a test for connection in variance
##
## The inputs are
##   r1, r2  - the two (residual) series
##         M - The maximum lag to consider
##   one.way - Whether to consider a one-way test or not
##    matrix - Do we include the matrix determinant test... as M grows it can
##             slow down the simulations
get.test.stats <- function(r1, r2, M, one.way="no", matrix=TRUE) {
  out <- non.matrix.test(r1, r2, lag.max=M, one.way=one.way)
  if(matrix)
    out <- rbind(out,
                 Determinant.Test(r1, r2, log.det=TRUE, lag.max=M, one.way=one.way),
                 Determinant.Test(r1*r1, r2*r2, log.det=TRUE, lag.max=M, one.way=one.way) );
  out
}  

#########################
## Find the bootstrapped p-value
## for a given test...
get.bootstrap.pval <- function(i, null.distro, test.stats) {
  (sum(null.distro[i,]>test.stats[i])+1)/(length(null.distro[i,])+1)
}


############################
## Given a simulated exmpiral distribution, this finds some summary stats
## simply a wrapper function...
get.null.distro.stats <- function(i, null.distro) {
  tab <- c(mean(null.distro[i,]),
           sd(null.distro[i,]),
           quantile(null.distro[i,], probs=c(0.90,0.95,0.99)) )
  names(tab) <- c("Mean", "Std.Dev.", "90th", "95th", "99th")
  tab
}
