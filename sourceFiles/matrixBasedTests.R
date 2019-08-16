###############################################################
# This file has the matrix based test
#
# The function to use is Determinant.test... see below
#
# Also note much of this code implements results in 
# the Robbins & Fisher (2015) paper, hence the function names
#
#
# RobbinsBlockMatrix builds the Block Cross-Correlation matrix
# we are using given series x & y, the maximum lag to construct the
# matrix lag.max. 
#
# The RobbinsBlockSingle and RobbinsBlockOneWay are special cases of the
# first function. In 'Single' the matrix is for a single lag (say lag 0)
# while in OneWay it builds a special block where the lower or upper off-diagonal
# is set to 0, for a one-directional test
#
# Determinant.Test - This is our new test. Given series x & y and maximum lag lag.max. It
#     will do the test based on the determinant. The boolean log.det tells the function 
#     to do a Mahdi McLeod type determinant test (log.det=TRUE) or a Pena Rodriguez type
#     determinant test (log.det=FALSE). The former tends to be more powerful but has a 
#     very slow convergence rate to the asymptotic distribution.
#
# Output for the function is a vector of length two, the first element
#     is the value of the test statistic, the second element is its corresponding p-value.

RobbinsBlockMatrix <- function (x, y, lag.max) 
{
  x <- as.matrix(x);
  y <- as.matrix(y);
  d.x <- NCOL(x);
  d.y <- NCOL(y);
  k <- d.x + d.y;
  n <- NROW(x);  # We can error check later
  m <- lag.max + 1
  out <- matrix(0, nrow = k*m, ncol = k*m)
  X <- cbind(x,y);
  Accmat <- acf(X, lag.max = lag.max, plot = FALSE, type = "covariance", na.action=na.pass)$acf
  C0 <- Accmat[1, ,];
  C0[(d.x+1):k, 1:(d.x)] <- 0;
  C0[1:(d.x), (d.x+1):k] <- 0;
  
  C0.inv <- solve(C0);
  
  L <- t(chol(C0.inv));
  for (i in 0:lag.max) for (j in i:lag.max) {
    tmp1 <- crossprod(t(crossprod(L, t(Accmat[j - i + 1, , ]))), L)
    tmp2 <- crossprod(t(crossprod(L, Accmat[j - i + 1, , ])), L)
    if(i!=j) {
      tmp1[1:d.x, 1:d.x] <- 0;
      tmp1[(d.x+1):k, (d.x+1):k] <- 0;
      tmp2[1:d.x, 1:d.x] <- 0;
      tmp2[(d.x+1):k, (d.x+1):k] <- 0;
    }
    
    out[(j * k + 1):(k * (j + 1)), (i * k + 1):(k * (i + 1))] <- tmp1;
    out[(i * k + 1):(k * (i + 1)), (j * k + 1):(k * (j + 1))] <- tmp2;
  }
  return(out);
}


RobbinsBlockSingle <- function(x, y, lag.spot )
{
  x <- as.matrix(x);
  y <- as.matrix(y);
  d.x <- NCOL(x);
  d.y <- NCOL(y);
  k <- d.x + d.y;
  n <- NROW(x);  # We can error check later

  out <- matrix(0, nrow = k, ncol = k)
  out <- diag(k);
  X <- cbind(x,y);
  Accmat <- acf(X, lag.max = abs(lag.spot), plot = FALSE, type = "covariance", na.action=na.pass)$acf
  C0 <- Accmat[1, ,];
  C01 <- C0[1:d.x, 1:d.x];
  C02 <- C0[(d.x+1):k, (d.x+1):k];
  C0[(d.x+1):k, 1:(d.x)] <- 0;
  C0[1:(d.x), (d.x+1):k] <- 0;
  L01 <- t(chol(solve(C01)));
  L02 <- t(chol(solve(C02)));
  C0.inv <- solve(C0);
  
  L <- t(chol(C0.inv));
  if(lag.spot >= 0)
    Ch <- Accmat[abs(lag.spot)+1, 1:d.x, (d.x+1):k]
  else
    Ch <- t(Accmat[abs(lag.spot)+1, (d.x+1):k, 1:d.x] );
  tmp <- sqrt(n/(n-abs(lag.spot)))*crossprod(t(crossprod(L01,Ch)), L02);
  out <- diag(k);
  out[1:d.x,(d.x+1):k]=tmp;
  out[(d.x+1):k,1:(d.x)]=t(tmp);
  return(out);
}


RobbinsBlockOneWay <- function(x, y, lag.max, right=TRUE)
{
  x <- as.matrix(x);
  y <- as.matrix(y);
  d.x <- NCOL(x);
  d.y <- NCOL(y);
  k <- d.x + d.y;
  n <- NROW(x);  # We can error check later
  m <- lag.max + 1
  out <- matrix(0, nrow = k*m, ncol = k*m)
  X <- cbind(x,y);
  Accmat <- acf(X, lag.max = lag.max, plot = FALSE, type = "covariance", na.action=na.pass)$acf
  C0 <- Accmat[1, ,];
  C0[(d.x+1):k, 1:(d.x)] <- 0;
  C0[1:(d.x), (d.x+1):k] <- 0;
  
  C0.inv <- solve(C0);  
  L <- t(chol(C0.inv));
  for (i in 0:lag.max) for (j in i:lag.max) {
    tmp1 <- crossprod(t(crossprod(L, t(Accmat[j - i + 1, , ]))), L)
    tmp2 <- crossprod(t(crossprod(L, Accmat[j - i + 1, , ])), L)
    if( (i==j) ) {
      tmp1[1:d.x, (d.x+1):k] <- 0;
      tmp1[(d.x+1):k, 1:d.x] <- 0;
      tmp2[1:d.x,(d.x+1):k] <- 0;
      tmp2[(d.x+1):k, 1:d.x] <- 0;
    }
    if( (i!=j && !right)  ) {
      tmp1[1:k, 1:d.x] <- 0;
      tmp1[(d.x+1):k, (d.x+1):k] <- 0;
      tmp2[1:d.x, 1:k] <- 0;
      tmp2[(d.x+1):k, (d.x+1):k] <- 0;      
      #tmp2[1:k, 1:k] <- 0;
    } else if(i!=j) {
      tmp1[1:d.x, 1:k] <- 0;
      tmp1[(d.x+1):k, (d.x+1):k] <- 0;
      tmp2[1:k, 1:d.x] <- 0;
      tmp2[(d.x+1):k, (d.x+1):k] <- 0; 
      #tmo1[1:k, 1:k] <- 0;
    }
    
    out[(j * k + 1):(k * (j + 1)), (i * k + 1):(k * (i + 1))] <- tmp1;
    out[(i * k + 1):(k * (i + 1)), (j * k + 1):(k * (j + 1))] <- tmp2;
  }
  return(out);
}


Determinant.Test <- function(x,y, log.det=TRUE, lag.max=5, single.lag=FALSE, one.way=c("no", "right", "left")) {
  one.way <- match.arg(one.way);
  if(single.lag && one.way!="no")
    stop("Cannot do single lag & one-way at same time");
  
  if(one.way=="right") {
    mat <- RobbinsBlockOneWay(x,y, lag.max=lag.max, right=TRUE);
  } else if(one.way=="left") {
    mat <- RobbinsBlockOneWay(x,y, lag.max=lag.max, right=FALSE);
  } else if(single.lag) {
    mat <- RobbinsBlockSingle(x,y, lag.spot=lag.max);
  } else {
    mat <- RobbinsBlockMatrix(x, y, lag.max=lag.max);
  }
  x <- as.matrix(x);
  y <- as.matrix(y);
  n <- NROW(x);
  d.x <- NCOL(x);
  d.y <- NCOL(y);
  if(lag.max==0 || single.lag) {
    weights <- 1;
    lag.max <- 0;
  } else if(one.way!="no") {
    weights <- (lag.max:1)/(lag.max+1)
  } else {
    weights <- c(1:lag.max/(lag.max+1), 1, lag.max:1/(lag.max+1) );
  }
  
  if(log.det) {
    stat <- -n*log(abs(det(mat)));
    weights <- weights*(lag.max+1);
  } else {
    stat <- n*(1-(abs(det(mat)))^(1/(lag.max+1) ) );
  }  
  
  cumulant1 <- d.x*d.y*sum(weights);
  cumulant2 <- 2*d.x*d.y*sum(weights^2);
  shape <- cumulant1^2/cumulant2;
  scale <- cumulant2/cumulant1;
  pval <- pgamma(stat, shape=shape, scale=scale, lower.tail=FALSE);
  c(stat, pval);
}

