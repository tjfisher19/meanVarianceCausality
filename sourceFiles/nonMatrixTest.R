#########################
## This function calculated all the non-matrix based tests
##
## This includes Q_{BP}, Q_{LB}, Z_m, Q_{FG}
##
## For each statistic it outputs the stat value and corresponding p-value
## first for the test for the mean, then the test for the variance.
##
## rbind(c(bp.stat1, bp.pval1),
##       c(bp.stat2, bp.pval2),
##       c(lb.stat1, lb.pval1),
##       c(lb.stat2, lb.pval2),
##       c(fg.stat1, fg.pval1),
##       c(fg.stat2, fg.pval2),
##       c(dan.stat1, dan.pval1),
##       c(dan.stat2, dan.pval2) )
##
## The inputs to the function are
##      x,y - the two (residual) series
##  lag.max - the maximum lag m to consider
##  one.way - If you wish to do a one-sided test, whether right (+) or left (-)
non.matrix.test <- function(x, y, lag.max=5, one.way=c("no", "right", "left") )
{
  one.way <- match.arg(one.way);
  if( (lag.max==0) )
    weight.type <-  "none"
  else
    lag.max <- abs(lag.max);
  
  #########################
  ## Get cross-correlations of series & squared series
  n <- length(x)
  r <- ccf(x,y, lag.max=lag.max, type="correlation", plot=FALSE, na.action=na.pass)$acf
  r1 <- r*r
  r <- ccf((x*x), (y*y), lag.max=lag.max, type="correlation", plot=FALSE, na.action=na.pass)$acf
  r2 <- r*r
  
  ###
  ## If one.way=="right" or "left" only need have the values
  if(one.way=="right") {
    ind <- 1:lag.max
    r1 <- r1[(lag.max+2):(lag.max*2 + 1)]
    r2 <- r2[(lag.max+2):(lag.max*2 + 1)]
  }
  else if(one.way=="left") {
    ind <- -lag.max:(-1)
    r1 <- r1[1:lag.max]
    r2 <- r2[1:lag.max]
  }
  else
    ind <- -lag.max:lag.max;
  
  #####################
  ## Daniell kernel... this can be defined in a few different ways
  daniell <- function(z) {
    if(z==0)
      1
    else
      sin(pi*z)/(pi*z);   #*(abs(z)<=1);
  }
  dan.weights <- sapply( (ind/lag.max), daniell);
  dan.weights <- dan.weights*dan.weights;    ### Weights for Hong(1996) stat

  ljung.weights <- n/(n-abs(ind) );    ### Weights for Ljung-Box(1978)
  
  ### Weights for Fisher-Gallagher type, see Robbins & Fisher (2015)
  if(one.way=="no")
    fg.weights1 <- c(1:lag.max/(lag.max+1), 1, lag.max:1/(lag.max+1) )
  else if(one.way=="right")
    fg.weights1 <- c(lag.max:1/(lag.max+1) )
  else if(one.way=="left")
    fg.weights1 <- c(1:lag.max/(lag.max+1) )
  fg.weights <- fg.weights1*n/(n-abs(ind) );
  
  
  ### Box-Pierce "weights", all equal 1
  bp.weights <- 1
  
  dan.stat1 <- n*sum(dan.weights*r1)
  bp.stat1 <- n*sum(bp.weights*r1)
  lb.stat1 <- n*sum(ljung.weights*r1)
  fg.stat1 <- n*sum(fg.weights*r1)
  
  dan.stat2 <- n*sum(dan.weights*r2)
  bp.stat2 <- n*sum(bp.weights*r2)
  lb.stat2 <- n*sum(ljung.weights*r2)
  fg.stat2 <- n*sum(fg.weights*r2)
  
  #######
  ## Hong(1996) standardized statistic
  Sn <- sum( (n-abs(ind))/n*dan.weights);
  Dn <- 2*sum( (n-abs(ind))/n*(n-(abs(ind)+1))/n*dan.weights^2);
  dan.stat1 <- (dan.stat1-Sn)/sqrt(Dn);
  dan.stat2 <- (dan.stat2-Sn)/sqrt(Dn);
  dan.pval1 <- pnorm(dan.stat1, lower.tail=FALSE);
  dan.pval2 <- pnorm(dan.stat2, lower.tail=FALSE);
  
  ###
  ## Now the BP and LB types, just chi square stats
  if(one.way=="no") {
    bp.pval1 <- pchisq(bp.stat1, df=(2*lag.max+1), lower.tail=FALSE )
    bp.pval2 <- pchisq(bp.stat2, df=(2*lag.max+1), lower.tail=FALSE )
    lb.pval1 <- pchisq(lb.stat1, df=(2*lag.max+1), lower.tail=FALSE );
    lb.pval2 <- pchisq(lb.stat2, df=(2*lag.max+1), lower.tail=FALSE );
  } else {
    bp.pval1 <- pchisq(bp.stat1, df=(lag.max), lower.tail=FALSE )
    bp.pval2 <- pchisq(bp.stat2, df=(lag.max), lower.tail=FALSE )
    lb.pval1 <- pchisq(lb.stat1, df=(lag.max), lower.tail=FALSE );
    lb.pval2 <- pchisq(lb.stat2, df=(lag.max), lower.tail=FALSE );    
  }
  
  ########################
  ## Now the FG type, can be approximated as a Gamma rv.
  k1 <- sum(fg.weights1);
  k2 <- 2*sum(fg.weights1^2);
  shape <- k1^2/k2;
  scale <- k2/k1;
  fg.pval1 <- pgamma(fg.stat1, shape=shape, scale=scale, lower.tail=FALSE);
  fg.pval2 <- pgamma(fg.stat2, shape=shape, scale=scale, lower.tail=FALSE);
  
  ################################
  ## report all the test stats and p-values
  rbind(c(bp.stat1, bp.pval1),
        c(bp.stat2, bp.pval2),
        c(lb.stat1, lb.pval1),
        c(lb.stat2, lb.pval2),
        c(fg.stat1, fg.pval1),
        c(fg.stat2, fg.pval2),
        c(dan.stat1, dan.pval1),
        c(dan.stat2, dan.pval2) )
}
