###############################
# A modification to varima.sim from portes library
# Now we can input an innovation sequence
# Giving greater control of any cross correlation
# 
# Note the code still requires a sigma matrix but it is not used when the 
# innovations are supplied.
new.varima.sim=function (phi = NULL, theta = NULL, d = NA, sigma, n, innovations=NA, constant = NA,
                         trend = NA, demean = NA, StableParameters = NA, Trunc.Series = NA)
{
  if (!is.null(phi) && class(phi) != "array" && class(phi) !=
        "numeric")
    stop("Phi must be entered as NULL or array with dimension (k*k*p) or numeric")
  if (!is.null(theta) && class(theta) != "array" && class(theta) !=
        "numeric")
    stop("Theta must be entered as NULL or array with dimension (k*k*q) or numeric")
  sigma <- as.matrix(sigma)
  k <- NCOL(sigma)
  if (all(is.na(d)) || all(d == 0))
    d <- rep(0, k)
  if (length(d) != k)
    stop("d must be entered as a vector with length equal to number of sigma rows")
  if (any(d < 0))
    stop("number of differences must be a nonnegative integer/integers")
  if (all(is.na(constant)))
    constant <- rep(0, k)
  if (length(constant) != k)
    stop("constant must be entered as a vector with length equal to number of sigma rows")
  if (all(is.na(trend)))
    trend <- rep(0, k)
  if (length(trend) != k)
    stop("trend must be entered as a vector with length equal to number of sigma rows")
  if (all(is.na(demean)))
    demean <- rep(0, k)
  if (length(demean) != k)
    stop("demean must be entered as a vector with length equal to number of sigma rows")
  if (class(phi) == "numeric")
    phi <- array(phi, dim = c(1, 1, length(phi)))
  if (class(theta) == "numeric")
    theta <- array(theta, dim = c(1, 1, length(theta)))
  if (all(phi == 0))
    phi <- NULL
  if (all(theta == 0))
    theta <- NULL
  p <- ifelse(is.null(phi), 0, dim(phi)[3])
  q <- ifelse(is.null(theta), 0, dim(theta)[3])
  if (p > 0 && ((dim(phi)[1] != dim(phi)[2]) || dim(phi)[2] !=
                  k))
    stop("Wrong dimensions of phi or/and sigma")
  if (q > 0 && ((dim(theta)[1] != dim(theta)[2]) || dim(theta)[2] !=
                  k))
    stop("Wrong dimensions of theta or/and sigma")
  StableQ <- all(!is.na(StableParameters))
  if (StableQ) {
    StableParameters <- matrix(StableParameters, nrow = k)
    stopifnot(NCOL(StableParameters) == 4)
    ALPHA <- StableParameters[, 1]
    BETA <- StableParameters[, 2]
    GAMMA <- StableParameters[, 3]
    DELTA <- StableParameters[, 4]
  }
  if (is.na(Trunc.Series))
  {Trunc.Series <- min(100, ceiling(n/3))}
  r <- max(p, q)
  if(!is.na(innovations[1]))
  {
    innovations <- as.matrix(innovations)
    if(nrow(innovations)<ncol(innovations))
    {
      innovations<-t(innovations)
      warning("Took transpose of innovations.")
    }
    if(p==0){dum<-0}
    else{dum<-Trunc.Series}
    
    if(n+dum+r > nrow(innovations))
    {
      cat("Num. innovations required n+dum+r with= ",n,dum,r," = ",n+dum+r,", Num. innovations provided = ",nrow(innovations),"\n",sep="",file=error.file, append=TRUE)
      stop("You need more innovations.  Error 1.")
    }
  }
  if (p == 0) {
    if (StableQ)
      epsilon <- rstable(n + q, ALPHA, BETA, GAMMA, DELTA)
    else {
      if(is.na(innovations[1])){epsilon <- t(crossprod(chol(sigma), matrix(rnorm(k * (n + q)), ncol = n + q)))}
      else{
        if(nrow(innovations)<(n+q)){stop("You need more innovations: Error 2")}
        epsilon <- innovations[1:(n+q),]
      }
    }
    if (q == 0) {
      SlopDrift <- t(matrix(trend * rep(1:NROW(epsilon),
                                        each = k), nrow = k, ncol = NROW(epsilon)))
      DriftTerm <- sweep(SlopDrift, 2L, -constant, check.margin = FALSE)
      CenterData <- scale(epsilon, center = -demean, scale = FALSE)
      Sim.Series <- DriftTerm + CenterData
      for (i in 1:length(d)) {
        if (d[i] > 0)
          Sim.Series[, i] <- as.matrix(diffinv(Sim.Series[,
                                                          i], differences = d[i]))[-(1:d[i]), ]
        else Sim.Series[, i] <- Sim.Series[, i]
      }
      return(ts(Sim.Series))
    }
    else InvertQ(theta)
    psi <- array(c(diag(k), -theta), dim = c(k, k, q + 1))
    Sim.VMA <- vma.sim(psi = psi, a = epsilon)
    SlopDrift <- t(matrix(trend * rep(1:NROW(Sim.VMA), each = k),
                          nrow = k, ncol = NROW(Sim.VMA)))
    DriftTerm <- sweep(SlopDrift, 2L, -constant, check.margin = FALSE)
    CenterData <- scale(Sim.VMA, center = -demean, scale = FALSE)
    Sim.Series <- DriftTerm + CenterData
    for (i in 1:length(d)) {
      if (d[i] > 0)
        Sim.Series[, i] <- as.matrix(diffinv(Sim.Series[,
                                                        i], differences = d[i]))[-(1:d[i]), ]
      else Sim.Series[, i] <- Sim.Series[, i]
    }
    return(ts(Sim.Series))
  }
  
  FirstSim.Series <- matrix(numeric(0), nrow = n, ncol = k)
  
  psi <- ImpulseVMA(phi = phi, theta = theta, Trunc.Series = Trunc.Series)
  if (StableQ)
  {epsilon <- rstable(Trunc.Series + r, ALPHA, BETA, GAMMA,
                      DELTA)}
  else {
    if(is.na(innovations[1])){epsilon <- t(crossprod(chol(sigma), matrix(rnorm(k * (Trunc.Series + r)), ncol = Trunc.Series + r)))}
    else{epsilon <- innovations[1:(Trunc.Series + r),]
         if(nrow(innovations)<(Trunc.Series + r)){stop("You need more innovations: Error 3")}
         innovations<-innovations[-(1:(Trunc.Series + r)),]}
  }
  FirstSim.Series[1:r, ] <- vma.sim(psi = psi, a = epsilon)
  a <- matrix(epsilon[1:r, ], nrow = r, ncol = k)
  if (StableQ)
  {epsilon <- rbind(a, rstable(n, ALPHA, BETA, GAMMA, DELTA))}
  else {
    if(is.na(innovations[1])){
      epsilon <- rbind(a, t(crossprod(chol(sigma), matrix(rnorm(k * n), ncol = n))))
    }
    else{
      if(nrow(innovations)<n){stop("You need more innovations: Error 4")}
      epsilon <- rbind(a,innovations[1:n,])
    }
  }
  if (q > 0) {
    extend.psi <- array(c(diag(k), -theta, rep(0, k * k *
                                                 (n - q))), dim = c(k, k, n + 1))
    u <- matrix(numeric(0), nrow = n, ncol = k)
    for (i in (q + 1):(n + q)) {
      out <- 0
      for (j in 0:q) {
        out = out + crossprod(t(extend.psi[, , j + 1]),
                              epsilon[i - j, ])
      }
      u[i - q, ] <- out
    }
  }
  else u <- epsilon
  for (i in (r + 1):n) {
    temp2 <- 0
    for (j in 1:p) temp2 <- temp2 + crossprod(t(phi[, , j]),
                                              FirstSim.Series[i - j, ])
    FirstSim.Series[i, ] <- temp2 + u[i, ]
  }
  SlopDrift <- t(matrix(trend * rep(1:NROW(FirstSim.Series),
                                    each = k), nrow = k, ncol = NROW(FirstSim.Series)))
  DriftTerm <- sweep(SlopDrift, 2L, -constant, check.margin = FALSE)
  CenterData <- scale(FirstSim.Series, center = -demean, scale = FALSE)
  Sim.Series <- DriftTerm + CenterData
  for (i in 1:length(d)) {
    if (d[i] > 0)
      Sim.Series[, i] <- as.matrix(diffinv(Sim.Series[,
                                                      i], differences = d[i]))[-(1:d[i]), ]
    else Sim.Series[, i] <- Sim.Series[, i]
  }
  return(ts(Sim.Series))
}