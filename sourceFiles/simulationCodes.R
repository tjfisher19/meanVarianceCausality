####################
## This is the main call function that performs the simulation for the paper
##
## simEmpiricalSize() performs the simulations for the ARMA(1,1)+GARCH(1,1) simulations
## in the paper. If we specify alternative != 0, it is technically computing power
##
##
## The second function is the call function for the VAR(p)+(marginal) GARCH(1,1) simulations.
##
##
## In each function we specify the sample size n
## the maximum lag M
## The number of replications, run.times
## The number of bootstrap samples, boot.runs
## alternative, which specifies if an alternative hypothesis is true or not (0)
## one.way, whether a one-sided test should be performed for causality.
##
## Lastly, we have the function sampleSizeDetermination() which calculares the emprical size
## as a function of the sample size n
##
simEmpiricalSize <- function(cluster=NULL, n=250, M=c(5), run.times=1000, boot.runs=1000, alternative=0, one.way="no") {
  cat("************************************\n", file=error.file)
  cat("  Error Log for parameters\n\n", file=error.file, append=TRUE)
  cat("          n: ", n, "\n", file=error.file, append=TRUE)
  cat("          M: ", M, "\n", file=error.file, append=TRUE)
  cat("  run.times: ", run.times, "\n", file=error.file, append=TRUE)
  cat("  boot.runs: ", boot.runs, "\n", file=error.file, append=TRUE)
  cat("alternative: ", alternative, "\n", file=error.file, append=TRUE)
  cat("************************************\n\n", file=error.file, append=TRUE)
  N <- rep(n, run.times)
  out <- sapply(N, do.one.simulation, cluster=cluster, M=M, alternative=alternative, boot.run=boot.runs, one.way=one.way)
  compare <- function(x, sig=0.05) { mean(x<sig) }
  tab <- matrix(nrow=20, ncol=4*length(M), 0)
  if(length(M)==1) {
    tab[,1] <- apply(out, 1, compare, sig=0.10)
    tab[,2] <- apply(out, 1, compare, sig=0.05)
    tab[,3] <- apply(out, 1, compare, sig=0.01)
    tab[,4] <- apply(out, 1, compare, sig=0.001)
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig")
  } else {
    tmp <- apply(out, 1, compare, sig=0.10)
    tab[,1] <- tmp[1:20]
    tab[,5] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.05)
    tab[,2] <- tmp[1:20]
    tab[,6] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.01)
    tab[,3] <- tmp[1:20]
    tab[,7] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.001)
    tab[,4] <- tmp[1:20]
    tab[,8] <- tmp[21:40]
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig", "10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig")
  }
  
  
  rownames(tab) <- c("Boot BP Mean", "Boot BP Var",
                     "Boot LB Mean", "Boot LB Var",
                     "Boot WLB Mean", "Boot WLB Var",
                     "Boot Dan Mean", "Boot Dan Var",
                     "Boot Mat Mean", "Boot Mat Var",
                     "Theo BP Mean", "Theo BP Var",
                     "Theo LB Mean", "Theo LB Var",
                     "Theo WLB Mean", "Theo WLB Var",
                     "Theo Dan Mean", "Theo Dan Var",
                     "Theo Mat Mean", "Theo Mat Var")
  tab
}



simEmpiricalVARSize <- function(cluster=NULL, n=250, M=c(5), run.times=1000, boot.runs=1000, alternative=0, one.way="no") {
  cat("************************************\n", file=error.file)
  cat("  Error Log for parameters\n\n", file=error.file, append=TRUE)
  cat("          n: ", n, "\n", file=error.file, append=TRUE)
  cat("          M: ", M, "\n", file=error.file, append=TRUE)
  cat("  run.times: ", run.times, "\n", file=error.file, append=TRUE)
  cat("  boot.runs: ", boot.runs, "\n", file=error.file, append=TRUE)
  cat("alternative: ", alternative, "\n", file=error.file, append=TRUE)
  cat("************************************\n\n", file=error.file, append=TRUE)
  N <- rep(n, run.times)
  out <- sapply(N, do.one.var.simulation, cluster=cluster, M=M, alternative=alternative, boot.run=boot.runs, one.way=one.way)
  compare <- function(x, sig=0.05) { mean(x<sig) }
  tab <- matrix(nrow=20, ncol=4*length(M), 0)
  if(length(M)==1) {
    tab[,1] <- apply(out, 1, compare, sig=0.10)
    tab[,2] <- apply(out, 1, compare, sig=0.05)
    tab[,3] <- apply(out, 1, compare, sig=0.01)
    tab[,4] <- apply(out, 1, compare, sig=0.001)
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig")
  } else {
    tmp <- apply(out, 1, compare, sig=0.10)
    tab[,1] <- tmp[1:20]
    tab[,5] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.05)
    tab[,2] <- tmp[1:20]
    tab[,6] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.01)
    tab[,3] <- tmp[1:20]
    tab[,7] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.001)
    tab[,4] <- tmp[1:20]
    tab[,8] <- tmp[21:40]
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig", "10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig")
  }
  
  
  rownames(tab) <- c("Boot BP Mean", "Boot BP Var",
                     "Boot LB Mean", "Boot LB Var",
                     "Boot WLB Mean", "Boot WLB Var",
                     "Boot Dan Mean", "Boot Dan Var",
                     "Boot Mat Mean", "Boot Mat Var",
                     "Theo BP Mean", "Theo BP Var",
                     "Theo LB Mean", "Theo LB Var",
                     "Theo WLB Mean", "Theo WLB Var",
                     "Theo Dan Mean", "Theo Dan Var",
                     "Theo Mat Mean", "Theo Mat Var")
  tab
}




simEmpiricalTrivialSize <- function(cluster=NULL, n=250, M=c(5), run.times=1000, boot.runs=1000, one.way="no") {
  cat("************************************\n", file=error.file)
  cat("  Error Log for parameters\n\n", file=error.file, append=TRUE)
  cat("          n: ", n, "\n", file=error.file, append=TRUE)
  cat("          M: ", M, "\n", file=error.file, append=TRUE)
  cat("  run.times: ", run.times, "\n", file=error.file, append=TRUE)
  cat("  boot.runs: ", boot.runs, "\n", file=error.file, append=TRUE)
  cat("************************************\n\n", file=error.file, append=TRUE)
  N <- rep(n, run.times)
  out <- sapply(N, do.one.trivial.simulation, cluster=cluster, M=M, boot.run=boot.runs, one.way=one.way)
  compare <- function(x, sig=0.05) { mean(x<sig) }
  tab <- matrix(nrow=20, ncol=4*length(M), 0)
  if(length(M)==1) {
    tab[,1] <- apply(out, 1, compare, sig=0.10)
    tab[,2] <- apply(out, 1, compare, sig=0.05)
    tab[,3] <- apply(out, 1, compare, sig=0.01)
    tab[,4] <- apply(out, 1, compare, sig=0.001)
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig")
  } else {
    tmp <- apply(out, 1, compare, sig=0.10)
    tab[,1] <- tmp[1:20]
    tab[,5] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.05)
    tab[,2] <- tmp[1:20]
    tab[,6] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.01)
    tab[,3] <- tmp[1:20]
    tab[,7] <- tmp[21:40]
    tmp <- apply(out, 1, compare, sig=0.001)
    tab[,4] <- tmp[1:20]
    tab[,8] <- tmp[21:40]
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig", "10% Sig", " 5% Sig", " 1% Sig", "0.1% Sig")
  }
  
  
  rownames(tab) <- c("Boot BP Mean", "Boot BP Var",
                     "Boot LB Mean", "Boot LB Var",
                     "Boot WLB Mean", "Boot WLB Var",
                     "Boot Dan Mean", "Boot Dan Var",
                     "Boot Mat Mean", "Boot Mat Var",
                     "Theo BP Mean", "Theo BP Var",
                     "Theo LB Mean", "Theo LB Var",
                     "Theo WLB Mean", "Theo WLB Var",
                     "Theo Dan Mean", "Theo Dan Var",
                     "Theo Mat Mean", "Theo Mat Var")
  tab
}

sampleSizeDetermination <- function(cluster=NULL, n=seq(100,500,100), M=NULL, run.times=1000, one.way="no") {

    compare <- function(x, sig=0.05) { mean(x<sig) }
    tab1 <- matrix(nrow=10, ncol=length(n), 0)
    tab2 <- matrix(nrow=10, ncol=length(n), 0)
    tab3 <- matrix(nrow=10, ncol=length(n), 0)
    tab4 <- matrix(nrow=10, ncol=length(n), 0)
    n.names <- rep(0, length(n) )
    for(i in 1:length(n) ) {
      if(is.null(M) )
        M <- round(log(n[i]))
      N <- rep(n[i], run.times)
      
      if(is.null(cluster)) {
        out <- sapply(N, do.one.trivial.noboot.simulation, M=M, one.way=one.way)
      } else {
        out <- parSapply(cl=cluster, X=N, FUN=do.one.trivial.noboot.simulation,  
                                 M=M, one.way=one.way )
      }

      tab1[,i] <- apply(out, 1, compare, sig=0.10)
      tab2[,i] <- apply(out, 1, compare, sig=0.05)
      tab3[,i] <- apply(out, 1, compare, sig=0.01)
      tab4[,i] <- apply(out, 1, compare, sig=0.001)
      n.names[i] <- paste("n=",n[i], sep="")
    }
    
    
    rownames(tab4) <- c("Theo BP Mean", "Theo BP Var",
                       "Theo LB Mean", "Theo LB Var",
                       "Theo WLB Mean", "Theo WLB Var",
                       "Theo Dan Mean", "Theo Dan Var",
                       "Theo Mat Mean", "Theo Mat Var")
    rownames(tab1) <- rownames(tab2) <- rownames(tab3) <- rownames(tab4)
    colnames(tab4) <- n.names
    colnames(tab1) <- colnames(tab2) <- colnames(tab3) <- colnames(tab4)
    list(tab1, tab2, tab3, tab4)
}
