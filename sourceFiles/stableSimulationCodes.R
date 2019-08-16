####################
## A modification of simulationCodes.R to perform the robustness studies.
##
##
simStableEmpiricalSize <- function(cluster=NULL, n=250, M=c(5), run.times=1000, boot.runs=1000, one.way="no") {
  cat("************************************\n", file=error.file)
  cat("  Error Log for parameters\n\n", file=error.file, append=TRUE)
  cat("          n: ", n, "\n", file=error.file, append=TRUE)
  cat("          M: ", M, "\n", file=error.file, append=TRUE)
  cat("  run.times: ", run.times, "\n", file=error.file, append=TRUE)
  cat("  boot.runs: ", boot.runs, "\n", file=error.file, append=TRUE)
  cat("************************************\n\n", file=error.file, append=TRUE)
  N <- rep(n, run.times)
  out <- sapply(N, do.one.var.stable.null.simulation, cluster=cluster, M=M, boot.run=boot.runs, one.way=one.way)
  compare <- function(x, sig=0.05) { mean(x<sig) }
  tab <- matrix(nrow=20, ncol=3*length(M), 0)
  if(length(M)==1) {
    tab[,1] <- apply(out, 1, compare, sig=0.10)
    tab[,2] <- apply(out, 1, compare, sig=0.05)
    tab[,3] <- apply(out, 1, compare, sig=0.01)
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig")
  } else {
    tmp <- apply(out, 1, compare, sig=0.10)
    tab[,1] <- tmp[1:20]
    tab[,4] <- tmp[21:40]
    tab[,7] <- tmp[41:60]
    tmp <- apply(out, 1, compare, sig=0.05)
    tab[,2] <- tmp[1:20]
    tab[,5] <- tmp[21:40]
    tab[,8] <- tmp[41:60]
    tmp <- apply(out, 1, compare, sig=0.01)
    tab[,3] <- tmp[1:20]
    tab[,6] <- tmp[21:40]
    tab[,9] <- tmp[41:60]
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig", "10% Sig", " 5% Sig", " 1% Sig", "10% Sig", " 5% Sig", " 1% Sig")
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



simStableVarEmpiricalSize <- function(cluster=NULL, n=250, M=c(5), run.times=1000, boot.runs=1000, one.way="no") {
  cat("************************************\n", file=error.file)
  cat("  Error Log for parameters\n\n", file=error.file, append=TRUE)
  cat("          n: ", n, "\n", file=error.file, append=TRUE)
  cat("          M: ", M, "\n", file=error.file, append=TRUE)
  cat("  run.times: ", run.times, "\n", file=error.file, append=TRUE)
  cat("  boot.runs: ", boot.runs, "\n", file=error.file, append=TRUE)
  cat("************************************\n\n", file=error.file, append=TRUE)
  N <- rep(n, run.times)
  out <- sapply(N, do.one.var.stable.var.simulation, cluster=cluster, M=M, boot.run=boot.runs, one.way=one.way)
  compare <- function(x, sig=0.05) { mean(x<sig) }
  tab <- matrix(nrow=20, ncol=3*length(M), 0)
  if(length(M)==1) {
    tab[,1] <- apply(out, 1, compare, sig=0.10)
    tab[,2] <- apply(out, 1, compare, sig=0.05)
    tab[,3] <- apply(out, 1, compare, sig=0.01)
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig")
  } else {
    tmp <- apply(out, 1, compare, sig=0.10)
    tab[,1] <- tmp[1:20]
    tab[,4] <- tmp[21:40]
    tab[,7] <- tmp[41:60]
    tmp <- apply(out, 1, compare, sig=0.05)
    tab[,2] <- tmp[1:20]
    tab[,5] <- tmp[21:40]
    tab[,8] <- tmp[41:60]
    tmp <- apply(out, 1, compare, sig=0.01)
    tab[,3] <- tmp[1:20]
    tab[,6] <- tmp[21:40]
    tab[,9] <- tmp[41:60]
    colnames(tab) <- c("10% Sig", " 5% Sig", " 1% Sig", "10% Sig", " 5% Sig", " 1% Sig", "10% Sig", " 5% Sig", " 1% Sig")
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
