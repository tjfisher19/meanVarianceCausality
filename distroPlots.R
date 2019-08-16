##################################################
## 
## This code generates the empirical density plots
##
## Fairly straightforward ggplot2 code here.
##
source("sourceFiles/testStatFunctions.R")
source("sourceFiles/bootstrapSimFunctions.R")
source("sourceFiles/simulationCodes.R")
source("sourceFiles/nonMatrixTest.R")
source("sourceFiles/matrixBasedTests.R")
source("sourceFiles/trivialSimFunctions.R")

quickPlots <- function(n=30, m=6, runs=1000, one.way="no", squares=FALSE) {
  get.stats <- function(n=30,m=6, one.way=one.way) {
    tmp <- generateNormalData(n=n)
    x1 <- tmp$x1
    x2 <- tmp$x2
    
    x1.res <- scale(x1)[,1]
    x2.res <- scale(x2)[,1]
    if(squares) {
      tmp <- get.test.stats(r1=x1.res^2, r2=x2.res^2, M=m, one.way=one.way, matrix=FALSE)
    } else {
      tmp <- get.test.stats(r1=x1.res, r2=x2.res, M=m, one.way=one.way, matrix=FALSE)
    }
    stats <- tmp[,1]
    stats
  }
  
  N <- rep(n, runs)
  stats <- sapply(N, get.stats, m=m, one.way=one.way)
  stats
  
  mean.test <- data.frame(type=c(rep("BP",runs),rep("LB",runs)), value=c(stats[1,],stats[3,]))
  
  if(squares) {
    plot.labels = c(expression(paste(italic(Q)[italic("BP")]^"\u204e")), 
                    expression(italic("Q")[italic("LB")]^"\u204e" ) )   
  } else {
    plot.labels = c(expression(paste(italic(Q)[italic("BP")])), 
                    expression(paste(italic(Q)[italic("LB")]) ) )
  }
  
  p <- ggplot(mean.test) +
    geom_density(aes(value, fill=type),alpha=0.88, size=0) +      ## empircal density
    stat_function(fun=dchisq, colour="black",                     ## chi square distro
                  args=list(df=(13)), linetype=2, size=1.25 )+
    geom_vline(xintercept=qchisq(0.90, df=13)) +                  ## 90% chi sq cutoff
    annotate("text", x=20.7, y=0.045, angle=75, label="90% cutoff", size=5) +
    geom_vline(xintercept=qchisq(0.95, df=13)) +                  ## 95% chi sq cutoff
    annotate("text", x=23.2, y=0.045, angle=75, label="95% cutoff", size=5) +
    geom_vline(xintercept=qchisq(0.99, df=13)) +                  ## 99% chi sq cutoff
    annotate("text", x=28.5, y=0.045, angle=75, label="99% cutoff", size=5) +
    
    annotate("text", x=25, y=0.005, angle=345,                    ## "Asymtotic chi^2"
             label="Asymptotic~chi[13]^2~distribution", parse=TRUE, size=5, color="white") +
    #annotate("text", x=25, y=0.0055, angle=337,                  ## "distribution"
    #         label="distribution", parse=TRUE, size=5, color="white") +
    
    coord_cartesian(xlim=c(18.85,32), ylim=c(0,0.05)) +           ## Limit coordinates
    labs(x=expression(paste("Statistic Value")), y="Density") +   ##  Axes
    scale_fill_manual("Legend Title\n",labels=plot.labels,        ##  colors of density
                      values = c("gray20", "gray50")) +
    theme_bw() +
    theme(axis.line=element_line(size=1, colour="black"),
          axis.text = element_text(colour="gray20"),
          legend.position = c(0.85, 0.5), 
          legend.text.align=0,
          legend.background = element_rect(colour = "black"),
          legend.title=element_blank(),
          legend.text=element_text(size=18, family="serif"),
          text = element_text(size=15, family="serif") );
  p
}

set.seed(12345)
p.mean <- quickPlots(runs=10000, n=30, m=6)
p.var <- quickPlots(runs=10000, n=30, m=6, squares = TRUE)

p.mean
p.var

cairo_ps(filename="simMeanStatsN30m6.eps", width=8, height=6)
print(p.mean)
dev.off()

cairo_ps(filename="simVarStatsN30m6.eps", width=8, height=6)
print(p.var)
dev.off()
