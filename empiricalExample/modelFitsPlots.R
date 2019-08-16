######################################
##
## This code performs the analysis of the S&P500 and EUR data
##
## We also generate the plots for the paper here.
##

library(fGarch)
library(forecast)
library(lubridate)
library(dplyr)
library(tidyr)
library(ggplot2)

#### Load data and process
data.sp500 <- read.csv("YAHOO-INDEX_GSPC.csv")
head(data.sp500)
data.sp500<-data.sp500[,c(1,7)]
data.sp500$Date<- ymd(data.sp500$Date)
colnames(data.sp500)<-c("Date","sp500")

data.EUR.2015 <- read.csv("EUR2015.csv")
data.EUR.2015$Date<- ymd(data.EUR.2015$Date)

### Combine
data.full<- right_join(data.sp500, data.EUR.2015)

data.full
head(data.full)
colnames(data.full)<-c("Date","S&P500 Index","Exchange Rate (EUR/USD)")
data.full$Closed <- as.factor(is.na(data.full$'S&P500 Index')*1)

rtn.obj<- gather(data.full,"variable","value",2:3)
rtn.obj$Date <- ymd(rtn.obj$Date)
head(rtn.obj)
str(rtn.obj)

gg.xts1 <- ggplot(rtn.obj, aes_string(x = "Date", y = "value", group = "variable", color="Closed" )) +
  facet_grid(variable~. , scales = "free_y", space = "fixed") +
  geom_line(data = subset(rtn.obj, variable == "S&P500 Index"),size=0.8) +
  geom_line(data = subset(rtn.obj, variable == "Exchange Rate (EUR/USD)"),size=0.8) +
  scale_colour_manual(values=c("black", "gray70")) + 
  ggtitle("") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        plot.title = element_text(size=20, face="bold",vjust=1),
        strip.background = element_rect(fill="white"),
        strip.text.y = element_text(size = 16, colour="black", family="serif",angle=270),
        legend.position="none") + 
  ylab("") +
  xlab("Date")
gg.xts1

cairo_ps(filename="sp500EuroDollars2015.eps", width=8, height=6)
print(gg.xts1)
dev.off()

### REMOVE MISSING VALUES
colnames(data.full)<-c("Date","sp500","EUR","Closed")
head(data.full)
data.full <- na.omit(data.full)
head(data.full)

## LOG DIFFERENCING
sp500<- diff(log(data.full$sp500))
EUR<- diff(log(data.full$EUR))


## Explore features about sp500
acf(sp500, lwd=2)
pacf(sp500, lwd=2)   ## Maybe an AR(4) or smaller ARMA
auto.arima(sp500)  ### Says ARMA(1,2)
arma.fit.sp500 <- Arima(sp500, order=c(1,0,2), include.mean=FALSE)

Weighted.Box.test(arma.fit.sp500$residuals, lag=30, fitdf=3)

acf(sp500^2)
pacf(sp500^2)
### GARCH(1,1)???
Weighted.Box.test(arma.fit.sp500$residuals, lag=30, fitdf=3, sqrd.res = TRUE)


sp500.fit <- garchFit(~arma(1,2)+garch(1,1), data=sp500, trace=FALSE, include.mean=FALSE, cond.dist="QMLE")
sp500.res <- sp500.fit@residuals/sqrt(sp500.fit@h.t) 
acf(sp500.res)
acf(sp500.res^2)
#####
## Marginall fit with an ARMA(1,2)+GARCH(1,1)



########################
# EUR
acf(EUR)
pacf(EUR)
auto.arima(EUR)  ## Says ARMA(1,3)
Arima(EUR, order=c(1,0,0), include.mean=FALSE)
Arima(EUR, order=c(1,0,3), include.mean=FALSE)
Arima(EUR, order=c(1,0,1), include.mean=FALSE)
Arima(EUR, order=c(0,0,1), include.mean=FALSE)  ## Good BIC, close to matching AIC...

arma.fit.EUR <- Arima(EUR, order=c(0,0,1), include.mean=FALSE)

Weighted.Box.test(arma.fit.EUR$residuals, lag=30, fitdf=1)
acf(arma.fit.EUR$residuals)
pacf(arma.fit.EUR$residual)

acf(arma.fit.EUR$residuals^2)
pacf(arma.fit.EUR$residuals^2)
acf(abs(arma.fit.EUR$residuals))  ### Nothing going on...
Weighted.Box.test(arma.fit.EUR$residuals, lag=30, fitdf=1, sqrd.res = TRUE)

#### No GARCH here...
EUR.res <- arma.fit.EUR$residuals
EUR.fit <- arma.fit.EUR
acf(EUR.res)
acf(EUR.res^2)
EUR.res <- scale(EUR.res)[,1]

########## Just an MA(1)...


mean.ccf <- ccf(sp500.res,EUR.res, main="", lwd=3, ci.col="gray20", lag.max=20 )
mean.ccfdf <- with(mean.ccf, data.frame(lag, acf, lo=qnorm(0.025, sd=1/sqrt(n.used)), up=qnorm(0.975, sd=1/sqrt(n.used))))

p.ccf.mean <- ggplot(data=mean.ccfdf, mapping=aes(x=lag, y=acf)) + 
        geom_hline(aes(yintercept=0)) +
        geom_segment(mapping=aes(xend=lag, yend=0), size=1.75) +
        geom_ribbon(aes(ymin=lo, ymax=up), alpha=0.2) +
        coord_cartesian(xlim=c(-16,16), ylim=c(-0.25,0.25)) +
        labs(x="Lag", y="Cross-correlation in mean") + 
        theme_bw() + 
        theme(axis.line=element_line(size=1, colour="black"),
              axis.text = element_text(colour="gray20"),
              text = element_text(size=15, family="serif") );
p.ccf.mean
cairo_ps(filename="sp500EuroResidCCFmean.eps", width=8, height=6)
print(p.ccf.mean)
dev.off()

var.ccf <- ccf(sp500.res^2,EUR.res^2, main="", lwd=3, ci.col="gray20", lag.max=20)
var.ccfdf <- with(var.ccf, data.frame(lag, acf, lo=qnorm(0.025, sd=1/sqrt(n.used)), up=qnorm(0.975, sd=1/sqrt(n.used))))
p.ccf.var <- ggplot(data=var.ccfdf, mapping=aes(x=lag, y=acf)) + 
  geom_hline(aes(yintercept=0)) +
  geom_segment(mapping=aes(xend=lag, yend=0), size=1.75) +
  geom_ribbon(aes(ymin=lo, ymax=up), alpha=0.2) +
  coord_cartesian(xlim=c(-16,16), ylim=c(-0.25,0.25)) +
  labs(x="Lag", y="Cross-correlation in variance") + 
  theme_bw() + 
  theme(axis.line=element_line(size=1, colour="black"),
        axis.text = element_text(colour="gray20"),
        text = element_text(size=15, family="serif") );
p.ccf.var
cairo_ps(filename="sp500EuroResidCCFvar.eps", width=8, height=6)
print(p.ccf.var)
dev.off()

sgedFit(sp500.res)$par
### > sgedFit(sp500.res)$par
### mean          sd          nu          xi 
### -0.01490927  0.99731078  1.58947645  0.91475410 
sp500.res.data <- data.frame(resid=sp500.res)
p.sp500.res <- ggplot(sp500.res.data, aes(resid)) +
  geom_density(size=0, fill="gray85", adjust=1) +      ## empircal density
  stat_function(fun=dnorm, colour="gray30",  
                linetype=3, size=1.25 )+
  stat_function(fun=dsged, colour="black", linetype=2, size=1.25,
                args=list(mean=-0.01490927, sd=0.99731078, nu=1.58947645, xi=0.91475410)) +
  scale_x_continuous(limits=c(-4,4)) + 
  labs(x="S&P 500 Redisuals", y="Density") +
  theme_bw() +
  theme(axis.line=element_line(size=1, colour="black"),
        axis.text = element_text(colour="gray20"),
        text = element_text(size=15, family="serif") );
p.sp500.res
cairo_ps(filename="sp500residDensity.eps", width=8, height=6)
print(p.sp500.res)
dev.off()


sgedFit(EUR.res)$par
### > sgedFit(EUR.res)$par
### mean           sd           nu           xi 
### -0.003652755  0.996272884  1.402415848  0.949465498 
eur.res.data <- data.frame(resid=EUR.res)
p.eur.res <- ggplot(eur.res.data, aes(resid)) +
  geom_density(size=0, fill="gray85", adjust=1) +      ## empircal density
  stat_function(fun=dnorm, colour="gray30",  
                linetype=3, size=1.25 )+
  stat_function(fun=dsged, colour="black", linetype=2, size=1.25,
                args=list(mean=-0.003652755, sd=0.996272884, nu=1.402415848, xi=0.949465498) ) +
  scale_x_continuous(limits=c(-4,4)) + 
  labs(x="EUR/USD Redisuals", y="Density") +
  theme_bw() +
  theme(axis.line=element_line(size=1, colour="black"),
        axis.text = element_text(colour="gray20"),
        text = element_text(size=15, family="serif") );
p.eur.res
cairo_ps(filename="EURUSDresidDensity.eps", width=8, height=6)
print(p.eur.res)
dev.off()


X <- cbind(sp500, EUR)
var.fit <- ar(X)
var.fit
resids <- var.fit$resid[-(1:4),]
acf(resids[,1])    ### Good
acf(resids[,1]^2)  
pacf(resids[,1]^2) ### Still a GARCH(1,1)??

acf(resids[,2])     ## Good
acf(resids[,2]^2)
pacf(resids[,2]^2)  ## Good

var.sp500.garch.fit <- garchFit(~garch(1,1), data=resids[,1], include.mean=FALSE, trace=FALSE, cond.dist="QMLE")
coef(var.sp500.garch.fit)
var.sp500.garch.res <- var.sp500.garch.fit@residuals/sqrt(var.sp500.garch.fit@h.t)
var(resids[,2])
var.eur.res <- scale(resids[,2])[,1]

ccf(var.sp500.garch.res, var.eur.res, main="", lwd=3, ci.col="gray20", lag.max = 16)
ccf(var.sp500.garch.res^2, var.eur.res^2, main="", lwd=3, ci.col="gray20", lag.max = 16)



var.mean.ccf <- ccf(var.sp500.garch.res,var.eur.res, main="", lwd=3, ci.col="gray20", lag.max=20 )
var.mean.ccfdf <- with(var.mean.ccf, data.frame(lag, acf, lo=qnorm(0.025, sd=1/sqrt(n.used)), up=qnorm(0.975, sd=1/sqrt(n.used))))

p.var.ccf.mean <- ggplot(data=var.mean.ccfdf, mapping=aes(x=lag, y=acf)) + 
  geom_hline(aes(yintercept=0)) +
  geom_segment(mapping=aes(xend=lag, yend=0), size=1.75) +
  geom_ribbon(aes(ymin=lo, ymax=up), alpha=0.2) +
  coord_cartesian(xlim=c(-16,16), ylim=c(-0.25,0.25)) +
  labs(x="Lag", y="Cross-correlation in mean") + 
  theme_bw() + 
  theme(axis.line=element_line(size=1, colour="black"),
        axis.text = element_text(colour="gray20"),
        text = element_text(size=15, family="serif") );
p.var.ccf.mean
cairo_ps(filename="sp500EuroVARResidCCFmean.eps", width=8, height=6)
print(p.var.ccf.mean)
dev.off()

var.var.ccf <- ccf(var.sp500.garch.res^2, var.eur.res^2, main="", lwd=3, ci.col="gray20", lag.max=20)
var.var.ccfdf <- with(var.var.ccf, data.frame(lag, acf, lo=qnorm(0.025, sd=1/sqrt(n.used)), up=qnorm(0.975, sd=1/sqrt(n.used))))
p.var.ccf.var <- ggplot(data=var.var.ccfdf, mapping=aes(x=lag, y=acf)) + 
  geom_hline(aes(yintercept=0)) +
  geom_segment(mapping=aes(xend=lag, yend=0), size=1.75) +
  geom_ribbon(aes(ymin=lo, ymax=up), alpha=0.2) +
  coord_cartesian(xlim=c(-16,16), ylim=c(-0.25,0.25)) +
  labs(x="Lag", y="Cross-correlation in variance") + 
  theme_bw() + 
  theme(axis.line=element_line(size=1, colour="black"),
        axis.text = element_text(colour="gray20"),
        text = element_text(size=15, family="serif") );
p.var.ccf.var
cairo_ps(filename="sp500EuroVARResidCCFvar.eps", width=8, height=6)
print(p.var.ccf.var)
dev.off()

