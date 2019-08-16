##################################################
## 
## This code generates the empirical size as a function
## of the sample size for Q_LB
##
## Fairly straightforward ggplot2 code here.
##

library(ggplot2)
library(dplyr)
library(tidyr)
tmp1 <- read.csv("resultCsvFiles/sizeDeterminationAlpha05MgrowsLogN.csv")
tmp2 <- read.csv("resultCsvFiles/sizeDeterminationAlpha05Mfixed5.csv")

my.data <- data.frame(n=seq(25,1500,1), LB1=as.numeric(tmp1[4,2:1477]), LB2=as.numeric(tmp2[4,2:1477]) )

names(tmp) <- c("Statistic", paste(seq(50,1500,10)))

# sizes.long <- tmp[c(4,6,8),] %>% gather(key, value, -Statistic)
sizes.long <- my.data %>% gather(n)
names(sizes.long) <- c("n", "Statistic", "value") 
sizes.long$n <- as.numeric(sizes.long$n)
sizes.long$lower <- (0.05 - 1.96*sqrt(0.05*0.95/10000) )*100
sizes.long$upper <- (0.05 + 1.96*sqrt(0.05*0.95/10000) )*100
sizes.long$value <- sizes.long$value*100
sizes.long$Tom <- as.factor(c("A", "B", "C"))


plot.labels = c(expression(italic("Q")[LB]^"\u2217"*"("*italic(m)*"=5)" ), 
                expression(italic("Q")[LB]^"\u2217"~"("*italic(m)*"=log("*italic(n)* "))" ) )


p.sizes <- ggplot(sizes.long) +
  geom_ribbon(aes(x=n, ymin=lower, ymax=upper), alpha=0.15) +
  geom_line(aes(x=n, y=value, linetype=Statistic, col=Statistic), alpha=0.35 ) +
  geom_smooth(aes(x=n, y=value, linetype=Statistic), col="black", se=FALSE, size=1.25) +
  scale_linetype_manual("Legend Title\n",labels=plot.labels, values = c(1, 6)) +
  scale_color_manual("Legend Title\n",labels=plot.labels, values = c("gray60", "gray30")) +
  labs(x=expression(paste("Sample size ", italic(n))), y="Empirical Size") +
  annotate("segment", x=50, y=5, xend=1500, yend=5, size=0.5) +
  theme_bw() +
  theme(axis.line=element_line(size=1, colour="black"),
        axis.text = element_text(colour="gray20"),
        legend.position = c(0.75, 0.85), 
        legend.text.align=0,
        legend.background = element_rect(colour = "black"),
        legend.title=element_blank(),
        legend.text=element_text(size=18, family="serif"),
        text = element_text(size=15, family="serif") );
p.sizes

cairo_ps(filename="sizeDeterminationLjungBox.eps", width=8, height=6)
print(p.sizes)
dev.off()


