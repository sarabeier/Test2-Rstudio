rm(list=ls()) #cleans up working space
library(lawstat)

setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")

dat <- read.csv('data.full.csv', row.names=1)
######################################################################################################
#ANOVAs functional data T5
res.aov1 <- aov((BP) ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov1)
TukeyHSD(res.aov1)
par(mfrow=c(2,2))
plot(res.aov1)
bartlett.test(dat$BP~interaction(dat$DOM,dat$source))

res.aov2 <- aov((resp) ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov2)
TukeyHSD(res.aov2)
par(mfrow=c(2,2))
plot(res.aov2)
bartlett.test(dat$resp~interaction(dat$DOM,dat$source))


######################################################################################################
#ANOVAs networks
#SARA: I have inspected variance and distribution of the residuals using the plot function
#SARA: based on visual inspection I have chosen the transformation of the data
res.aov1 <- aov(log(Mean.Degree) ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov1)
TukeyHSD(res.aov1)
par(mfrow=c(2,2))
plot(res.aov1)
#test for homogeneity of variances (ok if p>0.05)
bartlett.test(log(dat$Mean.Degree)~interaction(dat$DOM,dat$source))

dat$prop.neg.logit<-(log(dat$prop.neg/(100-dat$prop.neg))) #logit transformation (http://strata.uga.edu/8370/rtips/proportions.html)
res.aov2 <- aov(prop.neg.logit ~ source*DOM, data = dat)
summary (res.aov2)
TukeyHSD(res.aov2)
plot(res.aov2)
bartlett.test(dat$prop.neg.logit~interaction(dat$DOM,dat$source))

res.aov5 <- aov(log(Diameter) ~ source*DOM, data = Results)
summary (res.aov5)
TukeyHSD(res.aov5)
plot(res.aov5)
bartlett.test(log(dat$Diameter)~interaction(dat$DOM,dat$source))

res.aov3 <- aov(log(Density) ~ source*DOM, data = Results)
summary (res.aov3)
TukeyHSD(res.aov3)
plot(res.aov3)
bartlett.test(log(dat$Density)~interaction(dat$DOM,dat$source))

#no significance
res.aov4 <- aov(Modularity ~ source*DOM, data = Results)
summary (res.aov4)
TukeyHSD(res.aov4)
plot(res.aov4)
bartlett.test(log(dat$Modularity)~interaction(dat$DOM,dat$source))
