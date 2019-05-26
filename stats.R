rm(list=ls()) #cleans up working space
library(lawstat)

setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")

dat <- read.csv('data.full.csv', row.names=1)
######################################################################################################
#ANOVAs functional data T5
res.aov1 <- aov(BP ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov1)
TukeyHSD(res.aov1)
par(mfrow=c(2,2))
plot(res.aov1)
bartlett.test(dat$BP~interaction(dat$DOM,dat$source))
boxplot(BP ~ interaction(source,DOM), data = dat[dat$time=='T5',])  

res.aov2 <- aov(resp ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov2)
TukeyHSD(res.aov2)
par(mfrow=c(2,2))
plot(res.aov2)
bartlett.test(dat$resp~interaction(dat$DOM,dat$source))
boxplot(resp ~ interaction(source,DOM), data = dat[dat$time=='T5',])  


######################################################################################################
#ANOVAs networks
#SARA: I have inspected variance and distribution of the residuals using the plot function
#SARA: based on visual inspection I have chosen the transformation of the data
res.aov4 <- aov(log(Mean.Degree) ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov4)
TukeyHSD(res.aov4)
par(mfrow=c(2,2))
plot(res.aov4)
#test for homogeneity of variances (ok if p>0.05)
bartlett.test(log(dat$Mean.Degree)~interaction(dat$DOM,dat$source))
boxplot(log(Mean.Degree) ~ interaction(source,DOM), data = dat[dat$time=='T5',]) 

dat$prop.neg.logit<-(log(dat$prop.neg/(100-dat$prop.neg))) #logit transformation (http://strata.uga.edu/8370/rtips/proportions.html)
res.aov5 <- aov(prop.neg.logit ~ source*DOM, data = dat)
summary (res.aov5)
TukeyHSD(res.aov5)
plot(res.aov5)
bartlett.test(dat$prop.neg.logit~interaction(dat$DOM,dat$source))
boxplot(prop.neg.logit ~ interaction(source,DOM), data = dat[dat$time=='T5',]) 

res.aov5 <- aov(log(Diameter) ~ source*DOM, data = Results)
summary (res.aov6)
TukeyHSD(res.aov6)
plot(res.aov6)
bartlett.test(log(dat$Diameter)~interaction(dat$DOM,dat$source))
boxplot(log(Diameter) ~ interaction(source,DOM), data = dat[dat$time=='T5',]) 

res.aov7 <- aov(log(Density) ~ source*DOM, data = Results)
summary (res.aov7)
TukeyHSD(res.aov7)
plot(res.aov7)
bartlett.test(log(dat$Density)~interaction(dat$DOM,dat$source))
boxplot(log(Density) ~ interaction(source,DOM), data = dat[dat$time=='T5',]) 

res.aov8 <- aov(richness.bins.netw ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov8)
TukeyHSD(res.aov8)
plot(res.aov8)
bartlett.test(dat$richness.bins.netw~interaction(dat$DOM,dat$source))
boxplot(Density ~ interaction(source,DOM), data = dat[dat$time=='T5',]) 

hist(dat[dat$time=='T5',16])

#no significance
res.aov4 <- aov(Modularity ~ source*DOM, data = Results)
summary (res.aov4)
TukeyHSD(res.aov4)
plot(res.aov4)
bartlett.test(log(dat$Modularity)~interaction(dat$DOM,dat$source))
boxplot(log(Modularity) ~ interaction(source,DOM), data = dat[dat$time=='T5',]) 
