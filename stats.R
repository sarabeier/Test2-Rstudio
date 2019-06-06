rm(list=ls()) #cleans up working space
library(lawstat)
library(car)

setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")

dat <- read.csv('data.full.csv', row.names=1) #load data
dat$prop.neg.logit<-(log(dat$prop.neg/(100-dat$prop.neg))) #logit transformation for precent data (http://strata.uga.edu/8370/rtips/proportions.html)
head(dat)
######################################################################################################
#ANOVAs functional data T5
res.aov1 <- aov(BP ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov1)
TukeyHSD(res.aov1)
par(mfrow=c(2,2))
plot(res.aov1)
#tests for homogeneity of variances (ok if p>0.05)
bartlett.test(dat$BP~interaction(dat$DOM,dat$source)) #bartlett test sensitive to deviation from normal distribution
leveneTest(BP ~ source*DOM, data = dat[dat$time=='T5',]) #levene test not sensitive to deviation from normal distribution

res.aov2 <- aov(resp ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov2)
TukeyHSD(res.aov2)
par(mfrow=c(2,2))
plot(res.aov2)
leveneTest(resp ~ source*DOM, data = dat[dat$time=='T5',])

res.aov3 <- aov(resp ~ source*DOM, data = dat[dat$time=='T3',])
summary (res.aov3)
TukeyHSD(res.aov3)
par(mfrow=c(2,2))
plot(res.aov3)
leveneTest(resp ~ source*DOM, data = dat[dat$time=='T3',])

res.aov4 <- aov(cellcounts ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov4)
TukeyHSD(res.aov4)
par(mfrow=c(2,2))
plot(res.aov4)
leveneTest(cellcounts ~ source*DOM, data = dat[dat$time=='T5',])


######################################################################################################
#ANOVAs networks
#SARA: I have inspected variance and distribution of the residuals using the plot function
#SARA: based on visual inspection I have chosen the transformation of the data
res.aov5 <- aov(log(Mean.Degree) ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov5)
TukeyHSD(res.aov5)
par(mfrow=c(2,2))
plot(res.aov5)
leveneTest(log(Mean.Degree) ~ source*DOM, data = dat[dat$time=='T5',])

res.aov6 <- aov(prop.neg.logit ~ source*DOM, data = dat)
summary (res.aov6)
TukeyHSD(res.aov6)
plot(res.aov6)
leveneTest(prop.neg.logit ~ source*DOM, data = dat[dat$time=='T5',])

res.aov7 <- aov(log(Diameter) ~ source*DOM, data = dat)
summary (res.aov7)
TukeyHSD(res.aov7)
plot(res.aov7)
leveneTest(log(Diameter) ~ source*DOM, data = dat[dat$time=='T5',])

res.aov8 <- aov(log(Density) ~ source*DOM, data = dat)
summary (res.aov8)
TukeyHSD(res.aov8)
plot(res.aov8)
leveneTest(log(Density) ~ source*DOM, data = dat[dat$time=='T5',])

res.aov9 <- aov(richness.bins.netw ~ source*DOM, data = dat[dat$time=='T5',])
summary (res.aov9)
TukeyHSD(res.aov9)
plot(res.aov9)
leveneTest(richness.bins.netw ~ source*DOM, data = dat[dat$time=='T5',])

#no significance
res.aov10 <- aov(Modularity ~ source*DOM, data = dat)
summary (res.aov10)
TukeyHSD(res.aov10)
plot(res.aov10)
leveneTest(log(Modularity) ~ source*DOM, data = dat[dat$time=='T5',])

######################################################################################################
#ANOVAs diversity
res.aov11 <- aov(log((H.ASV)) ~ source*DOM, data = dat)
summary (res.aov11)
TukeyHSD(res.aov11)
plot(res.aov11)
leveneTest(log(H.ASV) ~ source*DOM, data = dat[dat$time=='T5',])

res.aov12 <- aov(((S.ASV)) ~ source*DOM, data = dat)
summary (res.aov12)
TukeyHSD(res.aov12)
plot(res.aov12)
leveneTest((S.ASV) ~ source*DOM, data = dat[dat$time=='T5',])

par(mfrow=c(4,3))

boxplot(BP ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='BP')  
boxplot(resp ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='res')
boxplot(cellcounts ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='cellcounts')
boxplot(log(Mean.Degree) ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='mean.degree') 
boxplot(prop.neg.logit ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='prop.neg') 
boxplot(Density ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='Density')
boxplot(log(Modularity) ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='Modularity') 
boxplot(richness.bins.netw ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='rich.bin') 
boxplot(S.ASV ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='rich.ASV') 
boxplot(H.ASV ~ interaction(source,DOM), data = dat[dat$time=='T5',], main='H.ASV') 

