rm(list=ls()) #cleans up working space
library(lawstat)
library(car)

setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")

dat <- read.csv('data.full.csv', row.names=1) #load data
dat$prop.neg.logit<-(log(dat$prop.neg/(100-dat$prop.neg))) #logit transformation for precent data (http://strata.uga.edu/8370/rtips/proportions.html)
head(dat)
traits <-read.table('picrust.traits.tab', header=T, row.names=1) #datafile with picrust traits
head(traits)
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

######################################################################################################
#one sample t.tests
#https://stats.stackexchange.com/questions/51242/statistical-difference-from-zero
C.traits <- traits[traits$source=='C',] #select for samples from Canet-community
C.traits$av.16sdiff <- C.traits$av.16s-traits[1,4] #substract 16s copy number for treatments from values of original inoculum
C.traits$av.genomesizediff <- C.traits$av.genomesize-traits[1,5] #substract genome size for treatments from values of original inoculum

S.traits <- traits[traits$source=='S',] #select for samples from SOLA-community
S.traits$av.16sdiff <- S.traits$av.16s-traits[1,4] #substract 16s copy number for treatments from values of original inoculum
S.traits$av.genomesizediff <- S.traits$av.genomesize-S.traits[7,5] #substract genome size for treatments from values of original inoculum

traits.diff <- rbind(C.traits[-1,],S.traits[-7,]) #join substracted values from the Canet and SOLA community into one dataframe
traits.diff #dataframe displaying differences of trait-values between inoculum and treatments

#function to extract p-value from one sample one-tailed p-test (hypothesis: 16s copies, genomesize increase compared to org)
ttestp <-function(values) {
  pvalue<-t.test(values, alternative='greater')[3] #this line codes for the one-tailed p-test (for two-tailed t-test replace 'greater' by 'two.sided')
  return(pvalue)
}


pvalues <-aggregate(. ~ source + DOM, data=traits.diff[-1,-c(3:5)], FUN=ttestp) #returns table with p-values
pvalues[with(pvalues, order(source, DOM)), ]
