#repeated measurement ANOVA
######################################################################################
#example
#https://www.r-bloggers.com/two-way-anova-with-repeated-measures/
set.seed(5250)

myData <- data.frame(PID = rep(seq(from = 1,
                                   to = 50, by = 1), 20),
                     stress = sample(x = 1:100,
                                     size = 1000,
                                     replace = TRUE),
                     image = sample(c("Happy", "Angry"),
                                    size = 1000,
                                    replace = TRUE),
                     music = sample(c("Disney", "Horror"),
                                    size = 1000,
                                    replace = TRUE)
)

myData <- within(myData, {
  PID   <- factor(PID)
  image <- factor(image)
  music <- factor(music)
})

myData <- myData[order(myData$PID), ]
head(myData)

#Extracting Condition Means
myData.mean <- aggregate(myData$stress,
                         by = list(myData$PID, myData$music,
                                   myData$image),
                         FUN = 'mean')

colnames(myData.mean) <- c("PID","music","image","stress")

myData.mean <- myData.mean[order(myData.mean$PID), ]
head(myData.mean)

#Building the ANOVA
stress.aov <- with(myData.mean,
                   aov(stress ~ music * image +
                         Error(PID / (music * image)))
)
summary(stress.aov)

######################################################################################
#rmANOVA with own data
#time = PID
#cellcounts = stress
dat <- read.csv('data.full.csv', row.names=1)
#####same format as above


#rm.aov <- with(dat,
#                aov(cellcounts ~ DOM*source +
#                      Error(time / (DOM*source)))
#)
#summary(rm.aov)

rm.aov <- aov(cellcounts ~ DOM*source +Error(time / (DOM*source)), data=dat) #rmANOVA output
summary(rm.aov)
TukeyHSD(rm.aov) #does not work!!

#individual post-hoc pairwise tests (p-values still need to be corrected manually for multiple testing)
rm.aov.p1 <- aov(cellcounts ~ source +Error(time / (source)), data=dat[dat$source!='S',]) #rmANOVA output
summary(rm.aov.p1)
rm.aov.p2 <- aov(cellcounts ~ source +Error(time / (source)), data=dat[dat$source!='C',]) #rmANOVA output
summary(rm.aov.p2)
rm.aov.p3 <- aov(cellcounts ~ source +Error(time / (source)), data=dat[dat$source!='CS',]) #rmANOVA output
summary(rm.aov.p3)

rm.aov.p4 <- aov(cellcounts ~ DOM +Error(time / (DOM)), data=dat[dat$source=='S',]) #rmANOVA output
summary(rm.aov.p4)
rm.aov.p5 <- aov(cellcounts ~ DOM +Error(time / (DOM)), data=dat[dat$source=='C',]) #rmANOVA output
summary(rm.aov.p5)
rm.aov.p6 <- aov(cellcounts ~ DOM +Error(time / (DOM)), data=dat[dat$source=='CS',]) #rmANOVA output
summary(rm.aov.p6)


#simple ANOVA at time T5:
aovT5 <- aov(cellcounts ~ DOM*source, data=dat[dat$time=='T5',])
summary(aovT5)
TukeyHSD(aovT5) 




summary(aov(cellcounts ~ DOM, data=dat[dat$time=='T5' & dat$source=='S',]))
summary(aov(cellcounts ~ DOM, data=dat[dat$time=='T5' & dat$source=='C',]))
summary(aov(cellcounts ~ DOM, data=dat[dat$time=='T5' & dat$source=='CS',]))
                 