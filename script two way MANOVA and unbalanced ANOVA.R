
rm(list=ls())
library(car)
library(dplyr)
library(ggplot2)



setwd("C:/Users/Gui/Documents/Doutorado/Colaboração Sara Beier/Experiment 14/Exp14_5")


totalD<-read_excel("day5.xlsx",
                   sheet ="Plan3")

as.factor(totalD$Community)
as.factor(totalD$DOM)
as.factor(totalD$Treatment)
as.numeric(totalD$Richness)
as.numeric(totalD$Shannon)

str(totalD)

ggplot(totalD, aes(x=Community, y=BP, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))  +
  xlab("Community") + ylab("Bacterial Production (ugC L-1 d-1)")

ggplot(totalD, aes(x=Community, y=BGE, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  xlab("Community") + ylab("Bacterial Growth Efficience (Percentage)")

ggplot(totalD, aes(x=Community, y=BR, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  #ggtitle("BR") +
  xlab("Community") + ylab("Bacterial Respiration (ugC L-1 d-1)")


ggplot(totalD, aes(x=Community, y=Cell_abund, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  #ggtitle("BR") +
  xlab("Community") + ylab("Cell abundance (cell ml-1)")



ggplot(data = na.omit(subset(totalD,select = c(Samples,Community,DOM, Richness))), aes(x=Community, y=Richness, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, na.rm = TRUE, outlier.size=4) + 
              stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
              #ggtitle("Richness") +
              xlab("Community") + ylab("Richness")


ggplot(data = na.omit(subset(totalD,select = c(Samples,Community,DOM, Shannon))), aes(x=Community, y=Shannon, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, na.rm = TRUE, outlier.size=4) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  #ggtitle("Shannon Diversity") +
  xlab("Community") + ylab("Shannon Index")

ggplot(data = na.omit(subset(totalD,select = c(Samples,Community,DOM, InvSimpson))), aes(x=Community, y=InvSimpson, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, na.rm = TRUE, outlier.size=4) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  #ggtitle("Shannon Diversity") +
  xlab("Community") + ylab("Inverse Simpson Index")


#test of homogeneity
#if p value is higher than 0.05, no evidence of that the variance across groups are
#statistically different. We assume homogeinity of variance in the different treatments groups
library(car)
leveneTest(Cell_abund ~ Community*DOM, data = totalD)
leveneTest(BP ~ Community*DOM, data = totalD)
leveneTest(BR ~ Community*DOM, data = totalD)
leveneTest(BGE ~ Community*DOM, data = totalD)
leveneTest(Richness ~ Community*DOM, data = totalD)
leveneTest(Shannon ~ Community*DOM, data = totalD)


#MANOVA
cell_abund<-totalD$Cell_abund
BP<-totalD$BP
BR<-totalD$BR
BGE<-totalD$BGE
Richness<-totalD$Richness
Shannon<-totalD$Shannon

res.man <- manova(cbind(cell_abund, BP, BR, BGE, Richness, Shannon) ~ Community*DOM, data = totalD)
summary(res.man)


#ANOVA
res.aovBP <- aov(BP ~ Community * DOM, data = totalD)
summary(res.aovBP)

res.aovBR <- aov(BR ~ Community * DOM, data = totalD)
summary(res.aovBR)

res.aovBGE <- aov(BGE ~ Community * DOM, data = totalD)
summary(res.aovBGE)

res.aovcellab <- aov(Cell_abund ~ Community * DOM, data = totalD)
summary(res.aovcellab)


#unbalanced ANOVA
res.aovRich <- aov(Richness ~ Community * DOM, data = totalD)
Anova(res.aovRich, type = "III")

res.aovSha <- aov(Shannon ~ Community * DOM, data = totalD)
Anova(res.aovSha, type = "III")


