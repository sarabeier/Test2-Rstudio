# from http://www.sthda.com/english/wiki/two-way-anova-test-in-r

rm(list=ls())
library(car)
library(dplyr)
library(ggplot2)
library(readxl)

#setwd("C:/Users/Gui/Documents/Doutorado/Colabora??o Sara Beier/Experiment 14/Exp14_5")
setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")

totalD<-read_excel("day5.xlsx",
                 sheet ="Plan3")

my_data<- read.csv("exp14_5.csv",sep=";" )

str(my_data)

as.factor(my_data$Time) 
as.factor(totalD$Community)
as.factor(totalD$DOM)
as.factor(totalD$Treatment)
as.numeric(totalD$Richness)
as.numeric(totalD$Shannon)

head(my_data)
str(my_data)
table(my_data$Community, my_data$DOM)

str(totalD)

day0<- my_data[my_data$Time==0,]
day1<-my_data[my_data$Time==1,]
day2<-my_data[my_data$Time==2,]
day3<-my_data[my_data$Time==3,]
day4<-my_data[my_data$Time==4,]
day5<-my_data[my_data$Time==5,]

ggplot(day0, aes(x=Community, y=Value, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  ggtitle("Initial Cell Abundance") +
  xlab("Community") + ylab("Cell Abundance")


ggplot(day1, aes(x=Community, y=Value, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  ggtitle("Cell Abundance Day1") +
  xlab("Community") + ylab("Cell Abundance")


ggplot(day2, aes(x=Community, y=Value, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  ggtitle("Cell Abundance Day2") +
  xlab("Community") + ylab("Cell Abundance")



ggplot(day3, aes(x=Community, y=Value, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  ggtitle("Cell abundance Day3") +
  xlab("Community") + ylab("Cell Abundance")


ggplot(day4, aes(x=Community, y=Value, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  ggtitle("Cell abundance Day4") +
  xlab("Community") + ylab("Cell Abundance")



ggplot(day5, aes(x=Community, y=Value, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  ggtitle("Final Cell Abundance Day 5") +
  xlab("Community") + ylab("Cell Abundance")


ggplot(totalD, aes(x=Community, y=BP, color=DOM)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  ggtitle("Final BP Day 5") +
  xlab("Community") + ylab("BP")


#if the test of homogeinity of variances fail, ANOVA is not the appropriate test.
#t.test( day3$Value[day3$Treatment=="S_S"], day5$Value[day3$Treatment=="S_M"])
#t.test( day3$Value[day3$Treatment=="C_S"], day5$Value[day3$Treatment=="C_M"])
#t.test( day3$Value[day3$Treatment=="SC_S"], day5$Value[day3$Treatment=="SC_M"])






#t.test( day5$Value[day5$Treatment=="S_S"], day5$Value[day5$Treatment=="S_M"])
#t.test( day5$Value[day5$Treatment=="C_S"], day5$Value[day5$Treatment=="C_M"])
#t.test( day5$Value[day5$Treatment=="SC_S"], day5$Value[day5$Treatment=="SC_M"])



# Two-way ANOVA with interaction effect
# These two calls are equivalent
#res.aov3 <- aov(Value ~ Community * DOM, data = day5)
#res.aov3 <- aov(len ~ supp + dose + supp:dose, data = my_data)

#If the interaction is not significant, use the additive model

#This call considers the factors additive 
res.aov5 <- aov(Value ~ Community + DOM, data = day5)
summary(res.aov5)

res.aov4 <- aov(Value ~ Community + DOM, data = day4)
summary(res.aov4)

res.aov3 <- aov(Value ~ Community + DOM, data = day3)
summary(res.aov3)


#test of homogeneity
#if p value is higher than 0.05, no evidence of that the variance across groups are
#statistically different. We assume homogeinity of variance in the different treatments groups
library(car)
leveneTest(Value ~ Community*DOM, data = day5)
leveneTest(Value ~ Community*DOM, data = day4)
leveneTest(Value ~ Community*DOM, data = day3)


#summary statistics
require("dplyr")
group_by(day4, Community, DOM) %>%
  summarise(
    count = n(),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE)
  )
T<-TukeyHSD(res.aov3, which = "Community")       
plot(T, las=1)

# 1. Homogeneity of variances
plot(res.aov3, 1)


# 2. Normality
plot(res.aov3, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov5)
# Run Shapiro-Wilk test for normality
shapiro.test(x = aov_residuals )
#if p greater than 0.05, the population is normally distributed