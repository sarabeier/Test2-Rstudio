
library(readxl)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(vegan)
library(tidyr)

setwd("C:/Users/Gui/Documents/Doutorado/Colaboração Sara Beier/Experiment 14/Exp14_5/Cytometry/FCS")

df<-read_excel("sites x bins.cyt2.exp14_5.xlsx",
               sheet = "sites x bins.cyt.exp14_5")
df<-as.data.frame(df)

df$Replicate<-as.factor(df$Replicate)
df$Time<-as.factor(df$Time)

str(df)
summary(df)
class(df)

time<-c("0","1","2","3","4","5")
repl<-c("1","2","3")

df.Canet<-df[df$Community=="Canet",] 

df.SS<-df[df$Treatment=="S_S",] 
df.SS.1<-df.SS[df.SS$Replicate==1,]
df.SS.2<-df.SS[df.SS$Replicate==2,]
df.SS.3<-df.SS[df.SS$Replicate==3,]

df.SM<-df[df$Treatment=="S_M",]
df.SM.1<-df.SM[df.SM$Replicate==1,]
df.SM.2<-df.SM[df.SM$Replicate==2,]
df.SM.3<-df.SM[df.SM$Replicate==3,]


df.CM<-df[df$Treatment=="C_M",] 
df.CM.1<-df.CM[df.CM$Replicate==1,]
df.CM.2<-df.CM[df.CM$Replicate==2,]
df.CM.3<-df.CM[df.CM$Replicate==3,]


df.CS<-df[df$Treatment=="C_S",] 
df.CS.1<-df.CS[df.CS$Replicate==1,]
df.CS.2<-df.CS[df.CS$Replicate==2,]
df.CS.3<-df.CS[df.CS$Replicate==3,]

df.SCS<-df[df$Treatment=="SC_S",]
df.SCS.1<-df.SCS[df.SCS$Replicate==1,]
df.SCS.2<-df.SCS[df.SCS$Replicate==2,]
df.SCS.3<-df.SCS[df.SCS$Replicate==3,]

df.SCM<-df[df$Treatment=="SC_M",]
df.SCM.1<-df.SCM[df.SCM$Replicate==1,]
df.SCM.2<-df.SCM[df.SCM$Replicate==2,]
df.SCM.3<-df.SCM[df.SCM$Replicate==3,]

write.csv2(df.Canet,"bin x canet.csv")
write.csv2(df.CM.1,"bin x CM.1.csv")
write.csv2(df.CM.2,"bin x CM.2.csv")
write.csv2(df.CM.3,"bin x CM.3.csv")

write.csv2(df.CS.1,"bin x CS.1.csv")
write.csv2(df.CS.2,"bin x CS.2.csv")
write.csv2(df.CS.3,"bin x CS.3.csv")

write.csv2(df.SM.1,"bin x SM.1.csv")
write.csv2(df.SM.2,"bin x SM.2.csv")
write.csv2(df.SM.3,"bin x SM.3.csv")

write.csv2(df.SS.1,"bin x SS.1.csv")
write.csv2(df.SS.2,"bin x SS.2.csv")
write.csv2(df.SS.3,"bin x SS.3.csv")

write.csv2(df.SCM.1,"bin x SCM.1.csv")
write.csv2(df.SCM.2,"bin x SCM.2.csv")
write.csv2(df.SCM.3,"bin x SCM.3.csv")

write.csv2(df.SCS.1,"bin x SCS.1.csv")
write.csv2(df.SCS.2,"bin x SCS.2.csv")
write.csv2(df.SCS.3,"bin x SCS.3.csv")
