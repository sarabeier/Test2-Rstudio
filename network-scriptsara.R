rm(list=ls()) #cleans up working space
library(Hmisc)
library(igraph)

#setwd("directory")
setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")
counts = list.files(pattern="*.csv")

#####################################################################################################
#loop with network analyses for all timeseries

datalist = list()
for (j in 2:length(counts)){
  #for (j in 10:11){
  df<-read.csv(counts[j], header = T, sep = ";")

  #df <-df[ , -which(names(df) %in% c("Bin1057"))] #SARA: removes outlier for dataset 11
  #df <-df[ , -which(names(df) %in% c("Bin628"))] #SARA: removes outlier for dataset 16

  #change the file name you are using here, to create a table with the network indices of the sub-counttable
  #filename<-"CM2"
  filename<-substring(counts[j],7,16)

  #finding zeros
  bincurve=matrix(NA,nrow=dim(df)[2],ncol=2)
  for (i in 1:dim(df)[2]){
    bincurve[i,]=c(i,length(which(df[,i]==0)))
  }

  #but also the bins with 50% of the samples in order to reduce 
  #complexity and make the correlation matrix
  counter=which(bincurve[,2]>3)
  FDclean=df[,-counter]
  dim(FDclean) #check the new dimensions

  #SARA: file 11 and 16 are problematic because they contain each one bin with 1 in for each time point (no variance), which messes up the network 
  out<-which(apply(FDclean[,-c(1:7)], 2, var)==0)+7 #SARA: identifies bins without variance; 7 needs to be added because the first 7 columns are not considered
  if (length(out)>=1) {
    FDclean <-FDclean[,-(which(apply(FDclean[,-c(1:7)], 2, var)==0)+7)] #removes bins without variance
    } else {
    FDclean <-FDclean
  }
  dim(FDclean) #check the new dimensions after removal of problematic bins
  MHC1=FDclean
  #Calculate the correlation values, filter by significance and by magnitude of correlation
  #HCcor1=rcorr(as.matrix(MHC1[,7:675]),type="spearman")
  HCcor1=rcorr(as.matrix(MHC1[,8:dim(FDclean)[2]]),type="spearman") #SARA: replaceing the latter number by dim(FDclean)[2] generalizes this command for all input, time needs to be excluded!!
  HCcor1$P=p.adjust(HCcor1$P,method="BH") # To control the false positive
  HCcor1$r[HCcor1$P>0.05]=0 #Onlysignificant values
  HCcor1$r[abs(HCcor1$r)<.7]=0 #Only correlations values higher than 0.7

  net1 <- graph_from_adjacency_matrix(HCcor1$r, weighted=T, mode="upper", diag=F) # to create the object
  net1$layout <-layout.fruchterman.reingold # the way to organizae the nodes in the plot
  V(net1)$label=NA # remove the rownames
  E(net1)$color=ifelse(E(net1)$weight> 0, "blue","red") #give color to edges
  egdes.net=E(net1)$weight #saving the edges in a new object
  E(net1)$weight=abs(E(net1)$weight) # replace negative egdes for positive (only for plotting)
  par(mfrow=c(2, 1)) #SARA: split panels to see both plots at once
  plot(net1,vertex.size=2,edge.width=.4) # figure
  assign(paste0(filename, '.net'), net1) #assignes new name to net1, which depends on filename (e.g. SM.3.csv.net)

  # Indexes
  mean.degree<-mean(degree(net1))
  median.degree<-median(degree(net1))
  density<-edge_density(net1, loops=F)
  diameter<-diameter(net1, directed=F)

  ## GUI: biological networks usually have long-tailed distribution of frequencies of degree
  hist.degree<-hist(degree(net1)) 

  #Describing edges (total, positive and negatives)

  Total.lenght<-length(which(!is.nan(egdes.net)))
  length(which(is.nan(egdes.net))) ## GUI: I believe that this one should always be zero
  Negative.length<-length(egdes.net[which(egdes.net<0)])
  Positive.length<-length(egdes.net[which(egdes.net>0)])

  #GUI: not sure yet how to explore betwenness, but we probably should
  bet=betweenness(net1, directed=F, weights=NA)
  ceb <- cluster_edge_betweenness(net1)
  Modularity<-modularity(ceb) 

  #partial.results<-c(filename,mean.degree,median.degree,density,diameter,Modularity,
  #                  Total.lenght,Positive.length,Negative.length)

  d<-data.frame(filename,mean.degree,median.degree,density,diameter,Modularity,
                   Total.lenght,Positive.length,Negative.length)
  datalist[[j]] <- d
}
Results = do.call(rbind, datalist)
colnames(Results)<-c("Sample","Mean.Degree", "Median.Degree", "Density", "Diameter", "Modularity", 
          "Total.Length", "Pos.Lenght", "Neg.Lenght" )

Results$prop.neg <- Results$Neg.Lenght/Results$Total.Length*100
Results$prop.pos <- Results$Pos.Lenght/Results$Total.Length*100

#create dataframe with schema for expreimental setup
source <- c(rep('C',6),rep('CS',6),rep('S',6) )
DOM <- c(rep('M',3),rep('S',3),rep('M',3),rep('S',3),rep('M',3),rep('S',3))
rep <- rep(c('r1','r2','r3'),6)
schema<- cbind.data.frame(source,DOM,rep)
schema$treat <-paste(source, DOM, sep='_')

Results <- cbind.data.frame(Results,schema)

write.csv(Results, 'Results.1.csv')

#####################################################################################################
#stats
res.aov1 <- aov(Mean.Degree ~ source*DOM, data = Results)
summary (res.aov1)
TukeyHSD(res.aov1)
par(mfrow=c(2,2))
plot(res.aov1)

res.aov2 <- aov(prop.neg ~ source*DOM, data = Results)
summary (res.aov2)
TukeyHSD(res.aov2)
par(mfrow=c(2,2))
plot(res.aov2)

Results$x<-(log(Results$prop.neg/(100-Results$prop.neg))) #logit transformation
res.aov2 <- aov(x ~ source*DOM, data = Results)
summary (res.aov2)
TukeyHSD(res.aov2)
plot(res.aov2)

res.aov3 <- aov(Density ~ source*DOM, data = Results)
summary (res.aov3)
TukeyHSD(res.aov3)
plot(res.aov3)

#no significance
res.aov <- aov(Modularity ~ source*DOM, data = Results)
summary (res.aov)
TukeyHSD(res.aov)


#####################################################################################################
#plots
par(mfrow=c(3,3))
plot(SS.1.csv.net,vertex.size=2,edge.width=.4, main=c('SS.1')) 
plot(SS.2.csv.net,vertex.size=2,edge.width=.4, main=c('SS.2')) 
plot(SS.3.csv.net,vertex.size=2,edge.width=.4, main=c('SS.3')) 
plot(CS.1.csv.net,vertex.size=2,edge.width=.4, main=c('CS.1')) 
plot(CS.2.csv.net,vertex.size=2,edge.width=.4, main=c('CS.2')) 
plot(CS.3.csv.net,vertex.size=2,edge.width=.4, main=c('CS.3')) 
plot(SCS.1.csv.net,vertex.size=2,edge.width=.4, main=c('SCS.1')) 
plot(SCS.2.csv.net,vertex.size=2,edge.width=.4, main=c('SCS.2')) 
plot(SCS.3.csv.net,vertex.size=2,edge.width=.4, main=c('SCS.3')) 

par(mfrow=c(3,3))
plot(SM.1.csv.net,vertex.size=2,edge.width=.4, main=c('SM.1')) 
plot(SM.2.csv.net,vertex.size=2,edge.width=.4, main=c('SM.2')) 
plot(SM.3.csv.net,vertex.size=2,edge.width=.4, main=c('SM.3')) 
plot(CM.1.csv.net,vertex.size=2,edge.width=.4, main=c('CM.1')) 
plot(CM.2.csv.net,vertex.size=2,edge.width=.4, main=c('CM.2')) 
plot(CM.3.csv.net,vertex.size=2,edge.width=.4, main=c('CM.3')) 
plot(SCM.1.csv.net,vertex.size=2,edge.width=.4, main=c('SCM.1')) 
plot(SCM.2.csv.net,vertex.size=2,edge.width=.4, main=c('SCM.2')) 
plot(SCM.3.csv.net,vertex.size=2,edge.width=.4, main=c('SCM.3')) 



