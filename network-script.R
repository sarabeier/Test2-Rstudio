rm(list=ls()) #cleans up working space
library(Hmisc)
library(igraph)

#setwd("directory")
setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")

#list with names of count-files
counts = list.files(pattern="*.csv") 

#####################################################################################################
#PART1) loop with network analyses for all timeseries

datalist = list()
for (j in 2:length(counts)){
  #for (j in 10:11){
  df<-read.csv(counts[j], header = T, sep = ";")
  
  #SARA: samples count[11] and counts[16] are problematice because of each one bin = 1 for all time points
  #SARA: the two line below remove the outlier bins manually (does not work if running loop)
  #SARA: I have added below an if..else statement to remove these outliers automatically based on their variance=0
  
  #df <-df[ , -which(names(df) %in% c("Bin1057"))] #SARA: removes outlier for dataset 11
  #df <-df[ , -which(names(df) %in% c("Bin628"))] #SARA: removes outlier for dataset 16
  
  #define the name of each time series based on the filename
  filename<-substring(counts[j],7,16) 
  
  #finding zeros
  bincurve=matrix(NA,nrow=dim(df)[2],ncol=2)
  for (i in 1:dim(df)[2]){
    bincurve[i,]=c(i,length(which(df[,i]==0)))
  }
  
  #exclude all bins which are present in less of half datapoints of each time series
  counter=which(bincurve[,2]>3)
  FDclean=df[,-counter]
  dim(FDclean) #check the new dimensions
  
  
  #SARA: file 11 and 16 are problematic because they contain each one bin with 1 in for each time point (no variance), which messes up the network 
  #identifies bins without variance; 7 needs to be added because the first 7 columns are not considered in order to identify the righ column number
  out<-which(apply(FDclean[,-c(1:7)], 2, var)==0)+7 
  #if..else code to remove bins with variance=0 if they occur
  if (length(out)>=1) {
    FDclean <-FDclean[,-(which(apply(FDclean[,-c(1:7)], 2, var)==0)+7)] #removes bins without variance
  } else {
    FDclean <-FDclean #in case no bins without variance are detected FDclean is not changed
  }
  dim(FDclean) #check the new dimensions after removal of problematic bins
  MHC1=FDclean
  rich <-dim(MHC1)[2]-7 #richness of OTU bins after filtering out by coverage
  #Calculate the correlation values, filter by significance and by magnitude of correlation
  HCcor1=rcorr(as.matrix(MHC1[,8:dim(MHC1)[2]]),type="spearman") #SARA: replaceing the latter number by dim(MHC1)[2] generalizes this command for all input, time needs to be excluded!!
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
  assign(paste0(filename, '.net'), net1) #assignes new name to net1, which depends on filename (e.g. SM.3.csv.net), e.g. to reproduce plots
  
  # Indexes (see also https://en.wikipedia.org/wiki/Network_science)
  mean.degree<-mean(degree(net1)) #mean degree indicates the average number of edges for each bin
  median.degree<-median(degree(net1))
  #The density D {\displaystyle D} D of a network is defined as 
  #a ratio of the number of edges E {\displaystyle E} E to the number of possible edges in a network with N {\displaystyle N} N nodes,
  density<-edge_density(net1, loops=F) 
  #the diametr is defined as the shortest distance between the two most distant nodes in the network
  diameter<-diameter(net1, directed=F)
  
  ## GUI: biological networks usually have long-tailed distribution of frequencies of degree
  #SARA: in the histogram I could see that counttables 11 and 16 had each one bin single bin with >400 interactions 
  #SARA: this made it impossible to calculate modularity, while all other parameters seemed to be +/- ok
  hist.degree<-hist(degree(net1))  #histogram allows to detect outliers
  
  #Describing edges (total, positive and negatives)
  Total.lenght<-length(which(!is.nan(egdes.net))) #total number of edges in the network
  length(which(is.nan(egdes.net))) ## GUI: I believe that this one should always be zero
  Negative.length<-length(egdes.net[which(egdes.net<0)]) #number of edges indicating negative interactions
  Positive.length<-length(egdes.net[which(egdes.net>0)]) #number of edges indicating positive interactions
  
  #GUI: not sure yet how to explore betwenness, but we probably should
  bet=betweenness(net1, directed=F, weights=NA)
  ceb <- cluster_edge_betweenness(net1)
  Modularity<-modularity(ceb) #measure the strength of division of a network into modules 
  
  #results from the loop are summarized in d
  d<-data.frame(filename,rich, mean.degree,median.degree,density,diameter,Modularity,
                Total.lenght,Positive.length,Negative.length)
  #results in d from each round in the loop are stored in the data list
  datalist[[j]] <- d
}
#results stored in the list are transformed into a dataframe
Results = do.call(rbind, datalist)
colnames(Results)<-c("Sample","Richness","Mean.Degree", "Median.Degree", "Density", "Diameter", "Modularity", 
                     "Total.Length", "Pos.Lenght", "Neg.Lenght" )


#adds columns with the proportion of positive and negative interactions
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
#SARA: I have inspected variance and distribution of the residuals using the plot function
#SARA: based on visual inspection I have chosen the transformation of the data
res.aov1 <- aov(log(Mean.Degree) ~ source*DOM, data = Results)
summary (res.aov1)
TukeyHSD(res.aov1)
par(mfrow=c(2,2))
plot(res.aov1)

Results$prop.neg.logit<-(log(Results$prop.neg/(100-Results$prop.neg))) #logit transformation (http://strata.uga.edu/8370/rtips/proportions.html)
res.aov2 <- aov(prop.neg.logit ~ source*DOM, data = Results)
summary (res.aov2)
TukeyHSD(res.aov2)
plot(res.aov2)

res.aov5 <- aov(log(Diameter) ~ source*DOM, data = Results)
summary (res.aov5)
TukeyHSD(res.aov5)
plot(res.aov5)

res.aov3 <- aov(log(Density) ~ source*DOM, data = Results)
summary (res.aov3)
TukeyHSD(res.aov3)
plot(res.aov3)

#no significance
res.aov4 <- aov(Modularity ~ source*DOM, data = Results)
summary (res.aov4)
TukeyHSD(res.aov4)
plot(res.aov4)

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



