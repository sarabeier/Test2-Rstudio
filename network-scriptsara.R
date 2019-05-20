rm(list=ls()) #cleans up working space
library(Hmisc)
library(igraph)

### Part 1 ###

#GUI: If you are using the subcounttables, you can skip Part 1 and go directly to Part 2 
#setwd("directory")
setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")

##another change
#SARA:another comment
  

#in the FlowJo analysis. Here, the gate "bact" was done with SCC-H vs FL1-H

FD<-flowDiv(myworkspaces= "workspace.wsp",gate.name="bact",dilutions=rep(10,n samples))

# 3 and 5 (SSC and FL1-H)
# we used 40

FD[["Matrices"]] #contain the count table

#SARA: before getting to the next setp please split the counttable into different 
#subcounttables, for each timeseries (each repliacte seperately!) create one counttable.
#SARA: to test the script maybe start for now only with one countable by subsetting 
#the full countttable




#### Part 2 ###


### GUI: my modifications on the code start from here


## GUI: creating a table to insert the results. Do it only once, at the first time
header<-c("Sample","Mean Degree", "Median Degree", "Density", "Diameter", "Modularity", 
          "Total Length", "Pos Lenght", "Neg Lenght" )
Results<-rbind(header)


#GUI: create an object with one replicate of a treatment, with all sampling days (from 0 to 6)
#GUI: I used only one of the subcounttables, but you just have to change the source .csv file to 
#analyse the other samples
counts = list.files(pattern="*.csv")

#########################################################################################################################
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
colnames(Results) <- header
Results

###########################
#create dataframe with schema for expreimental setup
source <- c(rep('C',6),rep('CS',6),rep('S',6) )
DOM <- c(rep('M',3),rep('S',3),rep('M',3),rep('S',3),rep('M',3),rep('S',3))
rep <- rep(c('r1','r2','r3'),6)
schema<- cbind.data.frame(source,DOM,rep)
schema$treat <-paste(source, DOM, sep='_')
Results <- cbind.data.frame(Results,schema)

#####################################################################################################
plot(SS.2.csv.net,vertex.size=2,edge.width=.4) 




write.csv(Results, 'Results.1.csv')


