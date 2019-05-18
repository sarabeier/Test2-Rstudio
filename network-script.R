### Part 1 ###

#GUI: If you are using the subcounttables, you can skip Part 1 and go directly to Part 2 
setwd("directory")

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


df<-read.csv("bin x CM.3.csv", header = T, sep = ";")

#change the file name you are using here, to create a table with the network indices of the sub-counttable
filename<-"CM3"

#finding zeros
bincurve=matrix(NA,nrow=dim(df)[2],ncol=2)
for (i in 1:dim(df)[2]){
  bincurve[i,]=c(i,length(which(df[,i]==0)))
}


#but also the bins with 50% of the samples in order to reduce 
#complexity and make the correlation matrix
#SARA: you have to replace the 'n sample' by the number of samples in you subset counttable 
#dived by 2.
#SARA: When you have six samples, each for one sampled day 
#(I remember we had 6 days sampled, right), the write '3'
counter=which(bincurve[,2]>3)
FDclean=df[,-counter]
dim(FDclean) #check the new dimensions

### Analisis per replicate
#High DOM Control
library(Hmisc)
library(igraph)


MHC1=FDclean



#GUI: I didn't look to see the meaning of the parameters values that Angel chose. I'll give a look at it


#Calculate the correlation values, filter by significance and by magnitude of correlation
## GUI: here, change the upper limit of column numbers accordingly to the dimension of the 
#'FDclean'. There is a way to get this automatically, but I can't remember now
HCcor1=rcorr(as.matrix(MHC1[,7:675]),type="spearman")
HCcor1$P=p.adjust(HCcor1$P,method="BH") # To control the false positive
HCcor1$r[HCcor1$P>0.05]=0 #Onlysignificant values
HCcor1$r[abs(HCcor1$r)<.7]=0 #Only correlations values higher than 0.7

net1 <- graph_from_adjacency_matrix(HCcor1$r, weighted=T, mode="upper", diag=F) # to create the object
net1$layout <-layout.fruchterman.reingold # the way to organizae the nodes in the plot
V(net1)$label=NA # remove the rownames
E(net1)$color=ifelse(E(net1)$weight> 0, "blue","red") #give color to edges
egdes.net=E(net1)$weight #saving the edges in a new object
E(net1)$weight=abs(E(net1)$weight) # replace negative egdes for positive (only for plotting)
plot(net1,vertex.size=2,edge.width=.4) # figure



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


partial.results<-c(filename,mean.degree,median.degree,density,diameter,Modularity,
                   Total.lenght,Positive.length,Negative.length)
Results<-rbind(Results,partial.results)

##GUI: start over from line #43 with another sample. Remember to change accordingly the 'filename' object

