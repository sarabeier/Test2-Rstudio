### Part 1 ###
setwd("directory")

##another change
#yet another

#in the FlowJo analysis. Here, the gate "bact" was done with SCC-H vs FL1-H

FD<-flowDiv(myworkspaces= "workspace.wsp",gate.name="bact",dilutions=rep(10,n samples))

# 3 and 5 (SSC and FL1-H)
# we used 40

FD[["Matrices"]] #contain the count table

#### Part 2 ###
#finding zeros
bincurve=matrix(NA,nrow=dim(FD[["Matrices"]])[2],ncol=2)
for (i in 1:dim(FD[["Matrices"]])[2]){
  bincurve[i,]=c(i,length(which(FD[["Matrices"]][,i]==0)))
}

#but also the bins with 50% of the samples in order to reduce complexity and make the correlation matrix
counter=which(bincurve[,2]>n samples)
FDclean=FD[["Matrices"]][,-counter]
dim(FDclean) #check the new dimensions

### Analisis per replicate
#High DOM Control
library(Hmisc)
library(igraph)
MHC1=FDclean[which(clases$Treat=="C" & clases$DOM=="H" & clases$Rep==1),]

#Calculate the correlation values, filter by significance and by magnitude of correlation
HCcor1=rcorr(MHC1,type="spearman")
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
mean(degree(net1))

median(degree(net1))

edge_density(net1, loops=F)

diameter(net1, directed=F)

#Describing edges (total, positive and negatives)

length(which(!is.nan(egdes.net)))
length(which(is.nan(egdes.net)))
length(egdes.net[which(egdes.net<0)])
length(egdes.net[which(egdes.net>0)])

bet=betweenness(net1, directed=F, weights=NA)
ceb <- cluster_edge_betweenness(net1)
modularity(ceb) 