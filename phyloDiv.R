rm(list=ls()) #cleans up working space
setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")
library(geiger)


#######################################################################################################
#PD (phylogenetic diversity) and PDAW (abundance weighted phylogenetic divesrity)
#######################################################################################################
#see: http://www.daijiang.name/en/2014/05/04/notes-func-phylo-book-1/
#PART1: example---------------------------------------------------------------------------------
my.sample = read.table("PD.example.sample.txt", sep = "\t", row.names = 1, header = T)
my.phylo = read.tree("PD.example.phylo.txt")
pruned.tree = treedata(my.phylo, data = t(my.sample[1, my.sample[1, ] > 0]), 
                       warnings = F)$phy #prunes tree for species that are present in first sample
#PD for one sample
PD1 <-sum(pruned.tree$edge.length)
PD1

#PD for all samples
PDall<-apply(my.sample, 1, function(x) {
  sum(treedata(my.phylo, x[x > 0])$phy$edge.length, warnings = F)
})
PDall

#PDAW
#prune tree
com.1.phylo = treedata(my.phylo, t(my.sample[1, my.sample[1, ] > 0]))$phy
plot.phylo(com.1.phylo)
nodelabels()
tiplabels()

branches = matrix(NA, nrow(com.1.phylo$edge), ncol = 4)
branches[, 1:2] = com.1.phylo$edge
branches[, 3] = com.1.phylo$edge.length
for (i in 1:nrow(branches)) {
  leaves.node = tips(com.1.phylo, branches[i, 2])
  branches[i, 4] = mean(t(my.sample[1, leaves.node]))
}
n.of.branches = nrow(com.1.phylo$edge)
denominator = sum(branches[, 4])
numerator = sum(branches[, 3] * branches[, 4])
weighted.fatith = n.of.branches * (numerator/denominator)
weighted.fatith #PDAW in one sample

#PART2: PDAW real data---------------------------------------------------------------------------------
tree = read.tree("dada.treeNJ")
samp = read.table("dada.rcounts.ASV.1.tab")

pruned.tree = treedata(tree, data = t(samp[1, samp[1, ] > 0]), 
                       warnings = F)$phy #prunes tree for species that are present in first sample
PD1 <-sum(pruned.tree$edge.length)
PD1 #phylogenetic distance in 1st sample

#loop for all samples
PDAW <- NULL
for (n in 1:nrow(samp)) {
  phylo = treedata(tree, t(samp[n, samp[n, ] > 0]))$phy
  plot.phylo(phylo)
  nodelabels()
  tiplabels()
  branches = matrix(NA, nrow(phylo$edge), ncol = 4)
  branches[, 1:2] = phylo$edge
  branches[, 3] = phylo$edge.length
  for (i in 1:nrow(branches)) {
    leaves.node = tips(phylo, branches[i, 2])
    branches[i, 4] = mean(t(samp[n, leaves.node]))
  }
  number.of.branches = nrow(phylo$edge)
  denominator = sum(branches[, 4])
  numerator = sum(branches[, 3] * branches[, 4])
  x = number.of.branches * (numerator/denominator)
  a <- as.list(x)
  if(n==1)PDAW<-a else PDAW=c(PDAW,a)
}

PDAWall <- as.data.frame(do.call(rbind, PDAW), row.names=rownames(samp))
PDAWall$treat <-substring(rownames(PDAWall),1,2)
PDAWall #output for PDAW in all samples

