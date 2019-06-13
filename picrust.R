rm(list=ls())

#@GUI: you don't need PART1 of this script, but can directly start working on part two

##################################
#PART1
#create input files for picrust (not more than two digits!!)
setwd("/Users/sara/Documents/DFG/cryofreezing/picrust/")
ref.stats<- read.csv("JGIuploads/stats.all.csv", row.names=1)
taxids <- read.table("16S.corr.txt", header=T)

ref.stats.taxids <- merge (taxids, ref.stats, by.x='assembly.corr', by.y=0)
ref.stats.taxids[,'assembly.corr']<-factor(ref.stats.taxids[,'assembly.corr'])#numeric to factor

colnames(ref.stats.taxids)[10] <- "genome.size"
colnames(ref.stats.taxids)[3] <- "x16S_rRNA_count.pic"
colnames(ref.stats.taxids)[7] <- "strain"
colnames(ref.stats.taxids)[11] <- "gene.count"
colnames(ref.stats.taxids)[12] <- "x16S_rRNA_count.jgi"
colnames(ref.stats.taxids)[13] <- "tRNA_count"
colnames(ref.stats.taxids)[14] <- "otherRNA_count"
colnames(ref.stats.taxids)[15] <- "HGT_count"
colnames(ref.stats.taxids)[16] <- "HGT_perc"

#compare 16s count from picrust file and jgi output
diff16S <-ref.stats.taxids[ref.stats.taxids$x16S_rRNA_count.pic!=ref.stats.taxids$x16S_rRNA_count.jgi,][,c(1,2,3,10,12,7)]

genomesize <-ref.stats.taxids[,c(2,10)]
genomesize.mod <-ref.stats.taxids[,c(2,10)]
genomesize.mod$genome.size <- genomesize.mod$genome.size/100000
genomesize.mod2 <-ref.stats.taxids[,c(2,10)]
genomesize.mod2$genome.size <- round(genomesize.mod2$genome.size/100000,0)
genomesize.mod3 <-ref.stats.taxids[,c(2,10)]
genomesize.mod3$genome.size <- round(genomesize.mod3$genome.size/10000,0)

#write.table(genomesize, 'genomesize.txt', row.names=F, sep ='\t',quote=FALSE)
#write.table(genomesize.mod, 'genomesize.mod1.txt', row.names=F, sep ='\t',quote=FALSE)
#write.table(genomesize.mod2, 'genomesize.mod2.txt', row.names=F, sep ='\t',quote=FALSE)
#write.table(genomesize.mod3, 'genomesize.mod3.txt', row.names=F, sep ='\t',quote=FALSE)
#write.table(x16s, '16Smod2.txt', row.names=F, sep ='\t',quote=FALSE)


###################
#PART2
#evaluation of the picrust output data
#read count file
setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio")
counts <- read.table("dada.rcounts.ASV.2.tab", header=T, row.names=1)

###########################
#create dataframe with schema for expreimental setup
source <- c(rep('C',7),rep('CS',4),rep('S',10) )
DOM <- c('org',rep('T.M',3),rep('T.S',3),rep('T.M',2),rep('T.S',2),rep('Sea',3),rep('T.M',3),'org',rep('T.S',3))
rep <- c('r1','r1','r2','r3','r1','r2','r3','r1','r2','r2','r3','r1','r2','r3','r1','r2','r3','r1','r1','r2','r3')
schema<- cbind.data.frame(source,DOM,rep)
row.names(schema)<- colnames(counts)
###################
#16s rRNA gene copy number
pic.16s <- read.table("marker_predicted_and_nsti.tsv", header=T) #load picrust output file for 16S
counts.s <- counts[row.names(counts) %in% pic.16s[pic.16s$metadata_NSTI<1,1], ] #extract ASVs with NSTI<1 (= ASVs with no clode relative in the picrust reference database)
colSums (counts.s) #check which proportion of sequences is left after removing ASVs with NSTI<1

counts.s.rel <- as.data.frame.matrix(prop.table(t(t(counts.s)),2)) #re-normlize remaining count-data
colSums (counts.s.rel) #should sum up again to 1

counts.16s <- merge(pic.16s,counts.s.rel, by.x='sequence', by.y=0) #create a column with predicted 16s copy numbers for each ASV
row.names(counts.16s) <- counts.16s[,1]
counts.16s <-counts.16s[,c(2,4:24)] #select releavnt samples

av.16s <- colSums(counts.16s[,1]*counts.16s[,2:22]) #average number 0f 16s rRNA gene copy per sample
av.16s

###################
#genome size
genomesize <- read.table("genomesize_predicted.mod2", header=T)  #load picrust output file for genome size
counts.s <- counts[row.names(counts) %in% genomesize[genomesize$metadata_NSTI<1,1], ] #extract ASVs with NSTI<1
colSums (counts.s)#check which proportion of sequences is left after removing ASVs with NSTI<1

counts.s.rel <- as.data.frame.matrix(prop.table(t(t(counts.s)),2)) #re-normlize remaining count-data
colSums (counts.s.rel) #should sum up again to 1

counts.genomesize <- merge(genomesize,counts.s.rel, by.x='sequence', by.y=0)#create a column with predicted genome size for each ASV
row.names(counts.genomesize) <- counts.genomesize[,1]
counts.genomesize <-counts.genomesize[,c(2,4:24)] #select releavnt samples

av.genomesize <- colSums(counts.genomesize[,1]*counts.genomesize[,2:22]/10) #average genome size (in MIObp)
av.genomesize
#summary
traits<-cbind(schema, av.16s,av.genomesize ) #data frame with trait data and sample schema
traits[,-c(1,2,3)]<-round(traits[,-c(1,2,3)],1) #round because precision of picrust is only two digits!

write.table(traits, 'picrust.traits.tab', sep='\t',col.names=NA)

traits.mean <- aggregate(. ~ source + DOM, data=traits[,-3], FUN=mean) #data frame with trait mean values for each treatment
traits.mean[with(traits.mean, order(source, DOM)), ]


