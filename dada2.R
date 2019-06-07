#@GUI: 
#I have replaced all spaces in outputfiles with '_'
#all outputfiles from this script now start with dada.*
#all output files are exported with tab-separator to avoid inconsistencies with csv files (,/;) in windows and macintosh (with exeption of the fastafile)


#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

#BiocManager::install("DECIPHER", version = "3.8")

#.cran_packages <- c("ggplot2", "gridExtra")
#.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
#.inst <- .cran_packages %in% installed.packages()
#if(any(!.inst)) {
#  install.packages(.cran_packages[!.inst])
#}
#.inst <- .bioc_packages %in% installed.packages()
#if(any(!.inst)) {
#  source("http://bioconductor.org/biocLite.R")
#  biocLite(.bioc_packages[!.inst], ask = F)
#}

# Load packages into session, and print package version
#sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

#install.packages("devtools")
#library(devtools)
#install_github("microbiome/microbiome")

###############################################
rm(list=ls()) #cleans up working space

library(microbiome)
library(dada2)
library(phyloseq)
library(ggplot2)
library(gridExtra)
library(DECIPHER)
library(plyr)
library(phangorn)

setwd("/Users/ccma_ufscar/Documents/Usuarios/Guilherme/Paired-End sequences/Sequences")
path<-("/Users/ccma_ufscar/Documents/Usuarios/Guilherme/Paired-End sequences/Sequences")
setwd("/Users/sara/Documents/DFG/cryofreezing/sequencing/PrimerClipped")
setwd("/Users/sara/Documents/R-scripts/cry/Test2-Rstudio/dada2.R")
path<-("/Users/sara/Documents/DFG/cryofreezing/sequencing/PrimerClipped")


list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[9:12])
plotQualityProfile(fnRs[9:12])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


## the truncation length (truncLen) is based on what you see with `plotQualityProfile`
#all other arguments are on their default levels
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,180),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


#Learn the error rates of each base
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Plot estimated errors 
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Remove replications
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)

dadaFs[[1]]
# 413 sequence variants were inferred from 8012 input unique sequences.
dadaRs[[1]]
# 420 sequence variants were inferred from 5737 input unique sequences.

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#[1]   21 5306

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#362  364  365  366  368  369  370  371  372  373  374  375  376  377  378  379  380  381 
#1    1    2    347  239    3  667  138 2369 1221   68  140   88    7    7    2    5    1  

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# 21 861
write.table(seqtab.nochim,file="dada.Seqtab.nochim.tab", sep = "\t")


sum(seqtab.nochim)/sum(seqtab)
#  0.8526282
#We lost a lot of sequences. But when we account for the abundances of 
#those variants we see they account for only about 15% of the merged sequence reads


#Track the number of reads that made it through so far
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.table(track,file="dada.Number_of_reads_per_sample.tab", sep = "\t") 

#Assign taxonomy. You have to download the `silva_nr_v128_train_set.fa.gz` file
#and put it inside the sae folder of your samples. See DADA2 documentation to download it
#remember to only use multithread=TRUE on Mac, not with Windows based computers


taxa <- assignTaxonomy(seqtab.nochim, "/Users/ccma_ufscar/Documents/Usuarios/Guilherme/Paired-End sequences/Sequences/silva_nr_v128_train_set.fa.gz",  multithread=TRUE)
taxa <- assignTaxonomy(seqtab.nochim, "/Users/sara/Documents/Silva/silva_nr_v128_train_set.fa.gz",  multithread=TRUE) #path on Sara's Mac
write.table(taxa, file="dada.taxonomic.table.tab",sep = "\t")


#Inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
taxa.print.names<-as.data.frame(taxa.print)

#library(seqinr) --> to export sequences for aligment in other programs

#seqnum <- paste0("Seq", seq(ncol(seqtab.nochim)))
#taxa.tree<-taxa[,2:4]
#uniqueSeqs <- as.list(colnames(seqtab.nochim))
#write.fasta(uniqueSeqs, taxa.tree, "/Users/ccma_ufscar/Documents/Usuarios/Guilherme/Paired-End sequences/Sequences/uniqueSeqs2.fasta")


#simple data.frame construction from the information encoded in the filenames
#Sample names are the 16S sequences -->ps0
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1) #usuful to see your samples and set the metadata (in my case Community, DOM and Treatment)
community <-c("C", "C", "C", "C", "C", "C", "C", "SC","SC", "SC", "SC", 
              "S", "S", "S", "S", "S", "S", "S","S", "S", "S") 
DOM<-c("C_initial", "M", "M", "M", "S", "S", "S", "M","M", "S", "S", 
       "SW", "SW", "SW", "M", "M", "M", "S_initial","S" ,"S", "S")
replicate<-c("1", "1" ,  "2"  , "3" ,  "1"  , "2" ,  "3"  , "1" , "2",  "2",  "3",
             "1",  "2",  "3",  "1",  "2"  , "3",   "1",  "1",   "2",   "3"  ) 
treatment<-c("C_initial", "CM", "CM", "CM", "CS", "CS", "CS", "SCM","SCM", "SCS", "SCS", 
             "SW", "SW", "SW", "SM", "SM", "SM", "S_initial","SS", "SS", "SS")

samdf <- data.frame(Community=community, DOM=DOM, Replicate=replicate, Treatment=treatment)
rownames(samdf) <- samples.out

##phylosq object from dada2 output
ps0 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                sample_data(samdf), 
                tax_table(taxa))
ps0 <- prune_samples(sample_names(ps0) != "Mock", ps0) # Remove mock sample (if you have one)
ps0

###

##construction of a phylogenetic tree
#alignment of sequences (if you exported, can be done in, e.g., the Geneious program)
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree. This could be a problem because 
#the default seqs names are the actual 16S sequence. Therefore, you might need to change these names
names(seqs)<- paste0("SV_", seq(ntaxa(ps0)), "_", taxa.print.names$Order)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)
seqs.tab<-as.data.frame(seqs)


phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
write.tree(treeNJ, file = "dada.treeNJ") #export tree

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR2, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

###

#phyloseq analysis
#simple data.frame construction from the information encoded in the filenames
#here I changed the sample names so I don't have the 16S sequences as names -->ps
seqtab.nochim2<-seqtab.nochim
colnames(seqtab.nochim2)<-rownames(seqs.tab)
taxa2<-taxa
rownames(taxa2)<-rownames(seqs.tab)

#I create another phyloseq object with the data with the new sample names.
#same procedure as before

#samples.out <- rownames(seqtab.nochim)
samples.out2 <- rownames(seqtab.nochim2)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1) #usuful to see your samples and set the metadata (in my case Community, DOM and Treatment)
community <-c("C", "C", "C", "C", "C", "C", "C", "SC","SC", "SC", "SC", 
              "S", "S", "S", "S", "S", "S", "S","S", "S", "S") 
DOM<-c("C_initial", "M", "M", "M", "S", "S", "S", "M","M", "S", "S", 
       "SW", "SW", "SW", "M", "M", "M", "S_initial","S" ,"S", "S")
replicate<-c("1", "1" ,  "2"  , "3" ,  "1"  , "2" ,  "3"  , "1" , "2",  "2",  "3",
             "1",  "2",  "3",  "1",  "2"  , "3",   "1",  "1",   "2",   "3"  ) 
treatment<-c("C_initial", "CM", "CM", "CM", "CS", "CS", "CS", "SCM","SCM", "SCS", "SCS", 
             "SW", "SW", "SW", "SM", "SM", "SM", "S_initial","SS", "SS", "SS")

samdf <- data.frame(Community=community, DOM=DOM, Replicate=replicate, Treatment=treatment)
rownames(samdf) <- samples.out2

#phylosq object from dada2 output
#with `seqtab.nochim2`, with the changed names (not the 16S sequences)
ps <- phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa2),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps


#extract count table (raw count data)
OTU1 = as(otu_table(ps), "matrix") #temporary file
# transpose if necessary
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf.ps = as.data.frame(OTU1)
row.names(OTUdf.ps) <-substring(row.names(OTUdf.ps),15,22) #simplifies sample names
write.table(OTUdf.ps, sep='\t', col.names=NA,quote=FALSE,"dada.counts.ASV.1.tab") ##abundance table @GUI: I have changed the file name to counts because it's the raw count data
write.table(tax_table(ps), file = "dada.Taxonomy.ASV.tab", sep = "\t") ## taxonomy table

#extract relative abundance table
ps.rel = transform_sample_counts(ps, function(x) x/sum(x))
rOTU1 = as(otu_table(ps.rel), "matrix") #temporary file
# transpose if necessary
if(taxa_are_rows(ps.rel)){rOTU1 <- t(rOTU1)}
# Coerce to data.frame
rOTUdf.ps = as.data.frame(rOTU1)
row.names(rOTUdf.ps) <-substring(row.names(rOTUdf.ps),15,22) #simplifies sample names
write.table(rOTUdf.ps, sep='\t', col.names=NA,quote=FALSE,"dada.rcounts.ASV.1.tab")
write.table(t(rOTUdf.ps), sep='\t', col.names=NA,quote=FALSE,"dada.rcounts.ASV.2.tab") #I expert also the transposed data in this case, beause I need it like this for the picrust

#extract fasta file
seqs <- getSequences(seqtab.nochim)
names(seqs)<- paste0(">SV_", seq(ntaxa(ps)), "_", taxa.print.names$Order)
write.table(seqs, "dada.seqs.nochim.fasta",quote=FALSE, col.names=FALSE)


save.image("dada2_r")
# Close R, Re-open R
load("dada2_r")

