library(dada2)
library(stringr)
path<-"~/Desktop/PhD/Metagenome/Bacterial_Seqs"
list.files(path)

fnFs<-sort(list.files(path, pattern="_R1_001.fastq.gz"))
fnRs<-sort(list.files(path, pattern="_R2_001.fastq.gz"))

pattern<-"GLU..."
sample.names<-str_match(fnFs, pattern)#works!!

sample.names

setwd(path)

plotQualityProfile(fnFs[4])
plotQualityProfile(fnRs[3])

filtFs<-file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs<-file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,150), trimLeft = c(18,19), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

out

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
rownames(seqtab.nochim)


uniquesToFasta(seqtab.nochim, fout='~/Desktop/PhD/Metagenome/rep-seqs.fna', ids=colnames(seqtab.nochim))

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE)

samdf<-read.csv("~/Desktop/PhD/Metagenome/GLUSEEN_HM2_metadata.csv")
head(samdf)
rownames(samdf)<-samdf$Sample_ID

library(ape)
#export sequences in a fasta file for QIIME tree building
uniquesToFasta(seqtab.nochim, fout='~/Desktop/PhD/Metagenome/rep-seqs.fna', ids=colnames(asv))
# alignment and tree run in QIIME
trefile = read.tree("phylo_tree.tre")# read tree from QIIME
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa),
               phy_tree(trefile))
saveRDS(ps, "~/Desktop/PhD/Metagenome/Bac16s_raw.RDS")

subsample<-round(50000*(sample_data(ps)$QPCR_16s/(mean(sample_data(ps)$QPCR_16s))))

ss<-sample_sums(ps)

QPCR16s<-rrarefy(otu_table(ps), subsample)
ps.rare<-ps
otu_table(ps.rare)<-QPCR16s

saveRDS(ps.rare, "~/Desktop/PhD/Metagenome/Bac16s_rarefied.RDS")
