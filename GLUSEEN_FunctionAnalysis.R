# import old Bac/Fun OTU table, load metadata
metadata<-read.csv("~/Desktop/PhD/Metagenome/GLUSEEN_HM2_metadata.csv")
metadata<-data.frame(metadata)
plot(metadata$Cd_tot~metadata$Codes)
par("mar")
par(mar=c(1,1,1,1))


# metadata stats ####

Zn.tot<-summary(aov(metadata$Zn_tot~metadata$Codes*metadata$Cities))
Ni.tot<-summary(aov(metadata$Ni_tot~metadata$Codes*metadata$Cities))
Co.tot<-summary(aov(metadata$Co_tot~metadata$Codes*metadata$Cities))
Cd.tot<-summary(aov(metadata$Cd_tot~metadata$Codes*metadata$Cities))

Zn.av<-summary(aov(metadata$Zn_avail~metadata$Codes*metadata$Cities))
Ni.av<-summary(aov(metadata$Ni_avail~metadata$Codes*metadata$Cities))
Co.av<-summary(aov(metadata$Co_avail~metadata$Codes*metadata$Cities))
Cd.av<-summary(aov(metadata$Cd_avail~metadata$Codes*metadata$Cities))

Zn.av
Ni.av
Co.av
Cd.av

Zn.tot
Ni.tot
Co.tot
Cd.tot

# p adjustment
Zn.tot$p.adj<-p.adjust(Zn.tot$Pr, method="bonferroni", n=4)
Ni.tot$p.adj<-p.adjust(Ni.tot$Pr, method="bonferroni", n=4)
Co.tot$p.adj<-p.adjust(Co.tot$Pr, method="bonferroni", n=4)
Cd.tot$p.adj<-p.adjust(Cd.tot$Pr, method="bonferroni", n=4)


#variance ? should do this ?
var(metadata$Cd_total, metadata$Codes)


# preprocess functional data ####

# template! and it works! ####

library(plyr)
e<-c("A", "A", "B", "B")
f<-c(1,1,1,1)
g<-c(2,2,4,4)
h<-c("a", "b", "c", "d")
test<-data.frame(e,f,g,h)

out<-ddply(test, .(a), summarize, # test is the input df, a is the column to summarize by
           b2=paste(mean(b)), #calculates the mean of column b
           c2=paste(mean(c)), #calculates the mean of column c
           d2=paste(unique(d), collapse=",")) # combines entries of column d as a vector
out

# dcast test conditions 

e<-c("A", "A", "A", "B", "B", "B", "B")
f<-c(2,2,2,2,2,2,2)
g<-c(2,2,2,4,4,4,4)
h<-c("a", "b", "c", "d","e","f","g")
test2<-data.frame(e,f,g,h)

l<-list(test, test2)

drp<-function(input){
  w<-ddply(input, .(e), summarize, 
           counts=mean(g),
           sample_id=unique(f),
           description=paste(unique(h), collapse=","))
  w
}

l.df<-ldply(l, drp)


drp<-function(input){
  w<-ddply(input, .(gene_id, Source), summarize, 
           counts=unique(gene_counts),
           misID=paste(unique(gene_counts), collapse=","),
           sample_id=unique(sample_id),
           norm_16sCounts=paste(unique(gene_16S_normalization),collapse=","),#??
           gene_category_id=paste(unique(gene_functional_category_id), collapse=","),
           gene_category=paste(unique(gene_functional_category), collapse=","),
           description=paste(unique(description), collapse=",")
           )
  o<-w$gene_id[grep(",",w$gene_16s_norm)]
  o
}
find<-ldply(list,derep)
finddoubles<-drp(sample_001)

sample_001<-read.csv("~/Desktop/PhD/Metagenome/Sample_001.csv")
sample_001<-data.frame(sample_001)
sample_001[sample_001]


sample_002<-read.csv("~/Desktop/PhD/Metagenome/Sample_002.csv")
sample_003<-read.csv("~/Desktop/PhD/Metagenome/Sample_003.csv")
sample_001<-data.frame(sample_001)
sample_002<-data.frame(sample_002)
sample_003<-data.frame(sample_003)
list.test2<-list("sample_001.csv","sample_002.csv","sample_003.csv")

test.tidy<-ldply(list.test, derep)
  
  #ddply(list, .(gene_id), summarize,
  #                    sample_id=paste(mean(sample_id)),
  #                    gene_counts=paste(mean(gene_counts)),
  #                    norm_16sCounts=paste(mean(gene_16S_normalization)),
  #                    function_ID=paste(unique(gene_functional_category_id), collapse=","),
  #                    description<-paste(unique(description), collapse=","))
# start Pipeline ####
derep<-function(input){
  a<-read.csv(paste0("~/Desktop/PhD/Metagenome/", input))
  a<-data.frame(a)
  a<-ddply(a, .(gene_id, Source), summarize, 
           #counts=mean(as.numeric(gene_counts), na.rm=TRUE),
           counts=unique(gene_counts),
           sample_id=unique(sample_id),
           #gene_16s_norm=paste(unique(gene_16S_normalization), collapse=","),
           gene_16s_norm=unique(gene_16S_normalization),
           gene_category_id=paste(unique(gene_functional_category_id), collapse=","),
           gene_category=paste(unique(gene_functional_category), collapse=","),
           description=paste(unique(description), collapse=",")) # works bc all calculations are done within a sample, and this works sample-wise through all samples

  return(a)
}#dereplicate counts in each sample
path<-"~/Desktop/PhD/Metagenome/"
#input<-"Sample_003.csv"
#Sample_test<-derep(path, input="Sample_003.csv")

make_tax_table<-function(input){
  
}

# automate over a list
library(plyr)
library(magrittr)
library(dplyr)
library(reshape2)

# make lists #
# A = Sequence run Number 1 #
# B = Sequence run Number 2 #
# list = samples run only in run 1 and not 2 #
listA<-as.list(list.files(path, pattern="A.csv"))
listB<-as.list(list.files(path, pattern="B.csv"))
list<-as.list(list.files(path, pattern="[^aAB].csv"))

#long.annotations<-ldply(list, derep) # combine sample datasets into one long dataframe
ListA.long.ann<-ldply(listA, derep)
ListB.long.ann<-ldply(listB, derep)
List.long.ann<-ldply(list, derep)

List.all<-rbind(ListA.long.ann,ListB.long.ann,List.long.ann)
nrow(ListA.long.ann)+nrow(ListB.long.ann)+nrow(List.long.ann)==nrow(List.all)

length(unique(List.all$gene_id)) # nice!!

# save progress #
saveRDS(List.long.ann, "~/Desktop/PhD/Metagenome/listC_longannotations.RDS")
saveRDS(ListA.long.ann, "~/Desktop/PhD/Metagenome/listA_longannotations.RDS")
saveRDS(ListB.long.ann, "~/Desktop/PhD/Metagenome/listB_longannotations.RDS")

write.csv(List.long.ann, "~/Desktop/PhD/Metagenome/listC_longannotations.csv")
write.csv(ListA.long.ann, "~/Desktop/PhD/Metagenome/listA_longannotations.csv")
write.csv(ListB.long.ann, "~/Desktop/PhD/Metagenome/listB_longannotations.csv")

List.long.ann<-read.csv("~/Desktop/PhD/Metagenome/listC_longannotations.csv", stringsAsFactors = F)
ListA.long.ann<-read.csv("~/Desktop/PhD/Metagenome/listA_longannotations.csv", stringsAsFactors = F)
ListB.long.ann<-read.csv("~/Desktop/PhD/Metagenome/listB_longannotations.csv", stringsAsFactors = F)

new.names1<-paste(List.all$gene_id[List.all$Source=="ARDB"], "_ARDB", sep="")
List.all$gene_id[List.all$Source=="ARDB"]<-new.names1

new.names2<-paste(List.all$gene_id[List.all$Source=="UniPkb"], "_UniPkb", sep="")
List.all$gene_id[List.all$Source=="UniPkb"]<-new.names2

new.names3<-paste(List.all$gene_id[List.all$Source=="BacMet"], "_BacMet", sep="")
List.all$gene_id[List.all$Source=="BacMet"]<-new.names3

new.names4<-paste(List.all$gene_id[List.all$Source=="Cazy"], "_Cazy", sep="")
List.all$gene_id[List.all$Source=="Cazy"]<-new.names4

saveRDS(List.all, "~/Desktop/PhD/Metagenome/All_annotations.RDS")
List.all<-readRDS("~/Desktop/PhD/Metagenome/All_annotations.RDS")

#datExpr<-dcast(long.annotations, sample_id~gene_id, value.var="counts", fill=0) # change from v to h
#troubleshooting####
err<-function(input){
  a<-read.csv(paste0("~/Desktop/PhD/Metagenome/", input))
  a<-data.frame(a)
  a<-paste(unique(a$sample_id),collapse=",")
  a 
}
sample_030<-read.csv("~/Desktop/PhD/Metagenome/Sample_030.csv")
sample_030$sample_id[which(sample_030$sample_id==NA)]
unique(sample_030$sample_id)
errors<-ldply(list, err)
# errors here; not sure why... it is a formating issue with one of the files. will have to seek and destroy 

# make ps tables ####
#TableA<-dcast(ListA.long.ann, sample_id~gene_id, value.var="counts", fill=0)
#TableB<-dcast(ListB.long.ann, sample_id~gene_id, value.var="counts", fill=0)
#Table<-dcast(List.long.ann, sample_id~gene_id, value.var="counts", fill=0)
reseq<-c(9,12,19:24,26,31:35,37:48,50,52:60,62:65,67,69:74,76,79:81,83,84,87,88,91,93,94,96:98)
seq1<-c(1:8,10,11,13:18,25,27:30,36,49,51,61,66,68,75,77,78,82,85,86,89,90,92,95,99,100)


Table<-data.frame(dcast(List.all, sample_id~gene_id, value.var="counts", fill=0))

# fixing some sh*t ####
table(List.all$Source)
unique(List.all$sample_id[List.all$Source=="ARBD"])
unique(List.all$sample_id[List.all$Source==""])
unique(List.all$sample_id[List.all$Source=="BacMat"])
unique(List.all$sample_id[List.all$Source=="BacMeta"])
unique(List.all$sample_id[List.all$Source=="MetaBac"])
unique(List.all$sample_id[List.all$Source=="MetaPkb"])
unique(List.all$sample_id[List.all$Source==" Cazy"])

List.all$Source[List.all$Source=="ARBD"]<-"ARDB"
List.all$Source[List.all$Source==""]<-"BacMet"
List.all$Source[List.all$Source=="BacMat"]<-"BacMet"
List.all$Source[List.all$Source=="BacMeta"]<-"BacMet"
List.all$Source[List.all$Source=="MetaBac"]<-"BacMet"
List.all$Source[List.all$Source=="MetaPkb"]<-"UniPkb"
List.all$Source[List.all$Source==" Cazy"]<-"Cazy"

#rownames(Table)<-Table$sample_id
#Table<-Table[,-1]
library(phyloseq)

otu<-data.frame(Table[62:161,])
otu$sample_id<-as.numeric(otu$sample_id)
otu<-otu[order(otu$sample_id),]
rownames(otu)<-otu$sample_id
otu<-otu[,-1]
otu<-otu_table(otu, taxa_are_rows = FALSE)

otu.reseq<-Table[1:61,]
otu.reseq$sample_id<-reseq
rownames(otu.reseq)<-otu.reseq$sample_id
otu.reseq<-otu.reseq[,-1]
otu.reseq<-otu_table(otu.reseq, taxa_are_rows = FALSE)
sample_names(otu)
sample_names(otu.reseq)

merged<-merge_phyloseq(otu, otu.reseq)

# make sample data ####
# I'm going to combine resequenced runs in phyloseq. First I need to make a matching sample data sheet to merge from
sam<-read.csv("~/Desktop/PhD/Metagenome/GLUSEEN_HM2_metadata.csv")
sam<-data.frame(sam)
sam$groups<-paste(sam$Cities, sam$Codes, sep=".")
rownames(sam)<-sam$site.number
sample.data<-sample_data(sam)
sample_names(sample.data)<-sam$site.number

# make tax table ####
tax<-ddply(List.all, .(gene_id), summarise, 
         gene_category_id=paste(unique(gene_category_id), collapse=","),
         gene_category=paste(unique(gene_category), collapse=","),
         description=paste(unique(description), collapse=","),
         source=unique(Source))

#needs to be done after correction from previous derep function is changed... description=description
#tax<-data.frame("Kingdom"=c(rep("Function", length(unique(List.all$gene_id)))), "Gene_ID"=unique(List.all$gene_id), "Function"=List.all$gene_category, "Function_ID"=List.all$gene_category_id) # not there yet
dim(tax)
names<-tax$gene_id
tax2<-tax[,c(5,3,4,1)]

tax<-tax_table(as.matrix(tax))
taxa_names(tax)<-names
length(taxa_names(tax))
head(taxa_names(tax))

tax2<-tax_table(as.matrix(tax2))
taxa_names(tax2)<-names

ps.meta<-phyloseq(merged, sample.data, tax)

psmSums<-sample_sums(ps.meta)
hist(sample_sums(ps.meta), breaks=30, xlab="Total Hits Per Sample", main = NULL)
sample_data(ps.meta)$metatotal<-psmSums


saveRDS(ps.meta, "~/Desktop/PhD/Metagenome/metagenome_ps.RDS")
ps.meta<-readRDS("~/Desktop/PhD/Metagenome/metagenome_ps.RDS")




#ps_bacarc<-read.RDS("")
#bsam<-sample_data(read.csv(""))
#ps_bacarc2<-read.RDS("")

#ps_bacarc<-merge_phyloseq(ps_bacarc, ps_bacarc2)
#
#ps_bacarc<-merge_samples()
#ps_bacarc<-phyloseq(ps_bacarc, bsam)

# sampling normalization functions ####

NormAdjmeta<-function(OTU, qpcr.val, t.hit1, t.hit2, set){ 
  # t.hit1=metagenomic 16s hits greengenes
  # t.hit2=total metagenomic anotations
  require(vegan)
  a<-sum(qpcr.val)/length(qpcr.val) # calculate average abundance
  b<-qpcr.val/a # calculate ratio of everage value to each sample
  ca<-t.hit2/t.hit1# calculate total metagenomic hits per 16s metagenomic hit
  val<-ca*b*set # set is the 16s rarefaction depth, ca*b is the adjustment for 16s abundance in the environment
  otu<-rrarefy(OTU, val)
  otu # return otu table
}

rarefy_uneven_meta<-function(ps, set){
  require(phyloseq)
  OTU<-otu_table(ps)
  qpcr.val<-sample_data(ps)$QPCR_16s
  t.hit1<-sample_data(ps)$GreenGenesHit
  t.hit2<-sample_data(ps)$metatotal
  otu<-NormAdjmeta(OTU, qpcr.val, t.hit1, t.hit2, set)
  otu_table(otu)
  otu
}

# Normalize data ####
otu.norm<-rarefy_uneven_meta(ps.meta, set=500)
otu.norm<-otu_table(otu.norm)
otu.norm2<-otu_table(otu.norm)
taxa_names(otu.norm)<-taxa_names(tax2)

ps.norm<-phyloseq(otu.norm, sample.data, tax2)
saveRDS(ps.norm, "~/Desktop/PhD/Metagenome/Rarefiedmetagenome_ps.RDS")

ps.norm2<-transform_sample_counts(ps.meta, function(x) x/sum(x))

# determine which samples to remove ####
a<-sum(sample_data(ps.norm)$QPCR_16s)/100 # calculate average abundance
b<-sample_data(ps.norm)$QPCR_16s/a # calculate ratio of everage value to each sample
ca<-sample_data(ps.meta)$metatotal/sample_data(ps.meta)$GreenGenesHit# calculate total metagenomic hits per 16s metagenomic hit
val<-ca*b*500
table(sample_sums(ps.meta)>val)

false<-c(9,20,24,33,38,41,59,83,93)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(Rmisc)
# make metal plots ####
se <- function(x) sqrt(var(x)/length(x))


arsenic<-subset_taxa(ps.norm, gene_category=="Arsenic")
arsenic<-tax_glom(arsenic, taxrank="gene_category")
merged.arsenic<-aggregate(otu_table(arsenic), list(sample_data(arsenic)$groups), mean)
SE<-aggregate(otu_table(arsenic), list(sample_data(arsenic)$groups), se)
merged.arsenic$SE<-SE$BAC0316_BacMet
merged.arsenic$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.arsenic$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.arsenic$Groups<-paste(merged.arsenic$Cities, merged.arsenic$Codes, sep=".")
merged.arsenic<-data.frame(merged.arsenic)

ggplot(data=merged.arsenic, mapping=aes(x=Groups,y=BAC0316_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0316_BacMet-SE, ymax=BAC0316_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Arsenic Resistance", y="Average Abundance")+coord_flip()



iron<-subset_taxa(ps.norm, gene_category=="Iron")
iron<-tax_glom(iron, taxrank="gene_category")
merged.iron<-aggregate(otu_table(iron), list(sample_data(iron)$groups), mean)
SE<-aggregate(otu_table(iron), list(sample_data(iron)$groups), se)
merged.iron$SE<-SE$BAC0003_BacMet
merged.iron$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.iron$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.iron$Groups<-paste(merged.iron$Cities, merged.iron$Codes, sep=".")
merged.iron<-data.frame(merged.iron)

ggplot(data=merged.iron, mapping=aes(x=Groups,y=BAC0003_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0003_BacMet-SE, ymax=BAC0003_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Iron Resistance", y="Average Abundance")+coord_flip()

Zinc<-subset_taxa(ps.norm, gene_category=="Zinc")
Zinc<-tax_glom(Zinc, taxrank="gene_category")
merged.Zinc<-aggregate(otu_table(Zinc), list(sample_data(Zinc)$groups), mean)
SE<-aggregate(otu_table(Zinc), list(sample_data(Zinc)$groups), se)
merged.Zinc$SE<-SE$BAC0646_BacMet
merged.Zinc$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Zinc$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Zinc$Groups<-paste(merged.Zinc$Cities, merged.Zinc$Codes, sep=".")
merged.Zinc<-data.frame(merged.Zinc)

ggplot(data=merged.Zinc, mapping=aes(x=Groups,y=BAC0646_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0646_BacMet-SE, ymax=BAC0646_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Zinc Resistance", y="Average Abundance")+coord_flip()

Nickel<-subset_taxa(ps.norm, gene_category=="Nickel")
Nickel<-tax_glom(Nickel, taxrank="gene_category")
merged.Nickel<-aggregate(otu_table(Nickel), list(sample_data(Nickel)$groups), mean)
SE<-aggregate(otu_table(Nickel), list(sample_data(Nickel)$groups), se)
merged.Nickel$SE<-SE$BAC0270_BacMet
merged.Nickel$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Nickel$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Nickel$Groups<-paste(merged.Nickel$Cities, merged.Nickel$Codes, sep=".")
merged.Nickel<-data.frame(merged.Nickel)

ggplot(data=merged.Nickel, mapping=aes(x=Groups,y=BAC0270_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0270_BacMet-SE, ymax=BAC0270_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Nickel Resistance", y="Average Abundance")+coord_flip()

Cadmium<-subset_taxa(ps.norm, gene_category=="Cadmium")
Cadmium<-tax_glom(Cadmium, taxrank="gene_category")
merged.Cadmium<-aggregate(otu_table(Cadmium), list(sample_data(Cadmium)$groups), mean)
SE<-aggregate(otu_table(Cadmium), list(sample_data(Cadmium)$groups), se)
merged.Cadmium$SE<-SE$BAC0058_BacMet
merged.Cadmium$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Cadmium$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Cadmium$Groups<-paste(merged.Cadmium$Cities, merged.Cadmium$Codes, sep=".")
merged.Cadmium<-data.frame(merged.Cadmium)

ggplot(data=merged.Cadmium, mapping=aes(x=Groups,y=BAC0058_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0058_BacMet-SE, ymax=BAC0058_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Cadmium Resistance", y="Average Abundance")+coord_flip()

Cobalt<-subset_taxa(ps.norm, gene_category=="Cobalt")
Cobalt<-tax_glom(Cobalt, taxrank="gene_category")
merged.Cobalt<-aggregate(otu_table(Cobalt), list(sample_data(Cobalt)$groups), mean)
SE<-aggregate(otu_table(Cobalt), list(sample_data(Cobalt)$groups), se)
merged.Cobalt$SE<-SE$BAC0092_BacMet
merged.Cobalt$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Cobalt$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Cobalt$Groups<-paste(merged.Cobalt$Cities, merged.Cobalt$Codes, sep=".")
merged.Cobalt<-data.frame(merged.Cobalt)

ggplot(data=merged.Cobalt, mapping=aes(x=Groups,y=BAC0092_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0092_BacMet-SE, ymax=BAC0092_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Cobalt Resistance", y="Average Abundance")+coord_flip()

Copper<-subset_taxa(ps.norm, gene_category=="Copper")
Copper<-tax_glom(Copper, taxrank="gene_category")
merged.Copper<-aggregate(otu_table(Copper), list(sample_data(Copper)$groups), mean)
SE<-aggregate(otu_table(Copper), list(sample_data(Copper)$groups), se)
merged.Copper$SE<-SE$BAC0133_BacMet
merged.Copper$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Copper$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Copper$Groups<-paste(merged.Copper$Cities, merged.Copper$Codes, sep=".")
merged.Copper<-data.frame(merged.Copper)

ggplot(data=merged.Copper, mapping=aes(x=Groups,y=BAC0133_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0133_BacMet-SE, ymax=BAC0133_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Copper Resistance", y="Average Abundance")+coord_flip()

#Manganese<-subset_taxa(ps.norm, gene_category=="Manganese")#no manganese
#Manganese<-tax_glom(Manganese, taxrank="gene_category")
#merged.Manganese<-aggregate(otu_table(Manganese), list(sample_data(Manganese)$groups), mean)
#SE<-aggregate(otu_table(Manganese), list(sample_data(Manganese)$groups), se)
#merged.Manganese$SE<-SE$BAC0058_BacMet
#merged.Manganese$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
#merged.Manganese$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
#merged.Manganese$Groups<-paste(merged.Manganese$Cities, merged.Manganese$Codes, sep=".")
#merged.Manganese<-data.frame(merged.Manganese)

#ggplot(data=merged.Manganese, mapping=aes(x=Groups,y=BAC0058_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0058_BacMet-SE, ymax=BAC0058_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Manganese Resistance", y="Average Abundance")+coord_flip()

Mercury<-subset_taxa(ps.norm, gene_category=="Mercury")
Mercury<-tax_glom(Mercury, taxrank="gene_category")
merged.Mercury<-aggregate(otu_table(Mercury), list(sample_data(Mercury)$groups), mean)
SE<-aggregate(otu_table(Mercury), list(sample_data(Mercury)$groups), se)
merged.Mercury$SE<-SE$BAC0652_BacMet
merged.Mercury$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Mercury$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Mercury$Groups<-paste(merged.Mercury$Cities, merged.Mercury$Codes, sep=".")
merged.Mercury<-data.frame(merged.Mercury)

ggplot(data=merged.Mercury, mapping=aes(x=Groups,y=BAC0652_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0652_BacMet-SE, ymax=BAC0652_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Mercury Resistance", y="Average Abundance")+coord_flip()

# make ardb plots ####

Multidrug<-subset_taxa(ps.norm, gene_category=="Resistance-nodulation-cell division transporter system. Multidrug resistance efflux pump.")
Multidrug<-tax_glom(Multidrug, taxrank="gene_category")
merged.Multidrug<-aggregate(otu_table(Multidrug), list(sample_data(Multidrug)$groups), mean)
SE<-aggregate(otu_table(Multidrug), list(sample_data(Multidrug)$groups), se)
merged.Multidrug$SE<-SE$ZP_02346647_ARDB
merged.Multidrug$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Multidrug$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Multidrug$Groups<-paste(merged.Multidrug$Cities, merged.Multidrug$Codes, sep=".")
merged.Multidrug<-data.frame(merged.Multidrug)

ggplot(data=merged.Multidrug, mapping=aes(x=Groups,y=ZP_02346647_ARDB, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=ZP_02346647_ARDB-SE, ymax=ZP_02346647_ARDB+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Multidrug Resistance", y="Average Abundance")+coord_flip()


ABC<-subset_taxa(ps.norm, gene_category=="ABC transporter system; bacitracin efflux pump.")
Mercury<-tax_glom(Mercury, taxrank="gene_category")
merged.Mercury<-aggregate(otu_table(Mercury), list(sample_data(Mercury)$groups), mean)
SE<-aggregate(otu_table(Mercury), list(sample_data(Mercury)$groups), se)
merged.Mercury$SE<-SE$BAC0652_BacMet
merged.Mercury$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Mercury$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Mercury$Groups<-paste(merged.Mercury$Cities, merged.Mercury$Codes, sep=".")
merged.Mercury<-data.frame(merged.Mercury)

ggplot(data=merged.Mercury, mapping=aes(x=Groups,y=BAC0652_BacMet, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=BAC0652_BacMet-SE, ymax=BAC0652_BacMet+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Mercury Resistance", y="Average Abundance")+coord_flip()

multidrug2<-subset_taxa(ps.norm, gene_category=="Resistance-nodulation-cell division transporter system. Multidrug resistance efflux pump.")
multidrug2<-tax_glom(multidrug2, taxrank="gene_category")
merged.multidrug2<-aggregate(otu_table(multidrug2), list(sample_data(multidrug2)$groups), mean)
SE<-aggregate(otu_table(multidrug2), list(sample_data(multidrug2)$groups), se)
merged.multidrug2$SE<-SE$ZP_02346647_ARDB
merged.multidrug2$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.multidrug2$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.multidrug2$Groups<-paste(merged.multidrug2$Cities, merged.multidrug2$Codes, sep=".")
merged.multidrug2<-data.frame(merged.multidrug2)

ggplot(data=merged.multidrug2, mapping=aes(x=Groups,y=ZP_02346647_ARDB, fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=ZP_02346647_ARDB-SE, ymax=ZP_02346647_ARDB+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Multidrug2 Resistance", y="Average Abundance")+coord_flip()

Bacitracin<-subset_taxa(ps.norm, description=="bacitracin")
Bacitracin2<-otu_table(Bacitracin)
merged.Bacitracin<-as.data.frame(rowSums(Bacitracin2))
merged.Bacitracin2<-aggregate(merged.Bacitracin, list(sample_data(Bacitracin)$groups), mean)
colnames(merged.Bacitracin2)<-c("Group1", "Sums")
SE<-aggregate(merged.Bacitracin, list(sample_data(Bacitracin)$groups), se)
colnames(SE)<-c("Group.1", "SE")
merged.Bacitracin2$SE<-SE$SE
merged.Bacitracin2$Codes<-c(rep(c("Reference", "Remnant", "Ruderal", "Turf"), 5))
merged.Bacitracin2$Cities<-c(rep("Baltimore", 4), rep("Budapest", 4), rep("Helsinki", 4), rep("Lahti", 4), rep("South Africa", 4))
merged.Bacitracin2$Groups<-paste(merged.Bacitracin2$Cities, merged.Bacitracin2$Codes, sep=".")
merged.Bacitracin2<-data.frame(merged.Bacitracin2)

ggplot(data=merged.Bacitracin2, mapping=aes(x=Groups,y=rowSums.Bacitracin2., fill=Cities))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=rowSums.Bacitracin2.-SE, ymax=rowSums.Bacitracin2.+SE))+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Bacitracin Resistance", y="Average Abundance")+coord_flip()


#statisticall test ####

# run deseq2 on dataset ####
diffDESeq<-function(ps){
  require(Deseq)
  require(Deseq2)
  require(phyloseq)
  L1 <-phyloseq_to_deseq2(ps, ~ Cities + Cities/Codes)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(L1), 1, gm_mean)
  diagdds = estimateSizeFactors(L1, geoMeans = geoMeans)
  diagdds = DESeq(diagdds, fitType="local")
  res = results(diagdds)
  return(res) 
}


# extract values, put in tax table ####
extractP<-function(ps, pvalue){
  require(phyloseq)
  if (identical(rownames(tax_table(ps)), rownames(pvalue))) {
    T1<-ps
    new_tax<-as.data.frame(as.matrix(tax_table(T1)))
    new_tax$DESeq_padj<-pvalue$padj
    tax_table(T1)<-tax_table(as.matrix(new_tax))
    T1
  }
  else {return("rownames out of order")}
}

Rare.deseq<-diffDESeq(ps.norm)

# add 

source->function_id/description->gene_id

# bacteria pipeline ####
BacDat<-data.frame(t(as.matrix(read.csv("/Volumes/Seagate/Backup/Thesis020717/Thesis/GithubFiles/Bacteria/Bacteria_OTU.csv"))))
dim(BacDat)
nrow(BacDat)
colnames(BacDat)
rownames(BacDat)
BacDat[1,]


FunDat<-read.csv("/Volumes/Seagate/Backup/Thesis020717/Thesis/GithubFiles/Fungi/Fungi_OTU_Main.csv")
dim(FunDat)

MetaDat<-read.csv("/Volumes/Seagate/Backup/Thesis020717/Thesis/GithubFiles/QPCR_Data_100615.csv")
dim(MetaDat)


# Limma-Voom ####
#https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

topTable(fit, coef=ncol(design))
limmaVoomPipe<-function(dat, design, contrasts){
  require(limma)
  require(edgeR)
  require(phyloseq)
  OTU<-as.matrix(t(otu_table(dat)))
  sData<-sample_data(dat)
  attach(sData)
  dge<-DGEList(counts=OTU)
  dge<-calcNormFactors(dge) #calculate the normalization factors via EdgeR (see EdgeR documentation)
  
  V <- voom(dge, design, plot=TRUE)
  Fit<-lmFit(V, design)
  Contr<-makeContrasts(categories , levels=colnames(coef(Fit)))
  Tmp<-contrasts.fit(Fit,Contr)
  Tmp<-eBayes(Tmp)
  Tmp2<-topTable(Tmp, coef=1, sort.by="P", n="INF")#add column for bayesian support
  Tmp2
  
}



# wgcna install ####
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html
library(WGCNA)
#install.packages("WGCNA")
#source("http://bioconductor.org/biocLite.R") 
#biocLite("impute")
#biocLite("preprocessCore")
#biocLite("GO.db")
# WGCNA ####

# filters genes based on total abundance in dataset (greater than 20), and total variance across samples 

# outputs a file for input into WGCNA functions

WGCNAfilt<-function(ps){
  require(phyloseq)
  
  #ps<-filter_taxa(ps, function(x) sum(x)>20, TRUE)
  datExpr<-otu_table(ps)
  
  allowWGCNAThreads()
  options(stringsAsFactors = FALSE)
  #colnames(datExpr)<-datExpr[,1]
  rownames(datExpr)<-datExpr$sample_id
  datExpr<-datExpr[,length(datExpr)]   # change to match length!!!
  # datExpr<-datExpr[,!(names(datExpr)%in%drops)]
  meanExpressionByArray=apply(datExpr, 1,mean, na.rm=T) # calculate mean sample-wise abundance
  barplot(meanExpressionByArray,
          xlab = "Sample", ylab = "Mean expression",
          main ="Mean expression across samples",
          cex.names = 0.7) # plot mean expression by sample
  # tutorial has a filter step of reducing the 
  #summary(NumberMissingByGene)
  variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T)) # 
  no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) ) #
  #table(no.presentdatExpr)
  KeepGenes= variancedatExpr>0 & no.presentdatExpr>=2 # set minimum variance for genes
  
  #table(KeepGenes)
  datExpr=datExpr[, KeepGenes]
  
  datExpr
}



library(WGCNA)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

# test WGCNA with/without filtering ####
datExpr<-otu_table(ps.meta)
datExpr2<-otu_table(ps.norm)
datExpr3<-otu_table(ps.norm2)
# Choose a set of soft-thresholding powers####
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

sft2 = pickSoftThreshold(datExpr3, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft2$fitIndices[,1], -sign(sft2$fitIndices[,3])*sft2$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence Rarefied"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft2$fitIndices[,1], sft2$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity Rarefied"))
text(sft2$fitIndices[,1], sft2$fitIndices[,5], labels=powers, cex=cex1,col="red")

# run blockwisemodules ####
datExprGOOD<-datExpr2[,goodSamplesGenes(datExpr2)$goodGenes==TRUE]
ncol(datExpr2)
ncol(datExprGOOD)

metanet<-blockwiseModules(datExprGOOD, weights=NULL,
                   checkMissingData = F,
                   maxBlockSize=10000,
                   corType = "p",
                   maxPOutliers = 1, 
                   quickCor = 0,
                   pearsonFallback = "individual",
                   cosineCorrelation = FALSE,
                   power = 9,
                   networkType = "signed",
                   replaceMissingAdjacencies = FALSE,
                   suppressTOMForZeroAdjacencies = FALSE,
                   TOMType = "signed",
                   #TOMDenom = "min",
                   deepSplit = 2,
                   detectCutHeight = 0.800, 
                   minModuleSize = 30,
                   reassignThreshold = 1e-6,
                   #minCoreKME = 0.5, 
                   #minCoreKMESize = minModuleSize/3,
                   #minKMEtoStay = 0.3,
                   mergeCutHeight = 0.25, 
                   impute = TRUE, 
                   trapErrors = FALSE 
                   )
table(metanet$colors)

datTraits<-as.matrix(sample_data(ps.norm))
datTraits<-datTraits[,c(5,8:22)]
names(datTraits)<-colnames(datTraits)
moduleColors=metanet$colors
ME=moduleEigengenes(datExprGOOD, metanet$colors)
MEs0 = ME$eigengenes #identify eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")# need to try with ancova!!


#moduleTraitAncova=lme(MEs~datTraits+sample_sums(ps.meta), random=~1|data.frame(datTraits)$Cities)
#moduleTraitAncova=aov(MEs~datTraits+sample_sums)


nSamples = nrow(datExprGOOD);
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               #ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

plotDendroAndColors(metanet$dendrograms[[1]], metanet$blockGenes[[1]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExprGOOD, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

City<-sample_data(ps.norm)$Cities
landuse<-sample_data(ps.norm)$Codes
an<-function(a){anova(lm(a~City+City/landuse+sample_sums(ps.meta)+sample_sums(ps.norm)))}
Ancova<-list()
Ancova<-lapply(MEs, an)

sink("~/Desktop/significant_modules_norm.txt")
lapply(MEs, an)
dev.off()
sigvectors<-c("MEdarkolivegreen","MElightpink3","MEturquoise","MEred,MEbrown","MEorangered1","MElightblue4","MEpalevioletred1","MEyellow2","MEpink3","MEtan4","MEskyblue3","MEskyblue1","MEnavajowhite","MEsalmon4","MEsteelblue","MEpink","MEblueviolet","MEsienna2","MEgrey")

names(datExpr3)[moduleColors=="MEturquoise"]
#anova(moduleTraitAncova)
geneTraitSignificance = as.data.frame(cor(datExpr3, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#moduleID2<-function(){} # correlate modules from one network to another; why do this if we can just combine networks?

WGCNApipe<-function(metanet){} # identify interesting modules (most positive and negative), then correlate to list of 
  #ADJ1=abs(cor(datExpr,use="p"))^6 # 
  #k=softConnectivity(datE=datExpr,power=6)
  #sizeGrWindow(10,5)
  #par(mfrow=c(1,2))
  #hist(k) # histogram of connectivity
  #scaleFreePlot(k, main="Check scale free topology\n") # plot connectivity
  
  # relate network to other factors ###
  
  # Define numbers of genes and samples
  nGenes = ncol(datExpr);
  nSamples = nrow(datExpr);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))

  # Define variable weight containing the weight column of datTrait
  weight = as.data.frame(datTraits$weight_g); #here we define what trait we are looking at correlating!!!
  names(weight) = "weight"
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
  names(GSPvalue) = paste("p.GS.", names(weight), sep="");
  
  module = "brown"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for body weight",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# WGCNA network pipeline ####
makeWGCNA<-function(name){
  a<-blockwiseModules(name) # figure out all of the blockwise function inputs!!
  saveRDS(a, paste0("path/", name))
}

l_ply(list, makeWGCNA)

list<- # look for all RDS with a particular pattern in folder
  
makeWcor<-function(datExpr, datTraits){
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p");
 
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr));
  #mtp<-modulTraitPvalue[]
  saveRDS(moduleTraitCor, "") #add path to save...
###### -> need to extract important correlations here !!! #########
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  saveRDS()
}

l_ply(list, makeWcorr) #make correlation heatmaps with all samples.

# Final analysis composition ####
# ID modules associated with Urbanization
# ID modules associated with City and the interaction w/ urbanization
# ID modules associated with HM abundance in soils
# Interesting Network associations... components of indicator modules...

# Test Run / Scratch space ####

library(WGCNA)
allowWGCNAThreads()
datMicroarrays=read.csv("~/Desktop/PhD/metaltable.csv")
datMicroarrays<-data.frame(datMicroarrays)
rownames(datMicroarrays)<-datMicroarrays$sample_id
datMicroarrays<-datMicroarrays[,-1]
traits<-read.csv("~/Desktop/PhD/GLUSEEN_HM_Parsed.csv")

datMicroarrays$sample_id==traits$sample_id

meanExpressionByArray=apply( datMicroarrays,1,mean, na.rm=T)
NumberMissingByArray=apply( is.na(data.frame(datMicroarrays)),1, sum)
sizeGrWindow(9, 5)
barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:50), cex.names = 0.7)

KeepArray= NumberMissingByArray<500
table(KeepArray)
datExpr=datMicroarrays[KeepArray,]
y=y[KeepArray]
ArrayName[KeepArray]

NumberMissingByGene =apply( is.na(data.frame(datMicroarrays)),2, sum)
summary(NumberMissingByGene)

variancedatExpr=as.vector(apply(as.matrix(datMicroarrays),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datMicroarrays)),2, sum) )

table(no.presentdatExpr)

KeepGenes= variancedatExpr>0 & no.presentdatExpr>=4

table(KeepGenes)
datExpr=datExpr[, KeepGenes]
GeneName=GeneName[KeepGenes]

sizeGrWindow(9, 5)
plotClusterTreeSamples(datExpr=datMicroarrays, y=traits$Landuse)


library(cluster)

#make a weighted gene co-expression network

ADJ1=abs(cor(datMicroarrays,use="p"))^6

k=as.vector(apply(ADJ1,2,sum, na.rm=T))
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

datExpr=datMicroarrays[, rank(-k,ties.method="first" )<=3600]#restrict number of genes
dissADJ=1-ADJ1

dissTOM=TOMdist(ADJ1)

pam4=pam(as.dist(dissADJ), 4)
pam5=pam(as.dist(dissADJ), 5)
pam6=pam(as.dist(dissADJ), 6)

colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))
hierADJ=hclust(as.dist(dissADJ), method="average" )
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ,colors=colorStaticADJ, dendroLabels = FALSE, hang = 0.03,
                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )

# tutorial ###


# statistical tests ####
eigengenes<-moduleEigengenes(datExpr)
summary(with(eigengenes, glm(eigengene~category*city+HM.Data))) # first make dataframe of eigengenes called eigengenes
# matR mg-rast import ####
#csv files load into phyloseq...

# plot HM data ####
plot(sam$Cd_tot~sam$Codes+sam$Cities, xlab="Location", ylab="Total Cadmium", main="Cadmium")
plot(sam$Ni_tot~sam$Codes+sam$Cities, xlab="Location", ylab="Total Nickel", main="Nickel")
plot(sam$Co_tot~sam$Codes+sam$Cities, xlab="Location", ylab="Total Cobalt", main="Cobalt")
plot(sam$Zn_tot~sam$Codes+sam$Cities, xlab="Location", ylab="Total Zinc", main="Zinc")

ggplot(data=sam, mapping=aes(x=groups,y=Cd_tot, fill=Cities))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Total Cadmium", y="Average Abundance")+coord_flip()

ggplot(data=sam, mapping=aes(x=groups,y=Co_tot, fill=Cities))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Total Cobalt", y="Average Abundance")+coord_flip()

ggplot(data=sam, mapping=aes(x=groups,y=Zn_tot, fill=Cities))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Total Zinc", y="Average Abundance")+coord_flip()

ggplot(data=sam, mapping=aes(x=groups,y=Ni_tot, fill=Cities))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Total Nickel", y="Average Abundance")+coord_flip()
#stats

summary(aov(sam$Cd_tot~sam$Codes*sam$Cities))
summary(aov(sam$Ni_tot~sam$Codes*sam$Cities))
summary(aov(sam$Co_tot~sam$Codes*sam$Cities))
summary(aov(sam$Zn_tot~sam$Codes*sam$Cities))


mgRAST<-read.csv("~/Desktop/PhD/Metagenome/testset/relative_function_pathways.csv")
mgRAST<-data.frame(mgRAST)
mgRAST$groups<-paste(mgRAST$CITY, mgRAST$TRT, sep=".")
names(mgRAST)
ggplot(data=mgRAST, mapping=aes(x=groups,y=Methanogenesis, fill=CITY))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Methanogenesis", y="Average Abundance")+coord_flip()

ggplot(data=mgRAST, mapping=aes(x=groups,y=Nitrogenfixation, fill=CITY))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Nitrogenfixation", y="Average Abundance")+coord_flip()

ggplot(data=mgRAST, mapping=aes(x=groups,y=TransportofZinc, fill=CITY))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="TransportofZinc", y="Average Abundance")+coord_flip()

ggplot(data=mgRAST, mapping=aes(x=groups,y=Arsenicresistance, fill=CITY))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Arsenicresistance", y="Average Abundance")+coord_flip()

ggplot(data=mgRAST, mapping=aes(x=groups,y=Zincresistance, fill=CITY))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Zincresistance", y="Average Abundance")+coord_flip()

ggplot(data=mgRAST, mapping=aes(x=groups,y=MultidrugResistanceEffluxPumps, fill=CITY))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="MultidrugResistanceEffluxPumps", y="Average Abundance")+coord_flip()

ggplot(data=mgRAST, mapping=aes(x=groups,y=Vibriopathogenicityisland, fill=CITY))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Vibriopathogenicityisland", y="Average Abundance")+coord_flip()

ggplot(data=mgRAST, mapping=aes(x=groups,y=Cobaltzinccadmiumresistance, fill=CITY))+geom_bar(stat="identity")+scale_fill_brewer(palette="Greens")+theme_classic()+labs(title="Cobaltzinccadmiumresistance", y="Average Abundance")+coord_flip()

summary(with(mgRAST, aov(Methanogenesis~CITY*TRT)))
summary(with(mgRAST, aov(Nitrogenfixation~CITY*TRT)))     
summary(with(mgRAST, aov(TransportofZinc~CITY*TRT)))  
summary(with(mgRAST, aov(Arsenicresistance~CITY*TRT)))  
summary(with(mgRAST, aov(Zincresistance~CITY*TRT)))  
summary(with(mgRAST, aov(MultidrugResistanceEffluxPumps~CITY*TRT)))  
summary(with(mgRAST, aov(Vibriopathogenicityisland~CITY*TRT)))  

metastorm.multidrug<-otu_table(Multidrug)[-c(1,2,4,5)]
plot(log(metastorm.multidrug)~mgRAST$MultidrugResistanceEffluxPumps)
