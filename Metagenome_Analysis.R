# function analysis pipeline 
 library(limma)
 library(Glimma)
 library(edgeR)
 library(mgrastr)
 library(phyloseq)
 library(roxygen2)
 library(vegan)
 library(igraph)
 library(ggplot2)
 library(stringr)
 library(Rmisc)
# import and format function data from MG-Rast/dada2 ####
O.key<-read.csv("~/Desktop/PhD/Metagenome/SampleKey.csv")
O.key<-data.frame(O.key)
O.key$Seq_ID<-as.character(O.key$Seq_ID) # unique sample IDs
O.key$Private_ID<-as.character(O.key$MG_ID) # unique accession codes in mgm... .3 format
 
l.sID<-as.character(O.key$MG_ID) # make list of accession numbers
names(l.sID)<-O.key$Seq_ID # name list according to unique sample IDs
sam<-read.csv("~/Desktop/PhD/Metagenome/GLUSEEN_HM2_metadata.csv")
rownames(sam)<-sam$Sample_ID2

url.15.4<-lapply(l.sID, mg.load, auth="gbdwiyCZqACamvn8fG59aTs3Z", ont="Subsystems", level="4", E="15", length="15", id="60")
dt15.l4<-downloadFunc(url.15.4, ont="Subsystems", level=4)
dt15.l4<-condense(dt15.l4, sample_key, sam)
ftax<-DFTax(url.15.4)
tax_table(dt15.l4)<-as.matrix(ftax)

# subsetting/normalizing data for analysis
saveRDS(dt15.l4, "~/Desktop/PhD/Metagenome/Function/dt15l4.RDS")
dt15.l4<-readRDS("~/Desktop/PhD/Metagenome/Function/dt15l4.RDS")

# rarefaction ####
level4OTU<-otu_table(dt15.l4)
sample_data(dt15.l4)<-sam
factor<-round(6000000*(sample_data(dt15.l4)$gperg_DNA/(mean(sample_data(dt15.l4)$gperg_DNA))))
table(sample_sums(dt15.l4)>factor)

normalizedOTU<-rrarefy(level4OTU, factor)

L4Norm<-dt15.l4
otu_table(L4Norm)<-otu_table(normalizedOTU, taxa_are_rows=F)
L4Norm<-prune_samples(sample_sums(L4Norm)==factor, L4Norm)


normseqs<-log10(sample_sums(AMOA))
log10<-log10(sample_data(L4Norm)$AmoA_copies)
plot(normseqs[normseqs!="-Inf"]~log10[normseqs!="-Inf"], col=sample_data(L4Norm)$Codes[normseqs!="-Inf"], pch=c(1,2,3,4)[as.numeric(sample_data(L4Norm)$Codes)[normseqs!="-Inf"]])
length(normseqs[normseqs!="-Inf"])
plot(normseqs[normseqs!="-Inf"]~log10[normseqs!="-Inf"])

cor(normseqs[normseqs!="-Inf"],log10[normseqs!="-Inf"])
L4Norm<-subset_samples(L4Norm, sample_names(L4Norm)!=c("Final_033", "Final_o20"))

# Norm rerun ####
level4OTU2<-otu_table(dt15.l4)
sample_data(dt15.l4)<-sam
factor<-round(6000000*(sample_data(dt15.l4)$gperg_DNA/(mean(sample_data(dt15.l4)$gperg_DNA))))
table(sample_sums(dt15.l4)>factor)

normalizedOTU2<-rrarefy(level4OTU2, factor)

L4Norm2<-dt15.l4
otu_table(L4Norm2)<-otu_table(normalizedOTU2, taxa_are_rows=F)
AMOA<-subset_taxa(L4Norm2, L4=="Ammonia monooxygenase")

normseqs<-log10(sample_sums(AMOA))
log10<-log10(sample_data(L4Norm)$AmoA_copies)

plot(sample_sums(AMOA), sample_data(AMOA)$Codes)

plot(normseqs[normseqs!="-Inf"]~log10[normseqs!="-Inf"], col=sample_data(L4Norm)$Codes[normseqs!="-Inf"], pch=c(1,2,3,4)[as.numeric(sample_data(L4Norm)$Codes)[normseqs!="-Inf"]])
length(normseqs[normseqs!="-Inf"])
plot(normseqs[normseqs!="-Inf"]~log10[normseqs!="-Inf"])

cor(normseqs[normseqs!="-Inf"],log10[normseqs!="-Inf"])
L4Norm<-subset_samples(L4Norm, sample_names(L4Norm)!=c("Final_033", "Final_o20"))

mg.amoa<-data.frame("mg"=sample_sums(AMOA), "QPCRAOA"=sample_data(AMOA)$AmoA_copies, "QPCRAOB"=sample_data(AMOA)$AmoB_copies,"Land"=sample_data(AMOA)$Codes, "City"=sample_data(AMOA)$Cities)

with(mg.amoa, cor(log10(mg), log10(QPCRAOA), use="pairwise.complete.obs"))
with(mg.amoa, plot(log10(mg)~log10(QPCRAOA), col=Land))
with(mg.amoa, plot(log10(mg)~log10(QPCRAOA), col=City, pch=c(1,2,3,4)[as.numeric(Land)]))

with(mg.amoa, plot(mg~log10(QPCRAOA), col=Land))
with(mg.amoa, plot(mg~QPCRAOB, col=Land))
with(mg.amoa, cor(mg, QPCRAOB))
with(mg.amoa, cor(mg, QPCRAO))


with(mg.amoa, plot(mg~Land))

with(mg.amoa, leveneTest(mg, Land))
plot(sample_sums(RDT.ammoa)~sample_data(RDT.ammoa)$Codes)

summary(aov(sample_sums(RDT.ammoa)~sample_data(RDT.ammoa)$Codes))
TukeyHSD(aov(sample_sums(RDT.ammoa)~sample_data(RDT.ammoa)$Codes))

####

sample_data(dt15.l4)<-sam
sample_data(dt15.l4)$samplecodes<-paste(as.character(sample_data(dt15.l4)$Cities), as.character(sample_data(dt15.l4)$Codes), sep="")
tax_table(dt15.l4)<-tax_table(dt15.l4)[,c(2,4,5,6)]

table(sample_sums(dt15.l4)<10000000)

rdt15.l4<-subset_samples(dt15.l4, sample_sums(dt15.l4)>10000000)
rdt15.l4  = transform_sample_counts(rdt15.l4, function(x) x / sum(x) ) # relative abundance of functions

sample_data(rdt15.l4)$BLK[sample_data(rdt15.l4)$Codes=="Reference"]<-1
sample_data(rdt15.l4)$BLK[sample_data(rdt15.l4)$Codes=="Remnant"]<-1
sample_data(rdt15.l4)$BLK[sample_data(rdt15.l4)$Codes=="Ruderal"]<-2
sample_data(rdt15.l4)$BLK[sample_data(rdt15.l4)$Codes=="Turf"]<-2

sample_data(rdt15.l4)$BLK2[sample_data(rdt15.l4)$Codes=="Reference"]<-1
sample_data(rdt15.l4)$BLK2[sample_data(rdt15.l4)$Codes=="Remnant"]<-1
sample_data(rdt15.l4)$BLK2[sample_data(rdt15.l4)$Codes=="Ruderal"]<-1
sample_data(rdt15.l4)$BLK2[sample_data(rdt15.l4)$Codes=="Turf"]<-2

sample_names(rdt15.l4)<-str_replace(sample_names(rdt15.l4), "Final_", "GLU") # v. important
sample_names(rdt15.l4)<-str_replace(sample_names(rdt15.l4), "South Africa", "SouthAfrica")
sample_names(rdt15.l4)<-str_replace(sample_names(rdt15.l4), "Lakti", "Lahti")

dt15.l1<-tax_glom(dt15.l4, "L1")
dt15.l3<-tax_glom(dt15.l4, "L3")

taxa_names(dt15.l3)<-tax_table(dt15.l3)[,2]

rdt15.l1<-tax_glom(rdt15.l4, "L1")
rdt15.l3<-tax_glom(rdt15.l4, "L3")


groups1<-c("BaltimoreRemnant","BaltimoreRuderal","BaltimoreTurf","BaltimoreReference","HelsinkiReference","HelsinkiRuderal","HelsinkiTurf","HelsinkiRemnant","LaktiReference","LaktiRuderal","LaktiTurf","LaktiRemnant","BudapestReference","BudapestRuderal","BudapestTurf","BudapestRemnant","South AfricaReference","South AfricaRemnant","South AfricaTurf","South AfricaRuderal")
groups2<-data.frame("Cities"=c(rep("Baltimore", 4), rep("Helsinki",4), rep("Lahti",4), rep("Budapest",4), rep("South Africa",4)), "Codes"=c("Remnant", "Ruderal", "Turf", "Reference","Reference", "Ruderal", "Turf", "Remnant", "Reference", "Ruderal", "Turf", "Remnant", "Reference","Ruderal", "Turf", "Remnant", "Reference", "Remnant", "Turf", "Ruderal"))

RDT.ammoa<-subset_taxa(L4Norm, L4=="Ammonia monooxygenase")

# subset functions ####
Nitrite.R.All<-subset_taxa(L4Norm, L4=="Copper-containing nitrite reductase (EC 1.7.2.1)"|
                             L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)"|
                             L4=="Nitric oxide reductase activation protein NorD"|
                             L4=="Nitric oxide reductase activation protein NorQ"|
                             L4=="Nitric-oxide reductase (EC 1.7.99.7), quinol-dependent"|
                             L4=="Nitric-oxide reductase subunit B (EC 1.7.99.7"|
                             L4=="Nitric-oxide reductase subunit C (EC 1.7.99.7)"|
                             L4=="Nitrous-oxide reductase (EC 1.7.99.6"|
                             L4=="Nitrous oxide reductase maturation protein NosD"|
                             L4=="Nitrous oxide reductase maturation protein NosF (ATPase)"|
                             L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)"|
                             L4=="Nitrite reductase associated c-type cytochorome NirN"|
                             L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)"|
                             L4=="Cytochrome c nitrite reductase, small subunit NrfH"|
                             L4=="Cytochrome c-type heme lyase subunit nrfE, nitrite reductase complex assembly"|
                             L4=="Cytochrome c-type heme lyase subunit nrfF, nitrite reductase complex assembly"|
                             L4=="Cytochrome c-type heme lyase subunit nrfG, nitrite reductase complex assembly"|
                             L4=="Ferredoxin--nitrite reductase (EC 1.7.7.1)"|
                             L4=="Nitrite reductase [NAD(P)H] large subunit (EC 1.7.1.4)"|
                             L4=="Nitrite reductase [NAD(P)H] small subunit (EC 1.7.1.4)"|
                             L4=="Nitrite reductase probable [NAD(P)H] subunit (EC 1.7.1.4)"|
                             L4=="Nitrous oxide reductase maturation protein NosR"|
                             L4=="Nitrous oxide reductase maturation protein, outer-membrane lipoprotein NosL"|
                             L4=="Nitrous oxide reductase maturation transmembrane protein NosY")

Nx.Reductase<-subset_taxa(L4Norm, L4=="Copper-containing nitrite reductase (EC 1.7.2.1)"|
                             L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)"|
                             L4=="Nitric-oxide reductase (EC 1.7.99.7), quinol-dependent"|
                             L4=="Nitric-oxide reductase subunit B (EC 1.7.99.7"|
                             L4=="Nitric-oxide reductase subunit C (EC 1.7.99.7)"|
                             L4=="Nitrous-oxide reductase (EC 1.7.99.6"|
                             L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)"|
                             L4=="Nitrite reductase associated c-type cytochorome NirN"|
                             L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)"|
                             L4=="Cytochrome c nitrite reductase, small subunit NrfH"|
                             L4=="Ferredoxin--nitrite reductase (EC 1.7.7.1)"|
                             L4=="Nitrite reductase [NAD(P)H] large subunit (EC 1.7.1.4)"|
                             L4=="Nitrite reductase [NAD(P)H] small subunit (EC 1.7.1.4)"|
                             L4=="Nitrite reductase probable [NAD(P)H] subunit (EC 1.7.1.4)")

Antibiotics<-subset_taxa(L4Norm, L4=="Multiple antibiotic resistance protein MarA"|
                           L4=="Multiple antibiotic resistance protein MarC"|
                           L4=="Multiple antibiotic resistance protein MarR"|
                           L4=="Multidrug resistance transporter, Bcr/CflA family"|
                           L4=="Multidrug efflux transporter MexF"|
                           L4=="Membrane fusion protein of RND family multidrug efflux pump"|
                           L4=="Multi antimicrobial extrusion protein (Na(+)/drug antiporter), MATE family of MDR efflux pumps"|
                           L4=="Multidrug efflux RND transporter MexD"|
                           L4=="Multidrug-efflux transporter, major facilitator superfamily (MFS) (TC 2.A.1"|
                           L4=="Multidrug resistance efflux pump PmrA"|
                           L4=="Multidrug transporter MdtB"|
                           L4=="Multidrug transporter MdtC"|
                           L4=="Multidrug transporter MdtD"|
                           L4=="MFS family multidrug efflux protein, similarity to bicyclomycin resistance protein Bcr"|
                           L4=="ABC transporter multidrug efflux pump, fused ATP-binding domains"|
                           L4=="ABC-type multidrug transport system, ATPase component"|
                           L4=="ABC-type multidrug transport system, permease component")

NickelTransport<-subset_taxa(L4Norm, L4=="Additional substrate-specific component NikN of nickel ECF transporter"|
                               L4=="ATPase component NikO of energizing module of nickel ECF transporter"|
                               L4=="Substrate-specific component NikM of nickel ECF transporter"|
                               L4=="Additional substrate-specific component NikN of nickel ECF transporter"|
                               L4=="ATPase component NikO of energizing module of nickel ECF transporter"|
                               L4=="Substrate-specific component NikM of nickel ECF transporter"|
                               L4=="Transmembrane component NikQ of energizing module of nickel ECF transporter"|
                               L4=="Nickel ABC transporter, periplasmic nickel-binding protein nikA2 (TC 3.A.1.5.3)"|
                               L4=="Nickel ABC transporter, periplasmic nickel-binding protein NikA (TC 3.A.1.5.3)"|
                               L4=="Nickel transport ATP-binding protein nikD2 (TC 3.A.1.5.3)"|
                               L4=="Nickel transport ATP-binding protein NikD (TC 3.A.1.5.3)"|
                               L4=="Nickel transport ATP-binding protein nikE2 (TC 3.A.1.5.3)"|
                               L4=="Nickel transport ATP-binding protein NikE (TC 3.A.1.5.3)"|
                               L4=="Nickel transport system permease protein nikB2 (TC 3.A.1.5.3"|
                               L4=="Nickel transport system permease protein NikB (TC 3.A.1.5.3)"|
                               L4=="Nickel transport system permease protein nikC2 (TC 3.A.1.5.3)"|
                               L4=="Nickel transport system permease protein NikC (TC 3.A.1.5.3)"|
                               L4=="HoxN/HupN/NixA family nickel/cobalt transporter")

CobaltTransport<-subset_taxa(L4Norm, L4=="ATPase component CbiO of energizing module of cobalt ECF transporter"|
                               L4=="Substrate-specific component CbiM of cobalt ECF transporter"|
                               L4=="Transmembrane component CbiQ of energizing module of cobalt ECF transporter"|
                               L4=="Additional substrate-specific component CbiN of cobalt ECF transporter"|
                               L4=="ATPase component CbiO of energizing module of cobalt ECF transporter"|
                               L4=="HoxN/HupN/NixA family cobalt transporter"|
                               L4=="HoxN/HupN/NixA family nickel/cobalt transporter"|
                               L4=="Predicted cobalt transporter CbtA"|
                               L4=="Substrate-specific component CbiM of cobalt ECF transporter"|
                               L4=="Transmembrane component CbiQ of energizing module of cobalt ECF transporter"|
                               L4=="Magnesium and cobalt efflux protein CorC"|
                               L4=="Magnesium and cobalt transport protein CorA"|
                               L4=="Magnesium and cobalt efflux protein CorC"|
                               L4=="Cobalt/zinc/cadmium efflux RND transporter, membrane fusion protein, CzcB family"|
                               L4=="Cobalt-zinc-cadmium resistance protein"|
                               L4=="Cobalt-zinc-cadmium resistance protein CzcA"|
                               L4=="Cobalt-zinc-cadmium resistance protein CzcD"|
                               L4=="Additional substrate-specific component CbiN of cobalt ECF transporter"|
                               L4=="HoxN/HupN/NixA family cobalt transporter"|
                               L4=="Predicted cobalt transporter CbtA"|
                               L4=="Substrate-specific component CbiM of cobalt ECF transporter"|
                               L4=="Transmembrane component CbiQ of energizing module of cobalt ECF transporter"|
                               L4=="Magnesium and cobalt efflux protein CorC"|
                               L4=="Additional substrate-specific component CbiN of cobalt ECF transporter"|
                               L4=="Probable Co/Zn/Cd efflux system membrane fusion protein2637")

ArsenicResistance<-subset_taxa(L4Norm, L4=="Arsenical-resistance protein ACR3"|L4=="Arsenic efflux pump protein"|L4=="Arsenic resistance protein ArsH")

CopperRegulation<-subset_taxa(L4Norm, L4=="Copper-translocating P-type ATPase (EC 3.6.3.4)"|
                                 L4=="Copper-translocating P-type ATPase (EC 3.6.3.4)"|
                                 L4=="Copper homeostasis protein CutE"|
                                 L4=="Copper homeostasis protein CutF precursor"|
                                 L4=="Cytoplasmic copper homeostasis protein cutC"|
                                 L4=="Membrane protein, suppressor for copper-sensitivity ScsB"|
                                 L4=="Membrane protein, suppressor for copper-sensitivity ScsD"|
                                 L4=="Secreted protein, suppressor for copper-sensitivity ScsC"|
                                 L4=="Suppression of copper sensitivity: putative copper binding protein ScsA"|
                                 L4=="Copper resistance protein B"|
                                 L4=="Copper resistance protein D"|
                                 L4=="Copper tolerance protein"|
                                 L4=="Copper-translocating P-type ATPase (EC 3.6.3.4)"|
                                 L4=="Copper homeostasis protein CutE")

CadmiumResistance<-subset_taxa(L4Norm, L4=="Cadmium efflux system accessory protein"|
                                 L4=="Cadmium-transporting ATPase (EC 3.6.3.3)"|
                                 L4=="Cobalt/zinc/cadmium efflux RND transporter, membrane fusion protein, CzcB family"|
                                 L4=="Cobalt-zinc-cadmium resistance protein"|
                                 L4=="Cobalt-zinc-cadmium resistance protein CzcA"|
                                 L4=="Cobalt-zinc-cadmium resistance protein CzcD"|
                                 L4=="Probable cadmium-transporting ATPase (EC 3.6.3.3)"|
                                 L4=="Probable Co/Zn/Cd efflux system membrane fusion protein"|
                                 L4=="Cadmium efflux system accessory protein"
                                 )


CZC<-subset_taxa(L4Norm, L4=="Cobalt-zinc-cadmium resistance protein"|
                   L4=="Cobalt-zinc-cadmium resistance protein CzcA"|
                   L4=="Cobalt-zinc-cadmium resistance protein CzcD"|
                   L4=="Cobalt/zinc/cadmium efflux RND transporter, membrane fusion protein, CzcB family"|
                   L4=="Probable Co/Zn/Cd efflux system membrane fusion protein"|
                   L4=="Putative metal chaperone, involved in Zn homeostasis, GTPase of COG0523 family")

Nitrogenase<-subset_taxa(L4Norm, L4=="Nitrogenase FeMo-cofactor scaffold and assembly protein NifE"|
                           L4=="Nitrogenase FeMo-cofactor scaffold and assembly protein NifN"|
                           L4=="Nitrogenase FeMo-cofactor synthesis FeS core scaffold and assembly protein NifB"|
                           L4=="Nitrogenase (molybdenum-iron) alpha chain (EC 1.18.6.1)"|
                           L4=="Nitrogenase (molybdenum-iron) beta chain (EC 1.18.6.1)"|
                           L4=="Nitrogenase (molybdenum-iron) reductase and maturation protein NifH"|
                           L4=="Nitrogenase (molybdenum-iron)-specific transcriptional regulator NifA")


ZincResistance<-subset_taxa(L4Norm, L4=="Cobalt-zinc-cadmium resistance protein"|
                              L4=="Cobalt-zinc-cadmium resistance protein CzcA"|
                              L4=="Cobalt-zinc-cadmium resistance protein CzcD"|
                              L4=="Zinc transporter ZitB"|
                              L4=="Zinc ABC transporter, ATP-binding protein ZnuC"|
                              L4=="Zinc ABC transporter, inner membrane permease protein ZnuB"|
                              L4=="Zinc ABC transporter, periplasmic-binding protein ZnuA"|
                              L4=="Zinc ABC transporter, ATP-binding protein ZnuC"|
                              L4=="Zinc ABC transporter, inner membrane permease protein ZnuB"|
                              L4=="Zinc ABC transporter, periplasmic-binding protein ZnuA"|
                              L4=="Zinc transporter, ZIP family"|
                              L4=="Zinc resistance-associated protein"|
                              L4=="Cobalt/zinc/cadmium efflux RND transporter, membrane fusion protein, CzcB family"|
                              L4=="Probable Co/Zn/Cd efflux system membrane fusion protein"|
                              L4=="Putative metal chaperone, involved in Zn homeostasis, GTPase of COG0523 family"
                   )


#Baltimore.F<-subset_samples(dt15.l3, Cities=="Baltimore")
#Helsinki.F<-subset_samples(dt15.l3, Cities=="Helsinki")
#Lahti.F<-subset_samples(dt15.l3, Cities=="Lahti")
#Budapest.F<-subset_samples(dt15.l3, Cities=="Budapest")
#SouthAfrica.F<-subset_samples(dt15.l3, Cities=="SouthAfrica")

# import bacteria
baclib<-readRDS("~/Desktop/PhD/Metagenome/Bac16s_rarefied.RDS") 
sample_data(baclib)$samplecodes<-paste(as.character(sample_data(baclib)$Cities), as.character(sample_data(baclib)$Codes), sep="")

baclib2<-readRDS("~/Desktop/PhD/Metagenome/Bac16s_raw.RDS") 
sample_data(baclib2)$samplecodes<-paste(as.character(sample_data(baclib2)$Cities), as.character(sample_data(baclib2)$Codes), sep="")
#sample_names(baclib)
#samdf<-read.csv("~/Desktop/PhD/Metagenome/GLUSEEN_HM2_metadata.csv")
#rownames(samdf)<-samdf$Sample_ID

baclib<-subset_samples(baclib, sample_sums(baclib)>100)

bac<-subset_taxa(baclib, Kingdom=="Bacteria")
arc<-subset_taxa(baclib, Kingdom=="Archaea")
AOB<-subset_taxa(baclib, Genus=="Nitrospira"|Genus=="Nitrospina"|Genus=="Nitrosomonas"|Genus=="Nitrobacter"|Genus=="Nitrococcus"|Genus=="Nitrosococcus"|Genus=="Nitrosovibreo")#|Family=="Nitrosomonadaceae")
AOA<-subset_taxa(baclib, Phylum=="Thaumarchaeota")
AOB<-subset_taxa(baclib, Genus=="Nitrosomonas"|Genus=="Nitrococcus"|Genus=="Nitrospira")

#
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Nitrogen Metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Miscellaneous"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Virulence, Disease and Defense"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Amino Acids and Derivatives"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Clustering-based subsystems"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="DNA Metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="RNA Metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Potassium metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Motility and Chemotaxis"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Potassium metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Iron acquisition and metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Sulfur Metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Cofactors, Vitamins, Prosthetic Groups, Pigments"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Stress Response"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Cofactors, Vitamins, Prosthetic Groups, Pigments"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Phosphorus Metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Fatty Acids, Lipids, and Isoprenoids"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Photosynthesis"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Nucleosides and Nucleotides"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Secondary Metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Carbohydrates"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Dormancy and Sporulation"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Phages, Prophages, Transposable elements, Plasmids"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Respiration"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Protein Metabolism"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Cell Division and Cell Cycle"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Membrane Transport"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Cell Wall and Capsule"))~sample_data(dt15.l1)$Codes)
plot(sample_sums(subset_taxa(rdt15.l1, L1=="Regulation and Cell signaling"))~sample_data(dt15.l1)$Codes)

clss<-subset_taxa(rdt15.l4, L1=="Clustering-based subsystems")
tax_table(clss)

# anova for soil param ####

env.param<-as.data.frame(sample_data(baclib))

summary(aov(env.param$Cd_tot~env.param$Codes*env.param$Cities))
summary(aov(env.param$Co_tot~env.param$Codes*env.param$Cities))
summary(aov(env.param$Ni_tot~env.param$Codes*env.param$Cities))
summary(aov(env.param$Zn_tot~env.param$Codes*env.param$Cities))

leveneTest(env.param$Cd_tot, group = env.param$Codes)
leveneTest(env.param$Co_tot, group = env.param$Codes)
leveneTest(env.param$Ni_tot, group = env.param$Codes)
leveneTest(env.param$Zn_tot, group = env.param$Codes)

summary(aov(env.param$Cd_avail~env.param$Codes*env.param$Cities))
summary(aov(env.param$Co_avail~env.param$Codes*env.param$Cities))
summary(aov(env.param$Ni_avail~env.param$Codes*env.param$Cities))
summary(aov(env.param$Zn_avail~env.param$Codes*env.param$Cities))

leveneTest(env.param$Cd_avail, group = env.param$Codes)
leveneTest(env.param$Co_avail, group = env.param$Codes)
leveneTest(env.param$Ni_avail, group = env.param$Codes)
leveneTest(env.param$Zn_avail, group = env.param$Codes)

TukeyHSD(aov(env.param$Ni_avail~env.param$Codes*env.param$Cities))
TukeyHSD(aov(env.param$Zn_avail~env.param$Codes*env.param$Cities))

Cd.sum<-summarySE(env.param, measurevar="Cd_tot", groupvars = "Cities", na.rm=T)
Co.sum<-summarySE(env.param, measurevar="Co_tot", groupvars = "Cities", na.rm=T)
Ni.sum<-summarySE(env.param, measurevar="Ni_tot", groupvars = "Cities", na.rm=T)
Zn.sum<-summarySE(env.param, measurevar="Zn_tot", groupvars = "Cities", na.rm=T)
Cdavail.sum<-summarySE(env.param, measurevar="Cd_avail", groupvars = "Codes", na.rm=T)
Coavail.sum<-summarySE(env.param, measurevar="Co_avail", groupvars = "Codes", na.rm=T)
Niavail.sum<-summarySE(env.param, measurevar="Ni_avail", groupvars = "Codes", na.rm=T)
Znavail.sum<-summarySE(env.param, measurevar="Zn_avail", groupvars = "Codes", na.rm=T)

AOA.sum<-summarySE(env.param, measurevar="AmoA_copies", groupvars = "Codes", na.rm=T)
AOB.sum<-summarySE(env.param, measurevar="AmoB_copies", groupvars = "Codes", na.rm=T)

colnames(AOA.sum)<-c("Codes","N", "Copies", "sd", "se", "ci")
colnames(AOB.sum)<-c("Codes", "N", "Copies", "sd", "se", "ci")
comb<-rbind(AOA.sum, AOB.sum)
comb$group<-c("A","A","A","A", "B","B","B","B")

ggplot(comb, aes(x=Codes, y=Copies)) +
  geom_errorbar(aes(ymin=Copies-se, ymax=Copies+se, color=group), 
    width=.1, size=1) +
  geom_point(aes(color=group))+
  scale_colour_grey()+
  theme_bw()

ggplot(AOB.sum, aes(x=Codes, y=Copies)) +
  geom_errorbar(aes(ymin=Copies-se, ymax=Copies+se), width=.1, color="black", size=1) +
  geom_point(aes())+theme_bw()

ggplot(Cd.sum, aes(x=Cities, y=Cd_tot)) +
  geom_errorbar(aes(ymin=Cd_tot-se, ymax=Cd_tot+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Co.sum, aes(x=Cities, y=Co_tot)) +
  geom_errorbar(aes(ymin=Co_tot-se, ymax=Co_tot+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()


ggplot(Ni.sum, aes(x=Cities, y=Ni_tot)) +
  geom_errorbar(aes(ymin=Ni_tot-se, ymax=Ni_tot+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Zn.sum, aes(x=Cities, y=Zn_tot)) +
  geom_errorbar(aes(ymin=Zn_tot-se, ymax=Zn_tot+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()


# run limma-voom 15 ####
# L4 City ####
dt.L4<-tax_glom(dt15.l4, "L4")
tax<-NULL
tax<-paste(tax_table(dt15.l4)[,3], c(1:12363), sep="")
f.L4<-dt.L4
taxa_names(f.L4)<-tax
counts<-t(as.data.frame(as.matrix(otu_table(f.L4))))
geneid <- rownames(counts)
categories<-data.frame("Cities"=as.factor(sample_data(f.L4)$Cities), "Landuse"=as.factor(sample_data(f.L4)$Codes))
L4_City.design<-model.matrix(~0+Cities, categories)
#colnames(design) <- gsub("group3", "", colnames(design))

dge.L4_city <- DGEList(counts=counts)


keep5 <- filterByExpr(dge.L4_city,L4_City.design)
dge.L4_city <- dge.L4_city[keep5,,keep.lib.sizes=FALSE] #error
dge.L4_city <- calcNormFactors(dge.L4_city)

dge.L4_city$genes<-geneid
lcpm<-cpm(dge.L4_city, log=TRUE)


v.L4_city <- voom(dge.L4_city, L4_City.design, plot=TRUE)
#dup<-duplicateCorrelation(v, design3, block=categories$Cities)

fitV.L4_city <- lmFit(v.L4_city, L4_City.design)

fitV.L4_city <- contrasts.fit(fitV.L4_city, contrasts=contr.matrix)
fitV.L4_city <- eBayes(fitV.L4_city, trend=TRUE)

length(fitV.L4_city$F.p.value)
length(fitV.L4_city$F.p.value[fitV.L4_city$F.p.value<0.05])

# L4 Landuse ####
dt.L4<-tax_glom(dt15.l4, "L4")
tax<-NULL
tax<-paste(tax_table(dt15.l4)[,3], c(1:12363), sep="")
f.L4<-dt.L4

taxa_names(f.L4)<-tax
counts<-t(as.data.frame(as.matrix(otu_table(f.L4))))
geneid <- rownames(counts)

categories<-data.frame("Cities"=as.factor(sample_data(f.L4)$Cities), "Landuse"=as.factor(sample_data(f.L4)$Codes))


L4.design<-model.matrix(~0+Landuse, categories)
#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(L4.design))
contr.matrix


dge.L4 <- DGEList(counts=counts)


keep3 <- filterByExpr(dge.L4,L4.design)
dge.L4 <- dge.L4[keep3,,keep.lib.sizes=FALSE] #error

dge.L4 <- calcNormFactors(dge.L4)

dge.L4$genes<-geneid
lcpm<-cpm(dge.L4, log=TRUE)
col.group <- categories$Cities
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
sym<-as.factor(categories$Landuse)
plotMDS(lcpm, col=col.group, pch=as.numeric(sym))

v.L4 <- voom(dge.L4, L4.design, plot=TRUE)
#dup<-duplicateCorrelation(v, design3, block=categories$Cities)

fitV.L4 <- lmFit(v.L4, L4.design) #what?

fitV.L4 <- contrasts.fit(fitV.L4, contrasts=contr.matrix)
fitV.L4 <- eBayes(fitV.L4, trend=TRUE)
topgenes<-topTable(fitV.L4, number=100, coef=3) #more significant results!!
siggenes<-topgenes[topgenes$adj.P.Val<0.05,]

summary(decideTests(fitV.L4))
L4sigsummary<-decideTests(fitV.L4)
L4.TurfP <- which(L4sigsummary[,1]==1)
L4.TurfN <- which(L4sigsummary[,1]==-1)
L4.RemP <- which(L4sigsummary[,2]==1)
L4.RemN <- which(L4sigsummary[,2]==-1)
L4.RudP <- which(L4sigsummary[,3]==1)
L4.RudN <- which(L4sigsummary[,3]==-1)

names(L4.TurfP)
names(L4.TurfN)
names(L4.RemP)
names(L4.RemN)
names(L4.RudP)
names(L4.RudN)

write.csv(names(which(L4sigsummary[,2]==0)), "~/Desktop/PhD/Metagenome/NonsigGenes.csv")

# on L3 ####
counts<-t(as.data.frame(as.matrix(otu_table(dt15.l3))))
geneid <- rownames(counts)

categories<-data.frame("Cities"=as.factor(sample_data(dt15.l3)$Cities), "Landuse"=as.factor(sample_data(dt15.l3)$Codes))


design3<-model.matrix(~0+Landuse, categories)
#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(design3))
contr.matrix


dge3 <- DGEList(counts=counts)


keep3 <- filterByExpr(dge3, design3)
dge3 <- dge3[keep3,,keep.lib.sizes=FALSE] #error

dge3 <- calcNormFactors(dge3)

dge3$genes<-geneid
lcpm<-cpm(dge3, log=TRUE)
col.group <- categories$Cities
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
sym<-as.factor(categories$Landuse)
plotMDS(lcpm, col=col.group, pch=as.numeric(sym))

v <- voom(dge3, design3, plot=TRUE)
#dup<-duplicateCorrelation(v, design3, block=categories$Cities)

fitV3 <- lmFit(v, design3) #what?

fitV3 <- contrasts.fit(fitV3, contrasts=contr.matrix)
fitV3 <- eBayes(fitV3, trend=TRUE)
topgenes<-topTable(fitV3, number=100, coef=3) #more significant results!!
siggenes<-topgenes[topgenes$adj.P.Val<0.05,]

summary(decideTests(fitV3))
sigsummary<-decideTests(fitV3)
de.TurfP <- which(sigsummary[,1]==1)
de.TurfN <- which(sigsummary[,1]==-1)
de.RemP <- which(sigsummary[,2]==1)
de.RemN <- which(sigsummary[,2]==-1)
de.RudP <- which(sigsummary[,3]==1)
de.RudN <- which(sigsummary[,3]==-1)

de.TurfP
de.TurfN
de.RemP
de.RemN
de.RudP
de.RudN

#plots that don't work...
tfit <- treat(fitV3, lfc=1)
dt <- decideTests(tfit)
basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)


glMDPlot(fitV3, coef=1, status=sigsummary, main=colnames(fitV3)[1],
         side.main="ENTREZID", counts=logCPM3, groups=categories$Landuse, launch=FALSE)

group<-as.factor(categories$Landuse)

order<-sample_data(dt15.l3)$Codes
names(order)<-sample_data(dt15.l3)$Sample_ID2
orderL<-sort(order, decreasing = F)
orderN<-names(orderL)

basal.vs.lp.topgenes <- basal.vs.lp$ID[1:30]
i <- which(v$genes %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
hm<-heatmap.2(lcpm[i,orderN], Colv=FALSE, scale="row", dendrogram = "none",
          labRow=v$genes[i], labCol=orderL, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10))

basal.vs.lp.topgenes <- rownames(topgenes)[1:30]
i <- which(v$genes %in% basal.vs.lp.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
hm<-heatmap.2(lcpm[i,orderN], Colv=FALSE, scale="row", dendrogram = "none",
              labRow=v$genes[i], labCol=orderL, 
              col=mycol, trace="none", density.info="none", 
              margin=c(8,6), lhei=c(2,10))


# Baltimore voom ####
Baltimorecounts<-t(as.data.frame(as.matrix(otu_table(Baltimore.F))))
geneid <- rownames(Baltimorecounts)

BaltimoreLanuse<-data.frame("Landuse"=sample_data(Baltimore.F)$Codes)


Baltimoredesign3<-model.matrix(~0+Landuse, BaltimoreLanuse)
#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(Baltimoredesign3))
contr.matrix


B1dge3 <- DGEList(counts=Baltimorecounts)


B1keep3 <- filterByExpr(B1dge3, Baltimoredesign3)
B1dge3 <- B1dge3[B1keep3,,keep.lib.sizes=FALSE] #error

B1dge3 <- calcNormFactors(B1dge3)

B1dge3$genes<-geneid
B1lcpm<-cpm(B1dge3, log=TRUE)

B1v <- voom(B1dge3, Baltimoredesign3, plot=TRUE)
#B1dup<-duplicateCorrelation(B1v, Baltimoredesign3)

B1fitV3 <- lmFit(B1v, Baltimoredesign3) #what?

B1fitV3 <- contrasts.fit(B1fitV3, contrasts=contr.matrix)
B1fitV3 <- eBayes(B1fitV3, trend=TRUE)
B1topgenes<-topTable(B1fitV3, number=100, coef=3) #more significant results!!

summary(decideTests(B1fitV3))
B1sigsummary<-decideTests(B1fitV3)
B1de.TurfP <- which(B1sigsummary[,1]==1)
B1de.TurfN <- which(B1sigsummary[,1]==-1)
B1de.RemP <- which(B1sigsummary[,1]==1)
B1de.RemN <- which(B1sigsummary[,1]==-1)
B1de.RudP <- which(B1sigsummary[,3]==1)
B1de.RudN <- which(B1sigsummary[,3]==-1)

B1de.TurfP
B1de.TurfN
B1de.RemP
B1de.RemN
B1de.RudP
B1de.RudN

# Helsinki
Helsinkicounts<-t(as.data.frame(as.matrix(otu_table(Helsinki.F))))
geneid <- rownames(Helsinkicounts)

HelsinkiLanuse<-data.frame("Landuse"=sample_data(Helsinki.F)$Codes)


Helsinkidesign3<-model.matrix(~0+Landuse, HelsinkiLanuse)
#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(Helsinkidesign3))
contr.matrix


Hdge3 <- DGEList(counts=Helsinkicounts)


Hkeep3 <- filterByExpr(Hdge3, Helsinkidesign3)
Hdge3 <- Hdge3[Hkeep3,,keep.lib.sizes=FALSE] #error

Hdge3 <- calcNormFactors(Hdge3)

Hdge3$genes<-geneid
Hlcpm<-cpm(Hdge3, log=TRUE)

Hv <- voom(Hdge3, Helsinkidesign3, plot=TRUE)
#Hdup<-duplicateCorrelation(Hv, Helsinkidesign3)

HfitV3 <- lmFit(Hv, Helsinkidesign3) #what?

HfitV3 <- contrasts.fit(HfitV3, contrasts=contr.matrix)
HfitV3 <- eBayes(HfitV3, trend=TRUE)
Htopgenes<-topTable(HfitV3, number=100, coef=3) #more significant results!!

summary(decideTests(HfitV3))
Hsigsummary<-decideTests(HfitV3)
Hde.TurfP <- which(Hsigsummary[,1]==1)
Hde.TurfN <- which(Hsigsummary[,1]==-1)
Hde.RemfP <- which(Hsigsummary[,1]==1)
Hde.RemfN <- which(Hsigsummary[,1]==-1)
Hde.RudP <- which(Hsigsummary[,3]==1)
Hde.RudN <- which(Hsigsummary[,3]==-1)

Hde.TurfP 
Hde.TurfN 
Hde.RemfP
Hde.RemfN
Hde.RudP
Hde.RudN

# Lahti
Lahticounts<-t(as.data.frame(as.matrix(otu_table(Lahti.F))))
geneid <- rownames(Lahticounts)

LahtiLanuse<-data.frame("Landuse"=sample_data(Lahti.F)$Codes)


Lahtidesign3<-model.matrix(~0+Landuse, LahtiLanuse)
#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(Lahtidesign3))
contr.matrix


Ldge3 <- DGEList(counts=Lahticounts)


Lkeep3 <- filterByExpr(Ldge3, Lahtidesign3)
Ldge3 <- Ldge3[Lkeep3,,keep.lib.sizes=FALSE] #error

Ldge3 <- calcNormFactors(Ldge3)

Ldge3$genes<-geneid
Llcpm<-cpm(Ldge3, log=TRUE)

Lv <- voom(Ldge3, Lahtidesign3, plot=TRUE)
#Ldup<-duplicateCorrelation(Lv, Lahtidesign3)

LfitV3 <- lmFit(Lv, Lahtidesign3) #what?

LfitV3 <- contrasts.fit(LfitV3, contrasts=contr.matrix)
LfitV3 <- eBayes(LfitV3, trend=TRUE)
Ltopgenes<-topTable(LfitV3, number=100, coef=3) #more significant results!!

summary(decideTests(LfitV3))
Lsigsummary<-decideTests(LfitV3)
Lde.TurfP <- which(Lsigsummary[,1]==1)
Lde.TurfN <- which(Lsigsummary[,1]==-1)

Lde.RudP <- which(Lsigsummary[,3]==1)
Lde.RudN <- which(Lsigsummary[,3]==-1)

# Budapest
Budapestcounts<-t(as.data.frame(as.matrix(otu_table(Budapest.F))))
geneid <- rownames(Budapestcounts)

BudapestLanuse<-data.frame("Landuse"=sample_data(Budapest.F)$Codes)


Budapestdesign3<-model.matrix(~0+Landuse, BudapestLanuse)
#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(Budapestdesign3))
contr.matrix


B2dge3 <- DGEList(counts=Budapestcounts)


B2keep3 <- filterByExpr(B2dge3, Budapestdesign3)
B2dge3 <- B2dge3[B2keep3,,keep.lib.sizes=FALSE] #error

B2dge3 <- calcNormFactors(B2dge3)

B2dge3$genes<-geneid
B2lcpm<-cpm(B2dge3, log=TRUE)

B2v <- voom(B2dge3, Budapestdesign3, plot=TRUE)
B2dup<-duplicateCorrelation(B2v, Budapestdesign3)

B2fitV3 <- lmFit(B2v, Budapestdesign3) #what?

B2fitV3 <- contrasts.fit(B2fitV3, contrasts=contr.matrix)
B2fitV3 <- eBayes(B2fitV3, trend=TRUE)
B2topgenes<-topTable(B2fitV3, number=100, coef=3) #more significant results!!

summary(decideTests(B2fitV3))
B2sigsummary<-decideTests(B2fitV3)
B2de.TurfP <- which(B2sigsummary[,1]==1)
B2de.TurfN <- which(B2sigsummary[,1]==-1)

B2de.RudP <- which(B2sigsummary[,3]==1)
B2de.RudN <- which(B2sigsummary[,3]==-1)

# South Africa

SAcounts<-t(as.data.frame(as.matrix(otu_table(SouthAfrica.F))))
geneid <- rownames(SAcounts)

SALanuse<-data.frame("Landuse"=sample_data(SouthAfrica.F)$Codes)


SAdesign3<-model.matrix(~0+Landuse, SALanuse)
#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(SAdesign3))
contr.matrix


SAdge3 <- DGEList(counts=SAcounts)


SAkeep3 <- filterByExpr(SAdge3, SAdesign3)
SAdge3 <- SAdge3[SAkeep3,,keep.lib.sizes=FALSE] #error

SAdge3 <- calcNormFactors(SAdge3)

SAdge3$genes<-geneid
SAlcpm<-cpm(SAdge3, log=TRUE)

SAv <- voom(SAdge3, SAdesign3, plot=TRUE)
#B1dup<-duplicateCorrelation(SAv, SAdesign3)

SAfitV3 <- lmFit(SAv, SAdesign3) #what?

SAfitV3 <- contrasts.fit(SAfitV3, contrasts=contr.matrix)
SAfitV3 <- eBayes(SAfitV3, trend=TRUE)
SAtopgenes<-topTable(SAfitV3, number=100, coef=3) #more significant results!!

summary(decideTests(SAfitV3))
SAsigsummary<-decideTests(SAfitV3)
SAde.TurfP <- which(SAsigsummary[,1]==1)
SAde.TurfN <- which(SAsigsummary[,1]==-1)

SAde.RudP <- which(SAsigsummary[,3]==1)
SAde.RudN <- which(SAsigsummary[,3]==-1)


# not sure this is worth the time...
sigtable<-matrix(nrow=6,ncol=6)
colnames(sigtable)<-c("Global", "Baltimore", "Budapest", "Helsinki", "Lahti", "Potchefstroom")
rownames(sigtable)<-c("Remnant Up", "Remnant Down", "Turf Up", "Turf Down", "Ruderal Up", "Ruderal Down")
sigtable[,1]<-c(list(names(de.RemP)),list(names(de.RemN)), list(names(de.TurfP)),list(names(de.TurfN)),list(names(de.RudP)),list(names(de.RudN)))


# scratch space

summary(anova(glm()))

rdt15.l3  = transform_sample_counts(dt15.l3, function(x) x / sum(x) )

# land use model on level 3
taxa_names(DT.Nitrogen3)<-as.matrix(tax_table(DT.Nitrogen3))[,2]

Ncounts3<-t(as.data.frame(as.matrix(otu_table(DT.Nitrogen3))))

categories<-data.frame("Cities"=as.factor(sample_data(DT.Nitrogen3)$Cities), "Landuse"=as.factor(sample_data(DT.Nitrogen3)$Codes))
design<-model.matrix(~Landuse, categories)
Ndge <- DGEList(counts=Ncounts3)
Nkeep <- filterByExpr(Ndge, design)
Ndge <- Ndge[Nkeep,,keep.lib.sizes=FALSE] #error
Ndge <- calcNormFactors(Ndge)

NlogCPM <- cpm(Ndge, log=TRUE, prior.count=3)

Nfit <- lmFit(NlogCPM, design)
Nfit <- eBayes(Nfit, trend=TRUE)
topTable(Nfit, number=100, coef=ncol(design))
saveRDS(fit,"~/Desktop/PhD/Metagenome/Function/limmaVoomFit.RDS")
fit.nested<-readRDS("~/Desktop/PhD/Metagenome/Function/limmaVoomFit.RDS")

# city-wise comparisons

design2<-model.matrix(~Cities, categories)
dge2 <- DGEList(counts=counts)
keep2 <- filterByExpr(dge2, design2)
dge2 <- dge2[keep2,,keep.lib.sizes=FALSE]
dge2 <- calcNormFactors(dge2)

logCPM2 <- cpm(dge2, log=TRUE, prior.count=3)

fit2 <- lmFit(logCPM2, design2)
fit2 <- eBayes(fit2, trend=TRUE)
topTable(fit2, coef=ncol(design2))
saveRDS(fit,"~/Desktop/PhD/Metagenome/Function/limmaVoomFit_citiesonly.RDS")
fit.cities<-readRDS("~/Desktop/PhD/Metagenome/Function/limmaVoomFit_citiesonly.RDS")

# landuse only

design3<-model.matrix(~Landuse, categories)
dge3 <- DGEList(counts=counts)
keep3 <- filterByExpr(dge3, design3)
dge3 <- dge3[keep3,,keep.lib.sizes=FALSE] #error
dge3 <- calcNormFactors(dge3)

logCPM3 <- cpm(dge3, log=TRUE, prior.count=3)

fit3 <- lmFit(logCPM3, design3)
fit3 <- eBayes(fit3, trend=TRUE)
topTable(fit3, number=100, coef=ncol(design3))
saveRDS(fit3,"~/Desktop/PhD/Metagenome/Function/limmaVoomFit_landuseonly.RDS")
fit.landuse<-readRDS("~/Desktop/PhD/Metagenome/Function/limmaVoomFitlanduseonly.RDS")
#<-data.frame(fit$stdev.unscaled)
coeffs<-data.frame(fit$coefficients)
rownames(coeffs)
coeffs$CitiesLakti.LanduseTurf


# two way anova

design4<-model.matrix(~Landuse*Cities, categories)
dge4 <- DGEList(counts=counts)
keep4 <- filterByExpr(dge4, design4)
dge4 <- dge4[keep4,,keep.lib.sizes=FALSE] #error
dge4 <- calcNormFactors(dge4)

logCPM4 <- cpm(dge4, log=TRUE, prior.count=3)

fit4 <- lmFit(logCPM4, design4)
fit4 <- eBayes(fit4, trend=TRUE)
topTable(fit4, number=100, coef=ncol(design4))
saveRDS(fit4,"~/Desktop/PhD/Metagenome/Function/limmaVoomFit_twowayfactorial.RDS")
fit.landuse<-readRDS("~/Desktop/PhD/Metagenome/Function/limmaVoomFitlanduseonly.RDS")

# limma trend ####

?subset_taxa
# correct statistical model
counts<-t(as.data.frame(as.matrix(otu_table(dt20))))

categories<-data.frame("Cities"=as.factor(sample_data(dt20)$Cities), "Landuse"=as.factor(sample_data(dt20)$Codes))
design<-model.matrix(~Cities/Landuse, categories)
dge <- DGEList(counts=counts)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE] #error
dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log=TRUE, prior.count=3)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, number=100, coef=ncol(design))
saveRDS(fit,"~/Desktop/PhD/Metagenome/Function/limmaVoomFit.RDS")
fit.nested<-readRDS("~/Desktop/PhD/Metagenome/Function/limmaVoomFit.RDS")

# city-wise comparisons

categories<-data.frame("Cities"=as.factor(sample_data(dt20)$Cities), "Landuse"=as.factor(sample_data(dt20)$Codes))
design2<-model.matrix(~Cities, categories)
dge2 <- DGEList(counts=counts)
keep2 <- filterByExpr(dge2, design2)
dge2 <- dge2[keep2,,keep.lib.sizes=FALSE]
dge2 <- calcNormFactors(dge2)

logCPM2 <- cpm(dge2, log=TRUE, prior.count=3)

fit2 <- lmFit(logCPM2, design2)
fit2 <- eBayes(fit2, trend=TRUE)
topTable(fit2, coef=ncol(design2))
saveRDS(fit,"~/Desktop/PhD/Metagenome/Function/limmaVoomFit_citiesonly.RDS")
fit.cities<-readRDS("~/Desktop/PhD/Metagenome/Function/limmaVoomFit_citiesonly.RDS")

# landuse only
counts<-t(as.data.frame(as.matrix(otu_table(dt20))))

categories<-data.frame("Cities"=as.factor(sample_data(dt20)$Cities), "Landuse"=as.factor(sample_data(dt20)$Codes))
design3<-model.matrix(~Landuse, categories)
dge3 <- DGEList(counts=counts)
keep3 <- filterByExpr(dge3, design3)
dge3 <- dge3[keep3,,keep.lib.sizes=FALSE] #error
dge3 <- calcNormFactors(dge3)

logCPM3 <- cpm(dge3, log=TRUE, prior.count=3)

fit3 <- lmFit(logCPM3, design3)
fit3 <- eBayes(fit3, trend=TRUE)
topTable(fit3, number=100, coef=ncol(design3))
saveRDS(fit3,"~/Desktop/PhD/Metagenome/Function/limmaVoomFit_landuseonly.RDS")
fit.landuse<-readRDS("~/Desktop/PhD/Metagenome/Function/limmaVoomFitlanduseonly.RDS")
#<-data.frame(fit$stdev.unscaled)
coeffs<-data.frame(fit$coefficients)
rownames(coeffs)
coeffs$CitiesLakti.LanduseTurf


# two way anova

design4<-model.matrix(~Landuse*Cities, categories)
dge4 <- DGEList(counts=counts)
keep4 <- filterByExpr(dge4, design4)
dge4 <- dge4[keep4,,keep.lib.sizes=FALSE] #error
dge4 <- calcNormFactors(dge4)

logCPM4 <- cpm(dge4, log=TRUE, prior.count=3)

fit4 <- lmFit(logCPM4, design4)
fit4 <- eBayes(fit4, trend=TRUE)
topTable(fit4, number=100, coef=ncol(design4))
saveRDS(fit4,"~/Desktop/PhD/Metagenome/Function/limmaVoomFit_twowayfactorial.RDS")
fit.landuse<-readRDS("~/Desktop/PhD/Metagenome/Function/limmaVoomFitlanduseonly.RDS")

# level 1 limma voom

# correct statistical model
countsL1<-t(as.data.frame(as.matrix(otu_table(dt20l1))))

design<-model.matrix(~Cities/Landuse, categories)
dgeL1 <- DGEList(counts=countsL1)
keepL1 <- filterByExpr(dgeL1, design)
dgeL1 <- dgeL1[keepL1,,keep.lib.sizes=FALSE] #error
dgeL1 <- calcNormFactors(dgeL1)

logCPML1 <- cpm(dgeL1, log=TRUE, prior.count=3)

fitL1 <- lmFit(logCPML1, design)
fitL1 <- eBayes(fitL1, trend=TRUE)
topTable(fitL1, number=14, coef=ncol(design))
saveRDS(fit,"~/Desktop/PhD/Metagenome/Function/limmaVoomFitL1.RDS")

# city-wise comparisons

categories<-data.frame("Cities"=as.factor(sample_data(dt20)$Cities), "Landuse"=as.factor(sample_data(dt20)$Codes))
design2<-model.matrix(~Cities, categories)
dge2L1 <- DGEList(counts=countsL1)
keep2L1 <- filterByExpr(dge2L1, design2)
dge2L1 <- dge2L1[keep2L1,,keep.lib.sizes=FALSE]
dge2L1 <- calcNormFactors(dge2L1)

logCPM2L1 <- cpm(dge2L1, log=TRUE, prior.count=3)

fit2L1 <- lmFit(logCPM2L1, design2)
fit2L1 <- eBayes(fit2L1, trend=TRUE)
topTable(fit2L1,number=20, coef=ncol(design2))
saveRDS(fit,"~/Desktop/PhD/Metagenome/Function/limmaVoomFitL1_citiesonly.RDS")


# landuse only
counts<-t(as.data.frame(as.matrix(otu_table(dt20))))

categories<-data.frame("Cities"=as.factor(sample_data(dt20)$Cities), "Landuse"=as.factor(sample_data(dt20)$Codes))
design3<-model.matrix(~Landuse, categories)
dge3L1 <- DGEList(counts=countsL1)
keep3L1 <- filterByExpr(dge3L1, design3)
dge3L1 <- dge3L1[keep3L1,,keep.lib.sizes=FALSE] #error
dge3L1 <- calcNormFactors(dge3L1)

logCPM3L1 <- cpm(dge3L1, log=TRUE, prior.count=3)

fit3L1 <- lmFit(logCPM3L1, design3)
fit3L1 <- eBayes(fit3L1, trend=TRUE)
topTable(fit3L1, number=14, coef=ncol(design3))
saveRDS(fit3L1,"~/Desktop/PhD/Metagenome/Function/limmaVoomFitL1_landuseonly.RDS")


# two way anova

design4<-model.matrix(~Landuse*Cities, categories)
dge4L1 <- DGEList(counts=countsL1)
keep4L1 <- filterByExpr(dge4L1, design4)
dge4L1 <- dge4L1[keep4L1,,keep.lib.sizes=FALSE] #error
dge4L1 <- calcNormFactors(dge4L1)

logCPM4L1 <- cpm(dge4L1, log=TRUE, prior.count=3)

fit4L1 <- lmFit(logCPM4L1, design4)
fit4L1 <- eBayes(fit4L1, trend=TRUE)
topTable(fit4L1, number=100, coef=ncol(design4))
saveRDS(fit4L1,"~/Desktop/PhD/Metagenome/Function/limmaVoomFitL1_twowayfactorial.RDS")


#<-data.frame(fit$stdev.unscaled)
coeffs<-data.frame(fit$coefficients)
rownames(coeffs)
coeffs$CitiesLakti.LanduseTurf
boxplot(coeffs[2:19],outcol="red")
library(plyr)
qt<-ldply(as.list(coeffs), quantile, probs=c(1,99)/100, names=TRUE)




fit$p.value
# parse p value table ####

# barplots ####

#methanogenesis




summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

dat<-sam
as.data.frame(dat)
NH4SE<-summarySE(dat, measurevar="NH4_N", groupvars = "Codes", na.rm=T)
NO3SE<-summarySE(dat, measurevar="NO3_N", groupvars="Codes", na.rm=T)
totSE<-summarySE(dat, measurevar="total_N", groupvars="Codes", na.rm=T)
AOASE<-summarySE(dat, measurevar="AmoA_copies", groupvars="Codes", na.rm=T)
AOBSE<-summarySE(dat, measurevar="AmoB_copies", groupvars="Codes", na.rm=T)
NH4SECity<-summarySE(dat, measurevar="NH4_N", groupvars = c("Cities"), na.rm=T)
Cobalt<-summarySE(dat, measurevar="Co_tot", groupvars=c("Cities"), na.rm=T)
Cadmium<-summarySE(dat, measurevar="Cd_tot", groupvars=c("Cities"), na.rm=T)
Nickel<-summarySE(dat, measurevar="Ni_tot", groupvars=c("Cities"), na.rm=T)
Zinc<-summarySE(dat, measurevar="Zn_tot", groupvars=c("Cities"), na.rm=T)
availCobalt<-summarySE(dat, measurevar="Co_avail", groupvars=c("Codes"), na.rm=T)
availNickel<-summarySE(dat, measurevar="Ni_avail", groupvars=c("Codes"), na.rm=T)

ggplot(availCobalt, aes(x=Codes, y=Co_avail)) +
  geom_errorbar(aes(ymin=Co_avail-se, ymax=Co_avail+se), width=.2, size=1) +
  geom_line(size=1) +
  labs(title="Available Co", xlab="Land-use Category", ylab="Co...")+
  geom_point()+theme_bw()

ggplot(availNickel, aes(x=Codes, y=Ni_avail)) +
  geom_errorbar(aes(ymin=Ni_avail-se, ymax=Ni_avail+se), width=.2, size=1) +
  labs(title="Available Ni", xlab="Land-use Category", ylab="Ni...")+
  geom_line(size=1) +
  geom_point()+theme_bw()


ggplot(Cobalt, aes(x=Cities, y=Co_tot)) +
  geom_errorbar(aes(ymin=Co_tot-se, ymax=Co_tot+se), width=.2, size=1) +
  geom_line(size=1) +
  labs(title="Total Co", xlab="Land-use Category", ylab="Co...")+
  geom_point()+theme_bw()

ggplot(Cadmium, aes(x=Cities, y=Cd_tot)) +
  geom_errorbar(aes(ymin=Cd_tot-se, ymax=Cd_tot+se), width=.1,size=1) +
  labs(title="Total Cd", xlab="Land-use Category", ylab="Cd...")+
  geom_line(size=1) +
  geom_point()+theme_bw()

ggplot(Nickel, aes(x=Cities, y=Ni_tot)) +
  geom_errorbar(aes(ymin=Ni_tot-se, ymax=Ni_tot+se), width=.1,size=1) +
  labs(title="Total Ni", xlab="Land-use Category", ylab="Ni ...")+
  geom_line(size=1) +
  geom_point()+theme_bw()

ggplot(Zinc, aes(x=Cities, y=Zn_tot)) +
  geom_errorbar(aes(ymin=Zn_tot-se, ymax=Zn_tot+se), width=.1, size=1) +
  labs(title="Total Zn", xlab="Land-use Category", ylab="Zn...")+
  geom_line(size=1) +
  geom_point()+theme_bw()

summary(aov(sample_data(dt15.l4)$Cd_tot~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))
summary(aov(sample_data(dt15.l4)$Co_tot~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))
summary(aov(sample_data(dt15.l4)$Zn_tot~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))
summary(aov(sample_data(dt15.l4)$Ni_tot~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))


summary(aov(sample_data(dt15.l4)$NH4_N~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))
summary(aov(sample_data(dt15.l4)$NO3_N~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))

TukeyHSD(aov(sample_data(dt15.l4)$Cd_tot~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))
TukeyHSD(aov(sample_data(dt15.l4)$Co_tot~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))
TukeyHSD(aov(sample_data(dt15.l4)$Ni_tot~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))

summary(aov(sample_data(dt15.l4)$AmoA_copies~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))
summary(aov(sample_data(dt15.l4)$AmoB_copies~sample_data(dt15.l4)$Cities*sample_data(dt15.l4)$Codes))


dat2<-data.frame("DissimilatoryNR"=as.numeric(sample_sums(Nx.Reductase)),
                 "N_reductase"=as.numeric(sample_sums(Nitrite.R.All)),
                 "AmmoniaOxidation"=as.numeric(sample_sums(RDT.ammoa)),
                "Nickel_Transport"=as.numeric(sample_sums(NickelTransport)),
                "Cobalt_transport"=as.numeric(sample_sums(CobaltTransport)),
                "ArsenicResistance"=as.numeric(sample_sums(ArsenicResistance)),
                "CopperResistance"=as.numeric(sample_sums(CopperRegulation)),
                "CadmiumResistance"=as.numeric(sample_sums(CadmiumResistance)),
                "CZC_Resistance"=as.numeric(sample_sums(CZC)),
                "ZincTransport"=as.numeric(sample_sums(ZincResistance)),
                "AntibioticResistance"=as.numeric(sample_sums(Antibiotics)),
                "Nitrogen_fixation"=as.numeric(sample_sums(Nitrogenase)))
                 
DNRA<-summarySE(dat2, measurevar="DissimilatoryNR", groupvars = "Codes", na.rm=T)
DNRA2<-summarySE(dat2, measurevar="N_reductase", groupvars = "Codes", na.rm=T)
Ammonia_oxidation<-summarySE(dat2, measurevar="AmmoniaOxidation", groupvars = "Codes", na.rm=T)
NickelT<-summarySE(dat2, measurevar="Nickel_Transport", groupvars = "Codes", na.rm=T)
CoT<-summarySE(dat2, measurevar="Cobalt_transport", groupvars = "Codes", na.rm=T)
ArR<-summarySE(dat2, measurevar="ArsenicResistance", groupvars = "Codes", na.rm=T)
CuR<-summarySE(dat2, measurevar="CopperResistance", groupvars="Codes", na.rm=T)
CdR<-summarySE(dat2, measurevar="CadmiumResistance", groupvars = "Codes", na.rm=T)
Nitrogen_fixation<-summarySE(dat2, measurevar="Nitrogen_fixation", groupvars = "Codes", na.rm=T)
CZCR<-summarySE(dat2, measurevar="CZC_Resistance", groupvars = "Codes", na.rm=T)
ZnT<-summarySE(dat2, measurevar="ZincTransport", groupvars = "Codes", na.rm=T)
AntibioR<-summarySE(dat2, measurevar="AntibioticResistance", groupvars = "Codes", na.rm=T)

dat4<-data.frame("Thaum"=sample_sums(subset_taxa(prokOrder, Family=="Nitrososphaeraceae")), "Codes"=sample_data(prokOrder)$Codes)
thaumarc<-summarySE(dat4, measurevar="Thaum", groupvars = "Codes", na.rm=T)
thP<-ggplot(thaumarc, aes(x=Codes, y=Thaum)) +
  geom_errorbar(aes(ymin=Thaum-se, ymax=Thaum+se), width=.1) +
  geom_line() +
  geom_point()+
  theme_bw()

thP+geom_point(data=b, aes(color="red"))


plot(sample_sums(RDT.ammoa)~sample_data(RDT.ammoa)$Codes)

ggplot(Nabund, aes(x=Codes, y=SS)) +
  geom_errorbar(aes(ymin=SS-se, ymax=SS+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Tabund, aes(x=Codes, y=SS2)) +
  geom_errorbar(aes(ymin=SS2-se, ymax=SS2+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(NH4SECity, aes(x=Codes, y=NH4_N, group=Cities, color=Cities)) +
  geom_errorbar(aes(ymin=NH4_N-se, ymax=NH4_N+se), width=.1) +
  geom_line() +
  geom_point()+
  scale_color_viridis_d()+theme_bw()

ggplot(TotalNCity, aes(x=Codes, y=total_N, group=Cities, color=Cities)) +
  geom_errorbar(aes(ymin=total_N-se, ymax=total_N+se), width=.1) +
  geom_line() +
  geom_point()+scale_color_viridis_d()+theme_bw()

ggplot(NH4SE, aes(x=Codes, y=NH4_N)) + 
  geom_errorbar(aes(ymin=NH4_N-se, ymax=NH4_N+se), width=.1, size=1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(NO3SE, aes(x=Codes, y=NO3_N)) + 
  geom_errorbar(aes(ymin=NO3_N-se, ymax=NO3_N+se), width=.1, size=1) +
  geom_line() +
  geom_point()+scale_color_viridis_d()+theme_bw()

ggplot(totSE, aes(x=Codes, y=total_N)) + 
  geom_errorbar(aes(ymin=total_N-se, ymax=total_N+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(AOASE, aes(x=Codes, y=AmoA_copies)) + 
  geom_errorbar(aes(ymin=AmoA_copies-se, ymax=AmoA_copies+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(AOBSE, aes(x=Codes, y=AmoB_copies)) + 
  geom_errorbar(aes(ymin=AmoB_copies-se, ymax=AmoB_copies+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(DNRA, aes(x=Codes, y=DissimilatoryNR)) + 
  geom_errorbar(aes(ymin=DissimilatoryNR-se, ymax=DissimilatoryNR+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Ammonia_oxidation, aes(x=Codes, y=AmmoniaOxidation)) + 
  geom_errorbar(aes(ymin=AmmoniaOxidation-se, ymax=AmmoniaOxidation+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()


ggplot(NickelT, aes(x=Codes, y=Nickel_Transport)) + 
  geom_errorbar(aes(ymin=Nickel_Transport-se, ymax=Nickel_Transport+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(DNRA2, aes(x=Codes, y=N_reductase)) + 
  geom_errorbar(aes(ymin=N_reductase-se, ymax=N_reductase+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Nitrate_and_nitrite_ammonification, aes(x=Codes, y=Nitrate_and_nitrite_ammonification)) + 
  geom_errorbar(aes(ymin=Nitrate_and_nitrite_ammonification-se, ymax=Nitrate_and_nitrite_ammonification+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Nitrate_and_nitrite_ammonificationC, aes(x=Codes, y=Nitrate_and_nitrite_ammonification, group=Cities, color=Cities)) + 
  geom_errorbar(aes(ymin=Nitrate_and_nitrite_ammonification-se, ymax=Nitrate_and_nitrite_ammonification+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Nitrilase, aes(x=Codes, y=Nitrilase)) + 
  geom_errorbar(aes(ymin=Nitrilase-se, ymax=Nitrilase+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Nitrogen_fixation, aes(x=Codes, y=Nitrogen_fixation)) + 
  geom_errorbar(aes(ymin=Nitrogen_fixation-se, ymax=Nitrogen_fixation+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Nitrosative_stress, aes(x=Codes, y=Nitrosative_stress)) + 
  geom_errorbar(aes(ymin=Nitrosative_stress-se, ymax=Nitrosative_stress+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(nitrile_hydratase_functions, aes(x=Codes, y=nitrile_hydratase_functions)) + 
  geom_errorbar(aes(ymin=nitrile_hydratase_functions-se, ymax=nitrile_hydratase_functions+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

ggplot(Nitric_oxide_synthase, aes(x=Codes, y=Nitric_oxide_synthase)) + 
  geom_errorbar(aes(ymin=Nitric_oxide_synthase-se, ymax=Nitric_oxide_synthase+se), width=.1) +
  geom_line() +
  geom_point()+theme_bw()

plot_bar(DT.Nitrogen, x="Codes", fill="L3", facet_grid="Cities")+geom_bar(stat="identity")+theme_bw()
plot_bar(RDT.Nitrogen, x="Codes", fill="L3")+geom_bar(stat="identity")+theme_bw()

dat3<-data.frame(otu_table(DT.Nitrogen3), sample_data(DT.Nitrogen3)$Codes, sample_data(DT.Nitrogen3)$Cities)

summarySE2<-function(x){
  a<-summarySE(dat3,measurevar=x,grouvars="Codes",na.rm=T)
  a}
allN<-ddply(dat3, dat3$Cities, summarySE2)

# bacterial barplots


plot_bar(transform_sample_counts(bac, function(x)x/sum(x)), x="Codes", fill="  Phylum", facet_grid="Cities")+geom_bar(stat="identity")+theme_bw()


plot_bar(bac, fill="Phylum", facet_grid="Codes")

# bac differential expression ####

tax<-paste(tax_table(bac)[,6], c(1:55700), sep="")
taxa_names(bac)<-tax

counts<-t(as.data.frame(as.matrix(otu_table(bac))))
geneid <- rownames(counts)

categories<-data.frame("Cities"=as.factor(sample_data(bac)$Cities), "Landuse"=as.factor(sample_data(bac)$Codes))


bacdesign<-model.matrix(~0+Landuse, categories)
#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(bacdesign))
contr.matrix


bacdge <- DGEList(counts=counts)


backeep <- filterByExpr(bacdge, bacdesign)
bacdge <- bacdge[backeep,,keep.lib.sizes=FALSE]

bacdge <- calcNormFactors(bacdge)

dge3$genes<-geneid
blcpm<-cpm(bacdge, log=TRUE)
col.group <- categories$Cities
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
sym<-as.factor(categories$Landuse)
plotMDS(blcpm, col=col.group, pch=as.numeric(sym))

v <- voom(bacdge, bacdesign, plot=TRUE)
#dup<-duplicateCorrelation(v, design3, block=categories$Cities)

bacfitV <- lmFit(v, bacdesign)

bacfitV <- contrasts.fit(bacfitV, contrasts=contr.matrix)
bacfitV <- eBayes(bacfitV, trend=TRUE)
topgenes<-topTable(bacfitV, number=100, coef=3)
siggenes<-topgenes[topgenes$adj.P.Val<0.05,]

summary(decideTests(bacfitV)) # run all together; naming scheme overwrites/not unique!!
bacsigsummary<-decideTests(bacfitV)
bac.TurfP <- which(bacsigsummary[,1]==1)
bac.TurfN <- which(bacsigsummary[,1]==-1)
bac.RemP <- which(bacsigsummary[,2]==1)
bac.RemN <- which(bacsigsummary[,2]==-1)
bac.RudP <- which(bacsigsummary[,3]==1)
bac.RudN <- which(bacsigsummary[,3]==-1)

names(bac.TurfP)
names(bac.TurfN)
bac.RemP
bac.RemN
names(bac.RudP)
names(bac.RudN)


t.bac<-as.matrix(tax_table(bac))
write.csv(t.bac[names(bac.TurfP),], "~/Desktop/PhD/Metagenome/TURFupbac.csv")
write.csv(t.bac[names(bac.TurfN),], "~/Desktop/PhD/Metagenome/TURFdownbac.csv")
write.csv(t.bac[names(bac.RudP),], "~/Desktop/PhD/Metagenome/RUDupbac.csv")
write.csv(t.bac[names(bac.RudN),], "~/Desktop/PhD/Metagenome/RUDdownbac.csv")

# archaeal differential abundance ####
# prokaryote indicators ####

prokOrder<-tax_glom(baclib, "Family")
tax_table(prokOrder)
proktax<-paste(tax_table(prokOrder)[,5], c(1:356), sep="")
taxa_names(prokOrder)<-proktax

counts<-t(as.data.frame(as.matrix(otu_table(prokOrder))))
geneid <- rownames(counts)

categories<-data.frame("Cities"=as.factor(sample_data(prokOrder)$Cities), "Landuse"=as.factor(sample_data(prokOrder)$Codes))
prokdesign<-model.matrix(~0+Landuse, categories)
prokdge <- DGEList(counts=counts)
prokeep <- filterByExpr(prokdge, prokdesign)
prokdge <- prokdge[prokeep,,keep.lib.sizes=FALSE] # doesn't work on relativised/proportion  data
prokdge <- calcNormFactors(prokdge)

logCPMprok <- cpm(prokdge, log=TRUE, prior.count=3)


#colnames(design) <- gsub("group3", "", colnames(design))
contr.matrix <- makeContrasts(
  TurfvRef = LanduseTurf-LanduseReference, 
  RemvRef = LanduseRemnant-LanduseReference, 
  RudvRef = LanduseRuderal-LanduseReference, 
  levels = colnames(prokdesign))
contr.matrix

#prokdge$genes<-geneid
# maybe skip this
proklcpm<-cpm(prokdge, log=TRUE)
col.group <- categories$Cities
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
sym<-as.factor(categories$Landuse)
plotMDS(blcpm, col=col.group, pch=as.numeric(sym))
## ok we're back
prokv <- voom(prokdge, prokdesign, plot=TRUE)
#dup<-duplicateCorrelation(v, design3, block=categories$Cities)

prokfitV <- lmFit(prokv, prokdesign)

prokfitV <- contrasts.fit(prokfitV, contrasts=contr.matrix)
prokfitV <- eBayes(prokfitV, trend=TRUE)
topgenes<-topTable(prokfitV, number=100, coef=3)
siggenes<-topgenes[topgenes$adj.P.Val<0.05,]

summary(decideTests(prokfitV)) # run all together; naming scheme overwrites/not unique!!
proksigsummary<-decideTests(prokfitV)
prok.TurfP <- which(proksigsummary[,1]==1)
prok.TurfN <- which(proksigsummary[,1]==-1)
prok.RemP <- which(proksigsummary[,2]==1)
prok.RemN <- which(proksigsummary[,2]==-1)
prok.RudP <- which(proksigsummary[,3]==1)
prok.RudN <- which(proksigsummary[,3]==-1)

names(prok.TurfP)
names(prok.TurfN)
bac.RemP
bac.RemN
names(prok.RudP)
names(prok.RudN)



# network analysis with bacteria
# import bacteria
library(phyloseq)
library(stringr)
library(vegan)
baclib<-readRDS("~/Desktop/PhD/Metagenome/Bac16s_rarefied.RDS") #need to change taxa names from sequences... impossible to interpret.

sample_names(dt20)<-str_replace(sample_names(dt20), "Final_", "GLU")
sample_names(bac)

orgfuncor<-cor(as.data.frame(as.matrix(otu_table(bac))), as.data.frame(as.matrix(otu_table(dt20))), method="spearman")
orgFuncordist<-vegdist(orgfuncor, method="bray") #bray chosen, but may highlight differences rather than similarities?
hc<-hclust(orgFuncordist, method="ward.D2") # ward.D2 aims for hight circular clusters in the network

# hclust does the module identification based on dist values; use bray-curtis dist values that embody the shared sites for each function/taxon 

#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

plot(as.phylo(hc),type = "unrooted", cex = 0.6, no.margin = TRUE)


# none of these color schemes have been tested yet...
dend <- as.dendrogram(hcl1)
col_aa_red <- ifelse(grepl("aa", labels(dend)), "red", "blue")#defines color scheme
dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = col_aa_red, edgePar = "col")
plot(dend2)

dend <- dist %>% hclust %>% as.phylo %>%
           set("branches_k_color", k=3) %>% set("branches_lwd", 1.2) %>%
           set("labels_colors") %>% set("labels_cex", c(.9,1.2)) %>% 
           set("leaves_pch", 19) %>% set("leaves_col", c("blue", "red"))
           # plot the dend in usual "base" plotting engine:
plot(dend)

#####

# Function ordination L 4 relative ####
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
l4nmds<-ordinate(rdt15.l4, method="NMDS", distance="bray")
pl4<-plot_ordination(rdt15.l4, l4nmds, type="samples", color="Cities", shape="Codes")+
  #scale_color_grey()+
  #scale_color_viridis_d()+
  scale_colour_manual(values=cbPalette)+
  geom_point(size=5)
pl4 + theme_bw()

funOTU<-as.matrix(otu_table(rdt15.l4))
funsam<-as.data.frame(sample_data(rdt15.l4))
funsam<-funsam[,c(4,5,7,10:24)]
funmeta<-metaMDS(funOTU, distance="bray")
plot(funmeta, display="sites", colour=funsam$Cities)
attach(funsam)
funfit<-envfit(funmeta~Codes+Cities+pH.H2O+C_org+CaCO3+AL_K2O+AL_P2O5+NH4_N+NO3_N+total_N+Cd_tot+Co_tot+Ni_tot+Zn_tot, na.rm=TRUE)
plot(funfit, col="black")

funfitarrows<-data.frame(labels=rownames(funfit$vectors$arrows), funfit$vectors$arrows)
funfitarrows$r<-funfit$vectors$r
funfitarrows$labels<-c("pH", "OM", "C03", "K", "P", "NH4", "NO3", "TN", "Cd", "Co", "Ni", "Zn")

arrow_map = aes(xend = NMDS1*r, yend = NMDS2*r, x = 0, y = 0, shape = NULL, color = NULL)
label_map = aes(x = 1.1*NMDS1*r, y = 1.1*NMDS2*r, shape = NULL, color = NULL, 
                label = labels)
arrowhead = arrow(length = unit(0.03, "npc"))
plot4 = ggplot() + geom_segment(arrow_map, data = funfitarrows, color = "black", 
                        arrow = arrowhead) + geom_text(label_map, size = 5, data = funfitarrows, color="black")
plot4+theme_bw()


# Function ordination L 4 Density ####
#density normalization

Dl4nmds<-ordinate(L4Norm, method="NMDS", distance="bray")
pl4D<-plot_ordination(L4Norm, Dl4nmds, type="samples", color="Cities", shape="Codes")+
  #scale_color_grey()+
  scale_colour_manual(values=cbPalette)+
  geom_point(size=5)
pl4D + theme_bw()

DfunOTU<-as.matrix(otu_table(L4Norm))
Dfunsam<-as.data.frame(sample_data(L4Norm))
Dfunsam<-Dfunsam[,c(4,5,7,10:24)]
Dfunmeta<-metaMDS(DfunOTU, distance="bray")
#plot(Dfunmeta, display="sites", colour=funsam$Cities)
attach(Dfunsam)
Dfunfit<-envfit(Dfunmeta~Codes+Cities+pH.H2O+C_org+CaCO3+AL_K2O+AL_P2O5+NH4_N+NO3_N+total_N+Cd_tot+Co_tot+Ni_tot+Zn_tot, na.rm=TRUE)
#plot(funfit, col="black")

Dfunfitarrows<-data.frame(labels=rownames(Dfunfit$vectors$arrows), Dfunfit$vectors$arrows)
Dfunfitarrows$r<-Dfunfit$vectors$r
Dfunfitarrows$labels<-c("pH", "OM", "Carbonate", "K", "Phosphate", "Ammonium", "Nitrate", "Total N", "Total Cd", "Total Co", "Total Ni", "Total Zn")

arrow_map = aes(xend = 1.8*NMDS1*r, yend = 1.8*NMDS2*r, x = 0, y = 0, shape = NULL, color = NULL)
label_map = aes(x = 2.1*NMDS1*r+c(-0.02, 0, 0, 0,0,0,-0.01,0,0.05, 0, 0.001, -0.025), y = 2.1*NMDS2*r+c(0.03, 0, 0.02, 0,0,0,0,0,-0.02, -0.035, -0.02, -0.03), shape = NULL, color = NULL, 
                label = labels)
arrowhead = arrow(length = unit(0.03, "npc"))
Dplot4 = pl4D + geom_segment(arrow_map, data = Dfunfitarrows, color = "grey", 
                        arrow = arrowhead) +
  
  geom_text(label_map, size = 3, data = Dfunfitarrows, color="black")
  #annotate("text", Dfunfitarrows$NMDS1[3], Dfunfitarrows$NMDS2[3], bquote("CaCO"[3]))
  #labs(title=bquote("CaCO"[3]))
Dplot4+theme_bw()

# bacterial mg-rast ordinations ####
mgbac<-readRDS("~/Desktop/PhD/Metagenome/Organism/ct_species.RDS")
mgbac2<-readRDS("~/Desktop/PhD/Metagenome/Organism/ct_species.RDS")

sample_data(mgbac)<-sam
#mgbac2<-subset_samples(mgbac2, sample_sums(mgbac)>1000000)
factor2<-round(1000000*(sample_data(mgbac)$gperg_DNA/(mean(sample_data(mgbac)$gperg_DNA))))

table(sample_sums(mgbac)>factor2)
mgbacOTU<-otu_table(mgbac)
mgbacnormalizedOTU<-rrarefy(mgbacOTU, factor2)
otu_table(mgbac)<-mgbacnormalizedOTU
mgbac<-prune_samples(sample_sums(mgbac)==factor2, mgbac)


mgbacnmds<-ordinate(mgbac, method="NMDS", distance="bray")
mgp1<-plot_ordination(mgbac, mgbacnmds, type="samples", color="Cities", shape="Codes")+
  #scale_color_grey()+
  scale_colour_manual(values=cbPalette)+
  geom_point(size=4)
mgp1+theme_bw()

mgbac2nmds<-ordinate(mgbac2, method="NMDS", distance="bray")
mgp2<-plot_ordination(mgbac2, mgbac2nmds, type="samples", color="Cities", shape="Codes")+
  #scale_color_grey()+
  scale_colour_manual(values=cbPalette)+
  geom_point(size=4)
mgp2+theme_bw()

mgbacOTU<-as.matrix(otu_table(mgbac))
mgbacsam<-as.data.frame(sample_data(mgbac))
mgbacsam<-mgbacsam[,c(4,5,7,10:20)]
mgbacmeta<-metaMDS(mgbacOTU, distance="bray")

mgbacOTU2<-as.matrix(otu_table(mgbac2))
mgbacsam2<-as.data.frame(sample_data(mgbac2))
mgbacsam2<-mgbacsam2[,c(4,5,7,10:20)]
mgbacmeta2<-metaMDS(mgbacOTU2, distance="bray")

attach(mgbacsam)

#bacfit<-envfit(bacmeta~bacsam$Codes+bacsam$Cities+bacsam$pH.H2O+bacsam$C_org+bacsam$CaCO3+bacsam$AL_K2O+bacsam$AL_P2O5+bacsam$NH4_N+bacsam$NO3_N+bacsam$total_N+bacsam$Cd_tot+bacsam$Co_tot+bacsam$Ni_tot, na.rm=TRUE)
mgbacfit2<-envfit(mgbacmeta2~Codes+Cities+pH.H2O+C_org+CaCO3+AL_K2O+AL_P2O5+NH4_N+NO3_N+total_N+Cd_tot+Co_tot+Ni_tot+Zn_tot, na.rm=TRUE)

mgbacfit<-envfit(mgbacmeta~Codes+Cities+pH.H2O+C_org+CaCO3+AL_K2O+AL_P2O5+NH4_N+NO3_N+total_N+Cd_tot+Co_tot+Ni_tot+Zn_tot, na.rm=TRUE)
plot(mgbacfit, col="red") # 

mgbacfitarrows<-data.frame(labels=rownames(mgbacfit$vectors$arrows), mgbacfit$vectors$arrows)
mgbacfitarrows$r<-mgbacfit$vectors$r
mgbacfitarrows$labels<-c("pH", "OM", "Carbonate", "K", "Phosphate", "Ammonium", "Nitrate", "TN", "Total Cd", "Total Co", "Total Ni", "Total Zn")
arrow_map = aes(xend = NMDS1*r, yend = NMDS2*r, x = 0, y = 0, shape = NULL, color = NULL)
label_map = aes(x = 1.1*NMDS1*r+c(0,0,0.02,0.02,0.02,0,0.02,0,0.01,-0.02,-0.02,-0.02), y = 1.1*NMDS2*r+c(0,0,0,0,0,0,0,0,0.04,0,0,0), shape = NULL, color = NULL, 
                label = labels)
arrowhead = arrow(length = unit(0.03, "npc"))
mgp3 = ggplot() + geom_segment(arrow_map, data = mgbacfitarrows, color = "black", 
                        arrow = arrowhead) + geom_text(label_map, size = 4, data = mgbacfitarrows)
mgp3+theme_bw()





# bacterial ordinations ####
bacnmds<-ordinate(bac, method="NMDS", distance="bray")
bp1<-plot_ordination(bac, bacnmds, type="samples", color="Cities", shape="Codes")+
  scale_colour_manual(values=cbPalette)+
  geom_point(size=4)
bp1+theme_bw()

bacOTU<-otu_table(bac)
bacsam<-as.data.frame(sample_data(bac))
bacsam<-bacsam[,c(4,5,7,10:24)]
bacmeta<-metaMDS(bacOTU, distance="bray")
#plot(bacmeta, display="sites") # this usually produces the same output as in phyloseq, make sure the plots are nearly identical before doing envfit, so that the vectors accurately correlate between plots


attach(bacsam)

#bacfit<-envfit(bacmeta~bacsam$Codes+bacsam$Cities+bacsam$pH.H2O+bacsam$C_org+bacsam$CaCO3+bacsam$AL_K2O+bacsam$AL_P2O5+bacsam$NH4_N+bacsam$NO3_N+bacsam$total_N+bacsam$Cd_tot+bacsam$Co_tot+bacsam$Ni_tot, na.rm=TRUE)
bacfit<-envfit(bacmeta~Codes+Cities+pH.H2O+C_org+CaCO3+AL_K2O+AL_P2O5+NH4_N+NO3_N+total_N+Cd_tot+Co_tot+Ni_tot+Zn_tot, na.rm=TRUE)
plot(bacfit, col="red") # 

bacfitarrows<-data.frame(labels=rownames(bacfit$vectors$arrows), bacfit$vectors$arrows)
bacfitarrows$labels<-c("pH", "OM", "C03", "K", "P", "NH4", "NO3", "TN", "Cd", "Co", "Ni", "Zn")
bacfitarrows$r<-bacfit$vectors$r
arrow_map = aes(xend = NMDS1*r, yend = NMDS2*r, x = 0, y = 0, shape = NULL, color = NULL)
label_map = aes(x = 1.1*NMDS1*r, y = 1.1*NMDS2*r, shape = NULL, color = NULL, 
                label = labels)
arrowhead = arrow(length = unit(0.03, "npc"))
p3 = ggplot() + geom_segment(arrow_map, data = bacfitarrows, color = "black", 
                       arrow = arrowhead) + geom_text(label_map, size = 3, data = bacfitarrows)

p3+theme_bw()

# Make a new graphic


arc<-prune_samples(sample_sums(arc)>5, arc)
arcnmds<-ordinate(arc, method="NMDS", distance="bray")
ap1<-plot_ordination(arc, arcnmds, type="samples", color="Cities", shape="Codes")+
  scale_colour_manual(values=cbPalette)+
  geom_point(size=4)
ap1+theme_bw()

arcOTU<-as.matrix(otu_table(arc))
arcsam<-as.data.frame(sample_data(arc))
arcsam<-arcsam[,c(4,5,7,10:20)]
arcmeta<-metaMDS(arcOTU, distance="bray")
plot(arcmeta, display="sites")

attach(arcsam)
arcfit<-envfit(arcmeta~Codes+Cities+pH.H2O+C_org+CaCO3+AL_K2O+AL_P2O5+NH4_N+NO3_N+total_N+Cd_tot+Co_tot+Ni_tot+Zn_tot, na.rm=TRUE)
plot(arcfit, col="red")

arcfitarrows<-data.frame(labels=rownames(arcfit$vectors$arrows), arcfit$vectors$arrows)
arcfitarrows$labels<-c("pH", "OM", "C03", "K", "P", "NH4", "NO3", "TN", "Cd", "Co", "Ni", "Zn")
arcfitarrows$r<-arcfit$vectors$r
arrow_map = aes(xend = 2.2*NMDS1*r, yend = 2.2*NMDS2*r, x = 0, y = 0, shape = NULL, color = NULL)
label_map = aes(x = 2.5*NMDS1*r, y = 2.5*NMDS2*r, shape = NULL, color = NULL, 
                label = labels)
arrowhead = arrow(length = unit(0.03, "npc"))
arcp3 = ap1 + geom_segment(arrow_map, data = arcfitarrows, color = "grey", 
                       arrow = arrowhead) + geom_text(label_map, size = 3, data = arcfitarrows)+
  scale_color_viridis_d()+
  geom_point(size=4)
arcp3+theme_bw()




# permanova and convergence ####
library(vegan)
# level 1 ontology
tab<-as.data.frame(as.matrix(otu_table(rdt15.l1)))
metadata<-data.frame(sample_data(dt15.l1))
#brayl1 <- adonis(tab~Cities/Codes, data=metadata)
brayl1 <- adonis(tab~Codes*Cities, data=metadata)
brayl1

# level 3 ontology
tabl3<-as.data.frame(as.matrix(otu_table(rdt15.l3)))
metadatal3<-data.frame(sample_data(rdt15.l3))
#brayl3 <- adonis(tabl3~Cities/Codes, data=metadatal3)
brayl3 <- adonis(tabl3~Cities*Codes, data=metadatal3)
brayl3

# level 4 ontology
tabl4<-as.data.frame(as.matrix(otu_table(rdt15.l4)))
metadatal4<-data.frame(sample_data(rdt15.l4))
brayl4 <- adonis(tabl4~Cities*Codes, data=metadatal4)
brayl4

tabl4<-as.data.frame(as.matrix(otu_table(L4Norm)))
metadatal4<-data.frame(sample_data(L4Norm))
brayl4 <- adonis(tabl4~Cities*Codes, data=metadatal4)
brayl4

# taxonomy function
tabl5<-as.data.frame(as.matrix(otu_table(mgbac)))
metadatal5<-data.frame(sample_data(mgbac))
brayl5 <- adonis(tabl5~Cities*Codes, data=metadatal5)
brayl5

mgbacOTU2
mgbacsam2<-data.frame(mgbacsam2)
brayl6 <- adonis(mgbacOTU2~Cities*Codes, data=mgbacsam2)
brayl6

# Sequence library

prokTab<-as.data.frame(as.matrix(otu_table(baclib)))
metabac<-data.frame(sample_data(baclib))

bactab1<-as.data.frame(as.matrix(otu_table(bac)))
metadata<-data.frame(sample_data(bac))
#brayl1 <- adonis(tab~Cities/Codes, data=metadata)
bacbrayl1 <- adonis(bactab1~Cities*Codes, data=metadata)
bacbrayl1  

raw.bac<-subset_taxa(baclib2, Kingdom=="Bacteria")
bactab2<-as.data.frame(otu_table(raw.bac))
metadata<-data.frame(sample_data(raw.bac))

bacbray2<-adonis(bactab2~Cities*Codes, data=metadata)
bacbray2

bootAdonis<-function(ps,n){
  require(phyloseq)
  require(vegan)
  require(plyr)
  l<-list(1:n)
  for(i in l){t<-as.data.frame(as.matrix(otu_table(ps)))
  o<-data.frame(sample_data(ps))}
  m<-matrix()
  pt<-ddply()
  p<-ddply(pt, summarise, "meanP"=mean(), "VP"=var())
  r<-ddply(pt, summarise, "meanP"=mean(), "VR"=var())
  out<-rbind(p,r)
}

# level 4 convergence 

distl4<-vegdist(otu_table(rdt15.l4), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
GroupC<-sample_data(rdt15.l4)$Codes
l4.C<-betadisper(distl4, GroupC, type=c("median"))


anova(l4.C)

permutest(l4.C)

l4.C

# function N L4/density normalized####

N.L1<-subset_taxa(rdt15.l4, L1=="Nitrogen Metabolism")
disNl4<-vegdist(otu_table(N.L1), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
GroupC<-sample_data(N.L1)$Codes
Nl4.C<-betadisper(disNl4, GroupC, type=c("median"))


anova(Nl4.C)

Nl4.C

# bacteria
GroupA=sample_data(bac)$BLK
GroupB=sample_data(bac)$BLK2
GroupC=sample_data(bac)$Codes

distbac<-vegdist(otu_table(bac), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
#bac.A<-betadisper(distbac, GroupA, type=c("median"))
#bac.B<-betadisper(distbac, GroupB, type=c("median"))
bac.C<-betadisper(distbac, GroupC, type=c("median"))

#anova(bac.A)
#anova(bac.B)
anova(bac.C)

#permutest(bac.A)
#permutest(bac.B)
permutest(bac.C)

#bac.A
#bac.B
bac.C
TukeyHSD(bac.C)

# archaea
GroupA=sample_data(bac)$BLK
GroupB=sample_data(baclib)$BLK2
GroupC=sample_data(arc)$Codes

distarc<-vegdist(otu_table(arc), method="bray", binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE)
#bac.A<-betadisper(distbac, GroupA, type=c("median"))
#bac.B<-betadisper(distbac, GroupB, type=c("median"))
arc.C<-betadisper(distarc, GroupC, type=c("median"))

#anova(bac.A)
#anova(bac.B)
anova(arc.C)

#permutest(bac.A)
#permutest(bac.B)
permutest(arc.C)

#bac.A
#bac.B
arc.C
TukeyHSD(arc.C)

# sig function convergence ####
plot_bar(dt15l1, x="Cities")

# ammonia monooxygenase

mg.amoa<-data.frame("mg"=sample_sums(RDT.ammoa), "QPCRAOA"=sample_data(RDT.ammoa)$AmoA_copies, "QPCRAOB"=sample_data(RDT.ammoa)$AmoB_copies,"Land"=sample_data(RDT.ammoa)$Codes, "City"=sample_data(RDT.ammoa)$Cities)

with(mg.amoa, cor(log10(mg), log10(QPCRAOA), use="pairwise.complete.obs"))
with(mg.amoa, plot(log10(mg)~log10(QPCRAOA), col=Land))
with(mg.amoa, plot(log10(mg)~log10(QPCRAOA), col=City, pch=c(1,2,3,4)[as.numeric(Land)]))

with(mg.amoa, plot(mg~log10(QPCRAOA), col=Land))
with(mg.amoa, plot(mg~QPCRAOB, col=Land))
with(mg.amoa, cor(mg, QPCRAOB))
with(mg.amoa, cor(mg, QPCRAOB))
with(mg.amoa, plot(mg~Land))

with(mg.amoa, leveneTest(mg, Land))
plot(sample_sums(RDT.ammoa)~sample_data(RDT.ammoa)$Codes)

summary(aov(sample_sums(RDT.ammoa)~sample_data(RDT.ammoa)$Codes))
TukeyHSD(aov(sample_sums(RDT.ammoa)~sample_data(RDT.ammoa)$Codes))

# DNRA

mg.dnra<-data.frame("mg"=sample_sums(Nx.Reductase), "mg2"=sample_sums(Nitrogenase), "Land"=sample_data(Nitrogenase)$Codes)

with(mg.dnra, leveneTest(mg, Land))
with(mg.dnra, plot(mg~Land))
with(mg.dnra, plot(mg2~Land))


plot(sample_sums(Nitrite.R.All)~sample_data(Nitrite.R.All)$Codes)

summary(aov(sample_sums(Nitrite.R.All)~sample_data(Nitrite.R.All)$Codes+sample_data(Nitrite.R.All)$Cities))
TukeyHSD(aov(sample_sums(Nitrite.R.All)~sample_data(Nitrite.R.All)$Codes))

summary(aov(sample_sums(Nitrite.Reductase)~sample_data(Nitrite.Reductase)$Codes))
TukeyHSD(aov(sample_sums(Nitrite.Reductase)~sample_data(Nitrite.Reductase)$Codes))


# heavy metal resistance

mg.Cd<-data.frame("mg"=sample_sums(CadmiumResistance), "Cd"=sample_data(CadmiumResistance)$Cd_tot, "Cdt"=sample_data(CadmiumResistance)$Cd_avail,"Land"=sample_data(CadmiumResistance)$Codes)

with(mg.Cd, cor(mg,Cd, use="pairwise.complete.obs"))
with(mg.Cd, cor(mg,Cdt, use="pairwise.complete.obs"))
with(mg.Cd, plot(mg~log10(Cd), col=Land))

mg.Co<-data.frame("mg"=sample_sums(CobaltTransport), "Co"=sample_data(CobaltTransport)$Co_tot, "Cdt"=sample_data(CobaltTransport)$Co_avail,"Land"=sample_data(CobaltTransport)$Codes)

with(mg.Co, cor(mg,Co, use="pairwise.complete.obs"))
with(mg.Co, cor(mg,Cdt, use="pairwise.complete.obs"))
with(mg.Co, plot(mg~Land))

mg.Ni<-data.frame("mg"=sample_sums(NickelTransport), "Ni"=sample_data(NickelTransport)$Ni_tot, "Nit"=sample_data(NickelTransport)$Ni_avail,"Land"=sample_data(NickelTransport)$Codes)

with(mg.Ni, cor(mg,Ni, use="pairwise.complete.obs"))
with(mg.Ni, cor(mg,Nit, use="pairwise.complete.obs"))

mg.Zn<-data.frame("mg"=sample_sums(NickelTransport), "Zn"=sample_data(NickelTransport)$Zn_tot, "Znt"=sample_data(NickelTransport)$Zn_avail,"Land"=sample_data(NickelTransport)$Codes)

with(mg.Ni, cor(mg,Ni, use="pairwise.complete.obs"))
with(mg.Ni, cor(mg,Nit, use="pairwise.complete.obs"))
# barcharts and p-values


# richness ####
plot_richness(dt15, x = "samplecodes", measures=c("Shannon")) + geom_boxplot()
plot_richness(dt15, x = "Codes", measures=c("Shannon")) + geom_boxplot()

plot_richness(dt15l1, x = "samplecodes", measures=c("Shannon")) + geom_boxplot()
plot_richness(dt15l1, x = "Codes", measures=c("Shannon")) + geom_boxplot()

summary(aov())

# network analysis with bacteria


orgfuncor<-cor(as.data.frame(as.matrix(otu_table(bac))), as.data.frame(as.matrix(otu_table(dt15))), method="spearman")
orgFuncordist<-vegdist(orgfuncor, method="bray") #bray chosen, but may highlight differences rather than similarities?
hc<-hclust(orgFuncordist, method="ward.D2") # ward.D2 aims for hight circular clusters in the network

# heatmap analysis ####
source("http://bioconductor.org/biocLite.R")
biocLite("Heatplus")
library(Heatplus)
library(vegan)
library(RColorBrewer)
disease.L1<-subset_taxa(dt15.l4, L1=="Virulence, Disease and Defense")
#disease.L1<-tax_glom(disease.L1, rank_names(disease.L1)[4])
disease.L1<-filter_taxa(disease.L1, function(x) max(x)>8000, TRUE)
table(disease.L1)
# make column dendrogram data
data.dist.g<-vegdist(t(otu_table(Frontiers.DT)), method="bray") # check inputs!!!
col.clus<-hclust(data.dist.g, "aver")
col.clus2<-as.dendrogram(col.clus)
col.clus2<-color_branches()
#make row dendrogram data
data.dist<-vegdist(otu_table(Frontiers.DT))
row.clus<-hclust(data.dist, "aver")
#make annotation data
ann.dat<-data.frame(sample_data(Frontiers.DT)$Cities, sample_data(Frontiers.DT)$Codes)
ann.dat2<-as.data.frame(as.matrix(tax_table(Frontiers.DT)))$L1
#plot the heatmap
plot(annHeatmap2(as.matrix(otu_table(disease.L1)),
                 col = colorRampPalette(c("lightblue", "red"), space = "rgb"),
                 breaks = 50,
                 dendrogram = list(Row = list(dendro = as.dendrogram(row.clus)), Col = list(dendro = as.dendrogram(col.clus))),
                 legend = 3,
                 #labels = list(Col = list(nrow = 1)),
                 ann = list(Row = list(data = ann.dat)),
                 cluster = list(Row = list(cuth = 0.4, col = brewer.pal(8, "Set2"))) # cuth gives the height at which the dedrogram should be cut to form clusters, and col specifies the colours for the clusters
))


# heatmap : https://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/
# subset permanova ####

Disease.PERM<-adonis(as.matrix(otu_table(RDT.Diseas3))~Cities*Codes, data=data.frame(sample_data(RDT.Diseas3)))
Disease.PERM

Nitrogen.PERM<-adonis(as.matrix(otu_table(RDT.Nitrogen3))~Cities*Codes, data=data.frame(sample_data(RDT.Nitrogen3)))
Nitrogen.PERM
Nitrogen.PERM$coefficients
# network ####

library(phyloseq)
library(vegan)
library(igraph)
# testing
#disease.L1N<-make_network(disease.L1, max.dist=0.2, type="taxa", dist.fun="bray")
#ig<-plot_net(Frontiers.DT, maxdist=0.2, color="Codes")
#ig

#d describes the edges of the network. Its first two columns are the IDs of the source and the target node for each edge. The following columns are edge attributes (weight, type, label, or anything else).
#vertices starts with a column of node IDs. Any following columns are interpreted as node attributes.
# http://kateto.net/networks-r-igraph
library(igraph)





List.Disease<-c("BaltimoreRemnant"=subset_samples(DT.Disease, samplecodes=="BaltimoreRemnant"),
                     "BaltimoreRuderal"=subset_samples(DT.Disease, samplecodes=="BaltimoreRuderal"),
                     "BaltimoreTurf"=subset_samples(DT.Disease, samplecodes=="BaltimoreTurf"),
                     "BaltimoreReference"=subset_samples(DT.Disease, samplecodes=="BaltimoreReference"),
                     "HelsinkiReference"=subset_samples(DT.Disease, samplecodes=="HelsinkiReference"),
                     "HelsinkiRuderal"=subset_samples(DT.Disease, samplecodes=="HelsinkiRuderal"),
                     "HelsinkiTurf"=subset_samples(DT.Disease, samplecodes=="HelsinkiTurf"),
                     "HelsinkiRemnant"=subset_samples(DT.Disease, samplecodes=="HelsinkiRemnant"),
                     "LaktiReference"=subset_samples(DT.Disease, samplecodes=="LaktiReference"),
                     "LaktiRuderal"=subset_samples(DT.Disease, samplecodes=="LaktiRuderal"),
                     "LaktiTurf"=subset_samples(DT.Disease, samplecodes=="LaktiTurf"),
                     "LaktiRemnant"=subset_samples(DT.Disease, samplecodes=="LaktiRemnant"),
                     "BudapestReference"=subset_samples(DT.Disease, samplecodes=="BudapestReference"),
                     "BudapestRuderal"=subset_samples(DT.Disease, samplecodes=="BudapestRuderal"),
                     "BudapestTurf"=subset_samples(DT.Disease, samplecodes=="BudapestTurf"),
                     "BudapestRemnant"=subset_samples(DT.Disease, samplecodes=="BudapestRemnant"),
                     "South AfricaReference"=subset_samples(DT.Disease, samplecodes=="South AfricaReference"),
                     "South AfricaRemnant"=subset_samples(DT.Disease, samplecodes=="South AfricaRemnant"),
                     "South AfricaTurf"=subset_samples(DT.Disease, samplecodes=="South AfricaTurf"),
                     "South AfricaRuderal"=subset_samples(DT.Disease, samplecodes=="South AfricaRuderal"))
#DiseaseNetwork60<-netStat(List.Disease, groups2, "Disease") # param set to 0.6
DiseaseNetwork70<-netStat(List.Disease, groups2, "Disease") # param set to 0.7
#DiseaseNetwork80<-netStat(List.Disease, groups2, "Disease") # param set to 0.8
#DiseaseNetwork90<-netStat(List.Disease, groups2, "Disease") # param set to 0.9

DiseaseNetwork500<-netStat(List.Disease, groups2, "Disease")
DiseaseNetwork1000<-netStat(List.Disease, groups2, "Disease")

List.Stress<-c("BaltimoreRemnant"=subset_samples(DT.Stress, samplecodes=="BaltimoreRemnant"),
                "BaltimoreRuderal"=subset_samples(DT.Stress, samplecodes=="BaltimoreRuderal"),
                "BaltimoreTurf"=subset_samples(DT.Stress, samplecodes=="BaltimoreTurf"),
                "BaltimoreReference"=subset_samples(DT.Stress, samplecodes=="BaltimoreReference"),
                "HelsinkiReference"=subset_samples(DT.Stress, samplecodes=="HelsinkiReference"),
                "HelsinkiRuderal"=subset_samples(DT.Stress, samplecodes=="HelsinkiRuderal"),
                "HelsinkiTurf"=subset_samples(DT.Stress, samplecodes=="HelsinkiTurf"),
                "HelsinkiRemnant"=subset_samples(DT.Stress, samplecodes=="HelsinkiRemnant"),
                "LaktiReference"=subset_samples(DT.Stress, samplecodes=="LaktiReference"),
                "LaktiRuderal"=subset_samples(DT.Stress, samplecodes=="LaktiRuderal"),
                "LaktiTurf"=subset_samples(DT.Stress, samplecodes=="LaktiTurf"),
                "LaktiRemnant"=subset_samples(DT.Stress, samplecodes=="LaktiRemnant"),
                "BudapestReference"=subset_samples(DT.Stress, samplecodes=="BudapestReference"),
                "BudapestRuderal"=subset_samples(DT.Stress, samplecodes=="BudapestRuderal"),
                "BudapestTurf"=subset_samples(DT.Stress, samplecodes=="BudapestTurf"),
                "BudapestRemnant"=subset_samples(DT.Stress, samplecodes=="BudapestRemnant"),
                "South AfricaReference"=subset_samples(DT.Stress, samplecodes=="South AfricaReference"),
                "South AfricaRemnant"=subset_samples(DT.Stress, samplecodes=="South AfricaRemnant"),
                "South AfricaTurf"=subset_samples(DT.Stress, samplecodes=="South AfricaTurf"),
                "South AfricaRuderal"=subset_samples(DT.Stress, samplecodes=="South AfricaRuderal"))
StressNetwork60<-netStat(List.Stress, groups2, "Stress") # param set to 0.6
StressNetwork70<-netStat(List.Stress, groups2, "Stress") # param set to 0.7
StressNetwork80<-netStat(List.Stress, groups2, "Disease") # param set to 0.8
StressNetwork90<-netStat(List.Stress, groups2, "Disease") # param set to 0.9

StressNetwork500<-netStat(List.Stress, groups2, "Stress")
StressNetwork2000<-netStat(List.Stress, groups2, "Stress")
StressNetwork3000<-netStat(List.Stress, groups2, "Stress")
StressNetwork6000<-netStat(List.Stress, groups2, "Stress")
StressNetwork10000<-netStat(List.Stress, groups2, "Stress")

List.Nitrogen<-c("BaltimoreRemnant"=subset_samples(DT.Nitrogen, samplecodes=="BaltimoreRemnant"),
               "BaltimoreRuderal"=subset_samples(DT.Nitrogen, samplecodes=="BaltimoreRuderal"),
               "BaltimoreTurf"=subset_samples(DT.Nitrogen, samplecodes=="BaltimoreTurf"),
               "BaltimoreReference"=subset_samples(DT.Nitrogen, samplecodes=="BaltimoreReference"),
               "HelsinkiReference"=subset_samples(DT.Nitrogen, samplecodes=="HelsinkiReference"),
               "HelsinkiRuderal"=subset_samples(DT.Nitrogen, samplecodes=="HelsinkiRuderal"),
               "HelsinkiTurf"=subset_samples(DT.Nitrogen, samplecodes=="HelsinkiTurf"),
               "HelsinkiRemnant"=subset_samples(DT.Nitrogen, samplecodes=="HelsinkiRemnant"),
               "LaktiReference"=subset_samples(DT.Nitrogen, samplecodes=="LaktiReference"),
               "LaktiRuderal"=subset_samples(DT.Nitrogen, samplecodes=="LaktiRuderal"),
               "LaktiTurf"=subset_samples(DT.Nitrogen, samplecodes=="LaktiTurf"),
               "LaktiRemnant"=subset_samples(DT.Nitrogen, samplecodes=="LaktiRemnant"),
               "BudapestReference"=subset_samples(DT.Nitrogen, samplecodes=="BudapestReference"),
               "BudapestRuderal"=subset_samples(DT.Nitrogen, samplecodes=="BudapestRuderal"),
               "BudapestTurf"=subset_samples(DT.Nitrogen, samplecodes=="BudapestTurf"),
               "BudapestRemnant"=subset_samples(DT.Nitrogen, samplecodes=="BudapestRemnant"),
               "South AfricaReference"=subset_samples(DT.Nitrogen, samplecodes=="South AfricaReference"),
               "South AfricaRemnant"=subset_samples(DT.Nitrogen, samplecodes=="South AfricaRemnant"),
               "South AfricaTurf"=subset_samples(DT.Nitrogen, samplecodes=="South AfricaTurf"),
               "South AfricaRuderal"=subset_samples(DT.Nitrogen, samplecodes=="South AfricaRuderal"))
#NitrogenNetwork60<-netStat(List.Nitrogen, groups2, "Nitrogen") # param set to 0.6
NitrogenNetwork70<-netStat(List.Nitrogen, groups2, "Nitrogen") # param set to 0.7
#NitrogenNetwork80<-netStat(List.Nitrogen, groups2, "Nitrogen") # param set to 0.8
#NitrogenNetwork90<-netStat(List.Nitrogen, groups2, "Nitrogen") # param set to 0.9

#NitrogenNetwork500<-netStat(List.Nitrogen, groups2, "Nitrogen")
#NitrogenNetwork800<-netStat(List.Nitrogen, groups2, "Nitrogen")

# make correlation between nitrogen metabolism genes and N 

df<-data.frame(sample_data(RDT.Nitrogen)$NH4_N,sample_data(RDT.Nitrogen)$NO3_N, sample_data(RDT.Nitrogen)$C_org, sample_data(RDT.Nitrogen)$total_N)

Ncordat<-cov(otu_table(RDT.Nitrogen), df, use="complete.obs")
Ncordat<-data.frame(Ncordat)

sort(Ncordat$sample_data.RDT.Nitrogen..NH4_N, decreasing=TRUE)
sort(Ncordat$sample_data.RDT.Nitrogen..NO3_N, decreasing=TRUE)
sort(Ncordat$sample_data.RDT.Nitrogen..total_N, decreasing=TRUE)

#RDT.Nitrogen2
Ncordat2<-cov(otu_table(DT.Nitrogen2), df, use="complete.obs")
Ncordat2<-data.frame(Ncordat2)

sort(Ncordat2$sample_data.RDT.Nitrogen..NH4_N, decreasing=TRUE)
sort(Ncordat2$sample_data.RDT.Nitrogen..NO3_N, decreasing=TRUE)
sort(Ncordat2$sample_data.RDT.Nitrogen..total_N, decreasing=TRUE)

library(plyr)
table<-as.data.frame(as.matrix(otu_table(DT.Nitrogen2)))
ddply(table, linmod, y=,data=)
# bacterial network analysis ####
List.Bac<-c("BaltimoreRemnant"=subset_samples(baclib, samplecodes=="BaltimoreRemnant"),
                "BaltimoreRuderal"=subset_samples(baclib, samplecodes=="BaltimoreRuderal"),
                "BaltimoreTurf"=subset_samples(baclib, samplecodes=="BaltimoreTurf"),
                "BaltimoreReference"=subset_samples(baclib, samplecodes=="BaltimoreReference"),
                "HelsinkiReference"=subset_samples(baclib, samplecodes=="HelsinkiReference"),
                "HelsinkiRuderal"=subset_samples(baclib, samplecodes=="HelsinkiRuderal"),
                "HelsinkiTurf"=subset_samples(baclib, samplecodes=="HelsinkiTurf"),
                "HelsinkiRemnant"=subset_samples(baclib, samplecodes=="HelsinkiRemnant"),
                "LaktiReference"=subset_samples(baclib, samplecodes=="LaktiReference"),
                "LaktiRuderal"=subset_samples(baclib, samplecodes=="LaktiRuderal"),
                "LaktiTurf"=subset_samples(baclib, samplecodes=="LaktiTurf"),
                "LaktiRemnant"=subset_samples(baclib, samplecodes=="LaktiRemnant"),
                "BudapestReference"=subset_samples(baclib, samplecodes=="BudapestReference"),
                "BudapestRuderal"=subset_samples(baclib, samplecodes=="BudapestRuderal"),
                "BudapestTurf"=subset_samples(baclib, samplecodes=="BudapestTurf"),
                "BudapestRemnant"=subset_samples(baclib, samplecodes=="BudapestRemnant"),
                "South AfricaReference"=subset_samples(baclib, samplecodes=="South AfricaReference"),
                "South AfricaRemnant"=subset_samples(baclib, samplecodes=="South AfricaRemnant"),
                "South AfricaTurf"=subset_samples(baclib, samplecodes=="South AfricaTurf"),
                "South AfricaRuderal"=subset_samples(baclib, samplecodes=="South AfricaRuderal"))
#DiseaseNetwork60<-netStat(List.Disease, groups2, "Disease") # param set to 0.6
BacNetwork70<-netStat(List.Bac, groups2, "Prok") # param set to 0.7
#DiseaseNetwork80<-netStat(List.Disease, groups2, "Disease") # param set to 0.8
#DiseaseNetwork90<-netStat(List.Disease, groups2, "Disease") # param set to 0.9
List.Bac2<-c("BaltimoreRemnant"=subset_samples(baclib2, samplecodes=="BaltimoreRemnant"),
            "BaltimoreRuderal"=subset_samples(baclib2, samplecodes=="BaltimoreRuderal"),
            "BaltimoreTurf"=subset_samples(baclib2, samplecodes=="BaltimoreTurf"),
            "BaltimoreReference"=subset_samples(baclib2, samplecodes=="BaltimoreReference"),
            "HelsinkiReference"=subset_samples(baclib2, samplecodes=="HelsinkiReference"),
            "HelsinkiRuderal"=subset_samples(baclib2, samplecodes=="HelsinkiRuderal"),
            "HelsinkiTurf"=subset_samples(baclib2, samplecodes=="HelsinkiTurf"),
            "HelsinkiRemnant"=subset_samples(baclib2, samplecodes=="HelsinkiRemnant"),
            "LaktiReference"=subset_samples(baclib2, samplecodes=="LaktiReference"),
            "LaktiRuderal"=subset_samples(baclib2, samplecodes=="LaktiRuderal"),
            "LaktiTurf"=subset_samples(baclib2, samplecodes=="LaktiTurf"),
            "LaktiRemnant"=subset_samples(baclib2, samplecodes=="LaktiRemnant"),
            "BudapestReference"=subset_samples(baclib2, samplecodes=="BudapestReference"),
            "BudapestRuderal"=subset_samples(baclib2, samplecodes=="BudapestRuderal"),
            "BudapestTurf"=subset_samples(baclib2, samplecodes=="BudapestTurf"),
            "BudapestRemnant"=subset_samples(baclib2, samplecodes=="BudapestRemnant"),
            "South AfricaReference"=subset_samples(baclib2, samplecodes=="South AfricaReference"),
            "South AfricaRemnant"=subset_samples(baclib2, samplecodes=="South AfricaRemnant"),
            "South AfricaTurf"=subset_samples(baclib2, samplecodes=="South AfricaTurf"),
            "South AfricaRuderal"=subset_samples(baclib2, samplecodes=="South AfricaRuderal"))
Bac2Network70<-netStat(List.Bac2, groups2, "Prok2") # param set to 0.7
# Subset prok. Phyloseq ####

#baclib<-readRDS("~/Desktop/PhD/Metagenome/Bac16s_rarefied.RDS")
#baclib2<-readRDS("~/Desktop/PhD/Metagenome/Bac16s_raw.RDS")
#baclib2<-transform_sample_counts(baclib2, transform)




# correlate function to taxa ####
sample_data(AOA)$SumsAOA<-sample_sums(AOA)
aoadat<-sample_data(AOA)



ggplot(aoadat, aes(x=log10(AmoA_copies), y=log10(SumsAOA)))+
  geom_point(shape=1,size=4) +    # Use hollow circles
    geom_smooth(method=lm,   # Add linear regression line
                se=FALSE, color="black")+
  theme_bw()

plot(log10(sample_data(AOA)$AmoA_copies)~log10(sample_sums(AOA)))
cor(sample_sums(AOA),sample_data(AOA)$AmoA_copies)
cor(log10(sample_data(AOA)$AmoA_copies),log10(sample_sums(AOA)))

plot(log10(sample_data(AOB)$AmoB_copies)~log10(sample_sums(AOB)))#come back to this after phylogenetic analysis!
plot(sample_data(AOB)$AmoB_copies~sample_sums(AOB))

cor(sample_sums(AOB),sample_data(AOB)$AmoB_copies)
cor(log10(sample_data(AOA)$AmoA_copies), log10(sample_sums(AOA)))

plot(sample_sums(AOA)~sample_data(AOA)$Codes)
plot(sample_data(AOA)$AmoA_copies~sample_data(AOA)$Codes)
cor(sample_sums(AOA),sample_data(AOA)$AmoA_copies)
log10(sample_sums(AOA))
cor(log10(sample_sums(AOA)),log10(sample_data(AOA)$AmoA_copies),use="all.obs")

#aob<-subset_taxa(baclib, Class=="Betaproteobacteria")

# correlate function to env. ####


CoCbZnAvail<-sample_data(baclib)$Co_avail+sample_data(baclib)$Cd_avail+sample_data(baclib)$Zn_avail
CoCbZnTot<-sample_data(baclib)$Co_tot+sample_data(baclib)$Cd_tot+sample_data(baclib)$Zn_tot
plot(log(sample_sums(CoZnCdFun))~log10(CoCbZnAvail))
plot(sample_sums(CoZnCdFun)~log10(sample_data(baclib)$Co_avail))
plot(sample_sums(CoZnCdFun)~log10(sample_data(baclib)$Cd_avail))
plot(log10(sample_sums(CoZnCdFun))~log10(sample_data(baclib)$Zn_avail))
plot(sample_sums(Choleratox)~sample_data(Choleratox)$Cities)

cor(as.numeric(sample_sums(CoZnCdFun)),CoCbZnAvail, use="pairwise.complete.obs")
cor(as.numeric(sample_sums(CoZnCdFun)),CoCbZnTot, use="pairwise.complete.obs")
cor(as.numeric(sample_sums(CoZnCdFun)),sample_data(CoZnCdFun)$Cd_tot, method="pearson", use="pairwise.complete.obs")
cor(as.numeric(sample_sums(CoZnCdFun)),sample_data(CoZnCdFun)$Cd_avail, method="pearson", use="pairwise.complete.obs")
cor(as.numeric(sample_sums(CoZnCdFun)),sample_data(CoZnCdFun)$Co_tot, method="pearson", use="pairwise.complete.obs")
cor(as.numeric(sample_sums(CoZnCdFun)),sample_data(CoZnCdFun)$Co_avail, method="pearson", use="pairwise.complete.obs")
cor(as.numeric(sample_sums(CoZnCdFun)),sample_data(CoZnCdFun)$Zn_tot, method="pearson", use="pairwise.complete.obs")
cor(as.numeric(sample_sums(CoZnCdFun)),sample_data(CoZnCdFun)$Zn_avail, method="pearson", use="pairwise.complete.obs")

cor(as.numeric(sample_sums(CdFun)),sample_data(CdFun)$Cd_tot, use="pairwise.complete.obs")
cor(as.numeric(sample_sums(CdFun)),sample_data(CdFun)$Cd_avail, use="pairwise.complete.obs")

# correlation summary ####
Nammon<-subset_taxa(RDT.Nitrogen, L3=="Nitrate_and_nitrite_ammonification")
NFT.cor<-cor(as.data.frame(as.matrix(otu_table(baclib))),as.data.frame(as.matrix(otu_table(Nammon))), method = "pearson")



length(NFT.cor[NFT.cor>0.6])
length(NFT.cor[NFT.cor>0.7])
length(NFT.cor[NFT.cor>0.8])
length(NFT.cor[NFT.cor>0.9])


F.L3<-tax_glom(RDT.Nitrogen, taxrank = "L3")
F.L1<-tax_glom(RDT.Nitrogen, taxrank = "L1")
T.g<-tax_glom(N.Prok, taxrank = "Genus")
T.f<-tax_glom(N.Prok, taxrank = "Family")
T.o<-tax_glom(N.Prok, taxrank = "Order")

NFT.corF4.a<-cor(as.data.frame(as.matrix(otu_table(N.Prok))),as.data.frame(as.matrix(otu_table(RDT.Nitrogen))), method = "pearson")
NFT.corF4.g<-cor(as.data.frame(as.matrix(otu_table(T.g))),as.data.frame(as.matrix(otu_table(RDT.Nitrogen))), method = "pearson")
NFT.corF4.f<-cor(as.data.frame(as.matrix(otu_table(T.f))),as.data.frame(as.matrix(otu_table(RDT.Nitrogen))), method = "pearson")
NFT.corF4.o<-cor(as.data.frame(as.matrix(otu_table(T.o))),as.data.frame(as.matrix(otu_table(RDT.Nitrogen))), method = "pearson")

NFT.corF3.a<-cor(as.data.frame(as.matrix(otu_table(N.Prok))),as.data.frame(as.matrix(otu_table(F.L3))), method = "pearson")
NFT.corF3.g<-cor(as.data.frame(as.matrix(otu_table(T.g))),as.data.frame(as.matrix(otu_table(F.L3))), method = "pearson")
NFT.corF3.f<-cor(as.data.frame(as.matrix(otu_table(T.f))),as.data.frame(as.matrix(otu_table(F.L3))), method = "pearson")
NFT.corF3.o<-cor(as.data.frame(as.matrix(otu_table(T.o))),as.data.frame(as.matrix(otu_table(F.L3))), method = "pearson")

NFT.corF1.a<-cor(as.data.frame(as.matrix(otu_table(N.Prok))),as.data.frame(as.matrix(otu_table(F.L1))), method = "pearson")
NFT.corF1.g<-cor(as.data.frame(as.matrix(otu_table(T.g))),as.data.frame(as.matrix(otu_table(F.L1))), method = "pearson")
NFT.corF1.f<-cor(as.data.frame(as.matrix(otu_table(T.f))),as.data.frame(as.matrix(otu_table(F.L1))), method = "pearson")
NFT.corF1.o<-cor(as.data.frame(as.matrix(otu_table(T.o))),as.data.frame(as.matrix(otu_table(F.L1))), method = "pearson")

corSummaries<-matrix(nrow=4, ncol=3)
colnames(corSummaries)<-c("L1", "L3", "L4")
rownames(corSummaries)<-c("ASV", "Genus", "Family", "Order")

#number of cor values above 0.8

corSummaries[1,1]<-length(NFT.corF1.a[abs(NFT.corF1.a)>0.8])
corSummaries[1,2]<-length(NFT.corF1.a[abs(NFT.corF3.a)>0.8])
corSummaries[1,3]<-length(NFT.corF1.a[abs(NFT.corF4.a)>0.8])

corSummaries[2,1]<-length(NFT.corF1.g[abs(NFT.corF1.g)>0.8])
corSummaries[2,2]<-length(NFT.corF3.g[abs(NFT.corF3.g)>0.8])
corSummaries[2,3]<-length(NFT.corF4.g[abs(NFT.corF4.g)>0.8])

corSummaries[3,1]<-length(NFT.corF1.f[abs(NFT.corF1.f)>0.8])
corSummaries[3,2]<-length(NFT.corF3.f[abs(NFT.corF3.f)>0.8])
corSummaries[3,3]<-length(NFT.corF4.f[abs(NFT.corF4.f)>0.8])

corSummaries[4,1]<-length(NFT.corF1.o[abs(NFT.corF1.o)>0.8])
corSummaries[4,2]<-length(NFT.corF3.o[abs(NFT.corF3.o)>0.8])
corSummaries[4,3]<-length(NFT.corF4.o[abs(NFT.corF4.o)>0.8])

corSummaries

length(NFT.cor[NFT.cor=1])
# map traits to phylogeny ####
library(phytools)
library(phylosignal)
library(phylobase)
library(ape)
library(adephylo)
library(adiv)
library(picante)

# RDT.Nitrogen == relativized N function data

# correlate all bac to the relativised nitrogen data
# subset to n function
# correlate bacteria
# network analysis of depauperate community
# plot tree of putative NO3 nitrifiers, colored by phylum/order?

sample_names(RDT.Nitrogen)
sample_names(ps.rare)<-sample_data(ps.rare)$Sample_ID2
N.Prok<-as.data.frame(as.matrix(otu_table(ps.rare)))
saveRDS(N.Prok, "~/Desktop/PhD/Metagenome/NProk.rds")
N.Prok<-readRDS("~/Desktop/PhD/Metagenome/NProk.rds")

dim(N.Prok)
N.func<-as.data.frame(as.matrix(otu_table(RDT.Nitrogen)))
NFT.cor<-cor(N.Prok[rownames(N.func),], N.func, method = "pearson")
#NFT.corlog<-cor(log10(N.Prok)[rownames(N.func),], log10(N.func), method = "pearson", use="complete.obs")
NFT.cor[is.na(NFT.cor)]<-0
length(NFT.cor[abs(NFT.cor)>0.99])
dim(NFT.cor)
NFT.cor<-as.data.frame(NFT.cor)

plot_tree(ps.rare, "treeonly", ladderize="left", color="Phylum")

sample_data(ps.rare)<-c("AOA"=sample_data(ps.rare)$AmoA_copies, "AOB"=sample_data(ps.rare)$AmoB_copies)
phylocor<-cor()

phy<-phy_tree(bac)
summary(phy)
names(NFT.cor$SS00448)<-names(NFT.cor)

setwd("~/Desktop/PhD/Metagenome/Bacterial_Seqs")
trefile = ape::read.tree("phylo_tree.tre")
trefile<-ape::collapse.singles()
nwk<-readNewick("phylo_tree.tre")

phy<-makeNodeLabel(trefile, method="number", prefix="Node")
names(phy$node.label)<-phy$node.label

saveRDS(NFT.cor, "~/Desktop/PhD/Metagenome/bacfunCorMat.rds")
g1 <- as(geospiza_raw$tree, "phylo4")

testfun<-NFT.cor$SS00448
names(testfun)<-rownames(NFT.cor)
testfun[1:10]

#Kstat<-K(trefile, NFT.cor$SS00448,nrep=999, alter="greater")
#Kcalc(NFT.cor$SS00448, phy, checkdata=TRUE)

physig<-phylosignal(testfun[phy$tip.label], phy, reps=999, checkdata=TRUE)

names(trefile$node.label)<-paste("a", c(1:56170), sep="")
names(trefile$tip.label)<-trefile$tip.label
g1<-as(trefile, "phylo4")# doesn't work, makes error

g2<-as(phy, "phylo4")
tips<-as.vector(taxa_names(ps.rare))
p4d$tree<-phy
p4d<-phylo4d(as(phy, "phylo4"), data.frame(NFT.cor))
p4d$dat<-
phylosig(trefile, NFT.cor$SS00448, method="lambda", test=TRUE)
phyloSignal(tree, methods="lambda", reps=999)


# convergence testing Networks ####
library(car)

leveneTest(StressNetwork70$metrics$Modularity, group=StressNetwork70$metrics$Codes)
leveneTest(StressNetwork70$metrics$Nodes, group=StressNetwork70$metrics$Codes)
leveneTest(StressNetwork70$metrics$Degree, group=StressNetwork70$metrics$Codes)
leveneTest(StressNetwork70$metrics$Closeness, group=StressNetwork70$metrics$Codes)
leveneTest(StressNetwork70$metrics$Edges, group=StressNetwork70$metrics$Codes)


leveneTest(NitrogenNetwork70$metrics$Modularity, group=NitrogenNetwork70$metrics$Codes)
leveneTest(NitrogenNetwork70$metrics$Nodes, group=NitrogenNetwork70$metrics$Codes)
leveneTest(NitrogenNetwork70$metrics$Degree, group=NitrogenNetwork70$metrics$Codes)
leveneTest(NitrogenNetwork70$metrics$Closeness, group=NitrogenNetwork70$metrics$Codes)
leveneTest(NitrogenNetwork70$metrics$Edges, group=NitrogenNetwork70$metrics$Codes)

leveneTest(DiseaseNetwork70$metrics$Modularity, group=DiseaseNetwork70$metrics$Codes)
leveneTest(DiseaseNetwork70$metrics$Nodes, group=DiseaseNetwork70$metrics$Codes)
leveneTest(DiseaseNetwork70$metrics$Degree, group=DiseaseNetwork70$metrics$Codes)
leveneTest(DiseaseNetwork70$metrics$Closeness, group=DiseaseNetwork70$metrics$Codes)
leveneTest(DiseaseNetwork70$metrics$Edges, group=DiseaseNetwork70$metrics$Codes)

leveneTest(BacNetwork70$metrics$Modularity, group=BacNetwork70$metrics$Codes)
leveneTest(BacNetwork70$metrics$Nodes, group=BacNetwork70$metrics$Codes)
leveneTest(BacNetwork70$metrics$Degree, group=BacNetwork70$metrics$Codes)
leveneTest(BacNetwork70$metrics$Closeness, group=BacNetwork70$metrics$Codes)
leveneTest(BacNetwork70$metrics$Edges, group=BacNetwork70$metrics$Codes)

# bac 2 uses the wrong normalization approach
leveneTest(Bac2Network70$metrics$Modularity, group=Bac2Network70$metrics$Codes)
leveneTest(Bac2Network70$metrics$Nodes, group=Bac2Network70$metrics$Codes)
leveneTest(Bac2Network70$metrics$Degree, group=Bac2Network70$metrics$Codes)
leveneTest(Bac2Network70$metrics$Closeness, group=Bac2Network70$metrics$Codes)
leveneTest(Bac2Network70$metrics$Edges, group=Bac2Network70$metrics$Codes)


# N cycle table calculations ####
#d<-data.frame("gene"=unlist(sample_sums(subset_taxa(rdt15.l4, L4=="Copper-containing nitrite reductase (EC 1.7.2.1)"))), "Land"=as.factor(sample_data(rdt15.l4)$Codes))
#ddply(d, .(Land), summarise, 
#      median=median(gene),
 #     mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Copper-containing nitrite reductase (EC 1.7.2.1)")), "Land"=as.factor(sample_data(rdt15.l4)$Codes)), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Cytochrome cd1 nitrite reductase (EC:1.7.2.1)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Cytochrome c nitrite reductase, small subunit NrfH")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Ferredoxin--nitrite reductase (EC 1.7.7.1)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrite reductase [NAD(P)H] large subunit (EC 1.7.1.4)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrite reductase [NAD(P)H] small subunit (EC 1.7.1.4)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrite reductase probable [NAD(P)H] subunit (EC 1.7.1.4)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitric-oxide reductase (EC 1.7.99.7), quinol-dependent")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitric-oxide reductase subunit B (EC 1.7.99.7)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise,
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitric-oxide reductase subunit C (EC 1.7.99.7)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise,
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitric oxide -responding transcriptional regulator NnrR (Crp/Fnr family)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrous-oxide reductase (EC 1.7.99.6)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrous oxide reductase maturation protein NosD")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrous oxide reductase maturation protein NosF (ATPase)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrous oxide reductase maturation protein NosR")), "Land"=sample_data(rdt15.l4)$Codes), .(Land),summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrous oxide reductase maturation protein, outer-membrane lipoprotein NosL")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrous oxide reductase maturation transmembrane protein NosY")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitric oxide reductase activation protein NorD")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitric oxide reductase activation protein NorQ")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Cytochrome c-type heme lyase subunit nrfE, nitrite reductase complex assembly")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Cytochrome c-type heme lyase subunit nrfF, nitrite reductase complex assembly")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Cytochrome c-type heme lyase subunit nrfG, nitrite reductase complex assembly")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrite reductase associated c-type cytochorome NirN")), "Land"=sample_data(rdt15.l4)$Codes), .(Land),summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Ammonia monooxygenase")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrogenase FeMo-cofactor scaffold and assembly protein NifE")), "Land"=sample_data(rdt15.l4)$Codes), .(Land),summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrogenase FeMo-cofactor scaffold and assembly protein NifN")), "Land"=sample_data(rdt15.l4)$Codes), .(Land),summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrogenase FeMo-cofactor synthesis FeS core scaffold and assembly protein NifB")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrogenase (molybdenum-iron) alpha chain (EC 1.18.6.1)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrogenase (molybdenum-iron) beta chain (EC 1.18.6.1)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrogenase (molybdenum-iron) reductase and maturation protein NifH")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nitrogenase (molybdenum-iron)-specific transcriptional regulator NifA")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))




ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Additional substrate-specific component NikN of nickel ECF transporter")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="ATPase component NikO of energizing module of nickel ECF transporter")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Substrate-specific component NikM of nickel ECF transporter")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Transmembrane component NikQ of energizing module of nickel ECF transporter")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Additional substrate-specific component NikN of nickel ECF transporter")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="ATPase component NikO of energizing module of nickel ECF transporter")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Substrate-specific component NikM of nickel ECF transporter")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Transmembrane component NikQ of energizing module of nickel ECF transporter")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))


ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel ABC transporter, periplasmic nickel-binding protein nikA2 (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel ABC transporter, periplasmic nickel-binding protein NikA (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel transport ATP-binding protein nikD2 (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel transport ATP-binding protein NikD (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel transport ATP-binding protein nikE2 (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel transport ATP-binding protein NikE (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel transport system permease protein nikB2 (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel transport system permease protein NikB (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel transport system permease protein nikC2 (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Nickel transport system permease protein NikC (TC 3.A.1.5.3)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))


ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Copper resistance protein B")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Copper resistance protein D")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Copper tolerance protein")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Membrane protein, suppressor for copper-sensitivity ScsB")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Membrane protein, suppressor for copper-sensitivity ScsD")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Secreted protein, suppressor for copper-sensitivity ScsC")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Suppression of copper sensitivity: putative copper binding protein ScsA")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))


ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multiple antibiotic resistance protein MarA")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multiple antibiotic resistance protein MarA")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multiple antibiotic resistance protein MarC")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multiple antibiotic resistance protein MarR")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multidrug resistance transporter, Bcr/CflA family")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multidrug efflux transporter MexF")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Membrane fusion protein of RND family multidrug efflux pump")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multi antimicrobial extrusion protein (Na(+)/drug antiporter), MATE family of MDR efflux pumps")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multidrug efflux RND transporter MexD")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multidrug-efflux transporter, major facilitator superfamily (MFS) (TC 2.A.1)")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multidrug resistance efflux pump PmrA")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multidrug transporter MdtB")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multidrug transporter MdtC")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Multidrug transporter MdtD")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="MFS family multidrug efflux protein, similarity to bicyclomycin resistance protein Bcr")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="ABC transporter multidrug efflux pump, fused ATP-binding domains")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="ABC-type multidrug transport system, ATPase component")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="ABC-type multidrug transport system, permease component")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(prokOrder, Family=="Nitrososphaeraceae")), "Land"=sample_data(prokOrder)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))

ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="Ni,Fe-hydrogenase I cytochrome b subunit")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
ddply(data.frame("gene"=sample_sums(subset_taxa(rdt15.l4, L4=="[NiFe] hydrogenase nickel incorporation protein HybF")), "Land"=sample_data(rdt15.l4)$Codes), .(Land), summarise, 
      median=median(gene),
      mean=mean(gene))
