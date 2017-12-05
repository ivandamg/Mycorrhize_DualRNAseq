

#GOseq analysis script.
# Do GENe ontology analysis
# need to import list of universe and differentially expressed genes.
####################################################################


#libraries
####################################################################

#install.packages(c("devtools","MatrixEQTL"))
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("Biobase","goseq","DESeq2"))

library(devtools)
library(Biobase)
library(goseq)
library(DESeq2)



#annotation info
####################################################################





#annotation info
####################################################################
#M esculenta
mesculenta_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/Mesculenta_305_v6.1.annotation_info.txt',h=T)
#R. irregularis
rirregularis_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/blast2go_annot_predicted_prot_hint_glomus_nu6_genome_masked.annot',h=F)
colnames(rirregularis_go)<-c('locusName','GO.ID')
####################################################################
GOterms<-mesculenta_go


### create transcript length
library(GenomicRanges)
library(rtracklayer)
GTFfile = "~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/mesculenta_305_v6.1.gene_exons.gff3"
GTF <- import.gff(GTFfile, format="gff3", genome="Manes.61",  feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$ID))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$ID <- rep(names(grl), elementLengths(grl))
elementMetadata(reducedGTF)$widths <- width(reducedGTF)
calc_length <- function(x) {sum(elementMetadata(x)$widths)}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$ID), calc_length))
colnames(output) <- c("Length")
output2<-as.vector(output)
output
#write.csv(output, "~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/temp_TranscriptLenths_new.csv")

# fin test.

tail(sign_YUCA)
sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
rownames(sign_YUCA)

ALL <- rownames(sign_YUCA)
DEG <- rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])
Transcriptlength <- output #read.csv("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/temp_TranscriptLenths_new.csv", header=FALSE)

DEG.vector <- c(t(DEG))
ALL.vector<-c(t(ALL))
gene.vector=as.integer(ALL.vector%in%DEG.vector)
names(gene.vector)=ALL.vector
head(gene.vector)
Length.vector <- c(t(Transcriptlength))
pwf=nullp(gene.vector,bias.data=Length.vector)


GO.wall=goseq(pwf,gene2cat=GOterms,use_genes_without_cat=TRUE)


GO.samp=goseq(pwf,gene2cat=GOterms,use_genes_without_cat=TRUE,method="Sampling",repcnt=1000)

plot(log10(GO.wall[,2]), log10(GO.samp[match(GO.samp[,1],GO.wall[,1]),2]),xlab="log10(Wallenius p-values)", ylab="log10(Sampling p-values)", xlim=c(-3,0), )
abline(0,1,col=3,lty=2)





GO.nobias=goseq(pwf,gene2cat=GOterms,use_genes_without_cat=TRUE,method="Hypergeometric")


enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method="BH")<.05]



write.csv(GO.wall,"GO.wall.csv")
write.csv(GO.samp,"GO.samp,csv")
write.csv(GO.nobias,"GO.nobias.csv")
write.csv(enriched.GO,"enriched.GO.csv")


source("http://bioconductor.org/biocLite.R")
bioclite("GO.db")
library(GO.db)
for(go in enriched.GO[1:10]){print(GOTERM[[go]] cat("-------—-\n")}

