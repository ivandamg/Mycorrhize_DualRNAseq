########################################################################################
#' Import transcript-level abundances and estimated counts for gene-level analysis packages
########################################################################################
#source("http://bioconductor.org/biocLite.R")
#biocLite("Rgraphviz")
#install.packages("pheatmap")
#biocLite("topGO")
#biocLite("limma")
#biocLite("edgeR")
library("gplots")
library("RColorBrewer")
library(Rgraphviz)
library(topGO)
library(edgeR)
library(limma)
library(heatmap2)
library(WGCNA)
library(pheatmap)
options("scipen"=3,"digits"=2)
#biocLite(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
#install.packages("~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/WGCNA_1.49.tgz",type = "source", repos = NULL, lib=.Library) 

########################################################################################
# ANALYSIS AT GENE LEVEL 
# input data: Kallisto pseudo-alignment. 

tximport <- function(files,
                     type=c("none","kallisto","salmon","sailfish","rsem"),
                     txIn=TRUE,
                     txOut=FALSE,
                     countsFromAbundance=c("no","scaledTPM","lengthScaledTPM"),
                     tx2gene=NULL,
                     reader=read.delim,
                     geneIdCol,
                     txIdCol,
                     abundanceCol,
                     countsCol,
                     lengthCol,
                     importer,
                     collatedFiles,
                     ignoreTxVersion=FALSE) {
  
  type <- match.arg(type, c("none","kallisto","salmon","sailfish","rsem"))
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  stopifnot(all(file.exists(files)))
  
  # kallisto presets
  if (type == "kallisto") {
    geneIdCol="gene_id"
    txIdCol <- "target_id"
    abundanceCol <- "tpm"
    countsCol <- "est_counts"
    lengthCol <- "eff_length"
    importer <- reader
  }
  
  # salmon/sailfish presets
  if (type %in% c("salmon","sailfish")) {
    geneIdCol="gene_id"
    txIdCol <- "Name"
    abundanceCol <- "TPM"
    countsCol <- "NumReads"
    lengthCol <- "EffectiveLength"
    importer <- function(x) reader(x, comment='#') 
  }
  
  # rsem presets
  if (type == "rsem") {
    txIn <- FALSE
    geneIdCol <- "gene_id"
    abundanceCol <- "FPKM"
    countsCol <- "expected_count"
    lengthCol <- "effective_length"
    importer <- reader
  }
  
  if (type == "cufflinks") {
    stop("reading from collated files not yet implemented")
  }
  
  # if input is tx-level, need to summarize abundances, counts and lengths to gene-level
  if (txIn) {
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      
      #####################################################################
      # some temporary code for detecting older fishes
      if ((i == 1) &
          (type %in% c("salmon","sailfish")) &
          !("EffectiveLength" %in% names(raw))) {
        lengthCol <- "Length" 
        # because the comment lines have the same comment character
        # as the header, need to name the column names
        importer <- function(x) {
          tmp <- reader(x, comment="#")
          names(tmp) <- c("Name","Length","TPM","NumReads")
          tmp
        }
        # re-read the first file
        raw <- as.data.frame(importer(files[i]))
      }
      #####################################################################
      
      # does the table contain gene association or was an external tx2gene table provided?
      if (is.null(tx2gene) & !txOut) {
        # e.g. Cufflinks includes the gene ID in the table
        stopifnot(all(c(geneIdCol, lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          geneId <- raw[[geneIdCol]]
        } else {
          stopifnot(all(geneId == raw[[geneIdCol]]))
        }
      } else {
        # e.g. Salmon and kallisto do not include the gene ID, need an external table
        stopifnot(all(c(lengthCol, abundanceCol) %in% names(raw)))
        if (i == 1) {
          txId <- raw[[txIdCol]]
        } else {
          stopifnot(all(txId == raw[[txIdCol]]))
        }
      }
      # create empty matrices
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[txIdCol]]
        colnames(mat) <- names(files)
        abundanceMatTx <- mat
        countsMatTx <- mat
        lengthMatTx <- mat
      }
      abundanceMatTx[,i] <- raw[[abundanceCol]]
      countsMatTx[,i] <- raw[[countsCol]]
      lengthMatTx[,i] <- raw[[lengthCol]]
    }
    message("")
    
    txi <- list(abundance=abundanceMatTx, counts=countsMatTx, length=lengthMatTx,
                countsFromAbundance="no")
    
    # if the user requested just the transcript-level data:
    if (txOut) {
      return(txi)
    }
    
    txi[["countsFromAbundance"]] <- NULL
    txiGene <- summarizeToGene(txi, tx2gene, ignoreTxVersion, countsFromAbundance)
    return(txiGene)  
    
    # e.g. RSEM already has gene-level summaries
    # just combine the gene-level summaries across files
  } else {
    # stating the obvious:
    if (txOut) stop("txOut only an option when transcript-level data is read in (txIn=TRUE)")
    
    message("reading in files")
    for (i in seq_along(files)) {
      message(i," ",appendLF=FALSE)
      raw <- as.data.frame(importer(files[i]))
      stopifnot(all(c(geneIdCol, abundanceCol, lengthCol) %in% names(raw)))
      if (i == 1) {
        mat <- matrix(nrow=nrow(raw),ncol=length(files))
        rownames(mat) <- raw[[geneIdCol]]
        colnames(mat) <- names(files)
        abundanceMat <- mat
        countsMat <- mat
        lengthMat <- mat
      }
      abundanceMat[,i] <- raw[[abundanceCol]]
      countsMat[,i] <- raw[[countsCol]]
      lengthMat[,i] <- raw[[lengthCol]]
    }
  } 
  message("")
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance="no"))
}

# summarizeToGene() splits out the summarization functions
# in tximport(), so it can be called by users to summarize
# transcript-level lists of matrices

#' @describeIn tximport Summarize tx-level matrices to gene-level
#' @export
summarizeToGene <- function(txi,
                            tx2gene,
                            ignoreTxVersion=FALSE,
                            countsFromAbundance=c("no","scaledTPM","lengthScaledTPM")
) {
  
  countsFromAbundance <- match.arg(countsFromAbundance, c("no","scaledTPM","lengthScaledTPM"))
  
  # unpack matrices from list for cleaner code
  abundanceMatTx <- txi$abundance
  countsMatTx <- txi$counts
  lengthMatTx <- txi$length
  
  txId <- rownames(abundanceMatTx)
  stopifnot(all(txId == rownames(countsMatTx)))
  stopifnot(all(txId == rownames(lengthMatTx)))
  
  # need to associate tx to genes
  # potentially remove unassociated transcript rows and warn user
  if (!is.null(tx2gene)) {
    colnames(tx2gene) <- c("tx","gene")
    if (ignoreTxVersion) {
      txId <- sapply(strsplit(as.character(txId), "\\."), "[[", 1)
    }
    tx2gene$gene <- factor(tx2gene$gene)
    tx2gene$tx <- factor(tx2gene$tx)
    # remove transcripts (and genes) not in the abundances
    tx2gene <- tx2gene[tx2gene$tx %in% txId,]
    tx2gene$gene <- droplevels(tx2gene$gene)
    ntxmissing <- sum(!txId %in% tx2gene$tx)
    if (ntxmissing > 0) message("transcripts missing genes: ", ntxmissing)
    sub.idx <- txId %in% tx2gene$tx
    abundanceMatTx <- abundanceMatTx[sub.idx,,drop=FALSE]
    countsMatTx <- countsMatTx[sub.idx,,drop=FALSE]
    lengthMatTx <- lengthMatTx[sub.idx,,drop=FALSE]
    txId <- txId[sub.idx]
    geneId <- tx2gene$gene[match(txId, tx2gene$tx)]
  }
  
  # summarize abundance and counts
  message("summarizing abundance")
  abundanceMat <- fastby(abundanceMatTx, geneId, colSums)
  message("summarizing counts")
  countsMat <- fastby(countsMatTx, geneId, colSums)
  message("summarizing length")
  
  # the next lines calculate a weighted average of transcript length, 
  # weighting by transcript abundance.
  # this can be used as an offset / normalization factor which removes length bias
  # for the differential analysis of estimated counts summarized at the gene level.
  weightedLength <- fastby(abundanceMatTx * lengthMatTx, geneId, colSums)
  lengthMat <- weightedLength / abundanceMat   
  
  # pre-calculate a simple average transcript length
  # for the case the abundances are all zero for all samples.
  # first, average the tx lengths over samples
  aveLengthSamp <- rowMeans(lengthMatTx)
  # then simple average of lengths within genes (not weighted by abundance)
  aveLengthSampGene <- tapply(aveLengthSamp, geneId, mean)
  
  stopifnot(all(names(aveLengthSampGene) == rownames(lengthMat)))
  
  # check for NaN and if possible replace these values with geometric mean of other samples.
  # (the geometic mean here implies an offset of 0 on the log scale)
  # NaN come from samples which have abundance of 0 for all isoforms of a gene, and 
  # so we cannot calculate the weighted average. our best guess is to use the average
  # transcript length from the other samples.
  lengthMat <- replaceMissingLength(lengthMat, aveLengthSampGene)
  
  if (countsFromAbundance != "no") {
    countsSum <- colSums(countsMat)
    if (countsFromAbundance == "lengthScaledTPM") {
      newCounts <- abundanceMat * rowMeans(lengthMat)
    } else {
      newCounts <- abundanceMat
    }
    newSum <- colSums(newCounts)
    countsMat <- t(t(newCounts) * (countsSum/newSum))
  }
  
  return(list(abundance=abundanceMat, counts=countsMat, length=lengthMat,
              countsFromAbundance=countsFromAbundance))
}

# this is much faster than by(), a bit slower than dplyr summarize_each()
fastby <- function(m, f, fun) {
  idx <- split(1:nrow(m), f)
  if (ncol(m) > 1) {
    t(sapply(idx, function(i) fun(m[i,,drop=FALSE])))
  } else {
    matrix(sapply(idx, function(i) fun(m[i,,drop=FALSE])),
           dimnames=list(levels(f), colnames(m)))
  }
}

# function for replacing missing average transcript length values
replaceMissingLength <- function(lengthMat, aveLengthSampGene) {
  nanRows <- which(apply(lengthMat, 1, function(row) any(is.nan(row))))
  if (length(nanRows) > 0) {
    for (i in nanRows) {
      if (all(is.nan(lengthMat[i,]))) {
        # if all samples have 0 abundances for all tx, use the simple average
        lengthMat[i,] <- aveLengthSampGene[i]
      } else {
        # otherwise use the geometric mean of the lengths from the other samples
        idx <- is.nan(lengthMat[i,])
        lengthMat[i,idx] <-  exp(mean(log(lengthMat[i,!idx]), na.rm=TRUE))
      }
    }
  }
  lengthMat
}

#####################################################################################################################################################
#####################################################################################################################################################
# transfrom transcripts to genes & import data
########################################################################################################################################
#Mesculenta
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta') 

filesToProcess <- dir(pattern = "*_abundance.tsv$")  #files to process.
names(filesToProcess)<-gsub('_abundance.tsv$','',filesToProcess,perl=TRUE)
samples<-read.table('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/Mesculenta_V1_10_B1_abundance.tsv',h=T)

tx2gene<-cbind.data.frame(samples$target_id,gsub('.[0-9].v6.1$','',samples$target_id,perl=TRUE))
colnames(tx2gene)<-c('TXNAME','GENEID')
txi_YUCA <- tximport(filesToProcess, type="kallisto", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length
FOUR_VARS_YUCA<-txi_YUCA[[2]]
colnames(FOUR_VARS_YUCA)<-gsub("Mesculenta_","",colnames(FOUR_VARS_YUCA),perl=T)

#Rirregularis
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Umap_Rirregularis')   
filesToProcess <- dir(pattern = "*_abundance.tsv$")  #files to process.
names(filesToProcess)<-gsub('_abundance.tsv$','',filesToProcess,perl=TRUE)
samples<-read.table('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Umap_Rirregularis/UmapAMF_V1_10_B1_abundance.tsv',h=T)

tx2gene<-cbind.data.frame(samples$target_id,gsub('.[0-9].v6.1$','',samples$target_id,perl=TRUE))
colnames(tx2gene)<-c('TXNAME','GENEID')

txi_AMF <- tximport(filesToProcess, type="kallisto", tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM") # normalized library size and transcript length
FOUR_VARS_AMF<-txi_AMF[[2]]
colnames(FOUR_VARS_AMF)<-gsub("UmapAMF_","",colnames(FOUR_VARS_AMF),perl=T)
colnames(FOUR_VARS_AMF)


########################################################################################################################################
########################################################################################################################################


########################################################################################################################################
#annotation infos
########################################################################################################################################
#M esculenta
mesculenta_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/Mesculenta_305_v6.1.annotation_info.txt',h=T)
mesculenta_go2<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/annotations_blast2go_mesculenta_FMv2/blast2go_FM_v3',h=T)
colnames(mesculenta_go2)<-c("locusName","transcriptName","GO","DEF")
#R. irregularis
rirregularis_go<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Bioinformatics/blast2go_annot_predicted_prot_hint_glomus_nu6_genome_masked.annot',h=F)
colnames(rirregularis_go)<-c('locusName','GO.ID')
head(mesculenta_go)
length(unique(mesculenta_go$locusName))
length(unique(rirregularis_go$locusName))

######################################################
#' 1) LOAD MERCATOR DATABASE
######################################################
# CASSAVA
#Mercator_Mesculenta<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Mesculenta_database_v2",h=F,quote="\'")
#colnames(Mercator_Mesculenta)<-c("Bin","Function","gene")
# with detailed info
Mercator_Mesculenta<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Mesculenta_database_v4.txt",h=F,quote="\'")
#replace 'NA'	'NA'	(\w+)\"(\b)   (\w+)\'(\b)   (\w+)\'\, (\w+)\'-  (\w+)\'\( (\w+)\'\) \.\'(\w+)   (\w+)\'\] \(\'(\w+)
colnames(Mercator_Mesculenta)<-c("Bin","Function","gene","details","Type")
head(Mercator_Mesculenta)
# AMF
#Mercator_Rirregularis<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Rirregularis_database2",h=F,quote="")
#colnames(Mercator_Rirregularis)<-c("Bin","Function","gene")

# with detailed info
Mercator_Rirregularis<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Rirregularis_database_v4",h=F,quote="\'")
colnames(Mercator_Rirregularis)<-c("Bin","Function","gene","details")



########################################################################################################################################
########################################################################################################################################
# SEQUENCING RESULTS ANALYSIS
########################################################################################################################################
#seq<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Sequencing_results_V1.txt",h=T)

#rownames(seq)<-seq$Samples
#colnames(FOUR_VARS_YUCA)
#rownames(seq)
# take out bad samples
#seq<-seq[grep(c("V1-4"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V1-8"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V1-12"), rownames(seq),invert =T),]
#
#seq<-seq[grep(c("V4-16"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V48"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V58"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V88b"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V8-8"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V8-6"), rownames(seq),invert =T),]
#
#seq<-seq[grep(c("V5-9"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V5-4"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V5-18"), rownames(seq),invert =T),]
#seq<-seq[grep(c("V518"), rownames(seq),invert =T),]
#
#seq<-seq[grep(c("V56"), rownames(seq),invert =T),]
#
#seq<-seq[seq$Variety!="V3",]
#seq<-seq[seq$Treatment!="CANB1",]
#
#seq<-seq[order(seq[,2],seq[,3]),]
#seq<-droplevels(seq)
#
#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Control1_TPM100.pdf")
#par(mfrow=c(2,2))
#  boxplot(seq$tpm10_M~seq$Treatment,ylab="Transcripts > 10 counts",las=2,ylim=c(0,14000))
#  boxplot(seq$tpm10_R~seq$Treatment,ylab="Transcripts > 10 counts",las=2)
#boxplot(seq$tpm10_M~seq$Cultivar,ylab="Transcripts > 10 counts",las=2,ylim=c(0,14000))
#boxplot(seq$tpm10_R~seq$Cultivar,ylab="Transcripts > 10 counts",las=2)
#dev.off()
#
#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Control2_STAR_KALLISTO_ALLIGN.pdf")
#par(mfrow=c(1,2))
#boxplot(seq$STAR_K_Mesc~seq$Treatment,ylab="Aligned reads by STAR / by kallisto",las=2,ylim=c(0,1.4) )
#boxplot(seq$STAR_K_Mesc~seq$Cultivar,ylab="Aligned reads by STAR / by kallisto",las=2,ylim=c(0,1.4))
#dev.off()
#
#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Control3_SNP_KB.pdf")
#par(mfrow=c(1,2))
#boxplot(seq$SNP_KB_Mesc~seq$Treatment,ylab="SNP / KB",las=2,ylim=c(0,6))
#boxplot(seq$SNP_KB_Mesc~seq$Cultivar,ylab="SNP / KB",las=2,ylim=c(0,6))
#dev.off()
#
#
#
#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Control4_TPM100_samples.pdf")
#colos<-as.vector(seq$Treatment)
#colos[colos=="CAN"]<-"lightsalmon3"
#colos[colos=="B1"]<-"dodgerblue3"
#colos[colos=="CTRL"]<-"azure3"
#par(mfrow=c(2,1))
# barplot(seq$tpm10_M, names=seq$Samples,las=2,col=colos,
#        ylab="Transcripts > 10 counts")
# barplot(seq$tpm10_R, names=seq$Samples,las=2,col=colos,
#        ylab="Transcripts > 10 counts")
#dev.off()
#
#boxplot(seq$tpm_100_R~interaction(seq$Treatment,seq$Cultivar),las=2,ylim=c(0,3000))
#boxplot(seq$tpm_100_R[seq$Treatment=="CAN"]~seq$Cultivar[seq$Treatment=="CAN"],las=2,cex.axis=0.6,ylim=c(0,3000))
#boxplot(seq$tpm_100_R[seq$Treatment=="B1"]~seq$Cultivar[seq$Treatment=="B1"],las=2,cex.axis=0.6,ylim=c(0,3000))
#boxplot(seq$tpm100_M~interaction(seq$Treatment,seq$Cultivar),las=2,ylim=c(0,3000))
#boxplot(seq$tpm100_M[seq$Treatment=="CAN"]~seq$Cultivar[seq$Treatment=="CAN"],las=2,cex.axis=0.6,ylim=c(0,3000))
#boxplot(seq$tpm100_M[seq$Treatment=="B1"]~seq$Cultivar[seq$Treatment=="B1"],las=2,cex.axis=0.6,ylim=c(0,3000))
#
#
#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Control5_reads_aligned_AMF.pdf")
#barplot(seq$pseudoaligned_R,names=seq$Samples,las=2,col=colos,
#        ylab="Reads aligned")
#dev.off()
#
#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Control6_reads_align_CASSAVA.pdf")
#barplot(seq$pseudoalign_M,names=seq$Samples,las=2,col=seq$Treatment,
#        ylab="Reads aligned")
#dev.off()
#
#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Control7_ratio_reads_AMF_CASSAVA.pdf")
#barplot((seq$pseudoaligned_R/seq$pseudoalign_M),names=seq$Samples,las=2,col=seq$Treatment,ylim=c(0,1),
#        ylab="Reads aligned AMF / cassava")
#dev.off()
#
#
##ACCUMULATION CURVES
#
#ratio_R_M<-seq$pseudoaligned_R/seq$pseudoalign_M
#plot(seq$tpm10_R~ratio_R_M,las=2,col=seq$Treatment,pch=16,
#        ylab="pseudoread_align")
#
#
#pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Control8_seq_effort.pdf")
#par(mfrow=c(2,1))
#options(scipen = 10)
#plot(seq$tpm100_M[seq$Treatment=="CTRL"]~seq$pseudoalign_M[seq$Treatment=="CTRL"], ylim=c(0,3000), pch=17,col="azure3",ylab="Transcripts > 10 counts",xlab="Reads",las=1)
#points(seq$tpm100_M[seq$Treatment=="CAN"]~seq$pseudoalign_M[seq$Treatment=="CAN"], ylim=c(0,3000), pch=16,col="lightsalmon3")
#points(seq$tpm100_M[seq$Treatment=="B1"]~seq$pseudoalign_M[seq$Treatment=="B1"], ylim=c(0,3000), pch=15,col="dodgerblue3")
#
#plot(seq$tpm_100_R[seq$Treatment=="CTRL"] ~seq$pseudoaligned_R[seq$Treatment=="CTRL"], ylim=c(0,2500),xlim=c(0,8000000),  pch=17,col="azure3",ylab="Transcripts > 10 counts",xlab="Reads",las=1)
#points(seq$tpm_100_R[seq$Treatment=="CAN"] ~seq$pseudoaligned_R[seq$Treatment=="CAN"],  pch=16,col="lightsalmon3")
#points(seq$tpm_100_R[seq$Treatment=="B1"] ~seq$pseudoaligned_R[seq$Treatment=="B1"], ylim=c(0,2000),  pch=15,col="dodgerblue3")
#dev.off()
#
########################################################################################################################################

########################################################################################################################################
########################################################################################################################################
# DATA SELECTION
########################################################################################################################################
# Selection all but co-inoculation
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("CANB1"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("CANB1"), colnames(FOUR_VARS_AMF),invert =T)]

# takw out bad variety
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V3"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3"), colnames(FOUR_VARS_AMF),invert =T)]


########################################################################################################################################
# Design DESEQUILIBRE (3 replicates each) V1 (2 ctrl),V4 (2ctrl),V5,V6, V8 2 cTRL
########################################################################################################################################

# take out bad samples
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V1_4"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V1_12"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V1_8"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V4_16b"), colnames(FOUR_VARS_YUCA),invert =T)]

FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V4_8"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_8"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8_8b"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8_8"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V8_6"), colnames(FOUR_VARS_YUCA),invert =T)]

FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_9"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_4"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_18"), colnames(FOUR_VARS_YUCA),invert =T)]
FOUR_VARS_YUCA<-FOUR_VARS_YUCA[,grep(c("V5_6"), colnames(FOUR_VARS_YUCA),invert =T)]


#take out controls for AMF
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("CTRL"), colnames(FOUR_VARS_AMF),invert =T)]
# take out samples

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_6"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_9"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_18"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_6"), colnames(FOUR_VARS_AMF),invert =T)]



colnames(FOUR_VARS_YUCA[,grep(c("V8"), colnames(FOUR_VARS_YUCA),invert =F)])
colnames(FOUR_VARS_AMF[,grep(c("V"), colnames(FOUR_VARS_AMF),invert =F)])

########################################################################################################################################




########################################################################################################################################
########################################################################################################################################
# DATA FILTERING AND NORMALIZATION
#################################################################################
# DATA FILTERING AND NORMALIZATION
DGE_YUCA <- DGEList(FOUR_VARS_YUCA)
DGE_AMF <- DGEList(FOUR_VARS_AMF)

#FILTERING  

keep_Y <- rowSums(DGE_YUCA$counts>100) >= 3
DGE_YUCA<- DGE_YUCA[keep_Y,]

keep_A <- rowSums(DGE_AMF$counts>100) >= 3
DGE_AMF<- DGE_AMF[keep_A,]

#NORMALIZATION
DGE_YUCA <- calcNormFactors(DGE_YUCA)
DGE_AMF <- calcNormFactors(DGE_AMF)


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################


###################################################################################################################################################
##SUPPLEMENTARY INFO ABOUT RNA-SEQ.
###################################################################################################################################################

########################################################################################################################################
#DATA VISUALIZATION NORMALIZED
########################################################################################################################################
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.SEQ2_MDS_2organism.pdf")
par(mfrow=c(2,2))
 
plotMDS(DGE_YUCA_N,col=TREAT_YUCA,pch=16 ) #by TREAT
#plotMDS(DGE_YUCA_N,pch=16,cex.lab=3,cex.axis=3,labels=colnames(DGE_YUCA_N) ) 
#by TREAT
#legend("topleft",c('CTRL','CAN','B1'),col=c(1,2,3),pch=16)
plotMDS(DGE_YUCA_N,col=as.numeric(gsub("V","",VAR,perl=TRUE)),pch=16) #by VAR
#legend("topleft",c('V1','V3','V4','V5','V6','V8'),col=c(1,3,4,5,6,8),pch=16)

plotMDS(DGE_AMF_N,col=TREAT_AMF,pch=16) #by TREAT
#legend("topleft",c('CAN','B1'),col=c(1,2),pch=16)
plotMDS(DGE_AMF_N,col=as.numeric(gsub("V","",VAR,perl=TRUE)),pch=16) #by VAR
#legend("topleft",c('V1','V3','V4','V5','V6','V8'),col=c(1,3,4,5,6,8),pch=16)
dev.off()
#
#plotMDS(DGE_AMF_N,col=TREAT_COL,pch=16,main="Rirregularis samples differentiation",cex=3) #by TREAT
#legend("topleft",c('CTRL','CAN','B1'),col=c(1,2,3),pch=16)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.SEQ3_HEATMAP_norm_samp.pdf")
for (vari in c("V1","V4","V5","V6","V8")) {
sampleDistMatrix <- as.matrix( dist(t(DGE_YUCA_N$E[,grep(vari, colnames(DGE_YUCA_N$E),invert =F)]) ))
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(sampleDistMatrix, fontsize=18, col=colors)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
sampleDistMatrix <- as.matrix( dist(t(DGE_AMF_N$E[,grep(vari, colnames(DGE_AMF_N$E),invert =F)]) ))
pheatmap(sampleDistMatrix, fontsize=18, col=colors)

}
dev.off()  
########################################################################################################################################
# Correlation replicate 1 vs 2. for each treatment.
#accumulaticon curve
###################################################################################################################################################
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.SEQ4_CORsamples.pdf")
par(mfrow=c(3,3))
for (vari in c("V1","V4","V5","V6","V8")){
  CORR_1<-DGE_YUCA_N$E[,grep(c(vari), colnames(DGE_YUCA_N$E),invert =F)]
  CORR_2<-CORR_1[,grep(c("CAN"),colnames(CORR_1))]
  CORR_3<-CORR_1[,grep(c("B1"),colnames(CORR_1))]
  CORR_4<-CORR_1[,grep(c("CTRL"),colnames(CORR_1))]
  
  plot(CORR_2[,c(1,2)],main=paste("CASSAVA",vari));plot(CORR_2[,c(1,3)],main=paste("CASSAVA",vari));plot(CORR_2[,c(2,3)],main=paste("CASSAVA",vari))
  plot(CORR_3[,c(1,2)],main=paste("CASSAVA",vari));plot(CORR_3[,c(1,3)],main=paste("CASSAVA",vari));plot(CORR_3[,c(2,3)],main=paste("CASSAVA",vari))
  plot(CORR_4[,c(1,2)],main=paste("CASSAVA",vari));tryCatch(plot(CORR_4[,c(1,3)],main=paste("CASSAVA",vari)),error=function(e) 1 );tryCatch(plot(CORR_4[,c(2,3)],main=paste("CASSAVA",vari)),error=function(e) 1 )

  CORR_1<-DGE_AMF_N$E[,grep(c(vari), colnames(DGE_AMF_N$E),invert =F)]
  CORR_2<-CORR_1[,grep(c("CAN"),colnames(CORR_1))]
  CORR_3<-CORR_1[,grep(c("B1"),colnames(CORR_1))]
  plot(CORR_2[,c(1,2)],main=paste("AMF",vari));plot(CORR_2[,c(1,3)],main=paste("AMF",vari));plot(CORR_2[,c(2,3)],main=paste("AMF",vari))
  plot(CORR_3[,c(1,2)],main=paste("AMF",vari));plot(CORR_3[,c(1,3)],main=paste("AMF",vari));plot(CORR_3[,c(2,3)],main=paste("AMF",vari))
  
}
dev.off()

########################################################################################################################################
# summary transcripts TPM -> GENE
# to publications
########################################################################################################################################
DTE_G_YUCA<-c('Total','More_than_0','more_than10','more_than100')
for (x in 1:dim(DGE_YUCA_N$E)[2]){
  samp<-c(length(DGE_YUCA_N$E[,x]),dim(DGE_YUCA_N$E[DGE_YUCA_N$E[,x]>0,])[1],dim(DGE_YUCA_N$E[DGE_YUCA_N$E[,x]>10,])[1],dim(DGE_YUCA_N$E[DGE_YUCA_N$E[,x]>100,])[1])
  DTE_G_YUCA<-cbind.data.frame(DTE_G_YUCA,samp)}

colnames(DTE_G_YUCA)[2:(dim(DGE_YUCA_N$E)[2]+1)]<-colnames(DGE_YUCA_N$E)
rownames(DTE_G_YUCA)<-c('Total','More_than_0','more_than10','more_than100')
colnames(DTE_G_YUCA)
max(as.numeric(DTE_G_YUCA[3,c(3,4,6,8,9,10,5,7,
              11,12,13,14,15,16,17,18,
              19,20,21,22,23,24,25,26,2,
              27,28,29,30,31,32,33,34,35,
              37,38,39,40,41,42,43,36
              )]))
rownames(seq)
seq<-cbind.data.frame(seq,as.numeric(DTE_G_YUCA[2,c(3,4,6,8,9,10,5,7,
                                               11,12,13,14,15,16,17,18,
                                               19,20,21,22,23,24,25,26,2,
                                               27,28,29,30,31,32,33,34,35,
                                               37,38,39,40,41,42,43,36)]))
colnames(seq)[dim(seq)[2]]<-"Gene_counts"


write.table(DTE_G_YUCA,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Genes_summary_CASSAVA_sample.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)
##############################################
# recreate normalized gene counts with control for AMF
# take out samples
FOUR_VARS_AMF<-txi_AMF[[2]]
colnames(FOUR_VARS_AMF)<-gsub("UmapAMF_","",colnames(FOUR_VARS_AMF),perl=T)
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("CANB1"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V3"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V8_6"), colnames(FOUR_VARS_AMF),invert =T)]

FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_9"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_18"), colnames(FOUR_VARS_AMF),invert =T)]
FOUR_VARS_AMF<-FOUR_VARS_AMF[,grep(c("V5_6"), colnames(FOUR_VARS_AMF),invert =T)]

DGE_AMF <- DGEList(FOUR_VARS_AMF)
#FILTERING  


keep_A <- rowSums(DGE_AMF$counts>100) >= 3
DGE_AMF<- DGE_AMF[keep_A,]

#NORMALIZATION
DGE_AMF <- calcNormFactors(DGE_AMF)

REPLICA<-unlist(lapply(strsplit(colnames(FOUR_VARS_AMF),"_"),function (x) x[[2]]))
VAR<-gsub( " ","",gsub( "m","",gsub( "m2","",gsub("_\\w+$","",colnames(FOUR_VARS_AMF),perl=T)    ) ))
TREAT<-gsub("^(\\w+)_","",colnames(FOUR_VARS_AMF),perl=T)
TREAT_AMF<-factor(TREAT,levels=c("CTRL","CAN","B1"))
Treatsss <- factor(paste(VAR,TREAT_AMF,sep="."))


design <- model.matrix(~0+Treatsss)
colnames(design) <- c(levels(Treatsss))

DGE_AMF_N <- voom(DGE_AMF, design,plot=TRUE)

DGE_AMF_N$E



DTE_G_AMF<-c('Total','More_than_0','more_than10','more_than100')
all_AMF_4var<-DGE_AMF_N$E
for (x in 1:dim(all_AMF_4var)[2]){
  samp<-c(length(all_AMF_4var[,x]),dim(all_AMF_4var[all_AMF_4var[,x]>0,])[1],dim(all_AMF_4var[all_AMF_4var[,x]>10,])[1],dim(all_AMF_4var[all_AMF_4var[,x]>100,])[1])
  DTE_G_AMF<-cbind.data.frame(DTE_G_AMF,samp)}
colnames(DTE_G_AMF)[2:(dim(all_AMF_4var)[2]+1)]<-colnames(all_AMF_4var)
DTE_G_AMF
colnames(DTE_G_AMF)
DGE_AMF_N$E

rownames(DTE_G_AMF)<-c('Total','More_than_0','more_than10','more_than100')
colnames(DTE_G_AMF)
rownames(seq)
head(DTE_G_AMF)

seq_A<-seq
seq_A<-seq_A[c(-7,-8,-10,-15,-19,-21,-24,-27,-30,-33,-37,-40),]
seq_A<-cbind.data.frame(seq_A,as.numeric(DTE_G_AMF[2,c(3,4,5,6,7,8,
                                                   9,10,11,12,13,14,
                                                   15,16,17,18,19,2,
                                                   20,21,22,23,24,25,
                                                   27,28,29,30,31,26)]))
colnames(seq_A)[dim(seq)[2]]<-"Gene_counts_A"

seq$Samples

write.table(DTE_G_AMF,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.Genes_summary_CASSAVA_sample.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/0.Control1_TPM100.pdf")
par(mfrow=c(2,2))
boxplot(seq$Gene_counts~seq$Treatment,ylab="Genes > 100 counts",las=2,ylim=c(0,27000))
kruskal.test(seq$Gene_counts~seq$Treatment)
boxplot(seq_A$Gene_counts_A~seq_A$Treatment,ylab="Genes > 100 counts",las=2)
kruskal.test(seq_A$Gene_counts_A~seq_A$Treatment)
boxplot(seq$Gene_counts~seq$Cultivar,ylab="Genes > 100 counts",las=2,ylim=c(0,27000))
kruskal.test(seq$Gene_counts~seq$Cultivar)
boxplot(seq_A$Gene_counts_A~seq_A$Cultivar,ylab="Genes > 100 counts",las=2)
kruskal.test(seq_A$Gene_counts_A~seq_A$Cultivar)
dev.off()


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/0.Control4_TPM100_samples.pdf")
colos<-as.vector(seq$Treatment)
colos[colos=="CAN"]<-"lightsalmon3"
colos[colos=="B1"]<-"dodgerblue3"
colos[colos=="CTRL"]<-"azure3"
par(mfrow=c(2,1))
barplot(seq$Gene_counts, names=seq$Samples,las=2,col=colos,
        ylab="genes > 100 counts")
barplot(seq$Gene_counts_A, names=seq$Samples,las=2,col=colos,
        ylab="genes > 100 counts")
dev.off()




#ACCUMULATION CURVES


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/0.Control8_seq_effort.pdf")
par(mfrow=c(2,1))
options(scipen = 10)
plot(seq$Gene_counts[seq$Treatment=="CTRL"]~seq$pseudoalign_M[seq$Treatment=="CTRL"], ylim=c(0,30000), pch=17,col="azure3",ylab=" Genes > 100 counts",xlab="Reads",las=1)
points(seq$Gene_counts[seq$Treatment=="CAN"]~seq$pseudoalign_M[seq$Treatment=="CAN"], ylim=c(0,25000), pch=16,col="lightsalmon3")
points(seq$Gene_counts[seq$Treatment=="B1"]~seq$pseudoalign_M[seq$Treatment=="B1"], ylim=c(0,25000), pch=15,col="dodgerblue3")

plot(seq$Gene_counts_A[seq$Treatment=="CTRL"] ~seq$pseudoaligned_R[seq$Treatment=="CTRL"], ylim=c(0,25000),xlim=c(0,8000000),  pch=17,col="azure3",ylab=" Genes > 100 counts",xlab="Reads",las=1)
points(seq$Gene_counts_A[seq$Treatment=="CAN"] ~seq$pseudoaligned_R[seq$Treatment=="CAN"],  pch=16,col="lightsalmon3")
points(seq$Gene_counts_A[seq$Treatment=="B1"] ~seq$pseudoaligned_R[seq$Treatment=="B1"], ylim=c(0,2000),  pch=15,col="dodgerblue3")
dev.off()





########################################################################################################################################
########################################################################################################################################

# DGE

########################################################################################################################################

########################################################################################################################################

# DGE
## Differential gene expression
# VERSION FINAL


REPLICA<-unlist(lapply(strsplit(colnames(FOUR_VARS_YUCA),"_"),function (x) x[[2]]))
VAR<-gsub( " ","",gsub( "m","",gsub( "m2","",gsub("_\\w+$","",colnames(FOUR_VARS_YUCA),perl=T)    ) ))
TREAT<-gsub("^(\\w+)_","",colnames(FOUR_VARS_YUCA),perl=T)
TREAT_YUCA<-factor(TREAT,levels=c("CTRL","CAN","B1"))

Treatsss <- factor(paste(VAR,TREAT_YUCA,sep="."))


design <- model.matrix(~0+Treatsss)
colnames(design) <- c(levels(Treatsss))

DGE_YUCA_N <- voom(DGE_YUCA, design,plot=TRUE)


#corfit$consensu
fit_YUCA <- lmFit(DGE_YUCA_N,design)
cm <- makeContrasts(V1_B1.CTRL= V1.B1-V1.CTRL,
                    V1_CAN.CTRL = V1.CAN-V1.CTRL,
                    V1_CAN.B1 = V1.CAN-V1.B1,
                    V4_B1.CTRL= V4.B1-V4.CTRL,
                    V4_CAN.CTRL = V4.CAN-V4.CTRL,
                    V4_CAN.B1 = V4.CAN-V4.B1,
                    V5_B1.CTRL= V5.B1-V5.CTRL,
                    V5_CAN.CTRL = V5.CAN-V5.CTRL,
                    V5_CAN.B1 = V5.CAN-V5.B1,
                    V6_B1.CTRL= V6.B1-V6.CTRL,
                    V6_CAN.CTRL = V6.CAN-V6.CTRL,
                    V6_CAN.B1 = V6.CAN-V6.B1,
                    V8_B1.CTRL= V8.B1-V8.CTRL,
                    V8_CAN.CTRL = V8.CAN-V8.CTRL,
                    V8_CAN.B1 = V8.CAN-V8.B1,
                    V1_V4_CAN= V1.CAN-V4.CAN,
                    V1_V5_CAN= V1.CAN-V5.CAN,
                    V1_V6_CAN= V1.CAN-V6.CAN,
                    V1_V8_CAN= V1.CAN-V8.CAN,
                    V4_V5_CAN= V4.CAN-V5.CAN,
                    V4_V6_CAN= V4.CAN-V6.CAN,
                    V4_V8_CAN= V4.CAN-V8.CAN,
                    V5_V6_CAN= V5.CAN-V6.CAN,
                    V5_V8_CAN= V5.CAN-V8.CAN,
                    V6_V8_CAN= V6.CAN-V8.CAN,
                    V1_V4_B1= V1.B1-V4.B1,
                    V1_V5_B1= V1.B1-V5.B1,
                    V1_V6_B1= V1.B1-V6.B1,
                    V1_V8_B1= V1.B1-V8.B1,
                    V4_V5_B1= V4.B1-V5.B1,
                    V4_V6_B1= V4.B1-V6.B1,
                    V4_V8_B1= V4.B1-V8.B1,
                    V5_V6_B1= V5.B1-V6.B1,
                    V5_V8_B1= V5.B1-V8.B1,
                    V6_V8_B1= V6.B1-V8.B1,
                    V1_V4_CTRL= V1.CTRL-V4.CTRL,
                    V1_V5_CTRL= V1.CTRL-V5.CTRL,
                    V1_V6_CTRL= V1.CTRL-V6.CTRL,
                    V1_V8_CTRL= V1.CTRL-V8.CTRL,
                    V4_V5_CTRL= V4.CTRL-V5.CTRL,
                    V4_V6_CTRL= V4.CTRL-V6.CTRL,
                    V4_V8_CTRL= V4.CTRL-V8.CTRL,
                    V5_V6_CTRL= V5.CTRL-V6.CTRL,
                    V5_V8_CTRL= V5.CTRL-V8.CTRL,
                    V6_V8_CTRL= V6.CTRL-V8.CTRL,
                    levels=design)

fit_YUCA<- eBayes(contrasts.fit(fit_YUCA, cm))

dim(fit_YUCA$t)
dim(fit_AMF$t)
###
#AMF treatment per variety
sign_YUCA_V1<-topTable(fit_YUCA, coef=c(1:3),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V4<-topTable(fit_YUCA, coef=c(4:6),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V5<-topTable(fit_YUCA, coef=c(7:9),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V6<-topTable(fit_YUCA, coef=c(10:12),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V8<-topTable(fit_YUCA, coef=c(13:15),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL

sign_YUCA_V1<-sign_YUCA_V1[sign_YUCA_V1$adj.P.Val<0.05,]
sign_YUCA_V4<-sign_YUCA_V4[sign_YUCA_V4$adj.P.Val<0.05,]
sign_YUCA_V5<-sign_YUCA_V5[sign_YUCA_V5$adj.P.Val<0.05,]
sign_YUCA_V6<-sign_YUCA_V6[sign_YUCA_V6$adj.P.Val<0.05,]
sign_YUCA_V8<-sign_YUCA_V8[sign_YUCA_V8$adj.P.Val<0.05,]
##
#  Different CAN-B1
dim(sign_YUCA_V1[sign_YUCA_V1[,3]>1| sign_YUCA_V1[,3]< (-1),])
dim(sign_YUCA_V4[sign_YUCA_V4[,3]>1| sign_YUCA_V4[,3]< (-1),])
dim(sign_YUCA_V5[sign_YUCA_V5[,3]>1| sign_YUCA_V5[,3]< (-1),])
dim(sign_YUCA_V6[sign_YUCA_V6[,3]>1| sign_YUCA_V6[,3]< (-1),])
dim(sign_YUCA_V8[sign_YUCA_V8[,3]>1| sign_YUCA_V8[,3]< (-1),])

########
# mock-b1 DW diff vs CAN mock in col 2215.
M_B1_2215<-topTable(fit_YUCA, coef=c(4),number=dim(DGE_YUCA_N)[1])
M_B1_2215[M_B1_2215$adj.P.Val<0.05,]
M_CAN_2215<-topTable(fit_YUCA, coef=c(5),number=dim(DGE_YUCA_N)[1])
M_CAN_2215[M_CAN_2215$adj.P.Val<0.05,]

# genes significant in M_B1 and not in M_CAN
M_B1_2215[M_B1_2215$adj.P.Val<0.05,][rownames(M_B1_2215[M_B1_2215$adj.P.Val<0.05,]) %in% rownames(M_CAN_2215[M_CAN_2215$adj.P.Val>0.05,]),]
M_B1_2215[M_B1_2215$adj.P.Val<0.05,][rownames(M_B1_2215[M_B1_2215$adj.P.Val<0.05,]) %in% rownames(M_CAN_2215[M_CAN_2215$P.Valuel>0.05,]),]


M_CAN_2215[grep("Manes.02G030800",rownames(M_CAN_2215)),]

# genes significant in M_CAN and not in M_B1
M_CAN_2215[M_CAN_2215$adj.P.Val<0.05,][rownames(M_CAN_2215[M_CAN_2215$adj.P.Val<0.05,]) %in% rownames(M_B1_2215[M_B1_2215$adj.P.Val>0.05,]),]
M_CAN_2215[M_CAN_2215$adj.P.Val<0.05,][rownames(M_CAN_2215[M_CAN_2215$adj.P.Val<0.05,]) %in% rownames(M_B1_2215[M_B1_2215$P.Value>0.05,]),]


M_B1_2215[grep("Manes.04G158800",rownames(M_B1_2215)),]

#install.packages('venneuler')
library(venneuler)

#########################################################################################
#YUCA by AMF diff 
YUCA_byAMF_diff<-unique(c(rownames(sign_YUCA_V1),
                   rownames(sign_YUCA_V4),
                   rownames(sign_YUCA_V5),
                   rownames(sign_YUCA_V6),
                   rownames(sign_YUCA_V8)))
sign_YUCA_V1<-topTable(fit_YUCA, coef=c(1:3),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V4<-topTable(fit_YUCA, coef=c(4:6),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V5<-topTable(fit_YUCA, coef=c(7:9),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V6<-topTable(fit_YUCA, coef=c(10:12),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V8<-topTable(fit_YUCA, coef=c(13:15),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL

sign_YUCA_V1<-cbind.data.frame(rownames(sign_YUCA_V1[rownames(sign_YUCA_V1) %in% YUCA_byAMF_diff ,]),sign_YUCA_V1[rownames(sign_YUCA_V1) %in% YUCA_byAMF_diff ,])
sign_YUCA_V4<-cbind.data.frame(rownames(sign_YUCA_V4[rownames(sign_YUCA_V4) %in% YUCA_byAMF_diff ,]),sign_YUCA_V4[rownames(sign_YUCA_V4) %in% YUCA_byAMF_diff ,])
sign_YUCA_V5<-cbind.data.frame(rownames(sign_YUCA_V5[rownames(sign_YUCA_V5) %in% YUCA_byAMF_diff ,]),sign_YUCA_V5[rownames(sign_YUCA_V5) %in% YUCA_byAMF_diff ,])
sign_YUCA_V6<-cbind.data.frame(rownames(sign_YUCA_V6[rownames(sign_YUCA_V6) %in% YUCA_byAMF_diff ,]),sign_YUCA_V6[rownames(sign_YUCA_V6) %in% YUCA_byAMF_diff ,])
sign_YUCA_V8<-cbind.data.frame(rownames(sign_YUCA_V8[rownames(sign_YUCA_V8) %in% YUCA_byAMF_diff ,]),sign_YUCA_V8[rownames(sign_YUCA_V8) %in% YUCA_byAMF_diff ,])
colnames(sign_YUCA_V1)[c(1,8)]<-c("gene","V1_adj.P.Val");colnames(sign_YUCA_V4)[c(1,8)]<-c("gene","V4_adj.P.Val");colnames(sign_YUCA_V5)[c(1,8)]<-c("gene","V5_adj.P.Val")
colnames(sign_YUCA_V6)[c(1,8)]<-c("gene","V6_adj.P.Val");colnames(sign_YUCA_V8)[c(1,8)]<-c("gene","V8_adj.P.Val")


V_YUCACANB1_pheatmap<-merge(sign_YUCA_V1[c(1,2,3,4,8)],merge(sign_YUCA_V4[c(1,2,3,4,8)],merge(sign_YUCA_V5[c(1,2,3,4,8)],merge(sign_YUCA_V6[c(1,2,3,4,8)],sign_YUCA_V8[c(1,2,3,4,8)],by="gene"),by="gene"),by="gene"),by="gene")

write.table(V_YUCACANB1_pheatmap,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/Supp.table3_YUCA_AMF1634.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)


#########################################################################################
# CAN-B1 _diff
CANB1_diff<-unique(c(rownames(sign_YUCA_V1[sign_YUCA_V1[,4]>2| sign_YUCA_V1[,4]< (-2),]),
                     rownames(sign_YUCA_V4[sign_YUCA_V4[,4]>2| sign_YUCA_V4[,4]< (-2),]),
                     rownames(sign_YUCA_V5[sign_YUCA_V5[,4]>2| sign_YUCA_V5[,4]< (-2),]),
                     rownames(sign_YUCA_V6[sign_YUCA_V6[,4]>2| sign_YUCA_V6[,4]< (-2),]),
                     rownames(sign_YUCA_V8[sign_YUCA_V8[,4]>2| sign_YUCA_V8[,4]< (-2),])))

unique(c(rownames(sign_YUCA_V1[sign_YUCA_V1[,4]>2 & sign_YUCA_V1$P.Value<0.05| sign_YUCA_V1[,4]< (-2) & sign_YUCA_V1$P.Value<0.05,]),
         rownames(sign_YUCA_V4[sign_YUCA_V4[,4]>2 & sign_YUCA_V4$P.Value<0.05| sign_YUCA_V4[,4]< (-2) & sign_YUCA_V4$P.Value<0.05,]),
         rownames(sign_YUCA_V5[sign_YUCA_V5[,4]>2 & sign_YUCA_V5$P.Value<0.05| sign_YUCA_V5[,4]< (-2) & sign_YUCA_V5$P.Value<0.05,]),
         rownames(sign_YUCA_V6[sign_YUCA_V6[,4]>2 & sign_YUCA_V6$P.Value<0.05| sign_YUCA_V6[,4]< (-2) & sign_YUCA_V6$P.Value<0.05, ]),
         rownames(sign_YUCA_V8[sign_YUCA_V8[,4]>2 & sign_YUCA_V8$P.Value<0.05| sign_YUCA_V8[,4]< (-2) & sign_YUCA_V8$P.Value<0.05,])))

unique(c(rownames(sign_YUCA_V1[sign_YUCA_V1[,4]>2 & sign_YUCA_V1$P.Value<0.05| sign_YUCA_V1[,4]< (-2) & sign_YUCA_V1$P.Value<0.05,]),
         rownames(sign_YUCA_V4[sign_YUCA_V4[,4]>2 & sign_YUCA_V4$P.Value<0.05| sign_YUCA_V4[,4]< (-2) & sign_YUCA_V4$P.Value<0.05,]),
         rownames(sign_YUCA_V5[sign_YUCA_V5[,4]>2 & sign_YUCA_V5$P.Value<0.05| sign_YUCA_V5[,4]< (-2) & sign_YUCA_V5$P.Value<0.05,]),
         rownames(sign_YUCA_V6[sign_YUCA_V6[,4]>2 & sign_YUCA_V6$P.Value<0.05| sign_YUCA_V6[,4]< (-2) & sign_YUCA_V6$P.Value<0.05, ]),
         rownames(sign_YUCA_V8[sign_YUCA_V8[,4]>2 & sign_YUCA_V8$P.Value<0.05| sign_YUCA_V8[,4]< (-2) & sign_YUCA_V8$P.Value<0.05,])))

sign_YUCA_all1<-topTable(fit_YUCA, coef=c(1:3),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_all4<-topTable(fit_YUCA, coef=c(4:6),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_all5<-topTable(fit_YUCA, coef=c(7:9),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_all6<-topTable(fit_YUCA, coef=c(10:12),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_all8<-topTable(fit_YUCA, coef=c(13:15),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL


dim(sign_YUCA_all4[rownames(sign_YUCA_all4) %in% CANB1_diff ,])

sign_YUCA_all1<-cbind.data.frame(rownames(sign_YUCA_all1[rownames(sign_YUCA_all1) %in% CANB1_diff ,]),sign_YUCA_all1[rownames(sign_YUCA_all1) %in% CANB1_diff ,])
sign_YUCA_all4<-cbind.data.frame(rownames(sign_YUCA_all4[rownames(sign_YUCA_all4) %in% CANB1_diff ,]),sign_YUCA_all4[rownames(sign_YUCA_all4) %in% CANB1_diff ,])
sign_YUCA_all5<-cbind.data.frame(rownames(sign_YUCA_all5[rownames(sign_YUCA_all5) %in% CANB1_diff ,]),sign_YUCA_all5[rownames(sign_YUCA_all5) %in% CANB1_diff ,])
sign_YUCA_all6<-cbind.data.frame(rownames(sign_YUCA_all6[rownames(sign_YUCA_all6) %in% CANB1_diff ,]),sign_YUCA_all6[rownames(sign_YUCA_all6) %in% CANB1_diff ,])
sign_YUCA_all8<-cbind.data.frame(rownames(sign_YUCA_all8[rownames(sign_YUCA_all8) %in% CANB1_diff ,]),sign_YUCA_all8[rownames(sign_YUCA_all8) %in% CANB1_diff ,])
colnames(sign_YUCA_all1)[c(1,8)]<-c("gene","V1_adj.P.Val");colnames(sign_YUCA_all4)[c(1,8)]<-c("gene","V4_adj.P.Val");colnames(sign_YUCA_all5)[c(1,8)]<-c("gene","V5_adj.P.Val")
colnames(sign_YUCA_all6)[c(1,8)]<-c("gene","V6_adj.P.Val");colnames(sign_YUCA_all8)[c(1,8)]<-c("gene","V8_adj.P.Val")


all_YUCACANB1_pheatmap<-merge(sign_YUCA_all1[c(1,4,8)],merge(sign_YUCA_all4[c(1,4,8)],merge(sign_YUCA_all5[c(1,4,8)],merge(sign_YUCA_all6[c(1,4,8)],sign_YUCA_all8[c(1,4,8)],by="gene"),by="gene"),by="gene"),by="gene")

write.table(all_YUCACANB1_pheatmap,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/Supp.table3_YUCA_CANB1FC2.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/figure1b_CANB1_cassava.pdf")
pheatmap(t(merge(sign_YUCA_all1[c(1,4)],merge(sign_YUCA_all4[c(1,4)],merge(sign_YUCA_all5[c(1,4)],merge(sign_YUCA_all6[c(1,4)],sign_YUCA_all8[c(1,4)],by="gene"),by="gene"),by="gene"),by="gene")
[,-1]),cellheight = 12,cellwidth =2 ,scale = 'column')
dev.off()
#########################################################################################

#VAR treatment per AMF
sign_YUCA_CAN<-topTable(fit_YUCA, coef=c(16:25),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_B1<-topTable(fit_YUCA, coef=c(26:35),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_CTRL<-topTable(fit_YUCA, coef=c(36:45),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL

sign_YUCA_CAN<-sign_YUCA_CAN[sign_YUCA_CAN$adj.P.Val<0.05,]
sign_YUCA_B1<-sign_YUCA_B1[sign_YUCA_B1$adj.P.Val<0.05,]
sign_YUCA_CTRL<-sign_YUCA_CTRL[sign_YUCA_CTRL$adj.P.Val<0.05,]

##
# difference VAR

VAR_diffe<-unique(c(rownames(sign_YUCA_CAN),rownames(sign_YUCA_B1),rownames(sign_YUCA_CTRL)))

#########################################################################################
# only VAR
VAR_diffe[!VAR_diffe %in% YUCA_byAMF_diff]

sign_YUCA_CAN<-topTable(fit_YUCA, coef=c(16:25),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_B1<-topTable(fit_YUCA, coef=c(26:35),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_CTRL<-topTable(fit_YUCA, coef=c(36:45),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL

sign_YUCA_CAN<-cbind.data.frame(rownames(sign_YUCA_CAN[rownames(sign_YUCA_CAN) %in% VAR_diffe[!VAR_diffe %in% YUCA_byAMF_diff] ,]),sign_YUCA_CAN[rownames(sign_YUCA_CAN) %in% VAR_diffe[!VAR_diffe %in% YUCA_byAMF_diff] ,])
sign_YUCA_B1<- cbind.data.frame(rownames(sign_YUCA_B1[rownames(sign_YUCA_B1) %in% VAR_diffe[!VAR_diffe %in% YUCA_byAMF_diff] ,]),sign_YUCA_B1[rownames(sign_YUCA_B1) %in% VAR_diffe[!VAR_diffe %in% YUCA_byAMF_diff] ,])
sign_YUCA_CTRL<-cbind.data.frame(rownames(sign_YUCA_CTRL[rownames(sign_YUCA_CTRL) %in% VAR_diffe[!VAR_diffe %in% YUCA_byAMF_diff] ,]),sign_YUCA_CTRL[rownames(sign_YUCA_CTRL) %in% VAR_diffe[!VAR_diffe %in% YUCA_byAMF_diff] ,])
colnames(sign_YUCA_CAN)[c(1,15)]<-c("gene","CAN_adj.P.Val");colnames(sign_YUCA_B1)[c(1,15)]<-c("gene","B1_adj.P.Val")
colnames(sign_YUCA_CTRL)[c(1,15)]<-c("gene","CTRL_adj.P.Val")
YUCAbyVAR<-merge(sign_YUCA_CAN[c(1:11,15)],merge(sign_YUCA_B1[c(1:11,15)],sign_YUCA_CTRL[c(1:11,15)],by="gene"),by="gene")

write.table(YUCAbyVAR,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/Supp.table2_YUCA_VAR_8559.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)





#########################################################################################
# only AMF and not VAR
YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe]

sign_YUCA_V1<-topTable(fit_YUCA, coef=c(1:3),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V4<-topTable(fit_YUCA, coef=c(4:6),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V5<-topTable(fit_YUCA, coef=c(7:9),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V6<-topTable(fit_YUCA, coef=c(10:12),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V8<-topTable(fit_YUCA, coef=c(13:15),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL

sign_YUCA_V1<-cbind.data.frame(rownames(sign_YUCA_V1[rownames(sign_YUCA_V1) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V1[rownames(sign_YUCA_V1) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,])
sign_YUCA_V4<-cbind.data.frame(rownames(sign_YUCA_V4[rownames(sign_YUCA_V4) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V4[rownames(sign_YUCA_V4) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,])
sign_YUCA_V5<-cbind.data.frame(rownames(sign_YUCA_V5[rownames(sign_YUCA_V5) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V5[rownames(sign_YUCA_V5) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,])
sign_YUCA_V6<-cbind.data.frame(rownames(sign_YUCA_V6[rownames(sign_YUCA_V6) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V6[rownames(sign_YUCA_V6) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,])
sign_YUCA_V8<-cbind.data.frame(rownames(sign_YUCA_V8[rownames(sign_YUCA_V8) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V8[rownames(sign_YUCA_V8) %in% YUCA_byAMF_diff[!YUCA_byAMF_diff %in% VAR_diffe] ,])
colnames(sign_YUCA_V1)[c(1,8)]<-c("gene","V1_adj.P.Val");colnames(sign_YUCA_V4)[c(1,8)]<-c("gene","V4_adj.P.Val");colnames(sign_YUCA_V5)[c(1,8)]<-c("gene","V5_adj.P.Val")
colnames(sign_YUCA_V6)[c(1,8)]<-c("gene","V6_adj.P.Val");colnames(sign_YUCA_V8)[c(1,8)]<-c("gene","V8_adj.P.Val")


YUCAbyAMF<-merge(sign_YUCA_V1[c(1:4,8)],merge(sign_YUCA_V4[c(1:4,8)],merge(sign_YUCA_V5[c(1:4,8)],merge(sign_YUCA_V6[c(1:4,8)],sign_YUCA_V8[c(1:4,8)],by="gene"),by="gene"),by="gene"),by="gene")

write.table(YUCAbyAMF,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/Supp.table2_YUCA_AMF_450.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)

#########################################################################################
# cassava genes influenced by VAR and AMF
YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe]

sign_YUCA_V1<-topTable(fit_YUCA, coef=c(1:3),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V4<-topTable(fit_YUCA, coef=c(4:6),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V5<-topTable(fit_YUCA, coef=c(7:9),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V6<-topTable(fit_YUCA, coef=c(10:12),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V8<-topTable(fit_YUCA, coef=c(13:15),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL

sign_YUCA_V1<-cbind.data.frame(rownames(sign_YUCA_V1[rownames(sign_YUCA_V1) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V1[rownames(sign_YUCA_V1) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,])
sign_YUCA_V4<-cbind.data.frame(rownames(sign_YUCA_V4[rownames(sign_YUCA_V4) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V4[rownames(sign_YUCA_V4) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,])
sign_YUCA_V5<-cbind.data.frame(rownames(sign_YUCA_V5[rownames(sign_YUCA_V5) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V5[rownames(sign_YUCA_V5) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,])
sign_YUCA_V6<-cbind.data.frame(rownames(sign_YUCA_V6[rownames(sign_YUCA_V6) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V6[rownames(sign_YUCA_V6) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,])
sign_YUCA_V8<-cbind.data.frame(rownames(sign_YUCA_V8[rownames(sign_YUCA_V8) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,]),sign_YUCA_V8[rownames(sign_YUCA_V8) %in% YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe] ,])
colnames(sign_YUCA_V1)[c(1,8)]<-c("gene","V1_adj.P.Val");colnames(sign_YUCA_V4)[c(1,8)]<-c("gene","V4_adj.P.Val");colnames(sign_YUCA_V5)[c(1,8)]<-c("gene","V5_adj.P.Val")
colnames(sign_YUCA_V6)[c(1,8)]<-c("gene","V6_adj.P.Val");colnames(sign_YUCA_V8)[c(1,8)]<-c("gene","V8_adj.P.Val")


YUCAbyAMF_VAR<-merge(sign_YUCA_V1[c(1:4,8)],merge(sign_YUCA_V4[c(1:4,8)],merge(sign_YUCA_V5[c(1:4,8)],merge(sign_YUCA_V6[c(1:4,8)],sign_YUCA_V8[c(1:4,8)],by="gene"),by="gene"),by="gene"),by="gene")

write.table(YUCAbyAMF_VAR,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/Supp.table2_YUCA_AMF&VAR.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)
#########################################################################################




VAR_diffe[!VAR_diffe %in% YUCA_byAMF_diff]

CANB1_diff[CANB1_diff %in% VAR_diffe]


YUCA_byAMF_diff[YUCA_byAMF_diff %in% VAR_diffe]

#########################################################################################

# FRED SYMBIOSIS GENES
Fred_YUCA_genes<-  c("Manes.06G156700","Manes.12G124500","Manes.12G124500","Manes.01G123300",
                     "Manes.04G002400","Manes.12G124300","Manes.05G053600","Manes.08G154600",
                     "Manes.18G034100","Manes.16G125700","Manes.06G143100","Manes.01G071000","Manes.01G193000")
sign_YUCA<-topTable(fit_YUCA, coef=c(1:15),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL

sign_YUCA[rownames(sign_YUCA) %in% Fred_YUCA_genes,]
pheatmap(sign_YUCA[rownames(sign_YUCA) %in% Fred_YUCA_genes,1:10])
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Fig2_symb_genes.pdf',height=10,width = 10)

symb<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Fred_YUCA_genes,]
par(las=2)
rownames(symb)<-gsub("Manes","M",rownames(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Fred_YUCA_genes,]))
pheatmap(symb,cellwidth =8 ,cellheight=8,)
dev.off()


#########################################################################################
#########################################################################################
# important genes that differ between CAN B1 and VAR.

###CASSAVA
sign_YUCA_V1<-topTable(fit_YUCA, coef=c(1:3),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V4<-topTable(fit_YUCA, coef=c(4:6),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V5<-topTable(fit_YUCA, coef=c(7:9),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V6<-topTable(fit_YUCA, coef=c(10:12),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
sign_YUCA_V8<-topTable(fit_YUCA, coef=c(13:15),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL

# cassava genes FC > 8x between CAN B1.
sign_YUCA_V1[sign_YUCA_V1$V1_CAN.B1>3 | sign_YUCA_V1$V1_CAN.B1<(-3),]
sign_YUCA_V4[sign_YUCA_V4$V4_CAN.B1>3 | sign_YUCA_V4$V4_CAN.B1<(-3),]
sign_YUCA_V5[sign_YUCA_V5$V5_CAN.B1>3 | sign_YUCA_V5$V5_CAN.B1<(-3),]
sign_YUCA_V6[sign_YUCA_V6$V6_CAN.B1>3 | sign_YUCA_V6$V6_CAN.B1<(-3),]
sign_YUCA_V8[sign_YUCA_V8$V8_CAN.B1>3 | sign_YUCA_V8$V8_CAN.B1<(-3),]

Y_FC8<-unique(c(rownames(sign_YUCA_V1[sign_YUCA_V1$V1_CAN.B1>3 | sign_YUCA_V1$V1_CAN.B1<(-3),]),
rownames(sign_YUCA_V4[sign_YUCA_V4$V4_CAN.B1>3 | sign_YUCA_V4$V4_CAN.B1<(-3),]),
rownames(sign_YUCA_V5[sign_YUCA_V5$V5_CAN.B1>3 | sign_YUCA_V5$V5_CAN.B1<(-3),]),
rownames(sign_YUCA_V6[sign_YUCA_V6$V6_CAN.B1>3 | sign_YUCA_V6$V6_CAN.B1<(-3),]),
rownames(sign_YUCA_V8[sign_YUCA_V8$V8_CAN.B1>3 | sign_YUCA_V8$V8_CAN.B1<(-3),]) ))


Y_FC8_V1<-sign_YUCA_V1[rownames(sign_YUCA_V1) %in% Y_FC8,]
Y_FC8_V4<-sign_YUCA_V4[rownames(sign_YUCA_V4) %in% Y_FC8,]
Y_FC8_V5<-sign_YUCA_V5[rownames(sign_YUCA_V5) %in% Y_FC8,]
Y_FC8_V6<-sign_YUCA_V6[rownames(sign_YUCA_V6) %in% Y_FC8,]
Y_FC8_V8<-sign_YUCA_V8[rownames(sign_YUCA_V8) %in% Y_FC8,]
Y_FC8_V1$gene<-rownames(Y_FC8_V1);Y_FC8_V4$gene<-rownames(Y_FC8_V4);Y_FC8_V5$gene<-rownames(Y_FC8_V5)
Y_FC8_V6$gene<-rownames(Y_FC8_V6);Y_FC8_V8$gene<-rownames(Y_FC8_V8)
Y_FC8_ALL<-merge(Y_FC8_V1,merge(Y_FC8_V4,merge(Y_FC8_V5,merge(Y_FC8_V6,Y_FC8_V8,by="gene"),by="gene"),by="gene"),by="gene")

Y_FC8_ALL_CANB1<-Y_FC8_ALL[,grep("CAN.B1",colnames(Y_FC8_ALL),invert=F)]
pheatmap(Y_FC8_ALL_CANB1)


Y_FC8_INFO<-merge(Y_FC8_ALL,Mercator_Mesculenta,by="gene")
Y_FC8_INFO<-Y_FC8_INFO[!duplicated(Y_FC8_INFO$gene),]
Y_FC8_INFO[,40]<- unlist(lapply(lapply(strsplit(as.character(Y_FC8_INFO[,38]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
Y_FC8_INFO[,40]<-as.factor(Y_FC8_INFO[,40])

table(Y_FC8_INFO[,40])

Y_FC8_function<-Y_FC8_INFO[grep("stress.biotic",Y_FC8_INFO[,40]),]
rownames(Y_FC8_function)<-Y_FC8_function$gene
pheatmap(Y_FC8_function[,grep("CAN.B1",colnames(Y_FC8_function),invert=F)])

###
#RAM2
sign_YUCA_V1[rownames(sign_YUCA_V1) %in% "Manes.01G193000",]
sign_YUCA_V4[rownames(sign_YUCA_V4) %in% "Manes.01G193000",]
sign_YUCA_V5[rownames(sign_YUCA_V5) %in% "Manes.01G193000",]
sign_YUCA_V6[rownames(sign_YUCA_V6) %in% "Manes.01G193000",]
sign_YUCA_V8[rownames(sign_YUCA_V8) %in% "Manes.01G193000",]

#########################################################################################
#########################################################################################
#########################################################################################

### AMF
# VERSION FINAL


REPLICA<-unlist(lapply(strsplit(colnames(FOUR_VARS_AMF),"_"),function (x) x[[2]]))
VAR<-gsub( " ","",gsub( "m","",gsub( "m2","",gsub("_\\w+$","",colnames(FOUR_VARS_AMF),perl=T)    ) ))
TREAT<-gsub("^(\\w+)_","",colnames(FOUR_VARS_AMF),perl=T)
TREAT_AMF<-factor(TREAT,levels=c("CAN","B1"))

Treatsss <- factor(paste(VAR,TREAT_AMF,sep="."))


design <- model.matrix(~0+Treatsss)
colnames(design) <- c(levels(Treatsss))

DGE_AMF_N <- voom(DGE_AMF, design,plot=TRUE)


#corfit$consensu
fit_AMF <- lmFit(DGE_AMF_N,design)
cm <- makeContrasts(V1_CAN.B1 = V1.CAN-V1.B1,
                    V4_CAN.B1 = V4.CAN-V4.B1,
                    V5_CAN.B1 = V5.CAN-V5.B1,
                    V6_CAN.B1 = V6.CAN-V6.B1,
                    V8_CAN.B1 = V8.CAN-V8.B1,
                    V1_V4_CAN= V1.CAN-V4.CAN,
                    V1_V5_CAN= V1.CAN-V5.CAN,
                    V1_V6_CAN= V1.CAN-V6.CAN,
                    V1_V8_CAN= V1.CAN-V8.CAN,
                    V4_V5_CAN= V4.CAN-V5.CAN,
                    V4_V6_CAN= V4.CAN-V6.CAN,
                    V4_V8_CAN= V4.CAN-V8.CAN,
                    V5_V6_CAN= V5.CAN-V6.CAN,
                    V5_V8_CAN= V5.CAN-V8.CAN,
                    V6_V8_CAN= V6.CAN-V8.CAN,
                    V1_V4_B1= V1.B1-V4.B1,
                    V1_V5_B1= V1.B1-V5.B1,
                    V1_V6_B1= V1.B1-V6.B1,
                    V1_V8_B1= V1.B1-V8.B1,
                    V4_V5_B1= V4.B1-V5.B1,
                    V4_V6_B1= V4.B1-V6.B1,
                    V4_V8_B1= V4.B1-V8.B1,
                    V5_V6_B1= V5.B1-V6.B1,
                    V5_V8_B1= V5.B1-V8.B1,
                    V6_V8_B1= V6.B1-V8.B1,
                    levels=design)

fit_AMF<- eBayes(contrasts.fit(fit_AMF, cm))

#AMF treatment per variety
sign_AMF_V1<-topTable(fit_AMF, coef=c(1),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V4<-topTable(fit_AMF, coef=c(2),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V5<-topTable(fit_AMF, coef=c(3),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V6<-topTable(fit_AMF, coef=c(4),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V8<-topTable(fit_AMF, coef=c(5),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL

sign_AMF_V1<-sign_AMF_V1[sign_AMF_V1$adj.P.Val<0.05,]
sign_AMF_V4<-sign_AMF_V4[sign_AMF_V4$adj.P.Val<0.05,]
sign_AMF_V5<-sign_AMF_V5[sign_AMF_V5$adj.P.Val<0.05,]
sign_AMF_V6<-sign_AMF_V6[sign_AMF_V6$adj.P.Val<0.05,]
sign_AMF_V8<-sign_AMF_V8[sign_AMF_V8$adj.P.Val<0.05,]


###################################################
# CAN-B1 FC 
CANB1_diff_A<-unique(c(rownames(sign_AMF_V1[sign_AMF_V1[,1]>4| sign_AMF_V1[,1]< (-4),]),
         rownames(sign_AMF_V4[sign_AMF_V4[,1]>4| sign_AMF_V4[,1]< (-4),]),
         rownames(sign_AMF_V5[sign_AMF_V5[,1]>4| sign_AMF_V5[,1]< (-4),]),
         rownames(sign_AMF_V6[sign_AMF_V6[,1]>4| sign_AMF_V6[,1]< (-4),]),
         rownames(sign_AMF_V8[sign_AMF_V8[,1]>4| sign_AMF_V8[,1]< (-4),])))

sign_AMF_V1<-topTable(fit_AMF, coef=c(1),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V4<-topTable(fit_AMF, coef=c(2),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V5<-topTable(fit_AMF, coef=c(3),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V6<-topTable(fit_AMF, coef=c(4),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V8<-topTable(fit_AMF, coef=c(5),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL

sign_AMF_V1<-cbind.data.frame(rownames(sign_AMF_V1[rownames(sign_AMF_V1) %in% CANB1_diff_A ,]),sign_AMF_V1[rownames(sign_AMF_V1) %in% CANB1_diff_A ,])
sign_AMF_V4<-cbind.data.frame(rownames(sign_AMF_V4[rownames(sign_AMF_V4) %in% CANB1_diff_A ,]),sign_AMF_V4[rownames(sign_AMF_V4) %in% CANB1_diff_A ,])
sign_AMF_V5<-cbind.data.frame(rownames(sign_AMF_V5[rownames(sign_AMF_V5) %in% CANB1_diff_A ,]),sign_AMF_V5[rownames(sign_AMF_V5) %in% CANB1_diff_A ,])
sign_AMF_V6<-cbind.data.frame(rownames(sign_AMF_V6[rownames(sign_AMF_V6) %in% CANB1_diff_A ,]),sign_AMF_V6[rownames(sign_AMF_V6) %in% CANB1_diff_A ,])
sign_AMF_V8<-cbind.data.frame(rownames(sign_AMF_V8[rownames(sign_AMF_V8) %in% CANB1_diff_A ,]),sign_AMF_V8[rownames(sign_AMF_V8) %in% CANB1_diff_A ,])
colnames(sign_AMF_V1)[1]<-"gene";colnames(sign_AMF_V4)[1]<-"gene";colnames(sign_AMF_V5)[1]<-"gene"
colnames(sign_AMF_V6)[1]<-"gene";colnames(sign_AMF_V8)[1]<-"gene"


all_AMFCANB1_pheatmap<-merge(sign_AMF_V1[c(1:2,6)],merge(sign_AMF_V4[c(1:2,6)],merge(sign_AMF_V5[c(1:2,6)],merge(sign_AMF_V6[c(1:2,6)],sign_AMF_V8[c(1:2,6)],by="gene"),by="gene"),by="gene"),by="gene")
colnames(all_AMFCANB1_pheatmap)[c(2,4,6,8,10)]<-c("V1_CANB1","V4_CANB1","V5_CANB1","V6_CANB1","V8_CANB1")
colnames(all_AMFCANB1_pheatmap)[c(3,5,7,9,11)]<-c("V1_adj.P.Val","V4_adj.P.Val","V5_adj.P.Val","V6_adj.P.Val","V8_adj.P.Val")

write.table(all_AMFCANB1_pheatmap,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/Supp.table5_AMF_CANB1FC4.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)
head(all_AMFCANB1_pheatmap)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/figure1c_CANB1_AMF.pdf")
pheatmap(t(all_AMFCANB1_pheatmap[,c(-1,-3,-5,-7,-9,-11)]),cellheight = 12,scale = 'column')
dev.off()


write.table(sign_AMF_V1,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/MAPman_AMF_V1.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)
write.table(sign_AMF_V4,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/MAPman_AMF_V4.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)
write.table(sign_AMF_V5,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/MAPman_AMF_V5.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)
write.table(sign_AMF_V6,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/MAPman_AMF_V6.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)
write.table(sign_AMF_V8,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/MAPman_AMF_V8.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)


###################################################

#AMF diff
sign_AMF_V1<-topTable(fit_AMF, coef=c(1),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V4<-topTable(fit_AMF, coef=c(2),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V5<-topTable(fit_AMF, coef=c(3),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V6<-topTable(fit_AMF, coef=c(4),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V8<-topTable(fit_AMF, coef=c(5),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL

sign_AMF_V1<-sign_AMF_V1[sign_AMF_V1$adj.P.Val<0.05,]
sign_AMF_V4<-sign_AMF_V4[sign_AMF_V4$adj.P.Val<0.05,]
sign_AMF_V5<-sign_AMF_V5[sign_AMF_V5$adj.P.Val<0.05,]
sign_AMF_V6<-sign_AMF_V6[sign_AMF_V6$adj.P.Val<0.05,]
sign_AMF_V8<-sign_AMF_V8[sign_AMF_V8$adj.P.Val<0.05,]

AMF_diff<-unique(c(rownames(sign_AMF_V1),
                   rownames(sign_AMF_V4),
                   rownames(sign_AMF_V5),
                   rownames(sign_AMF_V6),
                   rownames(sign_AMF_V8)))


#VAR treatment per AMF
sign_AMF_CAN<-topTable(fit_AMF, coef=c(6:15),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_B1<-topTable(fit_AMF, coef=c(16:25),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL

sign_AMF_CAN<-sign_AMF_CAN[sign_AMF_CAN$adj.P.Val<0.05,]
sign_AMF_B1<-sign_AMF_B1[sign_AMF_B1$adj.P.Val<0.05,]

##
# difference VAR

VAR_diffe<-unique(c(rownames(sign_AMF_CAN),rownames(sign_AMF_B1)))


########################################################################################
# only AMF and not VAR
AMF_diff[!AMF_diff %in% VAR_diffe]

sign_AMF_V1<-topTable(fit_AMF, coef=c(1),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V4<-topTable(fit_AMF, coef=c(2),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V5<-topTable(fit_AMF, coef=c(3),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V6<-topTable(fit_AMF, coef=c(4),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V8<-topTable(fit_AMF, coef=c(5),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V1<-cbind.data.frame(rownames(sign_AMF_V1[rownames(sign_AMF_V1) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,]),sign_AMF_V1[rownames(sign_AMF_V1) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,])
sign_AMF_V4<-cbind.data.frame(rownames(sign_AMF_V4[rownames(sign_AMF_V4) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,]),sign_AMF_V4[rownames(sign_AMF_V4) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,])
sign_AMF_V5<-cbind.data.frame(rownames(sign_AMF_V5[rownames(sign_AMF_V5) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,]),sign_AMF_V5[rownames(sign_AMF_V5) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,])
sign_AMF_V6<-cbind.data.frame(rownames(sign_AMF_V6[rownames(sign_AMF_V6) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,]),sign_AMF_V6[rownames(sign_AMF_V6) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,])
sign_AMF_V8<-cbind.data.frame(rownames(sign_AMF_V8[rownames(sign_AMF_V8) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,]),sign_AMF_V8[rownames(sign_AMF_V8) %in% AMF_diff[!AMF_diff %in% VAR_diffe] ,])
colnames(sign_AMF_V1)[1]<-"gene";colnames(sign_AMF_V4)[1]<-"gene";colnames(sign_AMF_V5)[1]<-"gene"
colnames(sign_AMF_V6)[1]<-"gene";colnames(sign_AMF_V8)[1]<-"gene"

AMF_AMF<-merge(sign_AMF_V1[c(1:2,6)],merge(sign_AMF_V4[c(1:2,6)],merge(sign_AMF_V5[c(1:2,6)],merge(sign_AMF_V6[c(1:2,6)],sign_AMF_V8[c(1:2,6)],by="gene"),by="gene"),by="gene"),by="gene")
colnames(AMF_AMF)[c(2,4,6,8,10)]<-c("V1_CANB1","V4_CANB1","V5_CANB1","V6_CANB1","V8_CANB1")
colnames(AMF_AMF)[c(3,5,7,9,11)]<-c("V1_adj.P.Val","V4_adj.P.Val","V5_adj.P.Val","V6_adj.P.Val","V8_adj.P.Val")


write.table(AMF_AMF,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/Supp.table5_AMF_AMF688.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)

########################################################################################

# only VAR and not AMF

VAR_diffe[!VAR_diffe %in% AMF_diff]


sign_AMF_CAN<-topTable(fit_AMF, coef=c(6:15),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_B1<-topTable(fit_AMF, coef=c(16:25),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_CAN<-cbind.data.frame(rownames(sign_AMF_CAN[rownames(sign_AMF_CAN) %in% VAR_diffe[!VAR_diffe %in% AMF_diff] ,]),sign_AMF_CAN[rownames(sign_AMF_CAN) %in% VAR_diffe[!VAR_diffe %in% AMF_diff] ,])
sign_AMF_B1<-cbind.data.frame(rownames(sign_AMF_B1[rownames(sign_AMF_B1) %in% VAR_diffe[!VAR_diffe %in% AMF_diff] ,]),sign_AMF_B1[rownames(sign_AMF_B1) %in% VAR_diffe[!VAR_diffe %in% AMF_diff] ,])
colnames(sign_AMF_CAN)[1]<-"gene";colnames(sign_AMF_B1)[1]<-"gene";
AMF_AMF<-merge(sign_AMF_CAN[c(1:11,15)],sign_AMF_B1[c(1:11,15)],by="gene")
colnames(AMF_AMF)[c(12,23)]<-c("CAN_adj.P.Val","B1_adj.P.Val")
write.table(AMF_AMF,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/Supp.table5_AMF_VAR1247.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)



########################################################################################
# Shared AMF and VAR
AMF_diff[AMF_diff %in% VAR_diffe]
VAR_diffe[VAR_diffe %in% AMF_diff]
sign_AMF_V1<-topTable(fit_AMF, coef=c(1),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V4<-topTable(fit_AMF, coef=c(2),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V5<-topTable(fit_AMF, coef=c(3),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V6<-topTable(fit_AMF, coef=c(4),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V8<-topTable(fit_AMF, coef=c(5),number=dim(DGE_AMF_N)[1]) # B1 vs CTRL
sign_AMF_V1<-cbind.data.frame(rownames(sign_AMF_V1[rownames(sign_AMF_V1) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,]),sign_AMF_V1[rownames(sign_AMF_V1) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,])
sign_AMF_V4<-cbind.data.frame(rownames(sign_AMF_V4[rownames(sign_AMF_V4) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,]),sign_AMF_V4[rownames(sign_AMF_V4) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,])
sign_AMF_V5<-cbind.data.frame(rownames(sign_AMF_V5[rownames(sign_AMF_V5) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,]),sign_AMF_V5[rownames(sign_AMF_V5) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,])
sign_AMF_V6<-cbind.data.frame(rownames(sign_AMF_V6[rownames(sign_AMF_V6) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,]),sign_AMF_V6[rownames(sign_AMF_V6) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,])
sign_AMF_V8<-cbind.data.frame(rownames(sign_AMF_V8[rownames(sign_AMF_V8) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,]),sign_AMF_V8[rownames(sign_AMF_V8) %in% VAR_diffe[VAR_diffe %in% AMF_diff] ,])
colnames(sign_AMF_V1)[1]<-"gene";colnames(sign_AMF_V4)[1]<-"gene";colnames(sign_AMF_V5)[1]<-"gene"
colnames(sign_AMF_V6)[1]<-"gene";colnames(sign_AMF_V8)[1]<-"gene"

AMF_AMFandVAR<-merge(sign_AMF_V1[c(1:2,6)],merge(sign_AMF_V4[c(1:2,6)],merge(sign_AMF_V5[c(1:2,6)],merge(sign_AMF_V6[c(1:2,6)],sign_AMF_V8[c(1:2,6)],by="gene"),by="gene"),by="gene"),by="gene")
colnames(AMF_AMFandVAR)[c(2,4,6,8,10)]<-c("V1_CANB1","V4_CANB1","V5_CANB1","V6_CANB1","V8_CANB1")
colnames(AMF_AMFandVAR)[c(3,5,7,9,11)]<-c("V1_adj.P.Val","V4_adj.P.Val","V5_adj.P.Val","V6_adj.P.Val","V8_adj.P.Val")


write.table(AMF_AMFandVAR,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/Supp.table5_AMF_AMF&VAR560.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T)



#########################################################################################
#########################################################################################
# important genes that differ between CAN B1 and VAR.

###AMF
sign_AMF_V1<-topTable(fit_AMF, coef=c(1),number=dim(DGE_AMF_N)[1]) # B1 vs CAN
sign_AMF_V4<-topTable(fit_AMF, coef=c(2),number=dim(DGE_AMF_N)[1]) # B1 vs CAN
sign_AMF_V5<-topTable(fit_AMF, coef=c(3),number=dim(DGE_AMF_N)[1]) # B1 vs CAN
sign_AMF_V6<-topTable(fit_AMF, coef=c(4),number=dim(DGE_AMF_N)[1]) # B1 vs CAN
sign_AMF_V8<-topTable(fit_AMF, coef=c(5),number=dim(DGE_AMF_N)[1]) # B1 vs CAN


# cassava genes FC > 8x between CAN B1.
sign_AMF_V1[sign_AMF_V1$logFC>3 | sign_AMF_V1$logFC<(-3),]
sign_AMF_V4[sign_AMF_V4$logFC>3 | sign_AMF_V4$logFC<(-3),]
sign_AMF_V5[sign_AMF_V5$logFC>3 | sign_AMF_V5$logFC<(-3),]
sign_AMF_V6[sign_AMF_V6$logFC>3 | sign_AMF_V6$logFC<(-3),]
sign_AMF_V8[sign_AMF_V8$logFC>3 | sign_AMF_V8$logFC<(-3),]

A_FC8<-unique(c(rownames(sign_AMF_V1[sign_AMF_V1$logFC>3 | sign_AMF_V1$logFC<(-3),]),
                rownames(sign_AMF_V4[sign_AMF_V4$logFC>3 | sign_AMF_V4$logFC<(-3),]),
                rownames(sign_AMF_V5[sign_AMF_V5$logFC>3 | sign_AMF_V5$logFC<(-3),]),
                rownames(sign_AMF_V6[sign_AMF_V6$logFC>3 | sign_AMF_V6$logFC<(-3),]),
                rownames(sign_AMF_V8[sign_AMF_V8$logFC>3 | sign_AMF_V8$logFC<(-3),]) ))


A_FC8_V1<-sign_AMF_V1[rownames(sign_AMF_V1) %in% A_FC8,]
A_FC8_V4<-sign_AMF_V4[rownames(sign_AMF_V4) %in% A_FC8,]
A_FC8_V5<-sign_AMF_V5[rownames(sign_AMF_V5) %in% A_FC8,]
A_FC8_V6<-sign_AMF_V6[rownames(sign_AMF_V6) %in% A_FC8,]
A_FC8_V8<-sign_AMF_V8[rownames(sign_AMF_V8) %in% A_FC8,]
A_FC8_V1$gene<-rownames(A_FC8_V1);A_FC8_V4$gene<-rownames(A_FC8_V4);A_FC8_V5$gene<-rownames(A_FC8_V5)
A_FC8_V6$gene<-rownames(A_FC8_V6);A_FC8_V8$gene<-rownames(A_FC8_V8)
A_FC8_ALL<-merge(A_FC8_V1,merge(A_FC8_V4,merge(A_FC8_V5,merge(A_FC8_V6,A_FC8_V8,by="gene"),by="gene"),by="gene"),by="gene")

A_FC8_ALL_CANB1<-A_FC8_ALL[,grep("logFC",colnames(A_FC8_ALL),invert=F)]
colnames(A_FC8_ALL_CANB1)<-c("V1_CAN.B1","V4_CAN.B1","V5_CAN.B1","V6_CAN.B1","V8_CAN.B1")
rownames(A_FC8_ALL_CANB1)<-A_FC8_ALL$gene
pheatmap(A_FC8_ALL_CANB1)


A_FC8_INFO<-merge(A_FC8_ALL,Mercator_Rirregularis,by="gene")
A_FC8_INFO<-A_FC8_INFO[!duplicated(A_FC8_INFO$gene),]
A_FC8_INFO[,35]<- unlist(lapply(lapply(strsplit(as.character(A_FC8_INFO[,33]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
A_FC8_INFO[,35]<-as.factor(A_FC8_INFO[,35])

table(A_FC8_INFO[,35])

# CHange function
A_FC8_function<-A_FC8_INFO[grep("stress.biotic",A_FC8_INFO[,35]),]
rownames(A_FC8_function)<-A_FC8_function$gene
pheatmap(A_FC8_function[,grep("logFC",colnames(A_FC8_function),invert=F)])

A_FC8_INFO[,c(1,34)]

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################


# test markers

genes_AMF<-sign_AMF[sign_AMF$adj.P.Val<0.05,]
genes_AMF_compe<-genes_AMF
genes_AMF_compe$gene<-rownames(genes_AMF_compe)
genes_AMF_compe<-merge(genes_AMF_compe,Mercator_Mesculenta,by="gene")
genes_AMF_compe[,9]<- unlist(lapply(lapply(strsplit(as.character(genes_AMF_compe[,9]),"\\."),function (x) x[1:2]),function (e) paste(e[1:2],collapse=".")))
genes_AMF_compe[,9]<-as.factor(genes_AMF_compe[,9])
#genes_AMF_compe<-genes_AMF_compe[16223,]

genes_AMF_compe2<-DGE_AMF$counts[rownames(DGE_AMF_N$E) %in% "Manes.S112700", ] 

genes_AMF_compe2<-DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% "Manes.S112700", ] 

rownames(Diff_CANB1[Diff_CANB1$adj.P.Val<0.05,])
genes_AMF_compe2<-DGE_AMF_N$E[rownames(DGE_AMF$counts) %in% rownames(sign_AMF[sign_AMF$V1_CAN.B1>2|sign_AMF$V4_CAN.B1>2|sign_AMF$V5_CAN.B1>2|
                                                                                    sign_AMF$V6_CAN.B1>2|sign_AMF$V8_CAN.B1>2,]), ] 

tratamiento<-unlist(lapply(strsplit(colnames(genes_AMF_compe2),"_"),function (x) x[3]))
cultivar<-gsub("V4","COL2215",
               gsub("V5","BRA337",
                    gsub("V6","CM4574-7",
                         gsub("m","",gsub("m2","",unlist(lapply(strsplit(colnames(genes_AMF_compe2),"_"),function (x) x[1]))   )))))
interacts<-interaction(tratamiento,cultivar)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/t_CAN_B1_diff_plant.pdf",width=10, height=10,useDingbats = F)
par(mfrow=c(3,3),mar=c(15,8,1,1),mgp=c(5, 1, 0))

for (i in 1:dim(genes_AMF_compe2)[1]) {
  stDevs <-tapply(genes_AMF_compe2[i,],interacts,sd)
  means<-tapply(genes_AMF_compe2[i,],interacts,mean)
  mp<-barplot(tapply(genes_AMF_compe2[i,],interacts,mean),ylim=c(min(means + stDevs)*2,max(means + stDevs)*2),
              beside=T,las =2,col=c("dodgerblue3","lightsalmon3","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
              ylab=paste("gene-expression\n",rownames(genes_AMF_compe2)[i]))
  segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
  segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
  segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
}
dev.off()

stDevs <-tapply(genes_AMF_compe2,interacts,sd)
means<-tapply(genes_AMF_compe2,interacts,mean)
mp<-barplot(tapply(genes_AMF_compe2,interacts,mean),
            beside=T,las =2,col=c("dodgerblue3","lightsalmon3","azure3"),cex.axis = 1.5,cex.names = 1.5,cex.lab=1.5,
            ylab=paste("gene-expression\n",rownames(genes_AMF_compe2)[i]))
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)






########################################################################################
########################################################################################
########################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
#################################################################################################################################################

#######################################################################################################################################

########################################################################################################################################



########################################################################################################################################
########################################    WGCNA MODULE COMPOSITIONS   ################################################################
########################################################################################################################################

## create eigene modules per sample. ##
###################################
#### M esculenta gene counts ######
###################################
YUCA_sign_TREAT<-unique(YUCA_byAMF_diff)
#normalized
sign_norm_genes_Y<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% YUCA_sign_TREAT,]
sign_norm_genes_Y<-sign_norm_genes_Y[,grep(c("CTRL"), colnames(sign_norm_genes_Y),invert =T)]
#pdf("~/Google Drive/Doctorat_shared_unil/Reports/FINAL_THESIS/Cassava_genes.pdf",useDingbats = F,width=10, height=10)
#pheatmap(sign_norm_genes_Y,show_rownames = F)
#dev.off()
#normalized data
dim(sign_norm_genes_Y) # 1634 genes.
colnames(t(FOUR_VARS_YUCA_S))
#create set of powers
powers=c(c(1:10),seq(from=12, to=30,by=2))
#network topology analysis fonction, show power of networks
sft_Y=pickSoftThreshold(t(sign_norm_genes_Y),powerVector = powers,verbose = 5)
par(mfrow=c(1,2))
plot(sft_Y$fitIndices[,1], -sign(sft_Y$fitIndices[,3])*sft_Y$fitIndices[,2],xlab="SoftThresholdpower)",ylab="ScaleFreeTopologyModelFit,signed^2",main =paste("Scaleindependence"),type="n")
text(sft_Y$fitIndices[,1], -sign(sft_Y$fitIndices[,3])*sft_Y$fitIndices[,2],labels=powers)

plot(sft_Y$fitIndices[,1], sft_Y$fitIndices[,5],type="n")
text(sft_Y$fitIndices[,1], sft_Y$fitIndices[,5],labels=powers)

####################################################
## clustering methods**
####################################################

# create adjacency,
softPower = 4
adjacency_Y = adjacency(t(sign_norm_genes_Y),power= softPower);
#creat topological overlap matrix
TOM_Y= TOMsimilarity(adjacency_Y)
dissTOM_Y=1-TOM_Y
#Clustering using TOM
geneTree_Y= hclust(as.dist(dissTOM_Y),method="average")
#sizeGrWindow(12,9)
#plot(geneTree_Y,main="Gene clustering based on TOM dissimilarity",labels=F,hang=0.04)
####################################################
#definition of module size.
####################################################
minModuleSize=50
#module identification using dynamic tree cut. module
dynmicMods_Y= cutreeDynamic(dendro=geneTree_Y, distM=dissTOM_Y,deepSplit=2,pamRespectsDendro=F,minClusterSize=minModuleSize)
table(dynmicMods_Y)

dynamicColors_Y=labels2colors(dynmicMods_Y)
table(dynamicColors_Y)
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree_Y,dynamicColors_Y,"Dynamic Tree Cut",dendroLabels= F,hang=0.03,addGuide = T,guideHang = 0.05,main="Gene dendrogram and module colors")
####################################################
#merging modules whose expression profiles are similar
####################################################
#calculate eignegens
MEList_Y=moduleEigengenes(t(sign_norm_genes_Y),colors=dynamicColors_Y)
MEs_Y=MEList_Y$eigengenes
#calculate dissimilarity od module eigengenes
MEDiss_Y= 1 -cor(MEs_Y)
#cluster module eigengenes
METree_Y=hclust(as.dist(MEDiss_Y),method="average")
sizeGrWindow(7,6)
plot(METree_Y, main="Clustering of module eigengenes") # choose cut for a correlation of 0.75

MEDissThres=0.1
abline(h=MEDissThres,col="red")
merge= mergeCloseModules(t(sign_norm_genes_Y),dynamicColors_Y,cutHeight = MEDissThres,verbose=3)
mergedColors=merge$colors
mergedMEs_Y=merge$newMEs
# to see effect of merging
sizeGrWindow(12,9)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA1_modules_CASSAVA.pdf")
par(las=1)
plotDendroAndColors(geneTree_Y,mergedColors,"Modules",dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05,
                    cex.colorLabels = 2, cex.dendroLabels = 2, 
                    cex.rowText = 2)
dev.off()
#change names (to work better)
moduleColors_Y=mergedColors
colorOrder=c("grey",standardColors(50))
moduleLabels_Y=match(moduleColors_Y,colorOrder)-1
MEs_Y=mergedMEs_Y

###########################################################################################################################
# PREPARE DATA YUCA FOR quantifying module-module associations 

#define nb genes and samples
nGenes=ncol(t(sign_norm_genes_Y))
nSamples=nrow(t(sign_norm_genes_Y))
#recalculate MEs with color labels
MEs0_Y=moduleEigengenes(t(sign_norm_genes_Y),moduleColors_Y)$eigengenes
MEs_Y=orderMEs(MEs0_Y)
names(MEs_Y)<-paste("Y",names(MEs_Y),sep="_")
rownames(MEs_Y)<-rownames(t(sign_norm_genes_Y))


###########################################################################################################################
###########################################################################################################################
###################################
#### R irregularis gene modules####
###################################

#normalized data
sign_norm_genes_A<-DGE_AMF_N$E
dim(sign_norm_genes_A)
#pdf("~/Google Drive/Doctorat_shared_unil/Reports/FINAL_THESIS/AMF_genes.pdf",useDingbats = F,width=10, height=10)
#pheatmap(sign_norm_genes_A,show_rownames = F)
#dev.off()

## create eigene modules per sample. 
#create set of powers
powers=c(c(1:10),seq(from=12, to=40,by=2))
#network topology analysis fx
sft_A=pickSoftThreshold(t(sign_norm_genes_A),powerVector = powers,verbose = 5) #FOUR_VARS_AMF_S
par(mfrow=c(1,2))
plot(sft_A$fitIndices[,1], -sign(sft_A$fitIndices[,3])*sft_A$fitIndices[,2],xlab="SoftThresholdpower)",ylab="ScaleFreeTopologyModelFit,signed^2",main =paste("Scaleindependence"),type="n")
text(sft_A$fitIndices[,1], -sign(sft_A$fitIndices[,3])*sft_A$fitIndices[,2],labels=powers)

plot(sft_A$fitIndices[,1], sft_A$fitIndices[,5],type="n")
text(sft_A$fitIndices[,1], sft_A$fitIndices[,5],labels=powers)

####################################################
## clustering methods**
####################################################
#adjacency
adjacency_A = adjacency(t(sign_norm_genes_A),power= 8); # change power FOUR_VARS_AMF_S
#topological overlap dissimilarity matrix 
dissTOM_A=1-TOMsimilarity(adjacency_A)

#Clustering using TOM
geneTree_A= hclust(as.dist(dissTOM_A),method="average")
#sizeGrWindow(12,9)
#plot(geneTree_A,main="Gene clustering based on TOM dissimilarity",labels=F,hang=0.04)

####################################################
#definition of module size.
####################################################
minModuleSize=50
#module identification using dynamic tree cut
dynmicMods_A= cutreeDynamic(dendro=geneTree_A, distM=dissTOM_A,deepSplit=2,pamRespectsDendro=F,minClusterSize=minModuleSize)
dynamicColors_A=labels2colors(dynmicMods_A)
table(dynamicColors_A)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree_A,dynamicColors_A,"Dynamic Tree Cut",dendroLabels= F,hang=0.03,addGuide = T,guideHang = 0.05,main="Gene dendrogram and module colors")

#merging modules whose expression profiles are similar
#calculate eignegens
MEList_AMF=moduleEigengenes(t(sign_norm_genes_A),colors=dynamicColors_A) #FOUR_VARS_AMF_S
MEs_AMF=MEList_AMF$eigengenes
#calculate dissimilarity od module eigengenes
MEDiss_AMF= 1 -cor(MEs_AMF)
#cluster module eigengenes
METree_AMF=hclust(as.dist(MEDiss_AMF),method="average")
sizeGrWindow(7,6)
plot(METree_AMF, main="Clustering of module eigengenes") # choose cut for a correlation of 0.75

MEDissThres=0.1
abline(h=MEDissThres,col="red")
merge_AMF= mergeCloseModules(t(sign_norm_genes_A),dynamicColors_A,cutHeight = MEDissThres,verbose=3)
mergedColors_AMF=merge_AMF$colors
mergedMEs_AMF=merge_AMF$newMEs
# to see effect of merging
sizeGrWindow(12,9)
plotDendroAndColors(geneTree_A,cbind(dynamicColors_A,mergedColors_AMF),c("Dynamic Tree cut","Merged Dynamic"),dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05)
table(mergedColors_AMF)
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA2_modules_AMF6.pdf")
par(las=1)
plotDendroAndColors(geneTree_A,mergedColors_AMF,"Modules",dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05,
                   cex.dendroLabels = 2, cex.colorLabels = 2,
                    cex.rowText = 2)
dev.off()


#change names (to work better)
moduleColors_AMF=mergedColors_AMF
colorOrder=c("grey",standardColors(50))
moduleLabels_AMF=match(moduleColors_AMF,colorOrder)-1
MEs_AMF=mergedMEs_AMF

###########################################################################################################################
# PREPARE DATA AMF FOR quantifying module-module associations 

#define nb genes and samples
nGenes=ncol(t(sign_norm_genes_A))
nSamples=nrow(t(sign_norm_genes_A))
#recalculate MEs with color labels
MEs0_AMF=moduleEigengenes(t(sign_norm_genes_A),moduleColors_AMF)$eigengenes

MEs_AMF=orderMEs(MEs0_AMF)
names(MEs_AMF)<-paste("A",names(MEs_AMF),sep="_")
rownames(MEs_AMF)<-rownames(t(sign_norm_genes_A))
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

#Modules expression in samples
##########################################################################################################################################################
# EIGENGENES
# MOdules expression in different treatments.
# Mesculenta
########??RELATE MODULE EIGENGENES APPARTENSHIP TO TREATMENTS.
modules_YUCA<-MEs_Y

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA3_module_eigenenes_CASSAVA.pdf")
#rownames(modules_YUCA)<-rownames(t(FOUR_VARS_YUCA_topheno))
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
heatmap.2( as.matrix(modules_YUCA), col=colors,scale="none",dendrogram = "both",key = F)
dev.off()
heatmap.2( as.matrix(modules_YUCA), col=colors,scale="column",dendrogram = "both",key = F)


modules_AMF<-MEs_AMF

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA4_module_eigenenes_AMF.pdf")
#rownames(modules_AMF)<-rownames(t(FOUR_VARS_AMF_topheno))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( as.matrix(modules_AMF), col=colors,scale="none",dendrogram = "both", key = F )
dev.off()



##################################################################################################
####correlation Modules eigengenes YUCA,AMF #####
##################################################################################################
#takeout samples from CTRL in cassava.



moduleTraitCor= cor(MEs_Y,MEs_AMF,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)


#Representation
sizeGrWindow(10,8)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",
                 signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))
#disply correlation values within heatmap plot
#X-axis co-expresion modules on AMF, Y-axis co-expresion modules CASSAVA
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA5_CORRELATION_Modules_CASSAVA_AMF.pdf",height = 8,width = 8)
col_HM<-c(colorRampPalette( rev(brewer.pal(9, "Greens")) )(20),rep("white",10),colorRampPalette( brewer.pal(9, "Blues") )(20)  )
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MEs_AMF),
               yLabels = names(MEs_Y),
               ySymbols = names(MEs_Y),
               xSymbols= names(MEs_AMF),
               colorLabels = FALSE,
               colors = col_HM,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,cex.lab=1,
               zlim = c(-1,1),
               main = paste("Module YUCA-module AMF relationships"))
dev.off()





###########################################################################################################################
########### CORRELATION MODULES CASSSAVA AMF

table(mergedColors_AMF) #AMF modules
table(mergedColors) #yuca modules

barplot(MEs_AMF[,4],col="black")
moduleTraitCor= cor(MEs_Y,MEs_AMF,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)
dim(moduleTraitPvalue)[2]
corModules2<-as.vector(as.matrix(moduleTraitPvalue))# TO CHANGE DEPENDING ON SIZE CORRELATION MATRIX
names(corModules2)<-c(paste(rep(colnames(moduleTraitPvalue)[1],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[2],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[3],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[4],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[5],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[6],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[7],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[8],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[9],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[10],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[11],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[12],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[13],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[14],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[15],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[16],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[17],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[18],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[19],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue)[20],dim(moduleTraitPvalue)[1]),rownames(moduleTraitPvalue),sep=" ")
                     
)
corModules_S<-corModules2[corModules2<0.001]
corModules_S
Cor_MOD_CASSAVA_2_AMF<-corModules_S
################################################################################################
################################################################################################
# expression of pair of modules in same plot.
Cor_MOD_CASSAVA_2_AMF

names(Cor_MOD_CASSAVA_2_AMF)[1]
sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[1]," "), "[[", 2)# YUCA
sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[1]," "), "[[", 1) # aMF

for (i in 1:length(names(Cor_MOD_CASSAVA_2_AMF))) {
cor_2_Modules_eigen<-data.frame(modules_YUCA[,grep(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[i]," "), "[[", 2),colnames(modules_YUCA))]
,modules_AMF[,grep(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[i]," "), "[[", 1),colnames(modules_AMF))]
)
rownames(cor_2_Modules_eigen)<-rownames(modules_YUCA)
colnames(cor_2_Modules_eigen)<-c(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[i]," "), "[[", 2),sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[i]," "), "[[", 1))

pdf(paste("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/COR2ORG_MODULES",names(Cor_MOD_CASSAVA_2_AMF)[i],".pdf",sep="_") )
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
heatmap.2( t(as.matrix(cor_2_Modules_eigen)),margins =c(8,8),
           notecol="black",col=colors,scale="row",trace ="none",notecex=.5,cexRow = 1,dendrogram = "col",key = F) #??change to row dendogram
dev.off()
}
# problem with greeping names that are included in others green greenyellow

cor_2_Modules_eigen<-data.frame(modules_YUCA[,grep(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[37]," "), "[[", 2),colnames(modules_YUCA))]
                                ,modules_AMF[,grep(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[37]," "), "[[", 1),colnames(modules_AMF))]
)
grep(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)[37]," "), "[[", 2),colnames(modules_YUCA))
colnames(modules_YUCA)

###########################################################################################################################
# Network of biological process inovolved in AMF-CASSAVA mOdules correlation
###########################################################################################################################
#install.packages("igraph")
library(igraph) 

nod<-cbind.data.frame(c(paste(rep("Y",length(table(merge$colors))),names(table(merge$colors)),sep="_"),
                        paste(rep("A",length(table(merge_AMF$colors))),names(table(merge_AMF$colors)) ,sep="_") ), 
                      c(names(table(merge$colors)),names(table(merge_AMF$colors))),
log(c(as.vector(table(merge$colors)),as.vector(table(merge_AMF$colors))) ),
c(as.vector(table(merge$colors)),as.vector(table(merge_AMF$colors))),
c(rep("CASSAVA",length(table(merge$colors))),rep("AMF",length(table(merge_AMF$colors))) ),
c(rep("1",length(table(merge$colors))),rep("2",length(table(merge_AMF$colors))) ))

colnames(nod)<-c("id","module","size","nsize","org","org.type")


lin<-cbind.data.frame(
gsub("ME","",sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)," "), "[[", 1)),
gsub("ME","",sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)," "), "[[", 2)),
-log10(as.vector(Cor_MOD_CASSAVA_2_AMF)   )   )

colnames(lin)<-c("from","to","weight")


net <- graph_from_data_frame(d=lin, vertices=nod, directed=F) 

# Generate colors based on media type:
V(net)$color <- c(rep("grey",length(table(merge$colors))),rep("brown",length(table(merge_AMF$colors))) )


# Set node size based on audience size:
V(net)$size <- V(net)$size*2
V(net)$label.cex = 1
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGNCA6_MODULES_NETWORK.pdf',useDingbats = F)

plot(net,vertex.label.dist=-0.4)
dev.off()


###################################################################################################
####################################################################

################################################################################################################################
# GO ANALYSIS IN MODULES
################################################################################################################################

#sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,]),]
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes_Y
top_gene_sign_intramod_conec_Y<-list()
GO_CASSAVA_Sign<-list()
for (correlated_Modula in   names(corModules_S))  {
  names(corModules_S)[1]
  gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)]#gsub("A_","",
  MEs_AMF$A_MEred
  ######################################################
  # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
  ######################################################
  #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
  geneModuleMembership_YUCA=as.data.frame(cor(t(sign_norm_genes),MEs_Y,use="p"))
  geneTraitSignificance_YUCA= as.data.frame(cor(t(sign_norm_genes), 
                                                MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)],use="p")) # define correlation with external variable gsub("A_","",
  ######################################################
  # 1) MODULE DEFINITION
  ######################################################
  #Cassava
  ## CHANGE FOR THE COMPARISON DESIRED
  modNames=substring(names(MEs_Y),5)
  module=gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))
  column= match(module,modNames)
  moduleGenes=moduleColors_Y==module
  ######################################################
  # 2) ANOTATION DATABASE
  ######################################################
  ## CASSAVA
  locusName<-rownames(sign_YUCA)
  sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
  sign_feat_go<-merge(sign_feat2,mesculenta_go2,by='locusName')
  sign_feat_go<-na.omit(sign_feat_go)
  gene_2_go<-sign_feat_go[,c(1,10)]
  head(sign_feat_go)
  write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
              col.names = F, row.names = F)
  geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
  ######################################################
  # 3) SELECT UNIVERSE
  ######################################################
  
  #TOTAL UNIVERSE ALL LOCI
  
  # TOTAL UNIVERSE CASSAVA
  locusName<-rownames(sign_YUCA)
  sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
  sign_feat_go<-merge(sign_feat2,mesculenta_go2[,c(1,3)],by='locusName') # do not select transcriptname. working at gene level
  sign_feat_go<-na.omit(sign_feat_go)
  sign_feat_go<-unique(sign_feat_go)
  length(unique(sign_feat_go$locusName))
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
  ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #CASSAVA
  GS.6_Y<-geneModuleMembership_YUCA[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0 & abs(geneModuleMembership_YUCA[moduleGenes, match(sapply(strsplit(correlated_Modula," "), "[[", 2), names(geneModuleMembership_YUCA))])>0,]
  dim(GS.6_Y)
  
  In_module<-rownames(GS.6_Y) # genes name in module CASSAVA
  
  most_correlated_genes<-geneModuleMembership_YUCA[moduleGenes,]
  
  na.omit(most_correlated_genes[,grep(sapply(strsplit(correlated_Modula," "), "[[", 2), colnames(most_correlated_genes)  )])
  
  
  
  
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  
  geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
  names(geneList) <- sign_feat_go$locusName
  
  sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5, ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
  
  
  sampleGOdata
  
  
  resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  RES_Y<-allRes
  
  GO_CASSAVA_Sign[[correlated_Modula]]<-RES_Y
  
  write.table(RES_Y,  paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/CROSSTALKS_AMF_CASSAVA/GO_ENRICH_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'GO_ENRICH.txt',sep="_") ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  write.table(RES_Y[RES_Y[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/CROSSTALKS_AMF_CASSAVA/GO_2REVIGO_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'GO_ENRICH.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = F, row.names = F)
  
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/CROSSTALKS_AMF_CASSAVA/BARPLOT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'.pdf',sep="_"), width=10, height=7)
  par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  tryCatch( barplot(RES_Y$Significant[RES_Y$classicFisher<0.05], names = paste(RES_Y$GO.ID[RES_Y$classicFisher<0.05],RES_Y$Term[RES_Y$classicFisher<0.05],sep=" "),
                    xlab = "Nb. of significant genes",horiz=T,las=1),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 2)) )
  dev.off()
  
  
  to_pheat_Y<-RES_Y[RES_Y$classicFisher<0.05,c(4)]
  names(to_pheat_Y)<-paste(RES_Y[RES_Y$classicFisher<0.05,c(1)],RES_Y[RES_Y$classicFisher<0.05,c(2)],sep=" ")
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/CROSSTALKS_AMF_CASSAVA/PHEAT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'.pdf',sep="_"), width=18, height=10)
  tryCatch(pheatmap(t(to_pheat_Y),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 2)))
  dev.off()
}

######################################################  ######################################################
# IN AMF
######################################################  ######################################################


sign_norm_genes<-DGE_AMF_N$E
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
top_gene_sign_intramod_conec_A<-list()
GO_AMF_Sign<-list()
for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_AMF)[c(1,101,119,173,19,55,37,165,83,47)])  {
  names(corModules_S)[1]
  gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)]
  
  
  ######################################################
  # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
  ######################################################
  #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
  geneModuleMembership_AMF=as.data.frame(cor(t(sign_norm_genes),MEs_AMF,use="p"))
  geneTraitSignificance_AMF= as.data.frame(cor(t(sign_norm_genes),
                                               MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)],use="p")) # define correlation with external variable
  ######################################################
  # 1) MODULE DEFINITION
  ######################################################
  #AMF
  modNames=substring(names(MEs_AMF),5)
  module=gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  column= match(module,modNames)
  moduleGenes=moduleColors_AMF==module
  ######################################################
  # 2) ANOTATION DATABASE
  ######################################################
  #AMF
  locusName<-rownames(sign_AMF)
  sign_feat2<-cbind.data.frame(sign_AMF,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
  sign_feat_go<-na.omit(sign_feat_go)
  gene_2_go<-sign_feat_go[,c(1,31)]
  
  write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
              col.names = F, row.names = F)
  geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
  ######################################################
  # 3) SELECT UNIVERSE
  ######################################################
  #TOTAL UNIVERSE AMF
  locusName<-rownames(sign_AMF)
  sign_feat2<-cbind.data.frame(sign_AMF,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName') # do not select transcriptname. working at gene level
  sign_feat_go<-na.omit(sign_feat_go)
  sign_feat_go<-unique(sign_feat_go)
  length(unique(sign_feat_go$locusName))
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
  ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #AMF
  # change module match
  GS.6_A<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_AMF[moduleGenes,1])>0 & abs(geneModuleMembership_AMF[moduleGenes,match(substring(sapply(strsplit(correlated_Modula," "), "[[", 1),1), names(geneModuleMembership_AMF))])>0,]
  dim(GS.6_A)
  In_module<-rownames(GS.6_A) # genes name in module AMF
  
  ######################################################
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  
  geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
  names(geneList) <- sign_feat_go$locusName
  
  sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5, ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
  
  
  sampleGOdata
  
  
  resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  
  RES_A<-allRes
  GO_AMF_Sign[[correlated_Modula]]<-RES_A
  
  
  ######################################################
  # 6) PRINT RESULTS
  ######################################################
  
  
  ################
  
  
  
  
  write.table(RES_A,    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/CROSSTALKS_AMF_CASSAVA/GO_ENRICH_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'GO_ENRICH.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  write.table(RES_A[RES_A[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/CROSSTALKS_AMF_CASSAVA/GO_2REVIGO_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'GO_ENRICH.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = F, row.names = F)
  
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/CROSSTALKS_AMF_CASSAVA/BARPLOT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.pdf',sep="_"), width=10, height=7)
  par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  tryCatch( barplot(RES_A$Significant[RES_A$classicFisher<0.05], names = paste(RES_A$GO.ID[RES_A$classicFisher<0.05],RES_A$Term[RES_A$classicFisher<0.05],sep=" "),
                    xlab = "Nb. of significant genes",horiz=T,las=1),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 1)) )
  dev.off()
  
  
  to_pheat_A<-RES_A[RES_A$classicFisher<0.05,c(4)]
  names(to_pheat_A)<-paste(RES_A[RES_A$classicFisher<0.05,c(1)],RES_A[RES_A$classicFisher<0.05,c(2)],sep=" ")
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/CROSSTALKS_AMF_CASSAVA/PHEAT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.pdf',sep="_"), width=18, height=10)
  tryCatch(pheatmap(t(to_pheat_A),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 1)))
  dev.off()
  
  
}


################################################################################################################################
# extract more representative and more significant genes, no go terms
#fastest
################################################################################################################################
# CASSAVA


sign_norm_genes<-sign_norm_genes_Y
top_gene_sign_intramod_conec_Y<-list()
Import_gene_M_Y<-list()
for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_AMF))  {
  
  ######################################################
  # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
  ######################################################
  #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
  geneModuleMembership_YUCA=as.data.frame(cor(t(sign_norm_genes_Y),MEs_Y,use="p"))
  geneTraitSignificance_YUCA= as.data.frame(cor(t(sign_norm_genes_Y), 
                                                MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)],use="p")) # define correlation with external variable gsub("A_","",
  ######################################################
  # 1) MODULE DEFINITION
  ######################################################
  #Cassava
  ## CHANGE FOR THE COMPARISON DESIRED
  modNames=substring(names(MEs_Y),5)
  module=gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))
  column= match(module,modNames)
  moduleGenes=moduleColors_Y==module
  ######################################################
  
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
  ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #CASSAVA
  GS.6_Y<-geneModuleMembership_YUCA[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0 & abs(geneModuleMembership_YUCA[moduleGenes, match(sapply(strsplit(correlated_Modula," "), "[[", 2), names(geneModuleMembership_YUCA))])>0,]
  dim(GS.6_Y)
  
  In_module<-rownames(GS.6_Y) # genes name in module CASSAVA
  most_correlated_genes<-geneModuleMembership_YUCA[moduleGenes,]
  na.omit(most_correlated_genes[,grep(sapply(strsplit(correlated_Modula," "), "[[", 2), colnames(most_correlated_genes)  )])
  
  ######################################################
  # 6) PRINT RESULTS
  ######################################################
  ## intramodular conectivity
  ADJ1=abs(cor(t(sign_norm_genes_Y),use="p"))^6
  Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors_Y)
  head(Alldegrees1)
  
  GS1=as.numeric(cor(MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)],t(sign_norm_genes_Y), use="p"))
  GeneSignificance=abs(GS1)
  colorlevels=unique(dynamicColors_Y)
  
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/ALL_MODULES_SIGN_Y',correlated_Modula,'.pdf',sep="_"),width=14,height=10,useDingbats = FALSE )
  par(mar=c(8,9,5,2),las=1,mgp=c(6, 2,0))
  whichmodule=gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2));
  restrict1 = (dynamicColors_Y==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=dynamicColors_Y[restrict1],
                     main=whichmodule,pch=16,cex=3, cex.axis=3,cex.lab = 3,cex.main = 3,lwd=10,
                     xlab = paste("Connectivity",whichmodule), ylab = paste(paste("Y",whichmodule,sep="_"),"correlation","to",gsub("A_ME","A_",sapply(strsplit(correlated_Modula," "), "[[", 1)),sep=" "), abline = TRUE)
  dev.off()
  
  datKME=signedKME(t(sign_norm_genes_Y), MEs_Y, outputColumnName="MM.")
  #FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub("Y_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 2)),colnames(datKME))]) >.8
  FilterGenes=  abs(GS1)[moduleGenes]> quantile(abs(GS1)[moduleGenes],.9) & abs(datKME[moduleGenes,grep(gsub("Y_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 2)),colnames(datKME))]) >.8
  
  top_sig_connect_high<-cbind.data.frame(rownames(datKME[moduleGenes,][FilterGenes,]),datKME[moduleGenes,grep(gsub("Y_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 2)),colnames(datKME))][FilterGenes])
 gene_cor<-data.frame(GS1)
  rownames(gene_cor)<- rownames(sign_norm_genes_Y)
  gene_cor<-cbind.data.frame(gene_cor,rownames(sign_norm_genes_Y))
  top_sig_correlated<-cbind.data.frame(rownames(gene_cor[moduleGenes,][FilterGenes,]),gene_cor[moduleGenes,][FilterGenes,])
  top_sig_correlated[,2:3]
  
  table(FilterGenes)
  top_gene_sign_intramod_conec_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-as.vector(top_sig_connect_high[,1]) #dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
 
  #sign_YUCA_V1_ALL<-topTable(fit_YUCA, coef=c(1:3),number=dim(DGE_YUCA_N)[1]) # B1 vs CTRL
  
  #Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-sign_YUCA_V1_ALL[rownames(sign_YUCA_V1_ALL) %in%   as.vector(unlist(top_gene_sign_intramod_conec_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]])),]
  
 #now create table swith top gene and top significancy
  #Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-cbind.data.frame(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]],rownames(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]))
  
  Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-cbind.data.frame(top_sig_correlated[,2:3], top_sig_connect_high[,2])
    
      

#  corr_list<-list()
 # pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/AMFgenes2pheno',correlated_Modula,'.pdf',sep=""),width=14,height=10,useDingbats = FALSE)
  #par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0),mfrow=c(3,3))
  #for( i in 1:length(top_gene_amf_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])){
   # NORM_EXP_TOPAMF<-t(sign_norm_genes_A[rownames(sign_norm_genes_A) %in% top_gene_amf_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],])
  #  corr_list[[i]]<-cor(phenoA[,c(6)],as.vector(NORM_EXP_TOPAMF[,i]),use="pairwise.complete.obs")
   # plot(NORM_EXP_TOPAMF[,i],phenoA[,c(6)],pch=16,cex=3,#col=phenoA$var,
    #     ylab=sapply(strsplit(correlated_Modula," "), "[[", 2),xlab=top_gene_amf_pheno[[1]][i],cex.lab=3 ,cex.axis=3)
    #abline(lm(phenoA[,c(6)]~as.vector(NORM_EXP_TOPAMF[,i])))
#  }
 # dev.off()
 # Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<- cbind.data.frame( unlist(corr_list),  Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])
#  colnames(Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[1]<-"Correlation"
  
  ####
  
  
  colnames(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]])<-c("gene_significance","gene","gene_correlation_to_module")
  Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-merge(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]],Mercator_Mesculenta,by="gene")
  
  
  toprint<-Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]
  toprint<-cbind( rep(sapply(strsplit(correlated_Modula," "), "[[", 1),dim(  toprint)[1]),
                  rep(sapply(strsplit(correlated_Modula," "), "[[", 2),dim(  toprint)[1]),
                  toprint[,c(1:6)])
  colnames(toprint)[1:2]<-c("AMF_module","YUCA_module")
  
  write.table( toprint,
               paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/TOP_GENES_FUNCTION',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 1),'.txt',sep="_"),
               ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  ################
  
  MEList_YUCA=moduleEigengenes(t(sign_norm_genes),colors=dynamicColors_Y)
  exp_and_modules<-cbind.data.frame(sign_norm_genes,MEList_YUCA$validColors)
  colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
  
  # extract top 20 genes
  exp_and_modules[exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2)),]
  cor_to_module<-as.data.frame(cor(t(sign_norm_genes_Y),
                                   MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)]
                                   ,use="p"))
  top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
  colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
  top_correlated2<-top_correlated[ top_correlated$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
  
  V_YUCACANB1_pheatmap
  #### all markers + p-value for MapMan
  COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
  colnames(COR_info_module)<-c("gene",paste("cor2",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="_"))
  FC_YUCA<-cbind.data.frame(V_YUCACANB1_pheatmap)
  colnames(FC_YUCA)[1]<-'gene'
  COR_info_module<-merge(merge(COR_info_module,Mercator_Mesculenta,by="gene"),FC_YUCA,by='gene')
  write.table(  COR_info_module,
                paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/GENES_2mapman',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 1),'.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  
  
  #do barplot
  exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
  Genes_module<-t(sign_norm_genes)[,colnames(t(sign_norm_genes)) %in%  rownames(exp_mod2)]
  write.table(t(Genes_module),    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/ALL_GENES_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'GO_ENRICH.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = T)
}

top_gene_sign_intramod_conec_Y
unlist(lapply(top_gene_sign_intramod_conec_Y,length))
table(mergedColors) #yuca modules
names(Import_gene_M_Y)
Import_gene_M_Y[2]
# all important genes in all modules
Interestingy<-sign_YUCA[rownames(sign_YUCA) %in%   as.vector(unlist(top_gene_sign_intramod_conec_Y)),]
Interestingy<-cbind.data.frame(Interestingy,rownames(Interestingy))
colnames(Interestingy)[dim(Interestingy)[2]]<-"gene"
Interestingy<-merge(Interestingy,Mercator_Mesculenta,by="gene")

top_gene_sign_intramod_conec_Y
######################################################  ######################################################
# IN AMF
######################################################  ######################################################




top_gene_sign_intramod_conec_A<-list()
Import_gene_M_A<-list()

for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_AMF))  {
 
  ######################################################
  # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
  ######################################################
  #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
  geneModuleMembership_AMF=as.data.frame(cor(t(sign_norm_genes_A),MEs_AMF,use="p"))
  geneTraitSignificance_AMF= as.data.frame(cor(t(sign_norm_genes_A),
                                               MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)],use="p")) # define correlation with external variable
  ######################################################
  # 1) MODULE DEFINITION
  ######################################################
  #AMF
  modNames=substring(names(MEs_AMF),5)
  module=gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  column= match(module,modNames)
  moduleGenes=moduleColors_AMF==module
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
  ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #AMF
  # change module match
  GS.6_A<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_AMF[moduleGenes,1])>0 & abs(geneModuleMembership_AMF[moduleGenes,match(substring(sapply(strsplit(correlated_Modula," "), "[[", 1),1), names(geneModuleMembership_AMF))])>0,]
  dim(GS.6_A)
  In_module<-rownames(GS.6_A) # genes name in module AMF
  
  ######################################################
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  ######################################################
  # 6) PRINT RESULTS
  ######################################################
  ## intramodular conectivity
  ADJ1=abs(cor(t(sign_norm_genes_A),use="p"))^6
  Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors_A)
  head(Alldegrees1)
  
  GS1=as.numeric(cor(MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)],t(sign_norm_genes_A), use="p"))
  GeneSignificance=abs(GS1)
  #GS.6_Y<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0]
  # quantile(abs(geneTraitSignificance_YUCA[moduleGenes,1]),.85)
  table( abs(GS1[moduleGenes])> quantile(abs(GS1[moduleGenes]),.9))
  colorlevels=unique(dynamicColors_A)
 
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/ALL_MODULES_SIGN',correlated_Modula,'.pdf',sep="_"),width=14,height=10,useDingbats = FALSE)
  par(mar=c(8,9,5,2),las=1,mgp=c(6, 2,0))
  
  whichmodule=gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1));
  restrict1 = (dynamicColors_A==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=dynamicColors_A[restrict1],
                     main=whichmodule,pch=16,cex=3, cex.axis=3,cex.lab = 3,cex.main = 3,
                     xlab = paste("Connectivity",whichmodule), ylab = paste(paste("A",whichmodule,sep="_"),"correlation","to",gsub("Y_ME","Y_",sapply(strsplit(correlated_Modula," "), "[[", 2)),sep=" "), abline = TRUE)
  
  dev.off()
  
  
  ####### # 
  #redo filter but with quantile and not absolute value
  #module black
  ##########
  datKME=signedKME(t(sign_norm_genes_A), MEs_AMF, outputColumnName="MM.")
  # FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
  FilterGenes=  abs(GS1)[moduleGenes]> quantile(abs(GS1)[moduleGenes],.9) & abs(datKME[moduleGenes,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
  
  
  
  top_sig_connect_high<-cbind.data.frame(rownames(datKME[moduleGenes,][FilterGenes,]),datKME[moduleGenes,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))][FilterGenes])
  gene_cor<-data.frame(GS1)
  rownames(gene_cor)<- rownames(sign_norm_genes_A)
  gene_cor<-cbind.data.frame(gene_cor,rownames(sign_norm_genes_A))
  top_sig_correlated<-cbind.data.frame(rownames(gene_cor[moduleGenes,][FilterGenes,]),gene_cor[moduleGenes,][FilterGenes,])
  top_sig_correlated[,2:3]
  
  table(FilterGenes)
  top_gene_sign_intramod_conec_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-as.vector(top_sig_connect_high[,1]) #dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
  
  Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-cbind.data.frame(top_sig_correlated[,2:3], top_sig_connect_high[,2])

  colnames(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])<-c("gene_significance","gene","gene_correlation_to_module")
Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-  merge(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],Mercator_Rirregularis,by="gene")
  


  toprint<-Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]
  toprint<-cbind( rep(sapply(strsplit(correlated_Modula," "), "[[", 1),dim(  toprint)[1]),
         rep(sapply(strsplit(correlated_Modula," "), "[[", 2),dim(  toprint)[1]),
         toprint[,c(1:6)])
  colnames(toprint)[1:2]<-c("AMF_module","YUCA_module")
  write.table( toprint,
               paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/TOP_GENES_FUNCTION',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_"),
               ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  summary( toprint)
  
  paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")
  
  ################
  
  MEList_AMF=moduleEigengenes(t(sign_norm_genes_A),colors=dynamicColors_A)
  exp_and_modules<-cbind.data.frame(sign_norm_genes_A,MEList_AMF$validColors)
  colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
  # extract top 20 genes
  exp_and_modules[exp_and_modules$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1)),]
  cor_to_module<-as.data.frame(cor(t(sign_norm_genes_A),
                                   MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)]
                                   ,use="p"))
  top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
  colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
  
  top_correlated2<-top_correlated[ top_correlated$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))  ,]
  
  
  #### all markers + p-value for MapMan
  COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
  colnames(COR_info_module)<-c("gene",paste("cor2",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="_"))
  FC_AMF<-cbind.data.frame(rownames(sign_AMF),sign_AMF)
  colnames(FC_AMF)[1]<-'gene'
  COR_info_module<-merge(merge(COR_info_module,Mercator_Rirregularis,by="gene"),FC_AMF,by='gene')
  
  write.table(  COR_info_module,
                paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/GENES_2mapman',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  
  #extract all genes in module
  exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
  
  Genes_module<-t(sign_norm_genes_A)[,colnames(t(sign_norm_genes_A)) %in%  rownames(exp_mod2)]
  
  write.table(t(Genes_module), paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_corA_Y/ALL_GENES_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = T)
}

names(Import_gene_M_A)
top_gene_sign_intramod_conec_A




#####################################################################
#nb of positive negative correlations between modules
names(top_gene_sign_intramod_conec_A)
names(top_gene_sign_intramod_conec_Y)
unlist(top_gene_sign_intramod_conec_A[[1]])
unlist(top_gene_sign_intramod_conec_Y[[1]])
nb_cor<-list()
nb_genes_A<-list()
nb_genes_Y<-list()
nb_genes_Ao<-list()
nb_genes_Yo<-list()
for (i in 1:length(top_gene_sign_intramod_conec_A)) {
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  unlist(top_gene_sign_intramod_conec_A[[i]]),])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  unlist(top_gene_sign_intramod_conec_Y[[i]]),])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
nb_cor[[i]]<-table(cor(NORM_EXP_TOPYUCA,NORM_EXP_TOPAMF)>0)
nb_genes_A[[i]]<-length(unlist(top_gene_sign_intramod_conec_A[[i]]))
  nb_genes_Y[[i]]<-length(unlist(top_gene_sign_intramod_conec_Y[[i]]))
  nb_genes_Yo[[i]]<- (length(cor(NORM_EXP_TOPYUCA))-dim(cor(NORM_EXP_TOPYUCA))) /2
  nb_genes_Ao[[i]]<- (length(cor(NORM_EXP_TOPAMF))-dim(cor(NORM_EXP_TOPAMF))) /2
}
nb_cor
unlist(nb_genes_Ao)
(length(cor(NORM_EXP_TOPYUCA))-dim(cor(NORM_EXP_TOPYUCA))) /2

################################################################################################################################
# plot of correlation between genes
######   ######  ######  ######  ###### ######
# manual correlation from two genes
t(sign_norm_genes_A)

NORM_EXP_TOPAMF<-t(sign_norm_genes_A[rownames(sign_norm_genes_A) %in%  "g2342.t1",])
NORM_EXP_TOPYUCA<-t(sign_norm_genes_Y[rownames(sign_norm_genes_Y) %in%  "Manes.18G142200",])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/remark_cor_genes_v1.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
plot(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA,pch=16,cex=3,ylab="Manes.09G039400 \n endopeptidase inhibitor activity",xlab="g2342.t1 peptidase S24/S26A-C",cex.lab=3,cex.axis=3)
abline(lm(NORM_EXP_TOPYUCA~as.vector(NORM_EXP_TOPAMF)),lwd=5)
dev.off()

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  "g2342.t1",])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.18G142200",])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/remark_cor_genes_v1.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
plot(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA,pch=16,cex=3,ylab="Manes.09G039400 \n endopeptidase inhibitor activity",xlab="g2342.t1 peptidase S24/S26A-C",cex.lab=3,cex.axis=3)
abline(lm(NORM_EXP_TOPYUCA~as.vector(NORM_EXP_TOPAMF)),lwd=5)
dev.off()

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  "g12728.t1",])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.09G039400",])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/remark_cor_genes_v1b.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
plot(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA,pch=16,cex=3,ylab="Manes.09G039400 \n endopeptidase inhibitor activity",xlab="g12728.t1 peptidase 20S proteasome beta subunit D1 (PBD1)",cex.lab=3,cex.axis=3)
abline(lm(NORM_EXP_TOPYUCA~as.vector(NORM_EXP_TOPAMF)),lwd=5)
dev.off()

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  "g1486.t1",])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.17G005400",])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/remark_cor_genes_v2.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
plot(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA,pch=16,cex=3,ylab="Manes.17G005400 \nSH3-domain protein (transport) ",xlab="g1486.t1 ALS3 (ABC-transporter)",cex.lab=3,cex.axis=3)
abline(lm(NORM_EXP_TOPYUCA~as.vector(NORM_EXP_TOPAMF)))
dev.off()




NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  "g12274.t1",])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.06G174400",])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/remark_cor_genes_v4.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
plot(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA,pch=16,cex=3,ylab="Manes.06G174400 \nKT1",xlab="g12274.t1 Glyoxal oxydase",cex.lab=3,cex.axis=3)
abline(lm(NORM_EXP_TOPYUCA~as.vector(NORM_EXP_TOPAMF)))
dev.off()

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  "g12274.t1",])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.12G031900",])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/remark_cor_genes_v3.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
plot(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA,pch=16,cex=3,ylab="Manes.12G031900 \nPT1",xlab="g12274.t1 Glyoxal oxydase",cex.lab=3,cex.axis=3)
abline(lm(NORM_EXP_TOPYUCA~as.vector(NORM_EXP_TOPAMF)))
dev.off()

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  "g12274.t1",])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.15G137300",])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/remark_cor_genes_v3b.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
plot(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA,pch=16,cex=3,ylab="Manes.15G137300 \nBarwin-related endoglucanase",xlab="g12274.t1 Glyoxal oxydase",cex.lab=3,cex.axis=3)
abline(lm(NORM_EXP_TOPYUCA~as.vector(NORM_EXP_TOPAMF)))
dev.off()


NORM_EXP_TOPYUCA2<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.12G031900",])
NORM_EXP_TOPYUCA2<-NORM_EXP_TOPYUCA2[,grep("CTRL",colnames(NORM_EXP_TOPYUCA2),invert=T)]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.06G174400",])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]

pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/remark_cor_genes_v5.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
options(scipen = 10)
plot(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA2,pch=16,cex=3,xlab="Manes.12G031900 PT1",ylab="Manes.06G174400 \nKT1",cex.lab=3,cex.axis=3)
abline(lm(as.vector(NORM_EXP_TOPYUCA2)~as.vector(NORM_EXP_TOPYUCA)))
dev.off()
#####  ######  ######  ######  ######  ######  ######
# all important correlations in a single plot correlogram

VIP_genes_Y<-c("Manes.12G031900","Manes.06G174400","Manes.15G137300",
               "Manes.09G039400","Manes.13G134300","Manes.09G142400","Manes.01G061100")
VIP_genes_A<-c("g12274.t1",
               "g2342.t1","g12728.t1","g1045.t1","g1360.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% VIP_genes_A,])

NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  VIP_genes_Y,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

correlacion<-cor(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)

colsc=c(rgb(241, 54, 23, maxColorValue=255),  "white" , rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space= 'Lab' )
colors = colramp(100)
labo<-data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/remark_cor_genes_vf.pdf',width=14,height=14,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
labo<-labo[,c(2,1,3:5,7,10,12,6,8:9,11)]
my.plotcorr(cex =2,cex.lab = 2,cor(labo), col=colors[((cor(labo) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )
dev.off()
#####  ######  ######  ######  ######  ######  ######
# function correlogram

my.plotcorr <- function (corr, outline = FALSE, col = "grey", upper.panel = c("ellipse", "number", "none"), lower.panel = c("ellipse", "number", "none"), diag = c("none", "ellipse", "number"), digits = 2, bty = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, cex.lab = par("cex.lab"), cex = 0.75 * par("cex"), mar = 0.1 + c(2, 2, 4, 2), ...)
{
  # this is a modified version of the plotcorr function from the ellipse package
  # this prints numbers and ellipses on the same plot but upper.panel and lower.panel changes what is displayed
  # diag now specifies what to put in the diagonal (numbers, ellipses, nothing)
  # digits specifies the number of digits after the . to round to
  # unlike the original, this function will always print x_i by x_i correlation rather than being able to drop it
  # modified by Esteban Buz
  if (!require('ellipse', quietly = TRUE, character = TRUE)) {
    stop("Need the ellipse library")
  }
  savepar <- par(pty = "s", mar = mar)
  on.exit(par(savepar))
  if (is.null(corr))
    return(invisible())
  if ((!is.matrix(corr)) || (round(min(corr, na.rm = TRUE), 6) < -1) || (round(max(corr, na.rm = TRUE), 6) > 1))
    stop("Need a correlation matrix")
  plot.new()
  par(new = TRUE)
  rowdim <- dim(corr)[1]
  coldim <- dim(corr)[2]
  rowlabs <- dimnames(corr)[[1]]
  collabs <- dimnames(corr)[[2]]
  if (is.null(rowlabs))
    rowlabs <- 1:rowdim
  if (is.null(collabs))
    collabs <- 1:coldim
  rowlabs <- as.character(rowlabs)
  collabs <- as.character(collabs)
  col <- rep(col, length = length(corr))
  dim(col) <- dim(corr)
  upper.panel <- match.arg(upper.panel)
  lower.panel <- match.arg(lower.panel)
  diag <- match.arg(diag)
  cols <- 1:coldim
  rows <- 1:rowdim
  maxdim <- max(length(rows), length(cols))
  plt <- par("plt")
  xlabwidth <- max(strwidth(rowlabs[rows], units = "figure", cex = cex.lab))/(plt[2] - plt[1])
  xlabwidth <- xlabwidth * maxdim/(1 - xlabwidth)
  ylabwidth <- max(strwidth(collabs[cols], units = "figure", cex = cex.lab))/(plt[4] - plt[3])
  ylabwidth <- ylabwidth * maxdim/(1 - ylabwidth)
  plot(c(-xlabwidth - 0.5, maxdim + 0.5), c(0.5, maxdim + 1 + ylabwidth), type = "n", bty = bty, axes = axes, xlab = "", ylab = "", asp = asp, cex.lab = cex.lab, ...)
  text(rep(0, length(rows)), length(rows):1, labels = rowlabs[rows], adj = 1, cex = cex.lab)
  text(cols, rep(length(rows) + 1, length(cols)), labels = collabs[cols], srt = 90, adj = 0, cex = cex.lab)
  mtext(xlab, 1, 0)
  mtext(ylab, 2, 0)
  mat <- diag(c(1, 1))
  plotcorrInternal <- function() {
    if (i == j){ #diag behavior
      if (diag == 'none'){
        return()
      } else if (diag == 'number'){
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else if (diag == 'ellipse') {
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      }
    } else if (i >= j){ #lower half of plot
      if (lower.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (lower.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    } else { #upper half of plot
      if (upper.panel == 'ellipse') { #check if ellipses should go here
        mat[1, 2] <- corr[i, j]
        mat[2, 1] <- mat[1, 2]
        ell <- ellipse(mat, t = 0.43)
        ell[, 1] <- ell[, 1] + j
        ell[, 2] <- ell[, 2] + length(rows) + 1 - i
        polygon(ell, col = col[i, j])
        if (outline)
          lines(ell)
      } else if (upper.panel == 'number') { #check if ellipses should go here
        text(j + 0.3, length(rows) + 1 - i, round(corr[i, j], digits=digits), adj = 1, cex = cex)
      } else {
        return()
      }
    }
  }
  for (i in 1:dim(corr)[1]) {
    for (j in 1:dim(corr)[2]) {
      plotcorrInternal()
    }
  }
  invisible()
}



names(top_gene_sign_intramod_conec_Y)[15]

names(top_gene_sign_intramod_conec_A)[15]

top_gene_sign_intramod_conec_A[[15]]

for (i in 1:length(names(top_gene_sign_intramod_conec_Y))) {
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% top_gene_sign_intramod_conec_A[[i]],])

NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  top_gene_sign_intramod_conec_Y[[i]],])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

correlacion<-cor(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)
colsc=c(rgb(241, 54, 23, maxColorValue=255),  "white" , rgb(0, 61, 104, maxColorValue=255))
colramp = colorRampPalette(colsc, space= 'Lab' )
colors = colramp(100)
labo<-data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)
pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/TOP_GENES/',names(top_gene_sign_intramod_conec_Y)[i],'cor_genes_vf.pdf',sep="_"),width=30,height=30,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
my.plotcorr(cex =2,cex.lab = 2,cor(labo), col=colors[((cor(labo) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )
dev.off()
}


###########################################################################################################################
# correlogram interesting genes in all modules


#cas<-unique(c("Manes.03G205000","Manes.11G078100","Manes.16G090600","Manes.01G208500","Manes.06G041300","Manes.07G034500","Manes.07G005600","Manes.16G024800","Manes.01G003900","Manes.09G041300","Manes.11G124400","Manes.13G071100","Manes.01G201800","Manes.03G044400","Manes.07G087500","Manes.09G120000","Manes.14G141400","Manes.18G060600","Manes.04G095800","Manes.09G182400","Manes.01G003900","Manes.16G020200","Manes.16G133800","Manes.03G194400","Manes.13G084400","Manes.16G133800","Manes.04G064300","Manes.10G116300","Manes.14G015500","Manes.17G053000","Manes.17G059000","Manes.01G205200","Manes.06G122900","Manes.09G129100","Manes.15G165000","Manes.17G059000","Manes.01G115400","Manes.01G265100","Manes.05G047400","Manes.06G122900","Manes.09G021000","Manes.10G121900","Manes.15G165000","Manes.17G059000","Manes.04G064300","Manes.09G129100","Manes.10G059900","Manes.16G110700","Manes.01G111700","Manes.01G139400","Manes.03G127800","Manes.03G171200","Manes.06G063500","Manes.06G155200","Manes.07G056100","Manes.14G100900","Manes.14G108900","Manes.06G148500","Manes.06G148600","Manes.07G094600","Manes.13G060300","Manes.16G052800","Manes.17G019200","Manes.08G141200","Manes.09G146200","Manes.09G156800","Manes.16G010100","Manes.17G018800"))
 #        amf<-unique(c("g1236.t1","g3417.t1","g8841.t1g11187.t1","g2022.t1","g8727.t1","g1108.t1","g6787.t1","g6816.t1","g981.t1","g10076.t1","g11821.t1","g30.t1","g5995.t1","g8763.t1g1004.t1","g12075.t1","g1651.t1","g3431.t1","g4807.t1","g4952.t1","g4976.t1","g5145.t1","g5801.t1","g7641.t1","g9899.t1g2998.t1","g7192.t1","g10703.t1","g2340.t1","g2594.t1","g5160.t1g11806.t1","g2305.t1","g3836.t1","g3987.t1","g4998.t1","g6035.t1","g646.t1","g6853.t1","g382.t1","g5472.t1","g687.t1","g7775.t1","g1045.t1","g11736.t1","g12727.t1","g1488.t1","g2417.t1","g2610.t1","g3371.t1","g4338.t1","g4376.t1","g5311.t1","g614.t1","g6404.t1","g65.t1","g6910.t1","g7358.t1","g9445.t1","g10703.t1","g1486.t1","g5160.t1","g5866.t1","g929.t1","g10995.t1","g1300.t1","g8590.t1","g8780.t1","g9390.t1","g7192.t1","g8161.t1","g8679.t1","g8727.t1","g14952.t1","g4959.t1","g7735.t1","g9089.t1","g10703.t1","g1207.t1","g2329.t1","g3510.t1","g5866.t1","g929.t1","g10995.t1","g3047.t1","g8590.t1"))

# RANDOM SELECTION OF 55 GENES per organism
         CasR<-unique(c("Manes.02G039800","Manes.02G059900","Manes.03G068200","Manes.03G068200","Manes.09G146200","Manes.09G146200","Manes.09G156800","Manes.09G156800","Manes.10G135100","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G126200","Manes.14G073700","Manes.17G018800","Manes.17G018800","Manes.02G039800","Manes.03G068200","Manes.03G068200","Manes.08G141200","Manes.09G146200","Manes.09G146200","Manes.09G156800","Manes.09G156800","Manes.13G126200","Manes.16G010100","Manes.17G018800","Manes.17G018800","Manes.02G039800","Manes.03G068200","Manes.03G068200","Manes.06G148500","Manes.06G148600","Manes.07G094600","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G126200","Manes.16G052800","Manes.17G019200","Manes.01G023200","Manes.01G111700","Manes.01G111700","Manes.01G123000","Manes.01G139400","Manes.02G069900","Manes.02G069900","Manes.02G069900","Manes.02G069900","Manes.03G127800","Manes.03G171200","Manes.03G202500","Manes.03G202500","Manes.05G047000","Manes.05G047000","Manes.06G063500","Manes.06G117800","Manes.06G155200","Manes.06G155200","Manes.06G164600","Manes.07G056100","Manes.09G059900","Manes.09G084500","Manes.09G087100","Manes.10G007300","Manes.11G004200","Manes.12G063500","Manes.12G063500","Manes.12G063500","Manes.14G007500","Manes.14G100900","Manes.14G108900","Manes.15G034400","Manes.15G094700","Manes.15G151100","Manes.15G151100","Manes.17G013000","Manes.01G023200"))
                       CasR2<- unique(c("Manes.01G111700","Manes.01G111700","Manes.01G139400","Manes.03G127800","Manes.03G171200","Manes.04G075000","Manes.04G075000","Manes.05G077100","Manes.05G077100","Manes.06G063500","Manes.06G117800","Manes.06G155200","Manes.06G155200","Manes.06G164600","Manes.07G056100","Manes.07G073600","Manes.09G059900","Manes.09G087100","Manes.10G007300","Manes.11G004200","Manes.12G035300","Manes.12G035300","Manes.12G063500","Manes.12G063500","Manes.12G063500","Manes.14G100900","Manes.14G108900","Manes.14G139600","Manes.15G094700","Manes.15G151100","Manes.15G151100","Manes.17G110900","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G115400","Manes.01G265100","Manes.03G206300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G141400","Manes.04G146400","Manes.05G047400","Manes.06G038200","Manes.06G122900","Manes.06G122900","Manes.08G049200","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G129100","Manes.09G129100","Manes.10G034400","Manes.10G034400","Manes.10G059900","Manes.11G114800","Manes.11G114800","Manes.12G013300","Manes.13G045600","Manes.14G015500","Manes.14G015500","Manes.15G165000","Manes.16G110700","Manes.17G059000","Manes.18G073000","Manes.01G115400","Manes.01G265100","Manes.02G119500","Manes.02G219500","Manes.03G206300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G146400","Manes.05G047400","Manes.06G038200","Manes.06G122900","Manes.06G122900","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.08G049200","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.10G121900","Manes.10G121900","Manes.14G015500","Manes.14G015500","Manes.15G165000","Manes.17G053000","Manes.17G059000","Manes.01G205200","Manes.01G205200","Manes.01G205200","Manes.04G146400","Manes.05G014200","Manes.05G136900","Manes.06G122900","Manes.06G122900","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.08G010900","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G129100","Manes.09G129100","Manes.10G059900","Manes.10G121900","Manes.10G121900","Manes.12G013300","Manes.15G165000","Manes.17G059000","Manes.18G031200","Manes.18G073000","Manes.18G094200","Manes.02G119500","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G090100","Manes.04G146400","Manes.04G161100","Manes.04G161100","Manes.04G161100","Manes.04G161100","Manes.10G116300","Manes.13G070400","Manes.14G015500","Manes.14G015500","Manes.17G053000","Manes.17G059000","Manes.01G208500","Manes.02G135300","Manes.06G041300","Manes.07G034500","Manes.14G122800","Manes.15G022900","Manes.01G201800","Manes.03G044400","Manes.07G087500","Manes.11G078100","Manes.13G008000","Manes.16G090600","Manes.16G090600","Manes.16G090600","Manes.03G205000","Manes.03G205000","Manes.04G095800","Manes.11G078100","Manes.14G061800","Manes.15G087400","Manes.16G090600","Manes.16G090600","Manes.16G090600","Manes.03G205000","Manes.03G205000","Manes.09G120000","Manes.09G120000","Manes.09G182400","Manes.14G141400","Manes.18G060600","Manes.18G060600","Manes.04G095800","Manes.09G182400","Manes.07G005600","Manes.07G092000","Manes.16G024800","Manes.01G003900","Manes.11G154700","Manes.16G020200","Manes.16G133800","Manes.16G133800","Manes.03G194400","Manes.11G154700","Manes.13G084400","Manes.16G133800","Manes.16G133800","Manes.01G003900","Manes.09G041300","Manes.11G124400","Manes.11G154700","Manes.13G071100"))
         CasR3<-unique(c(CasR,CasR2))
                cas<-CasR3      
         AMF_R<-unique(c("g1008.t1","g11187.t1","g11634.t1","g2998.t1","g3197.t1","g339.t1","g3641.t1","g4193.t1","g4624.t1","g5528.t1","g7192.t1","g7415.t1","g8161.t1","g8671.t1","g8671.t1","g8680.t1","g10995.t1","g1988.t1","g2737.t1","g3047.t1","g568.t1","g8590.t1","g105.t1","g10703.t1","g1207.t1","g1757.t1","g21.t1","g2329.t1","g3510.t1","g3844.t1","g5111.t1","g5456.t1","g5866.t1","g652.t1","g8553.t1","g8944.t1","g929.t1","g9749.t1","g10683.t1","g10723.t1","g10723.t1","g10723.t1","g11991.t1","g12075.t1","g12890.t1","g1385.t1","g1651.t1","g2394.t1","g2427.t1","g2715.t1","g2879.t1","g3065.t1","g3508.t1","g4871.t1","g4952.t1","g5145.t1","g5801.t1","g5821.t1","g5979.t1","g7641.t1","g8429.t1","g9899.t1","g14952.t1","g14952.t1","g4959.t1","g4978.t1","g5557.t1","g7735.t1","g9089.t1","g1008.t1","g11187.t1","g1565.t1","g1659.t1","g339.t1","g4193.t1","g7192.t1","g8161.t1","g8679.t1","g8727.t1","g8738.t1","g8981.t1","g10995.t1","g1300.t1","g1988.t1","g6205.t1","g8590.t1","g8780.t1","g9390.t1","g105.t1","g10703.t1","g1072.t1","g1486.t1","g1486.t1","g1486.t1","g177.t1","g2594.t1","g4372.t1","g5160.t1","g5456.t1","g5866.t1","g7772.t1","g8553.t1","g8944.t1","g9013.t1","g929.t1","g1012.t1","g10122.t1","g1022.t1","g1045.t1","g11736.t1","g12309.t1","g12727.t1","g12728.t1","g1339.t1","g1488.t1","g1507.t1","g1999.t1","g205.t1","g208.t1","g2331.t1","g2417.t1","g2421.t1","g2610.t1","g2610.t1","g2657.t1","g2662.t1","g2808.t1","g2910.t1","g3270.t1","g3327.t1","g3371.t1","g3377.t1","g3830.t1","g3904.t1","g4338.t1","g4338.t1","g4338.t1","g4376.t1","g441.t1","g4851.t1","g5093.t1","g526.t1","g5288.t1","g5311.t1","g548.t1","g5691.t1","g5697.t1","g614.t1","g623.t1","g629.t1","g638.t1","g6404.t1","g6408.t1","g6408.t1","g642.t1","g65.t1","g6748.t1","g6910.t1","g7241.t1","g7356.t1","g7358.t1","g7358.t1","g7654.t1","g8517.t1","g8613.t1","g8748.t1","g881.t1","g9174.t1","g9445.t1","g973.t1","g9739.t1","g9840.t1","g992.t1","g1008.t1","g11187.t1","g1565.t1","g2022.t1","g8679.t1","g8727.t1","g1004.t1","g1004.t1","g11991.t1","g12075.t1","g1651.t1","g2394.t1","g2879.t1","g3431.t1","g3431.t1","g3508.t1","g3613.t1","g4807.t1","g4871.t1","g4952.t1","g4976.t1","g5145.t1","g5801.t1","g5979.t1","g7641.t1","g9899.t1","g1236.t1","g1236.t1","g3417.t1","g407.t1","g4188.t1","g5983.t1","g8841.t1","g11187.t1","g11634.t1","g1565.t1","g2998.t1","g339.t1","g4193.t1","g7192.t1","g7415.t1","g8738.t1","g8952.t1","g8981.t1","g10703.t1","g1072.t1","g13844.t1","g177.t1","g2329.t1","g2340.t1","g2594.t1","g5160.t1","g5395.t1","g5456.t1","g7698.t1","g8944.t1","g9013.t1","g9388.t1","g1108.t1","g11163.t1","g11668.t1","g2567.t1","g3041.t1","g3094.t1","g3228.t1","g4191.t1","g4191.t1","g4191.t1","g5338.t1","g5338.t1","g6538.t1","g6787.t1","g6816.t1","g7673.t1","g8455.t1","g928.t1","g981.t1","g9834.t1","g10019.t1","g10098.t1","g10803.t1","g1165.t1","g11806.t1","g12383.t1","g1242.t1","g13288.t1","g1570.t1","g2056.t1","g2303.t1","g2303.t1","g2305.t1","g2345.t1","g2349.t1","g2350.t1","g2968.t1","g3399.t1","g372.t1","g3836.t1","g3987.t1","g4051.t1","g4435.t1","g4435.t1","g4998.t1","g5163.t1","g5253.t1","g5665.t1","g5802.t1","g6035.t1","g6165.t1","g6165.t1","g646.t1","g6853.t1","g7670.t1","g7867.t1","g8801.t1","g10271.t1","g11031.t1","g382.t1","g5472.t1","g6255.t1","g6598.t1","g687.t1","g7716.t1","g7775.t1","g8837.t1","g10076.t1","g10292.t1","g11821.t1","g12154.t1","g12914.t1","g13967.t1","g1518.t1","g2043.t1","g2337.t1","g2510.t1","g2748.t1","g284.t1","g296.t1","g30.t1","g3468.t1","g410.t1","g4281.t1","g4767.t1","g4790.t1","g5164.t1","g5223.t1","g5465.t1","g5995.t1","g6405.t1","g6750.t1","g8058.t1","g8061.t1","g8188.t1","g8763.t1","g9477.t1","g9871.t1","g9894.t1"))
         amf<-sample(AMF_R,120)

         NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])

NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/ALL_genes_random55choice.pdf',width=15,height=15,useDingbats = FALSE)
pheatmap(cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) )
dev.off()
#####################
# Select some genes for visualisation Secretory pathway AMF
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1")
                               ,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas
                                   ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot<-cor2plot[grep("Manes",colnames(cor2plot)),]
cor2plot<-cor2plot[,grep(".t1",colnames(cor2plot))]
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/AMF_secretory_plant.pdf',width=15,height=15,useDingbats = FALSE)
pheatmap(t(cor2plot) ,cellwidth = 8,cellheight = 8)
dev.off()
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas
                                 ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]

NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  c("Manes.01G023200","Manes.06G155200","Manes.15G034400","Manes.18G060600","Manes.01G139400","Manes.04G095800","Manes.10G007300")
                                 ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot<-cor2plot[colnames(cor2plot) %in%
                     c("Manes.01G023200","Manes.06G155200","Manes.15G034400","Manes.18G060600","Manes.01G139400","Manes.04G095800","Manes.10G007300") ,]
pheatmap(cor2plot)
# in AMF
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf
                               ,])

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1")
                               ,])
cor2plot<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot<-cor2plot[colnames(cor2plot) %in%
                     c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1") ,]
cor2plot<-cor2plot[,!colnames(cor2plot) %in%
                     c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1") ]
cor2plot<-cor2plot[,grep("\\.1",colnames(cor2plot),invert=T)]
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/AMF_secretory_amf.pdf',width=15,height=15,useDingbats = FALSE)
pheatmap(cor2plot ,cellwidth = 8,cellheight = 8)
dev.off()

pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/colors_cretory2.pdf',width=15,height=15,useDingbats = FALSE)
colfinal<-c("red","red","red","red","red","red","red","red","red","red","green","green","red","red","red","red","red","red","red","green","red","green","green","green","green","green","green","green","grey60","grey60","grey60","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","black","magenta","lightcyan","magenta","lightcyan","lightyellow","lightcyan","lightcyan","lightcyan","lightcyan","lightcyan","lightcyan","lightcyan","black","black","turquoise","turquoise","turquoise","lightyellow","lightyellow","turquoise","turquoise","turquoise","turquoise","black","black","black","black","black","black","black","black","black","black","black","royalblue","royalblue","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","magenta","lightcyan","cyan","cyan","cyan","red","green","green","green","grey60","red","grey60","green")
colsec<-c("grey60","red","turquoise","turquoise","turquoise","lightyellow","black","black","lightyellow","lightcyan")
colCASSAVA<-c("brown","brown","brown","brown","brown","brown","green","green","brown","black","brown","brown","brown","brown","brown","brown","brown","brown","brown","brown","brown","magenta","brown","brown","yellow","green","green","brown","brown","blue","brown","brown","yellow","yellow","yellow","yellow","yellow","yellow","brown","yellow","yellow","pink","pink","pink","magenta","magenta","magenta","magenta","magenta","blue","magenta","blue","magenta","magenta","magenta","magenta","magenta","black","brown","black","black","green","black","brown","black","black","black","brown","black","brown","brown","black","black","black","brown","black","black","green","brown","black","black","brown","black","brown","brown","brown","brown","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","magenta","magenta","blue","blue","blue","blue","blue","blue","blue","blue","blue","magenta","blue")
barplot(rep(5,10),col=colsec,border=NA)
#names
rev(c("Prot. Degr.","Vac. AA transp. 5","Iron-sulfur prot.","Heat-shock 12A","Methyltransferase","NA","F-box/WD 4","glucosyl transf.","KxDL ","Prot. Degr.","Yippee-like","NA","ABC transp.","NudC 2","Ubiquitin","NA","Methyltransferase","Sugar transport","uricase-2","RNA pol II transc.","Prot. Degr.","Endonuclease","NA","Secretory pathway","NA","NA","Intra-golgi transp.","NA","RNA binding","RNA binding","Cell organisation","RNA binding","Prot. Degr.","Signalling calcium","NA","NA","NA","Krebs cycle","Transc. Factor","postranslational mod.","Proton transport","Lipid metabo.","Glycolysis","Sucrose synthase","ATPase transporter","Cell organisation","Prot. Synthesis","Prot. Degr.","Secretory pathway","Hormone metabo.","Signalling","Secretory pathway","Secretory pathway","NA","DNAJ heat shock","Cell wall organisation","Signalling","Prot. Degr.","NA","Secretory pathway","Secretory pathway","Org. Dev.","NA","Cell organisation","Cell organisation","Prot. Synthesis","Prot. Degr.","Aldolase","Lipid metabo.","NA","Aldolase","Lipid Degr.","Alpha/beta hydrolase","Secretory pathway","Glycolysis","Glycolysis","Lipid metabo.","NA","Sugar transport","ATP synthesis","Transport at membrane","Sugar transport","Glycolysis","Prot. Degr.","Nucleotide metabo.","Peptide amidase","Alpha/beta hydrolase","Cytocrhome P450","Signalling calcium","NA","Cell organisation","Nucleoporin","ABC transp.","Proteine kinase","Cell organisation","RNA binding","NA","Phosphatase","Serine protease","NA","Secretory pathway","NA","postranslational mod.","NA","NA","NA","postranslational mod.","postranslational mod.","NA","Serine protease")
 
    dev.off()
    
# similar genes to secretory pathway.
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% unique(as.vector(unlist(top_gene_sign_intramod_conec_A)))
,])

cor2plot<-cor(data.frame(NORM_EXP_TOPAMF))
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/ALL_AMF.pdf',width=15,height=15,useDingbats = FALSE)

pheatmap(cor2plot,fontsize = 5)
dev.off()

####
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/colors_random55choice.pdf',width=15,height=15,useDingbats = FALSE)
color_fig3<-rev(c("yellow","green","yellow","grey60","grey60","pink","pink","pink","green","brown","brown","blue","green","brown","brown","brown","brown","brown","brown","brown","brown","magenta","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","turquoise","cyan","cyan","blue","cyan","cyan","magenta","magenta","magenta","magenta","black","brown","black","black","black","black","green","lightyellow","black","black","black","blue","brown","black","black","lightcyan","black","red","magenta","magenta","magenta","magenta","black","magenta","blue","blue","blue","black","black","black","blue","blue","blue","blue","blue","magenta","blue","blue","blue","blue","blue","blue","blue","royalblue","turquoise","turquoise","turquoise","turquoise","turquoise","green","green","turquoise","red","red","green","green","grey60","red","red","red","red","red","red","red","red","red"))
barplot(rep(5,110),col=color_fig3,border=NA)
#names
rev(c("Prot. Degr.","Vac. AA transp. 5","Iron-sulfur prot.","Heat-shock 12A","Methyltransferase","NA","F-box/WD 4","glucosyl transf.","KxDL ","Prot. Degr.","Yippee-like","NA","ABC transp.","NudC 2","Ubiquitin","NA","Methyltransferase","Sugar transport","uricase-2","RNA pol II transc.","Prot. Degr.","Endonuclease","NA","Secretory pathway","NA","NA","Intra-golgi transp.","NA","RNA binding","RNA binding","Cell organisation","RNA binding","Prot. Degr.","Signalling calcium","NA","NA","NA","Krebs cycle","Transc. Factor","postranslational mod.","Proton transport","Lipid metabo.","Glycolysis","Sucrose synthase","ATPase transporter","Cell organisation","Prot. Synthesis","Prot. Degr.","Secretory pathway","Hormone metabo.","Signalling","Secretory pathway","Secretory pathway","NA","DNAJ heat shock","Cell wall organisation","Signalling","Prot. Degr.","NA","Secretory pathway","Secretory pathway","Org. Dev.","NA","Cell organisation","Cell organisation","Prot. Synthesis","Prot. Degr.","Aldolase","Lipid metabo.","NA","Aldolase","Lipid Degr.","Alpha/beta hydrolase","Secretory pathway","Glycolysis","Glycolysis","Lipid metabo.","NA","Sugar transport","ATP synthesis","Transport at membrane","Sugar transport","Glycolysis","Prot. Degr.","Nucleotide metabo.","Peptide amidase","Alpha/beta hydrolase","Cytocrhome P450","Signalling calcium","NA","Cell organisation","Nucleoporin","ABC transp.","Proteine kinase","Cell organisation","RNA binding","NA","Phosphatase","Serine protease","NA","Secretory pathway","NA","postranslational mod.","NA","NA","NA","postranslational mod.","postranslational mod.","NA","Serine protease")
dev.off()


pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/correlo_interest.pdf',width=15,height=15,useDingbats = FALSE)
my.plotcorr(cex =1,cex.lab = 1,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )
dev.off()
rownames(cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)))
####
##############################################################################
# Visualization of functional genes correlation against AMF and PLANT

# all top genes
CasR<-unique(c("Manes.02G039800","Manes.02G059900","Manes.03G068200","Manes.03G068200","Manes.09G146200","Manes.09G146200","Manes.09G156800","Manes.09G156800","Manes.10G135100","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G126200","Manes.14G073700","Manes.17G018800","Manes.17G018800","Manes.02G039800","Manes.03G068200","Manes.03G068200","Manes.08G141200","Manes.09G146200","Manes.09G146200","Manes.09G156800","Manes.09G156800","Manes.13G126200","Manes.16G010100","Manes.17G018800","Manes.17G018800","Manes.02G039800","Manes.03G068200","Manes.03G068200","Manes.06G148500","Manes.06G148600","Manes.07G094600","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G126200","Manes.16G052800","Manes.17G019200","Manes.01G023200","Manes.01G111700","Manes.01G111700","Manes.01G123000","Manes.01G139400","Manes.02G069900","Manes.02G069900","Manes.02G069900","Manes.02G069900","Manes.03G127800","Manes.03G171200","Manes.03G202500","Manes.03G202500","Manes.05G047000","Manes.05G047000","Manes.06G063500","Manes.06G117800","Manes.06G155200","Manes.06G155200","Manes.06G164600","Manes.07G056100","Manes.09G059900","Manes.09G084500","Manes.09G087100","Manes.10G007300","Manes.11G004200","Manes.12G063500","Manes.12G063500","Manes.12G063500","Manes.14G007500","Manes.14G100900","Manes.14G108900","Manes.15G034400","Manes.15G094700","Manes.15G151100","Manes.15G151100","Manes.17G013000","Manes.01G023200"))
CasR2<- unique(c("Manes.01G111700","Manes.01G111700","Manes.01G139400","Manes.03G127800","Manes.03G171200","Manes.04G075000","Manes.04G075000","Manes.05G077100","Manes.05G077100","Manes.06G063500","Manes.06G117800","Manes.06G155200","Manes.06G155200","Manes.06G164600","Manes.07G056100","Manes.07G073600","Manes.09G059900","Manes.09G087100","Manes.10G007300","Manes.11G004200","Manes.12G035300","Manes.12G035300","Manes.12G063500","Manes.12G063500","Manes.12G063500","Manes.14G100900","Manes.14G108900","Manes.14G139600","Manes.15G094700","Manes.15G151100","Manes.15G151100","Manes.17G110900","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G115400","Manes.01G265100","Manes.03G206300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G141400","Manes.04G146400","Manes.05G047400","Manes.06G038200","Manes.06G122900","Manes.06G122900","Manes.08G049200","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G129100","Manes.09G129100","Manes.10G034400","Manes.10G034400","Manes.10G059900","Manes.11G114800","Manes.11G114800","Manes.12G013300","Manes.13G045600","Manes.14G015500","Manes.14G015500","Manes.15G165000","Manes.16G110700","Manes.17G059000","Manes.18G073000","Manes.01G115400","Manes.01G265100","Manes.02G119500","Manes.02G219500","Manes.03G206300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G146400","Manes.05G047400","Manes.06G038200","Manes.06G122900","Manes.06G122900","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.08G049200","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.10G121900","Manes.10G121900","Manes.14G015500","Manes.14G015500","Manes.15G165000","Manes.17G053000","Manes.17G059000","Manes.01G205200","Manes.01G205200","Manes.01G205200","Manes.04G146400","Manes.05G014200","Manes.05G136900","Manes.06G122900","Manes.06G122900","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.08G010900","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G129100","Manes.09G129100","Manes.10G059900","Manes.10G121900","Manes.10G121900","Manes.12G013300","Manes.15G165000","Manes.17G059000","Manes.18G031200","Manes.18G073000","Manes.18G094200","Manes.02G119500","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G090100","Manes.04G146400","Manes.04G161100","Manes.04G161100","Manes.04G161100","Manes.04G161100","Manes.10G116300","Manes.13G070400","Manes.14G015500","Manes.14G015500","Manes.17G053000","Manes.17G059000","Manes.01G208500","Manes.02G135300","Manes.06G041300","Manes.07G034500","Manes.14G122800","Manes.15G022900","Manes.01G201800","Manes.03G044400","Manes.07G087500","Manes.11G078100","Manes.13G008000","Manes.16G090600","Manes.16G090600","Manes.16G090600","Manes.03G205000","Manes.03G205000","Manes.04G095800","Manes.11G078100","Manes.14G061800","Manes.15G087400","Manes.16G090600","Manes.16G090600","Manes.16G090600","Manes.03G205000","Manes.03G205000","Manes.09G120000","Manes.09G120000","Manes.09G182400","Manes.14G141400","Manes.18G060600","Manes.18G060600","Manes.04G095800","Manes.09G182400","Manes.07G005600","Manes.07G092000","Manes.16G024800","Manes.01G003900","Manes.11G154700","Manes.16G020200","Manes.16G133800","Manes.16G133800","Manes.03G194400","Manes.11G154700","Manes.13G084400","Manes.16G133800","Manes.16G133800","Manes.01G003900","Manes.09G041300","Manes.11G124400","Manes.11G154700","Manes.13G071100"))
CasR3<-unique(c(CasR,CasR2))
cas<-CasR3      
amf<-unique(c("g1008.t1","g11187.t1","g11634.t1","g2998.t1","g3197.t1","g339.t1","g3641.t1","g4193.t1","g4624.t1","g5528.t1","g7192.t1","g7415.t1","g8161.t1","g8671.t1","g8671.t1","g8680.t1","g10995.t1","g1988.t1","g2737.t1","g3047.t1","g568.t1","g8590.t1","g105.t1","g10703.t1","g1207.t1","g1757.t1","g21.t1","g2329.t1","g3510.t1","g3844.t1","g5111.t1","g5456.t1","g5866.t1","g652.t1","g8553.t1","g8944.t1","g929.t1","g9749.t1","g10683.t1","g10723.t1","g10723.t1","g10723.t1","g11991.t1","g12075.t1","g12890.t1","g1385.t1","g1651.t1","g2394.t1","g2427.t1","g2715.t1","g2879.t1","g3065.t1","g3508.t1","g4871.t1","g4952.t1","g5145.t1","g5801.t1","g5821.t1","g5979.t1","g7641.t1","g8429.t1","g9899.t1","g14952.t1","g14952.t1","g4959.t1","g4978.t1","g5557.t1","g7735.t1","g9089.t1","g1008.t1","g11187.t1","g1565.t1","g1659.t1","g339.t1","g4193.t1","g7192.t1","g8161.t1","g8679.t1","g8727.t1","g8738.t1","g8981.t1","g10995.t1","g1300.t1","g1988.t1","g6205.t1","g8590.t1","g8780.t1","g9390.t1","g105.t1","g10703.t1","g1072.t1","g1486.t1","g1486.t1","g1486.t1","g177.t1","g2594.t1","g4372.t1","g5160.t1","g5456.t1","g5866.t1","g7772.t1","g8553.t1","g8944.t1","g9013.t1","g929.t1","g1012.t1","g10122.t1","g1022.t1","g1045.t1","g11736.t1","g12309.t1","g12727.t1","g12728.t1","g1339.t1","g1488.t1","g1507.t1","g1999.t1","g205.t1","g208.t1","g2331.t1","g2417.t1","g2421.t1","g2610.t1","g2610.t1","g2657.t1","g2662.t1","g2808.t1","g2910.t1","g3270.t1","g3327.t1","g3371.t1","g3377.t1","g3830.t1","g3904.t1","g4338.t1","g4338.t1","g4338.t1","g4376.t1","g441.t1","g4851.t1","g5093.t1","g526.t1","g5288.t1","g5311.t1","g548.t1","g5691.t1","g5697.t1","g614.t1","g623.t1","g629.t1","g638.t1","g6404.t1","g6408.t1","g6408.t1","g642.t1","g65.t1","g6748.t1","g6910.t1","g7241.t1","g7356.t1","g7358.t1","g7358.t1","g7654.t1","g8517.t1","g8613.t1","g8748.t1","g881.t1","g9174.t1","g9445.t1","g973.t1","g9739.t1","g9840.t1","g992.t1","g1008.t1","g11187.t1","g1565.t1","g2022.t1","g8679.t1","g8727.t1","g1004.t1","g1004.t1","g11991.t1","g12075.t1","g1651.t1","g2394.t1","g2879.t1","g3431.t1","g3431.t1","g3508.t1","g3613.t1","g4807.t1","g4871.t1","g4952.t1","g4976.t1","g5145.t1","g5801.t1","g5979.t1","g7641.t1","g9899.t1","g1236.t1","g1236.t1","g3417.t1","g407.t1","g4188.t1","g5983.t1","g8841.t1","g11187.t1","g11634.t1","g1565.t1","g2998.t1","g339.t1","g4193.t1","g7192.t1","g7415.t1","g8738.t1","g8952.t1","g8981.t1","g10703.t1","g1072.t1","g13844.t1","g177.t1","g2329.t1","g2340.t1","g2594.t1","g5160.t1","g5395.t1","g5456.t1","g7698.t1","g8944.t1","g9013.t1","g9388.t1","g1108.t1","g11163.t1","g11668.t1","g2567.t1","g3041.t1","g3094.t1","g3228.t1","g4191.t1","g4191.t1","g4191.t1","g5338.t1","g5338.t1","g6538.t1","g6787.t1","g6816.t1","g7673.t1","g8455.t1","g928.t1","g981.t1","g9834.t1","g10019.t1","g10098.t1","g10803.t1","g1165.t1","g11806.t1","g12383.t1","g1242.t1","g13288.t1","g1570.t1","g2056.t1","g2303.t1","g2303.t1","g2305.t1","g2345.t1","g2349.t1","g2350.t1","g2968.t1","g3399.t1","g372.t1","g3836.t1","g3987.t1","g4051.t1","g4435.t1","g4435.t1","g4998.t1","g5163.t1","g5253.t1","g5665.t1","g5802.t1","g6035.t1","g6165.t1","g6165.t1","g646.t1","g6853.t1","g7670.t1","g7867.t1","g8801.t1","g10271.t1","g11031.t1","g382.t1","g5472.t1","g6255.t1","g6598.t1","g687.t1","g7716.t1","g7775.t1","g8837.t1","g10076.t1","g10292.t1","g11821.t1","g12154.t1","g12914.t1","g13967.t1","g1518.t1","g2043.t1","g2337.t1","g2510.t1","g2748.t1","g284.t1","g296.t1","g30.t1","g3468.t1","g410.t1","g4281.t1","g4767.t1","g4790.t1","g5164.t1","g5223.t1","g5465.t1","g5995.t1","g6405.t1","g6750.t1","g8058.t1","g8061.t1","g8188.t1","g8763.t1","g9477.t1","g9871.t1","g9894.t1"))


# Select some genes for visualisation Secretory pathway AMF
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1")
                               ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas
                                 ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot<-cor2plot[grep("Manes",colnames(cor2plot)),]
cor2plot<-cor2plot[,grep(".t1",colnames(cor2plot))]
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/AMF_secretory_plant.pdf',width=15,height=15,useDingbats = FALSE)
pheatmap(t(cor2plot) ,cellwidth = 8,cellheight = 8)
dev.off()

# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1")
                               ,])
cor2plot<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot<-cor2plot[colnames(cor2plot) %in%
                     c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1") ,]
cor2plot<-cor2plot[,!colnames(cor2plot) %in%
                     c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1") ]
cor2plot<-cor2plot[,grep("\\.1",colnames(cor2plot),invert=T)]
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/AMF_secretory_amf.pdf',width=15,height=15,useDingbats = FALSE)
pheatmap(cor2plot ,cellwidth = 8,cellheight = 8)
dev.off()












##############################################################################
############# PAIRWISE GENE COMPARISONS
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopA.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.03G205000","Manes.11G078100","Manes.16G090600")
AMF_mod_toptop<-c("g1236.t1","g3417.t1","g8841.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopB.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.01G208500","Manes.06G041300","Manes.07G034500")
AMF_mod_toptop<-c("g11187.t1","g2022.t1","g8727.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )
dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopC.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.07G005600","Manes.16G024800")
AMF_mod_toptop<-c("g1108.t1","g11163.t1","g6787.t1","g6816.t1","g981.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopD.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.01G003900","Manes.09G041300","Manes.11G124400","Manes.11G154700","Manes.13G071100")
AMF_mod_toptop<-c("g10076.t1","g11821.t1","g30.t1","g5995.t1","g8763.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopE.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.01G201800","Manes.03G044400","Manes.07G087500","Manes.11G078100","Manes.16G090600")
AMF_mod_toptop<-c("g1004.t1","g11991.t1","g12075.t1","g1651.t1","g3431.t1","g4807.t1","g4952.t1","g4976.t1","g5145.t1","g5801.t1","g7641.t1","g9899.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopF.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.09G120000","Manes.14G141400","Manes.18G060600")
AMF_mod_toptop<-c("g2998.t1","g7192.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopG.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.04G095800","Manes.09G182400")
AMF_mod_toptop<-c("g10703.t1","g2329.t1","g2340.t1","g2594.t1","g5160.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopH.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.01G003900","Manes.11G154700","Manes.16G020200","Manes.16G133800")
AMF_mod_toptop<-c("g11806.t1","g12383.t1","g2305.t1","g3836.t1","g3987.t1","g4998.t1","g6035.t1","g646.t1","g6853.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )



dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopI.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.03G194400","Manes.13G084400","Manes.16G133800")
AMF_mod_toptop<-c("g10271.t1","g382.t1","g5472.t1","g687.t1","g7775.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopJ.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.04G064300","Manes.10G116300","Manes.14G015500","Manes.17G053000","Manes.17G059000")
AMF_mod_toptop<-c("g10122.t1","g1045.t1","g11736.t1","g12309.t1","g12727.t1","g1339.t1","g1488.t1","g205.t1","g2417.t1","g2610.t1","g3371.t1","g4338.t1","g4376.t1","g4851.t1","g5311.t1","g614.t1","g6404.t1","g65.t1","g6910.t1","g7358.t1","g9445.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopK.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.01G205200","Manes.06G122900","Manes.09G129100","Manes.10G121900","Manes.15G165000","Manes.17G059000")
AMF_mod_toptop<-c("g10703.t1","g1486.t1","g2594.t1","g5160.t1","g5866.t1","g929.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopL.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.01G115400","Manes.01G265100","Manes.03G206300","Manes.04G064300","Manes.05G047400","Manes.06G122900","Manes.09G021000","Manes.10G121900","Manes.15G165000","Manes.17G059000")
AMF_mod_toptop<-c("g10995.t1","g1300.t1","g8590.t1","g8780.t1","g9390.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopM.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.01G115400","Manes.04G064300","Manes.05G047400","Manes.06G122900","Manes.09G021000","Manes.09G129100","Manes.10G059900","Manes.16G110700")
AMF_mod_toptop<-c("g1565.t1","g4193.t1","g7192.t1","g8161.t1","g8679.t1","g8727.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopN.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.01G023200","Manes.01G111700","Manes.01G139400","Manes.03G127800","Manes.03G171200","Manes.06G063500","Manes.06G155200","Manes.07G056100","Manes.10G007300","Manes.14G100900","Manes.14G108900")
AMF_mod_toptop<-c("g14952.t1","g4959.t1","g7735.t1","g9089.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopO.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.06G148500","Manes.06G148600","Manes.07G094600","Manes.13G060300","Manes.16G052800","Manes.17G019200")
AMF_mod_toptop<-c("g10703.t1","g1207.t1","g2329.t1","g3510.t1","g3844.t1","g5866.t1","g652.t1","g929.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )

dev.off() 
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/fig3toptopP.pdf',width=8,height=8,useDingbats = FALSE)
cassava_mod_toptop<-c("Manes.03G068200","Manes.08G141200","Manes.09G146200","Manes.09G156800","Manes.16G010100","Manes.17G018800")
AMF_mod_toptop<-c("g10995.t1","g3047.t1","g8590.t1")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% AMF_mod_toptop,])
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cassava_mod_toptop,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
my.plotcorr(cex=1.5,cex.lab = 2,cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)), col=colors[((cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA)) + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0) )
dev.off() 
###########################################################################################################################
#### circos plot
library('migest')
library("plyr")
library("circlize")

df1<-cbind.data.frame(
c(unique(sapply(strsplit(gsub("Y_ME","Cassava.",gsub("A_ME","AMF.",names(top_gene_sign_intramod_conec_A))),"_"), "[[", 1)),
        unique(sapply(strsplit(gsub("Y_ME","Cassava.",gsub("A_ME","AMF.",names(top_gene_sign_intramod_conec_A))),"_"), "[[", 3))),
  c("248,129,153","232,120,141","217,111,130","202,102,119","186,93,107","171,84,96","156,75,85","140,66,73","125,57,62","110,48,51",
    "208,251,233","185,228,206","162,205,179","139,182,152","116,159,125","93,136,98","71,113,71"),
c( unique(sapply(strsplit(gsub("Y_ME","Cassava.",gsub("A_ME","AMF.",names(top_gene_sign_intramod_conec_A))),"_"), "[[", 1)),
  unique(sapply(strsplit(gsub("Y_ME","Cassava.",gsub("A_ME","AMF.",names(top_gene_sign_intramod_conec_A))),"_"), "[[", 3)))
)
colnames(df1)<-c("order","rgb","region")

df1$rgb<-as.character(df1$rgb)
df2a<-cbind.data.frame(
sapply(strsplit(gsub("Y_ME","Cassava.",gsub("A_ME","AMF.",names(top_gene_sign_intramod_conec_A))),"_"), "[[", 1),
sapply(strsplit(gsub("Y_ME","Cassava.",gsub("A_ME","AMF.",names(top_gene_sign_intramod_conec_A))),"_"), "[[", 3),
unlist(lapply(top_gene_sign_intramod_conec_A,function (x) length(x)))         )
colnames(df2a)<-c("orig","dest","m")
df2b<-cbind.data.frame(
sapply(strsplit(gsub("Y_ME","Cassava.",gsub("A_ME","AMF.",names(top_gene_sign_intramod_conec_Y))),"_"), "[[", 1),
sapply(strsplit(gsub("Y_ME","Cassava.",gsub("A_ME","AMF.",names(top_gene_sign_intramod_conec_Y))),"_"), "[[", 3),
unlist(lapply(top_gene_sign_intramod_conec_Y,function (x) length(x)))            )
colnames(df2b)<-c("orig","dest","m")
df1
df2<-rbind.data.frame(df2a,df2b)

top_gene_sign_intramod_conec_Y
sapply(strsplit(correlated_Modula," "), "[[", 1)


df1$xmin<-0
df1$xmax<-255
n<-nrow(df1)
df1<-cbind(df1, matrix(as.numeric(unlist(strsplit(df1$rgb,","))),nrow=n,byrow=TRUE))
names(df1)[ncol(df1)-2:0]<-c("r","g","b")
df1$rcol<-rgb(df1$r,df1$g,df1$b,max=255)
df1$rcol<-rgb(df1$r,df1$g,df1$b,alpha=200,max=255)

par(mar=rep(0,4))
circos.clear()
circos.par(cell.padding=c(0,0,0,0),
           track.margin=c(0,0.15),start.degree=90,gap.degree=4)
circos.initialize(factors=df1$region,
                  xlim=cbind(df1$xmin, df1$xmax))
circos.trackPlotRegion(ylim = c(0, 1),
                       factors = df1$region, track.height=0.1,
                       panel.fun = function(x, y) {
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270,
                                      "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa =c(0, 0.5)
                         circos.text(x=mean(xlim), y=1.7,
                                     labels=name, facing = dd, cex=0.6,  adj = aa)
                         circos.rect(xleft=xlim[1], ybottom=ylim[1],
                                     xright=xlim[2], ytop=ylim[2],
                                     col = df1$rcol[i], border=df1$rcol[i])
                         circos.rect(xleft=xlim[1], ybottom=ylim[1],
                                     xright=xlim[2]-rowSums(m)[i], ytop=ylim[1]+0.3,
                                     col = "white", border = "white")
                         circos.rect(xleft=xlim[1], ybottom=0.3,
                                     xright=xlim[2], ytop=0.32, col = "white",
                                     border = "white")
                         circos.axis(labels.cex=0.6, direction =
                                       "outside", major.at=seq(from=0,to=floor(df1$xmax)[i],by=5),
                                     minor.ticks=1,
                                     labels.away.percentage = 0.15)})

###########################################################################################################################
# Correlation of genes to phenotype.

pheno2<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Phenodata_vf/Harvest_VF_v2.txt",h=T)


# for AMF dataset
head(pheno2)
phenoA<-pheno2[pheno2$SAMPLE %in% colnames(FOUR_VARS_AMF),]
rownames(phenoA)<-phenoA$SAMPLE

head(phenoA[,c(6,11,13,15,14,27)])

#unique(as.vector(unlist(top_gene_sign_intramod_conec_A)))
#unique(as.vector(unlist(top_gene_sign_intramod_conec_Y)))
#NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  unique(as.vector(unlist(top_gene_sign_intramod_conec_Y))),])
#NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  unique(as.vector(unlist(top_gene_sign_intramod_conec_A))),])



pheatmap(t(cor(NORM_EXP_TOPAMF,phenoA[,c(11,13,15,14,27)])),cellheight = 12)
pheatmap(t(cor(NORM_EXP_TOPYUCA,phenoA[,c(11,13,15,14,27)])),cellheight = 12)
yuca_cor2pheno<-cor(NORM_EXP_TOPYUCA,phenoA[,c(11,13,15,14,27)])
yuca_cor2pheno[yuca_cor2pheno[,4]>0.2,]

amf_cor2pheno<-cor(NORM_EXP_TOPAMF,phenoA[,c(11,13,15,14,27)])
amf_cor2pheno[amf_cor2pheno[,4]>0.4,]
plot(t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  'g5922.t1',]),phenoA[,c(15)],pch=16)



NORM_EXP_TOPYUCA2<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.16G010100",])
NORM_EXP_TOPYUCA2<-NORM_EXP_TOPYUCA2[,grep("CTRL",colnames(NORM_EXP_TOPYUCA2),invert=T)]

plot(phenoA[,c(15)],NORM_EXP_TOPYUCA2,pch=16,cex=3,xlab="Manes.12G031900 PT1",ylab="Manes.06G174400 \nKT1",cex.lab=3,cex.axis=3)
abline(lm(as.vector(NORM_EXP_TOPYUCA2)~as.vector(phenoA[,c(15)])))
Mercator_Rirregularis[Mercator_Rirregularis$gene=="g5922.t1",]

###########################################################################################################################
###########################################################################################################################
##################RANDOM MODULEs CONTROLS AND NETWORKS

##CASSAVA

table(dynamicColors_Y)

Random_NET_cassava<-rep("a",length(dynamicColors_Y))
names(table(Random_NET_cassava))
table(Random_NET_cassava)

Random_NET_cassava[sample(length(Random_NET_cassava),table(dynamicColors_Y)[1],replace = F)]<-"Random1"

Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)][sample(length(Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)]),table(dynamicColors_Y)[2],replace = F)]<-"Random2"
Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)][sample(length(Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)]),table(dynamicColors_Y)[3],replace = F)]<-"Random3"
Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)][sample(length(Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)]),table(dynamicColors_Y)[4],replace = F)]<-"Random4"
Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)][sample(length(Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)]),table(dynamicColors_Y)[5],replace = F)]<-"Random5"
Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)][sample(length(Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)]),table(dynamicColors_Y)[6],replace = F)]<-"Random6"
Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)][sample(length(Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)]),table(dynamicColors_Y)[7],replace = F)]<-"Random7"
Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)][sample(length(Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)]),table(dynamicColors_Y)[8],replace = F)]<-"Random8"
Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)][sample(length(Random_NET_cassava[grep("Random",Random_NET_cassava,invert=T)]),table(dynamicColors_Y)[9],replace = F)]<-"Random9"


MEList_Y_RANDOM=moduleEigengenes(t(sign_norm_genes_Y),colors=Random_NET_cassava)
MEs_Y_RANDOM=MEList_Y_RANDOM$eigengenes


##AMF

table(dynamicColors_A)

Random_NET_AMF<-rep("a",length(dynamicColors_A))
names(table(Random_NET_AMF))
table(Random_NET_AMF)

Random_NET_AMF[sample(length(Random_NET_AMF),table(dynamicColors_A)[1],replace = F)]<-"Random1"

Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[2],replace = F)]<-"Random2"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[3],replace = F)]<-"Random3"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[4],replace = F)]<-"Random4"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[5],replace = F)]<-"Random5"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[6],replace = F)]<-"Random6"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[7],replace = F)]<-"Random7"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[8],replace = F)]<-"Random8"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[9],replace = F)]<-"Random9"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[10],replace = F)]<-"Random10"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[11],replace = F)]<-"Random11"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[12],replace = F)]<-"Random12"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[13],replace = F)]<-"Random13"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[14],replace = F)]<-"Random14"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[15],replace = F)]<-"Random15"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[16],replace = F)]<-"Random16"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[17],replace = F)]<-"Random17"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[18],replace = F)]<-"Random18"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[19],replace = F)]<-"Random19"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[20],replace = F)]<-"Random20"
Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)][sample(length(Random_NET_AMF[grep("Random",Random_NET_AMF,invert=T)]),table(dynamicColors_A)[21],replace = F)]<-"Random21"


MEList_A_RANDOM=moduleEigengenes(t(sign_norm_genes_A),colors=Random_NET_AMF)
MEs_A_RANDOM=MEList_A_RANDOM$eigengenes



moduleTraitCor_R= cor(MEs_Y_RANDOM,MEs_A_RANDOM,use="p")
moduleTraitPvalue_R=corPvalueStudent(moduleTraitCor_R,nSamples)
dim(moduleTraitPvalue_R)[2]
corModules2_R<-as.vector(as.matrix(moduleTraitPvalue_R))# TO CHANGE DEPENDING ON SIZE CORRELATION MATRIX
names(corModules2_R)<-c(paste(rep(colnames(moduleTraitPvalue_R)[1],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[2],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[3],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[4],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[5],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[6],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[7],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[8],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[9],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[10],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[11],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[12],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[13],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[14],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[15],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[16],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[17],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[18],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[19],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[20],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" "),
                      paste(rep(colnames(moduleTraitPvalue_R)[21],dim(moduleTraitPvalue_R)[1]),rownames(moduleTraitPvalue_R),sep=" ")
)
corModules_S_R<-corModules2_R[corModules2_R<0.001]
corModules_S_R

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/RANDOM_MODULES_CORRELATION_Y_A.pdf",height = 10,width = 10)

textMatrix=paste(signif(moduleTraitCor_R,2),"\n(",
                 signif(moduleTraitPvalue_R,1),")",sep="")

col_HM<-c(colorRampPalette( rev(brewer.pal(9, "Greens")) )(20),rep("white",10),colorRampPalette( brewer.pal(9, "Blues") )(20)  )
labeledHeatmap(Matrix = moduleTraitCor_R,
               xLabels = names(MEs_A_RANDOM),
               yLabels = names(MEs_Y_RANDOM),
               ySymbols = names(MEs_Y_RANDOM),
               xSymbols= names(MEs_A_RANDOM),
               colorLabels = FALSE,
               colors = col_HM,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,cex.lab=1,
               zlim = c(-1,1),
               main = paste("RANDOM Module YUCA- RANDOM module AMF relationships"))
dev.off()

# network graph RANDOM MODULES

library(igraph) 

nod<-cbind.data.frame(c(paste(rep("Y",length(table(Random_NET_cassava))),names(table(Random_NET_cassava)),sep="_"),
                        paste(rep("A",length(table(Random_NET_AMF))),names(table(Random_NET_AMF)) ,sep="_") ), 
                      c(names(table(Random_NET_cassava)),names(table(Random_NET_AMF))),
                      log(c(as.vector(table(Random_NET_cassava)),as.vector(table(Random_NET_AMF))) ),
                      c(as.vector(table(Random_NET_cassava)),as.vector(table(Random_NET_AMF))),
                      c(rep("CASSAVA",length(table(Random_NET_cassava))),rep("AMF",length(table(Random_NET_AMF))) ),
                      c(rep("1",length(table(Random_NET_cassava))),rep("2",length(table(Random_NET_AMF))) ))

colnames(nod)<-c("id","module","size","nsize","org","org.type")


lin<-cbind.data.frame(
  gsub("ME","",sapply(strsplit(names(corModules_S_R)," "), "[[", 1)),
  gsub("ME","",sapply(strsplit(names(corModules_S_R)," "), "[[", 2)),
  -log10(as.vector(corModules_S_R)   )   )

colnames(lin)<-c("from","to","weight")


net <- graph_from_data_frame(d=lin, vertices=nod, directed=F) 

# Generate colors based on media type:
V(net)$color <- c(rep("grey",length(table(Random_NET_cassava))),rep("brown",length(table(Random_NET_AMF))) )


# Set node size based on audience size:
V(net)$size <- V(net)$size*2
V(net)$label.cex = 1
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/RANDOM_MODULES_NETWORK_Y_A.pdf.pdf',useDingbats = F)

plot(net,vertex.label.dist=-0.4)
dev.off()

################################################################################################
# RANDOM MODULES
################################################################################################
# IN CASSAVA
################################################################################################
#module cyan CASSAVA

MEList_YUCA=moduleEigengenes(t(sign_norm_genes_Y),colors=dynamicColors_Y)
MEList_YUCA$validColors
exp_and_modules<-cbind.data.frame(sign_norm_genes_Y,MEList_YUCA$validColors)
exp_mod2<-exp_and_modules[sample(1:dim(exp_and_modules)[1],100),]


##GEne enrichment

### TOP GO SOLUTION

# 2) ANOTATION DATABASE
######################################################
## CASSAVA
locusName<-rownames(sign_YUCA)
sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
sign_feat_go<-merge(sign_feat2,mesculenta_go2,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,10)]
head(sign_feat_go)
write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
            col.names = F, row.names = F)
geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
######################################################
# 3) SELECT UNIVERSE
######################################################

#TOTAL UNIVERSE ALL LOCI

# TOTAL UNIVERSE CASSAVA
locusName<-rownames(sign_YUCA)
sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
sign_feat_go<-merge(sign_feat2,mesculenta_go2[,c(1,3)],by='locusName') # do not select transcriptname. working at gene level
sign_feat_go<-na.omit(sign_feat_go)
sign_feat_go<-unique(sign_feat_go)
length(unique(sign_feat_go$locusName))



In_module<-rownames(exp_mod2) # genes name in module
geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
names(geneList) <- sign_feat_go$locusName

sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5, ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)


resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)

write.table(allRes,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Modules_RANDOM100_Cassava.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = F)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Modules_RANDOM100_Cassava.pdf", width=14, height=10)
par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
barplot(allRes$Significant[allRes$classicFisher<0.05], names = paste(allRes$GO.ID[allRes$classicFisher<0.05],allRes$Term[allRes$classicFisher<0.05],sep=" "),
        xlab = "Nb. of significant genes",horiz=T)
dev.off()

RANDOM_CASSAVA<-allRes[allRes$classicFisher<0.05,]

dim(RANDOM_CASSAVA[RANDOM_CASSAVA$Significant>0,])[1]
dim(RANDOM_CASSAVA[RANDOM_CASSAVA$Significant>3,])[1]

############################################################################################################################
# IN AMF
############################################################################################################################
#module Magenta AMF
modules_AMF$A_MEmagenta
MEList_AMF=moduleEigengenes(t(sign_norm_genes_A),colors=dynamicColors_A)
MEList_AMF$validColors
exp_and_modules<-cbind.data.frame(sign_norm_genes_A,MEList_AMF$validColors)

dim(exp_and_modules)[2]
exp_mod2<-exp_and_modules[sample(1:dim(exp_and_modules)[1],100),] # random modules
colnames(exp_mod2)



##GEne enrichment

### TOP GO SOLUTION
#AMF
locusName<-rownames(sign_AMF)
sign_feat2<-cbind.data.frame(sign_AMF,locusName)
sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,31)]

write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
            col.names = F, row.names = F)
geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
######################################################
# 3) SELECT UNIVERSE
######################################################
#TOTAL UNIVERSE AMF
locusName<-rownames(sign_AMF)
sign_feat2<-cbind.data.frame(sign_AMF,locusName)
sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName') # do not select transcriptname. working at gene level
sign_feat_go<-na.omit(sign_feat_go)
sign_feat_go<-unique(sign_feat_go)
length(unique(sign_feat_go$locusName))



In_module<-rownames(exp_mod2) # genes name in module
geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
names(geneList) <- sign_feat_go$locusName

sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5, ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)


resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)

write.table(allRes,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Modules_RANDOM100_AMF.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = F)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Modules_RANDOM100_AMF.pdf", width=14, height=10)
par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
barplot(allRes$Significant[allRes$classicFisher<0.05], names = paste(allRes$GO.ID[allRes$classicFisher<0.05],allRes$Term[allRes$classicFisher<0.05],sep=" "),
        xlab = "Nb. of significant genes",horiz=T,cex.lab=2,,cex.axis=2)
dev.off()


RANDOM_AMF<-allRes[allRes$classicFisher<0.05,]

dim(RANDOM_AMF[RANDOM_AMF$Significant>0,])[1]
dim(RANDOM_AMF[RANDOM_AMF$Significant>3,])[1]

# figures Significant go terms.

# No gene number threshold

dim(RANDOM_AMF[RANDOM_AMF$Significant>0,])[1]
dim(RANDOM_CASSAVA[RANDOM_CASSAVA$Significant>0,])[1]

go_sig_Y<-list()
for (i in 1:length(GO_CASSAVA_Sign) ) {
  go_sig_Y[[i]]<-dim(GO_CASSAVA_Sign[[i]][GO_CASSAVA_Sign[[i]][,4]>0&GO_CASSAVA_Sign[[i]][,6]<0.05,])[1]
}

unlist(go_sig_Y)

go_sig_A<-list()
for (i in 1:length(GO_AMF_Sign) ) {
  go_sig_A[[i]]<-dim(GO_AMF_Sign[[i]][GO_AMF_Sign[[i]][,4]>0&GO_AMF_Sign[[i]][,6]<0.05,])[1]
  
}

unlist(go_sig_A)

namma<-c("Random Module AMF","Random Module Cassava",sapply(strsplit(names(GO_CASSAVA_Sign)," "), "[[", 1),sapply(strsplit(names(GO_CASSAVA_Sign)," "), "[[", 2))
namma<-gsub("A_","AMF_",namma)
namma<-gsub("Y_","Cassava_",namma)

go_0<-c(dim(RANDOM_AMF[RANDOM_AMF$Significant>0,])[1],dim(RANDOM_CASSAVA[RANDOM_CASSAVA$Significant>0,])[1]
,unlist(go_sig_A),unlist(go_sig_Y))
 
# With gene number threshold
dim(RANDOM_AMF[RANDOM_AMF$Significant>3,])[1]
dim(RANDOM_CASSAVA[RANDOM_CASSAVA$Significant>3,])[1]

go_sig_Y3<-list()
for (i in 1:length(GO_CASSAVA_Sign) ) {
  go_sig_Y3[[i]]<-dim(GO_CASSAVA_Sign[[i]][GO_CASSAVA_Sign[[i]][,4]>3&GO_CASSAVA_Sign[[i]][,6]<0.05,])[1]
  
}

unlist(go_sig_Y3)


go_sig_A3<-list()
for (i in 1:length(GO_AMF_Sign) ) {
  go_sig_A3[[i]]<-dim(GO_AMF_Sign[[i]][GO_AMF_Sign[[i]][,4]>3&GO_AMF_Sign[[i]][,6]<0.05,])[1]
  
}

unlist(go_sig_A3)

go_3<-c(dim(RANDOM_AMF[RANDOM_AMF$Significant>3,])[1],dim(RANDOM_CASSAVA[RANDOM_CASSAVA$Significant>03,])[1]
        ,unlist(go_sig_A3),unlist(go_sig_Y3))

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/GO_modules_and_random.pdf", width=10, height=8)
barplot(go_3[c(-4,-5,-6,-9,-15,-16,-18,-20,-32,-37,-27,-31,-33,-29,-34,-38,-30,-35,-28)],names.arg = namma[c(-4,-5,-6,-9,-15,-16,-18,-20,-32,-37,-27,-31,-33,-29,-34,-38,-30,-35,-28)],las=2)
dev.off()
############################################################################################
############################################################################################

####################################################################################################################

#sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,]),]
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes_Y
GO_CASSAVA_Sign_R<-list()
for (correlated_Modula in   names(corModules2_R))  {
  names(corModules_S)[1]
  gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)]#gsub("A_","",
  MEs_AMF$A_MEred
  ######################################################
  # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
  ######################################################
  #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
  geneModuleMembership_YUCA=as.data.frame(cor(t(sign_norm_genes),MEs_Y_RANDOM,use="p"))
  geneTraitSignificance_YUCA= as.data.frame(cor(t(sign_norm_genes), 
                                                MEs_A_RANDOM[,names(MEs_A_RANDOM)==sapply(strsplit(correlated_Modula," "), "[[", 1)],use="p")) # define correlation with external variable gsub("A_","",
  ######################################################
  # 1) MODULE DEFINITION
  ######################################################
  #Cassava
  ## CHANGE FOR THE COMPARISON DESIRED
  modNames=names(MEs_Y_RANDOM)
  module=gsub("ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))
  column= match(module,modNames)
  moduleGenes=Random_NET_cassava==module
  ######################################################
  # 2) ANOTATION DATABASE
  ######################################################
  ## CASSAVA
  locusName<-rownames(sign_YUCA)
  sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
  sign_feat_go<-merge(sign_feat2,mesculenta_go2,by='locusName')
  sign_feat_go<-na.omit(sign_feat_go)
  gene_2_go<-sign_feat_go[,c(1,10)]
  head(sign_feat_go)
  write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
              col.names = F, row.names = F)
  geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
  ######################################################
  # 3) SELECT UNIVERSE
  ######################################################
  
  #TOTAL UNIVERSE ALL LOCI
  
  # TOTAL UNIVERSE CASSAVA
  locusName<-rownames(sign_YUCA)
  sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
  sign_feat_go<-merge(sign_feat2,mesculenta_go2[,c(1,3)],by='locusName') # do not select transcriptname. working at gene level
  sign_feat_go<-na.omit(sign_feat_go)
  sign_feat_go<-unique(sign_feat_go)
  length(unique(sign_feat_go$locusName))
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
  ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #CASSAVA
  GS.6_Y<-geneModuleMembership_YUCA[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0 & abs(geneModuleMembership_YUCA[moduleGenes, match(sapply(strsplit(correlated_Modula," "), "[[", 2), names(geneModuleMembership_YUCA))])>0,]
  dim(GS.6_Y)
  GS.6_Y
  In_module<-rownames(GS.6_Y) # genes name in module CASSAVA
  
  most_correlated_genes<-geneModuleMembership_YUCA[moduleGenes,]
  
  na.omit(most_correlated_genes[,grep(sapply(strsplit(correlated_Modula," "), "[[", 2), colnames(most_correlated_genes)  )])
  
  
  
  
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  
  geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
  names(geneList) <- sign_feat_go$locusName
  
  sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5, ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
  
  
  sampleGOdata
  
  
  resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  RES_Y<-allRes
  
  GO_CASSAVA_Sign_R[[correlated_Modula]]<-RES_Y
  
   
}

######################################################  ######################################################
# IN AMF
######################################################  ######################################################

cassava9_randomM<-c(dim(GO_CASSAVA_Sign_R[[1]][GO_CASSAVA_Sign_R[[1]][,6]<0.05&GO_CASSAVA_Sign_R[[1]][,4]>3,])[1],
    dim(GO_CASSAVA_Sign_R[[2]][GO_CASSAVA_Sign_R[[2]][,6]<0.05&GO_CASSAVA_Sign_R[[2]][,4]>3,])[1],
        dim(GO_CASSAVA_Sign_R[[3]][GO_CASSAVA_Sign_R[[3]][,6]<0.05&GO_CASSAVA_Sign_R[[3]][,4]>3,])[1],
            dim(GO_CASSAVA_Sign_R[[4]][GO_CASSAVA_Sign_R[[4]][,6]<0.05&GO_CASSAVA_Sign_R[[4]][,4]>3,])[1],
                dim(GO_CASSAVA_Sign_R[[5]][GO_CASSAVA_Sign_R[[5]][,6]<0.05&GO_CASSAVA_Sign_R[[5]][,4]>3,])[1],
                    dim(GO_CASSAVA_Sign_R[[6]][GO_CASSAVA_Sign_R[[6]][,6]<0.05&GO_CASSAVA_Sign_R[[6]][,4]>3,])[1],
                        dim(GO_CASSAVA_Sign_R[[7]][GO_CASSAVA_Sign_R[[7]][,6]<0.05&GO_CASSAVA_Sign_R[[7]][,4]>3,])[1],
                            dim(GO_CASSAVA_Sign_R[[8]][GO_CASSAVA_Sign_R[[8]][,6]<0.05&GO_CASSAVA_Sign_R[[8]][,4]>3,])[1],
                                dim(GO_CASSAVA_Sign_R[[9]][GO_CASSAVA_Sign_R[[9]][,6]<0.05&GO_CASSAVA_Sign_R[[9]][,4]>3,])[1])

names(GO_CASSAVA_Sign_R)[c(1,2,3,4,5,6,7,8,9)]
cassava9_randomM_NAMES<-sapply(strsplit(names(GO_CASSAVA_Sign_R)[c(1,2,3,4,5,6,7,8,9)]," "), "[[", 2)

sign_norm_genes<-DGE_AMF_N$E
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
top_gene_sign_intramod_conec_A<-list()
GO_AMF_Sign_R<-list()
for (correlated_Modula in   names(corModules2_R)[c(1,100,127,136,145,154,163,172,181,10,19,28,37,46,55,64,73,82,91,109,118)])  {
  names(corModules_S)[1]
  gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)]
  
  
  ######################################################
  # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
  ######################################################
  #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
  geneModuleMembership_AMF=as.data.frame(cor(t(sign_norm_genes),MEs_A_RANDOM,use="p"))
  geneTraitSignificance_AMF= as.data.frame(cor(t(sign_norm_genes),
                                               MEs_Y[,names(MEs_Y_RANDOM)==sapply(strsplit(correlated_Modula," "), "[[", 2)],use="p")) # define correlation with external variable
  ######################################################
  # 1) MODULE DEFINITION
  ######################################################
  #AMF
  modNames=substring(names(MEs_A_RANDOM),5)
  module=gsub("ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  column= match(module,modNames)
  moduleGenes=Random_NET_AMF==module
  ######################################################
  # 2) ANOTATION DATABASE
  ######################################################
  #AMF
  locusName<-rownames(sign_AMF)
  sign_feat2<-cbind.data.frame(sign_AMF,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
  sign_feat_go<-na.omit(sign_feat_go)
  gene_2_go<-sign_feat_go[,c(1,31)]
  
  write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
              col.names = F, row.names = F)
  geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
  ######################################################
  # 3) SELECT UNIVERSE
  ######################################################
  #TOTAL UNIVERSE AMF
  locusName<-rownames(sign_AMF)
  sign_feat2<-cbind.data.frame(sign_AMF,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName') # do not select transcriptname. working at gene level
  sign_feat_go<-na.omit(sign_feat_go)
  sign_feat_go<-unique(sign_feat_go)
  length(unique(sign_feat_go$locusName))
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
  ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #AMF
  # change module match
  GS.6_A<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_AMF[moduleGenes,1])>0 & abs(geneModuleMembership_AMF[moduleGenes,match(substring(sapply(strsplit(correlated_Modula," "), "[[", 1),1), names(geneModuleMembership_AMF))])>0,]
  dim(GS.6_A)
  In_module<-rownames(GS.6_A) # genes name in module AMF
  
  ######################################################
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  
  geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
  names(geneList) <- sign_feat_go$locusName
  
  sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5, ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
  
  
  sampleGOdata
  
  
  resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  
  RES_A<-allRes
  GO_AMF_Sign_R[[correlated_Modula]]<-RES_A
  
  
  ######################################################
  # 6) PRINT RESULTS
  ######################################################
  ################
}

names(GO_AMF_Sign_R)
AMF21_randomM<-c(
dim(GO_AMF_Sign_R[[1]][GO_AMF_Sign_R[[1]][,6]<0.05&GO_AMF_Sign_R[[1]][,4]>3,])[1],
dim(GO_AMF_Sign_R[[2]][GO_AMF_Sign_R[[2]][,6]<0.05&GO_AMF_Sign_R[[2]][,4]>3,])[1],
    dim(GO_AMF_Sign_R[[3]][GO_AMF_Sign_R[[3]][,6]<0.05&GO_AMF_Sign_R[[3]][,4]>3,])[1],
        dim(GO_AMF_Sign_R[[4]][GO_AMF_Sign_R[[4]][,6]<0.05&GO_AMF_Sign_R[[4]][,4]>3,])[1],
            dim(GO_AMF_Sign_R[[5]][GO_AMF_Sign_R[[5]][,6]<0.05&GO_AMF_Sign_R[[5]][,4]>3,])[1],
                dim(GO_AMF_Sign_R[[6]][GO_AMF_Sign_R[[6]][,6]<0.05&GO_AMF_Sign_R[[6]][,4]>3,])[1],
                    dim(GO_AMF_Sign_R[[7]][GO_AMF_Sign_R[[7]][,6]<0.05&GO_AMF_Sign_R[[7]][,4]>3,])[1],
                        dim(GO_AMF_Sign_R[[8]][GO_AMF_Sign_R[[8]][,6]<0.05&GO_AMF_Sign_R[[8]][,4]>3,])[1],
dim(GO_AMF_Sign_R[[9]][GO_AMF_Sign_R[[9]][,6]<0.05&GO_AMF_Sign_R[[9]][,4]>3,])[1],
    dim(GO_AMF_Sign_R[[10]][GO_AMF_Sign_R[[10]][,6]<0.05&GO_AMF_Sign_R[[10]][,4]>3,])[1],
        dim(GO_AMF_Sign_R[[11]][GO_AMF_Sign_R[[11]][,6]<0.05&GO_AMF_Sign_R[[11]][,4]>3,])[1],
            dim(GO_AMF_Sign_R[[12]][GO_AMF_Sign_R[[12]][,6]<0.05&GO_AMF_Sign_R[[12]][,4]>3,])[1],
                dim(GO_AMF_Sign_R[[13]][GO_AMF_Sign_R[[13]][,6]<0.05&GO_AMF_Sign_R[[13]][,4]>3,])[1],
                    dim(GO_AMF_Sign_R[[14]][GO_AMF_Sign_R[[14]][,6]<0.05&GO_AMF_Sign_R[[14]][,4]>3,])[1],
                        dim(GO_AMF_Sign_R[[15]][GO_AMF_Sign_R[[15]][,6]<0.05&GO_AMF_Sign_R[[15]][,4]>3,])[1],
                            dim(GO_AMF_Sign_R[[16]][GO_AMF_Sign_R[[16]][,6]<0.05&GO_AMF_Sign_R[[16]][,4]>3,])[1],
 dim(GO_AMF_Sign_R[[17]][GO_AMF_Sign_R[[17]][,6]<0.05&GO_AMF_Sign_R[[17]][,4]>3,])[1],
     dim(GO_AMF_Sign_R[[18]][GO_AMF_Sign_R[[18]][,6]<0.05&GO_AMF_Sign_R[[18]][,4]>3,])[1],
         dim(GO_AMF_Sign_R[[19]][GO_AMF_Sign_R[[19]][,6]<0.05&GO_AMF_Sign_R[[19]][,4]>3,])[1],
             dim(GO_AMF_Sign_R[[20]][GO_AMF_Sign_R[[20]][,6]<0.05&GO_AMF_Sign_R[[20]][,4]>3,])[1],
                 dim(GO_AMF_Sign_R[[21]][GO_AMF_Sign_R[[21]][,6]<0.05&GO_AMF_Sign_R[[21]][,4]>3,])[1])

AMF21_randomM_NAMES<-sapply(strsplit(names(GO_AMF_Sign_R)," "), "[[", 1)


barplot(AMF21_randomM,names.arg = AMF21_randomM_NAMES, las=2)
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################

###########################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################

############################# ############################# ############################# 
############################# CONTROL WGCNA WITH 1 VAR 1 AMF ############################# 
############################# ############################# ############################# 

##### test var 5
########################################################################################################################################

FOUR_VARS_YUCA2<-FOUR_VARS_YUCA[,grep(c("V5"), colnames(FOUR_VARS_YUCA),invert =F)]

#take out controls for AMF
FOUR_VARS_AMF2<-FOUR_VARS_AMF[,grep(c("V5"), colnames(FOUR_VARS_AMF),invert =F)]
FOUR_VARS_AMF2<-FOUR_VARS_AMF[,grep(c("V5"), colnames(FOUR_VARS_AMF),invert =F)]
FOUR_VARS_AMF2<-FOUR_VARS_AMF2[,grep(c("CTRL"), colnames(FOUR_VARS_AMF2),invert =T)]


#######################################################################################################################################
# select significant compared to reference
#######################################################################################################################################
YUCA_sign_TREAT<-rownames(sign_YUCA_V1)
#normalized
sign_norm_genes_Y<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% YUCA_sign_TREAT,]



########################################################################################################################################



########################################################################################################################################
########################################    WGCNA MODULE COMPOSITIONS   ################################################################
########################################################################################################################################

## create eigene modules per sample. ##
###################################
#### M esculenta gene counts ######
###################################
#takeout controls
sign_norm_genes_Y<-sign_norm_genes_Y[,grep(c("CTRL"), colnames(sign_norm_genes_Y),invert =T)]
sign_norm_genes_Y<-sign_norm_genes_Y[,grep(c("V5"), colnames(sign_norm_genes_Y),invert =F)]


#normalized data
dim(sign_norm_genes_Y) # 1634 genes.
colnames(t(FOUR_VARS_YUCA_S))
#create set of powers
powers=c(c(1:10),seq(from=12, to=30,by=2))
#network topology analysis fonction, show power of networks
sft_Y=pickSoftThreshold(t(sign_norm_genes_Y),powerVector = powers,verbose = 5)
par(mfrow=c(1,2))
plot(sft_Y$fitIndices[,1], -sign(sft_Y$fitIndices[,3])*sft_Y$fitIndices[,2],xlab="SoftThresholdpower)",ylab="ScaleFreeTopologyModelFit,signed^2",main =paste("Scaleindependence"),type="n")
text(sft_Y$fitIndices[,1], -sign(sft_Y$fitIndices[,3])*sft_Y$fitIndices[,2],labels=powers)

plot(sft_Y$fitIndices[,1], sft_Y$fitIndices[,5],type="n")
text(sft_Y$fitIndices[,1], sft_Y$fitIndices[,5],labels=powers)

####################################################
## clustering methods**
####################################################

# create adjacency,
softPower = 18
adjacency_Y = adjacency(t(sign_norm_genes_Y),power= softPower);
#creat topological overlap matrix
TOM_Y= TOMsimilarity(adjacency_Y)
dissTOM_Y=1-TOM_Y
#Clustering using TOM
geneTree_Y= hclust(as.dist(dissTOM_Y),method="average")
#sizeGrWindow(12,9)
#plot(geneTree_Y,main="Gene clustering based on TOM dissimilarity",labels=F,hang=0.04)
####################################################
#definition of module size.
####################################################
minModuleSize=50
#module identification using dynamic tree cut. module
dynmicMods_Y= cutreeDynamic(dendro=geneTree_Y, distM=dissTOM_Y,deepSplit=2,pamRespectsDendro=F,minClusterSize=minModuleSize)
table(dynmicMods_Y)

dynamicColors_Y=labels2colors(dynmicMods_Y)
table(dynamicColors_Y)
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree_Y,dynamicColors_Y,"Dynamic Tree Cut",dendroLabels= F,hang=0.03,addGuide = T,guideHang = 0.05,main="Gene dendrogram and module colors")
####################################################
#merging modules whose expression profiles are similar
####################################################
#calculate eignegens
MEList_Y=moduleEigengenes(t(sign_norm_genes_Y),colors=dynamicColors_Y)
MEs_Y=MEList_Y$eigengenes
#calculate dissimilarity od module eigengenes
MEDiss_Y= 1 -cor(MEs_Y)
#cluster module eigengenes
METree_Y=hclust(as.dist(MEDiss_Y),method="average")
sizeGrWindow(7,6)
plot(METree_Y, main="Clustering of module eigengenes") # choose cut for a correlation of 0.75

MEDissThres=0.1
abline(h=MEDissThres,col="red")
merge= mergeCloseModules(t(sign_norm_genes_Y),dynamicColors_Y,cutHeight = MEDissThres,verbose=3)
mergedColors=merge$colors
mergedMEs_Y=merge$newMEs
# to see effect of merging
sizeGrWindow(12,9)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA1_modules_CASSAVA.pdf")
par(las=1)
plotDendroAndColors(geneTree_Y,mergedColors,"Modules",dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05,
                    cex.colorLabels = 2, cex.dendroLabels = 2, 
                    cex.rowText = 2)
dev.off()
#change names (to work better)
moduleColors_Y=mergedColors
colorOrder=c("grey",standardColors(50))
moduleLabels_Y=match(moduleColors_Y,colorOrder)-1
MEs_Y=mergedMEs_Y

###########################################################################################################################
# PREPARE DATA YUCA FOR quantifying module-module associations 

#define nb genes and samples
nGenes=ncol(t(sign_norm_genes_Y))
nSamples=nrow(t(sign_norm_genes_Y))
#recalculate MEs with color labels
MEs0_Y=moduleEigengenes(t(sign_norm_genes_Y),moduleColors_Y)$eigengenes
MEs_Y=orderMEs(MEs0_Y)
names(MEs_Y)<-paste("Y",names(MEs_Y),sep="_")
rownames(MEs_Y)<-rownames(t(sign_norm_genes_Y))



###########################################################################################################################
###########################################################################################################################
###################################
#### R irregularis gene modules####
###################################

#normalized data
sign_norm_genes_A<-DGE_AMF_N$E
dim(sign_norm_genes_A)
colnames(sign_norm_genes_A)
sign_norm_genes_A<-sign_norm_genes_A[,grep(c("V5"), colnames(sign_norm_genes_A),invert =F)]


## create eigene modules per sample. 
#create set of powers
powers=c(c(1:10),seq(from=12, to=40,by=2))
#network topology analysis fx
sft_A=pickSoftThreshold(t(sign_norm_genes_A),powerVector = powers,verbose = 5) #FOUR_VARS_AMF_S
par(mfrow=c(1,2))
plot(sft_A$fitIndices[,1], -sign(sft_A$fitIndices[,3])*sft_A$fitIndices[,2],xlab="SoftThresholdpower)",ylab="ScaleFreeTopologyModelFit,signed^2",main =paste("Scaleindependence"),type="n")
text(sft_A$fitIndices[,1], -sign(sft_A$fitIndices[,3])*sft_A$fitIndices[,2],labels=powers)

plot(sft_A$fitIndices[,1], sft_A$fitIndices[,5],type="n")
text(sft_A$fitIndices[,1], sft_A$fitIndices[,5],labels=powers)

####################################################
## clustering methods**
####################################################
#adjacency
adjacency_A = adjacency(t(sign_norm_genes_A),power= 8); # change power FOUR_VARS_AMF_S
#topological overlap dissimilarity matrix 
dissTOM_A=1-TOMsimilarity(adjacency_A)

#Clustering using TOM
geneTree_A= hclust(as.dist(dissTOM_A),method="average")
#sizeGrWindow(12,9)
#plot(geneTree_A,main="Gene clustering based on TOM dissimilarity",labels=F,hang=0.04)

####################################################
#definition of module size.
####################################################
minModuleSize=50
#module identification using dynamic tree cut
dynmicMods_A= cutreeDynamic(dendro=geneTree_A, distM=dissTOM_A,deepSplit=2,pamRespectsDendro=F,minClusterSize=minModuleSize)
dynamicColors_A=labels2colors(dynmicMods_A)
table(dynamicColors_A)

sizeGrWindow(8,6)
plotDendroAndColors(geneTree_A,dynamicColors_A,"Dynamic Tree Cut",dendroLabels= F,hang=0.03,addGuide = T,guideHang = 0.05,main="Gene dendrogram and module colors")

#merging modules whose expression profiles are similar
#calculate eignegens
MEList_AMF=moduleEigengenes(t(sign_norm_genes_A),colors=dynamicColors_A) #FOUR_VARS_AMF_S
MEs_AMF=MEList_AMF$eigengenes
#calculate dissimilarity od module eigengenes
MEDiss_AMF= 1 -cor(MEs_AMF)
#cluster module eigengenes
METree_AMF=hclust(as.dist(MEDiss_AMF),method="average")
sizeGrWindow(7,6)
plot(METree_AMF, main="Clustering of module eigengenes") # choose cut for a correlation of 0.75

MEDissThres=0.1
abline(h=MEDissThres,col="red")
merge_AMF= mergeCloseModules(t(sign_norm_genes_A),dynamicColors_A,cutHeight = MEDissThres,verbose=3)
mergedColors_AMF=merge_AMF$colors
mergedMEs_AMF=merge_AMF$newMEs
# to see effect of merging
sizeGrWindow(12,9)
plotDendroAndColors(geneTree_A,cbind(dynamicColors_A,mergedColors_AMF),c("Dynamic Tree cut","Merged Dynamic"),dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05)
table(mergedColors_AMF)
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA2_modules_AMF6.pdf")
par(las=1)
plotDendroAndColors(geneTree_A,mergedColors_AMF,"Modules",dendroLabels = F,hang=0.03,addGuide = T,guideHang = 0.05,
                    cex.dendroLabels = 2, cex.colorLabels = 2,
                    cex.rowText = 2)
dev.off()


#change names (to work better)
moduleColors_AMF=mergedColors_AMF
colorOrder=c("grey",standardColors(50))
moduleLabels_AMF=match(moduleColors_AMF,colorOrder)-1
MEs_AMF=mergedMEs_AMF

###########################################################################################################################
# PREPARE DATA AMF FOR quantifying module-module associations 

#define nb genes and samples
nGenes=ncol(t(sign_norm_genes_A))
nSamples=nrow(t(sign_norm_genes_A))
#recalculate MEs with color labels
MEs0_AMF=moduleEigengenes(t(sign_norm_genes_A),moduleColors_AMF)$eigengenes

MEs_AMF=orderMEs(MEs0_AMF)
names(MEs_AMF)<-paste("A",names(MEs_AMF),sep="_")
rownames(MEs_AMF)<-rownames(t(sign_norm_genes_A))
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

#Modules expression in samples
##########################################################################################################################################################
# EIGENGENES
# MOdules expression in different treatments.
# Mesculenta
########??RELATE MODULE EIGENGENES APPARTENSHIP TO TREATMENTS.
modules_YUCA<-MEs_Y

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA3_module_eigenenes_CASSAVA.pdf")
#rownames(modules_YUCA)<-rownames(t(FOUR_VARS_YUCA_topheno))
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
heatmap.2( as.matrix(modules_YUCA), col=colors,scale="none",dendrogram = "both",key = F)
dev.off()
heatmap.2( as.matrix(modules_YUCA), col=colors,scale="column",dendrogram = "both",key = F)


modules_AMF<-MEs_AMF

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA4_module_eigenenes_AMF.pdf")
#rownames(modules_AMF)<-rownames(t(FOUR_VARS_AMF_topheno))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( as.matrix(modules_AMF), col=colors,scale="none",dendrogram = "both", key = F )
dev.off()


################################################################################################
# RANDOM MODULES
################################################################################################
# IN CASSAVA
################################################################################################
#module cyan CASSAVA

MEList_YUCA=moduleEigengenes(t(sign_norm_genes_Y),colors=dynamicColors_Y)
MEList_YUCA$validColors
exp_and_modules<-cbind.data.frame(sign_norm_genes_Y,MEList_YUCA$validColors)
exp_mod2<-exp_and_modules[sample(1:dim(exp_and_modules)[1],100),]


##GEne enrichment

### TOP GO SOLUTION
locusName<-rownames(DGE_YUCA_N$E)
sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
sign_feat_go<-merge(sign_feat2,mesculenta_go,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,17)]

write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
            col.names = F, row.names = F)
geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")

In_module<-rownames(exp_mod2) # genes name in module
geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
names(geneList) <- sign_feat_go$locusName

sampleGOdata <-  tryCatch(new("topGOdata", ontology = "MF", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)


resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)

write.table(allRes,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Modules_RANDOM100_Cassava.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = F)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Modules_RANDOM100_Cassava.pdf", width=14, height=10)
par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
barplot(allRes$Significant[allRes$classicFisher<0.05], names = paste(allRes$GO.ID[allRes$classicFisher<0.05],allRes$Term[allRes$classicFisher<0.05],sep=" "),
        xlab = "Nb. of significant genes",horiz=T)
dev.off()


############################################################################################################################
# IN AMF
############################################################################################################################
#module Magenta AMF
modules_AMF$A_MEmagenta
MEList_AMF=moduleEigengenes(t(sign_norm_genes_A),colors=dynamicColors_A)
MEList_AMF$validColors
exp_and_modules<-cbind.data.frame(sign_norm_genes_A,MEList_AMF$validColors)

dim(exp_and_modules)[2]
exp_mod2<-exp_and_modules[sample(1:dim(exp_and_modules)[1],100),] # random modules
colnames(exp_mod2)



##GEne enrichment

### TOP GO SOLUTION
locusName<-rownames(DGE_AMF_N$E)
sign_feat2<-cbind.data.frame(sign_AMF,locusName)
sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,31)]

write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
            col.names = F, row.names = F)
geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")

In_module<-rownames(exp_mod2) # genes name in module
geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
names(geneList) <- sign_feat_go$locusName

sampleGOdata <-  tryCatch(new("topGOdata", ontology = "MF", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)


resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200)

write.table(allRes,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Modules_RANDOM100_AMF.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = F)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/Modules_RANDOM100_AMF.pdf", width=14, height=10)
par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
barplot(allRes$Significant[allRes$classicFisher<0.05], names = paste(allRes$GO.ID[allRes$classicFisher<0.05],allRes$Term[allRes$classicFisher<0.05],sep=" "),
        xlab = "Nb. of significant genes",horiz=T,cex.lab=2,,cex.axis=2)
dev.off()




############################################################################################
############################################################################################

##################################################################################################
####correlation Modules eigengenes YUCA,AMF #####
##################################################################################################
#takeout samples from CTRL in cassava.



moduleTraitCor= cor(MEs_Y,MEs_AMF,use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)


#Representation
sizeGrWindow(10,8)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",
                 signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)
par(mar=c(6,8.5,3,3))
#disply correlation values within heatmap plot
#X-axis co-expresion modules on AMF, Y-axis co-expresion modules CASSAVA
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V3/4.WGCNA5_CORRELATION_Modules_CASSAVA_AMF.pdf",height = 8,width = 8)
col_HM<-c(colorRampPalette( rev(brewer.pal(9, "Greens")) )(20),rep("white",10),colorRampPalette( brewer.pal(9, "Blues") )(20)  )
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(MEs_AMF),
               yLabels = names(MEs_Y),
               ySymbols = names(MEs_Y),
               xSymbols= names(MEs_AMF),
               colorLabels = FALSE,
               colors = col_HM,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,cex.lab=1,
               zlim = c(-1,1),
               main = paste("Module YUCA-module AMF relationships"))
dev.off()

##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
################## CORRELATION MODULES TO PHENOTYPE.

##################################################################################################

# Correlation of genes to phenotype.

pheno2<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Phenodata_vf/Harvest_VF_v2.txt",h=T)


# for AMF dataset
head(pheno2)
phenoA<-pheno2[pheno2$SAMPLE %in% colnames(FOUR_VARS_AMF),]
rownames(phenoA)<-phenoA$SAMPLE

head(phenoA[,c(6,11,13,15,14,27)])


### Correlation YUCA
moduleTraitCor=cor(MEs_Y,phenoA[,c(6,11:15,27)],use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",
                 signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)

    pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/info_for figures_v4/WGCNA7_Modules_CASSAVA_PHENOTYPE.pdf")
col_HM<-c(colorRampPalette( rev(brewer.pal(9, "Greens")) )(20),rep("white",10),colorRampPalette( brewer.pal(9, "Blues") )(20)  )
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels =  colnames(phenoA[,c(6,11:15,27)]),
               yLabels = names(MEs_Y),
               ySymbols = names(MEs_Y),
               colorLabels = FALSE,
               colors = col_HM,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,cex.lab=1,
               zlim = c(-1,1),
               main = paste("Module YUCA-Phenotype relationships"))
dev.off()

###############

dim(moduleTraitPvalue)[2]
corModules2_Y<-as.vector(as.matrix(moduleTraitPvalue))

names(corModules2_Y)<-c(paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[1],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[2],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[3],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[4],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[5],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[6],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[7],dim(moduleTraitPvalue)[2]),sep=" ")  )
corModules2_Y_S<-corModules2_Y[corModules2_Y<0.01]


top_gene_YUCA_pheno<-list()
Import_gene_YUCA_pheno<-list()

for (correlated_Modula in   names(corModules2_Y_S))  {
  
  ######################################################
  # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
  ######################################################
  #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
  geneModuleMembership_YUCA=as.data.frame(cor(t(sign_norm_genes_Y),MEs_Y,use="p"))
  geneTraitSignificance_YUCA= as.data.frame(cor(t(sign_norm_genes_Y),
                                               phenoA[,c(6)],use="p")) # define correlation with external variable
  ######################################################
  # 1) MODULE DEFINITION
  ######################################################
  #AMF
  modNames=substring(names(MEs_Y),5)
  module=gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  column= match(module,modNames)
  moduleGenes=moduleColors_Y==module
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
  ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #AMF
  # change module match
  GS.6_Y<-geneModuleMembership_YUCA[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0 & abs(geneModuleMembership_YUCA[moduleGenes,match(substring(sapply(strsplit(correlated_Modula," "), "[[", 1),1), names(geneModuleMembership_YUCA))])>0,]
  dim(GS.6_Y)
  In_module<-rownames(GS.6_Y) # genes name in module AMF
  
  ######################################################
  ## intramodular conectivity
  ADJ1=abs(cor(t(sign_norm_genes_Y),use="p"))^6
  Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors_Y)
  head(Alldegrees1)
  
  GS1=as.numeric(cor(phenoA[,c(6)],t(sign_norm_genes_Y),use="pairwise.complete.obs"))
  GeneSignificance=abs(GS1)
  #GS.6_Y<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0]
  # quantile(abs(geneTraitSignificance_YUCA[moduleGenes,1]),.85)
  table( abs(GS1[moduleGenes])> quantile(abs(GS1[moduleGenes]),.9))
  colorlevels=unique(dynamicColors_Y)
  
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/ALL_MODULES_SIGN',correlated_Modula,'.pdf',sep="_"),width=14,height=10,useDingbats = FALSE)
  par(mar=c(8,9,5,2),las=1,mgp=c(6, 2,0))
  
  whichmodule=gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1));
  restrict1 = (dynamicColors_Y==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=dynamicColors_Y[restrict1],
                     main=whichmodule,pch=16,cex=3, cex.axis=3,cex.lab = 3,cex.main = 3,
                     xlab = paste("Connectivity",whichmodule), ylab = paste(paste("Y",whichmodule,sep="_"),"correlation","to",gsub("Pheno","Y_",sapply(strsplit(correlated_Modula," "), "[[", 2)),sep=" "), abline = TRUE)
  
  dev.off()
  
  
  ####### # 
  #redo filter but with quantile and not absolute value
  #module black
  ##########
  datKME=signedKME(t(sign_norm_genes_Y), MEs_Y, outputColumnName="MM.")
  # FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
  FilterGenes=  abs(GS1)[moduleGenes]> quantile(abs(GS1)[moduleGenes],.9) & abs(datKME[moduleGenes,grep(gsub("Y_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
  top_sig_connect_high<-cbind.data.frame(rownames(datKME[moduleGenes,][FilterGenes,]),datKME[moduleGenes,grep(gsub("Y_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))][FilterGenes])
  as.vector(top_sig_connect_high[,1])
  
  table(FilterGenes)
  
  top_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-as.vector(top_sig_connect_high[,1]) #dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
  
  sign_YUCA_ALL<-topTable(fit_YUCA, coef=c(6:15),number=dim(DGE_YUCA_N)[1])
  Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<- sign_YUCA_ALL[rownames(sign_YUCA_ALL) %in%   as.vector(unlist(top_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])),]
  
  corr_list<-list()
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/YUCAgenes2pheno',correlated_Modula,'.pdf',sep=""),width=14,height=10,useDingbats = FALSE)
  par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0),mfrow=c(3,3))
  for( i in 1:length(top_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])){
    NORM_EXP_TOPYUCA<-t(sign_norm_genes_Y[rownames(sign_norm_genes_Y) %in% top_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],])
    tryCatch(corr_list[[i]]<-cor(phenoA[,c(6)],as.vector(NORM_EXP_TOPYUCA[,i]),use="pairwise.complete.obs"),error=function(e) 1)
    tryCatch(plot(NORM_EXP_TOPYUCA[,i],phenoA[,c(6)],pch=16,cex=3,#col=phenoA$var,
         ylab=sapply(strsplit(correlated_Modula," "), "[[", 2),xlab=top_gene_YUCA_pheno[[1]][i],cex.lab=3 ,cex.axis=3),error=function(e) 1)
    tryCatch(abline(lm(phenoA[,c(6)]~as.vector(NORM_EXP_TOPYUCA[,i]))),error=function(e) 1)
  }
  dev.off()
  Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<- cbind.data.frame( unlist(corr_list),  Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])
  colnames(Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[1]<-"Correlation"
  
  
  
  Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-cbind.data.frame(Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],rownames(Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]))
  colnames(Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[dim(Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[2]]<-"gene"
  Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-merge(Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],Mercator_Mesculenta,by="gene")
  
  toprint<-Import_gene_YUCA_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]
  
  toprint<-cbind( rep(sapply(strsplit(correlated_Modula," "), "[[", 1),dim(  toprint)[1]),
                  rep(sapply(strsplit(correlated_Modula," "), "[[", 2),dim(  toprint)[1]),
                  toprint[,c(1:11,15,16:19)])
  colnames(toprint)[1:2]<-c("YUCA_module","PHENO")
  write.table( toprint,
               paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/TOP_GENES_FUNCTION',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_"),
               ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  summary( toprint)
  
  paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")
  
  ################
  
  MEList_YUCA=moduleEigengenes(t(sign_norm_genes_Y),colors=dynamicColors_Y)
  exp_and_modules<-cbind.data.frame(sign_norm_genes_Y,MEList_YUCA$validColors)
  colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
  # extract top 20 genes
  exp_and_modules[exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1)),]
  cor_to_module<-as.data.frame(cor(t(sign_norm_genes_Y),
                                   phenoA[,c(6)]
                                   ,use="p"))
  top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
  colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
  
  top_correlated2<-top_correlated[ top_correlated$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))  ,]
  
  
  #### all markers + p-value for MapMan
  COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
  colnames(COR_info_module)<-c("gene",paste("cor2",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="_"))
  FC_YUCA<-cbind.data.frame(rownames(sign_YUCA),sign_YUCA)
  colnames(FC_YUCA)[1]<-'gene'
  COR_info_module<-merge(merge(COR_info_module,Mercator_Mesculenta,by="gene"),FC_YUCA,by='gene')
  
  write.table(  COR_info_module,
                paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/GENES_2mapman',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  
  #extract all genes in module
  exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))  ,]

  Genes_module<-t(sign_norm_genes_Y)[,colnames(t(sign_norm_genes_Y)) %in%  rownames(exp_mod2)]
  
  write.table(t(Genes_module), paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/ALL_GENES_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = T)
}

names(Import_gene_YUCA_pheno)
top_gene_YUCA_pheno[[3]]


NORM_EXP_TOPYUCA<-t(sign_norm_genes_Y[rownames(sign_norm_genes_Y) %in%  top_gene_YUCA_pheno[[1]],])
pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/YUCAgenes2pheno.pdf',width=14,height=10,useDingbats = FALSE)
par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0),mfrow=c(3,3))
for( i in 1:length(top_gene_amf_pheno[[1]])){
  plot(NORM_EXP_TOPYUCA[,i],phenoA[,c(6)],pch=16,cex=3,#col=phenoA$var,
       ylab="percentage colonization",xlab=top_gene_YUCA_pheno[[1]][i],cex.lab=3 ,cex.axis=3)
  abline(lm(phenoA[,c(6)]~as.vector(NORM_EXP_TOPYUCA[,i])))
}
dev.off()

cor(phenoA[,c(6)],as.vector(NORM_EXP_TOPYUCA[,2]),use="pairwise.complete.obs")
cor(phenoA[,c(6)],t(sign_norm_genes_Y))

Mercator_Mesculenta[Mercator_Mesculenta$gene %in% top_gene_YUCA_pheno[[3]],]

last_blastprotF<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/blast2go_annot_20160920_UniProtFungi.GOslim.annot',h=F)
colnames(last_blastprotF)<-c("gene","GO","Prot")
head(last_blastprotF)
last_blastprotF[last_blastprotF$gene %in%  top_gene_YUCA_pheno[[1]],]



########################################################################################
### Correlation AMF
moduleTraitCor=cor(MEs_AMF,phenoA[,c(6,11:15,27)],use="p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor,nSamples)
textMatrix=paste(signif(moduleTraitCor,2),"\n(",
                 signif(moduleTraitPvalue,1),")",sep="")
dim(textMatrix)=dim(moduleTraitCor)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/WGCNA8_Modules_AMF_PHENOTYPE.pdf")
col_HM<-c(colorRampPalette( rev(brewer.pal(9, "Greens")) )(20),rep("white",10),colorRampPalette( brewer.pal(9, "Blues") )(20)  )
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels =  colnames(phenoA[,c(6,11:15,27)]),
               yLabels = names(MEs_AMF),
               ySymbols = names(MEs_AMF),
               colorLabels = FALSE,
               colors = col_HM,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,cex.lab=1,
               zlim = c(-1,1),
               main = paste("Module YUCA-Phenotype relationships"))
dev.off()

###############

dim(moduleTraitPvalue)[2]
corModules2<-as.vector(as.matrix(moduleTraitPvalue))

names(corModules2)<-c(paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[1],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[2],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[3],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[4],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[5],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[6],dim(moduleTraitPvalue)[2]),sep=" "),
                      paste(rownames(moduleTraitPvalue),rep(colnames(moduleTraitPvalue)[7],dim(moduleTraitPvalue)[2]),sep=" ")  )



corModules2_S<-corModules2[corModules2<0.01]



top_gene_amf_pheno<-list()
Import_gene_AMF_pheno<-list()

for (correlated_Modula in   names(corModules2_S))  {
  
  ######################################################
  # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
  ######################################################
  #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
  geneModuleMembership_AMF=as.data.frame(cor(t(sign_norm_genes_A),MEs_AMF,use="p"))
  geneTraitSignificance_AMF= as.data.frame(cor(t(sign_norm_genes_A),
                                               phenoA[,c(6)],use="p")) # define correlation with external variable
  ######################################################
  # 1) MODULE DEFINITION
  ######################################################
  #AMF
  modNames=substring(names(MEs_AMF),5)
  module=gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))
  column= match(module,modNames)
  moduleGenes=moduleColors_AMF==module
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
  ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #AMF
  # change module match
  GS.6_A<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_AMF[moduleGenes,1])>0 & abs(geneModuleMembership_AMF[moduleGenes,match(substring(sapply(strsplit(correlated_Modula," "), "[[", 1),1), names(geneModuleMembership_AMF))])>0,]
  dim(GS.6_A)
  In_module<-rownames(GS.6_A) # genes name in module AMF
  
   ######################################################
  ## intramodular conectivity
  ADJ1=abs(cor(t(sign_norm_genes_A),use="p"))^6
  Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors_A)
  head(Alldegrees1)
  
  GS1=as.numeric(cor(phenoA[,c(6)],t(sign_norm_genes_A),use="pairwise.complete.obs"))
  GeneSignificance=abs(GS1)
  #GS.6_Y<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0]
  # quantile(abs(geneTraitSignificance_YUCA[moduleGenes,1]),.85)
  table( abs(GS1[moduleGenes])> quantile(abs(GS1[moduleGenes]),.9))
  colorlevels=unique(dynamicColors_A)
  
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/ALL_MODULES_SIGN',correlated_Modula,'.pdf',sep="_"),width=14,height=10,useDingbats = FALSE)
  par(mar=c(8,9,5,2),las=1,mgp=c(6, 2,0))
  
  whichmodule=gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1));
  restrict1 = (dynamicColors_A==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=dynamicColors_A[restrict1],
                     main=whichmodule,pch=16,cex=3, cex.axis=3,cex.lab = 3,cex.main = 3,
                     xlab = paste("Connectivity",whichmodule), ylab = paste(paste("A",whichmodule,sep="_"),"correlation","to",gsub("Y_ME","Y_",sapply(strsplit(correlated_Modula," "), "[[", 2)),sep=" "), abline = TRUE)
  
  dev.off()
  
  
  ####### # 
  #redo filter but with quantile and not absolute value
  #module black
  ##########
  datKME=signedKME(t(sign_norm_genes_A), MEs_AMF, outputColumnName="MM.")
  # FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
  FilterGenes=  abs(GS1)[moduleGenes]> quantile(abs(GS1)[moduleGenes],.9) & abs(datKME[moduleGenes,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
  
  top_sig_connect_high<-cbind.data.frame(rownames(datKME[moduleGenes,][FilterGenes,]),datKME[moduleGenes,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))][FilterGenes])
  as.vector(top_sig_connect_high[,1])
  
 
  table(FilterGenes)
  
  top_gene_amf_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-as.vector(top_sig_connect_high[,1]) #dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
  
  sign_AMF_ALL<-topTable(fit_AMF, coef=c(6:15),number=dim(DGE_AMF_N)[1])
  Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<- sign_AMF_ALL[rownames(sign_AMF_ALL) %in%   as.vector(unlist(top_gene_amf_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])),]
  
  
  corr_list<-list()
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/AMFgenes2pheno',correlated_Modula,'.pdf',sep=""),width=14,height=10,useDingbats = FALSE)
  par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0),mfrow=c(3,3))
  for( i in 1:length(top_gene_amf_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])){
    NORM_EXP_TOPAMF<-t(sign_norm_genes_A[rownames(sign_norm_genes_A) %in% top_gene_amf_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],])
    corr_list[[i]]<-cor(phenoA[,c(6)],as.vector(NORM_EXP_TOPAMF[,i]),use="pairwise.complete.obs")
    plot(NORM_EXP_TOPAMF[,i],phenoA[,c(6)],pch=16,cex=3,#col=phenoA$var,
         ylab=sapply(strsplit(correlated_Modula," "), "[[", 2),xlab=top_gene_amf_pheno[[1]][i],cex.lab=3 ,cex.axis=3)
    abline(lm(phenoA[,c(6)]~as.vector(NORM_EXP_TOPAMF[,i])))
  }
  dev.off()
  Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<- cbind.data.frame( unlist(corr_list),  Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])
  colnames(Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[1]<-"Correlation"
    
  Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-cbind.data.frame(Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],rownames(Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]))
  colnames(Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[dim(Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[2]]<-"gene"
  Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-merge(Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],Mercator_Rirregularis,by="gene")
  
  toprint<-Import_gene_AMF_pheno[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]
  
  toprint<-cbind( rep(sapply(strsplit(correlated_Modula," "), "[[", 1),dim(  toprint)[1]),
                  rep(sapply(strsplit(correlated_Modula," "), "[[", 2),dim(  toprint)[1]),
                  toprint[,c(1:11,15,16:19)])
  colnames(toprint)[1:2]<-c("AMF_module","YUCA_module")
  write.table( toprint,
               paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/TOP_GENES_FUNCTION',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_"),
               ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  summary( toprint)
  
  paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")
  
  ################
  
  MEList_AMF=moduleEigengenes(t(sign_norm_genes_A),colors=dynamicColors_A)
  exp_and_modules<-cbind.data.frame(sign_norm_genes_A,MEList_AMF$validColors)
  colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
  # extract top 20 genes
  exp_and_modules[exp_and_modules$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1)),]
  cor_to_module<-as.data.frame(cor(t(sign_norm_genes_A),
                                   phenoA[,c(6)]
                                   ,use="p"))
  top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
  colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
  
  top_correlated2<-top_correlated[ top_correlated$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))  ,]
  
  
  #### all markers + p-value for MapMan
  COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
  colnames(COR_info_module)<-c("gene",paste("cor2",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="_"))
  FC_AMF<-cbind.data.frame(rownames(sign_AMF),sign_AMF)
  colnames(FC_AMF)[1]<-'gene'
  COR_info_module<-merge(merge(COR_info_module,Mercator_Rirregularis,by="gene"),FC_AMF,by='gene')
  
  write.table(  COR_info_module,
                paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/GENES_2mapman',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  
  #extract all genes in module
  exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1))  ,]
  
  Genes_module<-t(sign_norm_genes_A)[,colnames(t(sign_norm_genes_A)) %in%  rownames(exp_mod2)]
  
  write.table(t(Genes_module), paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/info_for figures_v5/ALL_GENES_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = T)
}

names(Import_gene_AMF_pheno)
top_gene_amf_pheno[[1]]
as.vector(unlist(top_gene_amf_pheno))


Mercator_Rirregularis[Mercator_Rirregularis$gene %in% top_gene_amf_pheno[[3]],]

last_blastprotF<-read.delim('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/blast2go_annot_20160920_UniProtFungi.GOslim.annot',h=F)
colnames(last_blastprotF)<-c("gene","GO","Prot")
head(last_blastprotF)
last_blastprotF[last_blastprotF$gene %in%  top_gene_amf_pheno[[1]],]
