########################################################################################
#' Import MR bayes tree and make correlation to Differentially expressed genes
########################################################################################
#source("http://bioconductor.org/biocLite.R")

options("scipen"=100,"digits"=2)
install.packages('phyloch')

###############################################################
# import concensus tree
library(ape)
setwd('/Users/admin/Google Drive/Doctorat_shared_unil/Scripts_analysis/RNA-seq/AMFxCASSAVA_paper/Mr_bayes')

Cassava_tre<-read.nexus("mrbayesCASSAVA_reduced.nex.noclock.mb.con.tre")

plot(Cassava_tre)


plot(drop.tip(Cassava_tre,tip=Cassava_tre$tip.label[c(3:5,8:10,12:39)]))

Cassava_2<-drop.tip(Cassava_tre,tip=Cassava_tre$tip.label[c(3:5,8:10,12:39)])

cophenetic.phylo(Cassava_2)
cophenetic.phylo(Cassava_tre)

###############################################################
# do distance matrix of differences between samples

VAR<-gsub("m2","",gsub("_\\w+$"," ",colnames(FOUR_VARS_AMF),perl=T)       )
VAR<-factor(VAR,levels=c("V5 ","V1 ","V4 ","V6 ","V8 "))
VAR<-factor(VAR,levels=c("V1 ","V5 ","V4 ","V6 ","V8 "))
VAR<-factor(VAR,levels=c("V4 ","V5 ","V1 ","V6 ","V8 "))
VAR<-factor(VAR,levels=c("V6 ","V5 ","V1 ","V4 ","V8 "))
VAR<-factor(VAR,levels=c("V8 ","V5 ","V1 ","V4 ","V6 "))

lane <-VAR
lane<-as.character(lane)
lane[lane=="V1 "]<-"lane1"
lane[lane=="V5 "]<-"lane1"
lane[lane=="V8 "]<-"lane1"
lane[lane=="V4 "]<-"lane2"
lane[lane=="V6 "]<-"lane2"
lane<-as.factor(lane)
TREAT<-gsub("^(\\w+)_","",colnames(FOUR_VARS_AMF),perl=T)
TREAT_AMF<-factor(TREAT,levels=c("CAN","B1"))

design_AMF <- model.matrix(~lane+VAR+TREAT_AMF)

DGE_AMF_N <- voom(DGE_AMF, design_AMF,plot=TRUE)

fit_AMF <- eBayes(lmFit(DGE_AMF_N,design_AMF))
sign_AMF_P<-topTable(fit_AMF,c(3,4,5,6),number=dim(DGE_AMF_N$E)[1]) # all genes
rslt_AMF_P <- decideTests(fit_AMF[,c(3,4,5,6)],adjust.method="BH")

vennDiagram(rslt_AMF_P,include=c("up","down"),counts.col=c("#00441B","#08306B"),cex=c(2,2,1.5),main="DGE R. irregularis influenced by Cassava",circle.col=1)
rslt_AMF_P

V5_vs<-c(0,
length(rslt_AMF_P[,1][rslt_AMF_P[,1]!=0] ),
length(rslt_AMF_P[,2][rslt_AMF_P[,2]!=0] ),
length(rslt_AMF_P[,3][rslt_AMF_P[,3]!=0] ),
length(rslt_AMF_P[,3][rslt_AMF_P[,4]!=0] )
)

V1_vs<-c(
         length(rslt_AMF_P[,1][rslt_AMF_P[,1]!=0] ),
         0,
         length(rslt_AMF_P[,2][rslt_AMF_P[,2]!=0] ),
         length(rslt_AMF_P[,3][rslt_AMF_P[,3]!=0] ),
         length(rslt_AMF_P[,3][rslt_AMF_P[,4]!=0] )
)

V4_vs<-c(
  length(rslt_AMF_P[,1][rslt_AMF_P[,1]!=0] ),
  length(rslt_AMF_P[,2][rslt_AMF_P[,2]!=0] ),
  0,
  length(rslt_AMF_P[,3][rslt_AMF_P[,3]!=0] ),
  length(rslt_AMF_P[,3][rslt_AMF_P[,4]!=0] )
)

V6_vs<-c(
  length(rslt_AMF_P[,1][rslt_AMF_P[,1]!=0] ),
  length(rslt_AMF_P[,2][rslt_AMF_P[,2]!=0] ),
  length(rslt_AMF_P[,3][rslt_AMF_P[,3]!=0] ),
  0,
  length(rslt_AMF_P[,3][rslt_AMF_P[,4]!=0] )
)

V8_vs<-c(
  length(rslt_AMF_P[,1][rslt_AMF_P[,1]!=0] ),
  length(rslt_AMF_P[,2][rslt_AMF_P[,2]!=0] ),
  length(rslt_AMF_P[,3][rslt_AMF_P[,3]!=0] ),
  length(rslt_AMF_P[,3][rslt_AMF_P[,4]!=0] ),
  0
)

expression_diff<-rbind.data.frame(V5_vs,V1_vs,V4_vs,V6_vs,V8_vs)
colnames(expression_diff)<-c("m2V513","V1_2","V4_9","V6_5","V8_5")
rownames(expression_diff)<-c("m2V513","V1_2","V4_9","V6_5","V8_5")


hclust(expression_diff)
###############################################################
cophenetic.phylo(Cassava_2)
expression_diff
library(ade4)
mantel.test(cophenetic.phylo(Cassava_2),expression_diff,
            graph=T) 

