########################################################################################
#' Script to compare different methods of Go enrichments
#' 
########################################################################################
#' Needs already in memorory: annotations, DEG genes names and p-values, Universe genes. modules info.
#' 
#' Run befor RNA_seq_gene_ll_AMF_HOST_v6.R
#' 
#' ########################################################################################################################################
################################### correspondance to go. using  PHytozome annotation        ########################################### 
########################################################################################################################################

######################################################
# 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
######################################################

#' MODULES CORRELATED BETWEEN CASSAVA AND AMF
corModules_S
# number genes in modules
table(dynamicColors_Y)
table(dynamicColors_A)
# modules
MEs_Y
MEs_AMF

corModules_S[1]
gsub("A_ME","",sapply(strsplit(names(corModules_S)[1]," "), "[[", 1))

ORIGIN_MODULE_CASSAVA = gsub("Y_ME","",sapply(strsplit(names(corModules_S)[1]," "), "[[", 2))
ORIGIN_MODULE_AMF = gsub("A_ME","",sapply(strsplit(names(corModules_S)[1]," "), "[[", 1))
  
  TARGET_MODULE_CASSAVA = sapply(strsplit(names(corModules_S)[1]," "), "[[", 1)
  TARGET_MODULE_AMF =  sapply(strsplit(names(corModules_S)[1]," "), "[[", 2)
  
#take out controls of CASSAVA
FOUR_VARS_YUCA_S<-FOUR_VARS_YUCA_S[,grep(c("CTRL"), colnames(FOUR_VARS_YUCA_S),invert =T)]

#' MODULES INFO MM GS
#' 
#' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
geneModuleMembership_YUCA=as.data.frame(cor(t(FOUR_VARS_YUCA_S),MEs_Y,use="p"))
geneTraitSignificance_YUCA= as.data.frame(cor(t(FOUR_VARS_YUCA_S),MEs_AMF$MEred,use="p")) # define correlation with external variable

geneModuleMembership_AMF=as.data.frame(cor(t(FOUR_VARS_AMF_S),MEs_AMF,use="p"))
geneTraitSignificance_AMF= as.data.frame(cor(t(FOUR_VARS_AMF_S),MEs_Y$Y_MEturquoise,use="p")) # define correlation with external variable



######################################################
# 1) MODULE DEFINITION
######################################################
#Cassava
## CHANGE FOR THE COMPARISON DESIRED
modNames=substring(names(MEs_Y),5)
module=gsub("Y_ME","",sapply(strsplit(names(corModules_S)[1]," "), "[[", 2))
column= match(module,modNames)
moduleGenes=moduleColors_Y==module
#AMF
modNames=substring(names(MEs_AMF),5)
module=gsub("A_ME","",sapply(strsplit(names(corModules_S)[1]," "), "[[", 1))
column= match(module,modNames)
moduleGenes=moduleColors_AMF==module
 

######################################################
# 2) ANOTATION DATABASE
######################################################
## CASSAVA
locusName<-rownames(sign_YUCA)
sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
sign_feat_go<-merge(sign_feat2,mesculenta_go2,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,9)]
head(sign_feat_go)
write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
            col.names = F, row.names = F)
geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")

#AMF
locusName<-rownames(sign_AMF)
sign_feat2<-cbind.data.frame(sign_AMF,locusName)
sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
sign_feat_go<-na.omit(sign_feat_go)
gene_2_go<-sign_feat_go[,c(1,8)]

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
  #TOTAL UNIVERSE AMF
  locusName<-rownames(sign_AMF)
  sign_feat2<-cbind.data.frame(sign_AMF,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName') # do not select transcriptname. working at gene level
  sign_feat_go<-na.omit(sign_feat_go)
  sign_feat_go<-unique(sign_feat_go)
  length(unique(sign_feat_go$locusName))
  

# OTHER UNIVERSE. NOT USED
 ## DEG SIGNFICANT 
#  locusName<-rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])
#  sign_feat2<-cbind.data.frame(sign_YUCA[sign_YUCA$adj.P.Val<0.05,],locusName)
#  sign_feat_go<-merge(sign_feat2,mesculenta_go2[,c(1,3)],by='locusName') # do not select transcriptname. working at gene level
#  sign_feat_go<-na.omit(sign_feat_go)
#  sign_feat_go<-unique(sign_feat_go)
#  length(unique(sign_feat_go$locusName))
  
  ## GENES IN MODULE
#    locusName<-rownames(na.omit(geneModuleMembership_YUCA[moduleGenes,])) #take out [] if not work
#  sign_feat2<-cbind.data.frame(locusName,na.omit(abs(geneTraitSignificance_YUCA[moduleGenes,1])))#take out [] if not work
#  sign_feat_go<-merge(sign_feat2,mesculenta_go2,by='locusName')
#  sign_feat_go<-na.omit(sign_feat_go)
#  sign_feat_go<-unique(sign_feat_go)
  
  
  
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
  ## apply filter of 0.6 to data gene sign and module membership
            ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
  #CASSAVA
  GS.6_Y<-geneModuleMembership_YUCA[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0.4 & abs(geneModuleMembership_YUCA[moduleGenes, match(sapply(strsplit(names(corModules_S)[1]," "), "[[", 2), names(geneModuleMembership_YUCA))])>0.4,]
  dim(GS.6_Y)
  In_module<-rownames(GS.6_Y) # genes name in module CASSAVA
  


  #AMF
  # change module match
  GS.6_A<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_AMF[moduleGenes,1])>0.4 & abs(geneModuleMembership_AMF[moduleGenes,match(substring(sapply(strsplit(names(corModules_S)[1]," "), "[[", 1),3), names(geneModuleMembership_AMF))])>0.4,]
  dim(GS.6_A)
  In_module<-rownames(GS.6_A) # genes name in module AMF
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  
  geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
  names(geneList) <- sign_feat_go$locusName
  
  sampleGOdata <-  tryCatch(new("topGOdata",nodeSize = 5, ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
 
  
  sampleGOdata
  
  
    resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  RES_Y<-allRes
  
  
  RES_A<-allRes
  
  
  ############## PRINT RESULTS
  
  write.table(RES_Y,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_TURQUOISE_CASSAVA_GO_ENRICH.txt',quote=FALSE,sep='\t',
              col.names = T, row.names = F)
  
  write.table(RES_A,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_RED_AMF_GO_ENRICH.txt',quote=FALSE,sep='\t',
              col.names = T, row.names = F)
  
  #png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_TURQUOISE_CASSAVA_GO_ENRICH.png", width=1400, height=1000, units="px")
  #par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  barplot(RES_Y$Significant[RES_Y$classicFisher<0.05], names = paste(RES_Y$GO.ID[RES_Y$classicFisher<0.05],RES_Y$Term[RES_Y$classicFisher<0.05],sep=" "),
          xlab = "Nb. of significant genes",horiz=T,cex.lab=2,las=1,cex.axis=2)
  #dev.off()
  
#  png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_RED_AMF_GO_ENRICH.png", width=1400, height=1000, units="px")
 # par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  barplot(RES_A$Significant[RES_A$classicFisher<0.05], names = paste(RES_A$GO.ID[RES_A$classicFisher<0.05],RES_A$Term[RES_A$classicFisher<0.05],sep=" "),
          xlab = "Nb. of significant genes",horiz=T,cex.lab=2,las=1,cex.axis=2)
  #dev.off()
  

  to_pheat_Y<-RES_Y[RES_Y$classicFisher<0.05,c(4)]
  names(to_pheat_Y)<-paste(RES_Y[RES_Y$classicFisher<0.05,c(1)],RES_Y[RES_Y$classicFisher<0.05,c(2)],sep=" ")
  pheatmap(t(to_pheat_Y),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10)
  
  
  to_pheat_A<-RES_A[RES_A$classicFisher<0.05,c(4)]
  names(to_pheat_A)<-paste(RES_A[RES_A$classicFisher<0.05,c(1)],RES_A[RES_A$classicFisher<0.05,c(2)],sep=" ")
  pheatmap(t(to_pheat_A),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10)
  
 
  
    ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  
  
  
  
  ################################################################################################################################
  # GO ANALYSIS FOR DIFFERENTIALLY EXPRESSED
  # CASSAVA
  ################################################################################################################################

  ######################################################
  # 2) ANOTATION DATABASE
  ######################################################
  ## CASSAVA
  locusName<-rownames(sign_YUCA)
  sign_feat2<-cbind.data.frame(sign_YUCA,locusName)
  sign_feat_go<-merge(sign_feat2,mesculenta_go2,by='locusName')
  sign_feat_go<-na.omit(sign_feat_go)
  gene_2_go<-sign_feat_go[,c(1,9)]
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
  ## DEG SIGNFICANT 
  Significant_GENES_Y<-rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])
  
##### diff CAN b1 pval < 0.005
Dif_YUCA_t<-cbind.data.frame(rownames(rslt_YUCA@.Data),rslt_YUCA@.Data)
Diff_YUCA_CAN<-rownames(Dif_YUCA_t[Dif_YUCA_t$CAN!=0&Dif_YUCA_t$B1==0,])# Different to ctrl specific for CAN
Diff_YUCA_B1<-rownames(Dif_YUCA_t[Dif_YUCA_t$B1!=0&Dif_YUCA_t$CAN==0,])# Different to ctrl specific for CAN
Significant_GENES_YCANB1spe<-c(Diff_YUCA_CAN,Diff_YUCA_B1)
  Significant_GENES_Y2<-rownames(sign_YUCA2[sign_YUCA2$adj.P.Val<0.05,])
  GENES_YCANB1spe<-sign_YUCA[rownames(sign_YUCA) %in% Significant_GENES_YCANB1spe,]
  GENES_YCANB1spe[GENES_YCANB1spe$adj.P.Val<0.005,]
#####
  
  rand_genes<-sample(sign_feat_go$locusName,dim(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])[1],replace=F )
  rand_genes<-sample(sign_feat_go$locusName,100,replace=F )
  
  #exp_mod2<-exp_and_modules[sample(1:dim(exp_and_modules)[1],100),]
  write.table(sign_YUCA2[sign_YUCA2$adj.P.Val<0.05,],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/CASSAVA_DEG.txt',quote=FALSE,sep='\t',
              col.names = T, row.names = T)

  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  
  geneList <- factor(as.integer(sign_feat_go$locusName %in% Significant_GENES_Y))
  #geneList <- factor(as.integer(sign_feat_go$locusName %in% rownames(GENES_YCANB1spe[GENES_YCANB1spe$adj.P.Val<0.01,])))
  
  #geneList <- factor(as.integer(sign_feat_go$locusName %in% Significant_GENES_Y2))
  
  #geneList <- factor(as.integer(sign_feat_go$locusName %in% rand_genes))
  
  names(geneList) <- sign_feat_go$locusName
  
  sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5,ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
  
  
  sampleGOdata
  
  
  resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  RES_Y<-allRes
  

  
  ############## PRINT RESULTS
  
  write.table(RES_Y,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/2.GO_CASSAVA_SPE_CANB1.txt',quote=FALSE,sep='\t',
              col.names = T, row.names = F)
  
  write.table(RES_Y[RES_Y[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/2.GO_CASSAVA_SPE_CANB1.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = F, row.names = F)
  
  #png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_TURQUOISE_CASSAVA_GO_ENRICH.png", width=1400, height=1000, units="px")
  #par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  barplot(RES_Y$Significant[RES_Y$classicFisher<0.05], names = paste(RES_Y$GO.ID[RES_Y$classicFisher<0.05],RES_Y$Term[RES_Y$classicFisher<0.05],sep=" "),
          xlab = "Nb. of significant genes",horiz=T,cex.lab=2,las=1,cex.axis=2)
  #dev.off()
  
  
  to_pheat_Y<-RES_Y[RES_Y$classicFisher<0.05,c(4)]
  names(to_pheat_Y)<-paste(RES_Y[RES_Y$classicFisher<0.05,c(1)],RES_Y[RES_Y$classicFisher<0.05,c(2)],sep=" ")
  pheatmap(t(to_pheat_Y),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10)
  
  ################################################################################################################################
  # GO ANALYSIS FOR DIFFERENTIALLY EXPRESSED
  # AMF
  ################################################################################################################################
  
  ######################################################
  # 2) ANOTATION DATABASE
  ######################################################
  #AMF
  locusName<-rownames(sign_AMF)
  sign_feat2<-cbind.data.frame(sign_AMF,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
  sign_feat_go<-na.omit(sign_feat_go)
  gene_2_go<-sign_feat_go[,c(1,8)]
  
  write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
              col.names = F, row.names = F)
  geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
  
  
  ######################################################
  # 3) SELECT UNIVERSE
  ######################################################
  
  #TOTAL UNIVERSE ALL LOCI
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
  ## DEG SIGNFICANT 
  Significant_GENES_A<-rownames(sign_AMF[sign_AMF$adj.P.Val<0.05,])

  Significant_GENES_BYPLANT<-rownames(sign_AMF_P[sign_AMF_P$adj.P.Val<0.05,])
  
  rand_genes<-sample(sign_feat_go$locusName,dim(sign_AMF[sign_AMF$adj.P.Val<0.05,])[1],replace=F )
  rand_genes<-sample(sign_feat_go$locusName,100,replace=F )
  
  #exp_mod2<-exp_and_modules[sample(1:dim(exp_and_modules)[1],100),]
  
  #write.table(sign_AMF_P[sign_AMF_P$adj.P.Val<0.05,],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/AMF_by_plant_DEG.txt',quote=FALSE,sep='\t',
              col.names = T, row.names = T)
  
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  
  geneList <- factor(as.integer(sign_feat_go$locusName %in% Significant_GENES_A))
  #geneList <- factor(as.integer(sign_feat_go$locusName %in% rand_genes))
  #geneList <- factor(as.integer(sign_feat_go$locusName %in% Significant_GENES_BYPLANT))
  
  names(geneList) <- sign_feat_go$locusName
  
  sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5,ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
  
  
  sampleGOdata
  
  
  resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  RES_A<-allRes
  
  
  
  ############## PRINT RESULTS
  
  write.table(RES_A,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/1.GO_AMF.txt',quote=FALSE,sep='\t',
              col.names = T, row.names = F)
  write.table(RES_A[RES_A[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/GO6_AMF.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = F, row.names = F)
  
  #png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_TURQUOISE_CASSAVA_GO_ENRICH.png", width=1400, height=1000, units="px")
  #par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  barplot(RES_A$Significant[RES_A$classicFisher<0.05], names = paste(RES_A$GO.ID[RES_A$classicFisher<0.05],RES_A$Term[RES_A$classicFisher<0.05],sep=" "),
          xlab = "Nb. of significant genes",horiz=T,cex.lab=2,las=1,cex.axis=2)
  #dev.off()
  
  
  to_pheat_A<-RES_A[RES_A$classicFisher<0.05,c(4)]
  names(to_pheat_A)<-paste(RES_A[RES_A$classicFisher<0.05,c(1)],RES_A[RES_A$classicFisher<0.05,c(2)],sep=" ")
  pheatmap(t(to_pheat_A),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10)

    ################################################################################################################################
  # GO ANALYSIS FOR DIFFERENTIALLY EXPRESSED
  # AMF BY PLANT
  ################################################################################################################################
  
  ######################################################
  # 2) ANOTATION DATABASE
  ######################################################
  #AMF
  locusName<-rownames(sign_AMF_P)
  sign_feat2<-cbind.data.frame(sign_AMF_P,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName')
  sign_feat_go<-na.omit(sign_feat_go)
  gene_2_go<-sign_feat_go[,c(1,10)]
  
  write.table(gene_2_go,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',quote=FALSE,sep='\t',
              col.names = F, row.names = F)
  geneID2GO <- readMappings('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/kallisto_Mesculenta/gene2go_test.txt',sep = "\t", IDsep = ",")
  
  
  ######################################################
  # 3) SELECT UNIVERSE
  ######################################################
  
  #TOTAL UNIVERSE ALL LOCI
  #TOTAL UNIVERSE AMF
  locusName<-rownames(sign_AMF_P)
  sign_feat2<-cbind.data.frame(sign_AMF_P,locusName)
  sign_feat_go<-merge(sign_feat2,rirregularis_go,by='locusName') # do not select transcriptname. working at gene level
  sign_feat_go<-na.omit(sign_feat_go)
  sign_feat_go<-unique(sign_feat_go)
  length(unique(sign_feat_go$locusName))
  ######################################################
  # 4) SELECT TOP GENES
  ######################################################
  ## DEG SIGNFICANT 

  Significant_GENES_BYPLANT<-rownames(sign_AMF_P[sign_AMF_P$adj.P.Val<0.05,])
  
  rand_genes<-sample(sign_feat_go$locusName,dim(sign_AMF[sign_AMF$adj.P.Val<0.05,])[1],replace=F )
  
  ######################################################
  # 5) ANNOTATION
  ######################################################
  #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
  
  #geneList <- factor(as.integer(sign_feat_go$locusName %in% rand_genes))
  geneList <- factor(as.integer(sign_feat_go$locusName %in% Significant_GENES_BYPLANT))
  
  names(geneList) <- sign_feat_go$locusName
  
  sampleGOdata <-  tryCatch(new("topGOdata", nodeSize = 5,ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
  
  
  
  sampleGOdata
  
  
  resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
  resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
  resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
  
  allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                             classicKS = resultKS, elimKS = resultKS.elim,
                             orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
  
  RES_A<-allRes
  
  
  
  ############## PRINT RESULTS
  
  write.table(RES_A,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/3.GO_AMF_byplant.txt',quote=FALSE,sep='\t',
              col.names = T, row.names = F)
  write.table(RES_A[RES_A[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/GO6_AMF.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = F, row.names = F)
  
  #png("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_MODULE_TURQUOISE_CASSAVA_GO_ENRICH.png", width=1400, height=1000, units="px")
  #par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  barplot(RES_A$Significant[RES_A$classicFisher<0.05], names = paste(RES_A$GO.ID[RES_A$classicFisher<0.05],RES_A$Term[RES_A$classicFisher<0.05],sep=" "),
          xlab = "Nb. of significant genes",horiz=T,cex.lab=2,las=1,cex.axis=2)
  #dev.off()
  
  
  to_pheat_A<-RES_A[RES_A$classicFisher<0.05,c(4)]
  names(to_pheat_A)<-paste(RES_A[RES_A$classicFisher<0.05,c(1)],RES_A[RES_A$classicFisher<0.05,c(2)],sep=" ")
  pheatmap(t(to_pheat_A),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10)
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  
  
  
  
  
  ################################################################################################################################
  # GO ANALYSIS IN MODULES
  ################################################################################################################################

sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,]),]
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
top_gene_sign_intramod_conec_Y<-list()
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
    gene_2_go<-sign_feat_go[,c(1,9)]
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
 
       
    write.table(RES_Y,    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/GO_ENRICH_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'GO_ENRICH.txt',sep="_")
      ,quote=FALSE,sep='\t',col.names = T, row.names = F)
    write.table(RES_Y[RES_Y[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/GO_2REVIGO_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'GO_ENRICH.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = F, row.names = F)
    
    
    pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/BARPLOT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'.pdf',sep="_"), width=10, height=7)
    par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
    tryCatch( barplot(RES_Y$Significant[RES_Y$classicFisher<0.05], names = paste(RES_Y$GO.ID[RES_Y$classicFisher<0.05],RES_Y$Term[RES_Y$classicFisher<0.05],sep=" "),
            xlab = "Nb. of significant genes",horiz=T,las=1),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 2)) )
    dev.off()

    
    to_pheat_Y<-RES_Y[RES_Y$classicFisher<0.05,c(4)]
    names(to_pheat_Y)<-paste(RES_Y[RES_Y$classicFisher<0.05,c(1)],RES_Y[RES_Y$classicFisher<0.05,c(2)],sep=" ")
    
    pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/PHEAT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'.pdf',sep="_"), width=18, height=10)
    tryCatch(pheatmap(t(to_pheat_Y),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 2)))
    dev.off()
  }
    
  ######################################################  ######################################################
  # IN AMF
    ######################################################  ######################################################
  names(corModules_S)[c(1,10,13,16,19,24,26,27,28,34,35)]
  unique(sapply(strsplit(names(corModules_S)," "), "[[", 1))
  
  
  sign_norm_genes<-DGE_AMF_N$E
  #take out controls of CASSAVA
  sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
  top_gene_sign_intramod_conec_A<-list()
    for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_AMF))  {
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
    gene_2_go<-sign_feat_go[,c(1,8)]
    
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

    
    ######################################################
    # 6) PRINT RESULTS
    ######################################################
 
    
    ################

    
    
    
    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/ENRICH_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.txt',sep="_")
    write.table(RES_A,    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/GO_ENRICH_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'GO_ENRICH.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
    write.table(RES_A[RES_A[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/GO_2REVIGO_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'GO_ENRICH.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = F, row.names = F)

    
    pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/BARPLOT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.pdf',sep="_"), width=10, height=7)
    par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
    tryCatch( barplot(RES_A$Significant[RES_A$classicFisher<0.05], names = paste(RES_A$GO.ID[RES_A$classicFisher<0.05],RES_A$Term[RES_A$classicFisher<0.05],sep=" "),
                      xlab = "Nb. of significant genes",horiz=T,las=1),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 1)) )
    dev.off()
    
    
    to_pheat_A<-RES_A[RES_A$classicFisher<0.05,c(4)]
    names(to_pheat_A)<-paste(RES_A[RES_A$classicFisher<0.05,c(1)],RES_A[RES_A$classicFisher<0.05,c(2)],sep=" ")
    
    pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/PHEAT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.pdf',sep="_"), width=18, height=10)
    tryCatch(pheatmap(t(to_pheat_A),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 1)))
    dev.off()
        
    
    }
    
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ######################################################
  #' 1) LOAD MERCATOR DATABASE
  ######################################################
  #' 1) LOAD DATABASE
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
  
  
  
##############################################################################################################################
  ############# PLOT INFO AND PHEATMAP OF DIFFERENT PARTS
  
  ### Cassava total trancriptome influenced CAN B1
  
  GENES_YUCA<-sign_YUCA[sign_YUCA$adj.P.Val<0.05,]
  colnames(GENES_YUCA)[1:2]<-c("logFC CAN-CTRL","logFC B1-CTRL")
  pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/2.Pheat_YUCA_vsCTRL_genes.pdf',width=10,height=10)
  DEG_A<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])]
  pheatmap(DEG_A,show_colnames = F)
  dim(  GENES_YUCA)
  dev.off()
  GENES_YUCA$gene<-rownames(GENES_YUCA)
  GENES_YUCA_2<-merge(GENES_YUCA,Mercator_Mesculenta,by="gene")
  head(GENES_YUCA_2[,c(1:3,7,9:10)])
  write.table( GENES_YUCA_2[,c(1:3,6,9:10)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/2.Genes_YUCAvsCTRL.txt',
               quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  
   table(GENES_YUCA_2[,9])[table(GENES_YUCA_2[,9])>3]
  
  head( GENES_YUCA)
  
   ### AMF trancriptome CAN B1

GENES_AMF<-sign_AMF[sign_AMF$adj.P.Val<0.05,]
  colnames(GENES_AMF)[1]<-c("logFC CAN-B1")
  pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/1.Pheat_AMF_genes.pdf',width=10,height=10)
  DEG_A<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  rownames(sign_AMF[sign_AMF$adj.P.Val<0.05,])]
  pheatmap(DEG_A,show_colnames = F)
  dim(  GENES_AMF)
  dev.off()
  GENES_AMF$gene<-rownames(GENES_AMF)
  GENES_AMF_2<-merge(GENES_AMF,Mercator_Rirregularis,by="gene")
  head(GENES_AMF_2[,c(1:2,6,9:10)])
  write.table( GENES_AMF_2[,c(1:2,6,9:10)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/1.Genes_AMF.txt',
               quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  
  table(GENES_AMF_2[,9])[table(GENES_AMF_2[,9])>10]
  
  
  ################################################################################################################################
  ### Cassava different genes between CAN B1
  Significant_GENES_YCANB1spe
  
  GENES_YCANB1spe<-sign_YUCA[rownames(sign_YUCA) %in% Significant_GENES_YCANB1spe,]
  colnames(GENES_YCANB1spe)[1:2]<-c("logFC CAN-CTRL","logFC B1-CTRL")
  pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/Fig1_b_CANB1_genes.pdf',width=6,height=6)
  pheatmap(t(GENES_YCANB1spe[GENES_YCANB1spe$adj.P.Val<0.01,1:2]),cellheight = 20,show_colnames = F)
  tail(  GENES_YCANB1spe[GENES_YCANB1spe$adj.P.Val<0.01,])
  dim(  GENES_YCANB1spe[GENES_YCANB1spe$adj.P.Val<0.01,])
  dev.off()
  GENES_YCANB1spe$gene<-rownames(GENES_YCANB1spe)
  GENES_YCANB1spe_2<-merge(GENES_YCANB1spe[GENES_YCANB1spe$adj.P.Val<0.01,],Mercator_Mesculenta,by="gene")
  head(GENES_YCANB1spe_2)
  write.table( GENES_YCANB1spe_2[,c(1:3,7,9:10)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/2.Genes_CASSAVA_SPE_CANB1.txt',
               quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  
  table(GENES_YCANB1spe_2[,9])[table(GENES_YCANB1spe_2[,9])>3]
  
  ################################################################################################################################
  # AMF BY PLANT
  GENES_AbyP<-sign_AMF_P[rownames(sign_AMF_P) %in% Significant_GENES_BYPLANT,]
  colnames(GENES_AbyP)[1:4]<-c("logFC Cultivar 5-Cultivar 8","logFC Cultivar 1-Cultivar 8","logFC Cultivar 4-Cultivar 8","logFC Cultivar 6-Cultivar 8")
  pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/3.AMF_by_plant.pdf',width=10,height=10)
  pheatmap(t(GENES_AbyP[GENES_AbyP$adj.P.Val<0.05,1:4]),cellheight = 20,show_colnames = F)
  tail(  GENES_AbyP[GENES_AbyP$adj.P.Val<0.01,])
  dim(  GENES_AbyP[GENES_AbyP$adj.P.Val<0.05,])
  dev.off()
  GENES_AbyP$gene<-rownames(GENES_AbyP)
  GENES_AbyP<-merge(GENES_AbyP[GENES_AbyP$adj.P.Val<0.01,],Mercator_Rirregularis,by="gene")
  head(GENES_AbyP)
  write.table( GENES_AbyP[,c(1:5,9,11:12)],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/3.Genes_AMF_by_Plant.txt',
               quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  
  
  ##############################################################################################################################
  ##############################################################################################################################
  ##############################################################################################################################
  
  ################################################################################################################################
  # extract more representative and more significant genes, no go terms
  #fastest
  ################################################################################################################################
  

  sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,]),]
  #take out controls of CASSAVA
  sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
  top_gene_sign_intramod_conec_Y<-list()
  Import_gene_M_Y<-list()
  for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_AMF)[10])  {

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
    ADJ1=abs(cor(t(sign_norm_genes),use="p"))^6
    Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors_Y)
    head(Alldegrees1)
    
    GS1=as.numeric(cor(MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)],t(sign_norm_genes), use="p"))
    GeneSignificance=abs(GS1)
    colorlevels=unique(dynamicColors_Y)
   
    pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/ALL_MODULES_SIGN_Y',correlated_Modula,'.pdf',sep="_"),width=14,height=10,useDingbats = FALSE )
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
    as.vector(top_sig_connect_high[,1])
    
    
    
    table(FilterGenes)
    top_gene_sign_intramod_conec_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-as.vector(top_sig_connect_high[,1]) #dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
    
    Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-sign_YUCA[rownames(sign_YUCA) %in%   as.vector(unlist(top_gene_sign_intramod_conec_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]])),]
    Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-cbind.data.frame(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]],rownames(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]))
    colnames(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]])[dim(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]])[2]]<-"gene"
    Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]<-merge(Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]],Mercator_Mesculenta,by="gene")
    
    toprint<-Import_gene_M_Y[[paste(sapply(strsplit(correlated_Modula," "), "[[", 2),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="")]]
     write.table( toprint[,c(1:3,7,9:10)],
                 paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/TOP_GENES_FUNCTION',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 1),'.txt',sep="_"),
                 ,quote=FALSE,sep='\t',col.names = T, row.names = T)
    
                 ################
    
    MEList_YUCA=moduleEigengenes(t(sign_norm_genes),colors=dynamicColors_Y)
    exp_and_modules<-cbind.data.frame(sign_norm_genes,MEList_YUCA$validColors)
    colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
    
    # extract top 20 genes
    exp_and_modules[exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2)),]
    cor_to_module<-as.data.frame(cor(t(sign_norm_genes),
                                     MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)]
                                     ,use="p"))
    top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
    colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
    top_correlated2<-top_correlated[ top_correlated$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
    
    
    #### all markers + p-value for MapMan
    COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
    colnames(COR_info_module)<-c("gene",paste("cor2",sapply(strsplit(correlated_Modula," "), "[[", 1),sep="_"))
    FC_YUCA<-cbind.data.frame(rownames(sign_YUCA),sign_YUCA)
    colnames(FC_YUCA)[1]<-'gene'
    COR_info_module<-merge(merge(COR_info_module,Mercator_Mesculenta,by="gene"),FC_YUCA,by='gene')
    write.table(  COR_info_module,
                  paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/GENES_2mapman',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 1),'.txt',sep="_")
                  ,quote=FALSE,sep='\t',col.names = T, row.names = F)
    
    
    #do barplot
    exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
    Genes_module<-t(sign_norm_genes)[,colnames(t(sign_norm_genes)) %in%  rownames(exp_mod2)]
    write.table(t(Genes_module),    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/ALL_GENES_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'GO_ENRICH.txt',sep="_")
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
  §unique(sapply(strsplit(names(corModules_S)," "), "[[", 1))
  
  
  sign_norm_genes<-DGE_AMF_N$E
  #take out controls of CASSAVA
  sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
  top_gene_sign_intramod_conec_A<-list()
  Import_gene_M_A<-list()
  
  for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_AMF)[10])  {
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
    ADJ1=abs(cor(t(sign_norm_genes),use="p"))^6
    Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors_A)
    head(Alldegrees1)
    
    GS1=as.numeric(cor(MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)],t(sign_norm_genes), use="p"))
    GeneSignificance=abs(GS1)
    #GS.6_Y<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0]
    # quantile(abs(geneTraitSignificance_YUCA[moduleGenes,1]),.85)
   table( abs(GS1[moduleGenes])> quantile(abs(GS1[moduleGenes]),.9))
    colorlevels=unique(dynamicColors_A)
    sizeGrWindow(9,6)
    par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
    par(mar = c(4,5,3,1))
    pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/ALL_MODULES_SIGN',correlated_Modula,'.pdf',sep="_"),width=14,height=10,useDingbats = FALSE)
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
     datKME=signedKME(t(sign_norm_genes), MEs_AMF, outputColumnName="MM.")
   # FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
     FilterGenes=  abs(GS1)[moduleGenes]> quantile(abs(GS1)[moduleGenes],.9) & abs(datKME[moduleGenes,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))]) >.8
     
     top_sig_connect_high<-cbind.data.frame(rownames(datKME[moduleGenes,][FilterGenes,]),datKME[moduleGenes,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))][FilterGenes])
     as.vector(top_sig_connect_high[,1])
     
     table(FilterGenes)
    
    top_gene_sign_intramod_conec_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-as.vector(top_sig_connect_high[,1]) #dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]

    Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-sign_AMF[rownames(sign_AMF) %in%   as.vector(unlist(top_gene_sign_intramod_conec_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])),]
    Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-cbind.data.frame(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],rownames(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]))
    colnames(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[dim(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]])[2]]<-"gene"
    Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-merge(Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]],Mercator_Rirregularis,by="gene")
    
    toprint<-Import_gene_M_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]
    write.table( toprint[,c(1:2,6,9:10)],
                 paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/TOP_GENES_FUNCTION',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_"),
                 ,quote=FALSE,sep='\t',col.names = T, row.names = T)
    
    summary( toprint)
    
    paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")
    
    ################
    
    MEList_AMF=moduleEigengenes(t(sign_norm_genes),colors=dynamicColors_A)
    exp_and_modules<-cbind.data.frame(sign_norm_genes,MEList_AMF$validColors)
    colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
    # extract top 20 genes
    exp_and_modules[exp_and_modules$MODULE==gsub("A_ME","",sapply(strsplit(correlated_Modula," "), "[[", 1)),]
    cor_to_module<-as.data.frame(cor(t(sign_norm_genes),
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
                  paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/GENES_2mapman',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'COR_TO',sapply(strsplit(correlated_Modula," "), "[[", 2),'.txt',sep="_")
                  ,quote=FALSE,sep='\t',col.names = T, row.names = F)
    
    #extract all genes in module
    exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub("Y_ME","",sapply(strsplit(correlated_Modula," "), "[[", 2))  ,]
    
    Genes_module<-t(sign_norm_genes)[,colnames(t(sign_norm_genes)) %in%  rownames(exp_mod2)]
    
    write.table(t(Genes_module), paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/ALL_GENES_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  }
  
  names(Import_gene_M_A)
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  
################################################################################################################################
 ################################################################################################################################
  # CORRELATION PLOT BETWEEN TOP GENES IN AMF CASSAVA MODULE
  ################################################################################################################################
  
    top_gene_sign_intramod_conec_A
  top_gene_sign_intramod_conec_Y
  
  pos_cor<-list()
  neg_cor<-list()
  1:length(names(top_gene_sign_intramod_conec_Y))
  for (ca in  1:length(names(top_gene_sign_intramod_conec_Y))   ) {
  top_AMF_gene<-top_gene_sign_intramod_conec_A[[ca]]
  top_YUCA_gene<-top_gene_sign_intramod_conec_Y[[ca]]
  
  names(top_gene_sign_intramod_conec_Y)
  names(top_gene_sign_intramod_conec_A)
  
  NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  top_AMF_gene,])
  
  NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  top_YUCA_gene,])
  NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
  
#########################  
  # table info about modules and top genes
  total_genes_AMF<-gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_A),"_"), "[[", 2))
  
  for ( i in 1:length(table(mergedColors_AMF))) {
    total_genes_AMF[total_genes_AMF==names(table(mergedColors_AMF))[i]]<- as.vector(table(mergedColors_AMF)[i])
  }
  
  total_genes_YUCA<-gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_Y),"_"), "[[", 2))
  for ( i in 1:length(table(mergedColors))) {
    total_genes_YUCA[total_genes_YUCA==names(table(mergedColors))[i]]<- as.vector(table(mergedColors)[i])
  }
  # AMF
  nanaA<-gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_A),"_"), "[[", 2))
  tnanaA<-table(gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_A),"_"), "[[", 2)))
  for ( i in 1:length(tnanaA)) {
    nanaA[nanaA==names(tnanaA)[i]]<-c("Carboxylic acid metabolism","Negative regulation of RNA metabolism",
                                      "Meiotic cell cycle / Epithellium development / Phospholipid metabolism",
                                      "Small molecule catabolism","Multicellular organism development/process",
                                      "Transcription from RNA polymerase II promoter",
                                      "Aspartate family amino acid biosynthesis","Cell wall organization/ cell migration",
                                      "Nucleocytoplasmic transport","Quinone metabolism",
                                      "Negative regulation of protein catabolism")[i]    }
  
  
  # yuca
  nanaY<-gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_Y),"_"), "[[", 2))
  tnanaY<-table(gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_Y),"_"), "[[", 2)))
  for ( i in 1:length(tnanaY)) {
    nanaY[nanaY==names(tnanaY)[i]]<-c("ER to Golgi vesicle-mediated transport","Cell response to extra-cellular stimulus",
                                      "Sulfur aminoacid metabolism","Anatomical structure formation involved in morphogenesis",
                                      "Nucleotide biosynthesis","Mismatch repair","Riboflavin metabolism","Porphyrin-containing compound biosynthesis",
                                      "Macromolecule localization","L-phenylalamine biosynthesis","Photorespiration/proton transport",
                                      "Dicarboxylic acid transport","Negative regulation of cellular component organization",
                                      "Cellular homeostasis","Carbohydrate derivative biosynthesis/ Glycosylation","Cell wall biogenesis",
                                      "Monocarboxylic acid biosynthesis")[i]  }
  
    
 
 
 correlacion<-tryCatch(cor(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA),error=function(e) 1)
 
 pos_cor[[ca]]<-table(correlacion>0 )[names(table(correlacion>0 ))=="TRUE"]
 
 
 neg_cor[[ca]]<-table(correlacion>0 )[names(table(correlacion>0 ))=="FALSE"]

 lalass_neg=rep(0,length(unlist(pos_cor)))
 for (fa in 1:length(unlist(pos_cor)) ) {
   tryCatch(lalass_neg[fa]<-neg_cor[[fa]],error=function(e) 0)
 }
 

  #########################  
   #install.packages('corrgram')
  #install.packages('ellipse')
  #library(corrgram)
  #library(ellipse)
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/CorTOPgenes_',names(top_gene_sign_intramod_conec_Y)[ca] ,'.pdf',sep="_"),width=18,height=14)
  
  correlacion<-tryCatch(cor(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA),error=function(e) 1)
  
  colsc=c(rgb(241, 54, 23, maxColorValue=255),  "white" , rgb(0, 61, 104, maxColorValue=255))
  colramp = colorRampPalette(colsc, space= 'Lab' )
  colors = colramp(100)
  # function
  #https://hlplab.wordpress.com/2012/03/20/correlation-plot-matrices-using-the-ellipse-library/
  tryCatch(labo<-data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA),error=function(e) 1)
  #[,(dim(NORM_EXP_TOPAMF)[2]-4):dim(NORM_EXP_TOPAMF)[2]],labo[,(dim(NORM_EXP_TOPYUCA)[2]-4):dim(NORM_EXP_TOPYUCA)[2]]
  #tryCatch(my.plotcorr(cex =2,cex.lab = 2,cor(labo)[grep('g',colnames(labo)),grep('M',colnames(labo))], col=colors[((cor(labo)[grep('g',colnames(labo)),grep('M',colnames(labo))] + 1)/2) * 100], diag= 'ellipse' , upper.panel="number", mar=c(0,2,0,0),main=paste("Modules",names(top_gene_sign_intramod_conec_A)[i]) ),error=function(e) plot(1,1, main=names(top_gene_sign_intramod_conec_A)[i]) )
  tryCatch(pheatmap(correlacion,col=colramp(100),fontsize = 15,cellwidth = 15,cellheight = 15,main=paste("Modules",names(top_gene_sign_intramod_conec_A)[i])),error=function(e) plot(1,1, main=names(top_gene_sign_intramod_conec_A)[i]) )
  
  #tryCatch( pairs(labo[,c((dim(NORM_EXP_TOPAMF)[2]-4):dim(NORM_EXP_TOPAMF)[2], (dim(NORM_EXP_TOPYUCA)[2]-4):dim(NORM_EXP_TOPYUCA)[2])], pch=16,col="black",main=paste("Modules",names(top_gene_sign_intramod_conec_A)[i]) ),error=function(e) plot(1,1, main=names(top_gene_sign_intramod_conec_A)[i]) )
  
  dev.off()
  
  }
  
  #pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/CorTOPgenes_',names(top_gene_sign_intramod_conec_Y)[ca] ,'.pdf',sep="_"),width=8,height=12)
  #pheatmap(correlacion,col=colramp(100),fontsize = 8,cellwidth = 8,cellheight = 8,main=paste("Modules",names(top_gene_sign_intramod_conec_A)[i]))
  #dev.off()
  # final table

  Crosstalks_df_vf<-cbind.data.frame(gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_A),"_"), "[[", 2)),
                   gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_Y),"_"), "[[", 2)),
                   total_genes_AMF,total_genes_YUCA, nanaA,nanaY,
                   as.vector(unlist(lapply(top_gene_sign_intramod_conec_A,length))),
                   as.vector(unlist(lapply(top_gene_sign_intramod_conec_Y,length))),
                   unlist(pos_cor),unlist(lalass_neg)
                   )
  colnames(Crosstalks_df_vf)<-c("Modules Cassava","Modules AMF","Total genes Cassava","Total genes AMF",
                                "Principal GO Cassava","Principal GO AMF","Top genes Cassava","Top genes AMF",
                                "Nb. positive correlations","Nb. negative correlations")
  
  write.table(Crosstalks_df_vf,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/4.Crosstalks_df_vf.txt' ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  
   ####### ####### ####### ####### #######
  ##### ALL modules info (not onyl correlated ones)
  corModules2
  # table info about modules and top genes
  total_genes_AMF<-gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_A),"_"), "[[", 2))
  
  for ( i in 1:length(table(mergedColors_AMF))) {
    total_genes_AMF[total_genes_AMF==names(table(mergedColors_AMF))[i]]<- as.vector(table(mergedColors_AMF)[i])
  }
  
  total_genes_YUCA<-gsub("ME","",sapply(strsplit(names(top_gene_sign_intramod_conec_Y),"_"), "[[", 2))
  for ( i in 1:length(table(mergedColors))) {
    total_genes_YUCA[total_genes_YUCA==names(table(mergedColors))[i]]<- as.vector(table(mergedColors)[i])
  }
  
  
  ####### ####### ####### ####### #######
  tryCatch(pheatmap(t(correlacion),col=colramp(100),fontsize = 10,cellwidth = 10,cellheight = 10,main=paste("Modules",names(top_gene_sign_intramod_conec_A)[i])),error=function(e) plot(1,1, main=names(top_gene_sign_intramod_conec_A)[i]) )
  
  ######   ######  ######  ######  ###### ######
  # manual correlation from two genes

    NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  "g2342.t1",])
  NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.09G039400",])
  NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
  pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/remark_cor_genes_v1.pdf',width=14,height=10,useDingbats = FALSE)
  par(mar=c(8,11,5,2),las=1,mgp=c(5.5, 2,0))
  plot(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA,pch=16,cex=3,ylab="Manes.09G039400 \n endopeptidase inhibitor activity",xlab="g2342.t1 peptidase S24/S26A-C",cex.lab=3,cex.axis=3)
  abline(lm(NORM_EXP_TOPYUCA~as.vector(NORM_EXP_TOPAMF)),lwd=5)
  dev.off()
  
  NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  "g12728.t1",])
  NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  "Manes.09G039400",])
  NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[,grep("CTRL",colnames(NORM_EXP_TOPYUCA),invert=T)]
  pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/remark_cor_genes_v1b.pdf',width=14,height=10,useDingbats = FALSE)
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
  
################################################################################################################################
  # how many times genes are found in different modules?
  top_gene_sign_intramod_conec_A
  top_gene_sign_intramod_conec_Y
  top_gene_sign_intramod_conec_A[[2]]
  
  table(unlist(top_gene_sign_intramod_conec_A))
  table(unlist(top_gene_sign_intramod_conec_Y))
  
  length(table(unlist(top_gene_sign_intramod_conec_Y)))
  length(table(unlist(top_gene_sign_intramod_conec_A)))
  
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/TOPGENES_hist.pdf',sep="_"),width=10,height=6)
  par(mfrow=c(1,2))
  hist(table(unlist(top_gene_sign_intramod_conec_A)),main="TOP AMF GENES IN MODULES",xlab="Number of modules found for gene")
  hist(table(unlist(top_gene_sign_intramod_conec_Y)),main="TOP CASSAVA GENES IN MODULES",xlab="Number of modules found for gene")
  dev.off()
  
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/TOPGENES_times in other comparison.pdf',sep="_"),width=10,height=6)
  par(mfrow=c(2,1))
  
  names(top_gene_sign_intramod_conec_A)
  pheatmap( t(as.matrix(table(unlist(top_gene_sign_intramod_conec_A)))),cluster_rows = F,cluster_cols = T,cellheight= 8)
  pheatmap( t(as.matrix(table(unlist(top_gene_sign_intramod_conec_Y)))),cluster_rows = F,cluster_cols = T,cellheight= 8)
  dev.off()
  
  
  ###############
  ## how many time a gene is seen in a different comparison to other organism module.
  # gene specificity to one comparison or multiple comparisons
  
all_top_genes_R<-list()
  for (i in unique(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)," "), "[[", 1) ) ) {
    pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/TYPE_TOP_GENES/TYPE_TOP_GENES',i,'.pdf',sep="_"),width=10,height=6)
    
    top_mod1<-top_gene_sign_intramod_conec_A[grep(i,names(top_gene_sign_intramod_conec_A))]
  length(top_mod1)
  tryCatch(pheatmap( t(as.matrix(table(unlist(top_mod1))))/length(top_mod1),cluster_rows = F,
            main=paste(i,'nb. comparisons=',length(top_mod1)),cluster_cols = T,cellheight= 8),error=function(e) pheatmap( t(c((t(as.matrix(table(unlist(top_mod1))))/length(top_mod1))[1],(t(as.matrix(table(unlist(top_mod1))))/length(top_mod1))[2]*.9)),cluster_rows = F,
             main=i,cluster_cols = T,cellheight= 8) )    
  all_top_genes_R[[i]]<-colnames(t(as.matrix(table(unlist(top_mod1)))) ) 
  dev.off()
  general_top_gen<-cbind.data.frame(
    as.vector(table(unlist(top_mod1))[table(unlist(top_mod1))>(length(top_mod1)-2)])/length(top_mod1),
    names(table(unlist(top_mod1))[table(unlist(top_mod1))>(length(top_mod1)-2)]))
  colnames(general_top_gen)<-c("times_found","gene")
  general_top_gen<-merge(general_top_gen,Mercator_Rirregularis,by="gene")
  general_top_gen
  write.table(general_top_gen,paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/TYPE_TOP_GENES/GENERALIST_TOP_GENES_info',i,'.txt') ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  specialized_top_gen<-cbind.data.frame(
    as.vector(table(unlist(top_mod1))[table(unlist(top_mod1))==1]),
    names(table(unlist(top_mod1))[table(unlist(top_mod1))==1]))
  colnames(specialized_top_gen)<-c("times_found","gene")
  specialized_top_gen<-merge(specialized_top_gen,Mercator_Rirregularis,by="gene")
  specialized_top_gen
  write.table(specialized_top_gen,paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/TYPE_TOP_GENES/SPECIALIZED_TOP_GENES_info',i,'.txt') ,quote=FALSE,sep='\t',col.names = T, row.names = T)
  
  
  }



all_top_genes_Y<-list()
  for (i in unique(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)," "), "[[", 2) ) ) {
    pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/TYPE_TOP_GENES/TYPE_TOP_GENES',i,'.pdf',sep="_"),width=10,height=6)
    
        top_mod1<-top_gene_sign_intramod_conec_Y[grep(i,names(top_gene_sign_intramod_conec_Y))]
    length(top_mod1)
    tryCatch(pheatmap( t(as.matrix(table(unlist(top_mod1))))/length(top_mod1),cluster_rows = F,
                       main=paste(i,'nb. comparisons=',length(top_mod1)),cluster_cols = T,cellheight= 8),error=function(e) tryCatch(pheatmap( t(c((t(as.matrix(table(unlist(top_mod1))))/length(top_mod1))[1],(t(as.matrix(table(unlist(top_mod1))))/length(top_mod1))[2]*.9)),cluster_rows = F,
                        main=i,cluster_cols = T,cellheight= 8) ),error=function(e) plot(1,1,main=i) )
    all_top_genes_Y[[i]]<-colnames(t(as.matrix(table(unlist(top_mod1)))) ) 
    
    general_top_gen<-cbind.data.frame(
      as.vector(table(unlist(top_mod1))[table(unlist(top_mod1))>(length(top_mod1)-2)])/length(top_mod1),
      names(table(unlist(top_mod1))[table(unlist(top_mod1))>(length(top_mod1)-2)]))
    tryCatch(colnames(general_top_gen)<-c("times_found","gene"),error=function(e) 1)
    tryCatch( general_top_gen<-merge(general_top_gen,Mercator_Mesculenta,by="gene"),error=function(e) 1)
    write.table(general_top_gen,paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/TYPE_TOP_GENES/GENERALIST_TOP_GENES_info',i,'.txt') ,quote=FALSE,sep='\t',col.names = T, row.names = T)
    general_top_gen
    
    specialized_top_gen<-cbind.data.frame(
      as.vector(table(unlist(top_mod1))[table(unlist(top_mod1))==1]),
      names(table(unlist(top_mod1))[table(unlist(top_mod1))==1]))
    tryCatch(colnames(specialized_top_gen)<-c("times_found","gene"),error=function(e) 1)
    tryCatch(specialized_top_gen<-merge(specialized_top_gen,Mercator_Mesculenta,by="gene"),error=function(e) 1)
    write.table(specialized_top_gen,paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/TYPE_TOP_GENES/SPECIALIZED_TOP_GENES_info',i,'.txt') ,quote=FALSE,sep='\t',col.names = T, row.names = T)
    
    dev.off()
  }

all_top_genes_Y[1]
all_top_genes_R[1]


SAMPLE_A<-2
SAMPLE_Y<-10
names(all_top_genes_R[SAMPLE_A])
names(all_top_genes_Y[SAMPLE_Y])

  top_AMF_gene<-all_top_genes_R[[SAMPLE_A]]
  top_YUCA_gene<-all_top_genes_Y[[SAMPLE_Y]]
  
  NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in%  top_AMF_gene,])
  
  NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  top_YUCA_gene,])
  NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

  correlacion<-tryCatch(cor(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA),error=function(e) 1)
  
  colsc=c(rgb(241, 54, 23, maxColorValue=255),  "white" , rgb(0, 61, 104, maxColorValue=255))
  colramp = colorRampPalette(colsc, space= 'Lab' )
  colors = colramp(100)
  bk2 = c(seq(-1,0,length=50),seq(0.02, 1, length=49))
  # function
  #https://hlplab.wordpress.com/2012/03/20/correlation-plot-matrices-using-the-ellipse-library/
  tryCatch(labo<-data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA),error=function(e) 1)
  tryCatch(pheatmap(correlacion,col=colramp(100),breaks = bk2,main=paste("Modules",names(all_top_genes_R[SAMPLE_A]),names(all_top_genes_R[SAMPLE_Y]))),error=function(e) plot(1,1, main=names(top_gene_sign_intramod_conec_A)[i]) )
  #fontsize = 15,cellwidth = 15,cellheight = 15,
  
  
  ################################################################################################################################
  # find pathogen, symbiosi and fred genes in modules
  Fred_YUCA_genes
  pathogen_genes
  symbiosis_genes
  
  top_gene_sign_intramod_conec_A
  unlist(top_gene_sign_intramod_conec_Y)[unlist(top_gene_sign_intramod_conec_Y) %in% Fred_YUCA_genes]
  unlist(top_gene_sign_intramod_conec_Y)[unlist(top_gene_sign_intramod_conec_Y) %in% pathogen_genes]
  unlist(top_gene_sign_intramod_conec_Y)[unlist(top_gene_sign_intramod_conec_Y) %in% symbiosis_genes]
  
  ################################################################################################################################
  # do network of correlation between genes...
 length(table(unlist(top_gene_sign_intramod_conec_Y)))
 length(table(unlist(top_gene_sign_intramod_conec_A)))
 
 # module 37 correction
 top_gene_sign_intramod_conec_Y[[37]]<- top_gene_sign_intramod_conec_Y[[37]][grep("Manes",top_gene_sign_intramod_conec_Y[[37]]) ]
 #
 
 as.vector(table(unlist(top_gene_sign_intramod_conec_Y)))
as.vector(table(unlist(top_gene_sign_intramod_conec_A)))
names(table(unlist(top_gene_sign_intramod_conec_Y)))
names(table(unlist(top_gene_sign_intramod_conec_A)))

unlist(top_gene_sign_intramod_conec_Y[1])
unlist(top_gene_sign_intramod_conec_A[1])

nod_Y<-cbind.data.frame(names(table(unlist(top_gene_sign_intramod_conec_Y))),
as.vector(table(unlist(top_gene_sign_intramod_conec_Y))),
rep("Cassava",length(table(unlist(top_gene_sign_intramod_conec_Y)))),
rep(1,length(table(unlist(top_gene_sign_intramod_conec_Y))))
)
colnames(nod_Y)<-c("id","size","org","org.type")
nod_Y<-nod_Y[grep("Manes",nod_Y[,1]),]

nod_A<-cbind.data.frame(names(table(unlist(top_gene_sign_intramod_conec_A))),
as.vector(table(unlist(top_gene_sign_intramod_conec_A))),
rep("AMF",length(table(unlist(top_gene_sign_intramod_conec_A)))),
rep(2,length(table(unlist(top_gene_sign_intramod_conec_A))))
)
colnames(nod_A)<-c("id","size","org","org.type")

nod_vf<-rbind.data.frame(nod_Y,nod_A)
nod_vf[grep()]

lin2<-list()
lin3<-NULL
for (i in 1:length(top_gene_sign_intramod_conec_Y)) {
lin2[[i]]<-expand.grid(unlist(top_gene_sign_intramod_conec_Y[i]),unlist(top_gene_sign_intramod_conec_A[i]))
lin3<-rbind.data.frame(lin3,lin2[[i]])
}

lin_vf<-cbind.data.frame(lin3,rep(1,dim(lin3)[1]))
colnames(lin_vf) <-c("from","to","weight") 
lin_vf<-lin_vf[grep("Manes",lin_vf[,1]),]
library(igraph) 
head(lin_vf)
lin_vf[lin_vf[,3]!=1]

lin_
tail(nod_vf)
top_gene_sign_intramod_conec_Y[[1]]
top_gene_sign_intramod_conec_A[[1]]

  
  net <- graph_from_data_frame(d=lin_vf, vertices=nod_vf, directed=F) 
  
  # Generate colors based on media type:
  V(net)$color <- c(rep("grey",length(table(merge$colors))),rep("brown",length(table(merge_AMF$colors))) )
  
  
  # Set node size based on audience size:
  V(net)$size <- V(net)$size
  V(net)$label.cex = 0
  pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/genes_network_temp.pdf',width= 16, height=16)
  l <- layout_on_sphere(net)
  plot(net,vertex.label.dist=-0.4,vertex.label=NA,layout=layout_with_kk)
  dev.off()
  
  
  

  ####
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  
  
  
  
  
  
  
  ################################################################################################################################
  # SAME ANALYSIS BUT WITH MODULES CORRELATED TO COLONIZATION INTENSITY
  ################################################################################################################################
  ###### RELOAD CASSAVA CORRELATION TO COLONIZATION INTESITY
  
  
  # CASSAVA
  
  sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,]),]
  #take out controls of CASSAVA
  sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
  top_gene_sign_intramod_conec_Y<-list()
  for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_SENSITIVE))  {
    names(corModules_S)[1]
    gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
    MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)]#gsub("A_","",
    MEs_AMF$A_MEred
    ######################################################
    # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
    ######################################################
    #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
    geneModuleMembership_YUCA=as.data.frame(cor(t(sign_norm_genes),MEs_Y,use="p"))
    geneTraitSignificance_YUCA= as.data.frame(cor(t(sign_norm_genes),as.numeric(COV_10_samples3[,1]),use="p")) # define correlation with external variable gsub("A_","",
    ######################################################
    # 1) MODULE DEFINITION
    ######################################################
    #Cassava
    ## CHANGE FOR THE COMPARISON DESIRED
    modNames=substring(names(MEs_Y),5)
    module=gsub(" Y_ME","",correlated_Modula)
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
    gene_2_go<-sign_feat_go[,c(1,9)]
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
    GS.6_Y<-geneModuleMembership_YUCA[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0 & abs(geneModuleMembership_YUCA[moduleGenes, match(gsub(" ","",correlated_Modula), names(geneModuleMembership_YUCA))])>0,]
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
    
    sampleGOdata <-  tryCatch(new("topGOdata", ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
    
    
    
    sampleGOdata
    
    
    resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
    resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
    resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
    
    allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                               classicKS = resultKS, elimKS = resultKS.elim,
                               orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
    
    RES_Y<-allRes
    
    
    write.table(RES_Y,    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Correlation_MODULES_COL_INTENSITY/GO_ENRICH_Cassava_col_intensity',correlated_Modula ,'GO_ENRICH.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
    
    pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Correlation_MODULES_COL_INTENSITY/BARPLOT_Cassava_col_intensity',correlated_Modula ,'.pdf',sep="_"), width=10, height=7)
    par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
    tryCatch( barplot(RES_Y$Significant[RES_Y$classicFisher<0.05], names = paste(RES_Y$GO.ID[RES_Y$classicFisher<0.05],RES_Y$Term[RES_Y$classicFisher<0.05],sep=" "),
                      xlab = "Nb. of significant genes",horiz=T,las=1),error=function(e) plot(1,1,main=correlated_Modula) )
    dev.off()
    
    
    to_pheat_Y<-RES_Y[RES_Y$classicFisher<0.05,c(4)]
    names(to_pheat_Y)<-paste(RES_Y[RES_Y$classicFisher<0.05,c(1)],RES_Y[RES_Y$classicFisher<0.05,c(2)],sep=" ")
    
    pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Correlation_MODULES_COL_INTENSITY/PHEAT_Cassava_col_intensity',correlated_Modula,'.pdf',sep="_"), width=18, height=10)
    tryCatch(pheatmap(t(to_pheat_Y),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10),error=function(e) plot(1,1,main=correlated_Modula))
    dev.off()
  }
  
  ######################################################  ######################################################
  # IN AMF
  ######################################################  ######################################################
  names(corModules_S)[c(1,10,13,16,19,24,26,27,28,34,35)]
  unique(sapply(strsplit(names(corModules_S)," "), "[[", 1))
  
  
  sign_norm_genes<-DGE_AMF_N$E
  #take out controls of CASSAVA
  sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
  top_gene_sign_intramod_conec_A<-list()
  for (correlated_Modula in   names(Cor_MOD_AMF_2_SENSITIVE))  {
    names(corModules_S)[1]
    gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
    MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)]
    
    
    ######################################################
    # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
    ######################################################
    #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
    geneModuleMembership_AMF=as.data.frame(cor(t(sign_norm_genes),MEs_AMF,use="p"))
    geneTraitSignificance_AMF= as.data.frame(cor(t(sign_norm_genes),
                                                 as.numeric(COV_10_samples3[,1]),use="p")) # define correlation with external variable
    ######################################################
    # 1) MODULE DEFINITION
    ######################################################
    #AMF
    modNames=substring(names(MEs_AMF),5)
    module=gsub(" A_ME","",correlated_Modula)
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
    gene_2_go<-sign_feat_go[,c(1,8)]
    
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
    GS.6_A<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_AMF[moduleGenes,1])>0 & abs(geneModuleMembership_AMF[moduleGenes,match(gsub(" ","",correlated_Modula), names(geneModuleMembership_AMF))])>0,]
    dim(GS.6_A)
    In_module<-rownames(GS.6_A) # genes name in module AMF
    
    ######################################################
    ######################################################
    # 5) ANNOTATION
    ######################################################
    #AFTER SELECTING FOR DIFFERENT ORGANISMS NOW DO SCRIPT FOR ANNOTATION
    
    geneList <- factor(as.integer(sign_feat_go$locusName %in% In_module))
    names(geneList) <- sign_feat_go$locusName
    
    sampleGOdata <-  tryCatch(new("topGOdata", ontology = "BP", allGenes = geneList , annot = annFUN.gene2GO, gene2GO = geneID2GO),error=function(e) 1)
    
    
    
    sampleGOdata
    
    
    resultFisher <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "fisher"),error=function(e) 1)
    resultKS <- tryCatch(runTest(sampleGOdata, algorithm = "classic", statistic = "ks"),error=function(e) 1)
    resultKS.elim <- tryCatch(runTest(sampleGOdata, algorithm = "elim", statistic = "ks"),error=function(e) 1)
    
    allRes<- tryCatch(GenTable(sampleGOdata, classicFisher = resultFisher,
                               classicKS = resultKS, elimKS = resultKS.elim,
                               orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 200),error=function(e) 1)
    
    
    RES_A<-allRes
    
    
    ######################################################
    # 6) PRINT RESULTS
    ######################################################
    
    
    ################
    
    
    
    
    write.table(RES_A,    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/GO_ENRICH_AMF_col_intensity',correlated_Modula ,'GO_ENRICH.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = F)
    
    pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/BARPLOT_AMF_col_intensity',correlated_Modula ,'.pdf',sep="_"), width=10, height=7)
    par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
    tryCatch( barplot(RES_A$Significant[RES_A$classicFisher<0.05], names = paste(RES_A$GO.ID[RES_A$classicFisher<0.05],RES_A$Term[RES_A$classicFisher<0.05],sep=" "),
                      xlab = "Nb. of significant genes",horiz=T,las=1),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 1)) )
    dev.off()
    
    
    to_pheat_A<-RES_A[RES_A$classicFisher<0.05,c(4)]
    names(to_pheat_A)<-paste(RES_A[RES_A$classicFisher<0.05,c(1)],RES_A[RES_A$classicFisher<0.05,c(2)],sep=" ")
    
    pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/PHEAT_AMF_col_intensity',correlated_Modula, '.pdf',sep="_"), width=18, height=10)
    tryCatch(pheatmap(t(to_pheat_A),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 1)))
    dev.off()
    
    
  }
  
  top_gene_sign_intramod_conec_A
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ################################################################################################################################
  ######################################################
  #' 1) LOAD MERCATOR DATABASE
  ######################################################
  # CASSAVA
  Mercator_Mesculenta<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Mesculenta_database_v2",h=F,quote="")
  head(Mercator_Mesculenta)
  colnames(Mercator_Mesculenta)<-c("Bin","Function","gene")
  # AMF
  Mercator_Rirregularis<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Rirregularis_database2",h=F,quote="")
  head(Mercator_Rirregularis)
  colnames(Mercator_Rirregularis)<-c("Bin","Function","gene")
  
  # 
  # extract more representative and more significant genes, no go terms
  #fastest
  
  
  
  sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,]),]
  #take out controls of CASSAVA
  sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
  top_gene_intramod_COLINTENSITY_conec_Y<-list()
  Import_gene_COLINTENSITY_Y<-list()
  for (correlated_Modula in   names(Cor_MOD_CASSAVA_2_SENSITIVE))  {
    names(corModules_S)[1]
    gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
    MEs_AMF[,names(MEs_AMF)==sapply(strsplit(correlated_Modula," "), "[[", 1)]#gsub("A_","",
    MEs_AMF$A_MEred
    ######################################################
    # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
    ######################################################
    #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
    geneModuleMembership_YUCA=as.data.frame(cor(t(sign_norm_genes),MEs_Y,use="p"))
    geneTraitSignificance_YUCA= as.data.frame(cor(t(sign_norm_genes),as.numeric(COV_10_samples3[,1]),use="p")) # define correlation with external variable gsub("A_","",
    ######################################################
    # 1) MODULE DEFINITION
    ######################################################
    #Cassava
    ## CHANGE FOR THE COMPARISON DESIRED
    modNames=substring(names(MEs_Y),5)
    module=gsub(" Y_ME","",correlated_Modula)
    column= match(module,modNames)
    moduleGenes=moduleColors_Y==module
   
    ######################################################
    # 4) SELECT TOP GENES
    ######################################################
    ## CHANGE SELECTED LOCI IN MODULE BY FILTER GS and MM
    ## apply filter of 0.6 to data gene sign and module membership
    ### CAHNGE MATCH VALUE FOR DIFFERENT COMPARISONS
    #CASSAVA
    GS.6_Y<-geneModuleMembership_YUCA[moduleGenes,][abs(geneTraitSignificance_YUCA[moduleGenes,1])>0 & abs(geneModuleMembership_YUCA[moduleGenes, match(gsub(" ","",correlated_Modula), names(geneModuleMembership_YUCA))])>0,]
    dim(GS.6_Y)
    
    In_module<-rownames(GS.6_Y) # genes name in module CASSAVA
    
    most_correlated_genes<-geneModuleMembership_YUCA[moduleGenes,]
    
    na.omit(most_correlated_genes[,grep(sapply(strsplit(correlated_Modula," "), "[[", 2), colnames(most_correlated_genes)  )])
    
    
    
    
    
    
    
    ######################################################
    # 6) PRINT RESULTS
    ######################################################
    ## intramodular conectivity
    ADJ1=abs(cor(t(sign_norm_genes),use="p"))^6
    Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors_Y)
    head(Alldegrees1)
    
    GS1=as.numeric(cor(as.numeric(COV_10_samples3[,1]),t(sign_norm_genes), use="p"))
    GeneSignificance=abs(GS1)
    
    colorlevels=unique(dynamicColors_Y)
    
    pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/ALL_MODULES_Y_col_intensity_',correlated_Modula,'.pdf',sep="_"), width=18, height=14)
    
    whichmodule=gsub(" Y_ME","",correlated_Modula);
    restrict1 = (dynamicColors_Y==whichmodule);
    verboseScatterplot(Alldegrees1$kWithin[restrict1],
                       GeneSignificance[restrict1], col=dynamicColors_Y[restrict1],
                       main=whichmodule,
                       xlab = paste("Connectivity",whichmodule), ylab = paste(whichmodule,"Gene Significance","\n to Colonization intensity",sep=" "), abline = TRUE)
    
    dev.off()
    
    
    
    datKME=signedKME(t(sign_norm_genes_Y), MEs_Y, outputColumnName="MM.")
    FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub(" Y_","MM.",correlated_Modula),colnames(datKME))]) >.8
    table(FilterGenes)
    top_gene_intramod_COLINTENSITY_conec_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]<-dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
    
    Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]<-sign_YUCA[rownames(sign_YUCA) %in%   as.vector(unlist(top_gene_intramod_COLINTENSITY_conec_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]])),]
    Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]<-cbind.data.frame(Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]],rownames(Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]))
    colnames(Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]])[dim(Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]])[2]]<-"gene"
    Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]<-merge(Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]],Mercator_Mesculenta,by="gene")
    
    
    write.table( Import_gene_COLINTENSITY_Y[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]],
                 paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/TOP_GENES_FUNCTION_COL_INTENSITY',correlated_Modula ,'COR_TO_COLINTENSITY','.txt',sep="_"),
                 quote=FALSE,sep='\t',col.names = T, row.names = T)
    
    ################
    
    MEList_YUCA=moduleEigengenes(t(sign_norm_genes),colors=dynamicColors_Y)
    exp_and_modules<-cbind.data.frame(sign_norm_genes,MEList_YUCA$validColors)
    colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
    
    
    # extract top 20 genes
    exp_and_modules[exp_and_modules$MODULE==gsub(" Y_ME","",correlated_Modula),]
    cor_to_module<-as.data.frame(cor(t(sign_norm_genes),
                                     as.numeric(COV_10_samples3[,1]),
                                     use="p"))
    top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
    colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
    
    top_correlated2<-top_correlated[ top_correlated$MODULE==gsub(" Y_ME","",correlated_Modula)  ,]
    
    
    #### all markers + p-value for MapMan
    COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
    colnames(COR_info_module)<-c("gene",paste("cor2","COL INTENSITY",sep="_"))
    write.table(  COR_info_module,
                  paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/GENES_2mapman',correlated_Modula ,'COR_TO','COLINTENSITY','.txt',sep="_")
                  ,quote=FALSE,sep='\t',col.names = T, row.names = F)
    
    
    #do barplot
    exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub(" Y_ME","",correlated_Modula)  ,]
    
    Genes_module<-t(sign_norm_genes)[,colnames(t(sign_norm_genes)) %in%  rownames(exp_mod2)]
    
    write.table(t(Genes_module),    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/ALL_GENES_COLINTENSITY',correlated_Modula,'GO_ENRICH.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = T)
    
  }
  
  unlist(lapply(top_gene_sign_intramod_conec_Y,length))
  table(mergedColors) #yuca modules

  # all important genes in all modules
  Interestingy<-sign_YUCA[rownames(sign_YUCA) %in%   as.vector(unlist(top_gene_sign_intramod_conec_Y)),]
  Interestingy<-cbind.data.frame(Interestingy,rownames(Interestingy))
  colnames(Interestingy)[dim(Interestingy)[2]]<-"gene"
  Interestingy<-merge(Interestingy,Mercator_Mesculenta,by="gene")
  
  
  ######################################################  ######################################################
  # IN AMF
  ######################################################  ######################################################

  
  sign_norm_genes<-DGE_AMF_N$E
  #take out controls of CASSAVA
  sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
  
  top_gene_intramod_COLINTENSITY_conec_A<-list()
  Import_gene_COLINTENSITY_A<-list()  
  
  for (correlated_Modula in   names(Cor_MOD_AMF_2_SENSITIVE))  {
    names(corModules_S)[1]
    gsub("A_","",sapply(strsplit(correlated_Modula," "), "[[", 1))
    MEs_Y[,names(MEs_Y)==sapply(strsplit(correlated_Modula," "), "[[", 2)]
    
    
    ######################################################
    # 0) LOAD MODULES DEFINITION, SAMPLES, AND MM GS values
    ######################################################
    #' ## DEFINE MODULES TO COMPARE, CHANGE FOR DIFFERENT COMPARISONS
    geneModuleMembership_AMF=as.data.frame(cor(t(sign_norm_genes),MEs_AMF,use="p"))
    geneTraitSignificance_AMF= as.data.frame(cor(t(sign_norm_genes),
                                                 as.numeric(COV_10_samples3[,1]),use="p")) # define correlation with external variable
    ######################################################
    # 1) MODULE DEFINITION
    ######################################################
    #AMF
    modNames=substring(names(MEs_AMF),5)
    module=gsub(" A_ME","",correlated_Modula)
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
    GS.6_A<-geneModuleMembership_AMF[moduleGenes,][abs(geneTraitSignificance_AMF[moduleGenes,1])>0 & abs(geneModuleMembership_AMF[moduleGenes,match(gsub(" ","",correlated_Modula), names(geneModuleMembership_AMF))])>0,]
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
    ADJ1=abs(cor(t(sign_norm_genes),use="p"))^6
    Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors_A)
    head(Alldegrees1)
    
    GS1=as.numeric(cor(as.numeric(COV_10_samples3[,1]),t(sign_norm_genes), use="p"))
    GeneSignificance=abs(GS1)
    
    colorlevels=unique(dynamicColors_A)
      pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/ALL_MODULES_A_col_intensity_',correlated_Modula,'.pdf',sep="_"), width=18, height=14)
    
    
    whichmodule=gsub(" A_ME","",correlated_Modula);
    restrict1 = (dynamicColors_A==whichmodule);
    verboseScatterplot(Alldegrees1$kWithin[restrict1],
                       GeneSignificance[restrict1], col=dynamicColors_A[restrict1],
                       main=whichmodule,
                       xlab = paste("Connectivity",whichmodule), ylab = paste(whichmodule,"Gene Significance","\n to Colonization intensity",sep=" "), abline = TRUE)
    
    dev.off()
    
    #datKME=signedKME(t(sign_norm_genes), MEs_AMF, outputColumnName="MM.")
    #FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub("A_","MM.",sapply(strsplit(correlated_Modula," "), "[[", 1)),colnames(datKME))])>.8
    #table(FilterGenes)
    #top_gene_sign_intramod_conec_A[[paste(sapply(strsplit(correlated_Modula," "), "[[", 1),"_to_",sapply(strsplit(correlated_Modula," "), "[[", 2),sep="")]]<-dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
    #top_gene_sign_intramod_conec_A
    
    datKME=signedKME(t(sign_norm_genes), MEs_AMF, outputColumnName="MM.")
    FilterGenes= abs(GS1)> .7 & abs( datKME[,grep(gsub(" A_","MM.",correlated_Modula),colnames(datKME))]) >.8
    table(FilterGenes)
    top_gene_intramod_COLINTENSITY_conec_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]<-dimnames(data.frame(t(sign_norm_genes)))[[2]][FilterGenes]
    
    Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]<-sign_AMF[rownames(sign_AMF) %in%   as.vector(unlist(top_gene_intramod_COLINTENSITY_conec_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]])),]
    Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]<-cbind.data.frame(Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]],rownames(Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]))
    colnames(Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]])[dim(Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]])[2]]<-"gene"
    Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]]<-merge(Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]],Mercator_Rirregularis,by="gene")
    
    
    write.table( Import_gene_COLINTENSITY_A[[paste(correlated_Modula,"_to_COL_INTENSITY",sep="")]],
                 paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/TOP_GENES_FUNCTION_COL_INTENSITY',correlated_Modula ,'COR_TO_COLINTENSITY','.txt',sep="_"),
                 quote=FALSE,sep='\t',col.names = T, row.names = T)
    
    
    

    ################
    
    MEList_AMF=moduleEigengenes(t(sign_norm_genes),colors=dynamicColors_A)
    exp_and_modules<-cbind.data.frame(sign_norm_genes,MEList_AMF$validColors)
    colnames(exp_and_modules)[dim(exp_and_modules)[2]]<-"MODULE"
    # extract top 20 genes
    exp_and_modules[exp_and_modules$MODULE==gsub(" A_ME","",correlated_Modula),]
    cor_to_module<-as.data.frame(cor(t(sign_norm_genes),
                                     as.numeric(COV_10_samples3[,1]),
                                     use="p"))
    top_correlated<-cbind.data.frame(exp_and_modules,cor_to_module)
    colnames(top_correlated)[dim(top_correlated)[2]]<-"CORRELADO"
    
    top_correlated2<-top_correlated[ top_correlated$MODULE==gsub(" A_ME","",correlated_Modula)  ,]
    
    
    #### all markers + p-value for MapMan
    COR_info_module<-cbind.data.frame(rownames(top_correlated2),top_correlated2[,dim(top_correlated2)[2]])
    colnames(COR_info_module)<-c("gene",paste("cor2","COL INTENSITY",sep="_"))
    write.table(  COR_info_module,
                  paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/GENES_2mapman',correlated_Modula ,'COR_TO','COLINTENSITY','.txt',sep="_")
                  ,quote=FALSE,sep='\t',col.names = T, row.names = F)
    
    
    
    #extract all genes in module
    exp_mod2<-exp_and_modules[ exp_and_modules$MODULE==gsub(" A_ME","",correlated_Modula) ,]
    
    Genes_module<-t(sign_norm_genes)[,colnames(t(sign_norm_genes)) %in%  rownames(exp_mod2)]
    
    write.table(t(Genes_module), paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/ALL_GENES_COLINTENSITY',correlated_Modula,'GO_ENRICH.txt',sep="_")
                ,quote=FALSE,sep='\t',col.names = T, row.names = T)
    
  }
  
  names(Import_gene_M_A)
  
  
  
  ################################################################################################################################
  ################################################################################################################################
  
  ################################################################################################################################
  # how many times genes are found in different modules?
  top_gene_intramod_COLINTENSITY_conec_A
  top_gene_intramod_COLINTENSITY_conec_Y

  length(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y)))
  length(table(unlist(top_gene_intramod_COLINTENSITY_conec_A)))
  
  table(unlist(top_gene_intramod_COLINTENSITY_conec_A))
  table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))
  
  pdf( paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/TOPGENES_COLONIZATIONSENSIVITY_hist.pdf',sep="_"),width=10,height=6)
  par(mfrow=c(1,2))
  hist(table(unlist(top_gene_intramod_COLINTENSITY_conec_A)),main="TOP AMF GENES IN MODULES",xlab="Number of modules found for gene")
  hist(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y)),main="TOP CASSAVA GENES IN MODULES",xlab="Number of modules found for gene")
  dev.off()
  
  #cosmopolitan genes
  cosmo_A<-table(unlist(top_gene_intramod_COLINTENSITY_conec_A))[table(unlist(top_gene_intramod_COLINTENSITY_conec_A))>2]
  cosmo_Y<-table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))[table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))>2]
  
  
  cosmo_A<-cbind.data.frame(cosmo_A)
  colnames(cosmo_A)[1]<-"gene"
  cosmo_Y<-cbind.data.frame(cosmo_Y)
  colnames(cosmo_Y)[1]<-"gene"
  Gene_mod_Mercator<-merge(cosmo_A,Mercator_Rirregularis,by="gene")
  Gene_mod_Mercator
  Gene_mod_Mercator<-merge(cosmo_Y,Mercator_Mesculenta,by="gene")
  Gene_mod_Mercator
  

  ##############################################################
  ##############################################################
  ##############################################################
  ##############################################################
  ##############################################################
  ##############################################################
  ##############################################################
  
  ### COMMON GENES CROSSTALKS TO COLONIZATION SENSITIVITY
  
  Cor_MOD_AMF_2_SENSITIVE
  Cor_MOD_CASSAVA_2_SENSITIVE
  Cor_MOD_CASSAVA_2_AMF
  
  
unique(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)," "), "[[", 1))
unique(sapply(strsplit(names(Cor_MOD_CASSAVA_2_AMF)," "), "[[", 2))


table(unlist(top_gene_intramod_COLINTENSITY_conec_A))
table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))


table(unlist(top_gene_sign_intramod_conec_A))
table(unlist(top_gene_sign_intramod_conec_Y))
## common AMF
table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_A))),names(table(unlist(top_gene_sign_intramod_conec_A)))))
table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_A))),names(table(unlist(top_gene_sign_intramod_conec_A)))))[table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_A))),names(table(unlist(top_gene_sign_intramod_conec_A))))) >1]
length(table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_A))),names(table(unlist(top_gene_sign_intramod_conec_A)))))[table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_A))),names(table(unlist(top_gene_sign_intramod_conec_A))))) >1])
common_A<-names(table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_A))),names(table(unlist(top_gene_sign_intramod_conec_A)))))[table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_A))),names(table(unlist(top_gene_sign_intramod_conec_A))))) >1])
# function
common_A<-cbind.data.frame(common_A)
colnames(common_A)[1]<-"gene"

Gene_common_mod_Mercator<-merge(common_A,Mercator_Rirregularis,by="gene")
Gene_common_mod_Mercator



### COMMON CASSAVA
table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))),names(table(unlist(top_gene_sign_intramod_conec_Y)))))
table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))),names(table(unlist(top_gene_sign_intramod_conec_Y)))))[table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))),names(table(unlist(top_gene_sign_intramod_conec_Y))))) >1]
length(table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))),names(table(unlist(top_gene_sign_intramod_conec_Y)))))[table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))),names(table(unlist(top_gene_sign_intramod_conec_Y))))) >1])
common_Y<-names(table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))),names(table(unlist(top_gene_sign_intramod_conec_Y)))))[table(c(names(table(unlist(top_gene_intramod_COLINTENSITY_conec_Y))),names(table(unlist(top_gene_sign_intramod_conec_Y))))) >1])
#function

common_Y<-cbind.data.frame(common_Y)
colnames(common_Y)[1]<-"gene"
Gene_common_mod_Mercator<-merge(common_Y,Mercator_Mesculenta,by="gene")
unique(Gene_common_mod_Mercator$Function)
