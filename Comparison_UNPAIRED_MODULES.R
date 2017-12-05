################################################################################################################################
# GO ANALYSIS IN MODULES
################################################################################################################################

sign_norm_genes<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,]),]
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
top_gene_sign_intramod_conec_Y<-list()
for (correlated_Modula in   names(corModules2)[c(459,471,434,472,462,385,382)])  {
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
  
  
  write.table(RES_Y,    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/UNPAIR_GO_ENRICH_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'GO_ENRICH.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  write.table(RES_Y[RES_Y[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/UNPAIR_GO_2REVIGO_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'GO_ENRICH.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = F, row.names = F)
  
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/UNPAIR_BARPLOT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'.pdf',sep="_"), width=10, height=7)
  par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  tryCatch( barplot(RES_Y$Significant[RES_Y$classicFisher<0.05], names = paste(RES_Y$GO.ID[RES_Y$classicFisher<0.05],RES_Y$Term[RES_Y$classicFisher<0.05],sep=" "),
                    xlab = "Nb. of significant genes",horiz=T,las=1),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 2)) )
  dev.off()
  
  
  to_pheat_Y<-RES_Y[RES_Y$classicFisher<0.05,c(4)]
  names(to_pheat_Y)<-paste(RES_Y[RES_Y$classicFisher<0.05,c(1)],RES_Y[RES_Y$classicFisher<0.05,c(2)],sep=" ")
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/UNPAIR_PHEAT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 2) ,'.pdf',sep="_"), width=18, height=10)
  tryCatch(pheatmap(t(to_pheat_Y),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 2)))
  dev.off()
}

###############################

######################################################  ######################################################
# IN AMF
######################################################  ######################################################
names(corModules_S)[c(1,10,13,16,19,24,26,27,28,34,35)]
unique(sapply(strsplit(names(corModules2)," "), "[[", 1))


sign_norm_genes<-DGE_AMF_N$E
#take out controls of CASSAVA
sign_norm_genes<-sign_norm_genes[,grep(c("CTRL"), colnames(sign_norm_genes),invert =T)]
top_gene_sign_intramod_conec_A<-list()
for (correlated_Modula in   names(corModules2)[c(289,459,135,113,217,1,73,435,55)])  {
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
  write.table(RES_A,    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/UNPAIR_GO_ENRICH_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'GO_ENRICH.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = T, row.names = F)
  write.table(RES_A[RES_A[,6]<0.05,c(1,6)],    paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/UNPAIR_GO_2REVIGO_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'GO_ENRICH.txt',sep="_")
              ,quote=FALSE,sep='\t',col.names = F, row.names = F)
  
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/UNPAIR_BARPLOT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.pdf',sep="_"), width=10, height=7)
  par(mar=c(5,25,3,1),las=1,mgp=c(3, 1, 0))
  tryCatch( barplot(RES_A$Significant[RES_A$classicFisher<0.05], names = paste(RES_A$GO.ID[RES_A$classicFisher<0.05],RES_A$Term[RES_A$classicFisher<0.05],sep=" "),
                    xlab = "Nb. of significant genes",horiz=T,las=1),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 1)) )
  dev.off()
  
  
  to_pheat_A<-RES_A[RES_A$classicFisher<0.05,c(4)]
  names(to_pheat_A)<-paste(RES_A[RES_A$classicFisher<0.05,c(1)],RES_A[RES_A$classicFisher<0.05,c(2)],sep=" ")
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/UNPAIR_PHEAT_CROSSTALKS',sapply(strsplit(correlated_Modula," "), "[[", 1) ,'.pdf',sep="_"), width=18, height=10)
  tryCatch(pheatmap(t(to_pheat_A),cluster_cols = T, cluster_rows=F, cellwidth=10, cellheight=10),error=function(e) plot(1,1,main=sapply(strsplit(correlated_Modula," "), "[[", 1)))
  dev.off()
  
  
}

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

write.table(Crosstalks_df_vf,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/4.UNPAIR_Crosstalks_df_vf.txt' ,quote=FALSE,sep='\t',col.names = T, row.names = F)
