########################################################################################
#' Script to compare list of genes to Mercator database
#' 
########################################################################################
#' 
#' 
#' 
#######################################################################################

######################################################
#' 1) LOAD DATABASE
######################################################
# CASSAVA
Mercator_Mesculenta<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Mesculenta_database_v2",h=F,quote="\'")
colnames(Mercator_Mesculenta)<-c("Bin","Function","gene")
# with detailed info
Mercator_Mesculenta<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Mesculenta_database_v3",h=F,quote="\'")
colnames(Mercator_Mesculenta)<-c("Bin","Function","gene","details")
# AMF
Mercator_Rirregularis<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Rirregularis_database2",h=F,quote="")
colnames(Mercator_Rirregularis)<-c("Bin","Function","gene")

# with detailed info
Mercator_Rirregularis<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Rirregularis_database_v4",h=F,quote="\'")
colnames(Mercator_Rirregularis)<-c("Bin","Function","gene","details")




######################################################
#' 2) LOAD GENES in module
######################################################

Genes_in_module<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/GENES_2mapman_Y_MEpurple_COR.txt",h=T,quote="")
head(Genes_in_module)


######################################################
#' 3) LOAD GENES data
######################################################
# CASSAVA
YUCA_TOTAL_NORM_GENES<-DGE_YUCA_N$E
YUCA_TOTAL_NORM_GENES<-FOUR_VARS_YUCA
YUCA_TOTAL_NORM_GENES<-sign_YUCA

#AMF
AMF_TOTAL_NORM_GENES<-DGE_AMF_N$E
AMF_TOTAL_NORM_GENES<-FOUR_VARS_AMF
AMF_TOTAL_NORM_GENES<-sign_AMF

######################################################
#' 4) Merge Genes in Module and Mercator
######################################################
# cassava
Gene_mod_Mercator<-merge(Genes_in_module,Mercator_Mesculenta,by="gene")
head(Gene_mod_Mercator)

select_top_gen_mod_Mer<-Gene_mod_Mercator[Gene_mod_Mercator[,2]<quantile(Gene_mod_Mercator[,2],.10) |
                                            Gene_mod_Mercator[,2] >quantile(Gene_mod_Mercator[,2],.90)  ,]
head(select_top_gen_mod_Mer)
# show correlated genes in module
select_top_gen_mod_Mer[,c(1,2,4)]
select_top_gen_mod_Mer[,sort(select_top_gen_mod_Mer[,2])]

# AMF
Gene_mod_Mercator<-merge(Genes_in_module,Mercator_Rirregularis,by="gene")
head(Gene_mod_Mercator)

select_top_gen_mod_Mer<-Gene_mod_Mercator[Gene_mod_Mercator[,2]<quantile(Gene_mod_Mercator[,2],.10) |
                                            Gene_mod_Mercator[,2] >quantile(Gene_mod_Mercator[,2],.90)  ,]
head(select_top_gen_mod_Mer)
# show correlated genes in module
select_top_gen_mod_Mer[,c(1,2,4)]
select_top_gen_mod_Mer[,sort(select_top_gen_mod_Mer[,2])]


#### merge  with significant genes
# AMF
head(sign_AMF_P)

SIGN_AMF_MERCATOR<-cbind.data.frame(rownames(sign_AMF_P),sign_AMF_P)
colnames(SIGN_AMF_MERCATOR)[1]<-"gene"

SIGN_AMF_MERCATOR2<-merge(SIGN_AMF_MERCATOR,Mercator_Rirregularis,by="gene")

head(SIGN_AMF_MERCATOR3)
SIGN_AMF_MERCATOR3<-SIGN_AMF_MERCATOR2[SIGN_AMF_MERCATOR2$Function!="not_assigned.unknown",]
rownames(SIGN_AMF_MERCATOR3)<-paste(SIGN_AMF_MERCATOR3[,1],SIGN_AMF_MERCATOR3[,11],sep=" ")
dim(SIGN_AMF_MERCATOR3[SIGN_AMF_MERCATOR3$adj.P.Val<0.05,2:5])
pheatmap(t(SIGN_AMF_MERCATOR3[SIGN_AMF_MERCATOR3$adj.P.Val<0.05,2:5]),cellheight = 8)
######################################################
#' 5) Plot important gene in module DEG_N
######################################################
# CASSAVA
rownames(YUCA_TOTAL_NORM_GENES)
Select_genes_mod_mercator<-YUCA_TOTAL_NORM_GENES[ rownames(YUCA_TOTAL_NORM_GENES) %in% select_top_gen_mod_Mer[,1],]
head(Select_genes_mod_mercator)
Select_genes_mod_mercator2<-Select_genes_mod_mercator[,    c(grep("CTRL",colnames(Select_genes_mod_mercator) ),grep("CAN",colnames(Select_genes_mod_mercator) ), grep("B1",colnames(Select_genes_mod_mercator) ) )]
Select_genes_mod_mercator3<-Select_genes_mod_mercator2[,c(grep("V1",colnames(Select_genes_mod_mercator2)),grep("V4",colnames(Select_genes_mod_mercator2)),grep("V5",colnames(Select_genes_mod_mercator2)),
                      grep("V6",colnames(Select_genes_mod_mercator2)),grep("V8",colnames(Select_genes_mod_mercator2)))]


pheatmap(Select_genes_mod_mercator2[,1:2],cellwidth = 8,cellheight = 8,main = "Top correlated genes to AMF module \n logFoldchange againts non-mycorrhizal plants")

# AMF
rownames(AMF_TOTAL_NORM_GENES)
Select_genes_mod_mercator<-AMF_TOTAL_NORM_GENES[ rownames(AMF_TOTAL_NORM_GENES) %in% select_top_gen_mod_Mer[,1],]
head(Select_genes_mod_mercator)
Select_genes_mod_mercator2<-Select_genes_mod_mercator[,    c(grep("CAN",colnames(Select_genes_mod_mercator) ), grep("B1",colnames(Select_genes_mod_mercator) ) )]
Select_genes_mod_mercator3<-Select_genes_mod_mercator2[,c(grep("V1",colnames(Select_genes_mod_mercator2)),grep("V4",colnames(Select_genes_mod_mercator2)),grep("V5",colnames(Select_genes_mod_mercator2)),
                                                          grep("V6",colnames(Select_genes_mod_mercator2)),grep("V8",colnames(Select_genes_mod_mercator2)))]


pheatmap(t(Select_genes_mod_mercator[,2:1]),cluster_rows =F ,cellwidth = 8,cellheight = 8,show_rownames = T,main = "Top correlated AMF genes to Cassava module \n B1 logFoldchange againts CAN")

