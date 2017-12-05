
########################################################################################
#' Script to SEE BEHABIOUR OF GO TERMS IN SAMPLES
#' HIGHLIGTHED GO TERMS
########################################################################################
#' Needs already in memorory: annotations, DEG genes names and p-values, Universe genes. modules info.
#' 
#' Run AFTER RNA_seq_gene_ll_AMF_HOST_v6.R
#' RUN AFTER COMPARISON GO METHODS.R
#' ########################################################################################################################################
################################### correspondance to go. using  PHytozome annotation        ########################################### 

# Specific genes in Cassava and AMF
###################################################################################################################################################
DGE_YUCA_N$E

#Defense response GO:0006952
###################################################################################################################################################
#find genes list correspondant to this GO CASSAVA
# use normalized...
DEf_go<-mesculenta_go2[mesculenta_go2$GO=='GO:0006952',]

Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                               gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}

#in r irregularis
DEf_go<-rirregularis_go[rirregularis_go$GO=='GO:0006952',]
Normalised_sign<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  colnames(t(FOUR_VARS_AMF_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}

#Carbohydrate metabolic process GO:0044723,GO:005975 GO:0044724 GO:0016052 GO:1901137
###################################################################################################################################################
#find genes list correspondant to this GO CASSAVA
DEf_go<-mesculenta_go2[mesculenta_go2$GO=='GO:0044723',]
DEf_go2<-mesculenta_go2[mesculenta_go2$GO=='GO:0005975',]
DEf_go3<-mesculenta_go2[mesculenta_go2$GO=='GO:0044724',]
DEf_go4<-mesculenta_go2[mesculenta_go2$GO=='GO:0016052',]
DEf_go5<-mesculenta_go2[mesculenta_go2$GO=='GO:1901137',]

DEf_go<-rbind.data.frame(DEf_go,DEf_go2,DEf_go3,DEf_go4,DEf_go5)
Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}


#in r irregularis
DEf_go<-rirregularis_go[rirregularis_go$GO=='GO:0044723',]
DEf_go2<-rirregularis_go[rirregularis_go$GO=='GO:0005975',]
DEf_go3<-rirregularis_go[rirregularis_go$GO=='GO:0044724',]
DEf_go4<-rirregularis_go[rirregularis_go$GO=='GO:0016052',]
DEf_go5<-rirregularis_go[rirregularis_go$GO=='GO:1901137',]

DEf_go<-rbind.data.frame(DEf_go,DEf_go2,DEf_go3,DEf_go4,DEf_go5)

Normalised_sign<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  colnames(t(FOUR_VARS_AMF_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}


#ATPase activity, transmembrane mouvement GO:0042626
###################################################################################################################################################
#find genes list correspondant to this GO CASSAVA
DEf_go<-mesculenta_go2[mesculenta_go2$GO=='GO:0042626',]
Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}

#in r irregularis
DEf_go<-rirregularis_go[rirregularis_go$GO=='GO:0042626',]
Normalised_sign<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  colnames(t(FOUR_VARS_AMF_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}


#transport GO:0006810
###################################################################################################################################################
#find genes list correspondant to this GO CASSAVA
DEf_go<-mesculenta_go2[mesculenta_go2$GO=='GO:0050896',]
Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}

#in r irregularis
DEf_go<-rirregularis_go[rirregularis_go$GO=='GO:0009605',]
Normalised_sign<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  colnames(t(FOUR_VARS_AMF_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}


#transmembrane transport GO:0055085
###################################################################################################################################################
#find genes list correspondant to this GO CASSAVA
DEf_go<-mesculenta_go2[mesculenta_go2$GO=='GO:0055085',]
Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}

#in r irregularis
DEf_go<-rirregularis_go[rirregularis_go$GO=='GO:0055085',]
Normalised_sign<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  colnames(t(FOUR_VARS_AMF_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}



#transport GO:0006810
###################################################################################################################################################
#find genes list correspondant to this GO CASSAVA
DEf_go<-mesculenta_go2[mesculenta_go2$GO=='GO:0006810',]
Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}

#in r irregularis
DEf_go<-rirregularis_go[rirregularis_go$GO=='GO:0007584',]
Normalised_sign<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  colnames(t(FOUR_VARS_AMF_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)}


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
#ploting all genes in each sign. go in all modules

################
# CASSAVA
###############

setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1') 
filesToProcess <- dir(pattern = "^GO_ENRICH_CROSSTALKS_Y_")  #files to process.

for (module in filesToProcess[grep("CROSSTALKS_Y",filesToProcess)]) {
  gsub("_GO_ENRICH.txt","" ,  gsub("CROSSTALKS_","",module)   )
  all_go<-read.table(module,h=T,sep="\t")
# import Go in modules
#all_go<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_Y_MEgrey60_GO_ENRICH.txt",h=T,sep="\t")

sig_go<-all_go[all_go$classicFisher<0.05,]


pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/GENES_IN_GO',gsub("GO_ENRICH","" ,  gsub("CROSSTALKS_","",module)   ),".pdf",sep="_"
    ) , width=16, height=10)

for (GO in tryCatch( sig_go$GO.ID,error=function(e) 1) ){
  
  tryCatch( sig_go$Term[sig_go$GO.ID==GO],error=function(e) 1)
  
#find genes list correspondant to this GO CASSAVA
###################################################################################################################################################
DEf_go<-mesculenta_go2[mesculenta_go2$GO==GO,]
Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )


colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")

tryCatch( pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T, main=paste(gsub("_GO_ENRICH.txt","" ,  gsub("CROSSTALKS_","",module)   ),GO,"\n",sig_go$Term[sig_go$GO.ID==GO]),
         cellwidth=8, cellheight=8),error=function(e) 1 )
for(var in c("V1","V4","V5","V6","V8")) {
  df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
  tryCatch( pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,main=paste(gsub("_GO_ENRICH.txt","" ,  gsub("CROSSTALKS_","",module)   ),GO,"\n",sig_go$Term[sig_go$GO.ID==GO] ),
           cellwidth=8, cellheight=8),error=function(e)  1 )
           
           }
}
dev.off()

}



################
# AMF
###############

setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1') 
filesToProcess <- dir(pattern = "^GO_ENRICH_CROSSTALKS_A_")  #files to process.

for (module in filesToProcess[grep("CROSSTALKS_A",filesToProcess)]) {
  gsub("_GO_ENRICH.txt","" ,  gsub("CROSSTALKS_","",module)   )
  all_go<-read.table(module,h=T,sep="\t")
  # import Go in modules
  #all_go<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_Y_MEgrey60_GO_ENRICH.txt",h=T,sep="\t")
  
  sig_go<-all_go[all_go$classicFisher<0.05,]
  
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/GENES_IN_GO',gsub("GO_ENRICH","" ,  gsub("CROSSTALKS_","",module)   ),".pdf",sep="_"
  ) , width=14, height=10)

  for (GO in tryCatch( sig_go$GO.ID,error=function(e) 1) ) {
    
    tryCatch( sig_go$Term[sig_go$GO.ID==GO],error=function(e) 1)
    
    #find genes list correspondant to this GO CASSAVA
    ###################################################################################################################################################
    DEf_go<-rirregularis_go[rirregularis_go$GO==GO,]
    Normalised_sign<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  colnames(t(FOUR_VARS_AMF))]
    df_DEF<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% DEf_go$locusName ],
                              gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                              gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )
    
    
    colnames(df_DEF)[(dim(df_DEF)[2]-1):dim(df_DEF)[2]]<- c("treat","cultivar")
    
    tryCatch( pheatmap(t(df_DEF[,1:(dim(df_DEF)[2]-2)]),cluster_cols = T, cluster_rows=T, main=paste(gsub("_GO_ENRICH.txt","" ,  gsub("CROSSTALKS_","",module)   ),GO,"\n",sig_go$Term[sig_go$GO.ID==GO]),
                       cellwidth=8, cellheight=8),error=function(e)  1 )
    for(var in c("V1","V4","V5","V6","V8")) {
      df_DEF_v<-df_DEF[df_DEF$cultivar==var,]
      tryCatch( pheatmap(t(df_DEF_v[,1:(dim(df_DEF_v)[2]-2)]),cluster_cols = T, cluster_rows=T,main=paste(gsub("_GO_ENRICH.txt","" ,  gsub("CROSSTALKS_","",module)   ),GO,"\n",sig_go$Term[sig_go$GO.ID==GO] ),
                         cellwidth=8, cellheight=8),error=function(e)  1)
      
    }
  }
  dev.off()
  
}




#####################################################################################
#####################################################################################
### TOP GENES CORRELATION TO OTHER ORGANISM MODULE AND FuNCTION


#' 1) LOAD DATABASE
######################################################
# CASSAVA
Mercator_Mesculenta<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Mesculenta_database_v2",h=F,quote="")
head(Mercator_Mesculenta)
colnames(Mercator_Mesculenta)<-c("Bin","Function","gene")
# AMF
Mercator_Rirregularis<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/Mercator_Rirregularis_database2",h=F,quote="")
head(Mercator_Rirregularis)
colnames(Mercator_Rirregularis)<-c("Bin","Function","gene")

######################################################
###CASSAVA
# find important genes in modules
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1') 
filesToProcess <- dir(pattern = "^GENES_2mapman_Y")  #files to process.

all_genes<-list()
symb_in_module_Y<-list()
path_in_module_Y<-list()
fred_in_module_Y<-list()
top_funct_in_module_Y<-list()
for (module in filesToProcess[grep("^GENES_2mapman_Y",filesToProcess)]) {
  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )  
  all_genes[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-read.table(module,h=T,sep="\t")
  all_genes[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-merge(all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]],Mercator_Mesculenta,by="gene")
  
  
  symb_in_module_Y[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,1] %in% symbiosis_genes,]
  path_in_module_Y[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,1] %in% pathogen_genes,]
  fred_in_module_Y[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,1] %in% Fred_YUCA_genes,]
  
  symb_in_module_Y[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-merge(symb_in_module_Y[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]],Mercator_Mesculenta,by="gene")
  path_in_module_Y[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-merge(path_in_module_Y[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]],Mercator_Mesculenta,by="gene")
  fred_in_module_Y[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-merge(fred_in_module_Y[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]],Mercator_Mesculenta,by="gene")
  top_funct_in_module_Y[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-rbind.data.frame(all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,2]<quantile(all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,2],.10),],
                                                                                                         all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,2]>quantile(all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,2],.90),])
  
}  
#################
###AMF
# find important genes in modules
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1') 
filesToProcess <- dir(pattern = "^GENES_2mapman_A")  #files to process.

all_genes<-list()
symb_in_module_A<-list()
path_in_module_A<-list()
fred_in_module_A<-list()
top_funct_in_module_A<-list()
for (module in filesToProcess[grep("^GENES_2mapman_A",filesToProcess)]) {
  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )  
  all_genes[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-read.table(module,h=T,sep="\t")
  all_genes[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-merge(all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]],Mercator_Rirregularis,by="gene")
  
  
  symb_in_module_A[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,1] %in% symbiosis_genes,]
  path_in_module_A[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,1] %in% pathogen_genes,]
  fred_in_module_A[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[ gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,1] %in% Fred_YUCA_genes,]
  
  symb_in_module_A[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-merge(symb_in_module_A[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]],Mercator_Rirregularis,by="gene")
  path_in_module_A[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-merge(path_in_module_A[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]],Mercator_Rirregularis,by="gene")
  fred_in_module_A[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-merge(fred_in_module_A[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]],Mercator_Rirregularis,by="gene")
  top_funct_in_module_A[[  gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   ) ]]<-rbind.data.frame(all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,2]<quantile(all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,2],.10),],
                                                                                                         all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,2]>quantile(all_genes[[gsub("_COR.txt","" ,  gsub("GENES_2mapman_","",module)   )]][,2],.90),])
  
}  
top_funct_in_module_A[[  1 ]][,1]
top_funct_in_module_A[  1 ]
top_funct_in_module_Y[  1 ]


cor()
Manes.02G117700

plot(
  DGE_NYUCA_NOCTRL[rownames(DGE_NYUCA_NOCTRL) %in% top_funct_in_module_Y[[  1 ]][1:5,1],],
  
  DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% top_funct_in_module_A[[  1 ]][50:55,1],]
)


  DGE_NYUCA_NOCTRL<-DGE_YUCA_N$E[,grep("CTRL",colnames(DGE_YUCA_N$E) ,invert=T)]
  plot(
  DGE_NYUCA_NOCTRL[rownames(DGE_NYUCA_NOCTRL) %in% "Manes.02G117700",],

DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% "g3406.t1",]
)
"g3406.t1"

to_Disc<-rbind.data.frame(
DGE_NYUCA_NOCTRL[rownames(DGE_NYUCA_NOCTRL) %in% top_funct_in_module_Y[[  1 ]][,1],],
DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% top_funct_in_module_A[[  1 ]][,1],]
)

labello<-c(
  rep("CASSAVA",dim(DGE_NYUCA_NOCTRL[rownames(DGE_NYUCA_NOCTRL) %in% top_funct_in_module_Y[[  1 ]][,1],])[1]),
      rep("AMF",dim(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% top_funct_in_module_A[[  1 ]][,1],])[1])
)

head(to_Disc)
# multivariate
library(ade4)
to_Disc
isol<-dat[,2]
pca1 <- dudi.pca(to_Disc, scale= TRUE, scannf = F)
pca1.dis <- discrimin(pca1,as.factor(labello) , scannf = F) #discriminant analysis
s.class(pca1.dis$li, c("bleu","red") )# nmot working
s.arrow(pca1.dis$li )
s.label(pca1$li )
s.corcircle(pca1$co )
scatter(pca1)
pca1$eig/sum(pca1$eig)

install.packages('FactoMineR')
library(FactoMineR)

res.pca = PCA(to_pca, scale.unit=TRUE, ncp=5, graph=T)
plot.PCA(res.pca, axes=c(1, 2), choix="ind")




#####################################################################################
#####################################################################################


######################



###################################################################################################################################################
###################################################################################################################################################
# gene expression modules
# Cassava
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1') 
filesToProcess <- dir(pattern = "^GENES_CROSSTALKS_Y")  #files to process.

filesToProcess [12]
for (module in filesToProcess[grep("CROSSTALKS_A",filesToProcess)]) {
  gsub("_GO_ENRICH.txt","" ,  gsub("CROSSTALKS_","",module)   )
  
  # import Go in modules
  all_go<-read.table(filesToProcess [14],h=T,sep="\t")
  # take significant in cassava
  
  sign_in_go<-all_go[rownames(all_go) %in% rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.01,]),]
  sig_go<-all_go[all_go$classicFisher<0.05,]
  
  pheatmap(sign_in_go)
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/GENES_IN_GO',gsub("GO_ENRICH","" ,  gsub("CROSSTALKS_","",module)   ),".pdf",sep="_"
  ) , width=14, height=10)
  


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
################# NOW USE OF IMPORTANT SELECTED GENES AND GO TERMS
#########??USE ON MAPMAN AND OTHER FEATURES


Interesting_genes<-read.table("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/interesting_genes_vf.txt",h=T,sep="\t")

head(Interesting_genes)
Interesting_genes$GENE

Interesting_exp_Y<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  Interesting_genes$GENE]
Interesting_exp_A<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  Interesting_genes$GENE]

tail(sign_YUCA)
Interesting_FC_Y<-sign_YUCA[rownames(sign_YUCA) %in%  Interesting_genes$GENE,]
Interesting_FC_A<-sign_AMF[rownames(sign_AMF) %in%  Interesting_genes$GENE,]
Interesting_FC_Y[,1:2]

# Write table to import into MAPMAN. needs FOLD CHANGE OF SIGNIFICANT LOCI
write.table(sign_YUCA[sign_YUCA$adj.P.Val<0.05,1:3],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/FC_DEG_Y.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T) # in cassava by AMF 
write.table(sign_AMF[sign_AMF$adj.P.Val<0.05,1:2],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/FC_DEG_A.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T) # in AMF
write.table(sign_AMF_P[sign_AMF_P$adj.P.Val<0.05,1:5],'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/FC_DEG_A_byplant.txt',quote=FALSE,sep='\t',
            col.names = T, row.names = T) # in AMF by plant


DEG_Y<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])]
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/DEG5_Cassava_by_AMF.pdf")
pheatmap(DEG_Y)
dev.off()

DEG_Yb<- DEG_Y[grep("CTRL",rownames(DEG_Y),invert=T),]
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/DEG5b_Cassava_by_AMF_noCTRL.pdf")
pheatmap(DEG_Yb)
dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/DEG6_AMF.pdf")
DEG_A<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  rownames(sign_AMF[sign_AMF$adj.P.Val<0.05,])]
pheatmap(DEG_A)
dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/DEG7_AMF_by_Cassava.pdf")
DEG_A_P<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  rownames(sign_AMF_P[sign_AMF_P$adj.P.Val<0.05,])]
pheatmap(DEG_A_P)
dev.off()
head(DEG_A_P[grep("CAN",rownames(DEG_A_P)),])

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/DEG8_AMF_CAN_by_Cassava.pdf")
DEG_A_P<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  rownames(sign_AMF_P[sign_AMF_P$adj.P.Val<0.05,])]
pheatmap(DEG_A_P[grep("CAN",rownames(DEG_A_P)),])
dev.off()

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/DEG9_AMF_B1_by_Cassava.pdf")
DEG_A_P<-t(DGE_AMF_N$E)[,colnames(t(DGE_AMF_N$E)) %in%  rownames(sign_AMF_P[sign_AMF_P$adj.P.Val<0.05,])]
pheatmap(DEG_A_P[grep("B1",rownames(DEG_A_P)),])
dev.off()

###################################################################################################################################################
###################################################################################################################################################
################# SIMBIOSIS GENES WITHIN SIGNIFICANT IN PLANT influenced by CAN B1


# from mapman, change name of files m->M g->G
symbiosis_genes<-c("Manes.02G110000","Manes.12G121400","Manes.12G029800","Manes.01G104400","Manes.17G122900",
                   "Manes.12G029800","Manes.17G122900","Manes.08G165600","Manes.08G165600","Manes.08G165600",
                   "Manes.08G168800","Manes.02G110000","Manes.09G119200","Manes.14G162000","Manes.17G122900",
                   "Manes.15G092500","Manes.15G092500","Manes.15G092500","Manes.15G102400","Manes.01G125200",
                   "Manes.14G174200","Manes.03G083200","Manes.03G083200","Manes.03G091900","Manes.17G122900",
                   "Manes.17G122900","Manes.12G029800","Manes.17G122900","Manes.04G029000","Manes.08G060800",
                   "Manes.02G057000","Manes.02G110000","Manes.13G117700","Manes.13G117700","Manes.09G004400",
                   "Manes.12G029800","Manes.06G151900","Manes.11G136800","Manes.11G136800","Manes.01G094800",
                   "Manes.01G102500","Manes.17G122900","Manes.05G114600","Manes.12G025300","Manes.12G025400",
                   "Manes.09G164400","Manes.09G164400","Manes.17G122900","Manes.12G122800","Manes.07G045600",
                   "Manes.07G045800","Manes.03G096300","Manes.17G122900","Manes.10G095000","Manes.10G094900",
                   "Manes.02G070600","Manes.02G070600","Manes.02G070600","Manes.05G132300","Manes.05G132300",
                   "Manes.13G135000","Manes.13G026000","Manes.13G026000","Manes.12G102700","Manes.02G110000",
                   "Manes.12G029800","Manes.11G042900","Manes.11G042200","Manes.17G122900","Manes.08G124400",
                   "Manes.02G167300","Manes.05G114600","Manes.02G110000","Manes.12G094800","Manes.18G084400",
                   "Manes.07G119900","Manes.17G122900","Manes.02G110000","Manes.02G110000","Manes.12G121400",
                   "Manes.12G029800","Manes.17G122900","Manes.10G079500")

unique(symbiosis_genes)

dim(sign_YUCA[sign_YUCA$adj.P.Val<0.05,][  rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])  %in% symbiosis_genes , ])

Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
symbiosis2<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% symbiosis_genes ],
                          gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                          gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )

colnames(symbiosis2)[(dim(symbiosis2)[2]-1):dim(symbiosis2)[2]]<- c("treat","cultivar")

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/fig1_a_symbiosisgenes.pdf",width=10,height=8)
pheatmap(t(symbiosis2[,1:(dim(symbiosis2)[2]-2)]),cluster_cols = T, cluster_rows=T,
         cellwidth=8, cellheight=8)
dev.off()

symbiosis3<- symbiosis2[grep("CTRL",rownames(symbiosis2),invert=T),]
  pheatmap(t(symbiosis3[,1:(dim(symbiosis3)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=8)
  for(var in c("V1","V4","V5","V6","V8")) {
    symbiosis4<-symbiosis3[symbiosis3$cultivar==var,]
    pheatmap(t(symbiosis4[,1:(dim(symbiosis4)[2]-2)]),cluster_cols = T, cluster_rows=T,
             cellwidth=8, cellheight=8)
    }
  
  symbiosis4<-sign_YUCA[rownames(sign_YUCA) %in%  symbiosis_genes,]
  
  pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/MAPMAN1_CASSAVA_SYMBIOSISGENES_FC_CANB1_vsCTRL.pdf",width=10,height=7)
  pheatmap(t(symbiosis4[,1:2]),cellwidth=8, cellheight=8,cluster_cols = T, cluster_rows=T)
  dev.off()
  
  ##################################################################################################################################################
  ###################################################################################################################################################
  ################# PATHOGEN DEFENSE GENES WITHIN SIGNIFICANT IN PLANT influenced by CAN B1
  
  
  # from mapman, change name of files m->M g->G
  
  
  pathogen_genes<-c("Manes.05G025500","Manes.05G025500","Manes.08G137700","Manes.02G042500","Manes.02G042400",
    "Manes.05G175800","Manes.13G033200","Manes.06G131200","Manes.06G131100","Manes.18G040400",
    "Manes.15G092200","Manes.07G114500","Manes.01G085200","Manes.14G039000","Manes.14G038900",
    "Manes.03G165100","Manes.17G048400","Manes.17G051600","Manes.03G019200","Manes.12G152600",
    "Manes.16G050900","Manes.14G156100","Manes.02G206500","Manes.02G206500","Manes.02G206800",
    "Manes.13G133200","Manes.18G141400","Manes.02G018000","Manes.05G132100","Manes.05G132100",
    "Manes.05G031200","Manes.05G031200","Manes.01G253700","Manes.04G100700","Manes.08G096900",
    "Manes.08G145400","Manes.02G013200","Manes.02G154400","Manes.02G049200","Manes.05G097400",
    "Manes.05G175800","Manes.05G026200","Manes.13G048400","Manes.13G048400","Manes.16G090900",
    "Manes.12G108900","Manes.12G047500","Manes.12G047500","Manes.06G071700","Manes.06G032500",
    "Manes.18G040400","Manes.18G069000","Manes.18G125700","Manes.07G101700","Manes.07G135300",
    "Manes.11G013700","Manes.11G013700","Manes.11G013700","Manes.11G035600","Manes.11G035600",
    "Manes.11G035600","Manes.01G046200","Manes.01G222600","Manes.01G189600","Manes.14G098000",
    "Manes.03G056100","Manes.03G147500","Manes.03G199300","Manes.03G098000","Manes.03G044100",
    "Manes.04G102600","Manes.02G189600","Manes.02G189600","Manes.09G112700","Manes.09G112700",
    "Manes.16G129000","Manes.16G128900","Manes.06G153600","Manes.18G098700","Manes.03G009500",
    "Manes.02G206500","Manes.02G206500","Manes.02G206800","Manes.04G099100","Manes.16G045200",
    "Manes.11G117200","Manes.01G074600","Manes.11G117200","Manes.11G117200","Manes.08G036300",
    "Manes.10G059800","Manes.02G018000","Manes.05G132100","Manes.05G132100","Manes.05G031200",
    "Manes.05G031200","Manes.01G253700")
  
  unique(pathogen_genes)
  
  dim(sign_YUCA[sign_YUCA$adj.P.Val<0.05,][  rownames(sign_YUCA[sign_YUCA$adj.P.Val<0.05,])  %in% pathogen_genes , ])
  
  
  Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
  pathogen2<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% pathogen_genes ],
                                gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                                gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )
  
  colnames(pathogen2)[(dim(pathogen2)[2]-1):dim(pathogen2)[2]]<- c("treat","cultivar")
  
  pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/MAPMAN2b_CASSAVA_PATHOGENGENES_EXPRESSION_CANB1_vsCTRL.pdf")
  pheatmap(t(pathogen2[,1:(dim(pathogen2)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=6)
  dev.off()
  
  
  pathogen3<- pathogen2[grep("CTRL",rownames(pathogen2),invert=T),]
  pheatmap(t(pathogen3[,1:(dim(pathogen3)[2]-2)]),cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=6)
  for(var in c("V1","V4","V5","V6","V8")) {
    pathogen4<-pathogen3[pathogen3$cultivar==var,]
    pheatmap(t(pathogen4[,1:(dim(pathogen4)[2]-2)]),cluster_cols = T, cluster_rows=T,
             cellwidth=8, cellheight=6)
  }
  
  pathogen4<-sign_YUCA[rownames(sign_YUCA) %in%  pathogen_genes,]
  
  pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/MAPMAN2_CASSAVA_PATHOGENGENES_FC_CANB1_vsCTRL.pdf")
  pheatmap(t(pathogen4[,1:2]),cellwidth=8, cellheight=8,cluster_cols = T, cluster_rows=T)
  dev.off()

  
  ##################################################################################################################################################
  ###################################################################################################################################################
  ################# PATHOGEN DEFENSE GENES WITHIN SIGNIFICANT IN PLANT influenced by CAN B1
  # from mapman, change name of files m->M g->G
  
  diff_CAN_B1_genes<-c("Manes.04G102600","Manes.07G031900","Manes.15G070800","Manes.15G070700","Manes.03G130800",
    "Manes.03G183500","Manes.12G013700","Manes.02G185300","Manes.02G185300","Manes.11G083600",
    "Manes.08G030800","Manes.08G154200","Manes.17G083100","Manes.08G066100","Manes.09G173100",
    "Manes.12G018000","Manes.14G081100","Manes.14G009400","Manes.18G065400","Manes.01G009100",
    "Manes.14G135100","Manes.14G135100","Manes.14G135100","Manes.s037200","Manes.03G117700",
    "Manes.18G078400","Manes.09G100800","Manes.06G094200","Manes.06G094200","Manes.03G100900",
    "Manes.05G184300","Manes.03G142400","Manes.03G201100","Manes.03G201100","Manes.s058900",
    "Manes.09G071800","Manes.09G071800","Manes.06G161500","Manes.06G161500","Manes.03G208900",
    "Manes.02G206300","Manes.02G206600","Manes.10G026800","Manes.10G026800","Manes.15G046300",
    "Manes.15G046300","Manes.14G080600","Manes.14G080600","Manes.18G095100","Manes.10G056000",
    "Manes.03G193100","Manes.02G105800","Manes.08G130900","Manes.09G086300","Manes.09G090900",
    "Manes.09G086100","Manes.12G039900","Manes.04G075800","Manes.04G075800","Manes.s015500",
    "Manes.s015500","Manes.09G151200","Manes.16G004100","Manes.01G047200","Manes.06G099100",
    "Manes.04G033900","Manes.04G026400","Manes.14G089300","Manes.12G155900","Manes.01G235800",
    "Manes.14G080600","Manes.02G206300","Manes.02G206600","Manes.10G026800","Manes.10G026800",
    "Manes.15G046300","Manes.15G046300","Manes.14G080600","Manes.14G080600","Manes.18G095100",
    "Manes.10G056000","Manes.03G193100","Manes.02G135300","Manes.03G142000","Manes.10G118400",
    "Manes.11G116400","Manes.16G017800","Manes.16G017800","Manes.16G014900","Manes.16G017700",
    "Manes.11G116400","Manes.14G105900","Manes.s040000","Manes.s040000","Manes.s040000",
    "Manes.16G017700","Manes.14G080600","Manes.08G055800","Manes.02G159200","Manes.05G196400",
    "Manes.05G196400","Manes.09G125300","Manes.11G078000","Manes.11G078000","Manes.05G058800",
    "Manes.06G130700","Manes.14G170500","Manes.17G018400","Manes.17G018400","Manes.17G018400",
    "Manes.15G012500","Manes.10G014400","Manes.03G189000")
  
  unique( diff_CAN_B1_genes)
  Normalised_sign<-t(DGE_YUCA_N$E)[,colnames(t(DGE_YUCA_N$E)) %in%  colnames(t(FOUR_VARS_YUCA_S))]
  diff_CAN_B12<- cbind.data.frame(Normalised_sign[,colnames(Normalised_sign) %in% diff_CAN_B1_genes ],
                               gsub("m2","",gsub( "V[1-8]_\\w+_","",rownames(Normalised_sign))  ),
                               gsub("m2","",gsub( "_\\w+_\\w+$","",rownames(Normalised_sign)))   )
  
  colnames(diff_CAN_B12)[(dim(diff_CAN_B12)[2]-1):dim(diff_CAN_B12)[2]]<- c("treat","cultivar")
  colnames(diff_CAN_B12)<-gsub("Manes","M",colnames(diff_CAN_B12))
  #pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/Fig1_b_CANB1_genes.pdf",width=10,height=7)
  pheatmap(diff_CAN_B12[,1:(dim(diff_CAN_B12)[2]-2)],cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=6)
  #dev.off()
  
  
  diff_CAN_B13<- diff_CAN_B12[grep("CTRL",rownames(diff_CAN_B12),invert=T),]
  pheatmap(diff_CAN_B13[,1:(dim(diff_CAN_B13)[2]-2)],cluster_cols = T, cluster_rows=T,
           cellwidth=8, cellheight=6)
  for(var in c("V1","V4","V5","V6","V8")) {
    diff_CAN_B14<-diff_CAN_B13[diff_CAN_B13$cultivar==var,]
    pheatmap(t(diff_CAN_B14[,1:(dim(diff_CAN_B14)[2]-2)]),cluster_cols = T, cluster_rows=T,
             cellwidth=8, cellheight=6)
  }
  diff_CAN_B14<-sign_YUCA[rownames(sign_YUCA) %in%  diff_CAN_B1_genes,]
  rownames(diff_CAN_B14)<-gsub("Manes","M",rownames(diff_CAN_B14))
  pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/Fig1_b_CANB1_genes.pdf",width=11,height=7)
  pheatmap(t(diff_CAN_B14[,1:2]),cellwidth=8, cellheight=8,cluster_cols = T, cluster_rows=T,scale ="none")
  dev.off()
  t(diff_CAN_B14[,1:2])
  
  
  
  ######################################################
  ####
  # FRED LIST OF IMPORTANT PLANT GENES IN THE SYMBIOSIS
  
  Manes.06G156700
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.06G156700",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.12G124500",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.01G123300",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.04G002400",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.12G124300",]
  
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.05G053600",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.08G154600",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.18G034100",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.16G125700",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.06G143100",]
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.06G143100",] # not SIG in B1
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.14G029600",]
  
  sign_YUCA[rownames(sign_YUCA) %in% "Manes.01G071000",]
  
  Fred_YUCA_genes<-c("Manes.06G156700","Manes.12G124500","Manes.01G123300",
                     "Manes.04G002400","Manes.12G124300","Manes.05G053600",
                     "Manes.08G154600","Manes.18G034100","Manes.16G125700",
                     "Manes.06G143100","Manes.06G143100","Manes.14G029600",
                     "Manes.01G071000","Manes.01G193000")
  Fred_YUCA_genes<-  c("Manes.06G156700","Manes.12G124500","Manes.12G124500","Manes.01G123300",
"Manes.04G002400","Manes.12G124300","Manes.05G053600","Manes.08G154600",
"Manes.18G034100","Manes.16G125700","Manes.06G143100","Manes.01G071000","Manes.01G193000")
  sign_YUCA[rownames(sign_YUCA) %in% Fred_YUCA_genes,]
  pheatmap(sign_YUCA[rownames(sign_YUCA) %in% Fred_YUCA_genes,1:2])
  pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/Fig1_a_symb_genes.pdf',height=10,width = 10)
  
  symb<-DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Fred_YUCA_genes,]
  rownames(symb)<-gsub("Manes","M",rownames(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Fred_YUCA_genes,]))
  pheatmap(symb,cellwidth =8 ,cellheight=8)
dev.off()
  