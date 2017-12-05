########################################################################################
#' Import subvcf and do analysis on them
########################################################################################
#source("http://bioconductor.org/biocLite.R")

options("scipen"=100,"digits"=2)
#install.packages('pegas')
#install.packages('pheatmap')
#install.packages("adegenet", dep=TRUE)

library("ape")
library("pegas")
library("seqinr")
library("ggplot2")

library("adegenet")
library(pegas)
library(pheatmap)
library(phangorn)
#install.packages("phangorn")
library(devtools)
#install_github("emmanuelparadis/pegas/pegas")
########################################################################################


########################################################################################
########################################################################################
# ANALYSIS ON R: IRREGULARIS
########################################################################################
########################################################################################

########################################################################################
#' Import simplified vcf files from AMF
########################################################################################
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/variant_data/R_irregularis_reads') 
filesToProcess <- dir(pattern = "simplified_VCF$")  #files to process.

#######################################################################################################
# ALLELE FREQUENCY distribution of AMF samples
#######################################################################################################

filesToProcess
for (filesT in filesToProcess[c(13,14,11)]) {
v1<-read.delim(filesT,h=T,sep="\t",fill=T)
v1_b<-v1[v1$General_coverage>10,]
v1_b<-v1_b[grep(",",v1_b$Occurrence_Alleles),]
v1_c<-strsplit(as.character(v1_b$Occurrence_Alleles),",")
v1_d<-lapply(v1_c, function (x)  unlist(x[sample(1:2,1,replace=T)]))
v1_e<-cbind.data.frame(v1_b,unlist(v1_d))
dim(v1_e)
colnames(v1_e)[14]<-"rand_allele"
v1_e$rand_allele<-as.numeric(levels(v1_e$rand_allele))[v1_e$rand_allele]

hist(v1_e$rand_allele/v1_e$General_coverage, breaks=100,col="darkgray",xlim=c(0,1),
     main=paste(gsub("_20160307.RNAsubVCF","",filesT),"Cov.",dim(v1_e)[1],sep="_"),
     xlab="Allele fq.",las=1)
}
##########################################################################################



##########################################################################################
### import and rename samples
##########################################################################################

samp<-list()
for (filesT in filesToProcess) {
 samp[[filesT]]<-read.delim(filesT,h=T,sep="\t",fill=T)
}
names(samp)
head(samp[[65]],100)

df <- as.data.frame(samp)
df_pre<-df
colnames(df_pre)

colnames(df_pre)<-gsub('\\.8_',"_8_CTRL_",colnames(df_pre))
colnames(df_pre)<-gsub('V88b_',"V8_8b_CTRL_",colnames(df_pre))

colnames(df_pre)<- gsub("\\.4_","_4_CTRL_",colnames(df_pre))
colnames(df_pre)<- gsub("\\.12_","_12_CTRL_",colnames(df_pre))
colnames(df_pre)<- gsub("\\.16_","_16_CTRL_",colnames(df_pre))
colnames(df_pre)<- gsub("\\.16b_","_16b_CTRL_",colnames(df_pre))

colnames(df_pre)<- gsub("\\.20_","_20_CTRL_",colnames(df_pre))

colnames(df_pre)<-gsub('\\.1_',"_1_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.1a_',"_1a_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.5_',"_5_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.9_',"_9_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.9b_',"_9b_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.9a_',"_9a_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.13_',"_13_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.13b_',"_13b_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.17_',"_17_CAN_",colnames(df_pre))

colnames(df_pre)<-gsub('\\.2_',"_2_B1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.2a_',"_2a_B1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.6_',"_6_B1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.10_',"_10_B1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.14_',"_14_B1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.14b_',"_14_B1_",colnames(df_pre))
colnames(df_pre)<-gsub('V810_',"V8_10_B1_",colnames(df_pre))

colnames(df_pre)<-gsub('\\.18_',"_18_B1_",colnames(df_pre))

colnames(df_pre)<-gsub('\\.3_',"_3_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.7_',"_7_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.11_',"_11_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.11_',"_11_CANB1_",colnames(df_pre))

colnames(df_pre)<-gsub('\\.15_',"_15_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.19_',"_19_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('\\.19A_',"_19A_CANB1_",colnames(df_pre))

colnames(df_pre)<-gsub('mV13_',"mV1_3_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('mV17_',"mV1_7_CANB1_",colnames(df_pre))

colnames(df_pre)<-gsub('m2V513_',"m2V5_13_CAN_",colnames(df_pre))
colnames(df_pre)<-gsub('m2V518_',"m2V5_18_B1_",colnames(df_pre))

colnames(df_pre)<-gsub('mV47_',"mV4_7_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('mV411_',"mV4_11_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('mV419_',"mV4_19_CANB1_",colnames(df_pre))

colnames(df_pre)<-gsub('mV615_',"mV6_15_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('mV611_',"mV6_11_CANB1_",colnames(df_pre))

colnames(df_pre)<-gsub('mV67_',"mV6_7_CANB1_",colnames(df_pre))

colnames(df_pre)<-gsub('mV811_',"mV8_11_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('mV819_',"mV8_19_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('mV87_',"mV8_7_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('V116_',"V1_16_CTRL_",colnames(df_pre))
colnames(df_pre)<-gsub('V120_',"V1_20_CTRL_",colnames(df_pre))
colnames(df_pre)<-gsub('V48_',"V4_8_CTRL_",colnames(df_pre))

colnames(df_pre)<-gsub('V56_',"V5_6_B1_",colnames(df_pre))

colnames(df_pre)<-gsub('V57_',"V5_7_CANB1_",colnames(df_pre))
colnames(df_pre)<-gsub('V58_',"V5_8_CTRL_",colnames(df_pre))


colnames(df_pre)<-gsub("_20160529.simplified_VCF","",colnames(df_pre))
colnames(df_pre)<-gsub("\\.20160307.RNAsubVCF","",colnames(df_pre))
colnames(df_pre)<-gsub("_20160307..1..RNAsubVCF","",colnames(df_pre))


df_n<-df_pre[,grep("\\.Alleles",colnames(df_pre))] # extract alleles
colnames(df_n)

#############################################################################################################################
##### COVERAGE ANALYISIS on R.irregularis
#############################################################################################################################

#############################################################################################################################
# Min coverage analysis on samples, cov >30
# Find how bad samples are sequenced.
# to identify controls and bad sequenced samples.
#  take info about general coverage

head(df_pre)

for (i in seq(0,1118,13)[-1]) {
df_pre[,i][df_pre[,i]>=30]<-rep("sup30",length( df_pre[,i][df_pre[,i]>=30]))
df_pre[,i][df_pre[,i]<30 & df_pre[,i]>=10]<-rep("from11to30",length( df_pre[,i][df_pre[,i]<30&df_pre[,i]>=10]))
df_pre[,i][df_pre[,i]==0  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==0  ]))
df_pre[,i][df_pre[,i]==1  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==1  ]))
df_pre[,i][df_pre[,i]==2  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==2  ]))
df_pre[,i][df_pre[,i]==3  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==3  ]))
df_pre[,i][df_pre[,i]==4  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==4  ]))
df_pre[,i][df_pre[,i]==5  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==5  ]))
df_pre[,i][df_pre[,i]==6  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==6  ]))
df_pre[,i][df_pre[,i]==7  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==7  ]))
df_pre[,i][df_pre[,i]==8  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==8  ]))
df_pre[,i][df_pre[,i]==9  ]<-rep("NO DATA",length( df_pre[,i][df_pre[,i]==9  ]))
}
head(df_pre[1:2,1:50])



df_sequen<-df_pre[,grep("\\.General_coverage",colnames(df_pre))] # extract coverage per sample

dim(df_sequen)
nodata<-list() 
for (i in 1:dim(df_sequen)[2]) {
nodata[[i]]<-length(df_sequen[,i][df_sequen[,i]=="NO DATA"])
}
nodata2<-unlist(nodata)

names(nodata2)<-gsub("\\.General_coverage","",colnames(df_sequen))

### aLL ISOLATES coverage
pheatmap(t(nodata2),cluster_rows = F,cluster_cols = T,cellwidth = 8,cellheight =8,show_colnames = T)



###############################################################################################
### only isolates in manuscript RNASEQ1
# take out samples
nodata2<-nodata2[grep(c("CANB1"), names(nodata2),invert =T)]

nodata2<-nodata2[grep(c("V3"), names(nodata2),invert =T)]

nodata2<-nodata2[grep(c("V1_4"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V1_12"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V1_8"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V4_16b"), names(nodata2),invert =T)]

nodata2<-nodata2[grep(c("V4_8"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V5_8"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V8_8b"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V8_8"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V8_6"), names(nodata2),invert =T)]

nodata2<-nodata2[grep(c("V5_9"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V5_4"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V5_18"), names(nodata2),invert =T)]
nodata2<-nodata2[grep(c("V5_6"), names(nodata2),invert =T)]


dim(df_sequen)[1]# total number of sites
dim(df_sequen)[1]-t(nodata2)

#sites with coverage more than 10
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/0.SNP1_cov10.pdf")
pheatmap(dim(df_sequen)[1]-t(nodata2),cluster_rows = F,cluster_cols = T,cellwidth = 8,cellheight =8,show_colnames = T)
dev.off()
as.vector(nodata2)
COV_10_samples<-dim(df_sequen)[1]-t(nodata2)
names(COV_10_samples)
write.table(COV_10_samples,'~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/CROSSTALKS_AMF_CASSAVA/COV_10_samples.txt',
            sep="\n",quote=F,row.names = F,col.names = F)
##########################################################################################

##### SNP ANALYSIS
##########################################################################################
# count how many NA in data.
# give info about bad sequenced  and bad controls

df_n<-df_pre[,grep("\\.Alleles",colnames(df_pre))] # extract alleles
colnames(df_pre)

head(df_n)
names(df_n)
nb.NA<-list()
for (sampleT in  1:length(filesToProcess) ) {
  name <- names(df_n)[sampleT]
nb.NA[[name]]<-table(df_n[,sampleT]=="na")
}
ff<-unlist(nb.NA)

NA.data<-cbind.data.frame(colnames(df_n),ff[grep("TRUE",names(ff))],rep(1,length(nb.NA)))
colnames(NA.data)<-c("sample","Nb.NA","seq")

rownames(NA.data)<-gsub("\\.Alleles\\.TRUE","",rownames(NA.data))
barplot(NA.data$Nb.NA,names=rownames(NA.data) ,las=2,)

pheatmap(t(NA.data[,2:3]),cluster_cols = T, cluster_rows=F,cellwidth=8,
         scale="row",cellheight=8)

##########################################################################################
# DO PCoA to see difference on samples
# take out data with a lot of NA

head(df_pre)

df2<-df_pre[,grep("CTRL",colnames(df_pre),invert=T)]

df2<-df2[,grep("V8_18_B1",colnames(df2),invert=T)]
df2<-df2[,grep("V3_19A",colnames(df2),invert=T)]
df2<-df2[,grep("V5_9a",colnames(df2),invert=T)]
df2<-df2[,grep("V5_18",colnames(df2),invert=T)]

#### not do coverage
head(df2[1:2,c(13,26)])
colnames(df2)
coverage <- 30
df2 <- df2[ which( df2[,13] > coverage | df2[,26] > coverage | df2[,39] > coverage | df2[,52] > coverage ) , ]
df2 <- df2[ which( df2[,65] > coverage | df2[,78] > coverage | df2[,91] > coverage | df2[,104] > coverage ) , ]
df2 <- df2[ which(  df2[,117] > coverage | df2[,130] > coverage | df2[,143] > coverage | df2[,156] > coverage ) , ]
df2 <- df2[ which(  df2[,169] > coverage | df2[,182] > coverage | df2[,195] > coverage | df2[,208] > coverage ) , ]
df2 <- df2[ which(  df2[,221] > coverage | df2[,234] > coverage | df2[,247] > coverage | df2[,260] > coverage ) , ]
df2 <- df2[ which(  df2[,273] > coverage | df2[,286] > coverage | df2[,299] > coverage | df2[,312] > coverage ) , ]
df2 <- df2[ which(  df2[,325] > coverage | df2[,338] > coverage | df2[,351] > coverage | df2[,364] > coverage ) , ]
df2 <- df2[ which(  df2[,377] > coverage | df2[,390] > coverage | df2[,403] > coverage | df2[,416] > coverage ) , ]
df2 <- df2[ which(  df2[,429] > coverage | df2[,442] > coverage | df2[,455] > coverage | df2[,468] > coverage ) , ]
df2 <- df2[ which(  df2[,481] > coverage | df2[,494] > coverage | df2[,507] > coverage ) , ]
df2 <- df2[ which(  df2[,520] > coverage | df2[,533] > coverage | df2[,546] > coverage | df2[,559] > coverage ) , ]
df2 <- df2[ which(  df2[,572] > coverage | df2[,585] > coverage | df2[,598] > coverage | df2[,611] > coverage ) , ]
df2 <- df2[ which(  df2[,624] > coverage | df2[,637] > coverage | df2[,650] > coverage | df2[,663] > coverage ) , ]
df2 <- df2[ which(  df2[,676] > coverage | df2[,689] > coverage | df2[,702] > coverage | df2[,715] > coverage ) , ]

df2 <- df2[ which(  df2[,728] > coverage ), ] #| df2[,690] > coverage | df2[,705] > coverage  ) , ]
####??continue
head(df2,5)


df3<-df2[,grep("\\.Alleles",colnames(df2))]
colnames(df3)<-gsub("\\.Alleles","",colnames(df3))
  for (cool in 1:dim(df3)[2]) {
  df3<-df3[grep("na",df3[,cool],invert=T),]
}
colnames(df3)<-gsub("_20160307.RNAsubVCF.Alleles","",colnames(df3)) ;colnames(df3)


df.1<-df2genind(t(df3),sep="\t")
df.2<-dist(df.1,method="euclidean")

attr(df.2, "Labels") <-  colnames(df3)
fit<-cmdscale(df.2, k = 2, eig = TRUE, add = FALSE, x.ret = FALSE)
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab=(fit$eig/sum(fit$eig))[1], ylab=(fit$eig/sum(fit$eig))[2],cex=2,
     type="p",pch=16,col=c(36,36,36,37,37,37,7,7,7,12,12,33,33,33,41),cex.axis=1.5, cex.lab=1.5)
legend(x=-2,y=2.8,c("sp1-1","sp1-2","sp1-3","sp1-32","sp1-4","sp1-5"),
       cex=0.9,pch=16, col=c(36,37,7,12,33,41))
fit$eig/sum(fit$eig)

tre <- nj(df.2)


plot(tre,type="fan", edge.w = 2)

plot(tre, edge.w = 2,cex=0.8)

##################
### SAMPLES USED IN MANUSCRIPT

# take out samples
df3<-df3[,grep(c("CANB1"), names(df3),invert =T)]
df3<-df3[,grep(c("V3"), names(df3),invert =T)]

df3<-df3[,grep(c("V1_4"), names(df3),invert =T)]
df3<-df3[,grep(c("V1_12"), names(df3),invert =T)]
df3<-df3[,grep(c("V1_8"), names(df3),invert =T)]
df3<-df3[,grep(c("V4_16b"), names(df3),invert =T)]

df3<-df3[,grep(c("V4_8"), names(df3),invert =T)]
df3<-df3[,grep(c("V5_8"), names(df3),invert =T)]
df3<-df3[,grep(c("V8_8b"), names(df3),invert =T)]
df3<-df3[,grep(c("V8_8"), names(df3),invert =T)]
df3<-df3[,grep(c("V8_6"), names(df3),invert =T)]

df3<-df3[,grep(c("V5_9"), names(df3),invert =T)]
df3<-df3[,grep(c("V5_4"), names(df3),invert =T)]
df3<-df3[,grep(c("V5_18"), names(df3),invert =T)]
df3<-df3[,grep(c("V5_6"), names(df3),invert =T)]

colnames(df3)

df.1<-df2genind(t(df3),sep="\t")
df.2<-dist(df.1,method="euclidean")
#df.2<-dist.ml(df.1)

df.1<-phyDat(t(df3), return.index=TRUE)
df.2<-dist.ml(df.1)

tre <- nj(df.2)
myBoots <- boot.phylo(tre, df.1@tab, function(e) root(nj(dist(e,method="euclidean")),1))

plot(tre, "unrooted", main="NJ")
parsimony(tre, df.1)
fit = pml(tre, data=df.1)
fitJC<- optim.pml(fit, TRUE)
bs = bootstrap.pml(fitJC, bs=200, optNni=TRUE,control = pml.control(trace = 0))
plotBS(midpoint(fitJC$tre), bs, p = 50, type="p")

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/SNP2_CAN_B1_SNP.pdf")
#plot(tre, main="NJ")

plotBS(midpoint(fitJC$tre), bs, p = 50, type="p")
#nodelabels(myBoots, cex=.6)
dev.off()
######################################################################################
######################################################################################
######################################################################################



########################################################################################
########################################################################################
# ANALYSIS ON CASSAVA
########################################################################################
########################################################################################
######################################################################################

########## VARIETIEs CLUSTERING.
#using simplified VCF.
#' Import vcf files
########################################################################################
setwd('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/variant_data/M_esculenta_reads') 

# make simplifiedVCF. 
#for i in $(ls *.simplified_VCF); do cat $i | awk '$12 >10' > $(echo $i | cut -d'.' -f1).filt_sVCF ;done



filesToProcess <- dir(pattern = ".filt_sVCF$")  #files to process.
filesToProcess


### CAHNGE FOR INFORMATIVE NAME
filesToProcess<-filesToProcess[grep("V3",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V57",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V13",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V17",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V411",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V419",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V47",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V611",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V615",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V67",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V811",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V819",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V87",filesToProcess,invert=T)]

filesToProcess<-filesToProcess[grep("V5-11",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V5-15",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V5-19",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V518",filesToProcess,invert=T)]
filesToProcess<-filesToProcess[grep("V5-3",filesToProcess,invert=T)]


filesToProcess[grep("V57",filesToProcess)]

filesToProcess


# take out samples not used in manuscript

filesToProcess<-filesToProcess[grep(c("V1-4"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V1-12"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V1-8"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V4-16b"), filesToProcess,invert =T)]

filesToProcess<-filesToProcess[grep(c("V48"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V58"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V88b"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V8-8"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V8-6"), filesToProcess,invert =T)]

filesToProcess<-filesToProcess[grep(c("V5-9"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V5-4"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V5-18"), filesToProcess,invert =T)]
filesToProcess<-filesToProcess[grep(c("V56"), filesToProcess,invert =T)]



samp_C<-list()
for (filesT in filesToProcess) {
  samp_C[[filesT]]<-read.delim(filesT,h=T,sep="\t",fill=T)
}

dim(samp_C3[[40]])

samp_C2<-lapply(samp_C, function (x) cbind.data.frame(x[,c(2,4,11,12,13)],paste(x$scaffold,x$position,sep="_")) )
samp_C3<-lapply(samp_C2, function (x) x[x$Alleles!="na",])
samp_C3<-lapply(samp_C3, setNames, nm = c("scaffold","position","Alleles","Occurrence_Alleles","General_coverage","ID"))


head(samp_C3[[2]])

temp<-merge(samp_C3[[1]],samp_C3[[2]][,3:6] , by="ID")
colnames(temp)<-gsub("\\.x",gsub("_20160529.filt_sVCF","",names(samp_C[1])),colnames(temp))
colnames(temp)<-gsub("\\.y",gsub("_20160529.filt_sVCF","",names(samp_C[2])),colnames(temp))

# select only 5 samples per variety


for  (fifi in c(3:length(samp_C3))) {
temp<-merge(temp,samp_C3[[fifi]][,3:6] , by="ID")
colnames(temp)[(length(colnames(temp))-2):length(colnames(temp))]<-paste(colnames(temp)[(length(colnames(temp))-2):length(colnames(temp))],  gsub("_20160514.filt_sVCF","",names(samp_C[fifi])) ,sep="")

}

temp2<-temp[,seq(4,length(colnames(temp)), by=3)]

colnames(temp2)<-gsub("_20160529.filt_sVCF","",colnames(temp2))


df.1<-df2genind(t(temp2),sep="\t")
df.2<-dist(df.1,method="euclidean")

tre <- nj(df.2)

attr(df.2, "Labels") <-  colnames(temp2)
myBoots <- boot.phylo(tre, df.1@tab, function(e) root(nj(dist(e,method="euclidean")),1))
myBoots
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/SNP3_CASSAVA.pdf")
plot(tre,edge.w = 2,cex=0.5)
nodelabels(myBoots, cex=.6)
dev.off()


df.1<-phyDat(t(df3), return.index=TRUE)
df.2<-dist.ml(df.1)

tre <- nj(df.2)
myBoots <- boot.phylo(tre, df.1@tab, function(e) root(nj(dist(e,method="euclidean")),1))

plot(tre, "unrooted", main="NJ")
parsimony(tre, df.1)
fit = pml(tre, data=df.1)
fitJC<- optim.pml(fit, TRUE)
bs = bootstrap.pml(fitJC, bs=200, optNni=TRUE,control = pml.control(trace = 0))
plotBS(midpoint(fitJC$tre), bs, p = 50, type="p")

