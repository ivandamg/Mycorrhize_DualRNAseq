#################################################################################
#setwd("~/PhD/GWAS-AMF/Global_cassava_response_rirregularis_M4/data_vf")
setwd("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/Phenodata_vf")
#################################################################################
#libraries
library('nlme')
#options("scipen"=100,digits=3) 
library('lme4')
library('ggplot2')
library(reshape)
library("multcomp")
par(las=2)
#################################################################################
####### 0. DIFFERENTIATION CULTURES CAN B1


RI30<-read.table("~/PhD/GWAS-AMF/Correlation_pheno_allelefq/Phylo_tree_ri30/DIAG_dist_20150420_BTSP_25_TableFormat.txt",h=T)
mat_ID<-t(RI30)
colnames(mat_ID)<-gsub("-(\\w+)$", "", colnames(mat_ID), perl=TRUE) 
rownames(mat_ID)<-gsub(".(\\w+)$", "", rownames(mat_ID), perl=TRUE) 
to_hclust_ID<-dist(mat_ID,method = "maximum")
tree_ID <- hclust(to_hclust_ID)
plot(tree_ID)

r_irr<-read.table("pheno_CAN_B1_R_irregularis_2013.txt",h=T)

boxplot(r_irr$count~r_irr$treat,ylab="spore production",cex.lab=1,cex.names=1,cex.axis=1)
wilcox.test(r_irr$count~r_irr$treat)


pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/coinoculation_manuscript/V5_coinoc/Supplementary_Figure4_V5.pdf", width=10, height=8, useDingbats = F)
par(mar = c(10, 11.5, 1.5, 0.5), mgp = c(8.5, 2.5, 0),las=1)
par(mfrow=c(1,2))
mean_sp<-tapply(r_irr$count,list(r_irr$dish,r_irr$treat),mean)
mean_sp2<-melt(mean_sp)
mean_sp2$X2<-relevel(mean_sp2$X2, ref = "CAN")
wilcox.test(mean_sp2$value~mean_sp2$X2)
boxplot(mean_sp2$value~mean_sp2$X2,cex.axis=3,cex.names=3,cex.lab=3,ylim=c(0,3000),ylab= "Spores per cm^2 within plate (n=5)")

#Analysis take into account random effecnt among mesures within replicate

pdf( '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/supplementary_material/1.Pheno_diff.pdf', width=10, height=10)

glm.1<-glmer(count~treat + (1 | dish), family=poisson,data=r_irr)
glm.1<-lmer(count~treat + (1 | dish),data=r_irr)
glm.2<-lmer(count~1+ (1 | dish),data=r_irr)

glm.2<-glmer(count~1 + (1 | dish), family=poisson,data=r_irr)
shapiro.test(residuals(glm.1))
anova(glm.2,glm.1,test="F")
anova(glm.1)
stDevs <-tapply(r_irr$count,list(r_irr$dish,r_irr$treat),sd)
means<-tapply(r_irr$count,list(r_irr$dish,r_irr$treat),mean)
mp<-barplot(tapply(r_irr$count,list(r_irr$dish,r_irr$treat),mean),beside=T,col = c(rep("burlywood4",5),rep("cadetblue4",5)),
            ylim=c(0,3500))
title(ylab= "Nb. spores in cm^2", xlab="Isolate")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

dev.off()

#par(mar = c(10, 11.5, 1.5, 0.5), mgp = c(8.5, 2.5, 0),las=1)
mean_sp<-tapply(r_irr$X100.,list(r_irr$dish,r_irr$treat),mean)
mean_sp2<-melt(mean_sp)
mean_sp2$X2<-relevel(mean_sp2$X2, ref = "CAN")
wilcox.test(mean_sp2$value~mean_sp2$X2)
boxplot(mean_sp2$value~mean_sp2$X2,ylab= "Spore clustering within a cm^2 (n=5)",cex.axis=3,cex.names=3,cex.lab=3,)
#dev.off()



################################
####### 1. COLONIZATION MEASUREMENTS
col<-read.table("colonization_vf.txt",h=T) #tablero
col<-col[col$treat!="CAN+B1",]

#take out non used data in RNASeq
#col<-col[col$treat!="A",]
col<-droplevels(col)

summary(col)
head(col)
col$treat<-droplevels(col$treat)
glm_col<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var*col$treat,family=binomial)
glm_col2<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var+col$treat,family=binomial)
glm_col3<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var,family=binomial)

anova(glm_col,glm_col2,test='Chisq')
anova(glm_col2,glm_col3,test='Chisq')
anova(glm_col)
summary(glm_col)

pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/manuscript_V1/PHENO1_Colonization_ALL.pdf")
stDevs <-tapply(col$perc,list(col$treat,col$var),sd)
means<-tapply(col$perc,list(col$treat,col$var),mean)
mp<-barplot(tapply(col$perc,list(col$treat,col$var),mean),beside=T,ylim=c(0,1.1),
         ylab="% of colonization (n=10)",las=2)
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
dev.off()

######################################################################################################
######################################################################################################


######################################################################################################
######################################################################################################
######################################################################################################
########### FINAL COLONIZATION ANALYSIS
### ANALYSIS WITHOUT CONTROL
col<-read.table("colonization_vf.txt",h=T) #tablero

col<-col[col$treat!="CANB1"&col$treat!="A",]
col<-col[col$var!="V2_COL1468",];col<-col[col$var!="V7_VEN329A",];col<-col[col$var!="V9_NGA-16",];col<-col[col$var!="V3_TAI-8",]

summary(col)
col$treat<-droplevels(col$treat)
col$var<-droplevels(col$var)

glm_col<-glm(cbind(Positive,(Obs-Positive))~var*treat,data= col,family=binomial)
glm_col2<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var+col$treat,family=binomial)
glm_col3<-glm(cbind(col$Positive,(col$Obs-col$Positive))~col$var,family=binomial)
anova(glm_col,glm_col2,test='Chisq')

anova(glm_col,test="Chisq")
summary(glm_col)


# multiple comparisons
glm_col <- glm(cbind(Positive,(Obs-Positive)) ~ var * treat, data= col,family=binomial)
library(lsmeans)
lsmeans(glm_col, pairwise ~ var |treat )



tmp <- expand.grid(treat = unique(col$treat), var = unique(col$var))
X <- model.matrix(~ var * treat, data = tmp)
summary(glht(glm_col, linfct = X))

predict(glm_col, newdata = tmp)


## differences between CAN B1
Tukey <- contrMat(table(col$treat), "Tukey")
K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*8))
rownames(K1) <- paste(levels(col$var)[1], rownames(K1), sep = ":")
K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey,matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*7))
rownames(K2) <- paste(levels(col$var)[2], rownames(K2), sep = ":")
K3 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*2), Tukey,matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*6))
rownames(K3) <- paste(levels(col$var)[3], rownames(K3), sep = ":")
K4 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*3), Tukey,matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*5))
rownames(K4) <- paste(levels(col$var)[4], rownames(K4), sep = ":")
K5 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*4), Tukey,matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*4))
rownames(K5) <- paste(levels(col$var)[5], rownames(K5), sep = ":")
K6 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*5), Tukey,matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*3))
rownames(K6) <- paste(levels(col$var)[6], rownames(K6), sep = ":")
K7 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*6), Tukey,matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*2))
rownames(K7) <- paste(levels(col$var)[7], rownames(K7), sep = ":")
K8 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*7), Tukey,matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
rownames(K8) <- paste(levels(col$var)[8], rownames(K8), sep = ":")
K9 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)*8), Tukey)
rownames(K9) <- paste(levels(col$var)[9], rownames(K9), sep = ":")
K <- rbind(K1, K2, K3, K4,K5,K6,K7,K8,K9)
colnames(K) <- rep(colnames(Tukey), 9)
#K <- rbind(K1, K2, K3, K4,K5)
#K<-K[,1:10]
#colnames(K) <- rep(colnames(Tukey), 5)

anova(glm_col)

summary(glht(glm_col, linfct = K %*% X,alternative="two.sided"),test=adjusted("Shaffer"))
K %*% predict(glm_col, newdata = tmp)
 
### conclusions: 1.variety effect by colonization, varieties were colonized differently by AMF.
#2. interaction var:treat, v7 different colonization of B1 and CAN

# plot
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V4/2.PHENO2_Colonization.pdf",width=12, height=10)
par(mfrow=c(1,1),mar=c(20,8,6,2),mgp=c(6, 2, 0))

boxplot(col$perc~interaction(col$treat,col$var),las=2,ylim=c(0,1),ylab="% of colonization (n=10)",
        col=c("azure3","dodgerblue3","lightsalmon3"),cex.lab=2,cex.axis=2)
dev.off()



######################################################################################################
######################################################################################################
######################################################################################################




################################
####### 3. FINAL WEIGHT and HEIGHT MEASUREMENTs
### TO MANUSCRIPT
################################

dry<-read.table("Harvest_VF.txt",h=T)
summary(dry)
dry$BLOCK<-as.factor(dry$BLOCK)
dry$VARIETY<-as.factor(dry$VARIETY)
dry$aerial<-as.numeric(dry$aerial)
dry$roots<-as.numeric(dry$roots)
dry$tuber<-as.numeric(dry$tuber)
dry$total<-as.numeric(dry$total)
dry$ratio<-as.numeric(dry$ratio)
dry$X125dai<-as.numeric(dry$X125dai)

dry$TREATMENT<-relevel(dry$TREATMENT, ref = "C") #order factors to put control as reference
dry$AMF<-relevel(dry$AMF, ref = "CTRL") #order factors to put control as reference

all_root<-dry$roots+dry$tuber
dry<-cbind(dry,all_root)
dry<-dry[dry$TREATMENT!="A+B",]

#dry<-dry[dry$TREATMENT!="C",]
dry<-dry[dry$VARIETY!=2,];dry<-dry[dry$VARIETY!=7,];dry<-dry[dry$VARIETY!=9,];dry<-dry[dry$VARIETY!=3,]
dry<-droplevels(dry)
##checking outlayers to filter 

barplot(dry

#filter block no.11 and no 14. lowest variances
#var(dry$total[dry$block==14])
#var(dry$total[dry$block==11])
#dry<-dry[dry$block!=14&dry$block!=11,]
pdf("~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/2.PHENO_harvest_v4.pdf",width=18, height=16)
par(mfrow=c(2,3),mar=c(20,8,6,2),mgp=c(6, 2, 0))

boxplot(dry$total~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Total dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3"))
boxplot(dry$aerial~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Aboveground dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3"))
boxplot(dry$roots~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Fine roots dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3"))
boxplot(dry$tuber~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Tuberized roots dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3"))
boxplot(dry$all_root~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Underground dry weight (n=14)",col=c("azure3","dodgerblue3","lightsalmon3"))
boxplot(dry$X125dai~interaction(dry$AMF,dry$CULTIVAR),las=2,cex.axis=2,cex.lab=2,
        ylab="Height (n=14)",col=c("azure3","dodgerblue3","lightsalmon3"))
dev.off()

par(mfrow=c(2,3))
### ALL VARIETYieties

lsmeans(ALLV_lme, pairwise ~ AMF | CULTIVAR)

######### total dry weight
ALLV_lme<-lme(total~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)
########################
ALLV_lme<-lme(aerial~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)

ALLV_lme<-lme(roots~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)


ALLV_lme<-lme(all_root~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)

ALLV_lme<-lme(tuber~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)


ALLV_lme<-lme(X125dai~AMF*CULTIVAR,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
summary(ALLV_lme)
#summary(glht(ALLV_lme, linfct=mcp(CULTIVAR="Tukey")) )
#summary(glht(ALLV_lme, linfct=mcp(AMF="Tukey")) )


### ALL VARIETYieties
# Total
ALLV_lme<-lme(total~VARIETY*TREATMENT,random=~1|BLOCK,data=dry)
anova(ALLV_lme)

stDevs <-tapply(dry$total,list(dry$TREATMENT,dry$VARIETY),sd)
means<-tapply(dry$total,list(dry$TREATMENT,dry$VARIETY),mean)
mp<-barplot(tapply(dry$total,list(dry$TREATMENT,dry$VARIETY),mean),beside=T,ylim=c(0,21),
            ylab="Total dry weight of plant (n=14)")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

#paste(as.character(anova(ALLV_lme)[2,]),collapse="/ ")


table(list(dry$TREATMENT,dry$VARIETY))
#legend(1,18,c("Control","CAN","B1"),pch=16)

# aerial

ALLV_lme<-lme(aerial~VARIETY*TREATMENT,random=~1|BLOCK,data=dry)
anova(ALLV_lme)

stDevs <-tapply(dry$aerial,list(dry$TREATMENT,dry$VARIETY),sd)
means<-tapply(dry$aerial,list(dry$TREATMENT,dry$VARIETY),mean)
mp<-barplot(tapply(dry$aerial,list(dry$TREATMENT,dry$VARIETY),mean),beside=T,ylim=c(0,5),
            ylab="aerial dry weight of plant (n=14)")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

 #legend(1,18,c("Control","CAN","B1"),pch=16)
# roots
ALLV_lme<-lme(roots~VARIETY*TREATMENT,random=~1|BLOCK,data=dry)
anova(ALLV_lme)

stDevs <-tapply(dry$roots,list(dry$TREATMENT,dry$VARIETY),sd)
means<-tapply(dry$roots,list(dry$TREATMENT,dry$VARIETY),mean)
mp<-barplot(tapply(dry$roots,list(dry$TREATMENT,dry$VARIETY),mean),beside=T,ylim=c(0,2),
            ylab="roots dry weight of plant (n=14)")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

#legend(1,4,c("Control","CAN","B1"),pch=16)
# tuber
ALLV_lme<-lme(tuber~VARIETY*TREATMENT,random=~1|BLOCK,data=dry)
anova(ALLV_lme)

stDevs <-tapply(dry$tuber,list(dry$TREATMENT,dry$VARIETY),sd)
means<-tapply(dry$tuber,list(dry$TREATMENT,dry$VARIETY),mean)
mp<-barplot(tapply(dry$tuber,list(dry$TREATMENT,dry$VARIETY),mean),beside=T,ylim=c(0,15),
            ylab="tuber dry weight of plant (n=14)")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
#legend(1,13,c("Control","CAN","B1"),pch=16)
# all_root
ALLV_lme<-lme(all_root~VARIETY*TREATMENT,random=~1|BLOCK,data=dry)
anova(ALLV_lme)

stDevs <-tapply(dry$all_root,list(dry$TREATMENT,dry$VARIETY),sd)
means<-tapply(dry$all_root,list(dry$TREATMENT,dry$VARIETY),mean)
mp<-barplot(tapply(dry$all_root,list(dry$TREATMENT,dry$VARIETY),mean),beside=T,ylim=c(0,15),
            ylab="all_root dry weight of plant (n=14)")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)
#legend(1,13,c("Control","CAN","B1"),pch=16)
# ratio
ALLV_lme<-lme(ratio~VARIETY*TREATMENT,random=~1|BLOCK,data=dry)
anova(ALLV_lme)

stDevs <-tapply(dry$ratio,list(dry$TREATMENT,dry$VARIETY),sd)
means<-tapply(dry$ratio,list(dry$TREATMENT,dry$VARIETY),mean)
mp<-barplot(tapply(dry$ratio,list(dry$TREATMENT,dry$VARIETY),mean),beside=T,ylim=c(0,15),
            ylab="ratio dry weight of plant (n=14)")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

# height
ALLV_lme<-lme(X125dai~VARIETY*TREATMENT,random=~1|BLOCK,data=dry)
anova(ALLV_lme)
boxplot(X125dai~interaction(VARIETY,TREATMENT),data=dry, col=c(rep(1,6),rep(2,6),rep(3,6)))
stDevs <-tapply(dry$X125dai,list(dry$TREATMENT,dry$VARIETY),sd)
means<-tapply(dry$X125dai,list(dry$TREATMENT,dry$VARIETY),mean)
mp<-barplot(tapply(dry$X125dai,list(dry$TREATMENT,dry$VARIETY),mean),beside=T,
            ylab="Height of plant (n=14)")
segments(mp, means - stDevs, mp, means + stDevs, lwd=2)
segments(mp - 0.1, means - stDevs, mp + 0.1, means - stDevs, lwd=2)
segments(mp - 0.1, means + stDevs, mp + 0.1, means + stDevs, lwd=2)

