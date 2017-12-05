# gene-gene interaction script between functional groups.
# 27.10.2016
# after script RNA_seq_gene..._v11.R

##############################################################################
# Visualization of functional genes correlation against AMF and PLANT
library(pheatmap)
# all top genes
CasR<-unique(c("Manes.02G039800","Manes.02G059900","Manes.03G068200","Manes.03G068200","Manes.09G146200","Manes.09G146200","Manes.09G156800","Manes.09G156800","Manes.10G135100","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G126200","Manes.14G073700","Manes.17G018800","Manes.17G018800","Manes.02G039800","Manes.03G068200","Manes.03G068200","Manes.08G141200","Manes.09G146200","Manes.09G146200","Manes.09G156800","Manes.09G156800","Manes.13G126200","Manes.16G010100","Manes.17G018800","Manes.17G018800","Manes.02G039800","Manes.03G068200","Manes.03G068200","Manes.06G148500","Manes.06G148600","Manes.07G094600","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G060300","Manes.13G126200","Manes.16G052800","Manes.17G019200","Manes.01G023200","Manes.01G111700","Manes.01G111700","Manes.01G123000","Manes.01G139400","Manes.02G069900","Manes.02G069900","Manes.02G069900","Manes.02G069900","Manes.03G127800","Manes.03G171200","Manes.03G202500","Manes.03G202500","Manes.05G047000","Manes.05G047000","Manes.06G063500","Manes.06G117800","Manes.06G155200","Manes.06G155200","Manes.06G164600","Manes.07G056100","Manes.09G059900","Manes.09G084500","Manes.09G087100","Manes.10G007300","Manes.11G004200","Manes.12G063500","Manes.12G063500","Manes.12G063500","Manes.14G007500","Manes.14G100900","Manes.14G108900","Manes.15G034400","Manes.15G094700","Manes.15G151100","Manes.15G151100","Manes.17G013000","Manes.01G023200"))
CasR2<- unique(c("Manes.01G111700","Manes.01G111700","Manes.01G139400","Manes.03G127800","Manes.03G171200","Manes.04G075000","Manes.04G075000","Manes.05G077100","Manes.05G077100","Manes.06G063500","Manes.06G117800","Manes.06G155200","Manes.06G155200","Manes.06G164600","Manes.07G056100","Manes.07G073600","Manes.09G059900","Manes.09G087100","Manes.10G007300","Manes.11G004200","Manes.12G035300","Manes.12G035300","Manes.12G063500","Manes.12G063500","Manes.12G063500","Manes.14G100900","Manes.14G108900","Manes.14G139600","Manes.15G094700","Manes.15G151100","Manes.15G151100","Manes.17G110900","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G036800","Manes.01G115400","Manes.01G265100","Manes.03G206300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G141400","Manes.04G146400","Manes.05G047400","Manes.06G038200","Manes.06G122900","Manes.06G122900","Manes.08G049200","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G129100","Manes.09G129100","Manes.10G034400","Manes.10G034400","Manes.10G059900","Manes.11G114800","Manes.11G114800","Manes.12G013300","Manes.13G045600","Manes.14G015500","Manes.14G015500","Manes.15G165000","Manes.16G110700","Manes.17G059000","Manes.18G073000","Manes.01G115400","Manes.01G265100","Manes.02G119500","Manes.02G219500","Manes.03G206300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G146400","Manes.05G047400","Manes.06G038200","Manes.06G122900","Manes.06G122900","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.08G049200","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.10G121900","Manes.10G121900","Manes.14G015500","Manes.14G015500","Manes.15G165000","Manes.17G053000","Manes.17G059000","Manes.01G205200","Manes.01G205200","Manes.01G205200","Manes.04G146400","Manes.05G014200","Manes.05G136900","Manes.06G122900","Manes.06G122900","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.06G126700","Manes.08G010900","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G021000","Manes.09G129100","Manes.09G129100","Manes.10G059900","Manes.10G121900","Manes.10G121900","Manes.12G013300","Manes.15G165000","Manes.17G059000","Manes.18G031200","Manes.18G073000","Manes.18G094200","Manes.02G119500","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G064300","Manes.04G090100","Manes.04G146400","Manes.04G161100","Manes.04G161100","Manes.04G161100","Manes.04G161100","Manes.10G116300","Manes.13G070400","Manes.14G015500","Manes.14G015500","Manes.17G053000","Manes.17G059000","Manes.01G208500","Manes.02G135300","Manes.06G041300","Manes.07G034500","Manes.14G122800","Manes.15G022900","Manes.01G201800","Manes.03G044400","Manes.07G087500","Manes.11G078100","Manes.13G008000","Manes.16G090600","Manes.16G090600","Manes.16G090600","Manes.03G205000","Manes.03G205000","Manes.04G095800","Manes.11G078100","Manes.14G061800","Manes.15G087400","Manes.16G090600","Manes.16G090600","Manes.16G090600","Manes.03G205000","Manes.03G205000","Manes.09G120000","Manes.09G120000","Manes.09G182400","Manes.14G141400","Manes.18G060600","Manes.18G060600","Manes.04G095800","Manes.09G182400","Manes.07G005600","Manes.07G092000","Manes.16G024800","Manes.01G003900","Manes.11G154700","Manes.16G020200","Manes.16G133800","Manes.16G133800","Manes.03G194400","Manes.11G154700","Manes.13G084400","Manes.16G133800","Manes.16G133800","Manes.01G003900","Manes.09G041300","Manes.11G124400","Manes.11G154700","Manes.13G071100"))
CasR3<-unique(c(CasR,CasR2))
cas<-CasR3      
amf<-unique(c("g1008.t1","g11187.t1","g11634.t1","g2998.t1","g3197.t1","g339.t1","g3641.t1","g4193.t1","g4624.t1","g5528.t1","g7192.t1","g7415.t1","g8161.t1","g8671.t1","g8671.t1","g8680.t1","g10995.t1","g1988.t1","g2737.t1","g3047.t1","g568.t1","g8590.t1","g105.t1","g10703.t1","g1207.t1","g1757.t1","g21.t1","g2329.t1","g3510.t1","g3844.t1","g5111.t1","g5456.t1","g5866.t1","g652.t1","g8553.t1","g8944.t1","g929.t1","g9749.t1","g10683.t1","g10723.t1","g10723.t1","g10723.t1","g11991.t1","g12075.t1","g12890.t1","g1385.t1","g1651.t1","g2394.t1","g2427.t1","g2715.t1","g2879.t1","g3065.t1","g3508.t1","g4871.t1","g4952.t1","g5145.t1","g5801.t1","g5821.t1","g5979.t1","g7641.t1","g8429.t1","g9899.t1","g14952.t1","g14952.t1","g4959.t1","g4978.t1","g5557.t1","g7735.t1","g9089.t1","g1008.t1","g11187.t1","g1565.t1","g1659.t1","g339.t1","g4193.t1","g7192.t1","g8161.t1","g8679.t1","g8727.t1","g8738.t1","g8981.t1","g10995.t1","g1300.t1","g1988.t1","g6205.t1","g8590.t1","g8780.t1","g9390.t1","g105.t1","g10703.t1","g1072.t1","g1486.t1","g1486.t1","g1486.t1","g177.t1","g2594.t1","g4372.t1","g5160.t1","g5456.t1","g5866.t1","g7772.t1","g8553.t1","g8944.t1","g9013.t1","g929.t1","g1012.t1","g10122.t1","g1022.t1","g1045.t1","g11736.t1","g12309.t1","g12727.t1","g12728.t1","g1339.t1","g1488.t1","g1507.t1","g1999.t1","g205.t1","g208.t1","g2331.t1","g2417.t1","g2421.t1","g2610.t1","g2610.t1","g2657.t1","g2662.t1","g2808.t1","g2910.t1","g3270.t1","g3327.t1","g3371.t1","g3377.t1","g3830.t1","g3904.t1","g4338.t1","g4338.t1","g4338.t1","g4376.t1","g441.t1","g4851.t1","g5093.t1","g526.t1","g5288.t1","g5311.t1","g548.t1","g5691.t1","g5697.t1","g614.t1","g623.t1","g629.t1","g638.t1","g6404.t1","g6408.t1","g6408.t1","g642.t1","g65.t1","g6748.t1","g6910.t1","g7241.t1","g7356.t1","g7358.t1","g7358.t1","g7654.t1","g8517.t1","g8613.t1","g8748.t1","g881.t1","g9174.t1","g9445.t1","g973.t1","g9739.t1","g9840.t1","g992.t1","g1008.t1","g11187.t1","g1565.t1","g2022.t1","g8679.t1","g8727.t1","g1004.t1","g1004.t1","g11991.t1","g12075.t1","g1651.t1","g2394.t1","g2879.t1","g3431.t1","g3431.t1","g3508.t1","g3613.t1","g4807.t1","g4871.t1","g4952.t1","g4976.t1","g5145.t1","g5801.t1","g5979.t1","g7641.t1","g9899.t1","g1236.t1","g1236.t1","g3417.t1","g407.t1","g4188.t1","g5983.t1","g8841.t1","g11187.t1","g11634.t1","g1565.t1","g2998.t1","g339.t1","g4193.t1","g7192.t1","g7415.t1","g8738.t1","g8952.t1","g8981.t1","g10703.t1","g1072.t1","g13844.t1","g177.t1","g2329.t1","g2340.t1","g2594.t1","g5160.t1","g5395.t1","g5456.t1","g7698.t1","g8944.t1","g9013.t1","g9388.t1","g1108.t1","g11163.t1","g11668.t1","g2567.t1","g3041.t1","g3094.t1","g3228.t1","g4191.t1","g4191.t1","g4191.t1","g5338.t1","g5338.t1","g6538.t1","g6787.t1","g6816.t1","g7673.t1","g8455.t1","g928.t1","g981.t1","g9834.t1","g10019.t1","g10098.t1","g10803.t1","g1165.t1","g11806.t1","g12383.t1","g1242.t1","g13288.t1","g1570.t1","g2056.t1","g2303.t1","g2303.t1","g2305.t1","g2345.t1","g2349.t1","g2350.t1","g2968.t1","g3399.t1","g372.t1","g3836.t1","g3987.t1","g4051.t1","g4435.t1","g4435.t1","g4998.t1","g5163.t1","g5253.t1","g5665.t1","g5802.t1","g6035.t1","g6165.t1","g6165.t1","g646.t1","g6853.t1","g7670.t1","g7867.t1","g8801.t1","g10271.t1","g11031.t1","g382.t1","g5472.t1","g6255.t1","g6598.t1","g687.t1","g7716.t1","g7775.t1","g8837.t1","g10076.t1","g10292.t1","g11821.t1","g12154.t1","g12914.t1","g13967.t1","g1518.t1","g2043.t1","g2337.t1","g2510.t1","g2748.t1","g284.t1","g296.t1","g30.t1","g3468.t1","g410.t1","g4281.t1","g4767.t1","g4790.t1","g5164.t1","g5223.t1","g5465.t1","g5995.t1","g6405.t1","g6750.t1","g8058.t1","g8061.t1","g8188.t1","g8763.t1","g9477.t1","g9871.t1","g9894.t1"))
#########

###########################################################
################## AMF FUNCTIONAL GROUPS ################## 
###########################################################
###### Secretory pathway AMF

target_genes<- c("g2715.t1","g3508.t1","g382.t1","g5528.t1","g568.t1","g7358.t1","g3371.t1","g8780.t1","g10076.t1","g6408.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_secretory.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

###### Signalling AMF
target_genes<-c("g6404.t1","g4376.t1","g6910.t1","g4338.t1","g3431.t1","g1236.t1","g8671.t1","g4191.t1","g10995.t1","g8680.t1","g9390.t1","g2340.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_signalling.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### ABC transporter AMF
target_genes<-c("g1045.t1","g5801.t1","g9899.t1","g1486.t1","g646.t1","g1486.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_ABCtransport.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

###### Protein synthesis AMF
target_genes<-c("g1207.t1","g2594.t1","g284.t1","g1339.t1","g10122.t1","g1488.t1","g623.t1","g12309.t1","g9840.t1","g638.t1","g5691.t1","g2808.t1","g5866.t1","g11736.t1","g8590.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_Protsynthesis.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

###### Protein post translational modification AMF
target_genes<-c("g1236.t1","g5995.t1","g3431.t1","g11806.t1","g7867.t1","g2998.t1","g8671.t1","g4191.t1","g4338.t1","g6408.t1","g7358.t1","g4191.t1","g5338.t1","g4338.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_Protposttransmod.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Protein degradation AMF
target_genes<-c("g7356.t1","g4998.t1","g3641.t1","g8161.t1","g12728.t1","g2417.t1","g2662.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_Protdegradation.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Lipid metabo AMF
target_genes<-c("g8763.t1","g8429.t1","g12727.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_Lipidmetabo.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Hormone metabo AMF
target_genes<-c("g14952.t1","g2610.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_Hormonemetabo.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Cell organization AMF
target_genes<-c("g2305.t1","g65.t1","g1004.t1","g10723.t1","g4807.t1","g12383.t1","g5253.t1","g7670.t1","g2329.t1","g3510.t1","g5395.t1","g11821.t1","g5311.t1","g9174.t1","g10723.t1","g10703.t1","g10683.t1","g7641.t1","g6853.t1","g4193.t1","g1004.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_CELLorg.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Cell vesicle transport AMF
target_genes<-c("g7192.t1","g1300.t1")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_Y<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_Y<-cor2plot_Y[grep("Manes",colnames(cor2plot_Y)),]
cor2plot_Y<-cor2plot_Y[,grep(".t1",colnames(cor2plot_Y))]
colnames(cor2plot_Y)
# amf genes
NORM_EXP_TOPAMF_ALL<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf,])
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% target_genes ,])
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPAMF_ALL))
cor2plot_A<-cor2plot_A[colnames(cor2plot_A) %in% target_genes,]
cor2plot_A<-cor2plot_A[,!colnames(cor2plot_A) %in% target_genes]
cor2plot_A<-cor2plot_A[,grep("\\.1",colnames(cor2plot_A),invert=T)]
colnames(t(cor2plot_A))
write.table( rbind.data.frame(cor2plot_Y,t(cor2plot_A)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/AMF_Cellvesicletransport.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)




###########################################################
################## CASSAVA FUNCTIONAL GROUPS ############## 
###########################################################

###### Cell organisation cassava
target_genes<-c("Manes.09G021000","Manes.05G047400","Manes.04G064300")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Cellorganisation.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Glycolysis cassava
target_genes<-c("Manes.01G023200","Manes.09G120000","Manes.01G123000","Manes.18G060600","Manes.14G007500","Manes.06G155200","Manes.15G034400","Manes.03G171200","Manes.06G155200","Manes.07G056100")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Glycolysis.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Lipid metabolism cassava
target_genes<-c("Manes.10G007300","Manes.04G095800","Manes.01G139400","Manes.03G202500","Manes.14G100900","Manes.06G155200")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Lipidmetabo.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

###### MAjorCHO metabolism cassava
target_genes<-c("Manes.02G069900","Manes.01G111700","Manes.02G069900","Manes.03G044400","Manes.16G090600")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_MAjorCHOmetabo.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)



###### Protein degradation cassava
target_genes<-c("Manes.17G013000","Manes.01G205200","Manes.08G141200","Manes.06G148500","Manes.06G148600","Manes.16G133800","Manes.14G073700","Manes.16G133800","Manes.02G219500","Manes.06G041300","Manes.09G041300","Manes.13G084400","Manes.06G122900")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Protdegradation.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

###### Protein fold_synthesi_modification cassava
target_genes<-c("Manes.03G206300","Manes.10G121900","Manes.16G010100","Manes.15G151100","Manes.09G146200")

NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Protfoldsynthmod.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

###### secretory pathway cassava
target_genes<-c("Manes.09G156800","Manes.17G019200")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Secretory.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Calvin cycle cassava
target_genes<-c("Manes.06G063500","Manes.14G108900","Manes.05G047000")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Calvin.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

###### Signalling cassava
target_genes<-c("Manes.01G265100","Manes.17G059000","Manes.02G135300","Manes.15G165000")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Signalling.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### Transport cassava
target_genes<-c("Manes.07G034500","Manes.03G127800","Manes.06G164600","Manes.11G078100","Manes.14G141400","Manes.01G111700","Manes.02G069900","Manes.17G018800","Manes.01G003900","Manes.01G208500")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_transport.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)


###### RNAregulation cassava
target_genes<-c("Manes.14G061800","Manes.01G115400","Manes.07G094600","Manes.11G124400","Manes.10G116300")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_RNAregulation.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

###### Secondary metabolism cassava
target_genes<-c("Manes.12G063500","Manes.09G084500","Manes.07G087500","Manes.13G071100")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot_A<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot_A<-cor2plot_A[grep("Manes",colnames(cor2plot_A)),]
cor2plot_A<-cor2plot_A[,grep(".t1",colnames(cor2plot_A))]

t(cor2plot_A)
# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot_Y<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot_Y<-cor2plot_Y[colnames(cor2plot_Y) %in% target_genes,]
cor2plot_Y<-cor2plot_Y[,grep("\\.1",colnames(cor2plot_Y),invert=T)]
t(cor2plot_Y)
write.table( rbind.data.frame(t(cor2plot_A),t(cor2plot_Y)),
             '~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_Secondarymetabo.txt',
             quote=FALSE,sep='\t',col.names = T, row.names = T)

######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################


###########################################################
################## MOdule module correlation ############## 
###########################################################


#####################################################################

### TO IMPROVE BECAUSE BAD LABELING 

#nb of positive negative correlations between modules
names(top_gene_sign_intramod_conec_A)
names(top_gene_sign_intramod_conec_Y)
unlist(top_gene_sign_intramod_conec_A[[1]])
unlist(top_gene_sign_intramod_conec_Y[[18]])
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
  
  cor2plot<-cor(NORM_EXP_TOPYUCA,NORM_EXP_TOPAMF)
  #pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_2ndmetabolism_amf.pdf',width=32,height=15,useDingbats = FALSE)
  target_function<-Mercator_Rirregularis[Mercator_Rirregularis$gene %in% unlist(top_gene_sign_intramod_conec_A[[i]]) ,]
  target_function<-target_function[!duplicated(target_function$gene),]
  target_function2<-Mercator_Mesculenta[Mercator_Mesculenta$gene %in% unlist(top_gene_sign_intramod_conec_Y[[i]]) ,]
  target_function2<-target_function2[!duplicated(target_function2$gene),]
  cor2mod<-pheatmap(cor2plot)
  rownames(cor2plot)[cor2mod$tree_row[[3]]]
  pheatmap(cor2plot ,cellwidth = 8,cellheight = 8,labels_col = paste (colnames(cor2plot),target_function$Function[order(match(target_function$gene,colnames(cor2plot)))],sep="___"),
           labels_row = paste (rownames(cor2plot),target_function2$Function[order(match(target_function2$gene,rownames(cor2plot)))],sep="___"))
  #dev.off()

}
nb_cor
unlist(nb_genes_Ao)
(length(cor(NORM_EXP_TOPYUCA))-dim(cor(NORM_EXP_TOPYUCA))) /2




#############

names(top_gene_sign_intramod_conec_A)

###### Secondary metabolism cassava
target_genes<-c("Manes.12G063500","Manes.09G084500","Manes.07G087500","Manes.13G071100")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# amf genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes  ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot<-cor2plot[grep("Manes",colnames(cor2plot)),]
cor2plot<-cor2plot[,grep(".t1",colnames(cor2plot))]


pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_2ndmetabolism_amf.pdf',width=32,height=15,useDingbats = FALSE)
target_function<-Mercator_Rirregularis[Mercator_Rirregularis$gene %in% amf ,]
target_function<-target_function[!duplicated(target_function$gene),]
target_function2<-Mercator_Mesculenta[Mercator_Mesculenta$gene %in% cas ,]
target_function2<-target_function2[!duplicated(target_function2$gene),]

pheatmap(cor2plot ,cellwidth = 8,cellheight = 8,labels_col = paste (colnames(cor2plot),target_function$Function[order(match(target_function$gene,colnames(cor2plot)))],sep="___"),
         labels_row = paste (colnames(t(cor2plot)),target_function2$Function[order(match(target_function2$gene,colnames(t(cor2plot))))],sep="___"))
dev.off()

# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  target_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot<-cor2plot[colnames(cor2plot) %in% target_genes,]
cor2plot<-cor2plot[,grep("\\.1",colnames(cor2plot),invert=T)]
pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_2ndmetabolism_cassava.pdf',width=10,height=10,useDingbats = FALSE)
target_function<-Mercator_Mesculenta[Mercator_Mesculenta$gene %in% cas ,]
target_function<-target_function[!duplicated(target_function$gene),]
target_function2<-Mercator_Mesculenta[Mercator_Mesculenta$gene %in% cas ,]
target_function2<-target_function2[!duplicated(target_function2$gene),]
pheatmap(cor2plot ,cellwidth = 8,cellheight = 8,labels_col = paste (colnames(cor2plot),target_function$Function[order(match(target_function$gene,colnames(cor2plot)))],sep="___"),
         labels_row = paste (colnames(t(cor2plot)),target_function2$Function[order(match(target_function2$gene,colnames(t(cor2plot))))],sep="___"))
dev.off()




##################
unlist(top_gene_sign_intramod_conec_Y[[8]])


# plot gene transcription of top genes of each module
# cassava
for (i in 1:length(top_gene_sign_intramod_conec_Y)) {
blue_Y_V1<-sign_YUCA_V1[rownames(sign_YUCA_V1) %in% unlist(top_gene_sign_intramod_conec_Y[[i]]),]
blue_Y_V4<-sign_YUCA_V4[rownames(sign_YUCA_V4) %in% unlist(top_gene_sign_intramod_conec_Y[[i]]),]
blue_Y_V5<-sign_YUCA_V5[rownames(sign_YUCA_V5) %in% unlist(top_gene_sign_intramod_conec_Y[[i]]),]
blue_Y_V6<-sign_YUCA_V6[rownames(sign_YUCA_V6) %in% unlist(top_gene_sign_intramod_conec_Y[[i]]),]
blue_Y_V8<-sign_YUCA_V8[rownames(sign_YUCA_V8) %in% unlist(top_gene_sign_intramod_conec_Y[[i]]),]

pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/topgenes_inmodules_expression/',
          sapply(strsplit(names(top_gene_sign_intramod_conec_Y[i]),"_to_"), "[[", 1),".pdf",sep=""),width=10,height=10,useDingbats = FALSE)
pheatmap(data.frame(blue_Y_V1[,1:3],blue_Y_V4[,1:3],blue_Y_V5[,1:3],
                    blue_Y_V6[,1:3],blue_Y_V8[,1:3]),
         main = sapply(strsplit(names(top_gene_sign_intramod_conec_Y[i]),"_to_"), "[[", 1) )
dev.off()
}

# AMF
for (i in 1:length(top_gene_sign_intramod_conec_A)) {
  blue_Y_V1<-sign_AMF_V1[rownames(sign_AMF_V1) %in% unlist(top_gene_sign_intramod_conec_A[[i]]),]
  blue_Y_V4<-sign_AMF_V4[rownames(sign_AMF_V4) %in% unlist(top_gene_sign_intramod_conec_A[[i]]),]
  blue_Y_V5<-sign_AMF_V5[rownames(sign_AMF_V5) %in% unlist(top_gene_sign_intramod_conec_A[[i]]),]
  blue_Y_V6<-sign_AMF_V6[rownames(sign_AMF_V6) %in% unlist(top_gene_sign_intramod_conec_A[[i]]),]
  blue_Y_V8<-sign_AMF_V8[rownames(sign_AMF_V8) %in% unlist(top_gene_sign_intramod_conec_A[[i]]),]
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/topgenes_inmodules_expression/',
            sapply(strsplit(names(top_gene_sign_intramod_conec_A[i]),"_to_"), "[[", 1),".pdf",sep=""),width=10,height=10,useDingbats = FALSE)
  a_mod<-data.frame(blue_Y_V1,blue_Y_V4,blue_Y_V5,blue_Y_V6,blue_Y_V8)[,c(1,7,13,19,25)]
  colnames(a_mod)<-c("CAN_B1_V1","CAN_B1_V4","CAN_B1_V5","CAN_B1_V6","CAN_B1_V8")
  pheatmap(a_mod,
           main = sapply(strsplit(names(top_gene_sign_intramod_conec_A[i]),"_to_"), "[[", 1) )
  dev.off()
}




########################
Secreted_prot_lin<-c("g15235.t1","g14285.t1","g6828.t1","g11848.t1","g4212.t1","g2806.t1","g15235.t1","g7946.t1","g12227.t1","g8989.t1","g10176.t1","g11432.t1","g8749.t1","g13285.t1","g9659.t1","g12948.t1","g11823.t1","g3750.t1","g1241.t1","g12185.t1","g5958.t1","g10809.t1","g15428.t1","g5943.t1","g14837.t1","g3797.t1","g12202.t1","g7237.t1","g9470.t1","g11593.t1","g14668.t1","g5587.t1","g5586.t1","g10957.t1","g7816.t1","g7817.t1","g13231.t1","g12114.t1","g13767.t1","g13250.t1","g14364.t1","g14466.t1","g4797.t1","g5990.t1","g3639.t1","g14647.t1","g5993.t1","g11554.t1","g4474.t1","g6328.t1","g6326.t1","g4129.t1","g7478.t1","g8040.t1","g6839.t1","g4253.t1","g1101.t1","g10533.t1","g10535.t1","g10536.t1","g12225.t1","g9301.t1","g3101.t1","g5274.t1","g6732.t1","g14036.t1","g11335.t1","g12490.t1","g6251.t1","g9130.t1","g10211.t1","g2231.t1","g11124.t1","g13656.t1","g13196.t1","g13196.t1","g10799.t1","g10799.t1","g10799.t1","g1802.t1","g1802.t1","g1802.t1","g12953.t1","g3134.t1","g7509.t1","g6154.t1","g6154.t1","g9609.t1","g3978.t1","g3978.t1","g10879.t1","g9868.t1","g6678.t1","g13952.t1","g11500.t1","g7109.t1","g6090.t1","g11414.t1","g1225.t1","g2085.t1","g5938.t1","g4403.t1","g8016.t1","g6501.t1","g2151.t1","g11979.t1","g11978.t1","g4922.t1","g13939.t1","g9018.t1","g77.t1","g6559.t1","g77.t1","g8009.t1","g8684.t1","g4746.t1","g11608.t1","g10209.t1","g272.t1","g353.t1","g11544.t1","g5583.t1","g3271.t1","g4837.t1","g11199.t1","g4837.t1","g9681.t1","g14888.t1","g6353.t1","g4696.t1","g4964.t1","g11386.t1","g5592.t1","g10859.t1","g11780.t1","g12321.t1","g13171.t1","g13041.t1","g9300.t1","g6365.t1","g8335.t1","g5910.t1","g3514.t1","g5893.t1","g12354.t1","g1016.t1","g1016.t1","g6151.t1","g13590.t1","g11327.t1","g13313.t1","g3619.t1","g1696.t1","g1696.t1","g8243.t1","g7469.t1","g8525.t1","g10567.t1","g10638.t1","g14082.t1","g7534.t1","g6805.t1","g6805.t1","g6805.t1","g8535.t1","g1471.t1","g1264.t1","g10651.t1","g10651.t1","g3293.t1","g859.t1","g1071.t1","g6954.t1","g6954.t1","g13912.t1","g2320.t1","g4330.t1","g1004.t1","g4704.t1","g4244.t1","g14198.t1","g15913.t1","g12132.t1","g15203.t1","g10227.t1","g10463.t1","g4460.t1","g1730.t1","g10363.t1","g10419.t1","g13396.t1","g6639.t1","g7260.t1","g7260.t1","g12427.t1","g4682.t1","g4682.t1","g10904.t1","g5964.t1","g918.t1","g4934.t1","g9873.t1","g6307.t1","g6308.t1","g2195.t1","g6204.t1","g3564.t1","g6225.t1","g11828.t1","g962.t1","g7582.t1","g5770.t1","g5771.t1","g2347.t1","g6060.t1","g9590.t1","g4197.t1","g677.t1","g8125.t1","g3208.t1","g4305.t1","g7521.t1","g1633.t1","g2745.t1","g1627.t1","g6662.t1","g193.t1","g4267.t1","g54.t1","g2953.t1","g1879.t1","g10945.t1","g1836.t1","g13690.t1","g11951.t1","g11411.t1","g1456.t1","g7857.t1","g13682.t1","g11948.t1","g10892.t1","g7810.t1","g7809.t1","g9646.t1","g5013.t1","g11476.t1","g8931.t1","g12764.t1","g12093.t1","g13574.t1","g11209.t1","g9414.t1","g14505.t1","g13788.t1","g5070.t1","g6511.t1","g3321.t1","g10163.t1","g4901.t1","g11498.t1","g11362.t1","g13093.t1","g5516.t1","g12407.t1","g8213.t1","g6827.t1","g6636.t1","g11157.t1","g11862.t1","g13208.t1","g5032.t1","g13237.t1","g5636.t1","g4890.t1","g8170.t1","g8573.t1","g431.t1","g7884.t1","g14311.t1","g285.t1","g917.t1","g924.t1","g974.t1","g12784.t1","g5208.t1","g7821.t1","g352.t1","g2113.t1","g9904.t1")

# N6
Secreted_prot_N6<-c("g53.t1","g940.t1","g1107.t1","g1800.t1","g2001.t1","g2342.t1","g2499.t1","g3246.t1","g3407.t1","g3897.t1","g4420.t1","g4483.t1","g5262.t1","g5274.t1","g5397.t1","g5495.t1","g5572.t1","g5591.t1","g5675.t1","g5718.t1","g5986.t1","g6092.t1","g6615.t1","g6632.t1","g6761.t1","g6976.t1","g7070.t1","g7210.t1","g7589.t1","g7936.t1","g7970.t1","g7999.t1","g8143.t1","g9852.t1","g10016.t1","g10032.t1","g10816.t1","g10921.t1","g11538.t1","g11588.t1","g12102.t1","g12555.t1","g12841.t1","g13513.t1","g13850.t1","g14020.t1","g14745.t1","g15061.t1","g15327.t1","g15521.t1","g15637.t1","g15940.t1")
Topgenes_AMF<-unique(as.vector(unlist(top_gene_sign_intramod_conec_A)))

Topgenes_AMF[Topgenes_AMF %in% Secreted_prot_lin]


# AMF
  blue_Y_V1<-sign_AMF_V1[rownames(sign_AMF_V1) %in% Secreted_prot_lin,]
  blue_Y_V4<-sign_AMF_V4[rownames(sign_AMF_V4) %in% Secreted_prot_lin,]
  blue_Y_V5<-sign_AMF_V5[rownames(sign_AMF_V5) %in% Secreted_prot_lin,]
  blue_Y_V6<-sign_AMF_V6[rownames(sign_AMF_V6) %in% Secreted_prot_lin,]
  blue_Y_V8<-sign_AMF_V8[rownames(sign_AMF_V8) %in% Secreted_prot_lin,]
  
  pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/topgenes_inmodules_expression/',
            sapply(strsplit(names(top_gene_sign_intramod_conec_A[i]),"_to_"), "[[", 1),".pdf",sep=""),width=10,height=10,useDingbats = FALSE)
  a_mod<-data.frame(blue_Y_V1,blue_Y_V4,blue_Y_V5,blue_Y_V6,blue_Y_V8)[,c(1,7,13,19,25)]
  colnames(a_mod)<-c("CAN_B1_V1","CAN_B1_V4","CAN_B1_V5","CAN_B1_V6","CAN_B1_V8")
  pheatmap(a_mod,
           main = sapply(strsplit(names(top_gene_sign_intramod_conec_A[i]),"_to_"), "[[", 1) )
  dev.off()
}


##########
# N6 Method of detection of secreted proteins by tang gigaspora rosa 2016
Secreted_prot_N6<-c(c("g11362.t1","g15507.t1","g13250.t1","g14466.t1","g13767.t1","g7260.t1","g7262.t1","g11437.t1","g11578.t1","g9257.t1","g10904.t1","g7858.t1","g12784.t1","g15330.t1","g714.t1","g15229.t1","g11374.t1","g12456.t1","g13394.t1","g13894.t1","g14544.t1","g4267.t1","g15125.t1","g14888.t1","g13396.t1","g10893.t1","g12723.t1","g14974.t1","g15300.t1","g2703.t1","g11498.t1","g11876.t1","g2878.t1","g8603.t1","g10372.t1","g14364.t1","g12556.t1","g14311.t1","g6636.t1","g14649.t1","g9927.t1","g14248.t1","g4325.t1","g11209.t1","g11157.t1","g960.t1","g13932.t1","g11111.t1","g13171.t1","g13869.t1","g9646.t1","g10343.t1","g10799.t1","g12782.t1","g2113.t1","g4296.t1","g14000.t1","g6033.t1","g6460.t1","g5208.t1","g11828.t1","g864.t1","g13208.t1","g10363.t1","g13857.t1","g14942.t1","g10419.t1","g1157.t1","g3750.t1","g8495.t1","g10782.t1","g14201.t1","g6637.t1","g9683.t1","g5944.t1","g9873.t1","g13985.t1","g6726.t1","g9904.t1","g12511.t1","g13656.t1","g12490.t1","g12751.t1","g10417.t1","g11964.t1","g10211.t1","g3001.t1","g13690.t1","g7006.t1","g8170.t1","g11892.t1","g1073.t1","g11335.t1","g5274.t1","g7442.t1","g9130.t1","g11930.t1","g5516.t1","g6251.t1","g11334.t1","g5583.t1","g13590.t1","g11780.t1","g4837.t1","g8009.t1","g1900.t1","g4244.t1","g11779.t1","g11199.t1","g6307.t1","g9806.t1","g13010.t1","g917.t1","g10651.t1","g14994.t1","g3907.t1","g7108.t1","g10494.t1","g14107.t1","g14989.t1","g1668.t1","g1678.t1","g4704.t1","g2172.t1","g12321.t1","g9074.t1","g12407.t1","g8921.t1","g9601.t1","g13788.t1","g15733.t1","g2195.t1","g2892.t1","g4305.t1","g5964.t1","g6308.t1","g5587.t1","g6204.t1","g9857.t1","g8989.t1","g12093.t1","g13285.t1","g918.t1","g431.t1","g3514.t1","g5907.t1","g14270.t1","g4253.t1","g4516.t1","g9129.t1","g10533.t1","g12356.t1","g7668.t1","g8821.t1","g3247.t1","g13682.t1","g77.t1","g1696.t1","g13947.t1","g3134.t1","g12427.t1","g1333.t1","g1513.t1","g10536.t1","g14947.t1","g5893.t1","g6225.t1","g14499.t1","g11168.t1","g10228.t1","g1555.t1","g1985.t1","g3639.t1","g5213.t1","g8684.t1","g12391.t1","g13083.t1","g4797.t1","g6559.t1","g11554.t1","g15913.t1","g2342.t1","g8865.t1","g12953.t1","g13754.t1","g5032.t1","g12864.t1","g5993.t1","g7589.t1","g915.t1","g2085.t1","g12792.t1","g13237.t1","g1101.t1","g8929.t1","g6326.t1","g2764.t1","g5990.t1","g3288.t1","g4474.t1"),
                  c("g6328.t1","g7857.t1","g10991.t1","g6131.t1","g6133.t1","g7478.t1","g12071.t1","g12809.t1","g8525.t1","g14198.t1","g8335.t1","g11981.t1","g14647.t1","g758.t1","g10425.t1","g12028.t1","g8647.t1","g12312.t1","g4475.t1","g12948.t1","g1848.t1","g2392.t1","g6035.t1","g9206.t1","g10577.t1","g11605.t1","g15725.t1","g7884.t1","g10504.t1","g4696.t1","g4964.t1","g5091.t1","g5585.t1","g6353.t1","g11386.t1","g4129.t1","g5591.t1","g3515.t1","g6365.t1","g6782.t1","g913.t1","g3321.t1","g14507.t1","g1633.t1","g4090.t1","g10059.t1","g11847.t1","g962.t1","g3885.t1","g7159.t1","g5572.t1","g6074.t1","g3290.t1","g4610.t1","g6639.t1","g1656.t1","g3767.t1","g7582.t1","g8213.t1","g9141.t1","g10133.t1","g1867.t1","g3564.t1","g5646.t1","g1268.t1","g7509.t1","g9868.t1","g13231.t1","g10308.t1","g11145.t1","g1142.t1","g3289.t1","g3910.t1","g3425.t1","g4002.t1","g5586.t1","g7816.t1","g9470.t1","g11593.t1","g4349.t1","g3822.t1","g8107.t1","g9041.t1","g7817.t1","g8475.t1","g12516.t1","g12939.t1","g3666.t1","g5910.t1","g2851.t1","g10879.t1","g2288.t1","g2719.t1","g3519.t1","g6932.t1","g10456.t1","g14668.t1","g1483.t1","g6529.t1","g6837.t1","g10520.t1","g10756.t1","g10926.t1","g6530.t1","g2883.t1","g866.t1","g2884.t1","g14029.t1","g309.t1","g1840.t1","g6634.t1","g9590.t1","g7230.t1","g9558.t1","g13960.t1","g4147.t1","g4702.t1","g6206.t1","g6526.t1","g1119.t1","g1225.t1","g10682.t1","g1318.t1","g6527.t1","g9815.t1","g1685.t1","g8557.t1","g8648.t1","g13445.t1","g11098.t1","g9799.t1","g10574.t1","g3606.t1","g3817.t1","g5498.t1","g5861.t1","g13041.t1","g12434.t1","g1471.t1","g10076.t1","g12114.t1","g2440.t1","g9323.t1","g7206.t1","g10798.t1","g1742.t1","g13286.t1","g6154.t1","g5010.t1","g14199.t1","g1730.t1","g4460.t1","g4849.t1","g7890.t1","g11039.t1","g7820.t1","g7286.t1","g6984.t1","g14410.t1","g4892.t1","g6151.t1","g6455.t1","g1456.t1","g1332.t1","g9681.t1","g5725.t1","g12169.t1","g2017.t1","g3599.t1","g5844.t1","g7076.t1","g662.t1","g7695.t1","g5673.t1","g1341.t1","g12967.t1","g14376.t1","g1836.t1","g2153.t1","g2943.t1","g9485.t1","g2774.t1","g829.t1","g9358.t1","g1264.t1","g6608.t1","g3447.t1","g7747.t1","g1631.t1","g451.t1","g2991.t1","g3786.t1","g6983.t1","g3307.t1","g9121.t1","g1574.t1","g7456.t1","g7484.t1","g14996.t1","g2374.t1","g7336.t1","g3592.t1","g8378.t1","g14082.t1"),
                  c("g3619.t1","g4403.t1","g6457.t1","g8532.t1","g1016.t1","g4546.t1","g7983.t1","g3053.t1","g4313.t1","g1473.t1","g3739.t1","g7905.t1","g11058.t1","g1104.t1","g11046.t1","g6990.t1","g9300.t1","g10307.t1","g12990.t1","g13244.t1","g3177.t1","g4712.t1","g4920.t1","g6295.t1","g10638.t1","g6389.t1","g15565.t1","g6908.t1","g7480.t1","g8349.t1","g13952.t1","g7109.t1","g9269.t1","g12073.t1","g4407.t1","g7754.t1","g9393.t1","g7584.t1","g577.t1","g13057.t1","g7471.t1","g9083.t1","g3292.t1","g10704.t2","g14159.t1","g3000.t1","g10575.t1","g10992.t1","g13864.t1","g12917.t1","g4359.t1","g1949.t1","g524.t1","g4420.t1","g7534.t1","g9096.t1","g10681.t1","g5651.t1","g15557.t1","g1334.t1","g9233.t1","g10895.t1","g8985.t1","g10680.t1","g3978.t1","g3381.t1","g3485.t1","g12709.t1","g3831.t1","g698.t1","g3293.t1","g5047.t1","g6090.t1","g8386.t1","g10628.t1","g1814.t1","g11622.t1","g2668.t1","g6465.t1","g1287.t1","g2230.t1","g10839.t1","g7950.t1","g1649.t1","g9082.t1","g5897.t1","g8204.t1","g10792.t1","g1193.t1","g9965.t1","g9661.t1","g2140.t1","g4330.t1","g7752.t1","g7987.t1","g12819.t1","g12924.t1","g10717.t1","g13912.t1","g6678.t1","g14243.t1","g3448.t1","g9798.t1","g272.t1","g1731.t1","g4353.t1","g5997.t1","g7157.t1","g789.t1","g5345.t1","g7971.t1","g4746.t1","g4351.t1","g11500.t1","g2765.t1","g10226.t1","g2557.t1","g10132.t1","g10699.t1","g10391.t1","g5831.t1","g11920.t1","g6864.t1","g14589.t1","g1733.t1","g6125.t1","g315.t1","g4597.t1","g8125.t1","g4931.t1","g12878.t1","g7521.t1","g9522.t1","g11054.t1","g3144.t1","g4650.t1","g6468.t1","g8554.t1","g10530.t1","g847.t2","g1654.t1","g1734.t1","g11608.t1","g2003.t1","g6964.t1","g9290.t1","g10280.t1","g12665.t1","g11670.t1","g2733.t1","g9632.t1","g10345.t1","g2009.t1","g15054.t1","g1729.t1","g9672.t1","g12238.t1","g7740.t1","g11520.t1","g8899.t1","g9912.t1","g8096.t1","g1136.t1","g8662.t1","g5015.t1","g9554.t1","g10809.t1","g11332.t1","g7382.t1","g6162.t1","g6390.t1","g12185.t1","g7719.t1","g12202.t1","g7237.t1","g10831.t1","g3960.t1","g7254.t1","g11728.t1","g11745.t1","g4917.t1","g6955.t1","g2008.t1","g12089.t1","g1241.t1","g3271.t1","g5705.t1","g9369.t1","g11052.t1","g11273.t1","g4616.t1","g4954.t1","g10164.t1","g7445.t1","g11026.t1","g12717.t1","g9220.t1","g14952.t1","g3060.t1","g9885.t1","g5958.t1","g11776.t1","g3781.t1"),
                  c("g4708.t1","g10453.t1","g2151.t1","g4918.t1","g6501.t1","g4926.t1","g12007.t1","g3041.t1","g6877.t1","g8016.t1","g3752.t1","g9018.t1","g14508.t1","g9017.t1","g9911.t1","g1922.t1","g2255.t1","g4089.t1","g893.t1","g2745.t1","g8594.t1","g3852.t1","g6960.t1","g7573.t1","g4922.t1","g13006.t1","g12315.t1","g7410.t1","g5665.t1","g3411.t1","g668.t1","g5612.t1","g2602.t1","g2894.t1","g15346.t1","g3246.t1","g11411.t1","g6505.t1","g469.t1","g470.t1","g6452.t1","g3372.t1","g11192.t1","g2605.t1","g3712.t1","g307.t1","g9828.t1","g11373.t1","g2763.t1","g4678.t1","g5601.t1","g5713.t1","g6893.t1","g3260.t1","g8931.t1","g10668.t1","g4007.t1","g4102.t1","g9855.t1","g11714.t1","g10440.t1","g2107.t1","g1612.t1","g2607.t1","g2523.t1","g9589.t1","g9609.t1","g497.t1","g2938.t1","g8866.t1","g13501.t1","g1871.t1","g4526.t1","g4472.t1","g2612.t1","g10776.t1","g14219.t1","g4484.t1","g5842.t1","g10167.t1","g11216.t1","g3198.t1","g10807.t1","g6662.t1","g536.t1","g3687.t1","g10790.t1","g13769.t1","g4411.t1","g9329.t1","g6519.t1","g8093.t1","g8584.t1","g2105.t1","g4485.t1","g3547.t1","g8999.t1","g3089.t1","g2347.t1","g2561.t1","g12274.t1","g5010.t2","g5592.t1","g6556.t1","g6271.t1","g6779.t1","g1934.t1","g966.t1","g5821.t1","g7602.t1","g12429.t1","g4707.t1","g4709.t1","g7200.t1","g11682.t1","g3645.t1","g9737.t1","g8498.t1","g9631.t1","g9856.t1","g11732.t1","g10627.t1","g7553.t1","g8114.t1","g4317.t1","g6884.t1","g7462.t1","g7516.t1","g10301.t1","g11414.t1","g555.t1","g5968.t1","g4383.t1","g4936.t1","g3977.t1","g4880.t1","g4841.t1","g828.t1","g11069.t1","g5626.t1","g2156.t1","g10493.t1","g11315.t1","g1683.t1","g2903.t1","g4921.t1","g12493.t1","g11521.t1","g12102.t1","g10778.t1","g541.t1","g5495.t1","g5016.t1","g12883.t1","g11478.t1","g2734.t1","g10571.t1","g10247.t1","g4682.t1","g816.t1","g2338.t1","g10369.t1","g10051.t1","g4197.t1","g2289.t1","g2158.t1","g1596.t1","g8279.t1","g12579.t1","g4268.t1","g7780.t1","g3099.t1","g7208.t1","g1489.t1","g7296.t1","g5078.t1","g6638.t1","g8955.t1","g7469.t1","g8243.t1","g10859.t1","g2603.t1","g12592.t1","g4934.t1","g12613.t1","g5718.t1","g11544.t1","g2895.t1","g7285.t1","g3640.t1","g2283.t1","g4796.t1","g3284.t1","g9823.t1","g6784.t1","g793.t1","g4527.t1","g3543.t1","g3294.t1","g7007.t1","g8659.t1","g1177.t1","g8412.t1","g5994.t1","g6325.t1"),
                  c("g12311.t1","g5992.t1","g5476.t1","g2351.t1","g5991.t1","g6327.t1","g8247.t1","g5937.t1","g2071.t1","g11381.t1","g6814.t1","g558.t1","g4130.t1","g1947.t1","g10853.t1","g1945.t1","g11542.t1","g1941.t1","g7187.t1","g8536.t1","g10258.t1","g11845.t1","g5972.t1","g1624.t1"))

Secreted_prot_N6_info<-Mercator_Rirregularis[Mercator_Rirregularis$gene %in% Secreted_prot_N6,]
table(Secreted_prot_N6_info$Function)

function_gene<-Secreted_prot_N6_info$gene[grep("Galactose",Secreted_prot_N6_info$details)]


blue_Y_V1<-sign_AMF_V1[rownames(sign_AMF_V1) %in% function_gene,]
blue_Y_V4<-sign_AMF_V4[rownames(sign_AMF_V4) %in% function_gene,]
blue_Y_V5<-sign_AMF_V5[rownames(sign_AMF_V5) %in% function_gene,]
blue_Y_V6<-sign_AMF_V6[rownames(sign_AMF_V6) %in% function_gene,]
blue_Y_V8<-sign_AMF_V8[rownames(sign_AMF_V8) %in% function_gene,]

pdf(paste('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/topgenes_inmodules_expression/',
          sapply(strsplit(names(top_gene_sign_intramod_conec_A[i]),"_to_"), "[[", 1),".pdf",sep=""),width=10,height=10,useDingbats = FALSE)
a_mod<-data.frame(blue_Y_V1,blue_Y_V4,blue_Y_V5,blue_Y_V6,blue_Y_V8)[,c(1,7,13,19,25)]
colnames(a_mod)<-c("CAN_B1_V1","CAN_B1_V4","CAN_B1_V5","CAN_B1_V6","CAN_B1_V8")
pheatmap(a_mod,
         main = sapply(strsplit(names(top_gene_sign_intramod_conec_A[i]),"_to_"), "[[", 1) )

Topgenes_AMF[Topgenes_AMF %in% Secreted_prot_N6]



###########################################################
################## TEST IMPORTANT GENES IN AMF ############## 
###########################################################

###### Cell organisation cassava
Fred_YUCA_genes<-  c("Manes.06G156700","Manes.12G124500","Manes.12G124500","Manes.01G123300",
                     "Manes.04G002400","Manes.12G124300","Manes.05G053600","Manes.08G154600",
                     "Manes.18G034100","Manes.16G125700","Manes.06G143100","Manes.01G071000","Manes.01G193000")
NORM_EXP_TOPAMF<-t(DGE_AMF_N$E[rownames(DGE_AMF_N$E) %in% amf ,])
# cassava genes
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in% Fred_YUCA_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]
cor2plot<-cor(data.frame(NORM_EXP_TOPAMF,NORM_EXP_TOPYUCA))
cor2plot<-cor2plot[grep("Manes",colnames(cor2plot)),]
cor2plot<-cor2plot[,grep(".t1",colnames(cor2plot))]
#pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_cellorg_amf.pdf',width=32,height=15,useDingbats = FALSE)
target_function<-Mercator_Rirregularis[Mercator_Rirregularis$gene %in% amf ,]
target_function<-target_function[!duplicated(target_function$gene),]
target_function2<-Mercator_Mesculenta[Mercator_Mesculenta$gene %in% cas ,]
target_function2<-target_function2[!duplicated(target_function2$gene),]
pheatmap(cor2plot)
pheatmap(cor2plot ,cellwidth = 8,cellheight = 8,labels_col = paste (colnames(cor2plot),target_function$Function[order(match(target_function$gene,colnames(cor2plot)))],sep="___"),
         labels_row = paste (colnames(t(cor2plot)),target_function2$Function[order(match(target_function2$gene,colnames(t(cor2plot))))],sep="___"))

dev.off()

# in cassava
NORM_EXP_TOPYUCA_ALL<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  cas ,])
NORM_EXP_TOPYUCA_ALL<-NORM_EXP_TOPYUCA_ALL[grep("CTRL",rownames(NORM_EXP_TOPYUCA_ALL),invert=T),]
NORM_EXP_TOPYUCA<-t(DGE_YUCA_N$E[rownames(DGE_YUCA_N$E) %in%  Fred_YUCA_genes ,])
NORM_EXP_TOPYUCA<-NORM_EXP_TOPYUCA[grep("CTRL",rownames(NORM_EXP_TOPYUCA),invert=T),]

cor2plot<-cor(data.frame(NORM_EXP_TOPYUCA,NORM_EXP_TOPYUCA_ALL))
cor2plot<-cor2plot[colnames(cor2plot) %in% Fred_YUCA_genes,]
cor2plot<-cor2plot[,grep("\\.1",colnames(cor2plot),invert=T)]
##pdf('~/Google Drive/Doctorat_shared_unil/RNA-seq/Results/MANUSCRIPT_NATURE_v1/V5/gene_gene_interaction_func/Cassava_cellorg_cassava.pdf',width=10,height=10,useDingbats = FALSE)
target_function<-Mercator_Mesculenta[Mercator_Mesculenta$gene %in% cas ,]
target_function<-target_function[!duplicated(target_function$gene),]
target_function2<-Mercator_Mesculenta[Mercator_Mesculenta$gene %in% cas ,]
target_function2<-target_function2[!duplicated(target_function2$gene),]
pheatmap(cor2plot ,cellwidth = 8,cellheight = 8,labels_col = paste (colnames(cor2plot),target_function$Function[order(match(target_function$gene,colnames(cor2plot)))],sep="___"),
         labels_row = paste (colnames(t(cor2plot)),target_function2$Function[order(match(target_function2$gene,colnames(t(cor2plot))))],sep="___"))
pheatmap(cor2plot )
cor2plot
dim(cor2plot)
for (i in 1:dim(cor2plot)[2]) { cor2plot2<-cor2plot[cor2plot[,i]>quantile(abs(cor2plot[,i]),.8),]
pheatmap( cor2plot2)
}
dev.off()



