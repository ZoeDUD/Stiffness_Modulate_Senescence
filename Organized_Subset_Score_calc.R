
# ==========================================================================================
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)

#save.image('fibro_Object.RData')
load("fibro_Object.RData")

### subset all fibroblast cell types 
unique(lung$celltype)
keep <-  c("Myofibroblasts",
           "Fibroblasts",
           "PLIN2+ Fibroblasts",
           'HAS1 High Fibroblasts')

fibro <- subset(lung, cells = rownames(lung@meta.data[lung@meta.data$celltype %in% keep, ]))

## check the data
unique(fibro@meta.data$celltype)
unique(fibro$Study)


### look for DEG 
# DEG bw control v. ILD
ILDvC_DEG=FindMarkers(fibro,assay = "SCT",
                      group.by = "Status",
                      ident.1 = "ILD",
                      ident.2 = "Control",
                      test.use = "LR",latent.vars ='Study',
                      min.pct=0.005, logfc.threshold = 0.1)
ILDvC_DEG$Gene_name=row.names(ILDvC_DEG)
write.table(ILDvC_DEG,file="AllFibro_ILDvC.txt",quote=F,sep="\t",row.names=F,col.names = T)
GSEA = ILDvC_DEG[,c(6,2)]
attach(GSEA)
GSEA = GSEA[order(-avg_log2FC),]
write.table(GSEA,file="AllFibro_ILDvC.rnk",quote=F,sep="\t",row.names=F,col.names = F)

### input gene sets (Supplementary table 6)
geneSet = read.csv ('/opt/home/buckcenter.org/hdu/SeneGeneSet.csv',header = T)

csGO = geneSet$GOBP_CELLULAR_SENESCENCE[2:92]
fridman = geneSet$FRIDMAN_SENESCENCE_UP[2:78]
SenMayo = geneSet$SenMayo[2:126]
SFS2017 = read.csv('/opt/home/buckcenter.org/hdu/SenescenceFibroblastSignature2017.csv', header=T) %>%
  arrange(-Meta_NB.GLM_FC)
SFS2017 = SFS2017[c(1:727),1]
IMRSASP = geneSet$Schilling_IMR90_SASP[2:151]
CoreSASP = geneSet$Schilling_Core_SASP[2:26]

Ionchan= geneSet$GOMF_MECHANOSENSITIVE_ION_CHANNEL_ACTIVITY[2:18]
CresM = geneSet$GOBP_CELLULAR_RESPONSE_TO_MECHANICAL_STIMULUS[2:71]
ECM = geneSet$GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY[2:47]
KG_FA = geneSet$KEGG_FOCAL_ADHESION[2:200]

X250PL = geneSet$X250PLgene[2:23]
X250Gel = geneSet$X250Gelgene[2:13]
MS250a2 = geneSet$MS_PLv2[2:31]
MS250a50 = geneSet$MS_PLv50[2:31]
MS_PO_250v2 = c('PLEC','COL6A3','FLNA','FLNC','FLNB','AHNAK','MYH9','FN1','TLN1','ACTN1','DSP','SPTAN1','IQGAP1','KRT10','COL12A1','VCL',
                'DYNC1H1','SPTBN1','HSPG2','COL1A2','KRT1','THBS1','LAMA1','MAP4','LMNA','COL1A1','LAMA4','LAMB1','CLTC','MAP1B')
MS50a2 = c('PTMA','AP3B1','DCBLD2','CHORDC1','YKT6','LGALS7','GPNMB','CCT8','SDC4','NLN','FUCA1','TXNDC12','TPT1','MAP1LC3A',
           'LAMP2','PROCR','SEC13','CMPK1','RPS7','IDI1','TARS1','CALML5','S100A16','CPNE3','NCSTN','NME2P1','TMSB4X','ECI1',
           'PAFAH1B2','PIK3IP1')
collgen = c('COL1A1','COL1A2','COL2A1','COL3A1','COL4A1','COL4A2','COL4A3','COL4A4','COL4A5','COL4A6','COL5A1',
            'COL5A2','COL5A3','COL6A1','COL6A2','COL6A3','COL6A4P1','COL6A4P2','COL6A5','COL6A6','COL7A1','COL8A1',
            'COL8A2','COL9A1','COL9A2','COL9A3','COL10A1','COL11A1','COL11A2','COL12A1','COL13A1','COL14A1','COL15A1',
            'COL16A1','COL17A1','COL18A1','COL19A1','COL20A1','COL21A1','COL22A1','COL23A1','COL24A1','COL25A1','COL26A1',
            'COL27A1','COL28A1')

##################
#### Add module for score 
DefaultAssay(object = fibro) <- "RNA"
fibro= NormalizeData(object = fibro)
fibro = AddModuleScore (fibro, features = list(csGO), 
                        name = 'GOBPsene')
fibro = AddModuleScore (fibro, features = list(fridman), 
                        name = 'Fridman_up')
fibro = AddModuleScore (fibro, features = list(SFS2017), 
                        name = 'seneCB2017')
fibro = AddModuleScore (fibro, features = list(SenMayo), 
                        name = 'SenMayo')
fibro = AddModuleScore (fibro, features = list(IMRSASP), 
                        name = 'SASPatlas_IMR')
fibro = AddModuleScore (fibro, features = list(CoreSASP), 
                        name = 'SASPatlas_core')

fibro = AddModuleScore (fibro, features = list(Ionchan), 
                        name = 'Ionchan')
fibro = AddModuleScore (fibro, features = list(CresM), 
                        name = 'CelltoMech')
fibro = AddModuleScore (fibro, features = list(ECM), 
                        name = 'ECMassmble')
fibro = AddModuleScore (fibro, features = list(KG_FA), 
                        name = 'KEGG_FocalA')

fibro = AddModuleScore (fibro, features = list(X250PLECM), 
                        name = 'X250nmPLecm')
fibro = AddModuleScore (fibro, features = list(X250PL), 
                        name = 'X250nmPL_DEG')
fibro = AddModuleScore (fibro, features = list(X250Gel), 
                        name = 'X250nm2kpa_DEG')
fibro = AddModuleScore (fibro, features = list(MS250a2), 
                        name = 'X250nmPLv2k_MS')
fibro = AddModuleScore (fibro, features = list(MS250a50), 
                        name = 'X250nmPLv50k_MS')
fibro = AddModuleScore (fibro, features = list(MS_PO_250v2),
                        name = 'PoMS_250nmPLv2k')
fibro = AddModuleScore (fibro, feature = list(MS50a2),name ='X250nm50v2k_MS')

fibro = AddModuleScore (fibro, features = list(collgen), 
                        name = 'CollagenFam')

## plot box plot
library(dplyr)
library(ggpubr)
#### senescent gene sets  
t1 = c(filter(fibro@meta.data,Status == 'Control')$GOBPsene1, filter(fibro@meta.data,Status == 'ILD')$GOBPsene1,
       filter(fibro@meta.data,Status == 'Control')$Fridman_up1, filter(fibro@meta.data,Status == 'ILD')$Fridman_up1,
       filter(fibro@meta.data,Status == 'Control')$seneCB20171, filter(fibro@meta.data,Status == 'ILD')$seneCB20171,
       filter(fibro@meta.data,Status == 'Control')$SenMayo1,filter(fibro@meta.data,Status == 'ILD')$SenMayo1,
       filter(fibro@meta.data,Status == 'Control')$SASPatlas_IMR1,filter(fibro@meta.data,Status == 'ILD')$SASPatlas_IMR1,
       filter(fibro@meta.data,Status == 'Control')$SASPatlas_core1,filter(fibro@meta.data,Status == 'ILD')$SASPatlas_core1)
t2 = c(rep(c(rep('Control',731),rep('ILD',6844)),6))
t3= c(rep('GOBP_cellular_senescence',7575),rep('Fridman_Tainsky_2008',7575),rep('Hernandez-Segura_et_al_2017',7575),
      rep('Saul_et_al_2022_SenMayo',7575),rep('Basisty_et_al_2020_SASPAtlas_IMR90',7575),rep('Basisty_et_al_2020_SASPAtlas_Core',7575))

data = data.frame(ModuleScore=t1, Group=t2, GeneSet=t3)
ggboxplot(data, x = "GeneSet", y = "ModuleScore",
          color = "Group", palette = "jco",outlier.shape = NA,ylim =c(-1, 15) ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("SeneSetScoreBOX.pdf", width = 14.8, height = 8)

#### Mechano gene sets  
t1 = c(filter(fibro@meta.data,Status == 'Control')$KEGG_FocalA1,filter(fibro@meta.data,Status == 'ILD')$KEGG_FocalA1,
       filter(fibro@meta.data,Status == 'Control')$Ionchan1,filter(fibro@meta.data,Status == 'ILD')$Ionchan1,
       filter(fibro@meta.data,Status == 'Control')$CelltoMech1,filter(fibro@meta.data,Status == 'ILD')$CelltoMech1,
       filter(fibro@meta.data,Status == 'Control')$ECMassmble1,filter(fibro@meta.data,Status == 'ILD')$ECMassmble1
       )

length(filter(fibro@meta.data,Status == 'Control')$KEGG_FocalA1) #731
length(filter(fibro@meta.data,Status == 'ILD')$KEGG_FocalA1) #6844
t2 = c(rep(c(rep('Control',731),rep('ILD',6844)),4))
t3= c(rep('KEGG_focal_adhesion', 7575),rep('GOMF_mechanosensitive_ion_channel_activity',7575), rep('GOBP_cellular_response_to_mechanical_stimulus',7575),rep('GOBP_extracellular_matrix_assembly
',7575))

data = data.frame(ModuleScore=t1, Group=t2, GeneSet=t3)
ggboxplot(data, x = "GeneSet", y = "ModuleScore",
          color = "Group", palette = "jco",outlier.shape = NA,ylim =c(-0.5, 2.5) ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("MechanoSetScoreBOX.pdf", width = 14.8, height = 8)


##plot PL V2 DEG score (ILD V CON)
p_df1 = c(filter(fibro@meta.data,Status == 'Control')$X250nmPL_DEG1,filter(fibro@meta.data,Status == 'ILD')$X250nmPL_DEG1)
p_df2 =c(rep('Control',731),rep('ILD',6844))
p_df3= c(rep('250nMPLv2kPa_DEG',7575))
p_df = data.frame(ModuleScore=p_df1, Group=p_df2, GeneSet=p_df3)
ggboxplot(p_df, x = "GeneSet", y = "ModuleScore",
          color = "Group", palette = "jco",outlier.shape = NA,ylim =c(-3, 18) ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('plv2DEGScoreBOX.pdf', width = 7, height = 8)

wilcox.test(fibro$X250nmPL_DEG1[grep('ILD',fibro$Status)],y= fibro$X250nmPL_DEG1[grep('Control',fibro$Status)])

## Plot pool MS score 
p_df1 = c(filter(fibro@meta.data,Status == 'Control')$PoMS_250nmPLv2k1,filter(fibro@meta.data,Status == 'ILD')$PoMS_250nmPLv2k1)
p_df2 =c(rep('Control',731),rep('ILD',6844))
p_df3= c(rep('Pooled250nMPLv2kPa_SASP',7575))
p_df = data.frame(ModuleScore=p_df1, Group=p_df2, GeneSet=p_df3)
ggboxplot(p_df, x = "GeneSet", y = "ModuleScore",
          color = "Group", palette = "jco",outlier.shape = NA,ylim =c(-3, 9.5) ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('PoolMsScoreBOX.pdf', width = 7, height = 8)

wilcox.test(fibro$PoMS_250nmPLv2k1[grep('ILD',fibro$Status)],y= fibro$PoMS_250nmPLv2k1[grep('Control',fibro$Status)])

## Plot 50 V 2 (ILD V CON)
p_df1 = c(filter(fibro@meta.data,Status == 'Control')$X250nm50v2k_MS1,filter(fibro@meta.data,Status == 'ILD')$X250nm50v2k_MS1)
p_df2 =c(rep('Control',731),rep('ILD',6844))
p_df3= c(rep('250nM50v2kPa_MS',7575))
p_df = data.frame(ModuleScore=p_df1, Group=p_df2, GeneSet=p_df3)
ggboxplot(p_df, x = "GeneSet", y = "ModuleScore",
          color = "Group", palette = "jco",outlier.shape = NA,ylim =c(-3, 11) ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('50v2_MsScoreBOX.pdf', width = 7, height = 8)

wilcox.test(fibro$X250nm50v2k_MS1[grep('ILD',fibro$Status)],y= fibro$X250nm50v2k_MS1[grep('Control',fibro$Status)])


# #### plot coorelation of module score
# library(psych)
# df= fibro@meta.data[,15:38]
# pdf("ModuleScore_corplotR.pdf",width=19, height=7)
# corPlot(df,stars = TRUE, gr = colorRampPalette(c('blue','white','red')))
# dev.off()
# ### plot the selected sene and mechno set coorelation
# df= fibro@meta.data[,c(15,20,16,18,19,21,22,25,27,32,35:38)]
# pdf("SeleModuleScore_corplotR.pdf",width=10, height=7)
# corPlot(df,stars = TRUE, gr = colorRampPalette(c('blue','white','red')),upper = FALSE,xsrt=30)
# dev.off()

#### This package looks better
library(corrplot)
#library(dplyr)
df= fibro@meta.data[,c(16,18,19,20,21,22,32,24,25,27,35,36,39,40,38)]
colnames(df)= c('GOBP_cellular_senescence','Fridman_Tainsky_2008','Hernandez-Segura_et_al_2017','Saul_et_al_2022_SenMayo',
                'Basisty_2020_SASPAtlasIMR90','Basisty_2020_SASPAtlasCore','KEGG_focal_adhesion','GOMF_mechanosensitive_ion_channel_activity','GOBP_cellular_response_to_mechanical_stimulus',
                'GOBP_extracellular_matrix_assembly','250nMdoxo_PL_mRNA_DEG','250nMdoxo_2kPa_mRNA_DEG','250nMdoxo_PLv2kPa_MS_Cb_SASP','250nMdoxo_50v2kPa_MS_SASP','250nMdoxo_PLv50kPa_MS_SASP')
M = cor(df)
testRes = cor.mtest(df, conf.level = 0.95)
pAdj <- p.adjust(c(testRes[[1]]), method = "BH")
resAdj <- matrix(pAdj, ncol = dim(testRes[[1]])[1])
colnames(resAdj) = colnames(testRes$p)
rownames(resAdj) = rownames(testRes$p)
# corrplot(M,type = 'lower', order = 'hclust', tl.col = 'black',
#          cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))
# corrplot(M)
# selected bigger plot
pdf("ModuleScore_corplot_0815.2.pdf",width=9, height=7.6)
corrplot(M, p.mat = resAdj, method = 'circle', type = 'lower', insig='blank', 
         addCoef.col ='black', number.cex = 0.8, diag=FALSE,tl.col='black',tl.srt = 45,tl.cex=0.6, col=colorRampPalette(c("darkslateblue","white","firebrick3"))(200))
dev.off()

# selected smaller plot
df= fibro@meta.data[,c(20,16,18,19,21,22,25,27,32)]
colnames(df)= c('Saul_et_al_2022_SenMayo','GO:0090398','Fridman_Tainsky_2008','Hernandez-Segura_et_al_2017',
                'Basisty_2020_SASPAtlasIMR90','Basisty_2020_SASPAtlasCore','GO:0071260',
                'GO:0085029','hsa04510')
pdf("SeleModuleSenmayo_corsmall_half.pdf",width=7, height=7)
corrplot(M, p.mat = resAdj, method = 'circle', type = 'lower', insig='blank', 
         addCoef.col ='black', number.cex = 0.8, diag=FALSE,tl.srt = 45)
dev.off()
## smaller 2
df2 = df [,c(1:7,9:10)]
M2 = cor(df2)
testRes2 = cor.mtest(df2, conf.level = 0.95)
pAdj2 <- p.adjust(c(testRes2[[1]]), method = "BH")
resAdj2 <- matrix(pAdj2, ncol = dim(testRes2[[1]])[1])
colnames(resAdj2) = colnames(testRes2$p)
rownames(resAdj2) = rownames(testRes2$p)
pdf("ModuleScore_corplot2.pdf",width=9, height=7.6)
corrplot(M2, p.mat = resAdj2, method = 'circle', type = 'lower', insig='blank', 
         addCoef.col ='black', number.cex = 0.8, diag=FALSE,tl.col='black',tl.srt = 45,tl.cex=0.6, col=colorRampPalette(c("darkslateblue","white","firebrick3"))(200))
dev.off()




#######T test
wilcox.test(fibro$KEGG_FocalA1[grep('ILD',fibro$Status)],y= fibro$KEGG_FocalA1[grep('Control',fibro$Status)])
wilcox.test(fibro$Ionchan1[grep('ILD',fibro$Status)],y= fibro$Ionchan1[grep('Control',fibro$Status)])
wilcox.test(fibro$CelltoMech1[grep('ILD',fibro$Status)],y= fibro$CelltoMech1[grep('Control',fibro$Status)])
wilcox.test(fibro$ECMassmble1[grep('ILD',fibro$Status)],y= fibro$ECMassmble1[grep('Control',fibro$Status)])
wilcox.test(fibro$GOBPsene1[grep('ILD',fibro$Status)],y= fibro$GOBPsene1[grep('Control',fibro$Status)])
wilcox.test(fibro$Fridman_up1[grep('ILD',fibro$Status)],y= fibro$Fridman_up1[grep('Control',fibro$Status)])
wilcox.test(fibro$seneCB20171[grep('ILD',fibro$Status)],y= fibro$seneCB20171[grep('Control',fibro$Status)])
wilcox.test(fibro$SenMayo1[grep('ILD',fibro$Status)],y= fibro$SenMayo1[grep('Control',fibro$Status)])
wilcox.test(fibro$SASPatlas_IMR1[grep('ILD',fibro$Status)],y= fibro$SASPatlas_IMR1[grep('Control',fibro$Status)])
wilcox.test(fibro$SASPatlas_core1[grep('ILD',fibro$Status)],y= fibro$SASPatlas_core1[grep('Control',fibro$Status)])
wilcox.test(fibro$[grep('ILD',fibro$Status)],y= fibro$[grep('Control',fibro$Status)])
 

##plot 2 V PL DEG score (ILD V CON)
p_df1 = c(filter(fibro@meta.data,Status == 'Control')$X250nm2kpa_DEG1,filter(fibro@meta.data,Status == 'ILD')$X250nm2kpa_DEG1)
p_df2 =c(rep('Control',731),rep('ILD',6844))
p_df3= c(rep('X250nm2kpa_DEG1',7575))
p_df = data.frame(ModuleScore=p_df1, Group=p_df2, GeneSet=p_df3)
ggboxplot(p_df, x = "GeneSet", y = "ModuleScore",
          color = "Group", palette = "jco",outlier.shape = NA,ylim =c(-3, 18) ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('2kpDEGScoreBOX.pdf', width = 7, height = 8)
wilcox.test(fibro$X250nm2kpa_DEG1[grep('ILD',fibro$Status)],y= fibro$X250nm2kpa_DEG1[grep('Control',fibro$Status)])

# 
# gelset = c('TREM1','MMP3','PHPT1','TFPI2','MMP1','SVIP','STC1')
# placset = c('SPARC','CAV1','IGFBP5','ACTA2','LBH','DCN','PDE5A','TGFB2-OT1','FGF7','NREP','COL5A2','THBS2','COL1A1','WNT16','COL3A1',
#             'SOX6','MKI67')
# 
# fibro = AddModuleScore (fibro, features = list(gelset), 
#                         name = 'gelset')
# fibro = AddModuleScore (fibro, features = list(placset), 
#                         name = 'placset')
# 
# p_df1 = c(filter(fibro@meta.data,Status == 'Control')$gelset1,filter(fibro@meta.data,Status == 'ILD')$gelset1)
# p_df2 =c(rep('Control',731),rep('ILD',6844))
# p_df3= c(rep('placset',7575))
# p_df = data.frame(ModuleScore=p_df1, Group=p_df2, GeneSet=p_df3)
# ggboxplot(p_df, x = "GeneSet", y = "ModuleScore",
#           color = "Group", palette = "jco",outlier.shape = NA,ylim =c(-3, 18) ) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# wilcox.test(fibro$gelset1[grep('ILD',fibro$Status)],y= fibro$gelset1[grep('Control',fibro$Status)])
# ggsave('plv2DEGScoreBOX.pdf', width = 7, height = 8)
# 

