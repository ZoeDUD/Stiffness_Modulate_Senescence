### GSE128033, GSE128169, GSE122960 import, follow by anchor transfer using GSE135893, then all combined into 1 seurat object
# ==============================================================================
library(Seurat)
library(Matrix)
library(tidyverse)
library(scCustomize)
library(dplyr)


### Arrange .gz files by different samples, each sample has one folder
### Change file name to barcodes, features, and matrix
#### Load GSE128033 data
setwd('/opt/home/buckcenter.org/hdu/GSE128033/G1280')
## create seurat object for all the mtx count in G1280 and G1281
# Keep all cells with at least 200 detected genes
dirs <- paste0('/opt/home/buckcenter.org/hdu/GSE128033/G1280/',list.dirs(path = '/opt/home/buckcenter.org/hdu/GSE128033/G1280', recursive = F, full.names = F))
dirs2 <- paste0('/opt/home/buckcenter.org/hdu/GSE128033/G1281/',list.dirs(path = '/opt/home/buckcenter.org/hdu/GSE128033/G1281', recursive = F, full.names = F))
h5dirs1 <- paste0('/opt/home/buckcenter.org/hdu/GSE128033/G1281/',list.files(path='/opt/home/buckcenter.org/hdu/GSE128033/G1281',pattern= '.h5', recursive = F, full.names = F))
h5dirs2 <- paste0('/opt/home/buckcenter.org/hdu/GSE122960/H5Filtered/',list.files(path='/opt/home/buckcenter.org/hdu/GSE122960/H5Filtered',pattern= 'filtered', recursive = F, full.names = F))

for(x in dirs){
  sp = str_split(x, '_')
  name = sp[[1]][2]
  cts <- ReadMtx(mtx = paste0(x,'/matrix.mtx.gz'),
                 features = paste0(x,'/features.tsv.gz'),
                 cells = paste0(x,'/barcodes.tsv.gz'))
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts,min.features = 200,
                                  names.field = 1,
                                  names.delim = "_"))
}
#setwd('/opt/home/buckcenter.org/hdu/GSE128033/G1281')
for(x in dirs2){
  sp = str_split(x, '_')
  name = sp[[1]][2]
  cts <- ReadMtx(mtx = paste0(x,'/matrix.mtx.gz'),
                 features = paste0(x,'/features.tsv.gz'),
                 cells = paste0(x,'/barcodes.tsv.gz'))
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts,min.features = 200,
                                  names.field = 1,
                                  names.delim = "_"))
}

for (x in h5dirs1){
  name = gsub ('raw_feature_bc_matrix.h5','',x)
  h5 <- Read10X_h5(filename = x,
                   use.names = TRUE,
                   unique.features = TRUE)
  assign(name, CreateSeuratObject(counts= h5$`Gene Expression`, min.features = 200))
}

for (x in h5dirs2){
  name = gsub ('_filtered_gene_bc_matrices_h5.h5','',x)
  h5 <- Read10X_h5(filename = x,
                   use.names = TRUE,
                   unique.features = TRUE)
  assign(name, CreateSeuratObject(counts= h5, min.genes = 200))
}


### form a list 
samples_list = mget(unique(c(ls(pattern='GSM'),ls(pattern='SC'))))


# check QC
# $$$$$$Filter out cells with less than 250 nFeature and more than 10% percent.mt
# Run SCTransform

SmCT<-function(l){
  l <- PercentageFeatureSet(l,pattern = "^MT-", col.name = "percent.mt")
  l <- subset(l, subset = nFeature_RNA > 250 & percent.mt < 10)
  l <- SCTransform(l, vst.flavor = "v2", verbose = TRUE)
} 

Sample.list.qc<-lapply(samples_list,SmCT)
#remove individual files 
rm(list=ls(pattern='GSM'))
rm(list=ls(pattern='SC'))

# Find variable features, PCA , UMAP
PC_M = function (l){
  l <- RunPCA(l, verbose = TRUE)
  l <- RunUMAP(l,dims = 1:15, verbose = TRUE)
}

Sample.list.PCA = lapply(Sample.list.qc, PC_M)

#####Merge all objects into one for anchor transfer
id_cell = names(Sample.list.PCA)
A_ild = merge(Sample.list.PCA$GSM3489182_Donor_01, y= c(Sample.list.PCA$GSM3489183_IPF_01,Sample.list.PCA$GSM3489184_IPF_02,
                                                        Sample.list.PCA$GSM3489185_Donor_02, Sample.list.PCA$GSM3489186_Cryobiopsy_01,
                                                        Sample.list.PCA$GSM3489187_Donor_03, Sample.list.PCA$GSM3489188_IPF_03,
                                                        Sample.list.PCA$GSM3489189_Donor_04, Sample.list.PCA$GSM3489190_IPF_04,
                                                        Sample.list.PCA$GSM3489191_Donor_05, Sample.list.PCA$GSM3489192_HP_01,
                                                        Sample.list.PCA$GSM3489193_Donor_06, Sample.list.PCA$`GSM3489194_SSc-ILD_01`,
                                                        Sample.list.PCA$GSM3489195_Donor_07, Sample.list.PCA$`GSM3489196_Myositis-ILD_01`,
                                                        Sample.list.PCA$GSM3489197_Donor_08, Sample.list.PCA$`GSM3489198_SSc-ILD_02`,
                                                        Sample.list.PCA$GSM3909673_SC277, Sample.list.PCA$GSM3909674_SC281,
                                                        Sample.list.PCA$GSM3909675_SC284, Sample.list.PCA$SC108SSCLOW,
                                                        Sample.list.PCA$SC109SSCUP, Sample.list.PCA$SC135SSCLOW,
                                                        Sample.list.PCA$SC136SSCUP, Sample.list.PCA$SC14NOR,
                                                        Sample.list.PCA$SC153IPFLOW, Sample.list.PCA$SC154IPFUP,
                                                        Sample.list.PCA$SC155NORLOW, Sample.list.PCA$SC156NORUP,
                                                        Sample.list.PCA$SC228NORbal, Sample.list.PCA$SC249NORbal,
                                                        Sample.list.PCA$SC31DNOR, Sample.list.PCA$SC31NOR,
                                                        Sample.list.PCA$SC45NOR, Sample.list.PCA$SC51SSCLOW,
                                                        Sample.list.PCA$SC52SSCUP, Sample.list.PCA$SC56NOR,
                                                        Sample.list.PCA$SC59NOR, Sample.list.PCA$SC63SSCLOW,
                                                        Sample.list.PCA$SC64SSCUP, Sample.list.PCA$SC87IPFLOW,
                                                        Sample.list.PCA$SC88IPFUP, Sample.list.PCA$SC89IPFLOW,
                                                        Sample.list.PCA$SC93IPFLOW, Sample.list.PCA$SC94IPFUP,
                                                        Sample.list.PCA$SC95IPFLOW), 
              add.cell.ids = id_cell,merge.data = TRUE,project = "SeuratProject")

# load the GSE135893 R file and add SCTModel  
GSE135893_ILD_annotated_fullsize <- readRDS("/opt/home/buckcenter.org/hdu/GSE135893_ILD_annotated_fullsize.rds")
GSE135893_ILD_annotated_fullsize@assays$SCT
GSE135893_ILD_annotated_fullsize = SCTransform(object =GSE135893_ILD_annotated_fullsize , vst.flavor = "v2",verbose = T)
# GSE135893_ILD_annotated_fullsize <- FindVariableFeatures(GSE135893_ILD_annotated_fullsize, verbose = T, nfeatures = 3000)
# GSE135893_ILD_annotated_fullsize <- ScaleData(GSE135893_ILD_annotated_fullsize, features = row.names(GSE135893_ILD_annotated_fullsize@assays$SCT@data))
# GSE135893_ILD_annotated_fullsize <- RunPCA(GSE135893_ILD_annotated_fullsize)
# ElbowPlot(GSE135893_ILD_annotated_fullsize)

# Transfer Anchors using GSE135893 reference 
anchors <- FindTransferAnchors(reference = GSE135893_ILD_annotated_fullsize,
                               query = A_ild,
                               normalization.method = 'SCT',
                               reference.assay='SCT',
                               query.assay='SCT',
                               dims = 1:15)
unique(GSE135893_ILD_annotated_fullsize@meta.data$celltype)

predictions <- TransferData(anchorset = anchors,
                            refdata = GSE135893_ILD_annotated_fullsize@meta.data$celltype,
                            verbose = T)
A_ild <- AddMetaData(A_ild, metadata = predictions)
A_ild <- subset(A_ild, percent.mt < 10)
#A_ild@meta.data[8:length(colnames(A_ild@meta.data))] <- NULL

# Filter low id-score cell 
S_ild = subset(A_ild, prediction.score.max > 0.5)
S_ild = subset(S_ild, percent.mt < 10)
S_ild@meta.data[8:length(colnames(S_ild@meta.data))] <- NULL

## add study, status and diagnosis label  
S_ild@meta.data$Study = 'G'
S_ild@meta.data$Status = 'ILD'
S_ild@meta.data$Diagnosis = 'Control'
#
S_ild@meta.data[grepl('GSM3489182',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489182',rownames(S_ild@meta.data)),9]='Control'
S_ild@meta.data[grepl('GSM3489182',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3489183',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489183',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489183',rownames(S_ild@meta.data)),10]='IPF' 
#
S_ild@meta.data[grepl('GSM3489184',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489184',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489184',rownames(S_ild@meta.data)),10]='IPF' 
#
S_ild@meta.data[grepl('GSM3489185',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489185',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('GSM3489185',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3489186',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489186',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489186',rownames(S_ild@meta.data)),10]='IPF' 
#
S_ild@meta.data[grepl('GSM3489187',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489187',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('GSM3489187',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3489188',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489188',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489188',rownames(S_ild@meta.data)),10]='IPF' 
#
S_ild@meta.data[grepl('GSM3489189',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489189',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('GSM3489189',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3489190',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489190',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489190',rownames(S_ild@meta.data)),10]='IPF' 
#
S_ild@meta.data[grepl('GSM3489191',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489191',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('GSM3489191',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3489192',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489192',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489192',rownames(S_ild@meta.data)),10]='HP' 
#
S_ild@meta.data[grepl('GSM3489193',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489193',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('GSM3489193',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3489194',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489194',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489194',rownames(S_ild@meta.data)),10]='SSc' 
#
S_ild@meta.data[grepl('GSM3489195',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489195',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('GSM3489195',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3489196',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489196',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489196',rownames(S_ild@meta.data)),10]='PM' 
#
S_ild@meta.data[grepl('GSM3489197',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489197',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('GSM3489197',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3489198',rownames(S_ild@meta.data)),8]='GSE122960' 
S_ild@meta.data[grepl('GSM3489198',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3489198',rownames(S_ild@meta.data)),10]='SSc' 
#
S_ild@meta.data[grepl('GSM3909673',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('GSM3909673',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('GSM3909673',rownames(S_ild@meta.data)),10]='Control' 
#
S_ild@meta.data[grepl('GSM3909674',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('GSM3909674',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3909674',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('GSM3909675',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('GSM3909675',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('GSM3909675',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC108SSCLOW',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('SC108SSCLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC108SSCLOW',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC109SSCUP',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('SC109SSCUP',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC109SSCUP',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC135SSCLOW',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('SC135SSCLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC135SSCLOW',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC136SSCUP',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('SC136SSCUP',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC136SSCUP',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC14NOR',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC14NOR',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC14NOR',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC153IPFLOW',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC153IPFLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC153IPFLOW',rownames(S_ild@meta.data)),10]='IPF'
#
S_ild@meta.data[grepl('SC154IPFUP',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC154IPFUP',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC154IPFUP',rownames(S_ild@meta.data)),10]='IPF'
#
S_ild@meta.data[grepl('SC155NORLOW',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC155NORLOW',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC155NORLOW',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC156NORUP',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC156NORUP',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC156NORUP',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC228NORbal',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC228NORbal',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC228NORbal',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC249NORbal',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC249NORbal',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC249NORbal',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC31DNOR',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC31DNOR',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC31DNOR',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC31NOR',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC31NOR',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC31NOR',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC45NOR',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC45NOR',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC45NOR',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC51SSCLOW',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('SC51SSCLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC51SSCLOW',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC52SSCUP',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('SC52SSCUP',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC52SSCUP',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC56NOR',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC56NOR',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC56NOR',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC59NOR',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC59NOR',rownames(S_ild@meta.data)),9]='Control' 
S_ild@meta.data[grepl('SC59NOR',rownames(S_ild@meta.data)),10]='Control'
#
S_ild@meta.data[grepl('SC63SSCLOW',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('SC63SSCLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC63SSCLOW',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC64SSCUP',rownames(S_ild@meta.data)),8]='GSE128169' 
S_ild@meta.data[grepl('SC64SSCUP',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC64SSCUP',rownames(S_ild@meta.data)),10]='SSc'
#
S_ild@meta.data[grepl('SC87IPFLOW',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC87IPFLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC87IPFLOW',rownames(S_ild@meta.data)),10]='IPF'
#
S_ild@meta.data[grepl('SC88IPFUP',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC88IPFUP',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC88IPFUP',rownames(S_ild@meta.data)),10]='IPF'
#
S_ild@meta.data[grepl('SC89IPFLOW',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC89IPFLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC89IPFLOW',rownames(S_ild@meta.data)),10]='IPF'
#
S_ild@meta.data[grepl('SC93IPFLOW',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC93IPFLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC93IPFLOW',rownames(S_ild@meta.data)),10]='IPF'
#
S_ild@meta.data[grepl('SC94IPFUP',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC94IPFUP',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC94IPFUP',rownames(S_ild@meta.data)),10]='IPF'
#
S_ild@meta.data[grepl('SC95IPFLOW',rownames(S_ild@meta.data)),8]='GSE128033' 
S_ild@meta.data[grepl('SC95IPFLOW',rownames(S_ild@meta.data)),9]='ILD' 
S_ild@meta.data[grepl('SC95IPFLOW',rownames(S_ild@meta.data)),10]='IPF'

S_ild@meta.data = S_ild@meta.data %>% rename (celltype=predicted.id)
GSE135893_ILD_annotated_fullsize $Study = 'GSE135893'

### merge A_ild with GSE135893 data
lung = merge(GSE135893_ILD_annotated_fullsize, y=S_ild)



# Graph for supplementary 
# ==============================================================================
lung = RunPCA(lung)
lung = RunUMAP(lung,dims = 1:20, verbose = T)
DimPlot(lung, reduction = "umap",group.by = 'celltype',label = T,label.size = 2.5)
ggsave("celltypeUMAP.pdf",width=15,height=9.5)

DimPlot(lung, reduction = "umap",group.by = 'Study',label = T,label.size = 2.5)
ggsave("StudyA_UMAP.pdf",width=13,height=9.5)


#### population portion table
tot_cell <- prop.table(table(lung@meta.data$celltype))
tot_cell = (tot_cell*100)

G1229 = filter (lung@meta.data, Study == 'GSE122960')
G1280 = filter (lung@meta.data, Study == 'GSE128033')
G1281 = filter (lung@meta.data, Study == 'GSE128169')
G1358 = filter (lung@meta.data, Study == 'GSE135893')

G1229=(prop.table(table(G1229$celltype))*100)
G1280=(prop.table(table(G1280$celltype))*100)
G1281=(prop.table(table(G1281$celltype))*100)
G1358=(prop.table(table(G1358$celltype))*100)

comparison <- cbind(tot_cell, G1229,G1280,G1281,G1358)
write.csv(comparison, "CellportionbyGSE.csv", quote = F)

