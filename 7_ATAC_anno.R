library(Seurat)
library(Signac)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("chromVAR")

library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
set.seed(2023)

library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)

library(motifmatchr)
library(Signac)
library(Seurat)
library(JASPAR2018)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(getMatrixSet)
#T2_K2_integrated3
scmulti3 = T2_K2_integrated3
# Get a list of motif position frequency matrices from the JASPAR database
# 使用getMatrixSet函数从JASPAR数据库中提取Motif的PFM矩阵信息
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
# 使用CreateMotifMatrix函数构建Motif矩阵对象
#scmulti4 = scmulti3
DefaultAssay(object =  scmulti3) <- 'peaks'
# add motif information
scmulti3 <- AddMotifs(
  object = scmulti3,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

pfmlist = as.matrix(pfm)
pfmlist = as.data.frame(pfmlist$V1$)

########
#Computing motif activities
scmulti3  <- RunChromVAR(
  object = scmulti3,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(scmulti3) <- 'chromvar'


library(rtracklayer)
# declare names and types for the extra BED columns
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

###H3K4me3
H3K4me3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF246IEW_H3K4me3.bed",extraCols=extraCols_narrowPeak)

head(H3K4me3)
H3K4me3@seqnames
write.table(H3K4me3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF246IEW_H3K4me3.txt",sep = "\t")

H3K4me3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF246IEW_H3K4me3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K4me3)){
  x=paste(H3K4me3[i,1],H3K4me3[i,2],H3K4me3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K4me3_features = list(list$V1)


ATAC_peak = as.data.frame(scmulti2@assays$peaks@data)
ATAC_peak = rownames(ATAC_peak)
ATAC_peak = as.data.frame(ATAC_peak)

x1 = vector()
y1 = vector()
z1 = vector()
for (i in 1:nrow(ATAC_peak)){
  x = strsplit(ATAC_peak[i,1],split='-')[[1]][1]
  y = strsplit(ATAC_peak[i,1],split='-')[[1]][2]
  z = strsplit(ATAC_peak[i,1],split='-')[[1]][3]
  x1 = c(x1,x)
  y1 = c(y1,y)
  z1 = c(z1,z)
}
ATAC_peak = cbind(ATAC_peak,x1,as.numeric(y1),as.numeric(z1))
colnames(ATAC_peak)[3]="start"
colnames(ATAC_peak)[4]="end"
write.table(ATAC_peak,file="/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_k562.txt",sep="\t")

ATAC_peak_H3K4me3_list = data.frame()
for(i in 1:nrow(H3K4me3)){
  ATAC_peak_H3K4me3 <- ATAC_peak[ATAC_peak$x1 == H3K4me3[i,1] & ((H3K4me3[i,2] > ATAC_peak$start & H3K4me3[i,2] < ATAC_peak$end) | (H3K4me3[i,3] > ATAC_peak$start & H3K4me3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K4me3) != 0){
    ATAC_peak_H3K4me3$V5 <- H3K4me3[i,1]
    ATAC_peak_H3K4me3$V6 <- H3K4me3[i,2]
    ATAC_peak_H3K4me3$V7 <- H3K4me3[i,3]
    ATAC_peak_H3K4me3$V8 <- H3K4me3[i,6]
    ATAC_peak_H3K4me3$V9 <- H3K4me3[i,7]
    ATAC_peak_H3K4me3_list = rbind(ATAC_peak_H3K4me3_list,ATAC_peak_H3K4me3)
  }
}
write.table(ATAC_peak_H3K4me3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K4me3_list.txt",sep="\t")

##
H3K4me3_features=list("H3K4me3" = ATAC_peak_H3K4me3_list$ATAC_peak)

scmulti4 = scmulti3 

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K4me3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K4me3'
)

write.table(ATAC_peak_H3K4me3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K4me3_list.txt",sep="\t")


##H3K27ac
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K27ac = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF532MMV_H3K27ac.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K27ac)
H3K27ac@seqnames
write.table(H3K27ac,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF532MMV_H3K27ac.txt",sep = "\t")

H3K27ac = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF532MMV_H3K27ac.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K27ac)){
  x=paste(H3K27ac[i,1],H3K27ac[i,2],H3K27ac[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K27ac_features = list(list$V1)

ATAC_peak_H3K27ac_list = data.frame()
for(i in 1:nrow(H3K27ac)){
  ATAC_peak_H3K27ac <- ATAC_peak[ATAC_peak$x1 == H3K27ac[i,1] & ((H3K27ac[i,2] > ATAC_peak$start & H3K27ac[i,2] < ATAC_peak$end) | (H3K27ac[i,3] > ATAC_peak$start & H3K27ac[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K27ac) != 0){
    ATAC_peak_H3K27ac$V5 <- H3K27ac[i,1]
    ATAC_peak_H3K27ac$V6 <- H3K27ac[i,2]
    ATAC_peak_H3K27ac$V7 <- H3K27ac[i,3]
    ATAC_peak_H3K27ac$V8 <- H3K27ac[i,6]
    ATAC_peak_H3K27ac$V9 <- H3K27ac[i,7]
    ATAC_peak_H3K27ac_list = rbind(ATAC_peak_H3K27ac_list,ATAC_peak_H3K27ac)
  }
}
write.table(ATAC_peak_H3K27ac_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K27ac_list.txt",sep="\t")

##
H3K27ac_features=list("H3K27ac" = ATAC_peak_H3K27ac_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K27ac_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K27ac'
)

####H3K27me3
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K27me3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF795ZOS_H3K27me3.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K27me3)
H3K27me3@seqnames
write.table(H3K27me3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF795ZOS_H3K27me3.txt",sep = "\t")

H3K27me3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF795ZOS_H3K27me3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K27me3)){
  x=paste(H3K27me3[i,1],H3K27me3[i,2],H3K27me3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K27me3_features = list(list$V1)

ATAC_peak_H3K27me3_list = data.frame()
for(i in 1:nrow(H3K27me3)){
  ATAC_peak_H3K27me3 <- ATAC_peak[ATAC_peak$x1 == H3K27me3[i,1] & ((H3K27me3[i,2] > ATAC_peak$start & H3K27me3[i,2] < ATAC_peak$end) | (H3K27me3[i,3] > ATAC_peak$start & H3K27me3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K27me3) != 0){
    ATAC_peak_H3K27me3$V5 <- H3K27me3[i,1]
    ATAC_peak_H3K27me3$V6 <- H3K27me3[i,2]
    ATAC_peak_H3K27me3$V7 <- H3K27me3[i,3]
    ATAC_peak_H3K27me3$V8 <- H3K27me3[i,6]
    ATAC_peak_H3K27me3$V9 <- H3K27me3[i,7]
    ATAC_peak_H3K27me3_list = rbind(ATAC_peak_H3K27me3_list,ATAC_peak_H3K27me3)
  }
}
write.table(ATAC_peak_H3K27me3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K27me3_list.txt",sep="\t")

##
H3K27me3_features=list("H3K27me3" = ATAC_peak_H3K27me3_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K27me3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K27me3'
)

##H3K4me1
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K4me1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF135ZLM_H3K4me1.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K4me1)
H3K4me1@seqnames
write.table(H3K4me1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF135ZLM_H3K4me1.txt",sep = "\t")

H3K4me1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF135ZLM_H3K4me1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K4me1)){
  x=paste(H3K4me1[i,1],H3K4me1[i,2],H3K4me1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K4me1_features = list(list$V1)

ATAC_peak_H3K4me1_list = data.frame()
for(i in 1:nrow(H3K4me1)){
  ATAC_peak_H3K4me1 <- ATAC_peak[ATAC_peak$x1 == H3K4me1[i,1] & ((H3K4me1[i,2] > ATAC_peak$start & H3K4me1[i,2] < ATAC_peak$end) | (H3K4me1[i,3] > ATAC_peak$start & H3K4me1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K4me1) != 0){
    ATAC_peak_H3K4me1$V5 <- H3K4me1[i,1]
    ATAC_peak_H3K4me1$V6 <- H3K4me1[i,2]
    ATAC_peak_H3K4me1$V7 <- H3K4me1[i,3]
    ATAC_peak_H3K4me1$V8 <- H3K4me1[i,6]
    ATAC_peak_H3K4me1$V9 <- H3K4me1[i,7]
    ATAC_peak_H3K4me1_list = rbind(ATAC_peak_H3K4me1_list,ATAC_peak_H3K4me1)
  }
}
write.table(ATAC_peak_H3K4me1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K4me1_list.txt",sep="\t")

##
H3K4me1_features=list("H3K4me1" = ATAC_peak_H3K4me1_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K4me1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K4me1'
)

##H3K36me3
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K36me3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF520LHY_H3K36me3.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K36me3)
H3K36me3@seqnames
write.table(H3K36me3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF520LHY_H3K36me3.txt",sep = "\t")

H3K36me3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF520LHY_H3K36me3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K36me3)){
  x=paste(H3K36me3[i,1],H3K36me3[i,2],H3K36me3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K36me3_features = list(list$V1)

ATAC_peak_H3K36me3_list = data.frame()
for(i in 1:nrow(H3K36me3)){
  ATAC_peak_H3K36me3 <- ATAC_peak[ATAC_peak$x1 == H3K36me3[i,1] & ((H3K36me3[i,2] > ATAC_peak$start & H3K36me3[i,2] < ATAC_peak$end) | (H3K36me3[i,3] > ATAC_peak$start & H3K36me3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K36me3) != 0){
    ATAC_peak_H3K36me3$V5 <- H3K36me3[i,1]
    ATAC_peak_H3K36me3$V6 <- H3K36me3[i,2]
    ATAC_peak_H3K36me3$V7 <- H3K36me3[i,3]
    ATAC_peak_H3K36me3$V8 <- H3K36me3[i,6]
    ATAC_peak_H3K36me3$V9 <- H3K36me3[i,7]
    ATAC_peak_H3K36me3_list = rbind(ATAC_peak_H3K36me3_list,ATAC_peak_H3K36me3)
  }
}
write.table(ATAC_peak_H3K36me3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K36me3_list.txt",sep="\t")

##
H3K36me3_features=list("H3K36me3" = ATAC_peak_H3K36me3_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K36me3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K36me3'
)

###H3K9me3
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K9me3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF963GZJ_H3K9me3.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K9me3)
H3K9me3@seqnames
write.table(H3K9me3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF963GZJ_H3K9me3.txt",sep = "\t")

H3K9me3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF963GZJ_H3K9me3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K9me3)){
  x=paste(H3K9me3[i,1],H3K9me3[i,2],H3K9me3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K9me3_features = list(list$V1)

ATAC_peak_H3K9me3_list = data.frame()
for(i in 1:nrow(H3K9me3)){
  ATAC_peak_H3K9me3 <- ATAC_peak[ATAC_peak$x1 == H3K9me3[i,1] & ((H3K9me3[i,2] > ATAC_peak$start & H3K9me3[i,2] < ATAC_peak$end) | (H3K9me3[i,3] > ATAC_peak$start & H3K9me3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K9me3) != 0){
    ATAC_peak_H3K9me3$V5 <- H3K9me3[i,1]
    ATAC_peak_H3K9me3$V6 <- H3K9me3[i,2]
    ATAC_peak_H3K9me3$V7 <- H3K9me3[i,3]
    ATAC_peak_H3K9me3$V8 <- H3K9me3[i,6]
    ATAC_peak_H3K9me3$V9 <- H3K9me3[i,7]
    ATAC_peak_H3K9me3_list = rbind(ATAC_peak_H3K9me3_list,ATAC_peak_H3K9me3)
  }
}
write.table(ATAC_peak_H3K9me3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K9me3_list.txt",sep="\t")

##
H3K9me3_features=list("H3K9me3" = ATAC_peak_H3K9me3_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K9me3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K9me3'
)

##H3K9ac
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K9ac = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF558JOB_H3K9ac.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K9ac)
H3K9ac@seqnames
write.table(H3K9ac,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF558JOB_H3K9ac.txt",sep = "\t")

H3K9ac = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF558JOB_H3K9ac.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K9ac)){
  x=paste(H3K9ac[i,1],H3K9ac[i,2],H3K9ac[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K9ac_features = list(list$V1)

ATAC_peak_H3K9ac_list = data.frame()
for(i in 1:nrow(H3K9ac)){
  ATAC_peak_H3K9ac <- ATAC_peak[ATAC_peak$x1 == H3K9ac[i,1] & ((H3K9ac[i,2] > ATAC_peak$start & H3K9ac[i,2] < ATAC_peak$end) | (H3K9ac[i,3] > ATAC_peak$start & H3K9ac[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K9ac) != 0){
    ATAC_peak_H3K9ac$V5 <- H3K9ac[i,1]
    ATAC_peak_H3K9ac$V6 <- H3K9ac[i,2]
    ATAC_peak_H3K9ac$V7 <- H3K9ac[i,3]
    ATAC_peak_H3K9ac$V8 <- H3K9ac[i,6]
    ATAC_peak_H3K9ac$V9 <- H3K9ac[i,7]
    ATAC_peak_H3K9ac_list = rbind(ATAC_peak_H3K9ac_list,ATAC_peak_H3K9ac)
  }
}
write.table(ATAC_peak_H3K9ac_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K9ac_list.txt",sep="\t")

##
H3K9ac_features=list("H3K9ac" = ATAC_peak_H3K9ac_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K9ac_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K9ac'
)

##H3K4me2
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K4me2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF256AQN_H3K4me2.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K4me2)
H3K4me2@seqnames
write.table(H3K4me2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF256AQN_H3K4me2.txt",sep = "\t")

H3K4me2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF256AQN_H3K4me2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K4me2)){
  x=paste(H3K4me2[i,1],H3K4me2[i,2],H3K4me2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K4me2_features = list(list$V1)

ATAC_peak_H3K4me2_list = data.frame()
for(i in 1:nrow(H3K4me2)){
  ATAC_peak_H3K4me2<- ATAC_peak[ATAC_peak$x1 == H3K4me2[i,1] & ((H3K4me2[i,2] > ATAC_peak$start & H3K4me2[i,2] < ATAC_peak$end) | (H3K4me2[i,3] > ATAC_peak$start & H3K4me2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K4me2) != 0){
    ATAC_peak_H3K4me2$V5 <- H3K4me2[i,1]
    ATAC_peak_H3K4me2$V6 <- H3K4me2[i,2]
    ATAC_peak_H3K4me2$V7 <- H3K4me2[i,3]
    ATAC_peak_H3K4me2$V8 <- H3K4me2[i,6]
    ATAC_peak_H3K4me2$V9 <- H3K4me2[i,7]
    ATAC_peak_H3K4me2_list = rbind(ATAC_peak_H3K4me2_list,ATAC_peak_H3K4me2)
  }
}
write.table(ATAC_peak_H3K4me2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K4me2_list.txt",sep="\t")

##
H3K4me2_features=list("H3K4me2" = ATAC_peak_H3K4me2_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K4me2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K4me2'
)

##H4K20me1
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H4K20me1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF909RKY_H4K20me1.bed.gz",extraCols=extraCols_narrowPeak)

head(H4K20me1)
H4K20me1@seqnames
write.table(H4K20me1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF909RKY_H4K20me1.txt",sep = "\t")

H4K20me1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF909RKY_H4K20me1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H4K20me1)){
  x=paste(H4K20me1[i,1],H4K20me1[i,2],H4K20me1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H4K20me1_features = list(list$V1)

ATAC_peak_H4K20me1_list = data.frame()
for(i in 1:nrow(H4K20me1)){
  ATAC_peak_H4K20me1<- ATAC_peak[ATAC_peak$x1 == H4K20me1[i,1] & ((H4K20me1[i,2] > ATAC_peak$start & H4K20me1[i,2] < ATAC_peak$end) | (H4K20me1[i,3] > ATAC_peak$start & H4K20me1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H4K20me1) != 0){
    ATAC_peak_H4K20me1$V5 <- H4K20me1[i,1]
    ATAC_peak_H4K20me1$V6 <- H4K20me1[i,2]
    ATAC_peak_H4K20me1$V7 <- H4K20me1[i,3]
    ATAC_peak_H4K20me1$V8 <- H4K20me1[i,6]
    ATAC_peak_H4K20me1$V9 <- H4K20me1[i,7]
    ATAC_peak_H4K20me1_list = rbind(ATAC_peak_H4K20me1_list,ATAC_peak_H4K20me1)
  }
}
write.table(ATAC_peak_H4K20me1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H4K20me1_list.txt",sep="\t")

##
H4K20me1_features=list("H4K20me1" = ATAC_peak_H4K20me1_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H4K20me1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H4K20me1'
)

####H2AFZ
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H2AFZ= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF591KZT_H2AFZ.bed.gz",extraCols=extraCols_narrowPeak)

head(H2AFZ)
H2AFZ@seqnames
write.table(H2AFZ,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF591KZT_H2AFZ.txt",sep = "\t")

H2AFZ= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF591KZT_H2AFZ.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H2AFZ)){
  x=paste(H2AFZ[i,1],H2AFZ[i,2],H2AFZ[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H2AFZ_features = list(list$V1)

ATAC_peak_H2AFZ_list = data.frame()
for(i in 1:nrow(H2AFZ)){
  ATAC_peak_H2AFZ<- ATAC_peak[ATAC_peak$x1 == H2AFZ[i,1] & ((H2AFZ[i,2] > ATAC_peak$start & H2AFZ[i,2] < ATAC_peak$end) | (H2AFZ[i,3] > ATAC_peak$start & H2AFZ[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H2AFZ) != 0){
    ATAC_peak_H2AFZ$V5 <- H2AFZ[i,1]
    ATAC_peak_H2AFZ$V6 <- H2AFZ[i,2]
    ATAC_peak_H2AFZ$V7 <- H2AFZ[i,3]
    ATAC_peak_H2AFZ$V8 <- H2AFZ[i,6]
    ATAC_peak_H2AFZ$V9 <- H2AFZ[i,7]
    ATAC_peak_H2AFZ_list = rbind(ATAC_peak_H2AFZ_list,ATAC_peak_H2AFZ)
  }
}
write.table(ATAC_peak_H2AFZ_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H2AFZ_list.txt",sep="\t")

##
H2AFZ_features=list("H2AFZ" = ATAC_peak_H2AFZ_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H2AFZ_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H2AFZ'
)

##H3K79me2
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K79me2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF766GTJ_H3K79me2.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K79me2)
H3K79me2@seqnames
write.table(H3K79me2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF766GTJ_H3K79me2.txt",sep = "\t")

H3K79me2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF766GTJ_H3K79me2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K79me2)){
  x=paste(H3K79me2[i,1],H3K79me2[i,2],H3K79me2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K79me2_features = list(list$V1)

ATAC_peak_H3K79me2_list = data.frame()
for(i in 1:nrow(H3K79me2)){
  ATAC_peak_H3K79me2<- ATAC_peak[ATAC_peak$x1 == H3K79me2[i,1] & ((H3K79me2[i,2] > ATAC_peak$start & H3K79me2[i,2] < ATAC_peak$end) | (H3K79me2[i,3] > ATAC_peak$start & H3K79me2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K79me2) != 0){
    ATAC_peak_H3K79me2$V5 <- H3K79me2[i,1]
    ATAC_peak_H3K79me2$V6 <- H3K79me2[i,2]
    ATAC_peak_H3K79me2$V7 <- H3K79me2[i,3]
    ATAC_peak_H3K79me2$V8 <- H3K79me2[i,6]
    ATAC_peak_H3K79me2$V9 <- H3K79me2[i,7]
    ATAC_peak_H3K79me2_list = rbind(ATAC_peak_H3K79me2_list,ATAC_peak_H3K79me2)
  }
}
write.table(ATAC_peak_H3K79me2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K79me2_list.txt",sep="\t")

##
H3K79me2_features=list("H3K79me2" = ATAC_peak_H3K79me2_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K79me2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K79me2'
)

##H3K9me1
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

H3K9me1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF092WJN_H3K9me1.bed.gz",extraCols=extraCols_narrowPeak)

head(H3K9me1)
H3K9me1@seqnames
write.table(H3K9me1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF092WJN_H3K9me1.txt",sep = "\t")

H3K9me1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF092WJN_H3K9me1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K9me1)){
  x=paste(H3K9me1[i,1],H3K9me1[i,2],H3K9me1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K9me1_features = list(list$V1)

ATAC_peak_H3K9me1_list = data.frame()
for(i in 1:nrow(H3K9me1)){
  ATAC_peak_H3K9me1<- ATAC_peak[ATAC_peak$x1 == H3K9me1[i,1] & ((H3K9me1[i,2] > ATAC_peak$start & H3K9me1[i,2] < ATAC_peak$end) | (H3K9me1[i,3] > ATAC_peak$start & H3K9me1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K9me1) != 0){
    ATAC_peak_H3K9me1$V5 <- H3K9me1[i,1]
    ATAC_peak_H3K9me1$V6 <- H3K9me1[i,2]
    ATAC_peak_H3K9me1$V7 <- H3K9me1[i,3]
    ATAC_peak_H3K9me1$V8 <- H3K9me1[i,6]
    ATAC_peak_H3K9me1$V9 <- H3K9me1[i,7]
    ATAC_peak_H3K9me1_list = rbind(ATAC_peak_H3K9me1_list,ATAC_peak_H3K9me1)
  }
}
write.table(ATAC_peak_H3K9me1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K9me1_list.txt",sep="\t")

##
H3K9me1_features=list("H3K9me1" = ATAC_peak_H3K9me1_list$ATAC_peak)

DefaultAssay(scmulti4) <- 'peaks'
scmulti4 = AddChromatinModule(
  object = scmulti4,
  features = H3K9me1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K9me1'
)



##
all_control_0.03 <- subset(scmulti4, cells = c("AAACCGCGTACTGAAT-1","AAACGGATCCTGAATA-1","AAAGCCGCATTAGGCC-1","AAAGCTTGTTTGAGCA-1","AAAGGTTAGTCATCCC-1","AACATCATCCTGGTCT-1","AACTAGTGTAGCAGCT-1","AACTCACAGTAACTCA-1","AAGCCTCCATGAGTTT-1","AAGCGAATCTACCTCA-1","AAGCTAGAGTAGCGCC-1","AAGGATCCAAACAACA-1","AAGGTATAGGCCGGAA-1","AAGTGTTGTGAGGTGA-1",
                                                            "AAGTTACGTGAGCAAG-1","AATCATCCACCTCAGG-1","AATCTCAAGTTTCCGC-1","AATCTTGAGGACTAAG-1","AATTGCTCAGGCCTTG-1","ACAGGTAAGCACTAAC-1","ACATTAGTCAATCTAG-1","ACCAAACTCCATAAGC-1","ACCAAACTCTTAGGGT-1","ACCACATAGTTAGTGC-1","ACCAGGCTCGTTAAGC-1","ACCCTGTTCGGCCAGT-1","ACCGGTTCACATTAAC-1","ACGAACAAGGAACGGT-1",
                                                            "ACGATTCAGCGCCTTT-1","ACGTTACAGGCCTGGT-1","ACTAATCCAGCAAATA-1","ACTATGTCAGCGCTTG-1","ACTCACCTCGCTTGCT-1","ACTCAGTAGCCATCAG-1","ACTCGCGCAACAGCCT-1","ACTGAAACAGTAGCCT-1","ACTGAATGTCAGTAAT-1","ACTTACTTCTTGTCCA-1","AGAACCAAGGTCCGTA-1","AGACTATGTGAGCACT-1","AGAGAAGCACGTAATT-1","AGATGAAGTATTTGCC-1",
                                                            "AGCAACAAGGGTGAAC-1","AGCCGCTAGCGCAATT-1","AGCCGCTAGTAGCGCC-1","AGCTTAATCATGTTTC-1","AGGAACGGTCAAGTAT-1","AGGAGCTAGGTAAGCA-1","AGGCTAGCATGTTTGG-1","AGTAACGAGAAAGCAT-1","AGTACGCGTTAGTGAT-1","AGTAGGATCAACAAGG-1","AGTAGGATCCGCCAAA-1","AGTCAAGAGGGCCATC-1",
                                                            "AGTCAATGTATTCGCT-1","AGTTGGCGTGGGAACA-1","ATAGGTACATTAAGTC-1","ATCGCTTGTTCCATTA-1","ATGAATGCAGCCGCTA-1","ATGACGAAGCATGTCG-1","ATGACTCAGCGAGTAA-1","ATGCAAACACAGCCTG-1","ATGGCCGGTCGCAAAC-1","ATTAACCCAGAGAGCC-1","ATTACCCGTTATGTGG-1","ATTATGGTCATAACTG-1","ATTGCTCGTAACCAGC-1",
                                                            "CAACAATGTACTGATG-1","CAAGCTAGTGACATAT-1","CAAGTGAAGTTAGAGG-1","CACCTGTTCATGAGCT-1","CAGGCTATCCAAACAC-1","CAGGGCTTCTGGTCCT-1","CAGGTGGAGCAGGTTT-1","CATAATGTCGATTATG-1","CATAGACTCTAACCTT-1","CATAGGTTCGGTCATG-1","CATGAGGCAAGGAATC-1","CCAACATAGCGAGGTG-1","CCAACCAAGGCTGGCT-1","CCACACAAGCTTCTCA-1",
                                                            "CCACACAAGGCCTAAT-1","CCACACAAGTAAACCC-1","CCATAGCCAAACTGCC-1","CCATTATTCACGCCAA-1","CCCGCAACAAACGGGC-1","CCCGCAACAAGCCACT-1","CCCTCACCACCTGCTC-1","CCCTTAATCAGCACCA-1","CCTAATCGTTGTTGGA-1","CCTAGTTGTCCGCTGT-1","CCTCAAACACGAATCC-1","CCTGCTCCACCTCGCT-1","CCTGGTAAGCAACAAG-1",
                                                            "CCTTCGTAGTTTAACG-1","CGAACCGGTCACAGAC-1","CGAAGAATCAAACTCA-1","CGCAATCCAACTGGCT-1","CGCTCAGCAAGCGAGC-1","CGGAGCAAGACAGGTA-1","CGGCCATAGGAGCATA-1","CGGTGAACACTTCATC-1","CGTAACCCATAGCGAG-1","CGTCATTGTGGGAACA-1","CGTTAGGTCTCACAAA-1","CGTTTGGAGAGGAGTC-1","CGTTTGTGTATTGTGG-1",
                                                            "CTAATCCGTAACCTAG-1","CTAGCTTGTTCACCAT-1","CTAGGCGGTGGACATC-1","CTAGTCGAGGCAGGTG-1","CTAGTGAGTGTTTGTC-1","CTATGACAGGCGCTAC-1","CTATGGCCAATGCGCT-1","CTATGTTTCTTAATGG-1","CTATTGAAGTTGTCTT-1","CTCCCTGAGAGCAAGC-1","CTCCGGACAGGCTTCG-1","CTCTAGCTCAATCTAG-1","CTGCTCCCACCTATAG-1","CTGGACCAGCTGTCAG-1",
                                                            "CTGTAACAGCTATATG-1","CTTAACAAGGAGCAAC-1","CTTAGTTTCCCATAAA-1","CTTGGACCAAGTTATC-1","GAAGTAAGTTTGCTGT-1","GACATAGAGTAAAGGT-1","GACCTCAAGCTTAGCG-1","GACCTGCAGGTCCAAT-1","GACGCCTAGTTAGTTG-1","GAGAGGCGTAAAGCAA-1","GAGCAAGGTTAGCAGC-1","GAGTAATAGCGGATAA-1","GAGTCATTCCACCTTA-1",
                                                            "GATCGAGCAGTCTAAT-1","GATGGACAGGATGATG-1","GATTACTCACCTGTAA-1","GCAATAGAGAGGATAT-1","GCACCTAAGACTCGCA-1","GCACGGTTCACGCGGT-1","GCATGAAAGAAGCTAG-1","GCATGAGCACCTGCTC-1","GCATGAGCATTAAGCT-1","GCCATTACATTGTGAT-1","GCCTGACAGCCTCTCG-1","GCCTTAGAGGACAACA-1",
                                                            "GCGAAGCCACACTAAT-1","GCGGTTATCTGCAAAC-1","GCGTGCTAGTAGCCAT-1","GCGTTTCTCTGTAAGC-1","GCTATAGGTTAGCAGC-1","GCTCACAAGTTCCCGT-1","GCTCGATCATCCCGCT-1","GCTCTGGCATGTTTGG-1","GCTCTGTTCTGGCATG-1","GCTGACCAGCGATAAG-1","GCTGATCCAGGCCTTG-1","GGAACAATCACCATTT-1",
                                                            "GGAACAATCCCGTTTA-1","GGAACTAAGTATTGTG-1","GGATATTGTTGCACGG-1","GGATTATGTAATCGTG-1","GGCATTGTCCACCTGT-1","GGCGTTATCACCGGTA-1","GGGTGAAGTTACGCGG-1","GGTAGGAGTTCACCCA-1","GGTATGTTCACCAATA-1","GGTCCGTAGGACCGCT-1","GGTGCTTCACGCAACT-1","GGTGTTGTCCAAGTTA-1","GGTTGGTGTAGTTACG-1","GTAATAGCAGCTTACA-1",
                                                            "GTACAATGTACCGTTT-1","GTACTTCGTGGATTGC-1","GTATGTGGTAATCGTG-1","GTCATTAAGTAACCAC-1","GTGAATCTCCCGCCTA-1","GTGCACGGTAAAGCAA-1","GTTAAGTGTCACTCGC-1","GTTGGAGCAACATAAG-1","GTTTCTAGTCATAAGT-1","GTTTCTAGTTGCTTCG-1","TAACCTAAGGTGAAAT-1","TAATGGTGTTTGGCGG-1","TACAACATCTAAGTCA-1","TACATCAAGTAGGCGG-1",
                                                            "TACCAAATCCACCTTA-1","TACCGAAGTCCGGTTC-1","TACTGGCCACAATACT-1","TAGGAGTCAAGGTGGC-1","TAGGAGTCACCTAAGC-1","TAGGAGTCACCTAATG-1","TAGTTGTCAAATTCGT-1","TATAGCCAGCTTGCTC-1","TATCCAGCATAATTGC-1","TATGAAGCACCTAATG-1","TATGGATGTTAGCATG-1","TATGTGATCTTAGGGT-1","TCAAGAACAGCTACGT-1","TCACCGGCATTGTGAT-1",
                                                            "TCAGTGAGTATTGAGT-1","TCAGTGAGTCTAACCT-1","TCATTGTTCCTAATAG-1","TCCAGGTCAGCATGGA-1","TCCCGGACACTAAGCC-1","TCGTTAGCAGGATGGC-1","TGAAGCAAGACAAACG-1","TGACTTCGTATTGGTG-1","TGAGTTTCAGATAGAC-1","TGATGAACATGGCCCA-1","TGCCATTGTCACGGAT-1",
                                                            "TGCTCAACACCTCAGG-1","TGGCGGTTCCTTAAGA-1","TGGTCAGTCTGTGAGT-1","TGGTTCCTCAGGATGA-1","TGGTTCTGTCATCAAG-1","TGGTTCTGTGATGATT-1","TGTGAAACAAACCTTG-1","TGTGCAAGTAACGAGG-1","TGTGTTAAGTGTTGTA-1","TGTTACTTCGCAGGCT-1",
                                                            "TGTTGCACACTAAGCC-1","TTAAAGGCAGGCCTTG-1","TTACCTGTCTTTAGGA-1","TTACGTTTCTCACAAA-1","TTAGCCTGTGGATGTC-1","TTAGCGGTCGGTTAGT-1","TTAGGCGTCACCGGTA-1","TTATTGCTCATTATGG-1","TTATTGCTCTGCAAGT-1","TTCAGGTAGGATCCGC-1","TTCGCAACAGAAACGT-1","TTGCCCGTCATTTAGG-1","TTGGCTTGTATTGTCC-1",
                                                            "TTGGTGAGTCTAACCT-1","TTGTGAGGTAGCTAAT-1","TTGTTCCCAAGGTAAC-1","TTGTTTGGTGCATTTC-1","TTTATGGAGGCTACTG-1","TTTGTGTTCTCAATAG-1"
                                                            ,"AGTATAGCAAGGCCAA-1","CCAGACTCACACTAAT-1","CGAAATGAGAGGCTAA-1","GAGTGAGGTGGAAGGC-1","GTTCCTGGTGCATCGG-1","TTGTCCGGTTTGCTGT-1"
                                                            
))

all_control_0.03_H3K4me3 = as.data.frame(all_control_0.03@meta.data$H3K4me3)


#NEAT1
cell_NEAT1_KD <- subset(scmulti4, cells = c("AAGGAAGCAAACGGGC-1","ATTACCGCACGTAATT-1","CAGGAAGGTAATCCCT-1","CATTTGTTCGGGCCAT-1","CGATATTCAATAGCAA-1","CGCCTGTGTAACCAGC-1","GAACCGCTCAGCTAAC-1","GCACCTAAGTTACTTC-1","GCATTGCCAGCATGGA-1","GCGCAATGTCGCATAA-1","GTGCTCAAGTTCCCAC-1","GTTACGTAGCAAACCT-1","TAGCTTAAGGCGCTAC-1","TTACCGTGTGTTGCAA-1",
                                                         "TTAGACTCATAGACTT-1","TTAGCAATCCAGGAAA-1","TTGCATTTCATTGTCT-1"
))

cell_NEAT1_KD_H3K4me3 = as.data.frame(cell_NEAT1_KD@meta.data$H3K4me3)

ja = as.data.frame(pfm@listData$MA0025.1)

#
NEAT1_NT0.03_H3K4me3 = wilcox.test(all_control_0.03_H3K4me3[,1],cell_NEAT1_KD_H3K4me3[,1],alternative = "two.sided")
NEAT1_NT0.03_H3K4me3_p = NEAT1_NT0.03_H3K4me3$p.value
NEAT1_NT0.03_H3K4me3_diff = mean(cell_NEAT1_KD_H3K4me3[,1])-mean(all_control_0.03_H3K4me3[,1])


#LINC00578
cell_LINC00578_KD <- subset(scmulti4, cells = c("AACGACAAGTTAGCTA-1","AATGGCGCAGGAACCA-1","ACACCTTGTTGTCATC-1","ACTAGGCGTTAAGCTG-1","ACTCACCTCATCCTCA-1","AGGATGTCATGGTTAT-1","AGGCTAGCAAATTCGT-1","AGTAACGAGACAGGCG-1","ATCGCCCGTTCACCAT-1","CAAGTGAAGGGTGAAC-1",
                                                             "CAATGACTCACCAATA-1","CACAATATCCAGGTCA-1","CACGCTAAGTGCTGTG-1","CAGCCTAAGGGTCCAC-1","CAGGGCTTCATTATGG-1","CCAACCAAGGCTAAGA-1",
                                                             "CCCAAATAGGTTAGAG-1","CCCAATTGTGGACCTG-1","CCCTCACCAGTTTCTC-1","CCGCAAGGTTAGTTGG-1","CCTGATGAGAATCTCA-1","CGACAAGCAATGAGGT-1","CGAGGCAAGTCAGTAC-1","CGCAAATTCTCCTCAA-1",
                                                             "CGGGTCTAGGAGTCGG-1","CGTAACCCAGCCTGCA-1","CGTAACTAGCTACTGG-1","CGTATTGCAGTTATCG-1","CTAGGACGTGACCTGG-1","CTAGTCGAGTCTCACC-1",
                                                             "CTCCCTGAGTAGCGGG-1","CTCGCTAGTAAGCTCA-1","CTCGCTAGTAGCCATA-1","CTGCTCCCACAGAAAC-1","CTGGCTAAGCTGTACG-1","CTTCAAGCACCAACCG-1","CTTTAGTTCAAACTCA-1","CTTTATCAGTTGCCTC-1",
                                                             "GAAAGGCTCTCGCCTG-1","GACATAGAGGCTGGCT-1","GATTGATGTAGCTGGT-1","GCATTGCCAAATGCCC-1","GCCAATTAGCGTGCAC-1","GCCCAAATCTCCTCTT-1",
                                                             "GCCTACTTCTTAGCGG-1","GCCTGCTGTATCTGGA-1","GCGCGATTCATGCCTC-1","GCTGGTTCAGCTTAAT-1","GGAAGCTAGTCACGAT-1","GGACGGATCAATCTCT-1","GGATTATGTCCTTTAA-1","GGGATTAAGTTAGTGC-1",
                                                             "GGTTAATGTCCTAACT-1","GTACTTCGTATTGTGG-1","GTCCTAGAGATGGACA-1","GTGGCTTCAAATGCCC-1","TACTAAGTCCCGCATT-1","TAGCTAGGTGCAATAT-1","TATATCCTCGTTATCT-1","TCGACAAGTTGCAGTA-1",
                                                             "TCTTCAAGTATCTGGA-1","TGCTCACTCAGCATTA-1","TGCTTAAAGTAGCTTA-1","TGTCCTGGTTTCCACG-1","TTAGCGGTCCTGAATA-1","TTGACGTAGTGCACGC-1","TTGGTGAGTTAGCCAA-1",
                                                             "TTGTGAGGTGTTGCAA-1","TTTCTTGCAGGAATCG-1","TTTGTCCCATTGTGCA-1","TTTGTGGCATAATGAG-1"))
cell_LINC00578_KD_H3K4me3 = as.data.frame(cell_LINC00578_KD@meta.data$H3K4me3)

LINC00578_NT0.03_H3K4me3 = wilcox.test(all_control_0.03_H3K4me3[,1],cell_LINC00578_KD_H3K4me3[,1],alternative = "two.sided")
LINC00578_NT0.03_H3K4me3_p = LINC00578_NT0.03_H3K4me3$p.value
LINC00578_NT0.03_H3K4me3_diff = mean(cell_LINC00578_KD_H3K4me3[,1])-mean(all_control_0.03_H3K4me3[,1])



library(rtracklayer)
# declare names and types for the extra BED columns
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

###Bld.50.CDKN1B.AllCell
Bld.50.CDKN1B.AllCell = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.50.CDKN1B.AllCell.bed")

head(Bld.50.CDKN1B.AllCell)
Bld.50.CDKN1B.AllCell@seqnames
write.table(Bld.50.CDKN1B.AllCell,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF246IEW_Bld.50.CDKN1B.AllCell.txt",sep = "\t")

Bld.50.CDKN1B.AllCell = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF246IEW_Bld.50.CDKN1B.AllCell.txt",sep = "\t")


AllAg.Lung_fibroblasts = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Lng.50.AllAg.Lung_fibroblasts.bed")

head(AllAg.Lung_fibroblasts)
AllAg.Lung_fibroblasts@seqnames
write.table(AllAg.Lung_fibroblasts,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/AllAg.Lung_fibroblasts.txt",sep = "\t")

AllAg.Lung_fibroblasts = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/AllAg.Lung_fibroblasts.txt",sep = "\t")



pfmlist = as.matrix(pfm@listData)
#pfmlist = as.data.frame(pfmlist)
#write.table(pfmlist,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/pfmlist.txt",sep = "\t")

write.table(pfmlist, file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/pfmlist.txt", sep = ",", row.names = FALSE, col.names = TRUE)

pfmlist2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/pfmlist.txt",sep = ",",fill = T)

pfmlist2 = as.data.frame(pfmlist2[pfmlist2$V1 == 'new(PFMatrix',])


###JASPARlist

JASPARlist = as.data.frame(cbind(pfmlist2$V2,pfmlist2$V3))
JASPARlist$V1 <- sub(" ID = ", "", JASPARlist$V1)
JASPARlist$V2 <- sub(" name = ", "", JASPARlist$V2)

write.table(JASPARlist,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/JASPARlist.txt",sep = "\t")
#JASPARlist2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/JASPARlist2.txt",sep = "\t")


##ENCODE
ENCODE = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/experiment_report_2023_5_28_7h_52m.txt",sep = "\t",header = T,quote = "")
ENCODE_ID = as.data.frame(ENCODE[,3])

##senegene
senegene = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/cellage3.tsv",sep = "\t",header = T,quote = "")
senegene_ID = as.data.frame(senegene[,2])

senesig = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/signatures1.txt",sep = "\t",header = T,quote = "")
senesig_ID = as.data.frame(senesig[,1])

##
ENCODE_senegene = intersect(as.vector(ENCODE_ID$`ENCODE[, 3]`),as.vector(senegene_ID$`senegene[, 2]`))
write.table(ENCODE_senegene,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCODE_senegene.txt",sep = '\t')

ENCODE_senesig = intersect(as.vector(ENCODE_ID$`ENCODE[, 3]`),as.vector(senesig_ID$`senesig[, 1]`))
write.table(ENCODE_senesig,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCODE_senesig.txt",sep = '\t')

JASPAR_senegene = intersect(as.vector(JASPARlist2$JASPARlist2),as.vector(senegene_ID$`senegene[, 2]`))
write.table(JASPAR_senegene ,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/JASPAR_senegene.txt",sep = '\t')

JASPAR_senesig = intersect(as.vector(JASPARlist2$JASPARlist2),as.vector(senesig_ID$`senesig[, 1]`))
write.table(JASPAR_senesig,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/JASPAR_senesig.txt",sep = '\t')

#####

pfm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)

# Scan the DNA sequence of each peak for the presence of each motif
# 使用CreateMotifMatrix函数构建Motif矩阵对象
#scmulti4 = scmulti3
DefaultAssay(object =  scmulti2) <- 'peaks'
# add motif information
scmulti2 <- AddMotifs(
  object = scmulti2,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)


########
#Computing motif activities
scmulti2  <- RunChromVAR(
  object = scmulti2,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(scmulti2) <- 'chromvar'


library(rtracklayer)
# declare names and types for the extra BED columns
extraCols_narrowPeak <- c(FoldChange="numeric", pVal="numeric",
                          qVal="numeric", summit="integer")

###H3K4me3
H3K4me3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF246IEW_H3K4me3.bed",extraCols=extraCols_narrowPeak)

head(H3K4me3)
H3K4me3@seqnames
write.table(H3K4me3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF246IEW_H3K4me3.txt",sep = "\t")

H3K4me3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ENCFF246IEW_H3K4me3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(H3K4me3)){
  x=paste(H3K4me3[i,1],H3K4me3[i,2],H3K4me3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

H3K4me3_features = list(list$V1)

#

ATAC_peak = as.data.frame(scmulti2@assays$peaks@data)
ATAC_peak = rownames(ATAC_peak)
ATAC_peak = as.data.frame(ATAC_peak)

x1 = vector()
y1 = vector()
z1 = vector()
for (i in 1:nrow(ATAC_peak)){
  x = strsplit(ATAC_peak[i,1],split='-')[[1]][1]
  y = strsplit(ATAC_peak[i,1],split='-')[[1]][2]
  z = strsplit(ATAC_peak[i,1],split='-')[[1]][3]
  x1 = c(x1,x)
  y1 = c(y1,y)
  z1 = c(z1,z)
}
ATAC_peak = cbind(ATAC_peak,x1,as.numeric(y1),as.numeric(z1))
colnames(ATAC_peak)[3]="start"
colnames(ATAC_peak)[4]="end"
write.table(ATAC_peak,file="/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_k562.txt",sep="\t")

ATAC_peak_H3K4me3_list = data.frame()
for(i in 1:nrow(H3K4me3)){
  ATAC_peak_H3K4me3 <- ATAC_peak[ATAC_peak$x1 == H3K4me3[i,1] & ((H3K4me3[i,2] > ATAC_peak$start & H3K4me3[i,2] < ATAC_peak$end) | (H3K4me3[i,3] > ATAC_peak$start & H3K4me3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_H3K4me3) != 0){
    ATAC_peak_H3K4me3$V5 <- H3K4me3[i,1]
    ATAC_peak_H3K4me3$V6 <- H3K4me3[i,2]
    ATAC_peak_H3K4me3$V7 <- H3K4me3[i,3]
    ATAC_peak_H3K4me3$V8 <- H3K4me3[i,6]
    ATAC_peak_H3K4me3$V9 <- H3K4me3[i,7]
    ATAC_peak_H3K4me3_list = rbind(ATAC_peak_H3K4me3_list,ATAC_peak_H3K4me3)
  }
}
write.table(ATAC_peak_H3K4me3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K4me3_list.txt",sep="\t")

##
H3K4me3_features=list("H3K4me3" = ATAC_peak_H3K4me3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K4me3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K4me3'
)

#kk = as.data.frame(scmulti4$H3K4me3)
write.table(ATAC_peak_H3K4me3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K4me3_list.txt",sep="\t")

###
##H3K27ac
ATAC_peak_H3K27ac_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K27ac_list.txt",sep="\t")
##
H3K27ac_features=list("H3K27ac" = ATAC_peak_H3K27ac_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K27ac_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K27ac'
)

####H3K27me3

ATAC_peak_H3K27me3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K27me3_list.txt",sep="\t")
##
H3K27me3_features=list("H3K27me3" = ATAC_peak_H3K27me3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K27me3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K27me3'
)

##H3K4me1
ATAC_peak_H3K4me1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K4me1_list.txt",sep="\t")

##
H3K4me1_features=list("H3K4me1" = ATAC_peak_H3K4me1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K4me1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K4me1'
)

##H3K36me3
ATAC_peak_H3K36me3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K36me3_list.txt",sep="\t")

##
H3K36me3_features=list("H3K36me3" = ATAC_peak_H3K36me3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K36me3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K36me3'
)

###H3K9me3
ATAC_peak_H3K9me3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K9me3_list.txt",sep="\t")

##
H3K9me3_features=list("H3K9me3" = ATAC_peak_H3K9me3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K9me3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K9me3'
)

##H3K9ac

ATAC_peak_H3K9ac_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K9ac_list.txt",sep="\t")

##
H3K9ac_features=list("H3K9ac" = ATAC_peak_H3K9ac_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K9ac_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K9ac'
)

##H3K4me2
ATAC_peak_H3K4me2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K4me2_list.txt",sep="\t")


##
H3K4me2_features=list("H3K4me2" = ATAC_peak_H3K4me2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K4me2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K4me2'
)

##H4K20me1

ATAC_peak_H4K20me1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H4K20me1_list.txt",sep="\t")

##
H4K20me1_features=list("H4K20me1" = ATAC_peak_H4K20me1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H4K20me1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H4K20me1'
)

####H2AFZ
ATAC_peak_H2AFZ_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H2AFZ_list.txt",sep="\t")

##
H2AFZ_features=list("H2AFZ" = ATAC_peak_H2AFZ_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H2AFZ_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H2AFZ'
)

##H3K79me2
ATAC_peak_H3K79me2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K79me2_list.txt",sep="\t")

##
H3K79me2_features=list("H3K79me2" = ATAC_peak_H3K79me2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K79me2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K79me2'
)

##H3K9me1
ATAC_peak_H3K9me1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/ATAC_peak_H3K9me1_list.txt",sep="\t")

##
H3K9me1_features=list("H3K9me1" = ATAC_peak_H3K9me1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = H3K9me1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'H3K9me1'
)
scmulti2 = scmulti4

##ZNF184

ZNF184= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF008VFR_ZNF184.bed.gz",extraCols=extraCols_narrowPeak)

head(ZNF184)
ZNF184@seqnames
write.table(ZNF184,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF256AQN_ZNF184.txt",sep = "\t")

ZNF184= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF256AQN_ZNF184.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF184)){
  x=paste(ZNF184[i,1],ZNF184[i,2],ZNF184[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF184_features = list(list$V1)

ATAC_peak_ZNF184_list = data.frame()
for(i in 1:nrow(ZNF184)){
  ATAC_peak_ZNF184<- ATAC_peak[ATAC_peak$x1 == ZNF184[i,1] & ((ZNF184[i,2] > ATAC_peak$start & ZNF184[i,2] < ATAC_peak$end) | (ZNF184[i,3] > ATAC_peak$start & ZNF184[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF184) != 0){
    ATAC_peak_ZNF184$V5 <- ZNF184[i,1]
    ATAC_peak_ZNF184$V6 <- ZNF184[i,2]
    ATAC_peak_ZNF184$V7 <- ZNF184[i,3]
    ATAC_peak_ZNF184$V8 <- ZNF184[i,6]
    ATAC_peak_ZNF184$V9 <- ZNF184[i,7]
    ATAC_peak_ZNF184_list = rbind(ATAC_peak_ZNF184_list,ATAC_peak_ZNF184)
  }
}
write.table(ATAC_peak_ZNF184_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF184_list.txt",sep="\t")

ATAC_peak_ZNF184_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF184_list.txt",sep="\t")

ZNF184_features=list("ZNF184" = ATAC_peak_ZNF184_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF184_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  #name = 'ZNF184'
)


##TAF7
TAF7= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF010YHS_TAF7.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TAF7,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF010YHS_TAF7.txt",sep = "\t")

TAF7= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF010YHS_TAF7.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TAF7)){
  x=paste(TAF7[i,1],TAF7[i,2],TAF7[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TAF7_features = list(list$V1)

ATAC_peak_TAF7_list = data.frame()
for(i in 1:nrow(TAF7)){
  ATAC_peak_TAF7<- ATAC_peak[ATAC_peak$x1 == TAF7[i,1] & ((TAF7[i,2] > ATAC_peak$start & TAF7[i,2] < ATAC_peak$end) | (TAF7[i,3] > ATAC_peak$start & TAF7[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TAF7) != 0){
    ATAC_peak_TAF7$V5 <- TAF7[i,1]
    ATAC_peak_TAF7$V6 <- TAF7[i,2]
    ATAC_peak_TAF7$V7 <- TAF7[i,3]
    ATAC_peak_TAF7$V8 <- TAF7[i,6]
    ATAC_peak_TAF7$V9 <- TAF7[i,7]
    ATAC_peak_TAF7_list = rbind(ATAC_peak_TAF7_list,ATAC_peak_TAF7)
  }
}
write.table(ATAC_peak_TAF7_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAF7_list.txt",sep="\t")

ATAC_peak_TAF7_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAF7_list.txt",sep="\t")

TAF7_features=list("TAF7" = ATAC_peak_TAF7_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TAF7_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##PTRF
PTRF= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF011CUC_PTRF.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PTRF,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF011CUC_PTRF.txt",sep = "\t")

PTRF= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF011CUC_PTRF.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PTRF)){
  x=paste(PTRF[i,1],PTRF[i,2],PTRF[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PTRF_features = list(list$V1)

ATAC_peak_PTRF_list = data.frame()
for(i in 1:nrow(PTRF)){
  ATAC_peak_PTRF<- ATAC_peak[ATAC_peak$x1 == PTRF[i,1] & ((PTRF[i,2] > ATAC_peak$start & PTRF[i,2] < ATAC_peak$end) | (PTRF[i,3] > ATAC_peak$start & PTRF[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PTRF) != 0){
    ATAC_peak_PTRF$V5 <- PTRF[i,1]
    ATAC_peak_PTRF$V6 <- PTRF[i,2]
    ATAC_peak_PTRF$V7 <- PTRF[i,3]
    ATAC_peak_PTRF$V8 <- PTRF[i,6]
    ATAC_peak_PTRF$V9 <- PTRF[i,7]
    ATAC_peak_PTRF_list = rbind(ATAC_peak_PTRF_list,ATAC_peak_PTRF)
  }
}
write.table(ATAC_peak_PTRF_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PTRF_list.txt",sep="\t")

ATAC_peak_PTRF_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PTRF_list.txt",sep="\t")

PTRF_features=list("PTRF" = ATAC_peak_PTRF_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PTRF_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CBFB
CBFB= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF012EVG_CBFB.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CBFB,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF012EVG_CBFB.txt",sep = "\t")

CBFB= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF012EVG_CBFB.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CBFB)){
  x=paste(CBFB[i,1],CBFB[i,2],CBFB[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CBFB_features = list(list$V1)

ATAC_peak_CBFB_list = data.frame()
for(i in 1:nrow(CBFB)){
  ATAC_peak_CBFB<- ATAC_peak[ATAC_peak$x1 == CBFB[i,1] & ((CBFB[i,2] > ATAC_peak$start & CBFB[i,2] < ATAC_peak$end) | (CBFB[i,3] > ATAC_peak$start & CBFB[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CBFB) != 0){
    ATAC_peak_CBFB$V5 <- CBFB[i,1]
    ATAC_peak_CBFB$V6 <- CBFB[i,2]
    ATAC_peak_CBFB$V7 <- CBFB[i,3]
    ATAC_peak_CBFB$V8 <- CBFB[i,6]
    ATAC_peak_CBFB$V9 <- CBFB[i,7]
    ATAC_peak_CBFB_list = rbind(ATAC_peak_CBFB_list,ATAC_peak_CBFB)
  }
}
write.table(ATAC_peak_CBFB_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBFB_list.txt",sep="\t")

ATAC_peak_CBFB_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBFB_list.txt",sep="\t")

CBFB_features=list("CBFB" = ATAC_peak_CBFB_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CBFB_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MAZ
MAZ= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF014UGE_MAZ.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MAZ,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF014UGE_MAZ.txt",sep = "\t")

MAZ= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF014UGE_MAZ.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MAZ)){
  x=paste(MAZ[i,1],MAZ[i,2],MAZ[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MAZ_features = list(list$V1)

ATAC_peak_MAZ_list = data.frame()
for(i in 1:nrow(MAZ)){
  ATAC_peak_MAZ<- ATAC_peak[ATAC_peak$x1 == MAZ[i,1] & ((MAZ[i,2] > ATAC_peak$start & MAZ[i,2] < ATAC_peak$end) | (MAZ[i,3] > ATAC_peak$start & MAZ[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MAZ) != 0){
    ATAC_peak_MAZ$V5 <- MAZ[i,1]
    ATAC_peak_MAZ$V6 <- MAZ[i,2]
    ATAC_peak_MAZ$V7 <- MAZ[i,3]
    ATAC_peak_MAZ$V8 <- MAZ[i,6]
    ATAC_peak_MAZ$V9 <- MAZ[i,7]
    ATAC_peak_MAZ_list = rbind(ATAC_peak_MAZ_list,ATAC_peak_MAZ)
  }
}
write.table(ATAC_peak_MAZ_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MAZ_list.txt",sep="\t")

ATAC_peak_MAZ_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MAZ_list.txt",sep="\t")

MAZ_features=list("MAZ" = ATAC_peak_MAZ_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MAZ_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RUNX1
RUNX1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF015ADG_RUNX1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RUNX1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF015ADG_RUNX1.txt",sep = "\t")

RUNX1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF015ADG_RUNX1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RUNX1)){
  x=paste(RUNX1[i,1],RUNX1[i,2],RUNX1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RUNX1_features = list(list$V1)

ATAC_peak_RUNX1_list = data.frame()
for(i in 1:nrow(RUNX1)){
  ATAC_peak_RUNX1<- ATAC_peak[ATAC_peak$x1 == RUNX1[i,1] & ((RUNX1[i,2] > ATAC_peak$start & RUNX1[i,2] < ATAC_peak$end) | (RUNX1[i,3] > ATAC_peak$start & RUNX1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RUNX1) != 0){
    ATAC_peak_RUNX1$V5 <- RUNX1[i,1]
    ATAC_peak_RUNX1$V6 <- RUNX1[i,2]
    ATAC_peak_RUNX1$V7 <- RUNX1[i,3]
    ATAC_peak_RUNX1$V8 <- RUNX1[i,6]
    ATAC_peak_RUNX1$V9 <- RUNX1[i,7]
    ATAC_peak_RUNX1_list = rbind(ATAC_peak_RUNX1_list,ATAC_peak_RUNX1)
  }
}
write.table(ATAC_peak_RUNX1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RUNX1_list.txt",sep="\t")

ATAC_peak_RUNX1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RUNX1_list.txt",sep="\t")

RUNX1_features=list("RUNX1" = ATAC_peak_RUNX1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RUNX1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ESRRA
ESRRA= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF015PUY_ESRRA.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ESRRA,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF015PUY_ESRRA.txt",sep = "\t")

ESRRA= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF015PUY_ESRRA.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ESRRA)){
  x=paste(ESRRA[i,1],ESRRA[i,2],ESRRA[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ESRRA_features = list(list$V1)

ATAC_peak_ESRRA_list = data.frame()
for(i in 1:nrow(ESRRA)){
  ATAC_peak_ESRRA<- ATAC_peak[ATAC_peak$x1 == ESRRA[i,1] & ((ESRRA[i,2] > ATAC_peak$start & ESRRA[i,2] < ATAC_peak$end) | (ESRRA[i,3] > ATAC_peak$start & ESRRA[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ESRRA) != 0){
    ATAC_peak_ESRRA$V5 <- ESRRA[i,1]
    ATAC_peak_ESRRA$V6 <- ESRRA[i,2]
    ATAC_peak_ESRRA$V7 <- ESRRA[i,3]
    ATAC_peak_ESRRA$V8 <- ESRRA[i,6]
    ATAC_peak_ESRRA$V9 <- ESRRA[i,7]
    ATAC_peak_ESRRA_list = rbind(ATAC_peak_ESRRA_list,ATAC_peak_ESRRA)
  }
}
write.table(ATAC_peak_ESRRA_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ESRRA_list.txt",sep="\t")

ATAC_peak_ESRRA_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ESRRA_list.txt",sep="\t")

ESRRA_features=list("ESRRA" = ATAC_peak_ESRRA_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ESRRA_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##BMI1
BMI1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF023GES_BMI1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BMI1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF023GES_BMI1_age.txt",sep = "\t")

BMI1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF023GES_BMI1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BMI1)){
  x=paste(BMI1[i,1],BMI1[i,2],BMI1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BMI1_features = list(list$V1)

ATAC_peak_BMI1_list = data.frame()
for(i in 1:nrow(BMI1)){
  ATAC_peak_BMI1<- ATAC_peak[ATAC_peak$x1 == BMI1[i,1] & ((BMI1[i,2] > ATAC_peak$start & BMI1[i,2] < ATAC_peak$end) | (BMI1[i,3] > ATAC_peak$start & BMI1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BMI1) != 0){
    ATAC_peak_BMI1$V5 <- BMI1[i,1]
    ATAC_peak_BMI1$V6 <- BMI1[i,2]
    ATAC_peak_BMI1$V7 <- BMI1[i,3]
    ATAC_peak_BMI1$V8 <- BMI1[i,6]
    ATAC_peak_BMI1$V9 <- BMI1[i,7]
    ATAC_peak_BMI1_list = rbind(ATAC_peak_BMI1_list,ATAC_peak_BMI1)
  }
}
write.table(ATAC_peak_BMI1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BMI1_list.txt",sep="\t")

ATAC_peak_BMI1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BMI1_list.txt",sep="\t")

BMI1_features=list("BMI1" = ATAC_peak_BMI1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BMI1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NELFE
NELFE= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF024VAA_NELFE.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NELFE,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF024VAA_NELFE.txt",sep = "\t")

NELFE= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF024VAA_NELFE.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NELFE)){
  x=paste(NELFE[i,1],NELFE[i,2],NELFE[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NELFE_features = list(list$V1)

ATAC_peak_NELFE_list = data.frame()
for(i in 1:nrow(NELFE)){
  ATAC_peak_NELFE<- ATAC_peak[ATAC_peak$x1 == NELFE[i,1] & ((NELFE[i,2] > ATAC_peak$start & NELFE[i,2] < ATAC_peak$end) | (NELFE[i,3] > ATAC_peak$start & NELFE[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NELFE) != 0){
    ATAC_peak_NELFE$V5 <- NELFE[i,1]
    ATAC_peak_NELFE$V6 <- NELFE[i,2]
    ATAC_peak_NELFE$V7 <- NELFE[i,3]
    ATAC_peak_NELFE$V8 <- NELFE[i,6]
    ATAC_peak_NELFE$V9 <- NELFE[i,7]
    ATAC_peak_NELFE_list = rbind(ATAC_peak_NELFE_list,ATAC_peak_NELFE)
  }
}
write.table(ATAC_peak_NELFE_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NELFE_list.txt",sep="\t")

ATAC_peak_NELFE_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NELFE_list.txt",sep="\t")

NELFE_features=list("NELFE" = ATAC_peak_NELFE_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NELFE_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##MAX
MAX= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF034DOZ_MAX_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MAX,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF034DOZ_MAX_age.txt",sep = "\t")

MAX= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF034DOZ_MAX_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MAX)){
  x=paste(MAX[i,1],MAX[i,2],MAX[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MAX_features = list(list$V1)

ATAC_peak_MAX_list = data.frame()
for(i in 1:nrow(MAX)){
  ATAC_peak_MAX<- ATAC_peak[ATAC_peak$x1 == MAX[i,1] & ((MAX[i,2] > ATAC_peak$start & MAX[i,2] < ATAC_peak$end) | (MAX[i,3] > ATAC_peak$start & MAX[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MAX) != 0){
    ATAC_peak_MAX$V5 <- MAX[i,1]
    ATAC_peak_MAX$V6 <- MAX[i,2]
    ATAC_peak_MAX$V7 <- MAX[i,3]
    ATAC_peak_MAX$V8 <- MAX[i,6]
    ATAC_peak_MAX$V9 <- MAX[i,7]
    ATAC_peak_MAX_list = rbind(ATAC_peak_MAX_list,ATAC_peak_MAX)
  }
}
write.table(ATAC_peak_MAX_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MAX_list.txt",sep="\t")

ATAC_peak_MAX_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MAX_list.txt",sep="\t")

MAX_features=list("MAX" = ATAC_peak_MAX_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MAX_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZC3H8
ZC3H8= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF039UIB_ZC3H8.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZC3H8,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF039UIB_ZC3H8.txt",sep = "\t")

ZC3H8= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF039UIB_ZC3H8.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZC3H8)){
  x=paste(ZC3H8[i,1],ZC3H8[i,2],ZC3H8[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZC3H8_features = list(list$V1)

ATAC_peak_ZC3H8_list = data.frame()
for(i in 1:nrow(ZC3H8)){
  ATAC_peak_ZC3H8<- ATAC_peak[ATAC_peak$x1 == ZC3H8[i,1] & ((ZC3H8[i,2] > ATAC_peak$start & ZC3H8[i,2] < ATAC_peak$end) | (ZC3H8[i,3] > ATAC_peak$start & ZC3H8[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZC3H8) != 0){
    ATAC_peak_ZC3H8$V5 <- ZC3H8[i,1]
    ATAC_peak_ZC3H8$V6 <- ZC3H8[i,2]
    ATAC_peak_ZC3H8$V7 <- ZC3H8[i,3]
    ATAC_peak_ZC3H8$V8 <- ZC3H8[i,6]
    ATAC_peak_ZC3H8$V9 <- ZC3H8[i,7]
    ATAC_peak_ZC3H8_list = rbind(ATAC_peak_ZC3H8_list,ATAC_peak_ZC3H8)
  }
}
write.table(ATAC_peak_ZC3H8_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZC3H8_list.txt",sep="\t")

ATAC_peak_ZC3H8_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZC3H8_list.txt",sep="\t")

ZC3H8_features=list("ZC3H8" = ATAC_peak_ZC3H8_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZC3H8_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SIN3B
SIN3B= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF040MIA_SIN3B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SIN3B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF040MIA_SIN3B.txt",sep = "\t")

SIN3B= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF040MIA_SIN3B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SIN3B)){
  x=paste(SIN3B[i,1],SIN3B[i,2],SIN3B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SIN3B_features = list(list$V1)

ATAC_peak_SIN3B_list = data.frame()
for(i in 1:nrow(SIN3B)){
  ATAC_peak_SIN3B<- ATAC_peak[ATAC_peak$x1 == SIN3B[i,1] & ((SIN3B[i,2] > ATAC_peak$start & SIN3B[i,2] < ATAC_peak$end) | (SIN3B[i,3] > ATAC_peak$start & SIN3B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SIN3B) != 0){
    ATAC_peak_SIN3B$V5 <- SIN3B[i,1]
    ATAC_peak_SIN3B$V6 <- SIN3B[i,2]
    ATAC_peak_SIN3B$V7 <- SIN3B[i,3]
    ATAC_peak_SIN3B$V8 <- SIN3B[i,6]
    ATAC_peak_SIN3B$V9 <- SIN3B[i,7]
    ATAC_peak_SIN3B_list = rbind(ATAC_peak_SIN3B_list,ATAC_peak_SIN3B)
  }
}
write.table(ATAC_peak_SIN3B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIN3B_list.txt",sep="\t")

ATAC_peak_SIN3B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIN3B_list.txt",sep="\t")

SIN3B_features=list("SIN3B" = ATAC_peak_SIN3B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SIN3B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##STAT5B
STAT5B= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF043KIK_STAT5B_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(STAT5B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF043KIK_STAT5B_age.txt",sep = "\t")

STAT5B= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF043KIK_STAT5B_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(STAT5B)){
  x=paste(STAT5B[i,1],STAT5B[i,2],STAT5B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

STAT5B_features = list(list$V1)

ATAC_peak_STAT5B_list = data.frame()
for(i in 1:nrow(STAT5B)){
  ATAC_peak_STAT5B<- ATAC_peak[ATAC_peak$x1 == STAT5B[i,1] & ((STAT5B[i,2] > ATAC_peak$start & STAT5B[i,2] < ATAC_peak$end) | (STAT5B[i,3] > ATAC_peak$start & STAT5B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_STAT5B) != 0){
    ATAC_peak_STAT5B$V5 <- STAT5B[i,1]
    ATAC_peak_STAT5B$V6 <- STAT5B[i,2]
    ATAC_peak_STAT5B$V7 <- STAT5B[i,3]
    ATAC_peak_STAT5B$V8 <- STAT5B[i,6]
    ATAC_peak_STAT5B$V9 <- STAT5B[i,7]
    ATAC_peak_STAT5B_list = rbind(ATAC_peak_STAT5B_list,ATAC_peak_STAT5B)
  }
}
write.table(ATAC_peak_STAT5B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_STAT5B_list.txt",sep="\t")

ATAC_peak_STAT5B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_STAT5B_list.txt",sep="\t")

STAT5B_features=list("STAT5B" = ATAC_peak_STAT5B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = STAT5B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##MYC
MYC= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF664GSV_MYC_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MYC,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF664GSV_MYC_age.txt",sep = "\t")

MYC= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF664GSV_MYC_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MYC)){
  x=paste(MYC[i,1],MYC[i,2],MYC[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MYC_features = list(list$V1)

ATAC_peak_MYC_list = data.frame()
for(i in 1:nrow(MYC)){
  ATAC_peak_MYC<- ATAC_peak[ATAC_peak$x1 == MYC[i,1] & ((MYC[i,2] > ATAC_peak$start & MYC[i,2] < ATAC_peak$end) | (MYC[i,3] > ATAC_peak$start & MYC[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MYC) != 0){
    ATAC_peak_MYC$V5 <- MYC[i,1]
    ATAC_peak_MYC$V6 <- MYC[i,2]
    ATAC_peak_MYC$V7 <- MYC[i,3]
    ATAC_peak_MYC$V8 <- MYC[i,6]
    ATAC_peak_MYC$V9 <- MYC[i,7]
    ATAC_peak_MYC_list = rbind(ATAC_peak_MYC_list,ATAC_peak_MYC)
  }
}
write.table(ATAC_peak_MYC_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MYC_list.txt",sep="\t")

ATAC_peak_MYC_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MYC_list.txt",sep="\t")

MYC_features=list("MYC" = ATAC_peak_MYC_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MYC_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RBBP5
RBBP5= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF044FPH_RBBP5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RBBP5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF044FPH_RBBP5.txt",sep = "\t")

RBBP5= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF044FPH_RBBP5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RBBP5)){
  x=paste(RBBP5[i,1],RBBP5[i,2],RBBP5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RBBP5_features = list(list$V1)

ATAC_peak_RBBP5_list = data.frame()
for(i in 1:nrow(RBBP5)){
  ATAC_peak_RBBP5<- ATAC_peak[ATAC_peak$x1 == RBBP5[i,1] & ((RBBP5[i,2] > ATAC_peak$start & RBBP5[i,2] < ATAC_peak$end) | (RBBP5[i,3] > ATAC_peak$start & RBBP5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RBBP5) != 0){
    ATAC_peak_RBBP5$V5 <- RBBP5[i,1]
    ATAC_peak_RBBP5$V6 <- RBBP5[i,2]
    ATAC_peak_RBBP5$V7 <- RBBP5[i,3]
    ATAC_peak_RBBP5$V8 <- RBBP5[i,6]
    ATAC_peak_RBBP5$V9 <- RBBP5[i,7]
    ATAC_peak_RBBP5_list = rbind(ATAC_peak_RBBP5_list,ATAC_peak_RBBP5)
  }
}
write.table(ATAC_peak_RBBP5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RBBP5_list.txt",sep="\t")

ATAC_peak_RBBP5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RBBP5_list.txt",sep="\t")

RBBP5_features=list("RBBP5" = ATAC_peak_RBBP5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RBBP5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ARHGAP35
ARHGAP35= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF050CYN_ARHGAP35.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ARHGAP35,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF050CYN_ARHGAP35.txt",sep = "\t")

ARHGAP35= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF050CYN_ARHGAP35.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ARHGAP35)){
  x=paste(ARHGAP35[i,1],ARHGAP35[i,2],ARHGAP35[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ARHGAP35_features = list(list$V1)

ATAC_peak_ARHGAP35_list = data.frame()
for(i in 1:nrow(ARHGAP35)){
  ATAC_peak_ARHGAP35<- ATAC_peak[ATAC_peak$x1 == ARHGAP35[i,1] & ((ARHGAP35[i,2] > ATAC_peak$start & ARHGAP35[i,2] < ATAC_peak$end) | (ARHGAP35[i,3] > ATAC_peak$start & ARHGAP35[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ARHGAP35) != 0){
    ATAC_peak_ARHGAP35$V5 <- ARHGAP35[i,1]
    ATAC_peak_ARHGAP35$V6 <- ARHGAP35[i,2]
    ATAC_peak_ARHGAP35$V7 <- ARHGAP35[i,3]
    ATAC_peak_ARHGAP35$V8 <- ARHGAP35[i,6]
    ATAC_peak_ARHGAP35$V9 <- ARHGAP35[i,7]
    ATAC_peak_ARHGAP35_list = rbind(ATAC_peak_ARHGAP35_list,ATAC_peak_ARHGAP35)
  }
}
write.table(ATAC_peak_ARHGAP35_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARHGAP35_list.txt",sep="\t")

ATAC_peak_ARHGAP35_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARHGAP35_list.txt",sep="\t")

ARHGAP35_features=list("ARHGAP35" = ATAC_peak_ARHGAP35_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ARHGAP35_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HDAC3
HDAC3= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF051DMP_HDAC3_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HDAC3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF051DMP_HDAC3_age.txt",sep = "\t")

HDAC3= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF051DMP_HDAC3_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HDAC3)){
  x=paste(HDAC3[i,1],HDAC3[i,2],HDAC3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HDAC3_features = list(list$V1)

ATAC_peak_HDAC3_list = data.frame()
for(i in 1:nrow(HDAC3)){
  ATAC_peak_HDAC3<- ATAC_peak[ATAC_peak$x1 == HDAC3[i,1] & ((HDAC3[i,2] > ATAC_peak$start & HDAC3[i,2] < ATAC_peak$end) | (HDAC3[i,3] > ATAC_peak$start & HDAC3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HDAC3) != 0){
    ATAC_peak_HDAC3$V5 <- HDAC3[i,1]
    ATAC_peak_HDAC3$V6 <- HDAC3[i,2]
    ATAC_peak_HDAC3$V7 <- HDAC3[i,3]
    ATAC_peak_HDAC3$V8 <- HDAC3[i,6]
    ATAC_peak_HDAC3$V9 <- HDAC3[i,7]
    ATAC_peak_HDAC3_list = rbind(ATAC_peak_HDAC3_list,ATAC_peak_HDAC3)
  }
}
write.table(ATAC_peak_HDAC3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC3_list.txt",sep="\t")

ATAC_peak_HDAC3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC3_list.txt",sep="\t")

HDAC3_features=list("HDAC3" = ATAC_peak_HDAC3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HDAC3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SOX6
SOX6= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF060HPS_SOX6.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SOX6,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF060HPS_SOX6.txt",sep = "\t")

SOX6= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF060HPS_SOX6.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SOX6)){
  x=paste(SOX6[i,1],SOX6[i,2],SOX6[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SOX6_features = list(list$V1)

ATAC_peak_SOX6_list = data.frame()
for(i in 1:nrow(SOX6)){
  ATAC_peak_SOX6<- ATAC_peak[ATAC_peak$x1 == SOX6[i,1] & ((SOX6[i,2] > ATAC_peak$start & SOX6[i,2] < ATAC_peak$end) | (SOX6[i,3] > ATAC_peak$start & SOX6[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SOX6) != 0){
    ATAC_peak_SOX6$V5 <- SOX6[i,1]
    ATAC_peak_SOX6$V6 <- SOX6[i,2]
    ATAC_peak_SOX6$V7 <- SOX6[i,3]
    ATAC_peak_SOX6$V8 <- SOX6[i,6]
    ATAC_peak_SOX6$V9 <- SOX6[i,7]
    ATAC_peak_SOX6_list = rbind(ATAC_peak_SOX6_list,ATAC_peak_SOX6)
  }
}
write.table(ATAC_peak_SOX6_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SOX6_list.txt",sep="\t")

ATAC_peak_SOX6_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SOX6_list.txt",sep="\t")

SOX6_features=list("SOX6" = ATAC_peak_SOX6_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SOX6_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##GTF3C2
GTF3C2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF068ETD_GTF3C2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GTF3C2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF068ETD_GTF3C2.txt",sep = "\t")

GTF3C2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF068ETD_GTF3C2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GTF3C2)){
  x=paste(GTF3C2[i,1],GTF3C2[i,2],GTF3C2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GTF3C2_features = list(list$V1)

ATAC_peak_GTF3C2_list = data.frame()
for(i in 1:nrow(GTF3C2)){
  ATAC_peak_GTF3C2<- ATAC_peak[ATAC_peak$x1 == GTF3C2[i,1] & ((GTF3C2[i,2] > ATAC_peak$start & GTF3C2[i,2] < ATAC_peak$end) | (GTF3C2[i,3] > ATAC_peak$start & GTF3C2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GTF3C2) != 0){
    ATAC_peak_GTF3C2$V5 <- GTF3C2[i,1]
    ATAC_peak_GTF3C2$V6 <- GTF3C2[i,2]
    ATAC_peak_GTF3C2$V7 <- GTF3C2[i,3]
    ATAC_peak_GTF3C2$V8 <- GTF3C2[i,6]
    ATAC_peak_GTF3C2$V9 <- GTF3C2[i,7]
    ATAC_peak_GTF3C2_list = rbind(ATAC_peak_GTF3C2_list,ATAC_peak_GTF3C2)
  }
}
write.table(ATAC_peak_GTF3C2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF3C2_list.txt",sep="\t")

ATAC_peak_GTF3C2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF3C2_list.txt",sep="\t")

GTF3C2_features=list("GTF3C2" = ATAC_peak_GTF3C2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GTF3C2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CBX3
CBX3= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF068OEJ_CBX3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CBX3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF068OEJ_CBX3.txt",sep = "\t")

CBX3= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF068OEJ_CBX3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CBX3)){
  x=paste(CBX3[i,1],CBX3[i,2],CBX3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CBX3_features = list(list$V1)

ATAC_peak_CBX3_list = data.frame()
for(i in 1:nrow(CBX3)){
  ATAC_peak_CBX3<- ATAC_peak[ATAC_peak$x1 == CBX3[i,1] & ((CBX3[i,2] > ATAC_peak$start & CBX3[i,2] < ATAC_peak$end) | (CBX3[i,3] > ATAC_peak$start & CBX3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CBX3) != 0){
    ATAC_peak_CBX3$V5 <- CBX3[i,1]
    ATAC_peak_CBX3$V6 <- CBX3[i,2]
    ATAC_peak_CBX3$V7 <- CBX3[i,3]
    ATAC_peak_CBX3$V8 <- CBX3[i,6]
    ATAC_peak_CBX3$V9 <- CBX3[i,7]
    ATAC_peak_CBX3_list = rbind(ATAC_peak_CBX3_list,ATAC_peak_CBX3)
  }
}
write.table(ATAC_peak_CBX3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX3_list.txt",sep="\t")

ATAC_peak_CBX3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX3_list.txt",sep="\t")

CBX3_features=list("CBX3" = ATAC_peak_CBX3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CBX3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZNF507
ZNF507= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF072JDK_ZNF507.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF507,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF072JDK_ZNF507.txt",sep = "\t")

ZNF507= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF072JDK_ZNF507.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF507)){
  x=paste(ZNF507[i,1],ZNF507[i,2],ZNF507[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF507_features = list(list$V1)

ATAC_peak_ZNF507_list = data.frame()
for(i in 1:nrow(ZNF507)){
  ATAC_peak_ZNF507<- ATAC_peak[ATAC_peak$x1 == ZNF507[i,1] & ((ZNF507[i,2] > ATAC_peak$start & ZNF507[i,2] < ATAC_peak$end) | (ZNF507[i,3] > ATAC_peak$start & ZNF507[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF507) != 0){
    ATAC_peak_ZNF507$V5 <- ZNF507[i,1]
    ATAC_peak_ZNF507$V6 <- ZNF507[i,2]
    ATAC_peak_ZNF507$V7 <- ZNF507[i,3]
    ATAC_peak_ZNF507$V8 <- ZNF507[i,6]
    ATAC_peak_ZNF507$V9 <- ZNF507[i,7]
    ATAC_peak_ZNF507_list = rbind(ATAC_peak_ZNF507_list,ATAC_peak_ZNF507)
  }
}
write.table(ATAC_peak_ZNF507_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF507_list.txt",sep="\t")

ATAC_peak_ZNF507_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF507_list.txt",sep="\t")

ZNF507_features=list("ZNF507" = ATAC_peak_ZNF507_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF507_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##EZH2
EZH2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF080JPV_EZH2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(EZH2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF080JPV_EZH2.txt",sep = "\t")

EZH2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF080JPV_EZH2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(EZH2)){
  x=paste(EZH2[i,1],EZH2[i,2],EZH2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

EZH2_features = list(list$V1)

ATAC_peak_EZH2_list = data.frame()
for(i in 1:nrow(EZH2)){
  ATAC_peak_EZH2<- ATAC_peak[ATAC_peak$x1 == EZH2[i,1] & ((EZH2[i,2] > ATAC_peak$start & EZH2[i,2] < ATAC_peak$end) | (EZH2[i,3] > ATAC_peak$start & EZH2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_EZH2) != 0){
    ATAC_peak_EZH2$V5 <- EZH2[i,1]
    ATAC_peak_EZH2$V6 <- EZH2[i,2]
    ATAC_peak_EZH2$V7 <- EZH2[i,3]
    ATAC_peak_EZH2$V8 <- EZH2[i,6]
    ATAC_peak_EZH2$V9 <- EZH2[i,7]
    ATAC_peak_EZH2_list = rbind(ATAC_peak_EZH2_list,ATAC_peak_EZH2)
  }
}
write.table(ATAC_peak_EZH2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EZH2_list.txt",sep="\t")

ATAC_peak_EZH2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EZH2_list.txt",sep="\t")

EZH2_features=list("EZH2" = ATAC_peak_EZH2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = EZH2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ADNP
ADNP= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF083UZC_ADNP.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ADNP,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF083UZC_ADNP.txt",sep = "\t")

ADNP= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF083UZC_ADNP.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ADNP)){
  x=paste(ADNP[i,1],ADNP[i,2],ADNP[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ADNP_features = list(list$V1)

ATAC_peak_ADNP_list = data.frame()
for(i in 1:nrow(ADNP)){
  ATAC_peak_ADNP<- ATAC_peak[ATAC_peak$x1 == ADNP[i,1] & ((ADNP[i,2] > ATAC_peak$start & ADNP[i,2] < ATAC_peak$end) | (ADNP[i,3] > ATAC_peak$start & ADNP[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ADNP) != 0){
    ATAC_peak_ADNP$V5 <- ADNP[i,1]
    ATAC_peak_ADNP$V6 <- ADNP[i,2]
    ATAC_peak_ADNP$V7 <- ADNP[i,3]
    ATAC_peak_ADNP$V8 <- ADNP[i,6]
    ATAC_peak_ADNP$V9 <- ADNP[i,7]
    ATAC_peak_ADNP_list = rbind(ATAC_peak_ADNP_list,ATAC_peak_ADNP)
  }
}
write.table(ATAC_peak_ADNP_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ADNP_list.txt",sep="\t")

ATAC_peak_ADNP_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ADNP_list.txt",sep="\t")

ADNP_features=list("ADNP" = ATAC_peak_ADNP_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ADNP_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HMGN3
HMGN3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF092RLD_HMGN3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HMGN3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF092RLD_HMGN3.txt",sep = "\t")

HMGN3= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF092RLD_HMGN3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HMGN3)){
  x=paste(HMGN3[i,1],HMGN3[i,2],HMGN3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HMGN3_features = list(list$V1)

ATAC_peak_HMGN3_list = data.frame()
for(i in 1:nrow(HMGN3)){
  ATAC_peak_HMGN3<- ATAC_peak[ATAC_peak$x1 == HMGN3[i,1] & ((HMGN3[i,2] > ATAC_peak$start & HMGN3[i,2] < ATAC_peak$end) | (HMGN3[i,3] > ATAC_peak$start & HMGN3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HMGN3) != 0){
    ATAC_peak_HMGN3$V5 <- HMGN3[i,1]
    ATAC_peak_HMGN3$V6 <- HMGN3[i,2]
    ATAC_peak_HMGN3$V7 <- HMGN3[i,3]
    ATAC_peak_HMGN3$V8 <- HMGN3[i,6]
    ATAC_peak_HMGN3$V9 <- HMGN3[i,7]
    ATAC_peak_HMGN3_list = rbind(ATAC_peak_HMGN3_list,ATAC_peak_HMGN3)
  }
}
write.table(ATAC_peak_HMGN3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HMGN3_list.txt",sep="\t")

ATAC_peak_HMGN3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HMGN3_list.txt",sep="\t")

HMGN3_features=list("HMGN3" = ATAC_peak_HMGN3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HMGN3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##FOXO4
FOXO4 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF094TEN_FOXO4_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(FOXO4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF094TEN_FOXO4_age.txt",sep = "\t")

FOXO4= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF094TEN_FOXO4_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(FOXO4)){
  x=paste(FOXO4[i,1],FOXO4[i,2],FOXO4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

FOXO4_features = list(list$V1)

ATAC_peak_FOXO4_list = data.frame()
for(i in 1:nrow(FOXO4)){
  ATAC_peak_FOXO4<- ATAC_peak[ATAC_peak$x1 == FOXO4[i,1] & ((FOXO4[i,2] > ATAC_peak$start & FOXO4[i,2] < ATAC_peak$end) | (FOXO4[i,3] > ATAC_peak$start & FOXO4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_FOXO4) != 0){
    ATAC_peak_FOXO4$V5 <- FOXO4[i,1]
    ATAC_peak_FOXO4$V6 <- FOXO4[i,2]
    ATAC_peak_FOXO4$V7 <- FOXO4[i,3]
    ATAC_peak_FOXO4$V8 <- FOXO4[i,6]
    ATAC_peak_FOXO4$V9 <- FOXO4[i,7]
    ATAC_peak_FOXO4_list = rbind(ATAC_peak_FOXO4_list,ATAC_peak_FOXO4)
  }
}
write.table(ATAC_peak_FOXO4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXO4_list.txt",sep="\t")

ATAC_peak_FOXO4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXO4_list.txt",sep="\t")

FOXO4_features=list("FOXO4" = ATAC_peak_FOXO4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = FOXO4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NCOA4
NCOA4 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF096ZZW_NCOA4.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NCOA4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF096ZZW_NCOA4.txt",sep = "\t")

NCOA4= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF096ZZW_NCOA4.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NCOA4)){
  x=paste(NCOA4[i,1],NCOA4[i,2],NCOA4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NCOA4_features = list(list$V1)

ATAC_peak_NCOA4_list = data.frame()
for(i in 1:nrow(NCOA4)){
  ATAC_peak_NCOA4<- ATAC_peak[ATAC_peak$x1 == NCOA4[i,1] & ((NCOA4[i,2] > ATAC_peak$start & NCOA4[i,2] < ATAC_peak$end) | (NCOA4[i,3] > ATAC_peak$start & NCOA4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NCOA4) != 0){
    ATAC_peak_NCOA4$V5 <- NCOA4[i,1]
    ATAC_peak_NCOA4$V6 <- NCOA4[i,2]
    ATAC_peak_NCOA4$V7 <- NCOA4[i,3]
    ATAC_peak_NCOA4$V8 <- NCOA4[i,6]
    ATAC_peak_NCOA4$V9 <- NCOA4[i,7]
    ATAC_peak_NCOA4_list = rbind(ATAC_peak_NCOA4_list,ATAC_peak_NCOA4)
  }
}
write.table(ATAC_peak_NCOA4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOA4_list.txt",sep="\t")

ATAC_peak_NCOA4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOA4_list.txt",sep="\t")

NCOA4_features=list("NCOA4" = ATAC_peak_NCOA4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NCOA4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NKRF
NKRF = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF097NRE_NKRF.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NKRF,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF097NRE_NKRF.txt",sep = "\t")

NKRF= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF097NRE_NKRF.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NKRF)){
  x=paste(NKRF[i,1],NKRF[i,2],NKRF[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NKRF_features = list(list$V1)

ATAC_peak_NKRF_list = data.frame()
for(i in 1:nrow(NKRF)){
  ATAC_peak_NKRF<- ATAC_peak[ATAC_peak$x1 == NKRF[i,1] & ((NKRF[i,2] > ATAC_peak$start & NKRF[i,2] < ATAC_peak$end) | (NKRF[i,3] > ATAC_peak$start & NKRF[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NKRF) != 0){
    ATAC_peak_NKRF$V5 <- NKRF[i,1]
    ATAC_peak_NKRF$V6 <- NKRF[i,2]
    ATAC_peak_NKRF$V7 <- NKRF[i,3]
    ATAC_peak_NKRF$V8 <- NKRF[i,6]
    ATAC_peak_NKRF$V9 <- NKRF[i,7]
    ATAC_peak_NKRF_list = rbind(ATAC_peak_NKRF_list,ATAC_peak_NKRF)
  }
}
write.table(ATAC_peak_NKRF_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NKRF_list.txt",sep="\t")

ATAC_peak_NKRF_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NKRF_list.txt",sep="\t")

NKRF_features=list("NKRF" = ATAC_peak_NKRF_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NKRF_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##KDM1A
KDM1A = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF100CTO_KDM1A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(KDM1A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF100CTO_KDM1A.txt",sep = "\t")

KDM1A= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF100CTO_KDM1A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(KDM1A)){
  x=paste(KDM1A[i,1],KDM1A[i,2],KDM1A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

KDM1A_features = list(list$V1)

ATAC_peak_KDM1A_list = data.frame()
for(i in 1:nrow(KDM1A)){
  ATAC_peak_KDM1A<- ATAC_peak[ATAC_peak$x1 == KDM1A[i,1] & ((KDM1A[i,2] > ATAC_peak$start & KDM1A[i,2] < ATAC_peak$end) | (KDM1A[i,3] > ATAC_peak$start & KDM1A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_KDM1A) != 0){
    ATAC_peak_KDM1A$V5 <- KDM1A[i,1]
    ATAC_peak_KDM1A$V6 <- KDM1A[i,2]
    ATAC_peak_KDM1A$V7 <- KDM1A[i,3]
    ATAC_peak_KDM1A$V8 <- KDM1A[i,6]
    ATAC_peak_KDM1A$V9 <- KDM1A[i,7]
    ATAC_peak_KDM1A_list = rbind(ATAC_peak_KDM1A_list,ATAC_peak_KDM1A)
  }
}
write.table(ATAC_peak_KDM1A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KDM1A_list.txt",sep="\t")

ATAC_peak_KDM1A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KDM1A_list.txt",sep="\t")

KDM1A_features=list("KDM1A" = ATAC_peak_KDM1A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = KDM1A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##PCBP1
PCBP1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF103HVZ_PCBP1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PCBP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF103HVZ_PCBP1.txt",sep = "\t")

PCBP1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF103HVZ_PCBP1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PCBP1)){
  x=paste(PCBP1[i,1],PCBP1[i,2],PCBP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PCBP1_features = list(list$V1)

ATAC_peak_PCBP1_list = data.frame()
for(i in 1:nrow(PCBP1)){
  ATAC_peak_PCBP1<- ATAC_peak[ATAC_peak$x1 == PCBP1[i,1] & ((PCBP1[i,2] > ATAC_peak$start & PCBP1[i,2] < ATAC_peak$end) | (PCBP1[i,3] > ATAC_peak$start & PCBP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PCBP1) != 0){
    ATAC_peak_PCBP1$V5 <- PCBP1[i,1]
    ATAC_peak_PCBP1$V6 <- PCBP1[i,2]
    ATAC_peak_PCBP1$V7 <- PCBP1[i,3]
    ATAC_peak_PCBP1$V8 <- PCBP1[i,6]
    ATAC_peak_PCBP1$V9 <- PCBP1[i,7]
    ATAC_peak_PCBP1_list = rbind(ATAC_peak_PCBP1_list,ATAC_peak_PCBP1)
  }
}
write.table(ATAC_peak_PCBP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PCBP1_list.txt",sep="\t")

ATAC_peak_PCBP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PCBP1_list.txt",sep="\t")

PCBP1_features=list("PCBP1" = ATAC_peak_PCBP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PCBP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##FIP1L1
FIP1L1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF105YIV_FIP1L1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(FIP1L1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF105YIV_FIP1L1.txt",sep = "\t")

FIP1L1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF105YIV_FIP1L1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(FIP1L1)){
  x=paste(FIP1L1[i,1],FIP1L1[i,2],FIP1L1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

FIP1L1_features = list(list$V1)

ATAC_peak_FIP1L1_list = data.frame()
for(i in 1:nrow(FIP1L1)){
  ATAC_peak_FIP1L1<- ATAC_peak[ATAC_peak$x1 == FIP1L1[i,1] & ((FIP1L1[i,2] > ATAC_peak$start & FIP1L1[i,2] < ATAC_peak$end) | (FIP1L1[i,3] > ATAC_peak$start & FIP1L1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_FIP1L1) != 0){
    ATAC_peak_FIP1L1$V5 <- FIP1L1[i,1]
    ATAC_peak_FIP1L1$V6 <- FIP1L1[i,2]
    ATAC_peak_FIP1L1$V7 <- FIP1L1[i,3]
    ATAC_peak_FIP1L1$V8 <- FIP1L1[i,6]
    ATAC_peak_FIP1L1$V9 <- FIP1L1[i,7]
    ATAC_peak_FIP1L1_list = rbind(ATAC_peak_FIP1L1_list,ATAC_peak_FIP1L1)
  }
}
write.table(ATAC_peak_FIP1L1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FIP1L1_list.txt",sep="\t")

ATAC_peak_FIP1L1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FIP1L1_list.txt",sep="\t")

FIP1L1_features=list("FIP1L1" = ATAC_peak_FIP1L1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = FIP1L1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TOE1
TOE1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF108UNO_TOE1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TOE1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF108UNO_TOE1.txt",sep = "\t")

TOE1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF108UNO_TOE1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TOE1)){
  x=paste(TOE1[i,1],TOE1[i,2],TOE1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TOE1_features = list(list$V1)

ATAC_peak_TOE1_list = data.frame()
for(i in 1:nrow(TOE1)){
  ATAC_peak_TOE1<- ATAC_peak[ATAC_peak$x1 == TOE1[i,1] & ((TOE1[i,2] > ATAC_peak$start & TOE1[i,2] < ATAC_peak$end) | (TOE1[i,3] > ATAC_peak$start & TOE1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TOE1) != 0){
    ATAC_peak_TOE1$V5 <- TOE1[i,1]
    ATAC_peak_TOE1$V6 <- TOE1[i,2]
    ATAC_peak_TOE1$V7 <- TOE1[i,3]
    ATAC_peak_TOE1$V8 <- TOE1[i,6]
    ATAC_peak_TOE1$V9 <- TOE1[i,7]
    ATAC_peak_TOE1_list = rbind(ATAC_peak_TOE1_list,ATAC_peak_TOE1)
  }
}
write.table(ATAC_peak_TOE1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TOE1_list.txt",sep="\t")

ATAC_peak_TOE1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TOE1_list.txt",sep="\t")

TOE1_features=list("TOE1" = ATAC_peak_TOE1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TOE1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZFP91
ZFP91 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF110DOA_ZFP91.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZFP91,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF110DOA_ZFP91.txt",sep = "\t")

ZFP91= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF110DOA_ZFP91.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZFP91)){
  x=paste(ZFP91[i,1],ZFP91[i,2],ZFP91[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZFP91_features = list(list$V1)

ATAC_peak_ZFP91_list = data.frame()
for(i in 1:nrow(ZFP91)){
  ATAC_peak_ZFP91<- ATAC_peak[ATAC_peak$x1 == ZFP91[i,1] & ((ZFP91[i,2] > ATAC_peak$start & ZFP91[i,2] < ATAC_peak$end) | (ZFP91[i,3] > ATAC_peak$start & ZFP91[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZFP91) != 0){
    ATAC_peak_ZFP91$V5 <- ZFP91[i,1]
    ATAC_peak_ZFP91$V6 <- ZFP91[i,2]
    ATAC_peak_ZFP91$V7 <- ZFP91[i,3]
    ATAC_peak_ZFP91$V8 <- ZFP91[i,6]
    ATAC_peak_ZFP91$V9 <- ZFP91[i,7]
    ATAC_peak_ZFP91_list = rbind(ATAC_peak_ZFP91_list,ATAC_peak_ZFP91)
  }
}
write.table(ATAC_peak_ZFP91_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZFP91_list.txt",sep="\t")

ATAC_peak_ZFP91_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZFP91_list.txt",sep="\t")

ZFP91_features=list("ZFP91" = ATAC_peak_ZFP91_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZFP91_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NUFIP1
NUFIP1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF130GVU_NUFIP1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NUFIP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF130GVU_NUFIP1.txt",sep = "\t")

NUFIP1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF130GVU_NUFIP1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NUFIP1)){
  x=paste(NUFIP1[i,1],NUFIP1[i,2],NUFIP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NUFIP1_features = list(list$V1)

ATAC_peak_NUFIP1_list = data.frame()
for(i in 1:nrow(NUFIP1)){
  ATAC_peak_NUFIP1<- ATAC_peak[ATAC_peak$x1 == NUFIP1[i,1] & ((NUFIP1[i,2] > ATAC_peak$start & NUFIP1[i,2] < ATAC_peak$end) | (NUFIP1[i,3] > ATAC_peak$start & NUFIP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NUFIP1) != 0){
    ATAC_peak_NUFIP1$V5 <- NUFIP1[i,1]
    ATAC_peak_NUFIP1$V6 <- NUFIP1[i,2]
    ATAC_peak_NUFIP1$V7 <- NUFIP1[i,3]
    ATAC_peak_NUFIP1$V8 <- NUFIP1[i,6]
    ATAC_peak_NUFIP1$V9 <- NUFIP1[i,7]
    ATAC_peak_NUFIP1_list = rbind(ATAC_peak_NUFIP1_list,ATAC_peak_NUFIP1)
  }
}
write.table(ATAC_peak_NUFIP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NUFIP1_list.txt",sep="\t")

ATAC_peak_NUFIP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NUFIP1_list.txt",sep="\t")

NUFIP1_features=list("NUFIP1" = ATAC_peak_NUFIP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NUFIP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NUFIP1
NUFIP1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF130GVU_NUFIP1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NUFIP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF130GVU_NUFIP1.txt",sep = "\t")

NUFIP1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF130GVU_NUFIP1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NUFIP1)){
  x=paste(NUFIP1[i,1],NUFIP1[i,2],NUFIP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NUFIP1_features = list(list$V1)

ATAC_peak_NUFIP1_list = data.frame()
for(i in 1:nrow(NUFIP1)){
  ATAC_peak_NUFIP1<- ATAC_peak[ATAC_peak$x1 == NUFIP1[i,1] & ((NUFIP1[i,2] > ATAC_peak$start & NUFIP1[i,2] < ATAC_peak$end) | (NUFIP1[i,3] > ATAC_peak$start & NUFIP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NUFIP1) != 0){
    ATAC_peak_NUFIP1$V5 <- NUFIP1[i,1]
    ATAC_peak_NUFIP1$V6 <- NUFIP1[i,2]
    ATAC_peak_NUFIP1$V7 <- NUFIP1[i,3]
    ATAC_peak_NUFIP1$V8 <- NUFIP1[i,6]
    ATAC_peak_NUFIP1$V9 <- NUFIP1[i,7]
    ATAC_peak_NUFIP1_list = rbind(ATAC_peak_NUFIP1_list,ATAC_peak_NUFIP1)
  }
}
write.table(ATAC_peak_NUFIP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NUFIP1_list.txt",sep="\t")

ATAC_peak_NUFIP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NUFIP1_list.txt",sep="\t")

NUFIP1_features=list("NUFIP1" = ATAC_peak_NUFIP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NUFIP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##AGO1
AGO1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF132NIO_AGO1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(AGO1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF132NIO_AGO1.txt",sep = "\t")

AGO1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF132NIO_AGO1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(AGO1)){
  x=paste(AGO1[i,1],AGO1[i,2],AGO1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

AGO1_features = list(list$V1)

ATAC_peak_AGO1_list = data.frame()
for(i in 1:nrow(AGO1)){
  ATAC_peak_AGO1<- ATAC_peak[ATAC_peak$x1 == AGO1[i,1] & ((AGO1[i,2] > ATAC_peak$start & AGO1[i,2] < ATAC_peak$end) | (AGO1[i,3] > ATAC_peak$start & AGO1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_AGO1) != 0){
    ATAC_peak_AGO1$V5 <- AGO1[i,1]
    ATAC_peak_AGO1$V6 <- AGO1[i,2]
    ATAC_peak_AGO1$V7 <- AGO1[i,3]
    ATAC_peak_AGO1$V8 <- AGO1[i,6]
    ATAC_peak_AGO1$V9 <- AGO1[i,7]
    ATAC_peak_AGO1_list = rbind(ATAC_peak_AGO1_list,ATAC_peak_AGO1)
  }
}
write.table(ATAC_peak_AGO1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_AGO1_list.txt",sep="\t")

ATAC_peak_AGO1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_AGO1_list.txt",sep="\t")

AGO1_features=list("AGO1" = ATAC_peak_AGO1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = AGO1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RBFOX2
RBFOX2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF140GXH_RBFOX2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RBFOX2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF140GXH_RBFOX2.txt",sep = "\t")

RBFOX2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF140GXH_RBFOX2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RBFOX2)){
  x=paste(RBFOX2[i,1],RBFOX2[i,2],RBFOX2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RBFOX2_features = list(list$V1)

ATAC_peak_RBFOX2_list = data.frame()
for(i in 1:nrow(RBFOX2)){
  ATAC_peak_RBFOX2<- ATAC_peak[ATAC_peak$x1 == RBFOX2[i,1] & ((RBFOX2[i,2] > ATAC_peak$start & RBFOX2[i,2] < ATAC_peak$end) | (RBFOX2[i,3] > ATAC_peak$start & RBFOX2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RBFOX2) != 0){
    ATAC_peak_RBFOX2$V5 <- RBFOX2[i,1]
    ATAC_peak_RBFOX2$V6 <- RBFOX2[i,2]
    ATAC_peak_RBFOX2$V7 <- RBFOX2[i,3]
    ATAC_peak_RBFOX2$V8 <- RBFOX2[i,6]
    ATAC_peak_RBFOX2$V9 <- RBFOX2[i,7]
    ATAC_peak_RBFOX2_list = rbind(ATAC_peak_RBFOX2_list,ATAC_peak_RBFOX2)
  }
}
write.table(ATAC_peak_RBFOX2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RBFOX2_list.txt",sep="\t")

ATAC_peak_RBFOX2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RBFOX2_list.txt",sep="\t")

RBFOX2_features=list("RBFOX2" = ATAC_peak_RBFOX2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RBFOX2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##PRDM10
PRDM10 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF141RMX_PRDM10.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PRDM10,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF141RMX_PRDM10.txt",sep = "\t")

PRDM10= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF141RMX_PRDM10.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PRDM10)){
  x=paste(PRDM10[i,1],PRDM10[i,2],PRDM10[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PRDM10_features = list(list$V1)

ATAC_peak_PRDM10_list = data.frame()
for(i in 1:nrow(PRDM10)){
  ATAC_peak_PRDM10<- ATAC_peak[ATAC_peak$x1 == PRDM10[i,1] & ((PRDM10[i,2] > ATAC_peak$start & PRDM10[i,2] < ATAC_peak$end) | (PRDM10[i,3] > ATAC_peak$start & PRDM10[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PRDM10) != 0){
    ATAC_peak_PRDM10$V5 <- PRDM10[i,1]
    ATAC_peak_PRDM10$V6 <- PRDM10[i,2]
    ATAC_peak_PRDM10$V7 <- PRDM10[i,3]
    ATAC_peak_PRDM10$V8 <- PRDM10[i,6]
    ATAC_peak_PRDM10$V9 <- PRDM10[i,7]
    ATAC_peak_PRDM10_list = rbind(ATAC_peak_PRDM10_list,ATAC_peak_PRDM10)
  }
}
write.table(ATAC_peak_PRDM10_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PRDM10_list.txt",sep="\t")

ATAC_peak_PRDM10_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PRDM10_list.txt",sep="\t")

PRDM10_features=list("PRDM10" = ATAC_peak_PRDM10_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PRDM10_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##THRAP3
THRAP3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF147MOG_THRAP3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(THRAP3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF147MOG_THRAP3.txt",sep = "\t")

THRAP3= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF147MOG_THRAP3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(THRAP3)){
  x=paste(THRAP3[i,1],THRAP3[i,2],THRAP3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

THRAP3_features = list(list$V1)

ATAC_peak_THRAP3_list = data.frame()
for(i in 1:nrow(THRAP3)){
  ATAC_peak_THRAP3<- ATAC_peak[ATAC_peak$x1 == THRAP3[i,1] & ((THRAP3[i,2] > ATAC_peak$start & THRAP3[i,2] < ATAC_peak$end) | (THRAP3[i,3] > ATAC_peak$start & THRAP3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_THRAP3) != 0){
    ATAC_peak_THRAP3$V5 <- THRAP3[i,1]
    ATAC_peak_THRAP3$V6 <- THRAP3[i,2]
    ATAC_peak_THRAP3$V7 <- THRAP3[i,3]
    ATAC_peak_THRAP3$V8 <- THRAP3[i,6]
    ATAC_peak_THRAP3$V9 <- THRAP3[i,7]
    ATAC_peak_THRAP3_list = rbind(ATAC_peak_THRAP3_list,ATAC_peak_THRAP3)
  }
}
write.table(ATAC_peak_THRAP3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_THRAP3_list.txt",sep="\t")

ATAC_peak_THRAP3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_THRAP3_list.txt",sep="\t")

THRAP3_features=list("THRAP3" = ATAC_peak_THRAP3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = THRAP3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NR2C1
NR2C1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF151VRY_NR2C1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NR2C1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF151VRY_NR2C1.txt",sep = "\t")

NR2C1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF151VRY_NR2C1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NR2C1)){
  x=paste(NR2C1[i,1],NR2C1[i,2],NR2C1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NR2C1_features = list(list$V1)

ATAC_peak_NR2C1_list = data.frame()
for(i in 1:nrow(NR2C1)){
  ATAC_peak_NR2C1<- ATAC_peak[ATAC_peak$x1 == NR2C1[i,1] & ((NR2C1[i,2] > ATAC_peak$start & NR2C1[i,2] < ATAC_peak$end) | (NR2C1[i,3] > ATAC_peak$start & NR2C1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NR2C1) != 0){
    ATAC_peak_NR2C1$V5 <- NR2C1[i,1]
    ATAC_peak_NR2C1$V6 <- NR2C1[i,2]
    ATAC_peak_NR2C1$V7 <- NR2C1[i,3]
    ATAC_peak_NR2C1$V8 <- NR2C1[i,6]
    ATAC_peak_NR2C1$V9 <- NR2C1[i,7]
    ATAC_peak_NR2C1_list = rbind(ATAC_peak_NR2C1_list,ATAC_peak_NR2C1)
  }
}
write.table(ATAC_peak_NR2C1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NR2C1_list.txt",sep="\t")

ATAC_peak_NR2C1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NR2C1_list.txt",sep="\t")

NR2C1_features=list("NR2C1" = ATAC_peak_NR2C1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NR2C1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CBX8
CBX8 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF157NCP_CBX8.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CBX8,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF157NCP_CBX8.txt",sep = "\t")

CBX8= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF157NCP_CBX8.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CBX8)){
  x=paste(CBX8[i,1],CBX8[i,2],CBX8[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CBX8_features = list(list$V1)

ATAC_peak_CBX8_list = data.frame()
for(i in 1:nrow(CBX8)){
  ATAC_peak_CBX8<- ATAC_peak[ATAC_peak$x1 == CBX8[i,1] & ((CBX8[i,2] > ATAC_peak$start & CBX8[i,2] < ATAC_peak$end) | (CBX8[i,3] > ATAC_peak$start & CBX8[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CBX8) != 0){
    ATAC_peak_CBX8$V5 <- CBX8[i,1]
    ATAC_peak_CBX8$V6 <- CBX8[i,2]
    ATAC_peak_CBX8$V7 <- CBX8[i,3]
    ATAC_peak_CBX8$V8 <- CBX8[i,6]
    ATAC_peak_CBX8$V9 <- CBX8[i,7]
    ATAC_peak_CBX8_list = rbind(ATAC_peak_CBX8_list,ATAC_peak_CBX8)
  }
}
write.table(ATAC_peak_CBX8_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX8_list.txt",sep="\t")

ATAC_peak_CBX8_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX8_list.txt",sep="\t")

CBX8_features=list("CBX8" = ATAC_peak_CBX8_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CBX8_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TBP
TBP = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF159SBO_TBP_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TBP,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF159SBO_TBP_age.txt",sep = "\t")

TBP= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF159SBO_TBP_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TBP)){
  x=paste(TBP[i,1],TBP[i,2],TBP[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TBP_features = list(list$V1)

ATAC_peak_TBP_list = data.frame()
for(i in 1:nrow(TBP)){
  ATAC_peak_TBP<- ATAC_peak[ATAC_peak$x1 == TBP[i,1] & ((TBP[i,2] > ATAC_peak$start & TBP[i,2] < ATAC_peak$end) | (TBP[i,3] > ATAC_peak$start & TBP[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TBP) != 0){
    ATAC_peak_TBP$V5 <- TBP[i,1]
    ATAC_peak_TBP$V6 <- TBP[i,2]
    ATAC_peak_TBP$V7 <- TBP[i,3]
    ATAC_peak_TBP$V8 <- TBP[i,6]
    ATAC_peak_TBP$V9 <- TBP[i,7]
    ATAC_peak_TBP_list = rbind(ATAC_peak_TBP_list,ATAC_peak_TBP)
  }
}
write.table(ATAC_peak_TBP_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TBP_list.txt",sep="\t")

ATAC_peak_TBP_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TBP_list.txt",sep="\t")

TBP_features=list("TBP" = ATAC_peak_TBP_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TBP_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CEBPZ
CEBPZ = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF161IQF_CEBPZ.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CEBPZ,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF161IQF_CEBPZ.txt",sep = "\t")

CEBPZ= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF161IQF_CEBPZ.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CEBPZ)){
  x=paste(CEBPZ[i,1],CEBPZ[i,2],CEBPZ[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CEBPZ_features = list(list$V1)

ATAC_peak_CEBPZ_list = data.frame()
for(i in 1:nrow(CEBPZ)){
  ATAC_peak_CEBPZ<- ATAC_peak[ATAC_peak$x1 == CEBPZ[i,1] & ((CEBPZ[i,2] > ATAC_peak$start & CEBPZ[i,2] < ATAC_peak$end) | (CEBPZ[i,3] > ATAC_peak$start & CEBPZ[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CEBPZ) != 0){
    ATAC_peak_CEBPZ$V5 <- CEBPZ[i,1]
    ATAC_peak_CEBPZ$V6 <- CEBPZ[i,2]
    ATAC_peak_CEBPZ$V7 <- CEBPZ[i,3]
    ATAC_peak_CEBPZ$V8 <- CEBPZ[i,6]
    ATAC_peak_CEBPZ$V9 <- CEBPZ[i,7]
    ATAC_peak_CEBPZ_list = rbind(ATAC_peak_CEBPZ_list,ATAC_peak_CEBPZ)
  }
}
write.table(ATAC_peak_CEBPZ_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CEBPZ_list.txt",sep="\t")

ATAC_peak_CEBPZ_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CEBPZ_list.txt",sep="\t")

CEBPZ_features=list("CEBPZ" = ATAC_peak_CEBPZ_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CEBPZ_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZC3H11A
ZC3H11A = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF163FGU_ZC3H11A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZC3H11A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF163FGU_ZC3H11A.txt",sep = "\t")

ZC3H11A= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF163FGU_ZC3H11A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZC3H11A)){
  x=paste(ZC3H11A[i,1],ZC3H11A[i,2],ZC3H11A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZC3H11A_features = list(list$V1)

ATAC_peak_ZC3H11A_list = data.frame()
for(i in 1:nrow(ZC3H11A)){
  ATAC_peak_ZC3H11A<- ATAC_peak[ATAC_peak$x1 == ZC3H11A[i,1] & ((ZC3H11A[i,2] > ATAC_peak$start & ZC3H11A[i,2] < ATAC_peak$end) | (ZC3H11A[i,3] > ATAC_peak$start & ZC3H11A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZC3H11A) != 0){
    ATAC_peak_ZC3H11A$V5 <- ZC3H11A[i,1]
    ATAC_peak_ZC3H11A$V6 <- ZC3H11A[i,2]
    ATAC_peak_ZC3H11A$V7 <- ZC3H11A[i,3]
    ATAC_peak_ZC3H11A$V8 <- ZC3H11A[i,6]
    ATAC_peak_ZC3H11A$V9 <- ZC3H11A[i,7]
    ATAC_peak_ZC3H11A_list = rbind(ATAC_peak_ZC3H11A_list,ATAC_peak_ZC3H11A)
  }
}
write.table(ATAC_peak_ZC3H11A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZC3H11A_list.txt",sep="\t")

ATAC_peak_ZC3H11A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZC3H11A_list.txt",sep="\t")

ZC3H11A_features=list("ZC3H11A" = ATAC_peak_ZC3H11A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZC3H11A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HDGF
HDGF = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF166CUH_HDGF.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HDGF,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF166CUH_HDGF.txt",sep = "\t")

HDGF= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF166CUH_HDGF.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HDGF)){
  x=paste(HDGF[i,1],HDGF[i,2],HDGF[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HDGF_features = list(list$V1)

ATAC_peak_HDGF_list = data.frame()
for(i in 1:nrow(HDGF)){
  ATAC_peak_HDGF<- ATAC_peak[ATAC_peak$x1 == HDGF[i,1] & ((HDGF[i,2] > ATAC_peak$start & HDGF[i,2] < ATAC_peak$end) | (HDGF[i,3] > ATAC_peak$start & HDGF[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HDGF) != 0){
    ATAC_peak_HDGF$V5 <- HDGF[i,1]
    ATAC_peak_HDGF$V6 <- HDGF[i,2]
    ATAC_peak_HDGF$V7 <- HDGF[i,3]
    ATAC_peak_HDGF$V8 <- HDGF[i,6]
    ATAC_peak_HDGF$V9 <- HDGF[i,7]
    ATAC_peak_HDGF_list = rbind(ATAC_peak_HDGF_list,ATAC_peak_HDGF)
  }
}
write.table(ATAC_peak_HDGF_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDGF_list.txt",sep="\t")

ATAC_peak_HDGF_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDGF_list.txt",sep="\t")

HDGF_features=list("HDGF" = ATAC_peak_HDGF_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HDGF_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##THRA
THRA = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF174DRC_THRA.bed.gz",extraCols=extraCols_narrowPeak)
write.table(THRA,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF174DRC_THRA.txt",sep = "\t")

THRA= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF174DRC_THRA.txt",sep = "\t")

list = vector()
for (i in  1:nrow(THRA)){
  x=paste(THRA[i,1],THRA[i,2],THRA[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

THRA_features = list(list$V1)

ATAC_peak_THRA_list = data.frame()
for(i in 1:nrow(THRA)){
  ATAC_peak_THRA<- ATAC_peak[ATAC_peak$x1 == THRA[i,1] & ((THRA[i,2] > ATAC_peak$start & THRA[i,2] < ATAC_peak$end) | (THRA[i,3] > ATAC_peak$start & THRA[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_THRA) != 0){
    ATAC_peak_THRA$V5 <- THRA[i,1]
    ATAC_peak_THRA$V6 <- THRA[i,2]
    ATAC_peak_THRA$V7 <- THRA[i,3]
    ATAC_peak_THRA$V8 <- THRA[i,6]
    ATAC_peak_THRA$V9 <- THRA[i,7]
    ATAC_peak_THRA_list = rbind(ATAC_peak_THRA_list,ATAC_peak_THRA)
  }
}
write.table(ATAC_peak_THRA_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_THRA_list.txt",sep="\t")

ATAC_peak_THRA_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_THRA_list.txt",sep="\t")

THRA_features=list("THRA" = ATAC_peak_THRA_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = THRA_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZHX1
ZHX1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF175YPM_ZHX1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZHX1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF175YPM_ZHX1.txt",sep = "\t")

ZHX1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF175YPM_ZHX1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZHX1)){
  x=paste(ZHX1[i,1],ZHX1[i,2],ZHX1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZHX1_features = list(list$V1)

ATAC_peak_ZHX1_list = data.frame()
for(i in 1:nrow(ZHX1)){
  ATAC_peak_ZHX1<- ATAC_peak[ATAC_peak$x1 == ZHX1[i,1] & ((ZHX1[i,2] > ATAC_peak$start & ZHX1[i,2] < ATAC_peak$end) | (ZHX1[i,3] > ATAC_peak$start & ZHX1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZHX1) != 0){
    ATAC_peak_ZHX1$V5 <- ZHX1[i,1]
    ATAC_peak_ZHX1$V6 <- ZHX1[i,2]
    ATAC_peak_ZHX1$V7 <- ZHX1[i,3]
    ATAC_peak_ZHX1$V8 <- ZHX1[i,6]
    ATAC_peak_ZHX1$V9 <- ZHX1[i,7]
    ATAC_peak_ZHX1_list = rbind(ATAC_peak_ZHX1_list,ATAC_peak_ZHX1)
  }
}
write.table(ATAC_peak_ZHX1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZHX1_list.txt",sep="\t")

ATAC_peak_ZHX1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZHX1_list.txt",sep="\t")

ZHX1_features=list("ZHX1" = ATAC_peak_ZHX1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZHX1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SMARCA4
SMARCA4 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF193CLQ_SMARCA4.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMARCA4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF193CLQ_SMARCA4.txt",sep = "\t")

SMARCA4= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF193CLQ_SMARCA4.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMARCA4)){
  x=paste(SMARCA4[i,1],SMARCA4[i,2],SMARCA4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMARCA4_features = list(list$V1)

ATAC_peak_SMARCA4_list = data.frame()
for(i in 1:nrow(SMARCA4)){
  ATAC_peak_SMARCA4<- ATAC_peak[ATAC_peak$x1 == SMARCA4[i,1] & ((SMARCA4[i,2] > ATAC_peak$start & SMARCA4[i,2] < ATAC_peak$end) | (SMARCA4[i,3] > ATAC_peak$start & SMARCA4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMARCA4) != 0){
    ATAC_peak_SMARCA4$V5 <- SMARCA4[i,1]
    ATAC_peak_SMARCA4$V6 <- SMARCA4[i,2]
    ATAC_peak_SMARCA4$V7 <- SMARCA4[i,3]
    ATAC_peak_SMARCA4$V8 <- SMARCA4[i,6]
    ATAC_peak_SMARCA4$V9 <- SMARCA4[i,7]
    ATAC_peak_SMARCA4_list = rbind(ATAC_peak_SMARCA4_list,ATAC_peak_SMARCA4)
  }
}
write.table(ATAC_peak_SMARCA4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMARCA4_list.txt",sep="\t")

ATAC_peak_SMARCA4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMARCA4_list.txt",sep="\t")

SMARCA4_features=list("SMARCA4" = ATAC_peak_SMARCA4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMARCA4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZNF280A
ZNF280A = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF184TZK_ZNF280A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF280A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF184TZK_ZNF280A.txt",sep = "\t")

ZNF280A= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF184TZK_ZNF280A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF280A)){
  x=paste(ZNF280A[i,1],ZNF280A[i,2],ZNF280A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF280A_features = list(list$V1)

ATAC_peak_ZNF280A_list = data.frame()
for(i in 1:nrow(ZNF280A)){
  ATAC_peak_ZNF280A<- ATAC_peak[ATAC_peak$x1 == ZNF280A[i,1] & ((ZNF280A[i,2] > ATAC_peak$start & ZNF280A[i,2] < ATAC_peak$end) | (ZNF280A[i,3] > ATAC_peak$start & ZNF280A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF280A) != 0){
    ATAC_peak_ZNF280A$V5 <- ZNF280A[i,1]
    ATAC_peak_ZNF280A$V6 <- ZNF280A[i,2]
    ATAC_peak_ZNF280A$V7 <- ZNF280A[i,3]
    ATAC_peak_ZNF280A$V8 <- ZNF280A[i,6]
    ATAC_peak_ZNF280A$V9 <- ZNF280A[i,7]
    ATAC_peak_ZNF280A_list = rbind(ATAC_peak_ZNF280A_list,ATAC_peak_ZNF280A)
  }
}
write.table(ATAC_peak_ZNF280A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF280A_list.txt",sep="\t")

ATAC_peak_ZNF280A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF280A_list.txt",sep="\t")

ZNF280A_features=list("ZNF280A" = ATAC_peak_ZNF280A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF280A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CTBP1
CTBP1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF200APZ_CTBP1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CTBP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF200APZ_CTBP1.txt",sep = "\t")

CTBP1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF200APZ_CTBP1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CTBP1)){
  x=paste(CTBP1[i,1],CTBP1[i,2],CTBP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CTBP1_features = list(list$V1)

ATAC_peak_CTBP1_list = data.frame()
for(i in 1:nrow(CTBP1)){
  ATAC_peak_CTBP1<- ATAC_peak[ATAC_peak$x1 == CTBP1[i,1] & ((CTBP1[i,2] > ATAC_peak$start & CTBP1[i,2] < ATAC_peak$end) | (CTBP1[i,3] > ATAC_peak$start & CTBP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CTBP1) != 0){
    ATAC_peak_CTBP1$V5 <- CTBP1[i,1]
    ATAC_peak_CTBP1$V6 <- CTBP1[i,2]
    ATAC_peak_CTBP1$V7 <- CTBP1[i,3]
    ATAC_peak_CTBP1$V8 <- CTBP1[i,6]
    ATAC_peak_CTBP1$V9 <- CTBP1[i,7]
    ATAC_peak_CTBP1_list = rbind(ATAC_peak_CTBP1_list,ATAC_peak_CTBP1)
  }
}
write.table(ATAC_peak_CTBP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CTBP1_list.txt",sep="\t")

ATAC_peak_CTBP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CTBP1_list.txt",sep="\t")

CTBP1_features=list("CTBP1" = ATAC_peak_CTBP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CTBP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ATF1
ATF1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF221MGJ_ATF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ATF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF221MGJ_ATF1.txt",sep = "\t")

ATF1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF221MGJ_ATF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ATF1)){
  x=paste(ATF1[i,1],ATF1[i,2],ATF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ATF1_features = list(list$V1)

ATAC_peak_ATF1_list = data.frame()
for(i in 1:nrow(ATF1)){
  ATAC_peak_ATF1<- ATAC_peak[ATAC_peak$x1 == ATF1[i,1] & ((ATF1[i,2] > ATAC_peak$start & ATF1[i,2] < ATAC_peak$end) | (ATF1[i,3] > ATAC_peak$start & ATF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ATF1) != 0){
    ATAC_peak_ATF1$V5 <- ATF1[i,1]
    ATAC_peak_ATF1$V6 <- ATF1[i,2]
    ATAC_peak_ATF1$V7 <- ATF1[i,3]
    ATAC_peak_ATF1$V8 <- ATF1[i,6]
    ATAC_peak_ATF1$V9 <- ATF1[i,7]
    ATAC_peak_ATF1_list = rbind(ATAC_peak_ATF1_list,ATAC_peak_ATF1)
  }
}
write.table(ATAC_peak_ATF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATF1_list.txt",sep="\t")

ATAC_peak_ATF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATF1_list.txt",sep="\t")

ATF1_features=list("ATF1" = ATAC_peak_ATF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ATF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZEB2
ZEB2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF225SIS_ZEB2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZEB2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF225SIS_ZEB2.txt",sep = "\t")

ZEB2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF225SIS_ZEB2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZEB2)){
  x=paste(ZEB2[i,1],ZEB2[i,2],ZEB2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZEB2_features = list(list$V1)

ATAC_peak_ZEB2_list = data.frame()
for(i in 1:nrow(ZEB2)){
  ATAC_peak_ZEB2<- ATAC_peak[ATAC_peak$x1 == ZEB2[i,1] & ((ZEB2[i,2] > ATAC_peak$start & ZEB2[i,2] < ATAC_peak$end) | (ZEB2[i,3] > ATAC_peak$start & ZEB2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZEB2) != 0){
    ATAC_peak_ZEB2$V5 <- ZEB2[i,1]
    ATAC_peak_ZEB2$V6 <- ZEB2[i,2]
    ATAC_peak_ZEB2$V7 <- ZEB2[i,3]
    ATAC_peak_ZEB2$V8 <- ZEB2[i,6]
    ATAC_peak_ZEB2$V9 <- ZEB2[i,7]
    ATAC_peak_ZEB2_list = rbind(ATAC_peak_ZEB2_list,ATAC_peak_ZEB2)
  }
}
write.table(ATAC_peak_ZEB2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZEB2_list.txt",sep="\t")

ATAC_peak_ZEB2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZEB2_list.txt",sep="\t")

ZEB2_features=list("ZEB2" = ATAC_peak_ZEB2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZEB2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NCOA1
NCOA1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF226JZE_NCOA1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NCOA1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF226JZE_NCOA1.txt",sep = "\t")

NCOA1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF226JZE_NCOA1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NCOA1)){
  x=paste(NCOA1[i,1],NCOA1[i,2],NCOA1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NCOA1_features = list(list$V1)

ATAC_peak_NCOA1_list = data.frame()
for(i in 1:nrow(NCOA1)){
  ATAC_peak_NCOA1<- ATAC_peak[ATAC_peak$x1 == NCOA1[i,1] & ((NCOA1[i,2] > ATAC_peak$start & NCOA1[i,2] < ATAC_peak$end) | (NCOA1[i,3] > ATAC_peak$start & NCOA1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NCOA1) != 0){
    ATAC_peak_NCOA1$V5 <- NCOA1[i,1]
    ATAC_peak_NCOA1$V6 <- NCOA1[i,2]
    ATAC_peak_NCOA1$V7 <- NCOA1[i,3]
    ATAC_peak_NCOA1$V8 <- NCOA1[i,6]
    ATAC_peak_NCOA1$V9 <- NCOA1[i,7]
    ATAC_peak_NCOA1_list = rbind(ATAC_peak_NCOA1_list,ATAC_peak_NCOA1)
  }
}
write.table(ATAC_peak_NCOA1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOA1_list.txt",sep="\t")

ATAC_peak_NCOA1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOA1_list.txt",sep="\t")

NCOA1_features=list("NCOA1" = ATAC_peak_NCOA1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NCOA1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##BRD9
BRD9 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF226VFC_BRD9.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BRD9,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF226VFC_BRD9.txt",sep = "\t")

BRD9= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF226VFC_BRD9.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BRD9)){
  x=paste(BRD9[i,1],BRD9[i,2],BRD9[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BRD9_features = list(list$V1)

ATAC_peak_BRD9_list = data.frame()
for(i in 1:nrow(BRD9)){
  ATAC_peak_BRD9<- ATAC_peak[ATAC_peak$x1 == BRD9[i,1] & ((BRD9[i,2] > ATAC_peak$start & BRD9[i,2] < ATAC_peak$end) | (BRD9[i,3] > ATAC_peak$start & BRD9[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BRD9) != 0){
    ATAC_peak_BRD9$V5 <- BRD9[i,1]
    ATAC_peak_BRD9$V6 <- BRD9[i,2]
    ATAC_peak_BRD9$V7 <- BRD9[i,3]
    ATAC_peak_BRD9$V8 <- BRD9[i,6]
    ATAC_peak_BRD9$V9 <- BRD9[i,7]
    ATAC_peak_BRD9_list = rbind(ATAC_peak_BRD9_list,ATAC_peak_BRD9)
  }
}
write.table(ATAC_peak_BRD9_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRD9_list.txt",sep="\t")

ATAC_peak_BRD9_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRD9_list.txt",sep="\t")

BRD9_features=list("BRD9" = ATAC_peak_BRD9_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BRD9_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CCNT2
CCNT2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232GPW_CCNT2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CCNT2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232GPW_CCNT2.txt",sep = "\t")

CCNT2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232GPW_CCNT2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CCNT2)){
  x=paste(CCNT2[i,1],CCNT2[i,2],CCNT2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CCNT2_features = list(list$V1)

ATAC_peak_CCNT2_list = data.frame()
for(i in 1:nrow(CCNT2)){
  ATAC_peak_CCNT2<- ATAC_peak[ATAC_peak$x1 == CCNT2[i,1] & ((CCNT2[i,2] > ATAC_peak$start & CCNT2[i,2] < ATAC_peak$end) | (CCNT2[i,3] > ATAC_peak$start & CCNT2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CCNT2) != 0){
    ATAC_peak_CCNT2$V5 <- CCNT2[i,1]
    ATAC_peak_CCNT2$V6 <- CCNT2[i,2]
    ATAC_peak_CCNT2$V7 <- CCNT2[i,3]
    ATAC_peak_CCNT2$V8 <- CCNT2[i,6]
    ATAC_peak_CCNT2$V9 <- CCNT2[i,7]
    ATAC_peak_CCNT2_list = rbind(ATAC_peak_CCNT2_list,ATAC_peak_CCNT2)
  }
}
write.table(ATAC_peak_CCNT2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CCNT2_list.txt",sep="\t")

ATAC_peak_CCNT2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CCNT2_list.txt",sep="\t")

CCNT2_features=list("CCNT2" = ATAC_peak_CCNT2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CCNT2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MED1
MED1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232WLA_MED1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MED1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232WLA_MED1.txt",sep = "\t")

MED1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232WLA_MED1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MED1)){
  x=paste(MED1[i,1],MED1[i,2],MED1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MED1_features = list(list$V1)

ATAC_peak_MED1_list = data.frame()
for(i in 1:nrow(MED1)){
  ATAC_peak_MED1<- ATAC_peak[ATAC_peak$x1 == MED1[i,1] & ((MED1[i,2] > ATAC_peak$start & MED1[i,2] < ATAC_peak$end) | (MED1[i,3] > ATAC_peak$start & MED1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MED1) != 0){
    ATAC_peak_MED1$V5 <- MED1[i,1]
    ATAC_peak_MED1$V6 <- MED1[i,2]
    ATAC_peak_MED1$V7 <- MED1[i,3]
    ATAC_peak_MED1$V8 <- MED1[i,6]
    ATAC_peak_MED1$V9 <- MED1[i,7]
    ATAC_peak_MED1_list = rbind(ATAC_peak_MED1_list,ATAC_peak_MED1)
  }
}
write.table(ATAC_peak_MED1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MED1_list.txt",sep="\t")

ATAC_peak_MED1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MED1_list.txt",sep="\t")

MED1_features=list("MED1" = ATAC_peak_MED1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MED1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SMAD4
SMAD4 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232ZHD_SMAD4.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMAD4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232ZHD_SMAD4.txt",sep = "\t")

SMAD4= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF232ZHD_SMAD4.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMAD4)){
  x=paste(SMAD4[i,1],SMAD4[i,2],SMAD4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMAD4_features = list(list$V1)

ATAC_peak_SMAD4_list = data.frame()
for(i in 1:nrow(SMAD4)){
  ATAC_peak_SMAD4<- ATAC_peak[ATAC_peak$x1 == SMAD4[i,1] & ((SMAD4[i,2] > ATAC_peak$start & SMAD4[i,2] < ATAC_peak$end) | (SMAD4[i,3] > ATAC_peak$start & SMAD4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMAD4) != 0){
    ATAC_peak_SMAD4$V5 <- SMAD4[i,1]
    ATAC_peak_SMAD4$V6 <- SMAD4[i,2]
    ATAC_peak_SMAD4$V7 <- SMAD4[i,3]
    ATAC_peak_SMAD4$V8 <- SMAD4[i,6]
    ATAC_peak_SMAD4$V9 <- SMAD4[i,7]
    ATAC_peak_SMAD4_list = rbind(ATAC_peak_SMAD4_list,ATAC_peak_SMAD4)
  }
}
write.table(ATAC_peak_SMAD4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMAD4_list.txt",sep="\t")

ATAC_peak_SMAD4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMAD4_list.txt",sep="\t")

SMAD4_features=list("SMAD4" = ATAC_peak_SMAD4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMAD4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##PATZ1
PATZ1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF241TKU_PATZ1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PATZ1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF241TKU_PATZ1.txt",sep = "\t")

PATZ1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF241TKU_PATZ1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PATZ1)){
  x=paste(PATZ1[i,1],PATZ1[i,2],PATZ1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PATZ1_features = list(list$V1)

ATAC_peak_PATZ1_list = data.frame()
for(i in 1:nrow(PATZ1)){
  ATAC_peak_PATZ1<- ATAC_peak[ATAC_peak$x1 == PATZ1[i,1] & ((PATZ1[i,2] > ATAC_peak$start & PATZ1[i,2] < ATAC_peak$end) | (PATZ1[i,3] > ATAC_peak$start & PATZ1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PATZ1) != 0){
    ATAC_peak_PATZ1$V5 <- PATZ1[i,1]
    ATAC_peak_PATZ1$V6 <- PATZ1[i,2]
    ATAC_peak_PATZ1$V7 <- PATZ1[i,3]
    ATAC_peak_PATZ1$V8 <- PATZ1[i,6]
    ATAC_peak_PATZ1$V9 <- PATZ1[i,7]
    ATAC_peak_PATZ1_list = rbind(ATAC_peak_PATZ1_list,ATAC_peak_PATZ1)
  }
}
write.table(ATAC_peak_PATZ1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PATZ1_list.txt",sep="\t")

ATAC_peak_PATZ1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PATZ1_list.txt",sep="\t")

PATZ1_features=list("PATZ1" = ATAC_peak_PATZ1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PATZ1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RAD51
RAD51 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF244WAB_RAD51_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RAD51,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF244WAB_RAD51_age.txt",sep = "\t")

RAD51 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF244WAB_RAD51_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RAD51)){
  x=paste(RAD51[i,1],RAD51[i,2],RAD51[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RAD51_features = list(list$V1)

ATAC_peak_RAD51_list = data.frame()
for(i in 1:nrow(RAD51)){
  ATAC_peak_RAD51<- ATAC_peak[ATAC_peak$x1 == RAD51[i,1] & ((RAD51[i,2] > ATAC_peak$start & RAD51[i,2] < ATAC_peak$end) | (RAD51[i,3] > ATAC_peak$start & RAD51[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RAD51) != 0){
    ATAC_peak_RAD51$V5 <- RAD51[i,1]
    ATAC_peak_RAD51$V6 <- RAD51[i,2]
    ATAC_peak_RAD51$V7 <- RAD51[i,3]
    ATAC_peak_RAD51$V8 <- RAD51[i,6]
    ATAC_peak_RAD51$V9 <- RAD51[i,7]
    ATAC_peak_RAD51_list = rbind(ATAC_peak_RAD51_list,ATAC_peak_RAD51)
  }
}
write.table(ATAC_peak_RAD51_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RAD51_list.txt",sep="\t")

ATAC_peak_RAD51_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RAD51_list.txt",sep="\t")

RAD51_features=list("RAD51" = ATAC_peak_RAD51_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RAD51_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ATF3
ATF3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF246USC_ATF3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ATF3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF246USC_ATF3.txt",sep = "\t")

ATF3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF246USC_ATF3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ATF3)){
  x=paste(ATF3[i,1],ATF3[i,2],ATF3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ATF3_features = list(list$V1)

ATAC_peak_ATF3_list = data.frame()
for(i in 1:nrow(ATF3)){
  ATAC_peak_ATF3<- ATAC_peak[ATAC_peak$x1 == ATF3[i,1] & ((ATF3[i,2] > ATAC_peak$start & ATF3[i,2] < ATAC_peak$end) | (ATF3[i,3] > ATAC_peak$start & ATF3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ATF3) != 0){
    ATAC_peak_ATF3$V5 <- ATF3[i,1]
    ATAC_peak_ATF3$V6 <- ATF3[i,2]
    ATAC_peak_ATF3$V7 <- ATF3[i,3]
    ATAC_peak_ATF3$V8 <- ATF3[i,6]
    ATAC_peak_ATF3$V9 <- ATF3[i,7]
    ATAC_peak_ATF3_list = rbind(ATAC_peak_ATF3_list,ATAC_peak_ATF3)
  }
}
write.table(ATAC_peak_ATF3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATF3_list.txt",sep="\t")

ATAC_peak_ATF3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATF3_list.txt",sep="\t")

ATF3_features=list("ATF3" = ATAC_peak_ATF3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ATF3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NR2F6
NR2F6 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF246VEZ_NR2F6.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NR2F6,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF246VEZ_NR2F6.txt",sep = "\t")

NR2F6 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF246VEZ_NR2F6.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NR2F6)){
  x=paste(NR2F6[i,1],NR2F6[i,2],NR2F6[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NR2F6_features = list(list$V1)

ATAC_peak_NR2F6_list = data.frame()
for(i in 1:nrow(NR2F6)){
  ATAC_peak_NR2F6<- ATAC_peak[ATAC_peak$x1 == NR2F6[i,1] & ((NR2F6[i,2] > ATAC_peak$start & NR2F6[i,2] < ATAC_peak$end) | (NR2F6[i,3] > ATAC_peak$start & NR2F6[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NR2F6) != 0){
    ATAC_peak_NR2F6$V5 <- NR2F6[i,1]
    ATAC_peak_NR2F6$V6 <- NR2F6[i,2]
    ATAC_peak_NR2F6$V7 <- NR2F6[i,3]
    ATAC_peak_NR2F6$V8 <- NR2F6[i,6]
    ATAC_peak_NR2F6$V9 <- NR2F6[i,7]
    ATAC_peak_NR2F6_list = rbind(ATAC_peak_NR2F6_list,ATAC_peak_NR2F6)
  }
}
write.table(ATAC_peak_NR2F6_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NR2F6_list.txt",sep="\t")

ATAC_peak_NR2F6_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NR2F6_list.txt",sep="\t")

NR2F6_features=list("NR2F6" = ATAC_peak_NR2F6_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NR2F6_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZZZ3
ZZZ3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF251IFL_ZZZ3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZZZ3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF251IFL_ZZZ3.txt",sep = "\t")

ZZZ3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF251IFL_ZZZ3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZZZ3)){
  x=paste(ZZZ3[i,1],ZZZ3[i,2],ZZZ3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZZZ3_features = list(list$V1)

ATAC_peak_ZZZ3_list = data.frame()
for(i in 1:nrow(ZZZ3)){
  ATAC_peak_ZZZ3<- ATAC_peak[ATAC_peak$x1 == ZZZ3[i,1] & ((ZZZ3[i,2] > ATAC_peak$start & ZZZ3[i,2] < ATAC_peak$end) | (ZZZ3[i,3] > ATAC_peak$start & ZZZ3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZZZ3) != 0){
    ATAC_peak_ZZZ3$V5 <- ZZZ3[i,1]
    ATAC_peak_ZZZ3$V6 <- ZZZ3[i,2]
    ATAC_peak_ZZZ3$V7 <- ZZZ3[i,3]
    ATAC_peak_ZZZ3$V8 <- ZZZ3[i,6]
    ATAC_peak_ZZZ3$V9 <- ZZZ3[i,7]
    ATAC_peak_ZZZ3_list = rbind(ATAC_peak_ZZZ3_list,ATAC_peak_ZZZ3)
  }
}
write.table(ATAC_peak_ZZZ3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZZZ3_list.txt",sep="\t")

ATAC_peak_ZZZ3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZZZ3_list.txt",sep="\t")

ZZZ3_features=list("ZZZ3" = ATAC_peak_ZZZ3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZZZ3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CBX2
CBX2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF258XBJ_CBX2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CBX2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF258XBJ_CBX2.txt",sep = "\t")

CBX2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF258XBJ_CBX2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CBX2)){
  x=paste(CBX2[i,1],CBX2[i,2],CBX2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CBX2_features = list(list$V1)

ATAC_peak_CBX2_list = data.frame()
for(i in 1:nrow(CBX2)){
  ATAC_peak_CBX2<- ATAC_peak[ATAC_peak$x1 == CBX2[i,1] & ((CBX2[i,2] > ATAC_peak$start & CBX2[i,2] < ATAC_peak$end) | (CBX2[i,3] > ATAC_peak$start & CBX2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CBX2) != 0){
    ATAC_peak_CBX2$V5 <- CBX2[i,1]
    ATAC_peak_CBX2$V6 <- CBX2[i,2]
    ATAC_peak_CBX2$V7 <- CBX2[i,3]
    ATAC_peak_CBX2$V8 <- CBX2[i,6]
    ATAC_peak_CBX2$V9 <- CBX2[i,7]
    ATAC_peak_CBX2_list = rbind(ATAC_peak_CBX2_list,ATAC_peak_CBX2)
  }
}
write.table(ATAC_peak_CBX2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX2_list.txt",sep="\t")

ATAC_peak_CBX2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX2_list.txt",sep="\t")

CBX2_features=list("CBX2" = ATAC_peak_CBX2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CBX2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##FOXM1
FOXM1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF262OXZ_FOXM1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(FOXM1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF262OXZ_FOXM1_age.txt",sep = "\t")

FOXM1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF262OXZ_FOXM1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(FOXM1)){
  x=paste(FOXM1[i,1],FOXM1[i,2],FOXM1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

FOXM1_features = list(list$V1)

ATAC_peak_FOXM1_list = data.frame()
for(i in 1:nrow(FOXM1)){
  ATAC_peak_FOXM1<- ATAC_peak[ATAC_peak$x1 == FOXM1[i,1] & ((FOXM1[i,2] > ATAC_peak$start & FOXM1[i,2] < ATAC_peak$end) | (FOXM1[i,3] > ATAC_peak$start & FOXM1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_FOXM1) != 0){
    ATAC_peak_FOXM1$V5 <- FOXM1[i,1]
    ATAC_peak_FOXM1$V6 <- FOXM1[i,2]
    ATAC_peak_FOXM1$V7 <- FOXM1[i,3]
    ATAC_peak_FOXM1$V8 <- FOXM1[i,6]
    ATAC_peak_FOXM1$V9 <- FOXM1[i,7]
    ATAC_peak_FOXM1_list = rbind(ATAC_peak_FOXM1_list,ATAC_peak_FOXM1)
  }
}
write.table(ATAC_peak_FOXM1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXM1_list.txt",sep="\t")

ATAC_peak_FOXM1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXM1_list.txt",sep="\t")

FOXM1_features=list("FOXM1" = ATAC_peak_FOXM1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = FOXM1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZBTB2
ZBTB2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF264HLZ_ZBTB2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZBTB2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF264HLZ_ZBTB2.txt",sep = "\t")

ZBTB2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF264HLZ_ZBTB2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZBTB2)){
  x=paste(ZBTB2[i,1],ZBTB2[i,2],ZBTB2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZBTB2_features = list(list$V1)

ATAC_peak_ZBTB2_list = data.frame()
for(i in 1:nrow(ZBTB2)){
  ATAC_peak_ZBTB2<- ATAC_peak[ATAC_peak$x1 == ZBTB2[i,1] & ((ZBTB2[i,2] > ATAC_peak$start & ZBTB2[i,2] < ATAC_peak$end) | (ZBTB2[i,3] > ATAC_peak$start & ZBTB2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZBTB2) != 0){
    ATAC_peak_ZBTB2$V5 <- ZBTB2[i,1]
    ATAC_peak_ZBTB2$V6 <- ZBTB2[i,2]
    ATAC_peak_ZBTB2$V7 <- ZBTB2[i,3]
    ATAC_peak_ZBTB2$V8 <- ZBTB2[i,6]
    ATAC_peak_ZBTB2$V9 <- ZBTB2[i,7]
    ATAC_peak_ZBTB2_list = rbind(ATAC_peak_ZBTB2_list,ATAC_peak_ZBTB2)
  }
}
write.table(ATAC_peak_ZBTB2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB2_list.txt",sep="\t")

ATAC_peak_ZBTB2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB2_list.txt",sep="\t")

ZBTB2_features=list("ZBTB2" = ATAC_peak_ZBTB2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZBTB2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##GATAD2B
GATAD2B = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF270LDV_GATAD2B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GATAD2B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF270LDV_GATAD2B.txt",sep = "\t")

GATAD2B = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF270LDV_GATAD2B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GATAD2B)){
  x=paste(GATAD2B[i,1],GATAD2B[i,2],GATAD2B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GATAD2B_features = list(list$V1)

ATAC_peak_GATAD2B_list = data.frame()
for(i in 1:nrow(GATAD2B)){
  ATAC_peak_GATAD2B<- ATAC_peak[ATAC_peak$x1 == GATAD2B[i,1] & ((GATAD2B[i,2] > ATAC_peak$start & GATAD2B[i,2] < ATAC_peak$end) | (GATAD2B[i,3] > ATAC_peak$start & GATAD2B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GATAD2B) != 0){
    ATAC_peak_GATAD2B$V5 <- GATAD2B[i,1]
    ATAC_peak_GATAD2B$V6 <- GATAD2B[i,2]
    ATAC_peak_GATAD2B$V7 <- GATAD2B[i,3]
    ATAC_peak_GATAD2B$V8 <- GATAD2B[i,6]
    ATAC_peak_GATAD2B$V9 <- GATAD2B[i,7]
    ATAC_peak_GATAD2B_list = rbind(ATAC_peak_GATAD2B_list,ATAC_peak_GATAD2B)
  }
}
write.table(ATAC_peak_GATAD2B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GATAD2B_list.txt",sep="\t")

ATAC_peak_GATAD2B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GATAD2B_list.txt",sep="\t")

GATAD2B_features=list("GATAD2B" = ATAC_peak_GATAD2B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GATAD2B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##MBD2
MBD2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF271DMO_MBD2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MBD2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF271DMO_MBD2.txt",sep = "\t")

MBD2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF271DMO_MBD2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MBD2)){
  x=paste(MBD2[i,1],MBD2[i,2],MBD2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MBD2_features = list(list$V1)

ATAC_peak_MBD2_list = data.frame()
for(i in 1:nrow(MBD2)){
  ATAC_peak_MBD2<- ATAC_peak[ATAC_peak$x1 == MBD2[i,1] & ((MBD2[i,2] > ATAC_peak$start & MBD2[i,2] < ATAC_peak$end) | (MBD2[i,3] > ATAC_peak$start & MBD2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MBD2) != 0){
    ATAC_peak_MBD2$V5 <- MBD2[i,1]
    ATAC_peak_MBD2$V6 <- MBD2[i,2]
    ATAC_peak_MBD2$V7 <- MBD2[i,3]
    ATAC_peak_MBD2$V8 <- MBD2[i,6]
    ATAC_peak_MBD2$V9 <- MBD2[i,7]
    ATAC_peak_MBD2_list = rbind(ATAC_peak_MBD2_list,ATAC_peak_MBD2)
  }
}
write.table(ATAC_peak_MBD2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MBD2_list.txt",sep="\t")

ATAC_peak_MBD2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MBD2_list.txt",sep="\t")

MBD2_features=list("MBD2" = ATAC_peak_MBD2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MBD2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##JUND
JUND = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF273KIA_JUND.bed.gz",extraCols=extraCols_narrowPeak)
write.table(JUND,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/NCFF273KIA_JUND.txt",sep = "\t")

JUND = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/NCFF273KIA_JUND.txt",sep = "\t")

list = vector()
for (i in  1:nrow(JUND)){
  x=paste(JUND[i,1],JUND[i,2],JUND[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

JUND_features = list(list$V1)

ATAC_peak_JUND_list = data.frame()
for(i in 1:nrow(JUND)){
  ATAC_peak_JUND<- ATAC_peak[ATAC_peak$x1 == JUND[i,1] & ((JUND[i,2] > ATAC_peak$start & JUND[i,2] < ATAC_peak$end) | (JUND[i,3] > ATAC_peak$start & JUND[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_JUND) != 0){
    ATAC_peak_JUND$V5 <- JUND[i,1]
    ATAC_peak_JUND$V6 <- JUND[i,2]
    ATAC_peak_JUND$V7 <- JUND[i,3]
    ATAC_peak_JUND$V8 <- JUND[i,6]
    ATAC_peak_JUND$V9 <- JUND[i,7]
    ATAC_peak_JUND_list = rbind(ATAC_peak_JUND_list,ATAC_peak_JUND)
  }
}
write.table(ATAC_peak_JUND_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_JUND_list.txt",sep="\t")

ATAC_peak_JUND_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_JUND_list.txt",sep="\t")

JUND_features=list("JUND" = ATAC_peak_JUND_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = JUND_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##POLR2A
POLR2A = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF286QTF_POLR2A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(POLR2A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF286QTF_POLR2A.txt",sep = "\t")

POLR2A = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF286QTF_POLR2A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(POLR2A)){
  x=paste(POLR2A[i,1],POLR2A[i,2],POLR2A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

POLR2A_features = list(list$V1)

ATAC_peak_POLR2A_list = data.frame()
for(i in 1:nrow(POLR2A)){
  ATAC_peak_POLR2A<- ATAC_peak[ATAC_peak$x1 == POLR2A[i,1] & ((POLR2A[i,2] > ATAC_peak$start & POLR2A[i,2] < ATAC_peak$end) | (POLR2A[i,3] > ATAC_peak$start & POLR2A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_POLR2A) != 0){
    ATAC_peak_POLR2A$V5 <- POLR2A[i,1]
    ATAC_peak_POLR2A$V6 <- POLR2A[i,2]
    ATAC_peak_POLR2A$V7 <- POLR2A[i,3]
    ATAC_peak_POLR2A$V8 <- POLR2A[i,6]
    ATAC_peak_POLR2A$V9 <- POLR2A[i,7]
    ATAC_peak_POLR2A_list = rbind(ATAC_peak_POLR2A_list,ATAC_peak_POLR2A)
  }
}
write.table(ATAC_peak_POLR2A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLR2A_list.txt",sep="\t")

ATAC_peak_POLR2A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLR2A_list.txt",sep="\t")

POLR2A_features=list("POLR2A" = ATAC_peak_POLR2A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = POLR2A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SMC3
SMC3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF289LLT_SMC3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMC3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF289LLT_SMC3.txt",sep = "\t")

SMC3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF289LLT_SMC3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMC3)){
  x=paste(SMC3[i,1],SMC3[i,2],SMC3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMC3_features = list(list$V1)

ATAC_peak_SMC3_list = data.frame()
for(i in 1:nrow(SMC3)){
  ATAC_peak_SMC3<- ATAC_peak[ATAC_peak$x1 == SMC3[i,1] & ((SMC3[i,2] > ATAC_peak$start & SMC3[i,2] < ATAC_peak$end) | (SMC3[i,3] > ATAC_peak$start & SMC3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMC3) != 0){
    ATAC_peak_SMC3$V5 <- SMC3[i,1]
    ATAC_peak_SMC3$V6 <- SMC3[i,2]
    ATAC_peak_SMC3$V7 <- SMC3[i,3]
    ATAC_peak_SMC3$V8 <- SMC3[i,6]
    ATAC_peak_SMC3$V9 <- SMC3[i,7]
    ATAC_peak_SMC3_list = rbind(ATAC_peak_SMC3_list,ATAC_peak_SMC3)
  }
}
write.table(ATAC_peak_SMC3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMC3_list.txt",sep="\t")

ATAC_peak_SMC3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMC3_list.txt",sep="\t")

SMC3_features=list("SMC3" = ATAC_peak_SMC3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMC3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##EP400
EP400 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF294MTS_EP400.bed.gz",extraCols=extraCols_narrowPeak)
write.table(EP400,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF294MTS_EP400.txt",sep = "\t")

EP400 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF294MTS_EP400.txt",sep = "\t")

list = vector()
for (i in  1:nrow(EP400)){
  x=paste(EP400[i,1],EP400[i,2],EP400[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

EP400_features = list(list$V1)

ATAC_peak_EP400_list = data.frame()
for(i in 1:nrow(EP400)){
  ATAC_peak_EP400<- ATAC_peak[ATAC_peak$x1 == EP400[i,1] & ((EP400[i,2] > ATAC_peak$start & EP400[i,2] < ATAC_peak$end) | (EP400[i,3] > ATAC_peak$start & EP400[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_EP400) != 0){
    ATAC_peak_EP400$V5 <- EP400[i,1]
    ATAC_peak_EP400$V6 <- EP400[i,2]
    ATAC_peak_EP400$V7 <- EP400[i,3]
    ATAC_peak_EP400$V8 <- EP400[i,6]
    ATAC_peak_EP400$V9 <- EP400[i,7]
    ATAC_peak_EP400_list = rbind(ATAC_peak_EP400_list,ATAC_peak_EP400)
  }
}
write.table(ATAC_peak_EP400_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EP400_list.txt",sep="\t")

ATAC_peak_EP400_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EP400_list.txt",sep="\t")

EP400_features=list("EP400" = ATAC_peak_EP400_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = EP400_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##LARP7
LARP7 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF296KCF_LARP7.bed.gz",extraCols=extraCols_narrowPeak)
write.table(LARP7,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF296KCF_LARP7.txt",sep = "\t")

LARP7 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF296KCF_LARP7.txt",sep = "\t")

list = vector()
for (i in  1:nrow(LARP7)){
  x=paste(LARP7[i,1],LARP7[i,2],LARP7[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

LARP7_features = list(list$V1)

ATAC_peak_LARP7_list = data.frame()
for(i in 1:nrow(LARP7)){
  ATAC_peak_LARP7<- ATAC_peak[ATAC_peak$x1 == LARP7[i,1] & ((LARP7[i,2] > ATAC_peak$start & LARP7[i,2] < ATAC_peak$end) | (LARP7[i,3] > ATAC_peak$start & LARP7[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_LARP7) != 0){
    ATAC_peak_LARP7$V5 <- LARP7[i,1]
    ATAC_peak_LARP7$V6 <- LARP7[i,2]
    ATAC_peak_LARP7$V7 <- LARP7[i,3]
    ATAC_peak_LARP7$V8 <- LARP7[i,6]
    ATAC_peak_LARP7$V9 <- LARP7[i,7]
    ATAC_peak_LARP7_list = rbind(ATAC_peak_LARP7_list,ATAC_peak_LARP7)
  }
}
write.table(ATAC_peak_LARP7_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_LARP7_list.txt",sep="\t")

ATAC_peak_LARP7_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_LARP7_list.txt",sep="\t")

LARP7_features=list("LARP7" = ATAC_peak_LARP7_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = LARP7_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##PHF20
PHF20 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF304KDK_PHF20.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PHF20,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF304KDK_PHF20.txt",sep = "\t")

PHF20 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF304KDK_PHF20.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PHF20)){
  x=paste(PHF20[i,1],PHF20[i,2],PHF20[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PHF20_features = list(list$V1)

ATAC_peak_PHF20_list = data.frame()
for(i in 1:nrow(PHF20)){
  ATAC_peak_PHF20<- ATAC_peak[ATAC_peak$x1 == PHF20[i,1] & ((PHF20[i,2] > ATAC_peak$start & PHF20[i,2] < ATAC_peak$end) | (PHF20[i,3] > ATAC_peak$start & PHF20[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PHF20) != 0){
    ATAC_peak_PHF20$V5 <- PHF20[i,1]
    ATAC_peak_PHF20$V6 <- PHF20[i,2]
    ATAC_peak_PHF20$V7 <- PHF20[i,3]
    ATAC_peak_PHF20$V8 <- PHF20[i,6]
    ATAC_peak_PHF20$V9 <- PHF20[i,7]
    ATAC_peak_PHF20_list = rbind(ATAC_peak_PHF20_list,ATAC_peak_PHF20)
  }
}
write.table(ATAC_peak_PHF20_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PHF20_list.txt",sep="\t")

ATAC_peak_PHF20_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PHF20_list.txt",sep="\t")

PHF20_features=list("PHF20" = ATAC_peak_PHF20_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PHF20_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZBTB5
ZBTB5 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF306SAP_ZBTB5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZBTB5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF306SAP_ZBTB5.txt",sep = "\t")

ZBTB5 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF306SAP_ZBTB5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZBTB5)){
  x=paste(ZBTB5[i,1],ZBTB5[i,2],ZBTB5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZBTB5_features = list(list$V1)

ATAC_peak_ZBTB5_list = data.frame()
for(i in 1:nrow(ZBTB5)){
  ATAC_peak_ZBTB5<- ATAC_peak[ATAC_peak$x1 == ZBTB5[i,1] & ((ZBTB5[i,2] > ATAC_peak$start & ZBTB5[i,2] < ATAC_peak$end) | (ZBTB5[i,3] > ATAC_peak$start & ZBTB5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZBTB5) != 0){
    ATAC_peak_ZBTB5$V5 <- ZBTB5[i,1]
    ATAC_peak_ZBTB5$V6 <- ZBTB5[i,2]
    ATAC_peak_ZBTB5$V7 <- ZBTB5[i,3]
    ATAC_peak_ZBTB5$V8 <- ZBTB5[i,6]
    ATAC_peak_ZBTB5$V9 <- ZBTB5[i,7]
    ATAC_peak_ZBTB5_list = rbind(ATAC_peak_ZBTB5_list,ATAC_peak_ZBTB5)
  }
}
write.table(ATAC_peak_ZBTB5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB5_list.txt",sep="\t")

ATAC_peak_ZBTB5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB5_list.txt",sep="\t")

ZBTB5_features=list("ZBTB5" = ATAC_peak_ZBTB5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZBTB5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HMGXB4
HMGXB4 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF312OVC_HMGXB4.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HMGXB4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF312OVC_HMGXB4.txt",sep = "\t")

HMGXB4 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF312OVC_HMGXB4.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HMGXB4)){
  x=paste(HMGXB4[i,1],HMGXB4[i,2],HMGXB4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HMGXB4_features = list(list$V1)

ATAC_peak_HMGXB4_list = data.frame()
for(i in 1:nrow(HMGXB4)){
  ATAC_peak_HMGXB4<- ATAC_peak[ATAC_peak$x1 == HMGXB4[i,1] & ((HMGXB4[i,2] > ATAC_peak$start & HMGXB4[i,2] < ATAC_peak$end) | (HMGXB4[i,3] > ATAC_peak$start & HMGXB4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HMGXB4) != 0){
    ATAC_peak_HMGXB4$V5 <- HMGXB4[i,1]
    ATAC_peak_HMGXB4$V6 <- HMGXB4[i,2]
    ATAC_peak_HMGXB4$V7 <- HMGXB4[i,3]
    ATAC_peak_HMGXB4$V8 <- HMGXB4[i,6]
    ATAC_peak_HMGXB4$V9 <- HMGXB4[i,7]
    ATAC_peak_HMGXB4_list = rbind(ATAC_peak_HMGXB4_list,ATAC_peak_HMGXB4)
  }
}
write.table(ATAC_peak_HMGXB4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HMGXB4_list.txt",sep="\t")

ATAC_peak_HMGXB4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HMGXB4_list.txt",sep="\t")

HMGXB4_features=list("HMGXB4" = ATAC_peak_HMGXB4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HMGXB4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CHD4
CHD4 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF319EST_CHD4.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CHD4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF319EST_CHD4.txt",sep = "\t")

CHD4 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF319EST_CHD4.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CHD4)){
  x=paste(CHD4[i,1],CHD4[i,2],CHD4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CHD4_features = list(list$V1)

ATAC_peak_CHD4_list = data.frame()
for(i in 1:nrow(CHD4)){
  ATAC_peak_CHD4<- ATAC_peak[ATAC_peak$x1 == CHD4[i,1] & ((CHD4[i,2] > ATAC_peak$start & CHD4[i,2] < ATAC_peak$end) | (CHD4[i,3] > ATAC_peak$start & CHD4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CHD4) != 0){
    ATAC_peak_CHD4$V5 <- CHD4[i,1]
    ATAC_peak_CHD4$V6 <- CHD4[i,2]
    ATAC_peak_CHD4$V7 <- CHD4[i,3]
    ATAC_peak_CHD4$V8 <- CHD4[i,6]
    ATAC_peak_CHD4$V9 <- CHD4[i,7]
    ATAC_peak_CHD4_list = rbind(ATAC_peak_CHD4_list,ATAC_peak_CHD4)
  }
}
write.table(ATAC_peak_CHD4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHD4_list.txt",sep="\t")

ATAC_peak_CHD4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHD4_list.txt",sep="\t")

CHD4_features=list("CHD4" = ATAC_peak_CHD4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CHD4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZBTB8A
ZBTB8A = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF335GXH_ZBTB8A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZBTB8A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF335GXH_ZBTB8A.txt",sep = "\t")

ZBTB8A = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF335GXH_ZBTB8A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZBTB8A)){
  x=paste(ZBTB8A[i,1],ZBTB8A[i,2],ZBTB8A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZBTB8A_features = list(list$V1)

ATAC_peak_ZBTB8A_list = data.frame()
for(i in 1:nrow(ZBTB8A)){
  ATAC_peak_ZBTB8A<- ATAC_peak[ATAC_peak$x1 == ZBTB8A[i,1] & ((ZBTB8A[i,2] > ATAC_peak$start & ZBTB8A[i,2] < ATAC_peak$end) | (ZBTB8A[i,3] > ATAC_peak$start & ZBTB8A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZBTB8A) != 0){
    ATAC_peak_ZBTB8A$V5 <- ZBTB8A[i,1]
    ATAC_peak_ZBTB8A$V6 <- ZBTB8A[i,2]
    ATAC_peak_ZBTB8A$V7 <- ZBTB8A[i,3]
    ATAC_peak_ZBTB8A$V8 <- ZBTB8A[i,6]
    ATAC_peak_ZBTB8A$V9 <- ZBTB8A[i,7]
    ATAC_peak_ZBTB8A_list = rbind(ATAC_peak_ZBTB8A_list,ATAC_peak_ZBTB8A)
  }
}
write.table(ATAC_peak_ZBTB8A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB8A_list.txt",sep="\t")

ATAC_peak_ZBTB8A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB8A_list.txt",sep="\t")

ZBTB8A_features=list("ZBTB8A" = ATAC_peak_ZBTB8A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZBTB8A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TAF15
TAF15 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF337IVP_TAF15.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TAF15,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF337IVP_TAF15.txt",sep = "\t")

TAF15 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF337IVP_TAF15.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TAF15)){
  x=paste(TAF15[i,1],TAF15[i,2],TAF15[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TAF15_features = list(list$V1)

ATAC_peak_TAF15_list = data.frame()
for(i in 1:nrow(TAF15)){
  ATAC_peak_TAF15<- ATAC_peak[ATAC_peak$x1 == TAF15[i,1] & ((TAF15[i,2] > ATAC_peak$start & TAF15[i,2] < ATAC_peak$end) | (TAF15[i,3] > ATAC_peak$start & TAF15[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TAF15) != 0){
    ATAC_peak_TAF15$V5 <- TAF15[i,1]
    ATAC_peak_TAF15$V6 <- TAF15[i,2]
    ATAC_peak_TAF15$V7 <- TAF15[i,3]
    ATAC_peak_TAF15$V8 <- TAF15[i,6]
    ATAC_peak_TAF15$V9 <- TAF15[i,7]
    ATAC_peak_TAF15_list = rbind(ATAC_peak_TAF15_list,ATAC_peak_TAF15)
  }
}
write.table(ATAC_peak_TAF15_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAF15_list.txt",sep="\t")

ATAC_peak_TAF15_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAF15_list.txt",sep="\t")

TAF15_features=list("TAF15" = ATAC_peak_TAF15_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TAF15_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HDAC2
HDAC2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF339LFS_HDAC2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HDAC2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF339LFS_HDAC2.txt",sep = "\t")

HDAC2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF339LFS_HDAC2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HDAC2)){
  x=paste(HDAC2[i,1],HDAC2[i,2],HDAC2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HDAC2_features = list(list$V1)

ATAC_peak_HDAC2_list = data.frame()
for(i in 1:nrow(HDAC2)){
  ATAC_peak_HDAC2<- ATAC_peak[ATAC_peak$x1 == HDAC2[i,1] & ((HDAC2[i,2] > ATAC_peak$start & HDAC2[i,2] < ATAC_peak$end) | (HDAC2[i,3] > ATAC_peak$start & HDAC2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HDAC2) != 0){
    ATAC_peak_HDAC2$V5 <- HDAC2[i,1]
    ATAC_peak_HDAC2$V6 <- HDAC2[i,2]
    ATAC_peak_HDAC2$V7 <- HDAC2[i,3]
    ATAC_peak_HDAC2$V8 <- HDAC2[i,6]
    ATAC_peak_HDAC2$V9 <- HDAC2[i,7]
    ATAC_peak_HDAC2_list = rbind(ATAC_peak_HDAC2_list,ATAC_peak_HDAC2)
  }
}
write.table(ATAC_peak_HDAC2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC2_list.txt",sep="\t")

ATAC_peak_HDAC2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC2_list.txt",sep="\t")

HDAC2_features=list("HDAC2" = ATAC_peak_HDAC2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HDAC2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##E2F5
E2F5 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF347EPZ_E2F5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(E2F5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF347EPZ_E2F5.txt",sep = "\t")

E2F5 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF347EPZ_E2F5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(E2F5)){
  x=paste(E2F5[i,1],E2F5[i,2],E2F5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

E2F5_features = list(list$V1)

ATAC_peak_E2F5_list = data.frame()
for(i in 1:nrow(E2F5)){
  ATAC_peak_E2F5<- ATAC_peak[ATAC_peak$x1 == E2F5[i,1] & ((E2F5[i,2] > ATAC_peak$start & E2F5[i,2] < ATAC_peak$end) | (E2F5[i,3] > ATAC_peak$start & E2F5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_E2F5) != 0){
    ATAC_peak_E2F5$V5 <- E2F5[i,1]
    ATAC_peak_E2F5$V6 <- E2F5[i,2]
    ATAC_peak_E2F5$V7 <- E2F5[i,3]
    ATAC_peak_E2F5$V8 <- E2F5[i,6]
    ATAC_peak_E2F5$V9 <- E2F5[i,7]
    ATAC_peak_E2F5_list = rbind(ATAC_peak_E2F5_list,ATAC_peak_E2F5)
  }
}
write.table(ATAC_peak_E2F5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_E2F5_list.txt",sep="\t")

ATAC_peak_E2F5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_E2F5_list.txt",sep="\t")

E2F5_features=list("E2F5" = ATAC_peak_E2F5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = E2F5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ATM
ATM = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF347WVQ_ATM_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ATM,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF347WVQ_ATM_age.txt",sep = "\t")

ATM = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF347WVQ_ATM_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ATM)){
  x=paste(ATM[i,1],ATM[i,2],ATM[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ATM_features = list(list$V1)

ATAC_peak_ATM_list = data.frame()
for(i in 1:nrow(ATM)){
  ATAC_peak_ATM<- ATAC_peak[ATAC_peak$x1 == ATM[i,1] & ((ATM[i,2] > ATAC_peak$start & ATM[i,2] < ATAC_peak$end) | (ATM[i,3] > ATAC_peak$start & ATM[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ATM) != 0){
    ATAC_peak_ATM$V5 <- ATM[i,1]
    ATAC_peak_ATM$V6 <- ATM[i,2]
    ATAC_peak_ATM$V7 <- ATM[i,3]
    ATAC_peak_ATM$V8 <- ATM[i,6]
    ATAC_peak_ATM$V9 <- ATM[i,7]
    ATAC_peak_ATM_list = rbind(ATAC_peak_ATM_list,ATAC_peak_ATM)
  }
}
write.table(ATAC_peak_ATM_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATM_list.txt",sep="\t")

ATAC_peak_ATM_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATM_list.txt",sep="\t")

ATM_features=list("ATM" = ATAC_peak_ATM_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ATM_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##KAT2B
KAT2B = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF349VSP_KAT2B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(KAT2B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF349VSP_KAT2B.txt",sep = "\t")

KAT2B = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF349VSP_KAT2B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(KAT2B)){
  x=paste(KAT2B[i,1],KAT2B[i,2],KAT2B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

KAT2B_features = list(list$V1)

ATAC_peak_KAT2B_list = data.frame()
for(i in 1:nrow(KAT2B)){
  ATAC_peak_KAT2B<- ATAC_peak[ATAC_peak$x1 == KAT2B[i,1] & ((KAT2B[i,2] > ATAC_peak$start & KAT2B[i,2] < ATAC_peak$end) | (KAT2B[i,3] > ATAC_peak$start & KAT2B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_KAT2B) != 0){
    ATAC_peak_KAT2B$V5 <- KAT2B[i,1]
    ATAC_peak_KAT2B$V6 <- KAT2B[i,2]
    ATAC_peak_KAT2B$V7 <- KAT2B[i,3]
    ATAC_peak_KAT2B$V8 <- KAT2B[i,6]
    ATAC_peak_KAT2B$V9 <- KAT2B[i,7]
    ATAC_peak_KAT2B_list = rbind(ATAC_peak_KAT2B_list,ATAC_peak_KAT2B)
  }
}
write.table(ATAC_peak_KAT2B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KAT2B_list.txt",sep="\t")

ATAC_peak_KAT2B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KAT2B_list.txt",sep="\t")

KAT2B_features=list("KAT2B" = ATAC_peak_KAT2B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = KAT2B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##YBX3
YBX3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF351GRZ_YBX3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(YBX3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF351GRZ_YBX3.txt",sep = "\t")

YBX3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF351GRZ_YBX3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(YBX3)){
  x=paste(YBX3[i,1],YBX3[i,2],YBX3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

YBX3_features = list(list$V1)

ATAC_peak_YBX3_list = data.frame()
for(i in 1:nrow(YBX3)){
  ATAC_peak_YBX3<- ATAC_peak[ATAC_peak$x1 == YBX3[i,1] & ((YBX3[i,2] > ATAC_peak$start & YBX3[i,2] < ATAC_peak$end) | (YBX3[i,3] > ATAC_peak$start & YBX3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_YBX3) != 0){
    ATAC_peak_YBX3$V5 <- YBX3[i,1]
    ATAC_peak_YBX3$V6 <- YBX3[i,2]
    ATAC_peak_YBX3$V7 <- YBX3[i,3]
    ATAC_peak_YBX3$V8 <- YBX3[i,6]
    ATAC_peak_YBX3$V9 <- YBX3[i,7]
    ATAC_peak_YBX3_list = rbind(ATAC_peak_YBX3_list,ATAC_peak_YBX3)
  }
}
write.table(ATAC_peak_YBX3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_YBX3_list.txt",sep="\t")

ATAC_peak_YBX3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_YBX3_list.txt",sep="\t")

YBX3_features=list("YBX3" = ATAC_peak_YBX3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = YBX3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MTA1
MTA1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF352BVR_MTA1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MTA1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF352BVR_MTA1.txt",sep = "\t")

MTA1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF352BVR_MTA1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MTA1)){
  x=paste(MTA1[i,1],MTA1[i,2],MTA1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MTA1_features = list(list$V1)

ATAC_peak_MTA1_list = data.frame()
for(i in 1:nrow(MTA1)){
  ATAC_peak_MTA1<- ATAC_peak[ATAC_peak$x1 == MTA1[i,1] & ((MTA1[i,2] > ATAC_peak$start & MTA1[i,2] < ATAC_peak$end) | (MTA1[i,3] > ATAC_peak$start & MTA1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MTA1) != 0){
    ATAC_peak_MTA1$V5 <- MTA1[i,1]
    ATAC_peak_MTA1$V6 <- MTA1[i,2]
    ATAC_peak_MTA1$V7 <- MTA1[i,3]
    ATAC_peak_MTA1$V8 <- MTA1[i,6]
    ATAC_peak_MTA1$V9 <- MTA1[i,7]
    ATAC_peak_MTA1_list = rbind(ATAC_peak_MTA1_list,ATAC_peak_MTA1)
  }
}
write.table(ATAC_peak_MTA1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MTA1_list.txt",sep="\t")

ATAC_peak_MTA1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MTA1_list.txt",sep="\t")

MTA1_features=list("MTA1" = ATAC_peak_MTA1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MTA1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SAFB2
SAFB2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF352KMZ_SAFB2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SAFB2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF352KMZ_SAFB2.txt",sep = "\t")

SAFB2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF352KMZ_SAFB2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SAFB2)){
  x=paste(SAFB2[i,1],SAFB2[i,2],SAFB2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SAFB2_features = list(list$V1)

ATAC_peak_SAFB2_list = data.frame()
for(i in 1:nrow(SAFB2)){
  ATAC_peak_SAFB2<- ATAC_peak[ATAC_peak$x1 == SAFB2[i,1] & ((SAFB2[i,2] > ATAC_peak$start & SAFB2[i,2] < ATAC_peak$end) | (SAFB2[i,3] > ATAC_peak$start & SAFB2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SAFB2) != 0){
    ATAC_peak_SAFB2$V5 <- SAFB2[i,1]
    ATAC_peak_SAFB2$V6 <- SAFB2[i,2]
    ATAC_peak_SAFB2$V7 <- SAFB2[i,3]
    ATAC_peak_SAFB2$V8 <- SAFB2[i,6]
    ATAC_peak_SAFB2$V9 <- SAFB2[i,7]
    ATAC_peak_SAFB2_list = rbind(ATAC_peak_SAFB2_list,ATAC_peak_SAFB2)
  }
}
write.table(ATAC_peak_SAFB2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SAFB2_list.txt",sep="\t")

ATAC_peak_SAFB2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SAFB2_list.txt",sep="\t")

SAFB2_features=list("SAFB2" = ATAC_peak_SAFB2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SAFB2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SMARCB1
SMARCB1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF359MKF_SMARCB1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMARCB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF359MKF_SMARCB1.txt",sep = "\t")

SMARCB1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF359MKF_SMARCB1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMARCB1)){
  x=paste(SMARCB1[i,1],SMARCB1[i,2],SMARCB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMARCB1_features = list(list$V1)

ATAC_peak_SMARCB1_list = data.frame()
for(i in 1:nrow(SMARCB1)){
  ATAC_peak_SMARCB1<- ATAC_peak[ATAC_peak$x1 == SMARCB1[i,1] & ((SMARCB1[i,2] > ATAC_peak$start & SMARCB1[i,2] < ATAC_peak$end) | (SMARCB1[i,3] > ATAC_peak$start & SMARCB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMARCB1) != 0){
    ATAC_peak_SMARCB1$V5 <- SMARCB1[i,1]
    ATAC_peak_SMARCB1$V6 <- SMARCB1[i,2]
    ATAC_peak_SMARCB1$V7 <- SMARCB1[i,3]
    ATAC_peak_SMARCB1$V8 <- SMARCB1[i,6]
    ATAC_peak_SMARCB1$V9 <- SMARCB1[i,7]
    ATAC_peak_SMARCB1_list = rbind(ATAC_peak_SMARCB1_list,ATAC_peak_SMARCB1)
  }
}
write.table(ATAC_peak_SMARCB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMARCB1_list.txt",sep="\t")

ATAC_peak_SMARCB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMARCB1_list.txt",sep="\t")

SMARCB1_features=list("SMARCB1" = ATAC_peak_SMARCB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMARCB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RNF2
RNF2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF359QAT_RNF2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RNF2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF359QAT_RNF2.txt",sep = "\t")

RNF2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF359QAT_RNF2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RNF2)){
  x=paste(RNF2[i,1],RNF2[i,2],RNF2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RNF2_features = list(list$V1)

ATAC_peak_RNF2_list = data.frame()
for(i in 1:nrow(RNF2)){
  ATAC_peak_RNF2<- ATAC_peak[ATAC_peak$x1 == RNF2[i,1] & ((RNF2[i,2] > ATAC_peak$start & RNF2[i,2] < ATAC_peak$end) | (RNF2[i,3] > ATAC_peak$start & RNF2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RNF2) != 0){
    ATAC_peak_RNF2$V5 <- RNF2[i,1]
    ATAC_peak_RNF2$V6 <- RNF2[i,2]
    ATAC_peak_RNF2$V7 <- RNF2[i,3]
    ATAC_peak_RNF2$V8 <- RNF2[i,6]
    ATAC_peak_RNF2$V9 <- RNF2[i,7]
    ATAC_peak_RNF2_list = rbind(ATAC_peak_RNF2_list,ATAC_peak_RNF2)
  }
}
write.table(ATAC_peak_RNF2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RNF2_list.txt",sep="\t")

ATAC_peak_RNF2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RNF2_list.txt",sep="\t")

RNF2_features=list("RNF2" = ATAC_peak_RNF2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RNF2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SMAD2
SMAD2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF360SUP_SMAD2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMAD2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF360SUP_SMAD2.txt",sep = "\t")

SMAD2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF360SUP_SMAD2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMAD2)){
  x=paste(SMAD2[i,1],SMAD2[i,2],SMAD2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMAD2_features = list(list$V1)

ATAC_peak_SMAD2_list = data.frame()
for(i in 1:nrow(SMAD2)){
  ATAC_peak_SMAD2<- ATAC_peak[ATAC_peak$x1 == SMAD2[i,1] & ((SMAD2[i,2] > ATAC_peak$start & SMAD2[i,2] < ATAC_peak$end) | (SMAD2[i,3] > ATAC_peak$start & SMAD2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMAD2) != 0){
    ATAC_peak_SMAD2$V5 <- SMAD2[i,1]
    ATAC_peak_SMAD2$V6 <- SMAD2[i,2]
    ATAC_peak_SMAD2$V7 <- SMAD2[i,3]
    ATAC_peak_SMAD2$V8 <- SMAD2[i,6]
    ATAC_peak_SMAD2$V9 <- SMAD2[i,7]
    ATAC_peak_SMAD2_list = rbind(ATAC_peak_SMAD2_list,ATAC_peak_SMAD2)
  }
}
write.table(ATAC_peak_SMAD2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMAD2_list.txt",sep="\t")

ATAC_peak_SMAD2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMAD2_list.txt",sep="\t")

SMAD2_features=list("SMAD2" = ATAC_peak_SMAD2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMAD2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZBTB1
ZBTB1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF364JZV_ZBTB1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZBTB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF364JZV_ZBTB1.txt",sep = "\t")

ZBTB1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF364JZV_ZBTB1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZBTB1)){
  x=paste(ZBTB1[i,1],ZBTB1[i,2],ZBTB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZBTB1_features = list(list$V1)

ATAC_peak_ZBTB1_list = data.frame()
for(i in 1:nrow(ZBTB1)){
  ATAC_peak_ZBTB1<- ATAC_peak[ATAC_peak$x1 == ZBTB1[i,1] & ((ZBTB1[i,2] > ATAC_peak$start & ZBTB1[i,2] < ATAC_peak$end) | (ZBTB1[i,3] > ATAC_peak$start & ZBTB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZBTB1) != 0){
    ATAC_peak_ZBTB1$V5 <- ZBTB1[i,1]
    ATAC_peak_ZBTB1$V6 <- ZBTB1[i,2]
    ATAC_peak_ZBTB1$V7 <- ZBTB1[i,3]
    ATAC_peak_ZBTB1$V8 <- ZBTB1[i,6]
    ATAC_peak_ZBTB1$V9 <- ZBTB1[i,7]
    ATAC_peak_ZBTB1_list = rbind(ATAC_peak_ZBTB1_list,ATAC_peak_ZBTB1)
  }
}
write.table(ATAC_peak_ZBTB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB1_list.txt",sep="\t")

ATAC_peak_ZBTB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB1_list.txt",sep="\t")

ZBTB1_features=list("ZBTB1" = ATAC_peak_ZBTB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZBTB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SMARCA5
SMARCA5 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF365CXF_SMARCA5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMARCA5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/EENCFF365CXF_SMARCA5.txt",sep = "\t")

SMARCA5 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/EENCFF365CXF_SMARCA5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMARCA5)){
  x=paste(SMARCA5[i,1],SMARCA5[i,2],SMARCA5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMARCA5_features = list(list$V1)

ATAC_peak_SMARCA5_list = data.frame()
for(i in 1:nrow(SMARCA5)){
  ATAC_peak_SMARCA5<- ATAC_peak[ATAC_peak$x1 == SMARCA5[i,1] & ((SMARCA5[i,2] > ATAC_peak$start & SMARCA5[i,2] < ATAC_peak$end) | (SMARCA5[i,3] > ATAC_peak$start & SMARCA5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMARCA5) != 0){
    ATAC_peak_SMARCA5$V5 <- SMARCA5[i,1]
    ATAC_peak_SMARCA5$V6 <- SMARCA5[i,2]
    ATAC_peak_SMARCA5$V7 <- SMARCA5[i,3]
    ATAC_peak_SMARCA5$V8 <- SMARCA5[i,6]
    ATAC_peak_SMARCA5$V9 <- SMARCA5[i,7]
    ATAC_peak_SMARCA5_list = rbind(ATAC_peak_SMARCA5_list,ATAC_peak_SMARCA5)
  }
}
write.table(ATAC_peak_SMARCA5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMARCA5_list.txt",sep="\t")

ATAC_peak_SMARCA5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMARCA5_list.txt",sep="\t")

SMARCA5_features=list("SMARCA5" = ATAC_peak_SMARCA5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMARCA5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##BRF1
BRF1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF368OLS_BRF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BRF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF368OLS_BRF1.txt",sep = "\t")

BRF1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF368OLS_BRF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BRF1)){
  x=paste(BRF1[i,1],BRF1[i,2],BRF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BRF1_features = list(list$V1)

ATAC_peak_BRF1_list = data.frame()
for(i in 1:nrow(BRF1)){
  ATAC_peak_BRF1<- ATAC_peak[ATAC_peak$x1 == BRF1[i,1] & ((BRF1[i,2] > ATAC_peak$start & BRF1[i,2] < ATAC_peak$end) | (BRF1[i,3] > ATAC_peak$start & BRF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BRF1) != 0){
    ATAC_peak_BRF1$V5 <- BRF1[i,1]
    ATAC_peak_BRF1$V6 <- BRF1[i,2]
    ATAC_peak_BRF1$V7 <- BRF1[i,3]
    ATAC_peak_BRF1$V8 <- BRF1[i,6]
    ATAC_peak_BRF1$V9 <- BRF1[i,7]
    ATAC_peak_BRF1_list = rbind(ATAC_peak_BRF1_list,ATAC_peak_BRF1)
  }
}
write.table(ATAC_peak_BRF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRF1_list.txt",sep="\t")

ATAC_peak_BRF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRF1_list.txt",sep="\t")

BRF1_features=list("BRF1" = ATAC_peak_BRF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BRF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##FOXJ3
FOXJ3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF370RQD_FOXJ3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(FOXJ3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF370RQD_FOXJ3.txt",sep = "\t")

FOXJ3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF370RQD_FOXJ3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(FOXJ3)){
  x=paste(FOXJ3[i,1],FOXJ3[i,2],FOXJ3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

FOXJ3_features = list(list$V1)

ATAC_peak_FOXJ3_list = data.frame()
for(i in 1:nrow(FOXJ3)){
  ATAC_peak_FOXJ3<- ATAC_peak[ATAC_peak$x1 == FOXJ3[i,1] & ((FOXJ3[i,2] > ATAC_peak$start & FOXJ3[i,2] < ATAC_peak$end) | (FOXJ3[i,3] > ATAC_peak$start & FOXJ3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_FOXJ3) != 0){
    ATAC_peak_FOXJ3$V5 <- FOXJ3[i,1]
    ATAC_peak_FOXJ3$V6 <- FOXJ3[i,2]
    ATAC_peak_FOXJ3$V7 <- FOXJ3[i,3]
    ATAC_peak_FOXJ3$V8 <- FOXJ3[i,6]
    ATAC_peak_FOXJ3$V9 <- FOXJ3[i,7]
    ATAC_peak_FOXJ3_list = rbind(ATAC_peak_FOXJ3_list,ATAC_peak_FOXJ3)
  }
}
write.table(ATAC_peak_FOXJ3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXJ3_list.txt",sep="\t")

ATAC_peak_FOXJ3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXJ3_list.txt",sep="\t")

FOXJ3_features=list("FOXJ3" = ATAC_peak_FOXJ3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = FOXJ3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##EHMT2
EHMT2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF373EMV_EHMT2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(EHMT2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF373EMV_EHMT2.txt",sep = "\t")

EHMT2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF373EMV_EHMT2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(EHMT2)){
  x=paste(EHMT2[i,1],EHMT2[i,2],EHMT2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

EHMT2_features = list(list$V1)

ATAC_peak_EHMT2_list = data.frame()
for(i in 1:nrow(EHMT2)){
  ATAC_peak_EHMT2<- ATAC_peak[ATAC_peak$x1 == EHMT2[i,1] & ((EHMT2[i,2] > ATAC_peak$start & EHMT2[i,2] < ATAC_peak$end) | (EHMT2[i,3] > ATAC_peak$start & EHMT2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_EHMT2) != 0){
    ATAC_peak_EHMT2$V5 <- EHMT2[i,1]
    ATAC_peak_EHMT2$V6 <- EHMT2[i,2]
    ATAC_peak_EHMT2$V7 <- EHMT2[i,3]
    ATAC_peak_EHMT2$V8 <- EHMT2[i,6]
    ATAC_peak_EHMT2$V9 <- EHMT2[i,7]
    ATAC_peak_EHMT2_list = rbind(ATAC_peak_EHMT2_list,ATAC_peak_EHMT2)
  }
}
write.table(ATAC_peak_EHMT2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EHMT2_list.txt",sep="\t")

ATAC_peak_EHMT2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EHMT2_list.txt",sep="\t")

EHMT2_features=list("EHMT2" = ATAC_peak_EHMT2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = EHMT2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##FOS
FOS = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF378UUI_FOS_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(FOS,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF378UUI_FOS_age.txt",sep = "\t")

FOS = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF378UUI_FOS_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(FOS)){
  x=paste(FOS[i,1],FOS[i,2],FOS[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

FOS_features = list(list$V1)

ATAC_peak_FOS_list = data.frame()
for(i in 1:nrow(FOS)){
  ATAC_peak_FOS<- ATAC_peak[ATAC_peak$x1 == FOS[i,1] & ((FOS[i,2] > ATAC_peak$start & FOS[i,2] < ATAC_peak$end) | (FOS[i,3] > ATAC_peak$start & FOS[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_FOS) != 0){
    ATAC_peak_FOS$V5 <- FOS[i,1]
    ATAC_peak_FOS$V6 <- FOS[i,2]
    ATAC_peak_FOS$V7 <- FOS[i,3]
    ATAC_peak_FOS$V8 <- FOS[i,6]
    ATAC_peak_FOS$V9 <- FOS[i,7]
    ATAC_peak_FOS_list = rbind(ATAC_peak_FOS_list,ATAC_peak_FOS)
  }
}
write.table(ATAC_peak_FOS_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOS_list.txt",sep="\t")

ATAC_peak_FOS_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOS_list.txt",sep="\t")

FOS_features=list("FOS" = ATAC_peak_FOS_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = FOS_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TRIM28
TRIM28 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF382UDO_TRIM28.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TRIM28,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF382UDO_TRIM28.txt",sep = "\t")

TRIM28 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF382UDO_TRIM28.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TRIM28)){
  x=paste(TRIM28[i,1],TRIM28[i,2],TRIM28[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TRIM28_features = list(list$V1)

ATAC_peak_TRIM28_list = data.frame()
for(i in 1:nrow(TRIM28)){
  ATAC_peak_TRIM28<- ATAC_peak[ATAC_peak$x1 == TRIM28[i,1] & ((TRIM28[i,2] > ATAC_peak$start & TRIM28[i,2] < ATAC_peak$end) | (TRIM28[i,3] > ATAC_peak$start & TRIM28[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TRIM28) != 0){
    ATAC_peak_TRIM28$V5 <- TRIM28[i,1]
    ATAC_peak_TRIM28$V6 <- TRIM28[i,2]
    ATAC_peak_TRIM28$V7 <- TRIM28[i,3]
    ATAC_peak_TRIM28$V8 <- TRIM28[i,6]
    ATAC_peak_TRIM28$V9 <- TRIM28[i,7]
    ATAC_peak_TRIM28_list = rbind(ATAC_peak_TRIM28_list,ATAC_peak_TRIM28)
  }
}
write.table(ATAC_peak_TRIM28_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TRIM28_list.txt",sep="\t")

ATAC_peak_TRIM28_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TRIM28_list.txt",sep="\t")

TRIM28_features=list("TRIM28" = ATAC_peak_TRIM28_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TRIM28_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MCM3
MCM3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF384CRI_MCM3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MCM3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF384CRI_MCM3.txt",sep = "\t")

MCM3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF384CRI_MCM3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MCM3)){
  x=paste(MCM3[i,1],MCM3[i,2],MCM3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MCM3_features = list(list$V1)

ATAC_peak_MCM3_list = data.frame()
for(i in 1:nrow(MCM3)){
  ATAC_peak_MCM3<- ATAC_peak[ATAC_peak$x1 == MCM3[i,1] & ((MCM3[i,2] > ATAC_peak$start & MCM3[i,2] < ATAC_peak$end) | (MCM3[i,3] > ATAC_peak$start & MCM3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MCM3) != 0){
    ATAC_peak_MCM3$V5 <- MCM3[i,1]
    ATAC_peak_MCM3$V6 <- MCM3[i,2]
    ATAC_peak_MCM3$V7 <- MCM3[i,3]
    ATAC_peak_MCM3$V8 <- MCM3[i,6]
    ATAC_peak_MCM3$V9 <- MCM3[i,7]
    ATAC_peak_MCM3_list = rbind(ATAC_peak_MCM3_list,ATAC_peak_MCM3)
  }
}
write.table(ATAC_peak_MCM3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MCM3_list.txt",sep="\t")

ATAC_peak_MCM3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MCM3_list.txt",sep="\t")

MCM3_features=list("MCM3" = ATAC_peak_MCM3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MCM3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZNF318
ZNF318 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF384XCG_ZNF318.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF318,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF384XCG_ZNF318.txt",sep = "\t")

ZNF318 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF384XCG_ZNF318.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF318)){
  x=paste(ZNF318[i,1],ZNF318[i,2],ZNF318[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF318_features = list(list$V1)

ATAC_peak_ZNF318_list = data.frame()
for(i in 1:nrow(ZNF318)){
  ATAC_peak_ZNF318<- ATAC_peak[ATAC_peak$x1 == ZNF318[i,1] & ((ZNF318[i,2] > ATAC_peak$start & ZNF318[i,2] < ATAC_peak$end) | (ZNF318[i,3] > ATAC_peak$start & ZNF318[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF318) != 0){
    ATAC_peak_ZNF318$V5 <- ZNF318[i,1]
    ATAC_peak_ZNF318$V6 <- ZNF318[i,2]
    ATAC_peak_ZNF318$V7 <- ZNF318[i,3]
    ATAC_peak_ZNF318$V8 <- ZNF318[i,6]
    ATAC_peak_ZNF318$V9 <- ZNF318[i,7]
    ATAC_peak_ZNF318_list = rbind(ATAC_peak_ZNF318_list,ATAC_peak_ZNF318)
  }
}
write.table(ATAC_peak_ZNF318_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF318_list.txt",sep="\t")

ATAC_peak_ZNF318_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF318_list.txt",sep="\t")

ZNF318_features=list("ZNF318" = ATAC_peak_ZNF318_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF318_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SETDB1
SETDB1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF397SBJ_SETDB1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SETDB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF397SBJ_SETDB1.txt",sep = "\t")

SETDB1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF397SBJ_SETDB1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SETDB1)){
  x=paste(SETDB1[i,1],SETDB1[i,2],SETDB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SETDB1_features = list(list$V1)

ATAC_peak_SETDB1_list = data.frame()
for(i in 1:nrow(SETDB1)){
  ATAC_peak_SETDB1<- ATAC_peak[ATAC_peak$x1 == SETDB1[i,1] & ((SETDB1[i,2] > ATAC_peak$start & SETDB1[i,2] < ATAC_peak$end) | (SETDB1[i,3] > ATAC_peak$start & SETDB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SETDB1) != 0){
    ATAC_peak_SETDB1$V5 <- SETDB1[i,1]
    ATAC_peak_SETDB1$V6 <- SETDB1[i,2]
    ATAC_peak_SETDB1$V7 <- SETDB1[i,3]
    ATAC_peak_SETDB1$V8 <- SETDB1[i,6]
    ATAC_peak_SETDB1$V9 <- SETDB1[i,7]
    ATAC_peak_SETDB1_list = rbind(ATAC_peak_SETDB1_list,ATAC_peak_SETDB1)
  }
}
write.table(ATAC_peak_SETDB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SETDB1_list.txt",sep="\t")

ATAC_peak_SETDB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SETDB1_list.txt",sep="\t")

SETDB1_features=list("SETDB1" = ATAC_peak_SETDB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SETDB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SIX5
SIX5 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF397SOD_SIX5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SIX5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF397SOD_SIX5.txt",sep = "\t")

SIX5 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF397SOD_SIX5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SIX5)){
  x=paste(SIX5[i,1],SIX5[i,2],SIX5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SIX5_features = list(list$V1)

ATAC_peak_SIX5_list = data.frame()
for(i in 1:nrow(SIX5)){
  ATAC_peak_SIX5<- ATAC_peak[ATAC_peak$x1 == SIX5[i,1] & ((SIX5[i,2] > ATAC_peak$start & SIX5[i,2] < ATAC_peak$end) | (SIX5[i,3] > ATAC_peak$start & SIX5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SIX5) != 0){
    ATAC_peak_SIX5$V5 <- SIX5[i,1]
    ATAC_peak_SIX5$V6 <- SIX5[i,2]
    ATAC_peak_SIX5$V7 <- SIX5[i,3]
    ATAC_peak_SIX5$V8 <- SIX5[i,6]
    ATAC_peak_SIX5$V9 <- SIX5[i,7]
    ATAC_peak_SIX5_list = rbind(ATAC_peak_SIX5_list,ATAC_peak_SIX5)
  }
}
write.table(ATAC_peak_SIX5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIX5_list.txt",sep="\t")

ATAC_peak_SIX5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIX5_list.txt",sep="\t")

SIX5_features=list("SIX5" = ATAC_peak_SIX5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SIX5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TCF7
TCF7 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF401GJM_TCF7.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TCF7,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF401GJM_TCF7.txt",sep = "\t")

TCF7 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF401GJM_TCF7.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TCF7)){
  x=paste(TCF7[i,1],TCF7[i,2],TCF7[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TCF7_features = list(list$V1)

ATAC_peak_TCF7_list = data.frame()
for(i in 1:nrow(TCF7)){
  ATAC_peak_TCF7<- ATAC_peak[ATAC_peak$x1 == TCF7[i,1] & ((TCF7[i,2] > ATAC_peak$start & TCF7[i,2] < ATAC_peak$end) | (TCF7[i,3] > ATAC_peak$start & TCF7[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TCF7) != 0){
    ATAC_peak_TCF7$V5 <- TCF7[i,1]
    ATAC_peak_TCF7$V6 <- TCF7[i,2]
    ATAC_peak_TCF7$V7 <- TCF7[i,3]
    ATAC_peak_TCF7$V8 <- TCF7[i,6]
    ATAC_peak_TCF7$V9 <- TCF7[i,7]
    ATAC_peak_TCF7_list = rbind(ATAC_peak_TCF7_list,ATAC_peak_TCF7)
  }
}
write.table(ATAC_peak_TCF7_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TCF7_list.txt",sep="\t")

ATAC_peak_TCF7_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TCF7_list.txt",sep="\t")

TCF7_features=list("TCF7" = ATAC_peak_TCF7_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TCF7_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##UBTF
UBTF = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF401MXT_UBTF.bed.gz",extraCols=extraCols_narrowPeak)
write.table(UBTF,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF401MXT_UBTF.txt",sep = "\t")

UBTF = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF401MXT_UBTF.txt",sep = "\t")

list = vector()
for (i in  1:nrow(UBTF)){
  x=paste(UBTF[i,1],UBTF[i,2],UBTF[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

UBTF_features = list(list$V1)

ATAC_peak_UBTF_list = data.frame()
for(i in 1:nrow(UBTF)){
  ATAC_peak_UBTF<- ATAC_peak[ATAC_peak$x1 == UBTF[i,1] & ((UBTF[i,2] > ATAC_peak$start & UBTF[i,2] < ATAC_peak$end) | (UBTF[i,3] > ATAC_peak$start & UBTF[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_UBTF) != 0){
    ATAC_peak_UBTF$V5 <- UBTF[i,1]
    ATAC_peak_UBTF$V6 <- UBTF[i,2]
    ATAC_peak_UBTF$V7 <- UBTF[i,3]
    ATAC_peak_UBTF$V8 <- UBTF[i,6]
    ATAC_peak_UBTF$V9 <- UBTF[i,7]
    ATAC_peak_UBTF_list = rbind(ATAC_peak_UBTF_list,ATAC_peak_UBTF)
  }
}
write.table(ATAC_peak_UBTF_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_UBTF_list.txt",sep="\t")

ATAC_peak_UBTF_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_UBTF_list.txt",sep="\t")

UBTF_features=list("UBTF" = ATAC_peak_UBTF_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = UBTF_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MLLT1
MLLT1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF404UXE_MLLT1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MLLT1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF404UXE_MLLT1.txt",sep = "\t")

MLLT1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF404UXE_MLLT1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MLLT1)){
  x=paste(MLLT1[i,1],MLLT1[i,2],MLLT1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MLLT1_features = list(list$V1)

ATAC_peak_MLLT1_list = data.frame()
for(i in 1:nrow(MLLT1)){
  ATAC_peak_MLLT1<- ATAC_peak[ATAC_peak$x1 == MLLT1[i,1] & ((MLLT1[i,2] > ATAC_peak$start & MLLT1[i,2] < ATAC_peak$end) | (MLLT1[i,3] > ATAC_peak$start & MLLT1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MLLT1) != 0){
    ATAC_peak_MLLT1$V5 <- MLLT1[i,1]
    ATAC_peak_MLLT1$V6 <- MLLT1[i,2]
    ATAC_peak_MLLT1$V7 <- MLLT1[i,3]
    ATAC_peak_MLLT1$V8 <- MLLT1[i,6]
    ATAC_peak_MLLT1$V9 <- MLLT1[i,7]
    ATAC_peak_MLLT1_list = rbind(ATAC_peak_MLLT1_list,ATAC_peak_MLLT1)
  }
}
write.table(ATAC_peak_MLLT1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MLLT1_list.txt",sep="\t")

ATAC_peak_MLLT1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MLLT1_list.txt",sep="\t")

MLLT1_features=list("MLLT1" = ATAC_peak_MLLT1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MLLT1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##POLR2G
POLR2G = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF407HBE_POLR2G.bed.gz",extraCols=extraCols_narrowPeak)
write.table(POLR2G,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF407HBE_POLR2G.txt",sep = "\t")

POLR2G = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF407HBE_POLR2G.txt",sep = "\t")

list = vector()
for (i in  1:nrow(POLR2G)){
  x=paste(POLR2G[i,1],POLR2G[i,2],POLR2G[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

POLR2G_features = list(list$V1)

ATAC_peak_POLR2G_list = data.frame()
for(i in 1:nrow(POLR2G)){
  ATAC_peak_POLR2G<- ATAC_peak[ATAC_peak$x1 == POLR2G[i,1] & ((POLR2G[i,2] > ATAC_peak$start & POLR2G[i,2] < ATAC_peak$end) | (POLR2G[i,3] > ATAC_peak$start & POLR2G[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_POLR2G) != 0){
    ATAC_peak_POLR2G$V5 <- POLR2G[i,1]
    ATAC_peak_POLR2G$V6 <- POLR2G[i,2]
    ATAC_peak_POLR2G$V7 <- POLR2G[i,3]
    ATAC_peak_POLR2G$V8 <- POLR2G[i,6]
    ATAC_peak_POLR2G$V9 <- POLR2G[i,7]
    ATAC_peak_POLR2G_list = rbind(ATAC_peak_POLR2G_list,ATAC_peak_POLR2G)
  }
}
write.table(ATAC_peak_POLR2G_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLR2G_list.txt",sep="\t")

ATAC_peak_POLR2G_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLR2G_list.txt",sep="\t")

POLR2G_features=list("POLR2G" = ATAC_peak_POLR2G_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = POLR2G_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CHD1
CHD1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF408NUX_CHD1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CHD1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF408NUX_CHD1.txt",sep = "\t")

CHD1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF408NUX_CHD1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CHD1)){
  x=paste(CHD1[i,1],CHD1[i,2],CHD1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CHD1_features = list(list$V1)

ATAC_peak_CHD1_list = data.frame()
for(i in 1:nrow(CHD1)){
  ATAC_peak_CHD1<- ATAC_peak[ATAC_peak$x1 == CHD1[i,1] & ((CHD1[i,2] > ATAC_peak$start & CHD1[i,2] < ATAC_peak$end) | (CHD1[i,3] > ATAC_peak$start & CHD1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CHD1) != 0){
    ATAC_peak_CHD1$V5 <- CHD1[i,1]
    ATAC_peak_CHD1$V6 <- CHD1[i,2]
    ATAC_peak_CHD1$V7 <- CHD1[i,3]
    ATAC_peak_CHD1$V8 <- CHD1[i,6]
    ATAC_peak_CHD1$V9 <- CHD1[i,7]
    ATAC_peak_CHD1_list = rbind(ATAC_peak_CHD1_list,ATAC_peak_CHD1)
  }
}
write.table(ATAC_peak_CHD1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHD1_list.txt",sep="\t")

ATAC_peak_CHD1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHD1_list.txt",sep="\t")

CHD1_features=list("CHD1" = ATAC_peak_CHD1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CHD1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZNF512
ZNF512 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF415LUE_ZNF512.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF512,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF415LUE_ZNF512.txt",sep = "\t")

ZNF512 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF415LUE_ZNF512.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF512)){
  x=paste(ZNF512[i,1],ZNF512[i,2],ZNF512[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF512_features = list(list$V1)

ATAC_peak_ZNF512_list = data.frame()
for(i in 1:nrow(ZNF512)){
  ATAC_peak_ZNF512<- ATAC_peak[ATAC_peak$x1 == ZNF512[i,1] & ((ZNF512[i,2] > ATAC_peak$start & ZNF512[i,2] < ATAC_peak$end) | (ZNF512[i,3] > ATAC_peak$start & ZNF512[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF512) != 0){
    ATAC_peak_ZNF512$V5 <- ZNF512[i,1]
    ATAC_peak_ZNF512$V6 <- ZNF512[i,2]
    ATAC_peak_ZNF512$V7 <- ZNF512[i,3]
    ATAC_peak_ZNF512$V8 <- ZNF512[i,6]
    ATAC_peak_ZNF512$V9 <- ZNF512[i,7]
    ATAC_peak_ZNF512_list = rbind(ATAC_peak_ZNF512_list,ATAC_peak_ZNF512)
  }
}
write.table(ATAC_peak_ZNF512_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF512_list.txt",sep="\t")

ATAC_peak_ZNF512_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF512_list.txt",sep="\t")

ZNF512_features=list("ZNF512" = ATAC_peak_ZNF512_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF512_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HNRNPL
HNRNPL = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF423RKR_HNRNPL.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HNRNPL,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF423RKR_HNRNPL.txt",sep = "\t")

HNRNPL = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF423RKR_HNRNPL.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HNRNPL)){
  x=paste(HNRNPL[i,1],HNRNPL[i,2],HNRNPL[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HNRNPL_features = list(list$V1)

ATAC_peak_HNRNPL_list = data.frame()
for(i in 1:nrow(HNRNPL)){
  ATAC_peak_HNRNPL<- ATAC_peak[ATAC_peak$x1 == HNRNPL[i,1] & ((HNRNPL[i,2] > ATAC_peak$start & HNRNPL[i,2] < ATAC_peak$end) | (HNRNPL[i,3] > ATAC_peak$start & HNRNPL[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HNRNPL) != 0){
    ATAC_peak_HNRNPL$V5 <- HNRNPL[i,1]
    ATAC_peak_HNRNPL$V6 <- HNRNPL[i,2]
    ATAC_peak_HNRNPL$V7 <- HNRNPL[i,3]
    ATAC_peak_HNRNPL$V8 <- HNRNPL[i,6]
    ATAC_peak_HNRNPL$V9 <- HNRNPL[i,7]
    ATAC_peak_HNRNPL_list = rbind(ATAC_peak_HNRNPL_list,ATAC_peak_HNRNPL)
  }
}
write.table(ATAC_peak_HNRNPL_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HNRNPL_list.txt",sep="\t")

ATAC_peak_HNRNPL_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HNRNPL_list.txt",sep="\t")

HNRNPL_features=list("HNRNPL" = ATAC_peak_HNRNPL_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HNRNPL_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZNF830
ZNF830 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF426LEV_ZNF830.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF830,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF426LEV_ZNF830.txt",sep = "\t")

ZNF830 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF426LEV_ZNF830.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF830)){
  x=paste(ZNF830[i,1],ZNF830[i,2],ZNF830[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF830_features = list(list$V1)

ATAC_peak_ZNF830_list = data.frame()
for(i in 1:nrow(ZNF830)){
  ATAC_peak_ZNF830<- ATAC_peak[ATAC_peak$x1 == ZNF830[i,1] & ((ZNF830[i,2] > ATAC_peak$start & ZNF830[i,2] < ATAC_peak$end) | (ZNF830[i,3] > ATAC_peak$start & ZNF830[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF830) != 0){
    ATAC_peak_ZNF830$V5 <- ZNF830[i,1]
    ATAC_peak_ZNF830$V6 <- ZNF830[i,2]
    ATAC_peak_ZNF830$V7 <- ZNF830[i,3]
    ATAC_peak_ZNF830$V8 <- ZNF830[i,6]
    ATAC_peak_ZNF830$V9 <- ZNF830[i,7]
    ATAC_peak_ZNF830_list = rbind(ATAC_peak_ZNF830_list,ATAC_peak_ZNF830)
  }
}
write.table(ATAC_peak_ZNF830_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF830_list.txt",sep="\t")

ATAC_peak_ZNF830_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF830_list.txt",sep="\t")

ZNF830_features=list("ZNF830" = ATAC_peak_ZNF830_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF830_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RAD21
RAD21 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF439DYW_RAD21.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RAD21,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF439DYW_RAD21.txt",sep = "\t")

RAD21 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF439DYW_RAD21.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RAD21)){
  x=paste(RAD21[i,1],RAD21[i,2],RAD21[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RAD21_features = list(list$V1)

ATAC_peak_RAD21_list = data.frame()
for(i in 1:nrow(RAD21)){
  ATAC_peak_RAD21<- ATAC_peak[ATAC_peak$x1 == RAD21[i,1] & ((RAD21[i,2] > ATAC_peak$start & RAD21[i,2] < ATAC_peak$end) | (RAD21[i,3] > ATAC_peak$start & RAD21[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RAD21) != 0){
    ATAC_peak_RAD21$V5 <- RAD21[i,1]
    ATAC_peak_RAD21$V6 <- RAD21[i,2]
    ATAC_peak_RAD21$V7 <- RAD21[i,3]
    ATAC_peak_RAD21$V8 <- RAD21[i,6]
    ATAC_peak_RAD21$V9 <- RAD21[i,7]
    ATAC_peak_RAD21_list = rbind(ATAC_peak_RAD21_list,ATAC_peak_RAD21)
  }
}
write.table(ATAC_peak_RAD21_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RAD21_list.txt",sep="\t")

ATAC_peak_RAD21_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RAD21_list.txt",sep="\t")

RAD21_features=list("RAD21" = ATAC_peak_RAD21_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RAD21_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##GTF2I
GTF2I = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF439JZB_GTF2I.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GTF2I,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF439JZB_GTF2I.txt",sep = "\t")

GTF2I = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF439JZB_GTF2I.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GTF2I)){
  x=paste(GTF2I[i,1],GTF2I[i,2],GTF2I[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GTF2I_features = list(list$V1)

ATAC_peak_GTF2I_list = data.frame()
for(i in 1:nrow(GTF2I)){
  ATAC_peak_GTF2I<- ATAC_peak[ATAC_peak$x1 == GTF2I[i,1] & ((GTF2I[i,2] > ATAC_peak$start & GTF2I[i,2] < ATAC_peak$end) | (GTF2I[i,3] > ATAC_peak$start & GTF2I[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GTF2I) != 0){
    ATAC_peak_GTF2I$V5 <- GTF2I[i,1]
    ATAC_peak_GTF2I$V6 <- GTF2I[i,2]
    ATAC_peak_GTF2I$V7 <- GTF2I[i,3]
    ATAC_peak_GTF2I$V8 <- GTF2I[i,6]
    ATAC_peak_GTF2I$V9 <- GTF2I[i,7]
    ATAC_peak_GTF2I_list = rbind(ATAC_peak_GTF2I_list,ATAC_peak_GTF2I)
  }
}
write.table(ATAC_peak_GTF2I_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF2I_list.txt",sep="\t")

ATAC_peak_GTF2I_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF2I_list.txt",sep="\t")

GTF2I_features=list("GTF2I" = ATAC_peak_GTF2I_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GTF2I_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZNF274
ZNF274 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF446CGL_ZNF274.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF274,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF446CGL_ZNF274.txt",sep = "\t")

ZNF274 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF446CGL_ZNF274.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF274)){
  x=paste(ZNF274[i,1],ZNF274[i,2],ZNF274[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF274_features = list(list$V1)

ATAC_peak_ZNF274_list = data.frame()
for(i in 1:nrow(ZNF274)){
  ATAC_peak_ZNF274<- ATAC_peak[ATAC_peak$x1 == ZNF274[i,1] & ((ZNF274[i,2] > ATAC_peak$start & ZNF274[i,2] < ATAC_peak$end) | (ZNF274[i,3] > ATAC_peak$start & ZNF274[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF274) != 0){
    ATAC_peak_ZNF274$V5 <- ZNF274[i,1]
    ATAC_peak_ZNF274$V6 <- ZNF274[i,2]
    ATAC_peak_ZNF274$V7 <- ZNF274[i,3]
    ATAC_peak_ZNF274$V8 <- ZNF274[i,6]
    ATAC_peak_ZNF274$V9 <- ZNF274[i,7]
    ATAC_peak_ZNF274_list = rbind(ATAC_peak_ZNF274_list,ATAC_peak_ZNF274)
  }
}
write.table(ATAC_peak_ZNF274_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF274_list.txt",sep="\t")

ATAC_peak_ZNF274_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF274_list.txt",sep="\t")

ZNF274_features=list("ZNF274" = ATAC_peak_ZNF274_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF274_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZFX
ZFX = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF451YXQ_ZFX.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZFX,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF451YXQ_ZFX.txt",sep = "\t")

ZFX = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF451YXQ_ZFX.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZFX)){
  x=paste(ZFX[i,1],ZFX[i,2],ZFX[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZFX_features = list(list$V1)

ATAC_peak_ZFX_list = data.frame()
for(i in 1:nrow(ZFX)){
  ATAC_peak_ZFX<- ATAC_peak[ATAC_peak$x1 == ZFX[i,1] & ((ZFX[i,2] > ATAC_peak$start & ZFX[i,2] < ATAC_peak$end) | (ZFX[i,3] > ATAC_peak$start & ZFX[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZFX) != 0){
    ATAC_peak_ZFX$V5 <- ZFX[i,1]
    ATAC_peak_ZFX$V6 <- ZFX[i,2]
    ATAC_peak_ZFX$V7 <- ZFX[i,3]
    ATAC_peak_ZFX$V8 <- ZFX[i,6]
    ATAC_peak_ZFX$V9 <- ZFX[i,7]
    ATAC_peak_ZFX_list = rbind(ATAC_peak_ZFX_list,ATAC_peak_ZFX)
  }
}
write.table(ATAC_peak_ZFX_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZFX_list.txt",sep="\t")

ATAC_peak_ZFX_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZFX_list.txt",sep="\t")

ZFX_features=list("ZFX" = ATAC_peak_ZFX_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZFX_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CHD2
CHD2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF454FVE_CHD2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CHD2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF454FVE_CHD2.txt",sep = "\t")

CHD2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF454FVE_CHD2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CHD2)){
  x=paste(CHD2[i,1],CHD2[i,2],CHD2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CHD2_features = list(list$V1)

ATAC_peak_CHD2_list = data.frame()
for(i in 1:nrow(CHD2)){
  ATAC_peak_CHD2<- ATAC_peak[ATAC_peak$x1 == CHD2[i,1] & ((CHD2[i,2] > ATAC_peak$start & CHD2[i,2] < ATAC_peak$end) | (CHD2[i,3] > ATAC_peak$start & CHD2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CHD2) != 0){
    ATAC_peak_CHD2$V5 <- CHD2[i,1]
    ATAC_peak_CHD2$V6 <- CHD2[i,2]
    ATAC_peak_CHD2$V7 <- CHD2[i,3]
    ATAC_peak_CHD2$V8 <- CHD2[i,6]
    ATAC_peak_CHD2$V9 <- CHD2[i,7]
    ATAC_peak_CHD2_list = rbind(ATAC_peak_CHD2_list,ATAC_peak_CHD2)
  }
}
write.table(ATAC_peak_CHD2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHD2_list.txt",sep="\t")

ATAC_peak_CHD2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHD2_list.txt",sep="\t")

CHD2_features=list("CHD2" = ATAC_peak_CHD2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CHD2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MXI1
MXI1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF457COE_MXI1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MXI1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF457COE_MXI1_age.txt",sep = "\t")

MXI1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF457COE_MXI1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MXI1)){
  x=paste(MXI1[i,1],MXI1[i,2],MXI1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MXI1_features = list(list$V1)

ATAC_peak_MXI1_list = data.frame()
for(i in 1:nrow(MXI1)){
  ATAC_peak_MXI1<- ATAC_peak[ATAC_peak$x1 == MXI1[i,1] & ((MXI1[i,2] > ATAC_peak$start & MXI1[i,2] < ATAC_peak$end) | (MXI1[i,3] > ATAC_peak$start & MXI1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MXI1) != 0){
    ATAC_peak_MXI1$V5 <- MXI1[i,1]
    ATAC_peak_MXI1$V6 <- MXI1[i,2]
    ATAC_peak_MXI1$V7 <- MXI1[i,3]
    ATAC_peak_MXI1$V8 <- MXI1[i,6]
    ATAC_peak_MXI1$V9 <- MXI1[i,7]
    ATAC_peak_MXI1_list = rbind(ATAC_peak_MXI1_list,ATAC_peak_MXI1)
  }
}
write.table(ATAC_peak_MXI1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MXI1_list.txt",sep="\t")

ATAC_peak_MXI1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MXI1_list.txt",sep="\t")

MXI1_features=list("MXI1" = ATAC_peak_MXI1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MXI1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##DEAF1
DEAF1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF457VMF_DEAF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(DEAF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF457VMF_DEAF1.txt",sep = "\t")

DEAF1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF457VMF_DEAF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(DEAF1)){
  x=paste(DEAF1[i,1],DEAF1[i,2],DEAF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

DEAF1_features = list(list$V1)

ATAC_peak_DEAF1_list = data.frame()
for(i in 1:nrow(DEAF1)){
  ATAC_peak_DEAF1<- ATAC_peak[ATAC_peak$x1 == DEAF1[i,1] & ((DEAF1[i,2] > ATAC_peak$start & DEAF1[i,2] < ATAC_peak$end) | (DEAF1[i,3] > ATAC_peak$start & DEAF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_DEAF1) != 0){
    ATAC_peak_DEAF1$V5 <- DEAF1[i,1]
    ATAC_peak_DEAF1$V6 <- DEAF1[i,2]
    ATAC_peak_DEAF1$V7 <- DEAF1[i,3]
    ATAC_peak_DEAF1$V8 <- DEAF1[i,6]
    ATAC_peak_DEAF1$V9 <- DEAF1[i,7]
    ATAC_peak_DEAF1_list = rbind(ATAC_peak_DEAF1_list,ATAC_peak_DEAF1)
  }
}
write.table(ATAC_peak_DEAF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DEAF1_list.txt",sep="\t")

ATAC_peak_DEAF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DEAF1_list.txt",sep="\t")

DEAF1_features=list("DEAF1" = ATAC_peak_DEAF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = DEAF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ID3
ID3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF458GQT_ID3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ID3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF458GQT_ID3.txt",sep = "\t")

ID3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF458GQT_ID3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ID3)){
  x=paste(ID3[i,1],ID3[i,2],ID3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ID3_features = list(list$V1)

ATAC_peak_ID3_list = data.frame()
for(i in 1:nrow(ID3)){
  ATAC_peak_ID3<- ATAC_peak[ATAC_peak$x1 ==ID3[i,1] & ((ID3[i,2] > ATAC_peak$start &ID3[i,2] < ATAC_peak$end) | (ID3[i,3] > ATAC_peak$start &ID3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ID3) != 0){
    ATAC_peak_ID3$V5 <-ID3[i,1]
    ATAC_peak_ID3$V6 <- ID3[i,2]
    ATAC_peak_ID3$V7 <- ID3[i,3]
    ATAC_peak_ID3$V8 <- ID3[i,6]
    ATAC_peak_ID3$V9 <- ID3[i,7]
    ATAC_peak_ID3_list = rbind(ATAC_peak_ID3_list,ATAC_peak_ID3)
  }
}
write.table(ATAC_peak_ID3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ID3_list.txt",sep="\t")

ATAC_peak_ID3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ID3_list.txt",sep="\t")

ID3_features=list("ID3" = ATAC_peak_ID3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ID3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MIER1
MIER1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF465EPI_MIER1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MIER1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF465EPI_MIER1.txt",sep = "\t")

MIER1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF465EPI_MIER1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MIER1)){
  x=paste(MIER1[i,1],MIER1[i,2],MIER1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MIER1_features = list(list$V1)

ATAC_peak_MIER1_list = data.frame()
for(i in 1:nrow(MIER1)){
  ATAC_peak_MIER1<- ATAC_peak[ATAC_peak$x1 ==MIER1[i,1] & ((MIER1[i,2] > ATAC_peak$start &MIER1[i,2] < ATAC_peak$end) | (MIER1[i,3] > ATAC_peak$start &MIER1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MIER1) != 0){
    ATAC_peak_MIER1$V5 <-MIER1[i,1]
    ATAC_peak_MIER1$V6 <- MIER1[i,2]
    ATAC_peak_MIER1$V7 <- MIER1[i,3]
    ATAC_peak_MIER1$V8 <- MIER1[i,6]
    ATAC_peak_MIER1$V9 <- MIER1[i,7]
    ATAC_peak_MIER1_list = rbind(ATAC_peak_MIER1_list,ATAC_peak_MIER1)
  }
}
write.table(ATAC_peak_MIER1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MIER1_list.txt",sep="\t")

ATAC_peak_MIER1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MIER1_list.txt",sep="\t")

MIER1_features=list("MIER1" = ATAC_peak_MIER1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MIER1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CLOCK
CLOCK = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF475CFL_CLOCK_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CLOCK,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF475CFL_CLOCK_age.txt",sep = "\t")

CLOCK = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF475CFL_CLOCK_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CLOCK)){
  x=paste(CLOCK[i,1],CLOCK[i,2],CLOCK[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CLOCK_features = list(list$V1)

ATAC_peak_CLOCK_list = data.frame()
for(i in 1:nrow(CLOCK)){
  ATAC_peak_CLOCK<- ATAC_peak[ATAC_peak$x1 ==CLOCK[i,1] & ((CLOCK[i,2] > ATAC_peak$start &CLOCK[i,2] < ATAC_peak$end) | (CLOCK[i,3] > ATAC_peak$start &CLOCK[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CLOCK) != 0){
    ATAC_peak_CLOCK$V5 <-CLOCK[i,1]
    ATAC_peak_CLOCK$V6 <- CLOCK[i,2]
    ATAC_peak_CLOCK$V7 <- CLOCK[i,3]
    ATAC_peak_CLOCK$V8 <- CLOCK[i,6]
    ATAC_peak_CLOCK$V9 <- CLOCK[i,7]
    ATAC_peak_CLOCK_list = rbind(ATAC_peak_CLOCK_list,ATAC_peak_CLOCK)
  }
}
write.table(ATAC_peak_CLOCK_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CLOCK_list.txt",sep="\t")

ATAC_peak_CLOCK_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CLOCK_list.txt",sep="\t")

CLOCK_features=list("CLOCK" = ATAC_peak_CLOCK_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CLOCK_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SRSF1
SRSF1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF475SHO_SRSF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SRSF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF475SHO_SRSF1.txt",sep = "\t")

SRSF1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF475SHO_SRSF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SRSF1)){
  x=paste(SRSF1[i,1],SRSF1[i,2],SRSF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SRSF1_features = list(list$V1)

ATAC_peak_SRSF1_list = data.frame()
for(i in 1:nrow(SRSF1)){
  ATAC_peak_SRSF1<- ATAC_peak[ATAC_peak$x1 ==SRSF1[i,1] & ((SRSF1[i,2] > ATAC_peak$start &SRSF1[i,2] < ATAC_peak$end) | (SRSF1[i,3] > ATAC_peak$start &SRSF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SRSF1) != 0){
    ATAC_peak_SRSF1$V5 <-SRSF1[i,1]
    ATAC_peak_SRSF1$V6 <- SRSF1[i,2]
    ATAC_peak_SRSF1$V7 <- SRSF1[i,3]
    ATAC_peak_SRSF1$V8 <- SRSF1[i,6]
    ATAC_peak_SRSF1$V9 <- SRSF1[i,7]
    ATAC_peak_SRSF1_list = rbind(ATAC_peak_SRSF1_list,ATAC_peak_SRSF1)
  }
}
write.table(ATAC_peak_SRSF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SRSF1_list.txt",sep="\t")

ATAC_peak_SRSF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SRSF1_list.txt",sep="\t")

SRSF1_features=list("SRSF1" = ATAC_peak_SRSF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SRSF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SKIL
SKIL = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF480QGO_SKIL.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SKIL,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF480QGO_SKIL.txt",sep = "\t")

SKIL = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF480QGO_SKIL.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SKIL)){
  x=paste(SKIL[i,1],SKIL[i,2],SKIL[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SKIL_features = list(list$V1)

ATAC_peak_SKIL_list = data.frame()
for(i in 1:nrow(SKIL)){
  ATAC_peak_SKIL<- ATAC_peak[ATAC_peak$x1 ==SKIL[i,1] & ((SKIL[i,2] > ATAC_peak$start &SKIL[i,2] < ATAC_peak$end) | (SKIL[i,3] > ATAC_peak$start &SKIL[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SKIL) != 0){
    ATAC_peak_SKIL$V5 <-SKIL[i,1]
    ATAC_peak_SKIL$V6 <- SKIL[i,2]
    ATAC_peak_SKIL$V7 <- SKIL[i,3]
    ATAC_peak_SKIL$V8 <- SKIL[i,6]
    ATAC_peak_SKIL$V9 <- SKIL[i,7]
    ATAC_peak_SKIL_list = rbind(ATAC_peak_SKIL_list,ATAC_peak_SKIL)
  }
}
write.table(ATAC_peak_SKIL_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SKIL_list.txt",sep="\t")

ATAC_peak_SKIL_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SKIL_list.txt",sep="\t")

SKIL_features=list("SKIL" = ATAC_peak_SKIL_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SKIL_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##KAT8
KAT8 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF488QFU_KAT8.bed.gz",extraCols=extraCols_narrowPeak)
write.table(KAT8,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF488QFU_KAT8.txt",sep = "\t")

KAT8 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF488QFU_KAT8.txt",sep = "\t")

list = vector()
for (i in  1:nrow(KAT8)){
  x=paste(KAT8[i,1],KAT8[i,2],KAT8[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

KAT8_features = list(list$V1)

ATAC_peak_KAT8_list = data.frame()
for(i in 1:nrow(KAT8)){
  ATAC_peak_KAT8<- ATAC_peak[ATAC_peak$x1 ==KAT8[i,1] & ((KAT8[i,2] > ATAC_peak$start &KAT8[i,2] < ATAC_peak$end) | (KAT8[i,3] > ATAC_peak$start &KAT8[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_KAT8) != 0){
    ATAC_peak_KAT8$V5 <-KAT8[i,1]
    ATAC_peak_KAT8$V6 <- KAT8[i,2]
    ATAC_peak_KAT8$V7 <- KAT8[i,3]
    ATAC_peak_KAT8$V8 <- KAT8[i,6]
    ATAC_peak_KAT8$V9 <- KAT8[i,7]
    ATAC_peak_KAT8_list = rbind(ATAC_peak_KAT8_list,ATAC_peak_KAT8)
  }
}
write.table(ATAC_peak_KAT8_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KAT8_list.txt",sep="\t")

ATAC_peak_KAT8_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KAT8_list.txt",sep="\t")

KAT8_features=list("KAT8" = ATAC_peak_KAT8_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = KAT8_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ASH1L
ASH1L = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF491RRP_ASH1L.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ASH1L,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF491RRP_ASH1L.txt",sep = "\t")

ASH1L = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF491RRP_ASH1L.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ASH1L)){
  x=paste(ASH1L[i,1],ASH1L[i,2],ASH1L[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ASH1L_features = list(list$V1)

ATAC_peak_ASH1L_list = data.frame()
for(i in 1:nrow(ASH1L)){
  ATAC_peak_ASH1L<- ATAC_peak[ATAC_peak$x1 ==ASH1L[i,1] & ((ASH1L[i,2] > ATAC_peak$start &ASH1L[i,2] < ATAC_peak$end) | (ASH1L[i,3] > ATAC_peak$start &ASH1L[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ASH1L) != 0){
    ATAC_peak_ASH1L$V5 <-ASH1L[i,1]
    ATAC_peak_ASH1L$V6 <- ASH1L[i,2]
    ATAC_peak_ASH1L$V7 <- ASH1L[i,3]
    ATAC_peak_ASH1L$V8 <- ASH1L[i,6]
    ATAC_peak_ASH1L$V9 <- ASH1L[i,7]
    ATAC_peak_ASH1L_list = rbind(ATAC_peak_ASH1L_list,ATAC_peak_ASH1L)
  }
}
write.table(ATAC_peak_ASH1L_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ASH1L_list.txt",sep="\t")

ATAC_peak_ASH1L_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ASH1L_list.txt",sep="\t")

ASH1L_features=list("ASH1L" = ATAC_peak_ASH1L_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ASH1L_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##POLR2AphosphoS5
POLR2AphosphoS5 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF493NUB_POLR2AphosphoS5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(POLR2AphosphoS5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF493NUB_POLR2AphosphoS5.txt",sep = "\t")

POLR2AphosphoS5 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF493NUB_POLR2AphosphoS5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(POLR2AphosphoS5)){
  x=paste(POLR2AphosphoS5[i,1],POLR2AphosphoS5[i,2],POLR2AphosphoS5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

POLR2AphosphoS5_features = list(list$V1)

ATAC_peak_POLR2AphosphoS5_list = data.frame()
for(i in 1:nrow(POLR2AphosphoS5)){
  ATAC_peak_POLR2AphosphoS5<- ATAC_peak[ATAC_peak$x1 ==POLR2AphosphoS5[i,1] & ((POLR2AphosphoS5[i,2] > ATAC_peak$start &POLR2AphosphoS5[i,2] < ATAC_peak$end) | (POLR2AphosphoS5[i,3] > ATAC_peak$start &POLR2AphosphoS5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_POLR2AphosphoS5) != 0){
    ATAC_peak_POLR2AphosphoS5$V5 <-POLR2AphosphoS5[i,1]
    ATAC_peak_POLR2AphosphoS5$V6 <- POLR2AphosphoS5[i,2]
    ATAC_peak_POLR2AphosphoS5$V7 <- POLR2AphosphoS5[i,3]
    ATAC_peak_POLR2AphosphoS5$V8 <- POLR2AphosphoS5[i,6]
    ATAC_peak_POLR2AphosphoS5$V9 <- POLR2AphosphoS5[i,7]
    ATAC_peak_POLR2AphosphoS5_list = rbind(ATAC_peak_POLR2AphosphoS5_list,ATAC_peak_POLR2AphosphoS5)
  }
}
write.table(ATAC_peak_POLR2AphosphoS5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLR2AphosphoS5_list.txt",sep="\t")

ATAC_peak_POLR2AphosphoS5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLR2AphosphoS5_list.txt",sep="\t")

POLR2AphosphoS5_features=list("POLR2AphosphoS5" = ATAC_peak_POLR2AphosphoS5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = POLR2AphosphoS5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ETS2
ETS2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF493PSP_ETS2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ETS2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF493PSP_ETS2.txt",sep = "\t")

ETS2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF493PSP_ETS2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ETS2)){
  x=paste(ETS2[i,1],ETS2[i,2],ETS2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ETS2_features = list(list$V1)

ATAC_peak_ETS2_list = data.frame()
for(i in 1:nrow(ETS2)){
  ATAC_peak_ETS2<- ATAC_peak[ATAC_peak$x1 ==ETS2[i,1] & ((ETS2[i,2] > ATAC_peak$start &ETS2[i,2] < ATAC_peak$end) | (ETS2[i,3] > ATAC_peak$start &ETS2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ETS2) != 0){
    ATAC_peak_ETS2$V5 <-ETS2[i,1]
    ATAC_peak_ETS2$V6 <- ETS2[i,2]
    ATAC_peak_ETS2$V7 <- ETS2[i,3]
    ATAC_peak_ETS2$V8 <- ETS2[i,6]
    ATAC_peak_ETS2$V9 <- ETS2[i,7]
    ATAC_peak_ETS2_list = rbind(ATAC_peak_ETS2_list,ATAC_peak_ETS2)
  }
}
write.table(ATAC_peak_ETS2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ETS2_list.txt",sep="\t")

ATAC_peak_ETS2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ETS2_list.txt",sep="\t")

ETS2_features=list("ETS2" = ATAC_peak_ETS2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ETS2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##BDP1
BDP1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF496QIZ_BDP1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BDP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF496QIZ_BDP1.txt",sep = "\t")

BDP1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF496QIZ_BDP1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BDP1)){
  x=paste(BDP1[i,1],BDP1[i,2],BDP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BDP1_features = list(list$V1)

ATAC_peak_BDP1_list = data.frame()
for(i in 1:nrow(BDP1)){
  ATAC_peak_BDP1<- ATAC_peak[ATAC_peak$x1 ==BDP1[i,1] & ((BDP1[i,2] > ATAC_peak$start &BDP1[i,2] < ATAC_peak$end) | (BDP1[i,3] > ATAC_peak$start &BDP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BDP1) != 0){
    ATAC_peak_BDP1$V5 <-BDP1[i,1]
    ATAC_peak_BDP1$V6 <- BDP1[i,2]
    ATAC_peak_BDP1$V7 <- BDP1[i,3]
    ATAC_peak_BDP1$V8 <- BDP1[i,6]
    ATAC_peak_BDP1$V9 <- BDP1[i,7]
    ATAC_peak_BDP1_list = rbind(ATAC_peak_BDP1_list,ATAC_peak_BDP1)
  }
}
write.table(ATAC_peak_BDP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BDP1_list.txt",sep="\t")

ATAC_peak_BDP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BDP1_list.txt",sep="\t")

BDP1_features=list("BDP1" = ATAC_peak_BDP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BDP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##JUN
JUN = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF498BRE_JUN_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(JUN,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF498BRE_JUN_age.txt",sep = "\t")

JUN = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF498BRE_JUN_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(JUN)){
  x=paste(JUN[i,1],JUN[i,2],JUN[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

JUN_features = list(list$V1)

ATAC_peak_JUN_list = data.frame()
for(i in 1:nrow(JUN)){
  ATAC_peak_JUN<- ATAC_peak[ATAC_peak$x1 ==JUN[i,1] & ((JUN[i,2] > ATAC_peak$start &JUN[i,2] < ATAC_peak$end) | (JUN[i,3] > ATAC_peak$start &JUN[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_JUN) != 0){
    ATAC_peak_JUN$V5 <-JUN[i,1]
    ATAC_peak_JUN$V6 <- JUN[i,2]
    ATAC_peak_JUN$V7 <- JUN[i,3]
    ATAC_peak_JUN$V8 <- JUN[i,6]
    ATAC_peak_JUN$V9 <- JUN[i,7]
    ATAC_peak_JUN_list = rbind(ATAC_peak_JUN_list,ATAC_peak_JUN)
  }
}
write.table(ATAC_peak_JUN_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_JUN_list.txt",sep="\t")

ATAC_peak_JUN_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_JUN_list.txt",sep="\t")

JUN_features=list("JUN" = ATAC_peak_JUN_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = JUN_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##HDAC6
HDAC6 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF499GKO_HDAC6.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HDAC6,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF499GKO_HDAC6.txt",sep = "\t")

HDAC6 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF499GKO_HDAC6.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HDAC6)){
  x=paste(HDAC6[i,1],HDAC6[i,2],HDAC6[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HDAC6_features = list(list$V1)

ATAC_peak_HDAC6_list = data.frame()
for(i in 1:nrow(HDAC6)){
  ATAC_peak_HDAC6<- ATAC_peak[ATAC_peak$x1 ==HDAC6[i,1] & ((HDAC6[i,2] > ATAC_peak$start &HDAC6[i,2] < ATAC_peak$end) | (HDAC6[i,3] > ATAC_peak$start &HDAC6[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HDAC6) != 0){
    ATAC_peak_HDAC6$V5 <-HDAC6[i,1]
    ATAC_peak_HDAC6$V6 <- HDAC6[i,2]
    ATAC_peak_HDAC6$V7 <- HDAC6[i,3]
    ATAC_peak_HDAC6$V8 <- HDAC6[i,6]
    ATAC_peak_HDAC6$V9 <- HDAC6[i,7]
    ATAC_peak_HDAC6_list = rbind(ATAC_peak_HDAC6_list,ATAC_peak_HDAC6)
  }
}
write.table(ATAC_peak_HDAC6_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC6_list.txt",sep="\t")

ATAC_peak_HDAC6_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC6_list.txt",sep="\t")

HDAC6_features=list("HDAC6" = ATAC_peak_HDAC6_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HDAC6_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

save(scmulti2, file = "/Users/ramzipit/Desktop/scmulti2_temp.Rdata")
##SIN3A
SIN3A = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF503JBM_SIN3A_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SIN3A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF503JBM_SIN3A_age.txt",sep = "\t")

SIN3A = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF503JBM_SIN3A_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SIN3A)){
  x=paste(SIN3A[i,1],SIN3A[i,2],SIN3A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SIN3A_features = list(list$V1)

ATAC_peak_SIN3A_list = data.frame()
for(i in 1:nrow(SIN3A)){
  ATAC_peak_SIN3A<- ATAC_peak[ATAC_peak$x1 ==SIN3A[i,1] & ((SIN3A[i,2] > ATAC_peak$start &SIN3A[i,2] < ATAC_peak$end) | (SIN3A[i,3] > ATAC_peak$start &SIN3A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SIN3A) != 0){
    ATAC_peak_SIN3A$V5 <-SIN3A[i,1]
    ATAC_peak_SIN3A$V6 <- SIN3A[i,2]
    ATAC_peak_SIN3A$V7 <- SIN3A[i,3]
    ATAC_peak_SIN3A$V8 <- SIN3A[i,6]
    ATAC_peak_SIN3A$V9 <- SIN3A[i,7]
    ATAC_peak_SIN3A_list = rbind(ATAC_peak_SIN3A_list,ATAC_peak_SIN3A)
  }
}
write.table(ATAC_peak_SIN3A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIN3A_list.txt",sep="\t")

ATAC_peak_SIN3A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIN3A_list.txt",sep="\t")

SIN3A_features=list("SIN3A" = ATAC_peak_SIN3A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SIN3A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##HNRNPK
HNRNPK = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF505RNR_HNRNPK.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HNRNPK,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF505RNR_HNRNPK.txt",sep = "\t")

HNRNPK = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF505RNR_HNRNPK.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HNRNPK)){
  x=paste(HNRNPK[i,1],HNRNPK[i,2],HNRNPK[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HNRNPK_features = list(list$V1)

ATAC_peak_HNRNPK_list = data.frame()
for(i in 1:nrow(HNRNPK)){
  ATAC_peak_HNRNPK<- ATAC_peak[ATAC_peak$x1 ==HNRNPK[i,1] & ((HNRNPK[i,2] > ATAC_peak$start &HNRNPK[i,2] < ATAC_peak$end) | (HNRNPK[i,3] > ATAC_peak$start &HNRNPK[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HNRNPK) != 0){
    ATAC_peak_HNRNPK$V5 <-HNRNPK[i,1]
    ATAC_peak_HNRNPK$V6 <- HNRNPK[i,2]
    ATAC_peak_HNRNPK$V7 <- HNRNPK[i,3]
    ATAC_peak_HNRNPK$V8 <- HNRNPK[i,6]
    ATAC_peak_HNRNPK$V9 <- HNRNPK[i,7]
    ATAC_peak_HNRNPK_list = rbind(ATAC_peak_HNRNPK_list,ATAC_peak_HNRNPK)
  }
}
write.table(ATAC_peak_HNRNPK_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HNRNPK_list.txt",sep="\t")

ATAC_peak_HNRNPK_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HNRNPK_list.txt",sep="\t")

HNRNPK_features=list("HNRNPK" = ATAC_peak_HNRNPK_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HNRNPK_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##YBX1
YBX1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF513QUL_YBX1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(YBX1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF513QUL_YBX1.txt",sep = "\t")

YBX1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF513QUL_YBX1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(YBX1)){
  x=paste(YBX1[i,1],YBX1[i,2],YBX1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

YBX1_features = list(list$V1)

ATAC_peak_YBX1_list = data.frame()
for(i in 1:nrow(YBX1)){
  ATAC_peak_YBX1<- ATAC_peak[ATAC_peak$x1 ==YBX1[i,1] & ((YBX1[i,2] > ATAC_peak$start &YBX1[i,2] < ATAC_peak$end) | (YBX1[i,3] > ATAC_peak$start &YBX1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_YBX1) != 0){
    ATAC_peak_YBX1$V5 <-YBX1[i,1]
    ATAC_peak_YBX1$V6 <- YBX1[i,2]
    ATAC_peak_YBX1$V7 <- YBX1[i,3]
    ATAC_peak_YBX1$V8 <- YBX1[i,6]
    ATAC_peak_YBX1$V9 <- YBX1[i,7]
    ATAC_peak_YBX1_list = rbind(ATAC_peak_YBX1_list,ATAC_peak_YBX1)
  }
}
write.table(ATAC_peak_YBX1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_YBX1_list.txt",sep="\t")

ATAC_peak_YBX1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_YBX1_list.txt",sep="\t")

YBX1_features=list("YBX1" = ATAC_peak_YBX1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = YBX1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NCOR1
NCOR1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF518LOT_NCOR1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NCOR1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF518LOT_NCOR1_age.txt",sep = "\t")

NCOR1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF518LOT_NCOR1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NCOR1)){
  x=paste(NCOR1[i,1],NCOR1[i,2],NCOR1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NCOR1_features = list(list$V1)

ATAC_peak_NCOR1_list = data.frame()
for(i in 1:nrow(NCOR1)){
  ATAC_peak_NCOR1<- ATAC_peak[ATAC_peak$x1 ==NCOR1[i,1] & ((NCOR1[i,2] > ATAC_peak$start &NCOR1[i,2] < ATAC_peak$end) | (NCOR1[i,3] > ATAC_peak$start &NCOR1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NCOR1) != 0){
    ATAC_peak_NCOR1$V5 <-NCOR1[i,1]
    ATAC_peak_NCOR1$V6 <- NCOR1[i,2]
    ATAC_peak_NCOR1$V7 <- NCOR1[i,3]
    ATAC_peak_NCOR1$V8 <- NCOR1[i,6]
    ATAC_peak_NCOR1$V9 <- NCOR1[i,7]
    ATAC_peak_NCOR1_list = rbind(ATAC_peak_NCOR1_list,ATAC_peak_NCOR1)
  }
}
write.table(ATAC_peak_NCOR1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOR1_list.txt",sep="\t")

ATAC_peak_NCOR1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOR1_list.txt",sep="\t")

NCOR1_features=list("NCOR1" = ATAC_peak_NCOR1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NCOR1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZNF592
ZNF592 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF520ESS_ZNF592.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF592,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF520ESS_ZNF592.txt",sep = "\t")

ZNF592 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF520ESS_ZNF592.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF592)){
  x=paste(ZNF592[i,1],ZNF592[i,2],ZNF592[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF592_features = list(list$V1)

ATAC_peak_ZNF592_list = data.frame()
for(i in 1:nrow(ZNF592)){
  ATAC_peak_ZNF592<- ATAC_peak[ATAC_peak$x1 ==ZNF592[i,1] & ((ZNF592[i,2] > ATAC_peak$start &ZNF592[i,2] < ATAC_peak$end) | (ZNF592[i,3] > ATAC_peak$start &ZNF592[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF592) != 0){
    ATAC_peak_ZNF592$V5 <-ZNF592[i,1]
    ATAC_peak_ZNF592$V6 <- ZNF592[i,2]
    ATAC_peak_ZNF592$V7 <- ZNF592[i,3]
    ATAC_peak_ZNF592$V8 <- ZNF592[i,6]
    ATAC_peak_ZNF592$V9 <- ZNF592[i,7]
    ATAC_peak_ZNF592_list = rbind(ATAC_peak_ZNF592_list,ATAC_peak_ZNF592)
  }
}
write.table(ATAC_peak_ZNF592_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF592_list.txt",sep="\t")

ATAC_peak_ZNF592_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF592_list.txt",sep="\t")

ZNF592_features=list("ZNF592" = ATAC_peak_ZNF592_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF592_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CREBBP
CREBBP = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF532VPN_CREBBP_AGE.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CREBBP,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF532VPN_CREBBP_AGE.txt",sep = "\t")

CREBBP = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF532VPN_CREBBP_AGE.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CREBBP)){
  x=paste(CREBBP[i,1],CREBBP[i,2],CREBBP[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CREBBP_features = list(list$V1)

ATAC_peak_CREBBP_list = data.frame()
for(i in 1:nrow(CREBBP)){
  ATAC_peak_CREBBP<- ATAC_peak[ATAC_peak$x1 ==CREBBP[i,1] & ((CREBBP[i,2] > ATAC_peak$start &CREBBP[i,2] < ATAC_peak$end) | (CREBBP[i,3] > ATAC_peak$start &CREBBP[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CREBBP) != 0){
    ATAC_peak_CREBBP$V5 <-CREBBP[i,1]
    ATAC_peak_CREBBP$V6 <- CREBBP[i,2]
    ATAC_peak_CREBBP$V7 <- CREBBP[i,3]
    ATAC_peak_CREBBP$V8 <- CREBBP[i,6]
    ATAC_peak_CREBBP$V9 <- CREBBP[i,7]
    ATAC_peak_CREBBP_list = rbind(ATAC_peak_CREBBP_list,ATAC_peak_CREBBP)
  }
}
write.table(ATAC_peak_CREBBP_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CREBBP_list.txt",sep="\t")

ATAC_peak_CREBBP_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CREBBP_list.txt",sep="\t")

CREBBP_features=list("CREBBP" = ATAC_peak_CREBBP_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CREBBP_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZSCAN29
ZSCAN29 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF538JGD_ZSCAN29.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZSCAN29,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF538JGD_ZSCAN29.txt",sep = "\t")

ZSCAN29 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF538JGD_ZSCAN29.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZSCAN29)){
  x=paste(ZSCAN29[i,1],ZSCAN29[i,2],ZSCAN29[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZSCAN29_features = list(list$V1)

ATAC_peak_ZSCAN29_list = data.frame()
for(i in 1:nrow(ZSCAN29)){
  ATAC_peak_ZSCAN29<- ATAC_peak[ATAC_peak$x1 ==ZSCAN29[i,1] & ((ZSCAN29[i,2] > ATAC_peak$start &ZSCAN29[i,2] < ATAC_peak$end) | (ZSCAN29[i,3] > ATAC_peak$start &ZSCAN29[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZSCAN29) != 0){
    ATAC_peak_ZSCAN29$V5 <-ZSCAN29[i,1]
    ATAC_peak_ZSCAN29$V6 <- ZSCAN29[i,2]
    ATAC_peak_ZSCAN29$V7 <- ZSCAN29[i,3]
    ATAC_peak_ZSCAN29$V8 <- ZSCAN29[i,6]
    ATAC_peak_ZSCAN29$V9 <- ZSCAN29[i,7]
    ATAC_peak_ZSCAN29_list = rbind(ATAC_peak_ZSCAN29_list,ATAC_peak_ZSCAN29)
  }
}
write.table(ATAC_peak_ZSCAN29_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZSCAN29_list.txt",sep="\t")

ATAC_peak_ZSCAN29_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZSCAN29_list.txt",sep="\t")

ZSCAN29_features=list("ZSCAN29" = ATAC_peak_ZSCAN29_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZSCAN29_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##KDM4B
KDM4B = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF539PFU_KDM4B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(KDM4B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF539PFU_KDM4B.txt",sep = "\t")

KDM4B = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF539PFU_KDM4B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(KDM4B)){
  x=paste(KDM4B[i,1],KDM4B[i,2],KDM4B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

KDM4B_features = list(list$V1)

ATAC_peak_KDM4B_list = data.frame()
for(i in 1:nrow(KDM4B)){
  ATAC_peak_KDM4B<- ATAC_peak[ATAC_peak$x1 ==KDM4B[i,1] & ((KDM4B[i,2] > ATAC_peak$start &KDM4B[i,2] < ATAC_peak$end) | (KDM4B[i,3] > ATAC_peak$start &KDM4B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_KDM4B) != 0){
    ATAC_peak_KDM4B$V5 <-KDM4B[i,1]
    ATAC_peak_KDM4B$V6 <- KDM4B[i,2]
    ATAC_peak_KDM4B$V7 <- KDM4B[i,3]
    ATAC_peak_KDM4B$V8 <- KDM4B[i,6]
    ATAC_peak_KDM4B$V9 <- KDM4B[i,7]
    ATAC_peak_KDM4B_list = rbind(ATAC_peak_KDM4B_list,ATAC_peak_KDM4B)
  }
}
write.table(ATAC_peak_KDM4B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KDM4B_list.txt",sep="\t")

ATAC_peak_KDM4B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KDM4B_list.txt",sep="\t")

KDM4B_features=list("KDM4B" = ATAC_peak_KDM4B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = KDM4B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##EP300
EP300 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF539ZQW_EP300_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(EP300,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF539ZQW_EP300_age.txt",sep = "\t")

EP300 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF539ZQW_EP300_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(EP300)){
  x=paste(EP300[i,1],EP300[i,2],EP300[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

EP300_features = list(list$V1)

ATAC_peak_EP300_list = data.frame()
for(i in 1:nrow(EP300)){
  ATAC_peak_EP300<- ATAC_peak[ATAC_peak$x1 ==EP300[i,1] & ((EP300[i,2] > ATAC_peak$start &EP300[i,2] < ATAC_peak$end) | (EP300[i,3] > ATAC_peak$start &EP300[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_EP300) != 0){
    ATAC_peak_EP300$V5 <-EP300[i,1]
    ATAC_peak_EP300$V6 <- EP300[i,2]
    ATAC_peak_EP300$V7 <- EP300[i,3]
    ATAC_peak_EP300$V8 <- EP300[i,6]
    ATAC_peak_EP300$V9 <- EP300[i,7]
    ATAC_peak_EP300_list = rbind(ATAC_peak_EP300_list,ATAC_peak_EP300)
  }
}
write.table(ATAC_peak_EP300_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EP300_list.txt",sep="\t")

ATAC_peak_EP300_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EP300_list.txt",sep="\t")

EP300_features=list("EP300" = ATAC_peak_EP300_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = EP300_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##XRCC4
XRCC4 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF540QPM_XRCC4.bed.gz",extraCols=extraCols_narrowPeak)
write.table(XRCC4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/EENCFF540QPM_XRCC4.txt",sep = "\t")

XRCC4 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/EENCFF540QPM_XRCC4.txt",sep = "\t")

list = vector()
for (i in  1:nrow(XRCC4)){
  x=paste(XRCC4[i,1],XRCC4[i,2],XRCC4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

XRCC4_features = list(list$V1)

ATAC_peak_XRCC4_list = data.frame()
for(i in 1:nrow(XRCC4)){
  ATAC_peak_XRCC4<- ATAC_peak[ATAC_peak$x1 ==XRCC4[i,1] & ((XRCC4[i,2] > ATAC_peak$start &XRCC4[i,2] < ATAC_peak$end) | (XRCC4[i,3] > ATAC_peak$start &XRCC4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_XRCC4) != 0){
    ATAC_peak_XRCC4$V5 <-XRCC4[i,1]
    ATAC_peak_XRCC4$V6 <- XRCC4[i,2]
    ATAC_peak_XRCC4$V7 <- XRCC4[i,3]
    ATAC_peak_XRCC4$V8 <- XRCC4[i,6]
    ATAC_peak_XRCC4$V9 <- XRCC4[i,7]
    ATAC_peak_XRCC4_list = rbind(ATAC_peak_XRCC4_list,ATAC_peak_XRCC4)
  }
}
write.table(ATAC_peak_XRCC4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_XRCC4_list.txt",sep="\t")

ATAC_peak_XRCC4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_XRCC4_list.txt",sep="\t")

XRCC4_features=list("XRCC4" = ATAC_peak_XRCC4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = XRCC4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CBFA2T2
CBFA2T2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF543GET_CBFA2T2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CBFA2T2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF543GET_CBFA2T2.txt",sep = "\t")

CBFA2T2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF543GET_CBFA2T2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CBFA2T2)){
  x=paste(CBFA2T2[i,1],CBFA2T2[i,2],CBFA2T2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CBFA2T2_features = list(list$V1)

ATAC_peak_CBFA2T2_list = data.frame()
for(i in 1:nrow(CBFA2T2)){
  ATAC_peak_CBFA2T2<- ATAC_peak[ATAC_peak$x1 ==CBFA2T2[i,1] & ((CBFA2T2[i,2] > ATAC_peak$start &CBFA2T2[i,2] < ATAC_peak$end) | (CBFA2T2[i,3] > ATAC_peak$start &CBFA2T2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CBFA2T2) != 0){
    ATAC_peak_CBFA2T2$V5 <-CBFA2T2[i,1]
    ATAC_peak_CBFA2T2$V6 <- CBFA2T2[i,2]
    ATAC_peak_CBFA2T2$V7 <- CBFA2T2[i,3]
    ATAC_peak_CBFA2T2$V8 <- CBFA2T2[i,6]
    ATAC_peak_CBFA2T2$V9 <- CBFA2T2[i,7]
    ATAC_peak_CBFA2T2_list = rbind(ATAC_peak_CBFA2T2_list,ATAC_peak_CBFA2T2)
  }
}
write.table(ATAC_peak_CBFA2T2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBFA2T2_list.txt",sep="\t")

ATAC_peak_CBFA2T2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBFA2T2_list.txt",sep="\t")

CBFA2T2_features=list("CBFA2T2" = ATAC_peak_CBFA2T2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CBFA2T2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MTA3
MTA3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF553ZOP_MTA3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MTA3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF553ZOP_MTA3.txt",sep = "\t")

MTA3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF553ZOP_MTA3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MTA3)){
  x=paste(MTA3[i,1],MTA3[i,2],MTA3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MTA3_features = list(list$V1)

ATAC_peak_MTA3_list = data.frame()
for(i in 1:nrow(MTA3)){
  ATAC_peak_MTA3<- ATAC_peak[ATAC_peak$x1 ==MTA3[i,1] & ((MTA3[i,2] > ATAC_peak$start &MTA3[i,2] < ATAC_peak$end) | (MTA3[i,3] > ATAC_peak$start &MTA3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MTA3) != 0){
    ATAC_peak_MTA3$V5 <-MTA3[i,1]
    ATAC_peak_MTA3$V6 <- MTA3[i,2]
    ATAC_peak_MTA3$V7 <- MTA3[i,3]
    ATAC_peak_MTA3$V8 <- MTA3[i,6]
    ATAC_peak_MTA3$V9 <- MTA3[i,7]
    ATAC_peak_MTA3_list = rbind(ATAC_peak_MTA3_list,ATAC_peak_MTA3)
  }
}
write.table(ATAC_peak_MTA3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MTA3_list.txt",sep="\t")

ATAC_peak_MTA3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MTA3_list.txt",sep="\t")

MTA3_features=list("MTA3" = ATAC_peak_MTA3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MTA3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TRIP13
TRIP13 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF558MHP_TRIP13.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TRIP13,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/NCFF558MHP_TRIP13.txt",sep = "\t")

TRIP13 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/NCFF558MHP_TRIP13.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TRIP13)){
  x=paste(TRIP13[i,1],TRIP13[i,2],TRIP13[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TRIP13_features = list(list$V1)

ATAC_peak_TRIP13_list = data.frame()
for(i in 1:nrow(TRIP13)){
  ATAC_peak_TRIP13<- ATAC_peak[ATAC_peak$x1 ==TRIP13[i,1] & ((TRIP13[i,2] > ATAC_peak$start &TRIP13[i,2] < ATAC_peak$end) | (TRIP13[i,3] > ATAC_peak$start &TRIP13[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TRIP13) != 0){
    ATAC_peak_TRIP13$V5 <-TRIP13[i,1]
    ATAC_peak_TRIP13$V6 <- TRIP13[i,2]
    ATAC_peak_TRIP13$V7 <- TRIP13[i,3]
    ATAC_peak_TRIP13$V8 <- TRIP13[i,6]
    ATAC_peak_TRIP13$V9 <- TRIP13[i,7]
    ATAC_peak_TRIP13_list = rbind(ATAC_peak_TRIP13_list,ATAC_peak_TRIP13)
  }
}
write.table(ATAC_peak_TRIP13_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TRIP13_list.txt",sep="\t")

ATAC_peak_TRIP13_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TRIP13_list.txt",sep="\t")

TRIP13_features=list("TRIP13" = ATAC_peak_TRIP13_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TRIP13_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ATF2
ATF2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF559LTZ_ATF2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ATF2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF559LTZ_ATF2.txt",sep = "\t")

ATF2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF559LTZ_ATF2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ATF2)){
  x=paste(ATF2[i,1],ATF2[i,2],ATF2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ATF2_features = list(list$V1)

ATAC_peak_ATF2_list = data.frame()
for(i in 1:nrow(ATF2)){
  ATAC_peak_ATF2<- ATAC_peak[ATAC_peak$x1 ==ATF2[i,1] & ((ATF2[i,2] > ATAC_peak$start &ATF2[i,2] < ATAC_peak$end) | (ATF2[i,3] > ATAC_peak$start &ATF2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ATF2) != 0){
    ATAC_peak_ATF2$V5 <-ATF2[i,1]
    ATAC_peak_ATF2$V6 <- ATF2[i,2]
    ATAC_peak_ATF2$V7 <- ATF2[i,3]
    ATAC_peak_ATF2$V8 <- ATF2[i,6]
    ATAC_peak_ATF2$V9 <- ATF2[i,7]
    ATAC_peak_ATF2_list = rbind(ATAC_peak_ATF2_list,ATAC_peak_ATF2)
  }
}
write.table(ATAC_peak_ATF2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATF2_list.txt",sep="\t")

ATAC_peak_ATF2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATF2_list.txt",sep="\t")

ATF2_features=list("ATF2" = ATAC_peak_ATF2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ATF2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ARID3A
ARID3A = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF563SWF_ARID3A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ARID3A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF563SWF_ARID3A.txt",sep = "\t")

ARID3A = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF563SWF_ARID3A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ARID3A)){
  x=paste(ARID3A[i,1],ARID3A[i,2],ARID3A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ARID3A_features = list(list$V1)

ATAC_peak_ARID3A_list = data.frame()
for(i in 1:nrow(ARID3A)){
  ATAC_peak_ARID3A<- ATAC_peak[ATAC_peak$x1 ==ARID3A[i,1] & ((ARID3A[i,2] > ATAC_peak$start &ARID3A[i,2] < ATAC_peak$end) | (ARID3A[i,3] > ATAC_peak$start &ARID3A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ARID3A) != 0){
    ATAC_peak_ARID3A$V5 <-ARID3A[i,1]
    ATAC_peak_ARID3A$V6 <- ARID3A[i,2]
    ATAC_peak_ARID3A$V7 <- ARID3A[i,3]
    ATAC_peak_ARID3A$V8 <- ARID3A[i,6]
    ATAC_peak_ARID3A$V9 <- ARID3A[i,7]
    ATAC_peak_ARID3A_list = rbind(ATAC_peak_ARID3A_list,ATAC_peak_ARID3A)
  }
}
write.table(ATAC_peak_ARID3A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARID3A_list.txt",sep="\t")

ATAC_peak_ARID3A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARID3A_list.txt",sep="\t")

ARID3A_features=list("ARID3A" = ATAC_peak_ARID3A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ARID3A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##NFRKB
NFRKB = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF568ETE_NFRKB.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NFRKB,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF568ETE_NFRKB.txt",sep = "\t")

NFRKB = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF568ETE_NFRKB.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NFRKB)){
  x=paste(NFRKB[i,1],NFRKB[i,2],NFRKB[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NFRKB_features = list(list$V1)

ATAC_peak_NFRKB_list = data.frame()
for(i in 1:nrow(NFRKB)){
  ATAC_peak_NFRKB<- ATAC_peak[ATAC_peak$x1 ==NFRKB[i,1] & ((NFRKB[i,2] > ATAC_peak$start &NFRKB[i,2] < ATAC_peak$end) | (NFRKB[i,3] > ATAC_peak$start &NFRKB[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NFRKB) != 0){
    ATAC_peak_NFRKB$V5 <-NFRKB[i,1]
    ATAC_peak_NFRKB$V6 <- NFRKB[i,2]
    ATAC_peak_NFRKB$V7 <- NFRKB[i,3]
    ATAC_peak_NFRKB$V8 <- NFRKB[i,6]
    ATAC_peak_NFRKB$V9 <- NFRKB[i,7]
    ATAC_peak_NFRKB_list = rbind(ATAC_peak_NFRKB_list,ATAC_peak_NFRKB)
  }
}
write.table(ATAC_peak_NFRKB_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFRKB_list.txt",sep="\t")

ATAC_peak_NFRKB_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFRKB_list.txt",sep="\t")

NFRKB_features=list("NFRKB" = ATAC_peak_NFRKB_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NFRKB_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZFP36
ZFP36 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF574TFC_ZFP36.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZFP36,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF574TFC_ZFP36.txt",sep = "\t")

ZFP36 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF574TFC_ZFP36.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZFP36)){
  x=paste(ZFP36[i,1],ZFP36[i,2],ZFP36[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZFP36_features = list(list$V1)

ATAC_peak_ZFP36_list = data.frame()
for(i in 1:nrow(ZFP36)){
  ATAC_peak_ZFP36<- ATAC_peak[ATAC_peak$x1 ==ZFP36[i,1] & ((ZFP36[i,2] > ATAC_peak$start &ZFP36[i,2] < ATAC_peak$end) | (ZFP36[i,3] > ATAC_peak$start &ZFP36[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZFP36) != 0){
    ATAC_peak_ZFP36$V5 <-ZFP36[i,1]
    ATAC_peak_ZFP36$V6 <- ZFP36[i,2]
    ATAC_peak_ZFP36$V7 <- ZFP36[i,3]
    ATAC_peak_ZFP36$V8 <- ZFP36[i,6]
    ATAC_peak_ZFP36$V9 <- ZFP36[i,7]
    ATAC_peak_ZFP36_list = rbind(ATAC_peak_ZFP36_list,ATAC_peak_ZFP36)
  }
}
write.table(ATAC_peak_ZFP36_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZFP36_list.txt",sep="\t")

ATAC_peak_ZFP36_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZFP36_list.txt",sep="\t")

ZFP36_features=list("ZFP36" = ATAC_peak_ZFP36_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZFP36_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CREB1
CREB1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF576BKY_CREB1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CREB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF576BKY_CREB1_age.txt",sep = "\t")

CREB1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF576BKY_CREB1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CREB1)){
  x=paste(CREB1[i,1],CREB1[i,2],CREB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CREB1_features = list(list$V1)

ATAC_peak_CREB1_list = data.frame()
for(i in 1:nrow(CREB1)){
  ATAC_peak_CREB1<- ATAC_peak[ATAC_peak$x1 ==CREB1[i,1] & ((CREB1[i,2] > ATAC_peak$start &CREB1[i,2] < ATAC_peak$end) | (CREB1[i,3] > ATAC_peak$start &CREB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CREB1) != 0){
    ATAC_peak_CREB1$V5 <-CREB1[i,1]
    ATAC_peak_CREB1$V6 <- CREB1[i,2]
    ATAC_peak_CREB1$V7 <- CREB1[i,3]
    ATAC_peak_CREB1$V8 <- CREB1[i,6]
    ATAC_peak_CREB1$V9 <- CREB1[i,7]
    ATAC_peak_CREB1_list = rbind(ATAC_peak_CREB1_list,ATAC_peak_CREB1)
  }
}
write.table(ATAC_peak_CREB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CREB1_list.txt",sep="\t")

ATAC_peak_CREB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CREB1_list.txt",sep="\t")

CREB1_features=list("CREB1" = ATAC_peak_CREB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CREB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##E2F1
E2F1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF576KKB_E2F1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(E2F1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF576KKB_E2F1_age.txt",sep = "\t")

E2F1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF576KKB_E2F1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(E2F1)){
  x=paste(E2F1[i,1],E2F1[i,2],E2F1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

E2F1_features = list(list$V1)

ATAC_peak_E2F1_list = data.frame()
for(i in 1:nrow(E2F1)){
  ATAC_peak_E2F1<- ATAC_peak[ATAC_peak$x1 ==E2F1[i,1] & ((E2F1[i,2] > ATAC_peak$start &E2F1[i,2] < ATAC_peak$end) | (E2F1[i,3] > ATAC_peak$start &E2F1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_E2F1) != 0){
    ATAC_peak_E2F1$V5 <-E2F1[i,1]
    ATAC_peak_E2F1$V6 <- E2F1[i,2]
    ATAC_peak_E2F1$V7 <- E2F1[i,3]
    ATAC_peak_E2F1$V8 <- E2F1[i,6]
    ATAC_peak_E2F1$V9 <- E2F1[i,7]
    ATAC_peak_E2F1_list = rbind(ATAC_peak_E2F1_list,ATAC_peak_E2F1)
  }
}
write.table(ATAC_peak_E2F1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_E2F1_list.txt",sep="\t")

ATAC_peak_E2F1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_E2F1_list.txt",sep="\t")

E2F1_features=list("E2F1" = ATAC_peak_E2F1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = E2F1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##XRCC3
XRCC3 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF577FBP_XRCC3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(XRCC3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF577FBP_XRCC3.txt",sep = "\t")

XRCC3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF577FBP_XRCC3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(XRCC3)){
  x=paste(XRCC3[i,1],XRCC3[i,2],XRCC3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

XRCC3_features = list(list$V1)

ATAC_peak_XRCC3_list = data.frame()
for(i in 1:nrow(XRCC3)){
  ATAC_peak_XRCC3<- ATAC_peak[ATAC_peak$x1 ==XRCC3[i,1] & ((XRCC3[i,2] > ATAC_peak$start &XRCC3[i,2] < ATAC_peak$end) | (XRCC3[i,3] > ATAC_peak$start &XRCC3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_XRCC3) != 0){
    ATAC_peak_XRCC3$V5 <-XRCC3[i,1]
    ATAC_peak_XRCC3$V6 <- XRCC3[i,2]
    ATAC_peak_XRCC3$V7 <- XRCC3[i,3]
    ATAC_peak_XRCC3$V8 <- XRCC3[i,6]
    ATAC_peak_XRCC3$V9 <- XRCC3[i,7]
    ATAC_peak_XRCC3_list = rbind(ATAC_peak_XRCC3_list,ATAC_peak_XRCC3)
  }
}
write.table(ATAC_peak_XRCC3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_XRCC3_list.txt",sep="\t")

ATAC_peak_XRCC3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_XRCC3_list.txt",sep="\t")

XRCC3_features=list("XRCC3" = ATAC_peak_XRCC3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = XRCC3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##GTF2B
GTF2B = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF579OUH_GTF2B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GTF2B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF579OUH_GTF2B.txt",sep = "\t")

GTF2B = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF579OUH_GTF2B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GTF2B)){
  x=paste(GTF2B[i,1],GTF2B[i,2],GTF2B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GTF2B_features = list(list$V1)

ATAC_peak_GTF2B_list = data.frame()
for(i in 1:nrow(GTF2B)){
  ATAC_peak_GTF2B<- ATAC_peak[ATAC_peak$x1 ==GTF2B[i,1] & ((GTF2B[i,2] > ATAC_peak$start &GTF2B[i,2] < ATAC_peak$end) | (GTF2B[i,3] > ATAC_peak$start &GTF2B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GTF2B) != 0){
    ATAC_peak_GTF2B$V5 <-GTF2B[i,1]
    ATAC_peak_GTF2B$V6 <- GTF2B[i,2]
    ATAC_peak_GTF2B$V7 <- GTF2B[i,3]
    ATAC_peak_GTF2B$V8 <- GTF2B[i,6]
    ATAC_peak_GTF2B$V9 <- GTF2B[i,7]
    ATAC_peak_GTF2B_list = rbind(ATAC_peak_GTF2B_list,ATAC_peak_GTF2B)
  }
}
write.table(ATAC_peak_GTF2B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF2B_list.txt",sep="\t")

ATAC_peak_GTF2B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF2B_list.txt",sep="\t")

GTF2B_features=list("GTF2B" = ATAC_peak_GTF2B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GTF2B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##FUS
FUS = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF581MLB_FUS.bed.gz",extraCols=extraCols_narrowPeak)
write.table(FUS,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF581MLB_FUS.txt",sep = "\t")

FUS = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF581MLB_FUS.txt",sep = "\t")

list = vector()
for (i in  1:nrow(FUS)){
  x=paste(FUS[i,1],FUS[i,2],FUS[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

FUS_features = list(list$V1)

ATAC_peak_FUS_list = data.frame()
for(i in 1:nrow(FUS)){
  ATAC_peak_FUS<- ATAC_peak[ATAC_peak$x1 ==FUS[i,1] & ((FUS[i,2] > ATAC_peak$start &FUS[i,2] < ATAC_peak$end) | (FUS[i,3] > ATAC_peak$start &FUS[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_FUS) != 0){
    ATAC_peak_FUS$V5 <-FUS[i,1]
    ATAC_peak_FUS$V6 <- FUS[i,2]
    ATAC_peak_FUS$V7 <- FUS[i,3]
    ATAC_peak_FUS$V8 <- FUS[i,6]
    ATAC_peak_FUS$V9 <- FUS[i,7]
    ATAC_peak_FUS_list = rbind(ATAC_peak_FUS_list,ATAC_peak_FUS)
  }
}
write.table(ATAC_peak_FUS_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FUS_list.txt",sep="\t")

ATAC_peak_FUS_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FUS_list.txt",sep="\t")

FUS_features=list("FUS" = ATAC_peak_FUS_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = FUS_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##KLF1
KLF1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF583PGH_KLF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(KLF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF583PGH_KLF1.txt",sep = "\t")

KLF1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF583PGH_KLF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(KLF1)){
  x=paste(KLF1[i,1],KLF1[i,2],KLF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

KLF1_features = list(list$V1)

ATAC_peak_KLF1_list = data.frame()
for(i in 1:nrow(KLF1)){
  ATAC_peak_KLF1<- ATAC_peak[ATAC_peak$x1 ==KLF1[i,1] & ((KLF1[i,2] > ATAC_peak$start &KLF1[i,2] < ATAC_peak$end) | (KLF1[i,3] > ATAC_peak$start &KLF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_KLF1) != 0){
    ATAC_peak_KLF1$V5 <-KLF1[i,1]
    ATAC_peak_KLF1$V6 <- KLF1[i,2]
    ATAC_peak_KLF1$V7 <- KLF1[i,3]
    ATAC_peak_KLF1$V8 <- KLF1[i,6]
    ATAC_peak_KLF1$V9 <- KLF1[i,7]
    ATAC_peak_KLF1_list = rbind(ATAC_peak_KLF1_list,ATAC_peak_KLF1)
  }
}
write.table(ATAC_peak_KLF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KLF1_list.txt",sep="\t")

ATAC_peak_KLF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KLF1_list.txt",sep="\t")

KLF1_features=list("KLF1" = ATAC_peak_KLF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = KLF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TAF9B
TAF9B = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF585DTZ_TAF9B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TAF9B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF585DTZ_TAF9B.txt",sep = "\t")

TAF9B = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF585DTZ_TAF9B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TAF9B)){
  x=paste(TAF9B[i,1],TAF9B[i,2],TAF9B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TAF9B_features = list(list$V1)

ATAC_peak_TAF9B_list = data.frame()
for(i in 1:nrow(TAF9B)){
  ATAC_peak_TAF9B<- ATAC_peak[ATAC_peak$x1 ==TAF9B[i,1] & ((TAF9B[i,2] > ATAC_peak$start &TAF9B[i,2] < ATAC_peak$end) | (TAF9B[i,3] > ATAC_peak$start &TAF9B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TAF9B) != 0){
    ATAC_peak_TAF9B$V5 <-TAF9B[i,1]
    ATAC_peak_TAF9B$V6 <- TAF9B[i,2]
    ATAC_peak_TAF9B$V7 <- TAF9B[i,3]
    ATAC_peak_TAF9B$V8 <- TAF9B[i,6]
    ATAC_peak_TAF9B$V9 <- TAF9B[i,7]
    ATAC_peak_TAF9B_list = rbind(ATAC_peak_TAF9B_list,ATAC_peak_TAF9B)
  }
}
write.table(ATAC_peak_TAF9B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAF9B_list.txt",sep="\t")

ATAC_peak_TAF9B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAF9B_list.txt",sep="\t")

TAF9B_features=list("TAF9B" = ATAC_peak_TAF9B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TAF9B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TAF1
TAF1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF598CIP_TAF1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TAF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF598CIP_TAF1_age.txt",sep = "\t")

TAF1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF598CIP_TAF1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TAF1)){
  x=paste(TAF1[i,1],TAF1[i,2],TAF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TAF1_features = list(list$V1)

ATAC_peak_TAF1_list = data.frame()
for(i in 1:nrow(TAF1)){
  ATAC_peak_TAF1<- ATAC_peak[ATAC_peak$x1 ==TAF1[i,1] & ((TAF1[i,2] > ATAC_peak$start &TAF1[i,2] < ATAC_peak$end) | (TAF1[i,3] > ATAC_peak$start &TAF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TAF1) != 0){
    ATAC_peak_TAF1$V5 <-TAF1[i,1]
    ATAC_peak_TAF1$V6 <- TAF1[i,2]
    ATAC_peak_TAF1$V7 <- TAF1[i,3]
    ATAC_peak_TAF1$V8 <- TAF1[i,6]
    ATAC_peak_TAF1$V9 <- TAF1[i,7]
    ATAC_peak_TAF1_list = rbind(ATAC_peak_TAF1_list,ATAC_peak_TAF1)
  }
}
write.table(ATAC_peak_TAF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAF1_list.txt",sep="\t")

ATAC_peak_TAF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAF1_list.txt",sep="\t")

TAF1_features=list("TAF1" = ATAC_peak_TAF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TAF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##DDX20
DDX20 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF602OHR_DDX20.bed.gz",extraCols=extraCols_narrowPeak)
write.table(DDX20,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF602OHR_DDX20.txt",sep = "\t")

DDX20 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF602OHR_DDX20.txt",sep = "\t")

list = vector()
for (i in  1:nrow(DDX20)){
  x=paste(DDX20[i,1],DDX20[i,2],DDX20[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

DDX20_features = list(list$V1)

ATAC_peak_DDX20_list = data.frame()
for(i in 1:nrow(DDX20)){
  ATAC_peak_DDX20<- ATAC_peak[ATAC_peak$x1 ==DDX20[i,1] & ((DDX20[i,2] > ATAC_peak$start &DDX20[i,2] < ATAC_peak$end) | (DDX20[i,3] > ATAC_peak$start &DDX20[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_DDX20) != 0){
    ATAC_peak_DDX20$V5 <-DDX20[i,1]
    ATAC_peak_DDX20$V6 <- DDX20[i,2]
    ATAC_peak_DDX20$V7 <- DDX20[i,3]
    ATAC_peak_DDX20$V8 <- DDX20[i,6]
    ATAC_peak_DDX20$V9 <- DDX20[i,7]
    ATAC_peak_DDX20_list = rbind(ATAC_peak_DDX20_list,ATAC_peak_DDX20)
  }
}
write.table(ATAC_peak_DDX20_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DDX20_list.txt",sep="\t")

ATAC_peak_DDX20_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DDX20_list.txt",sep="\t")

DDX20_features=list("DDX20" = ATAC_peak_DDX20_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = DDX20_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##GMEB1
GMEB1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF606HUJ_GMEB1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GMEB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF606HUJ_GMEB1.txt",sep = "\t")

GMEB1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF606HUJ_GMEB1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GMEB1)){
  x=paste(GMEB1[i,1],GMEB1[i,2],GMEB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GMEB1_features = list(list$V1)

ATAC_peak_GMEB1_list = data.frame()
for(i in 1:nrow(GMEB1)){
  ATAC_peak_GMEB1<- ATAC_peak[ATAC_peak$x1 ==GMEB1[i,1] & ((GMEB1[i,2] > ATAC_peak$start &GMEB1[i,2] < ATAC_peak$end) | (GMEB1[i,3] > ATAC_peak$start &GMEB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GMEB1) != 0){
    ATAC_peak_GMEB1$V5 <-GMEB1[i,1]
    ATAC_peak_GMEB1$V6 <- GMEB1[i,2]
    ATAC_peak_GMEB1$V7 <- GMEB1[i,3]
    ATAC_peak_GMEB1$V8 <- GMEB1[i,6]
    ATAC_peak_GMEB1$V9 <- GMEB1[i,7]
    ATAC_peak_GMEB1_list = rbind(ATAC_peak_GMEB1_list,ATAC_peak_GMEB1)
  }
}
write.table(ATAC_peak_GMEB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GMEB1_list.txt",sep="\t")

ATAC_peak_GMEB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GMEB1_list.txt",sep="\t")

GMEB1_features=list("GMEB1" = ATAC_peak_GMEB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GMEB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##TBP
TBP = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF609UTS_TBP.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TBP,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF609UTS_TBP.txt",sep = "\t")

TBP = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF609UTS_TBP.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TBP)){
  x=paste(TBP[i,1],TBP[i,2],TBP[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TBP_features = list(list$V1)

ATAC_peak_TBP_list = data.frame()
for(i in 1:nrow(TBP)){
  ATAC_peak_TBP<- ATAC_peak[ATAC_peak$x1 ==TBP[i,1] & ((TBP[i,2] > ATAC_peak$start &TBP[i,2] < ATAC_peak$end) | (TBP[i,3] > ATAC_peak$start &TBP[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TBP) != 0){
    ATAC_peak_TBP$V5 <-TBP[i,1]
    ATAC_peak_TBP$V6 <- TBP[i,2]
    ATAC_peak_TBP$V7 <- TBP[i,3]
    ATAC_peak_TBP$V8 <- TBP[i,6]
    ATAC_peak_TBP$V9 <- TBP[i,7]
    ATAC_peak_TBP_list = rbind(ATAC_peak_TBP_list,ATAC_peak_TBP)
  }
}
write.table(ATAC_peak_TBP_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TBP_list.txt",sep="\t")

ATAC_peak_TBP_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TBP_list.txt",sep="\t")

TBP_features=list("TBP" = ATAC_peak_TBP_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TBP_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SUPT5H
SUPT5H = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF612PGJ_SUPT5H.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SUPT5H,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF612PGJ_SUPT5H.txt",sep = "\t")

SUPT5H = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF612PGJ_SUPT5H.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SUPT5H)){
  x=paste(SUPT5H[i,1],SUPT5H[i,2],SUPT5H[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SUPT5H_features = list(list$V1)

ATAC_peak_SUPT5H_list = data.frame()
for(i in 1:nrow(SUPT5H)){
  ATAC_peak_SUPT5H<- ATAC_peak[ATAC_peak$x1 ==SUPT5H[i,1] & ((SUPT5H[i,2] > ATAC_peak$start &SUPT5H[i,2] < ATAC_peak$end) | (SUPT5H[i,3] > ATAC_peak$start &SUPT5H[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SUPT5H) != 0){
    ATAC_peak_SUPT5H$V5 <-SUPT5H[i,1]
    ATAC_peak_SUPT5H$V6 <- SUPT5H[i,2]
    ATAC_peak_SUPT5H$V7 <- SUPT5H[i,3]
    ATAC_peak_SUPT5H$V8 <- SUPT5H[i,6]
    ATAC_peak_SUPT5H$V9 <- SUPT5H[i,7]
    ATAC_peak_SUPT5H_list = rbind(ATAC_peak_SUPT5H_list,ATAC_peak_SUPT5H)
  }
}
write.table(ATAC_peak_SUPT5H_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SUPT5H_list.txt",sep="\t")

ATAC_peak_SUPT5H_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SUPT5H_list.txt",sep="\t")

SUPT5H_features=list("SUPT5H" = ATAC_peak_SUPT5H_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SUPT5H_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##HBP1
HBP1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF617ZLR_HBP1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HBP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF617ZLR_HBP1_age.txt",sep = "\t")

HBP1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF617ZLR_HBP1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HBP1)){
  x=paste(HBP1[i,1],HBP1[i,2],HBP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HBP1_features = list(list$V1)

ATAC_peak_HBP1_list = data.frame()
for(i in 1:nrow(HBP1)){
  ATAC_peak_HBP1<- ATAC_peak[ATAC_peak$x1 ==HBP1[i,1] & ((HBP1[i,2] > ATAC_peak$start &HBP1[i,2] < ATAC_peak$end) | (HBP1[i,3] > ATAC_peak$start &HBP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HBP1) != 0){
    ATAC_peak_HBP1$V5 <-HBP1[i,1]
    ATAC_peak_HBP1$V6 <- HBP1[i,2]
    ATAC_peak_HBP1$V7 <- HBP1[i,3]
    ATAC_peak_HBP1$V8 <- HBP1[i,6]
    ATAC_peak_HBP1$V9 <- HBP1[i,7]
    ATAC_peak_HBP1_list = rbind(ATAC_peak_HBP1_list,ATAC_peak_HBP1)
  }
}
write.table(ATAC_peak_HBP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HBP1_list.txt",sep="\t")

ATAC_peak_HBP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HBP1_list.txt",sep="\t")

HBP1_features=list("HBP1" = ATAC_peak_HBP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HBP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CBX1
CBX1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF622IAZ_CBX1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CBX1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF622IAZ_CBX1.txt",sep = "\t")

CBX1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF622IAZ_CBX1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CBX1)){
  x=paste(CBX1[i,1],CBX1[i,2],CBX1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CBX1_features = list(list$V1)

ATAC_peak_CBX1_list = data.frame()
for(i in 1:nrow(CBX1)){
  ATAC_peak_CBX1<- ATAC_peak[ATAC_peak$x1 ==CBX1[i,1] & ((CBX1[i,2] > ATAC_peak$start &CBX1[i,2] < ATAC_peak$end) | (CBX1[i,3] > ATAC_peak$start &CBX1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CBX1) != 0){
    ATAC_peak_CBX1$V5 <-CBX1[i,1]
    ATAC_peak_CBX1$V6 <- CBX1[i,2]
    ATAC_peak_CBX1$V7 <- CBX1[i,3]
    ATAC_peak_CBX1$V8 <- CBX1[i,6]
    ATAC_peak_CBX1$V9 <- CBX1[i,7]
    ATAC_peak_CBX1_list = rbind(ATAC_peak_CBX1_list,ATAC_peak_CBX1)
  }
}
write.table(ATAC_peak_CBX1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX1_list.txt",sep="\t")

ATAC_peak_CBX1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX1_list.txt",sep="\t")

CBX1_features=list("CBX1" = ATAC_peak_CBX1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CBX1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CHAMP1
CHAMP1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF626FDW_CHAMP1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CHAMP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF626FDW_CHAMP1.txt",sep = "\t")

CHAMP1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF626FDW_CHAMP1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CHAMP1)){
  x=paste(CHAMP1[i,1],CHAMP1[i,2],CHAMP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CHAMP1_features = list(list$V1)

ATAC_peak_CHAMP1_list = data.frame()
for(i in 1:nrow(CHAMP1)){
  ATAC_peak_CHAMP1<- ATAC_peak[ATAC_peak$x1 ==CHAMP1[i,1] & ((CHAMP1[i,2] > ATAC_peak$start &CHAMP1[i,2] < ATAC_peak$end) | (CHAMP1[i,3] > ATAC_peak$start &CHAMP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CHAMP1) != 0){
    ATAC_peak_CHAMP1$V5 <-CHAMP1[i,1]
    ATAC_peak_CHAMP1$V6 <- CHAMP1[i,2]
    ATAC_peak_CHAMP1$V7 <- CHAMP1[i,3]
    ATAC_peak_CHAMP1$V8 <- CHAMP1[i,6]
    ATAC_peak_CHAMP1$V9 <- CHAMP1[i,7]
    ATAC_peak_CHAMP1_list = rbind(ATAC_peak_CHAMP1_list,ATAC_peak_CHAMP1)
  }
}
write.table(ATAC_peak_CHAMP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHAMP1_list.txt",sep="\t")

ATAC_peak_CHAMP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHAMP1_list.txt",sep="\t")

CHAMP1_features=list("CHAMP1" = ATAC_peak_CHAMP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CHAMP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##GTF2A2
GTF2A2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF632HFK_GTF2A2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GTF2A2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF632HFK_GTF2A2.txt",sep = "\t")

GTF2A2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF632HFK_GTF2A2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GTF2A2)){
  x=paste(GTF2A2[i,1],GTF2A2[i,2],GTF2A2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GTF2A2_features = list(list$V1)

ATAC_peak_GTF2A2_list = data.frame()
for(i in 1:nrow(GTF2A2)){
  ATAC_peak_GTF2A2<- ATAC_peak[ATAC_peak$x1 ==GTF2A2[i,1] & ((GTF2A2[i,2] > ATAC_peak$start &GTF2A2[i,2] < ATAC_peak$end) | (GTF2A2[i,3] > ATAC_peak$start &GTF2A2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GTF2A2) != 0){
    ATAC_peak_GTF2A2$V5 <-GTF2A2[i,1]
    ATAC_peak_GTF2A2$V6 <- GTF2A2[i,2]
    ATAC_peak_GTF2A2$V7 <- GTF2A2[i,3]
    ATAC_peak_GTF2A2$V8 <- GTF2A2[i,6]
    ATAC_peak_GTF2A2$V9 <- GTF2A2[i,7]
    ATAC_peak_GTF2A2_list = rbind(ATAC_peak_GTF2A2_list,ATAC_peak_GTF2A2)
  }
}
write.table(ATAC_peak_GTF2A2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF2A2_list.txt",sep="\t")

ATAC_peak_GTF2A2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF2A2_list.txt",sep="\t")

GTF2A2_features=list("GTF2A2" = ATAC_peak_GTF2A2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GTF2A2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##BRD4
BRD4 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF632ZRY_BRD4.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BRD4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF632ZRY_BRD4.txt",sep = "\t")

BRD4 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF632ZRY_BRD4.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BRD4)){
  x=paste(BRD4[i,1],BRD4[i,2],BRD4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BRD4_features = list(list$V1)

ATAC_peak_BRD4_list = data.frame()
for(i in 1:nrow(BRD4)){
  ATAC_peak_BRD4<- ATAC_peak[ATAC_peak$x1 ==BRD4[i,1] & ((BRD4[i,2] > ATAC_peak$start &BRD4[i,2] < ATAC_peak$end) | (BRD4[i,3] > ATAC_peak$start &BRD4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BRD4) != 0){
    ATAC_peak_BRD4$V5 <-BRD4[i,1]
    ATAC_peak_BRD4$V6 <- BRD4[i,2]
    ATAC_peak_BRD4$V7 <- BRD4[i,3]
    ATAC_peak_BRD4$V8 <- BRD4[i,6]
    ATAC_peak_BRD4$V9 <- BRD4[i,7]
    ATAC_peak_BRD4_list = rbind(ATAC_peak_BRD4_list,ATAC_peak_BRD4)
  }
}
write.table(ATAC_peak_BRD4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRD4_list.txt",sep="\t")

ATAC_peak_BRD4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRD4_list.txt",sep="\t")

BRD4_features=list("BRD4" = ATAC_peak_BRD4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BRD4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NCOA2
NCOA2 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF635KPD_NCOA2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NCOA2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/NCFF635KPD_NCOA2.txt",sep = "\t")

NCOA2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/NCFF635KPD_NCOA2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NCOA2)){
  x=paste(NCOA2[i,1],NCOA2[i,2],NCOA2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NCOA2_features = list(list$V1)

ATAC_peak_NCOA2_list = data.frame()
for(i in 1:nrow(NCOA2)){
  ATAC_peak_NCOA2<- ATAC_peak[ATAC_peak$x1 ==NCOA2[i,1] & ((NCOA2[i,2] > ATAC_peak$start &NCOA2[i,2] < ATAC_peak$end) | (NCOA2[i,3] > ATAC_peak$start &NCOA2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NCOA2) != 0){
    ATAC_peak_NCOA2$V5 <-NCOA2[i,1]
    ATAC_peak_NCOA2$V6 <- NCOA2[i,2]
    ATAC_peak_NCOA2$V7 <- NCOA2[i,3]
    ATAC_peak_NCOA2$V8 <- NCOA2[i,6]
    ATAC_peak_NCOA2$V9 <- NCOA2[i,7]
    ATAC_peak_NCOA2_list = rbind(ATAC_peak_NCOA2_list,ATAC_peak_NCOA2)
  }
}
write.table(ATAC_peak_NCOA2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOA2_list.txt",sep="\t")

ATAC_peak_NCOA2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOA2_list.txt",sep="\t")

NCOA2_features=list("NCOA2" = ATAC_peak_NCOA2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NCOA2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##HDAC1
HDAC1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF638PDM_HDAC1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HDAC1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF638PDM_HDAC1_age.txt",sep = "\t")

HDAC1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF638PDM_HDAC1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HDAC1)){
  x=paste(HDAC1[i,1],HDAC1[i,2],HDAC1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HDAC1_features = list(list$V1)

ATAC_peak_HDAC1_list = data.frame()
for(i in 1:nrow(HDAC1)){
  ATAC_peak_HDAC1<- ATAC_peak[ATAC_peak$x1 ==HDAC1[i,1] & ((HDAC1[i,2] > ATAC_peak$start &HDAC1[i,2] < ATAC_peak$end) | (HDAC1[i,3] > ATAC_peak$start &HDAC1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HDAC1) != 0){
    ATAC_peak_HDAC1$V5 <-HDAC1[i,1]
    ATAC_peak_HDAC1$V6 <- HDAC1[i,2]
    ATAC_peak_HDAC1$V7 <- HDAC1[i,3]
    ATAC_peak_HDAC1$V8 <- HDAC1[i,6]
    ATAC_peak_HDAC1$V9 <- HDAC1[i,7]
    ATAC_peak_HDAC1_list = rbind(ATAC_peak_HDAC1_list,ATAC_peak_HDAC1)
  }
}
write.table(ATAC_peak_HDAC1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC1_list.txt",sep="\t")

ATAC_peak_HDAC1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC1_list.txt",sep="\t")

HDAC1_features=list("HDAC1" = ATAC_peak_HDAC1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HDAC1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##EGR1
EGR1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF640AKG_EGR1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(EGR1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF640AKG_EGR1_age.txt",sep = "\t")

EGR1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF640AKG_EGR1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(EGR1)){
  x=paste(EGR1[i,1],EGR1[i,2],EGR1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

EGR1_features = list(list$V1)

ATAC_peak_EGR1_list = data.frame()
for(i in 1:nrow(EGR1)){
  ATAC_peak_EGR1<- ATAC_peak[ATAC_peak$x1 ==EGR1[i,1] & ((EGR1[i,2] > ATAC_peak$start &EGR1[i,2] < ATAC_peak$end) | (EGR1[i,3] > ATAC_peak$start &EGR1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_EGR1) != 0){
    ATAC_peak_EGR1$V5 <-EGR1[i,1]
    ATAC_peak_EGR1$V6 <- EGR1[i,2]
    ATAC_peak_EGR1$V7 <- EGR1[i,3]
    ATAC_peak_EGR1$V8 <- EGR1[i,6]
    ATAC_peak_EGR1$V9 <- EGR1[i,7]
    ATAC_peak_EGR1_list = rbind(ATAC_peak_EGR1_list,ATAC_peak_EGR1)
  }
}
write.table(ATAC_peak_EGR1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EGR1_list.txt",sep="\t")

ATAC_peak_EGR1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_EGR1_list.txt",sep="\t")

EGR1_features=list("EGR1" = ATAC_peak_EGR1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = EGR1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##C11orf30
C11orf30 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF640MQS_C11orf30.bed.gz",extraCols=extraCols_narrowPeak)
write.table(C11orf30,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF640MQS_C11orf30.txt",sep = "\t")

C11orf30 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF640MQS_C11orf30.txt",sep = "\t")

list = vector()
for (i in  1:nrow(C11orf30)){
  x=paste(C11orf30[i,1],C11orf30[i,2],C11orf30[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

C11orf30_features = list(list$V1)

ATAC_peak_C11orf30_list = data.frame()
for(i in 1:nrow(C11orf30)){
  ATAC_peak_C11orf30<- ATAC_peak[ATAC_peak$x1 ==C11orf30[i,1] & ((C11orf30[i,2] > ATAC_peak$start &C11orf30[i,2] < ATAC_peak$end) | (C11orf30[i,3] > ATAC_peak$start &C11orf30[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_C11orf30) != 0){
    ATAC_peak_C11orf30$V5 <-C11orf30[i,1]
    ATAC_peak_C11orf30$V6 <- C11orf30[i,2]
    ATAC_peak_C11orf30$V7 <- C11orf30[i,3]
    ATAC_peak_C11orf30$V8 <- C11orf30[i,6]
    ATAC_peak_C11orf30$V9 <- C11orf30[i,7]
    ATAC_peak_C11orf30_list = rbind(ATAC_peak_C11orf30_list,ATAC_peak_C11orf30)
  }
}
write.table(ATAC_peak_C11orf30_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_C11orf30_list.txt",sep="\t")

ATAC_peak_C11orf30_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_C11orf30_list.txt",sep="\t")

C11orf30_features=list("C11orf30" = ATAC_peak_C11orf30_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = C11orf30_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZNF121
ZNF121 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF650DWZ_ZNF121.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF121,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF650DWZ_ZNF121.txt",sep = "\t")

ZNF121 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF650DWZ_ZNF121.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF121)){
  x=paste(ZNF121[i,1],ZNF121[i,2],ZNF121[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF121_features = list(list$V1)

ATAC_peak_ZNF121_list = data.frame()
for(i in 1:nrow(ZNF121)){
  ATAC_peak_ZNF121<- ATAC_peak[ATAC_peak$x1 ==ZNF121[i,1] & ((ZNF121[i,2] > ATAC_peak$start &ZNF121[i,2] < ATAC_peak$end) | (ZNF121[i,3] > ATAC_peak$start &ZNF121[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF121) != 0){
    ATAC_peak_ZNF121$V5 <-ZNF121[i,1]
    ATAC_peak_ZNF121$V6 <- ZNF121[i,2]
    ATAC_peak_ZNF121$V7 <- ZNF121[i,3]
    ATAC_peak_ZNF121$V8 <- ZNF121[i,6]
    ATAC_peak_ZNF121$V9 <- ZNF121[i,7]
    ATAC_peak_ZNF121_list = rbind(ATAC_peak_ZNF121_list,ATAC_peak_ZNF121)
  }
}
write.table(ATAC_peak_ZNF121_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF121_list.txt",sep="\t")

ATAC_peak_ZNF121_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF121_list.txt",sep="\t")

ZNF121_features=list("ZNF121" = ATAC_peak_ZNF121_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF121_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TBL1XR1
TBL1XR1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF650YZH_TBL1XR1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TBL1XR1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF650YZH_TBL1XR1.txt",sep = "\t")

TBL1XR1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF650YZH_TBL1XR1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TBL1XR1)){
  x=paste(TBL1XR1[i,1],TBL1XR1[i,2],TBL1XR1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TBL1XR1_features = list(list$V1)

ATAC_peak_TBL1XR1_list = data.frame()
for(i in 1:nrow(TBL1XR1)){
  ATAC_peak_TBL1XR1<- ATAC_peak[ATAC_peak$x1 ==TBL1XR1[i,1] & ((TBL1XR1[i,2] > ATAC_peak$start &TBL1XR1[i,2] < ATAC_peak$end) | (TBL1XR1[i,3] > ATAC_peak$start &TBL1XR1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TBL1XR1) != 0){
    ATAC_peak_TBL1XR1$V5 <-TBL1XR1[i,1]
    ATAC_peak_TBL1XR1$V6 <- TBL1XR1[i,2]
    ATAC_peak_TBL1XR1$V7 <- TBL1XR1[i,3]
    ATAC_peak_TBL1XR1$V8 <- TBL1XR1[i,6]
    ATAC_peak_TBL1XR1$V9 <- TBL1XR1[i,7]
    ATAC_peak_TBL1XR1_list = rbind(ATAC_peak_TBL1XR1_list,ATAC_peak_TBL1XR1)
  }
}
write.table(ATAC_peak_TBL1XR1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TBL1XR1_list.txt",sep="\t")

ATAC_peak_TBL1XR1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TBL1XR1_list.txt",sep="\t")

TBL1XR1_features=list("TBL1XR1" = ATAC_peak_TBL1XR1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TBL1XR1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##PHF21A
PHF21A = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF663JIA_PHF21A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PHF21A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF663JIA_PHF21A.txt",sep = "\t")

PHF21A = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF663JIA_PHF21A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PHF21A)){
  x=paste(PHF21A[i,1],PHF21A[i,2],PHF21A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PHF21A_features = list(list$V1)

ATAC_peak_PHF21A_list = data.frame()
for(i in 1:nrow(PHF21A)){
  ATAC_peak_PHF21A<- ATAC_peak[ATAC_peak$x1 ==PHF21A[i,1] & ((PHF21A[i,2] > ATAC_peak$start &PHF21A[i,2] < ATAC_peak$end) | (PHF21A[i,3] > ATAC_peak$start &PHF21A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PHF21A) != 0){
    ATAC_peak_PHF21A$V5 <-PHF21A[i,1]
    ATAC_peak_PHF21A$V6 <- PHF21A[i,2]
    ATAC_peak_PHF21A$V7 <- PHF21A[i,3]
    ATAC_peak_PHF21A$V8 <- PHF21A[i,6]
    ATAC_peak_PHF21A$V9 <- PHF21A[i,7]
    ATAC_peak_PHF21A_list = rbind(ATAC_peak_PHF21A_list,ATAC_peak_PHF21A)
  }
}
write.table(ATAC_peak_PHF21A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PHF21A_list.txt",sep="\t")

ATAC_peak_PHF21A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PHF21A_list.txt",sep="\t")

PHF21A_features=list("PHF21A" = ATAC_peak_PHF21A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PHF21A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TAL1
TAL1 = import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF665HHH_TAL1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TAL1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF665HHH_TAL1.txt",sep = "\t")

TAL1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF665HHH_TAL1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TAL1)){
  x=paste(TAL1[i,1],TAL1[i,2],TAL1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TAL1_features = list(list$V1)

ATAC_peak_TAL1_list = data.frame()
for(i in 1:nrow(TAL1)){
  ATAC_peak_TAL1<- ATAC_peak[ATAC_peak$x1 ==TAL1[i,1] & ((TAL1[i,2] > ATAC_peak$start &TAL1[i,2] < ATAC_peak$end) | (TAL1[i,3] > ATAC_peak$start &TAL1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TAL1) != 0){
    ATAC_peak_TAL1$V5 <-TAL1[i,1]
    ATAC_peak_TAL1$V6 <- TAL1[i,2]
    ATAC_peak_TAL1$V7 <- TAL1[i,3]
    ATAC_peak_TAL1$V8 <- TAL1[i,6]
    ATAC_peak_TAL1$V9 <- TAL1[i,7]
    ATAC_peak_TAL1_list = rbind(ATAC_peak_TAL1_list,ATAC_peak_TAL1)
  }
}
write.table(ATAC_peak_TAL1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAL1_list.txt",sep="\t")

ATAC_peak_TAL1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TAL1_list.txt",sep="\t")

TAL1_features=list("TAL1" = ATAC_peak_TAL1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TAL1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##SMARCC2
SMARCC2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF669FKJ_SMARCC2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMARCC2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF669FKJ_SMARCC2.txt",sep = "\t")

SMARCC2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF669FKJ_SMARCC2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMARCC2)){
  x=paste(SMARCC2[i,1],SMARCC2[i,2],SMARCC2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMARCC2_features = list(list$V1)

ATAC_peak_SMARCC2_list = data.frame()
for(i in 1:nrow(SMARCC2)){
  ATAC_peak_SMARCC2<- ATAC_peak[ATAC_peak$x1 ==SMARCC2[i,1] & ((SMARCC2[i,2] > ATAC_peak$start &SMARCC2[i,2] < ATAC_peak$end) | (SMARCC2[i,3] > ATAC_peak$start &SMARCC2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMARCC2) != 0){
    ATAC_peak_SMARCC2$V5 <-SMARCC2[i,1]
    ATAC_peak_SMARCC2$V6 <- SMARCC2[i,2]
    ATAC_peak_SMARCC2$V7 <- SMARCC2[i,3]
    ATAC_peak_SMARCC2$V8 <- SMARCC2[i,6]
    ATAC_peak_SMARCC2$V9 <- SMARCC2[i,7]
    ATAC_peak_SMARCC2_list = rbind(ATAC_peak_SMARCC2_list,ATAC_peak_SMARCC2)
  }
}
write.table(ATAC_peak_SMARCC2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMARCC2_list.txt",sep="\t")

ATAC_peak_SMARCC2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMARCC2_list.txt",sep="\t")

SMARCC2_features=list("SMARCC2" = ATAC_peak_SMARCC2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMARCC2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##AFF1
AFF1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF671NPK_AFF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(AFF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF671NPK_AFF1.txt",sep = "\t")

AFF1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF671NPK_AFF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(AFF1)){
  x=paste(AFF1[i,1],AFF1[i,2],AFF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

AFF1_features = list(list$V1)

ATAC_peak_AFF1_list = data.frame()
for(i in 1:nrow(AFF1)){
  ATAC_peak_AFF1<- ATAC_peak[ATAC_peak$x1 ==AFF1[i,1] & ((AFF1[i,2] > ATAC_peak$start &AFF1[i,2] < ATAC_peak$end) | (AFF1[i,3] > ATAC_peak$start &AFF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_AFF1) != 0){
    ATAC_peak_AFF1$V5 <-AFF1[i,1]
    ATAC_peak_AFF1$V6 <- AFF1[i,2]
    ATAC_peak_AFF1$V7 <- AFF1[i,3]
    ATAC_peak_AFF1$V8 <- AFF1[i,6]
    ATAC_peak_AFF1$V9 <- AFF1[i,7]
    ATAC_peak_AFF1_list = rbind(ATAC_peak_AFF1_list,ATAC_peak_AFF1)
  }
}
write.table(ATAC_peak_AFF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_AFF1_list.txt",sep="\t")

ATAC_peak_AFF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_AFF1_list.txt",sep="\t")

AFF1_features=list("AFF1" = ATAC_peak_AFF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = AFF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ATF2
ATF2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF673TEA_ATF2_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ATF2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF673TEA_ATF2_age.txt",sep = "\t")

ATF2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF673TEA_ATF2_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ATF2)){
  x=paste(ATF2[i,1],ATF2[i,2],ATF2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ATF2_features = list(list$V1)

ATAC_peak_ATF2_list = data.frame()
for(i in 1:nrow(ATF2)){
  ATAC_peak_ATF2<- ATAC_peak[ATAC_peak$x1 ==ATF2[i,1] & ((ATF2[i,2] > ATAC_peak$start &ATF2[i,2] < ATAC_peak$end) | (ATF2[i,3] > ATAC_peak$start &ATF2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ATF2) != 0){
    ATAC_peak_ATF2$V5 <-ATF2[i,1]
    ATAC_peak_ATF2$V6 <- ATF2[i,2]
    ATAC_peak_ATF2$V7 <- ATF2[i,3]
    ATAC_peak_ATF2$V8 <- ATF2[i,6]
    ATAC_peak_ATF2$V9 <- ATF2[i,7]
    ATAC_peak_ATF2_list = rbind(ATAC_peak_ATF2_list,ATAC_peak_ATF2)
  }
}
write.table(ATAC_peak_ATF2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATF2_list.txt",sep="\t")

ATAC_peak_ATF2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ATF2_list.txt",sep="\t")

ATF2_features=list("ATF2" = ATAC_peak_ATF2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ATF2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##PHB2
PHB2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF674BWC_PHB2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PHB2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF674BWC_PHB2.txt",sep = "\t")

PHB2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF674BWC_PHB2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PHB2)){
  x=paste(PHB2[i,1],PHB2[i,2],PHB2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PHB2_features = list(list$V1)

ATAC_peak_PHB2_list = data.frame()
for(i in 1:nrow(PHB2)){
  ATAC_peak_PHB2<- ATAC_peak[ATAC_peak$x1 ==PHB2[i,1] & ((PHB2[i,2] > ATAC_peak$start &PHB2[i,2] < ATAC_peak$end) | (PHB2[i,3] > ATAC_peak$start &PHB2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PHB2) != 0){
    ATAC_peak_PHB2$V5 <-PHB2[i,1]
    ATAC_peak_PHB2$V6 <- PHB2[i,2]
    ATAC_peak_PHB2$V7 <- PHB2[i,3]
    ATAC_peak_PHB2$V8 <- PHB2[i,6]
    ATAC_peak_PHB2$V9 <- PHB2[i,7]
    ATAC_peak_PHB2_list = rbind(ATAC_peak_PHB2_list,ATAC_peak_PHB2)
  }
}
write.table(ATAC_peak_PHB2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PHB2_list.txt",sep="\t")

ATAC_peak_PHB2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PHB2_list.txt",sep="\t")

PHB2_features=list("PHB2" = ATAC_peak_PHB2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PHB2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RFX1
RFX1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF678NAO_RFX1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RFX1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF678NAO_RFX1.txt",sep = "\t")

RFX1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF678NAO_RFX1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RFX1)){
  x=paste(RFX1[i,1],RFX1[i,2],RFX1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RFX1_features = list(list$V1)

ATAC_peak_RFX1_list = data.frame()
for(i in 1:nrow(RFX1)){
  ATAC_peak_RFX1<- ATAC_peak[ATAC_peak$x1 ==RFX1[i,1] & ((RFX1[i,2] > ATAC_peak$start &RFX1[i,2] < ATAC_peak$end) | (RFX1[i,3] > ATAC_peak$start &RFX1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RFX1) != 0){
    ATAC_peak_RFX1$V5 <-RFX1[i,1]
    ATAC_peak_RFX1$V6 <- RFX1[i,2]
    ATAC_peak_RFX1$V7 <- RFX1[i,3]
    ATAC_peak_RFX1$V8 <- RFX1[i,6]
    ATAC_peak_RFX1$V9 <- RFX1[i,7]
    ATAC_peak_RFX1_list = rbind(ATAC_peak_RFX1_list,ATAC_peak_RFX1)
  }
}
write.table(ATAC_peak_RFX1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RFX1_list.txt",sep="\t")

ATAC_peak_RFX1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RFX1_list.txt",sep="\t")

RFX1_features=list("RFX1" = ATAC_peak_RFX1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RFX1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CDC5L
CDC5L= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF680NJD_CDC5L.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CDC5L,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF680NJD_CDC5L.txt",sep = "\t")

CDC5L= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF680NJD_CDC5L.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CDC5L)){
  x=paste(CDC5L[i,1],CDC5L[i,2],CDC5L[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CDC5L_features = list(list$V1)

ATAC_peak_CDC5L_list = data.frame()
for(i in 1:nrow(CDC5L)){
  ATAC_peak_CDC5L<- ATAC_peak[ATAC_peak$x1 ==CDC5L[i,1] & ((CDC5L[i,2] > ATAC_peak$start &CDC5L[i,2] < ATAC_peak$end) | (CDC5L[i,3] > ATAC_peak$start &CDC5L[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CDC5L) != 0){
    ATAC_peak_CDC5L$V5 <-CDC5L[i,1]
    ATAC_peak_CDC5L$V6 <- CDC5L[i,2]
    ATAC_peak_CDC5L$V7 <- CDC5L[i,3]
    ATAC_peak_CDC5L$V8 <- CDC5L[i,6]
    ATAC_peak_CDC5L$V9 <- CDC5L[i,7]
    ATAC_peak_CDC5L_list = rbind(ATAC_peak_CDC5L_list,ATAC_peak_CDC5L)
  }
}
write.table(ATAC_peak_CDC5L_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CDC5L_list.txt",sep="\t")

ATAC_peak_CDC5L_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CDC5L_list.txt",sep="\t")

CDC5L_features=list("CDC5L" = ATAC_peak_CDC5L_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CDC5L_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MXD1
MXD1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF684AHJ_MXD1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MXD1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF684AHJ_MXD1_age.txt",sep = "\t")

MXD1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF684AHJ_MXD1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MXD1)){
  x=paste(MXD1[i,1],MXD1[i,2],MXD1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MXD1_features = list(list$V1)

ATAC_peak_MXD1_list = data.frame()
for(i in 1:nrow(MXD1)){
  ATAC_peak_MXD1<- ATAC_peak[ATAC_peak$x1 ==MXD1[i,1] & ((MXD1[i,2] > ATAC_peak$start &MXD1[i,2] < ATAC_peak$end) | (MXD1[i,3] > ATAC_peak$start &MXD1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MXD1) != 0){
    ATAC_peak_MXD1$V5 <-MXD1[i,1]
    ATAC_peak_MXD1$V6 <- MXD1[i,2]
    ATAC_peak_MXD1$V7 <- MXD1[i,3]
    ATAC_peak_MXD1$V8 <- MXD1[i,6]
    ATAC_peak_MXD1$V9 <- MXD1[i,7]
    ATAC_peak_MXD1_list = rbind(ATAC_peak_MXD1_list,ATAC_peak_MXD1)
  }
}
write.table(ATAC_peak_MXD1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MXD1_list.txt",sep="\t")

ATAC_peak_MXD1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MXD1_list.txt",sep="\t")

MXD1_features=list("MXD1" = ATAC_peak_MXD1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MXD1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##RELA
RELA= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF687VTH_RELA_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RELA,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF687VTH_RELA_age.txt",sep = "\t")

RELA= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF687VTH_RELA_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RELA)){
  x=paste(RELA[i,1],RELA[i,2],RELA[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RELA_features = list(list$V1)

ATAC_peak_RELA_list = data.frame()
for(i in 1:nrow(RELA)){
  ATAC_peak_RELA<- ATAC_peak[ATAC_peak$x1 ==RELA[i,1] & ((RELA[i,2] > ATAC_peak$start &RELA[i,2] < ATAC_peak$end) | (RELA[i,3] > ATAC_peak$start &RELA[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RELA) != 0){
    ATAC_peak_RELA$V5 <-RELA[i,1]
    ATAC_peak_RELA$V6 <- RELA[i,2]
    ATAC_peak_RELA$V7 <- RELA[i,3]
    ATAC_peak_RELA$V8 <- RELA[i,6]
    ATAC_peak_RELA$V9 <- RELA[i,7]
    ATAC_peak_RELA_list = rbind(ATAC_peak_RELA_list,ATAC_peak_RELA)
  }
}
write.table(ATAC_peak_RELA_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RELA_list.txt",sep="\t")

ATAC_peak_RELA_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RELA_list.txt",sep="\t")

RELA_features=list("RELA" = ATAC_peak_RELA_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RELA_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZMYM3
ZMYM3= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF700JAQ_ZMYM3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZMYM3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF700JAQ_ZMYM3.txt",sep = "\t")

ZMYM3= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF700JAQ_ZMYM3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZMYM3)){
  x=paste(ZMYM3[i,1],ZMYM3[i,2],ZMYM3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZMYM3_features = list(list$V1)

ATAC_peak_ZMYM3_list = data.frame()
for(i in 1:nrow(ZMYM3)){
  ATAC_peak_ZMYM3<- ATAC_peak[ATAC_peak$x1 ==ZMYM3[i,1] & ((ZMYM3[i,2] > ATAC_peak$start &ZMYM3[i,2] < ATAC_peak$end) | (ZMYM3[i,3] > ATAC_peak$start &ZMYM3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZMYM3) != 0){
    ATAC_peak_ZMYM3$V5 <-ZMYM3[i,1]
    ATAC_peak_ZMYM3$V6 <- ZMYM3[i,2]
    ATAC_peak_ZMYM3$V7 <- ZMYM3[i,3]
    ATAC_peak_ZMYM3$V8 <- ZMYM3[i,6]
    ATAC_peak_ZMYM3$V9 <- ZMYM3[i,7]
    ATAC_peak_ZMYM3_list = rbind(ATAC_peak_ZMYM3_list,ATAC_peak_ZMYM3)
  }
}
write.table(ATAC_peak_ZMYM3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZMYM3_list.txt",sep="\t")

ATAC_peak_ZMYM3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZMYM3_list.txt",sep="\t")

ZMYM3_features=list("ZMYM3" = ATAC_peak_ZMYM3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZMYM3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##NFE2L1
NFE2L1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF700WYH_NFE2L1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NFE2L1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF700WYH_NFE2L1_age.txt",sep = "\t")

NFE2L1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF700WYH_NFE2L1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NFE2L1)){
  x=paste(NFE2L1[i,1],NFE2L1[i,2],NFE2L1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NFE2L1_features = list(list$V1)

ATAC_peak_NFE2L1_list = data.frame()
for(i in 1:nrow(NFE2L1)){
  ATAC_peak_NFE2L1<- ATAC_peak[ATAC_peak$x1 ==NFE2L1[i,1] & ((NFE2L1[i,2] > ATAC_peak$start &NFE2L1[i,2] < ATAC_peak$end) | (NFE2L1[i,3] > ATAC_peak$start &NFE2L1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NFE2L1) != 0){
    ATAC_peak_NFE2L1$V5 <-NFE2L1[i,1]
    ATAC_peak_NFE2L1$V6 <- NFE2L1[i,2]
    ATAC_peak_NFE2L1$V7 <- NFE2L1[i,3]
    ATAC_peak_NFE2L1$V8 <- NFE2L1[i,6]
    ATAC_peak_NFE2L1$V9 <- NFE2L1[i,7]
    ATAC_peak_NFE2L1_list = rbind(ATAC_peak_NFE2L1_list,ATAC_peak_NFE2L1)
  }
}
write.table(ATAC_peak_NFE2L1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFE2L1_list.txt",sep="\t")

ATAC_peak_NFE2L1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFE2L1_list.txt",sep="\t")

NFE2L1_features=list("NFE2L1" = ATAC_peak_NFE2L1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NFE2L1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##TFDP1
TFDP1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF701EYT_TFDP1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TFDP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF701EYT_TFDP1_age.txt",sep = "\t")

TFDP1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF701EYT_TFDP1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TFDP1)){
  x=paste(TFDP1[i,1],TFDP1[i,2],TFDP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TFDP1_features = list(list$V1)

ATAC_peak_TFDP1_list = data.frame()
for(i in 1:nrow(TFDP1)){
  ATAC_peak_TFDP1<- ATAC_peak[ATAC_peak$x1 ==TFDP1[i,1] & ((TFDP1[i,2] > ATAC_peak$start &TFDP1[i,2] < ATAC_peak$end) | (TFDP1[i,3] > ATAC_peak$start &TFDP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TFDP1) != 0){
    ATAC_peak_TFDP1$V5 <-TFDP1[i,1]
    ATAC_peak_TFDP1$V6 <- TFDP1[i,2]
    ATAC_peak_TFDP1$V7 <- TFDP1[i,3]
    ATAC_peak_TFDP1$V8 <- TFDP1[i,6]
    ATAC_peak_TFDP1$V9 <- TFDP1[i,7]
    ATAC_peak_TFDP1_list = rbind(ATAC_peak_TFDP1_list,ATAC_peak_TFDP1)
  }
}
write.table(ATAC_peak_TFDP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TFDP1_list.txt",sep="\t")

ATAC_peak_TFDP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TFDP1_list.txt",sep="\t")

TFDP1_features=list("TFDP1" = ATAC_peak_TFDP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TFDP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##RCOR1
RCOR1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF703TCX_RCOR1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RCOR1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF703TCX_RCOR1.txt",sep = "\t")

RCOR1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF703TCX_RCOR1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RCOR1)){
  x=paste(RCOR1[i,1],RCOR1[i,2],RCOR1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RCOR1_features = list(list$V1)

ATAC_peak_RCOR1_list = data.frame()
for(i in 1:nrow(RCOR1)){
  ATAC_peak_RCOR1<- ATAC_peak[ATAC_peak$x1 ==RCOR1[i,1] & ((RCOR1[i,2] > ATAC_peak$start &RCOR1[i,2] < ATAC_peak$end) | (RCOR1[i,3] > ATAC_peak$start &RCOR1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RCOR1) != 0){
    ATAC_peak_RCOR1$V5 <-RCOR1[i,1]
    ATAC_peak_RCOR1$V6 <- RCOR1[i,2]
    ATAC_peak_RCOR1$V7 <- RCOR1[i,3]
    ATAC_peak_RCOR1$V8 <- RCOR1[i,6]
    ATAC_peak_RCOR1$V9 <- RCOR1[i,7]
    ATAC_peak_RCOR1_list = rbind(ATAC_peak_RCOR1_list,ATAC_peak_RCOR1)
  }
}
write.table(ATAC_peak_RCOR1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RCOR1_list.txt",sep="\t")

ATAC_peak_RCOR1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RCOR1_list.txt",sep="\t")

RCOR1_features=list("RCOR1" = ATAC_peak_RCOR1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RCOR1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)




##XRCC5
XRCC5= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF704NMU_XRCC5_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(XRCC5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF704NMU_XRCC5_age.txt",sep = "\t")

XRCC5= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF704NMU_XRCC5_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(XRCC5)){
  x=paste(XRCC5[i,1],XRCC5[i,2],XRCC5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

XRCC5_features = list(list$V1)

ATAC_peak_XRCC5_list = data.frame()
for(i in 1:nrow(XRCC5)){
  ATAC_peak_XRCC5<- ATAC_peak[ATAC_peak$x1 ==XRCC5[i,1] & ((XRCC5[i,2] > ATAC_peak$start &XRCC5[i,2] < ATAC_peak$end) | (XRCC5[i,3] > ATAC_peak$start &XRCC5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_XRCC5) != 0){
    ATAC_peak_XRCC5$V5 <-XRCC5[i,1]
    ATAC_peak_XRCC5$V6 <- XRCC5[i,2]
    ATAC_peak_XRCC5$V7 <- XRCC5[i,3]
    ATAC_peak_XRCC5$V8 <- XRCC5[i,6]
    ATAC_peak_XRCC5$V9 <- XRCC5[i,7]
    ATAC_peak_XRCC5_list = rbind(ATAC_peak_XRCC5_list,ATAC_peak_XRCC5)
  }
}
write.table(ATAC_peak_XRCC5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_XRCC5_list.txt",sep="\t")

ATAC_peak_XRCC5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_XRCC5_list.txt",sep="\t")

XRCC5_features=list("XRCC5" = ATAC_peak_XRCC5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = XRCC5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)




##NONO
NONO= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF708JRU_NONO.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NONO,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF708JRU_NONO.txt",sep = "\t")

NONO= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF708JRU_NONO.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NONO)){
  x=paste(NONO[i,1],NONO[i,2],NONO[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NONO_features = list(list$V1)

ATAC_peak_NONO_list = data.frame()
for(i in 1:nrow(NONO)){
  ATAC_peak_NONO<- ATAC_peak[ATAC_peak$x1 ==NONO[i,1] & ((NONO[i,2] > ATAC_peak$start &NONO[i,2] < ATAC_peak$end) | (NONO[i,3] > ATAC_peak$start &NONO[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NONO) != 0){
    ATAC_peak_NONO$V5 <-NONO[i,1]
    ATAC_peak_NONO$V6 <- NONO[i,2]
    ATAC_peak_NONO$V7 <- NONO[i,3]
    ATAC_peak_NONO$V8 <- NONO[i,6]
    ATAC_peak_NONO$V9 <- NONO[i,7]
    ATAC_peak_NONO_list = rbind(ATAC_peak_NONO_list,ATAC_peak_NONO)
  }
}
write.table(ATAC_peak_NONO_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NONO_list.txt",sep="\t")

ATAC_peak_NONO_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NONO_list.txt",sep="\t")

NONO_features=list("NONO" = ATAC_peak_NONO_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NONO_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##HCFC1
HCFC1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF709HJC_HCFC1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HCFC1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF709HJC_HCFC1.txt",sep = "\t")

HCFC1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF709HJC_HCFC1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HCFC1)){
  x=paste(HCFC1[i,1],HCFC1[i,2],HCFC1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HCFC1_features = list(list$V1)

ATAC_peak_HCFC1_list = data.frame()
for(i in 1:nrow(HCFC1)){
  ATAC_peak_HCFC1<- ATAC_peak[ATAC_peak$x1 ==HCFC1[i,1] & ((HCFC1[i,2] > ATAC_peak$start &HCFC1[i,2] < ATAC_peak$end) | (HCFC1[i,3] > ATAC_peak$start &HCFC1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HCFC1) != 0){
    ATAC_peak_HCFC1$V5 <-HCFC1[i,1]
    ATAC_peak_HCFC1$V6 <- HCFC1[i,2]
    ATAC_peak_HCFC1$V7 <- HCFC1[i,3]
    ATAC_peak_HCFC1$V8 <- HCFC1[i,6]
    ATAC_peak_HCFC1$V9 <- HCFC1[i,7]
    ATAC_peak_HCFC1_list = rbind(ATAC_peak_HCFC1_list,ATAC_peak_HCFC1)
  }
}
write.table(ATAC_peak_HCFC1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HCFC1_list.txt",sep="\t")

ATAC_peak_HCFC1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HCFC1_list.txt",sep="\t")

HCFC1_features=list("HCFC1" = ATAC_peak_HCFC1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HCFC1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##HDAC8
HDAC8= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF711GFV_HDAC8.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HDAC8,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF711GFV_HDAC8.txt",sep = "\t")

HDAC8= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF711GFV_HDAC8.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HDAC8)){
  x=paste(HDAC8[i,1],HDAC8[i,2],HDAC8[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HDAC8_features = list(list$V1)

ATAC_peak_HDAC8_list = data.frame()
for(i in 1:nrow(HDAC8)){
  ATAC_peak_HDAC8<- ATAC_peak[ATAC_peak$x1 ==HDAC8[i,1] & ((HDAC8[i,2] > ATAC_peak$start &HDAC8[i,2] < ATAC_peak$end) | (HDAC8[i,3] > ATAC_peak$start &HDAC8[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HDAC8) != 0){
    ATAC_peak_HDAC8$V5 <-HDAC8[i,1]
    ATAC_peak_HDAC8$V6 <- HDAC8[i,2]
    ATAC_peak_HDAC8$V7 <- HDAC8[i,3]
    ATAC_peak_HDAC8$V8 <- HDAC8[i,6]
    ATAC_peak_HDAC8$V9 <- HDAC8[i,7]
    ATAC_peak_HDAC8_list = rbind(ATAC_peak_HDAC8_list,ATAC_peak_HDAC8)
  }
}
write.table(ATAC_peak_HDAC8_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC8_list.txt",sep="\t")

ATAC_peak_HDAC8_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HDAC8_list.txt",sep="\t")

HDAC8_features=list("HDAC8" = ATAC_peak_HDAC8_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HDAC8_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##MCM7
MCM7= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF716AEV_MCM7.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MCM7,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF716AEV_MCM7.txt",sep = "\t")

MCM7= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF716AEV_MCM7.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MCM7)){
  x=paste(MCM7[i,1],MCM7[i,2],MCM7[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MCM7_features = list(list$V1)

ATAC_peak_MCM7_list = data.frame()
for(i in 1:nrow(MCM7)){
  ATAC_peak_MCM7<- ATAC_peak[ATAC_peak$x1 ==MCM7[i,1] & ((MCM7[i,2] > ATAC_peak$start &MCM7[i,2] < ATAC_peak$end) | (MCM7[i,3] > ATAC_peak$start &MCM7[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MCM7) != 0){
    ATAC_peak_MCM7$V5 <-MCM7[i,1]
    ATAC_peak_MCM7$V6 <- MCM7[i,2]
    ATAC_peak_MCM7$V7 <- MCM7[i,3]
    ATAC_peak_MCM7$V8 <- MCM7[i,6]
    ATAC_peak_MCM7$V9 <- MCM7[i,7]
    ATAC_peak_MCM7_list = rbind(ATAC_peak_MCM7_list,ATAC_peak_MCM7)
  }
}
write.table(ATAC_peak_MCM7_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MCM7_list.txt",sep="\t")

ATAC_peak_MCM7_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MCM7_list.txt",sep="\t")

MCM7_features=list("MCM7" = ATAC_peak_MCM7_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MCM7_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CHD7
CHD7= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF722UJW_CHD7.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CHD7,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF722UJW_CHD7.txt",sep = "\t")

CHD7= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF722UJW_CHD7.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CHD7)){
  x=paste(CHD7[i,1],CHD7[i,2],CHD7[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CHD7_features = list(list$V1)

ATAC_peak_CHD7_list = data.frame()
for(i in 1:nrow(CHD7)){
  ATAC_peak_CHD7<- ATAC_peak[ATAC_peak$x1 ==CHD7[i,1] & ((CHD7[i,2] > ATAC_peak$start &CHD7[i,2] < ATAC_peak$end) | (CHD7[i,3] > ATAC_peak$start &CHD7[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CHD7) != 0){
    ATAC_peak_CHD7$V5 <-CHD7[i,1]
    ATAC_peak_CHD7$V6 <- CHD7[i,2]
    ATAC_peak_CHD7$V7 <- CHD7[i,3]
    ATAC_peak_CHD7$V8 <- CHD7[i,6]
    ATAC_peak_CHD7$V9 <- CHD7[i,7]
    ATAC_peak_CHD7_list = rbind(ATAC_peak_CHD7_list,ATAC_peak_CHD7)
  }
}
write.table(ATAC_peak_CHD7_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHD7_list.txt",sep="\t")

ATAC_peak_CHD7_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CHD7_list.txt",sep="\t")

CHD7_features=list("CHD7" = ATAC_peak_CHD7_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CHD7_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##E4F1
E4F1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF724SFS_E4F1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(E4F1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF724SFS_E4F1.txt",sep = "\t")

E4F1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF724SFS_E4F1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(E4F1)){
  x=paste(E4F1[i,1],E4F1[i,2],E4F1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

E4F1_features = list(list$V1)

ATAC_peak_E4F1_list = data.frame()
for(i in 1:nrow(E4F1)){
  ATAC_peak_E4F1<- ATAC_peak[ATAC_peak$x1 ==E4F1[i,1] & ((E4F1[i,2] > ATAC_peak$start &E4F1[i,2] < ATAC_peak$end) | (E4F1[i,3] > ATAC_peak$start &E4F1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_E4F1) != 0){
    ATAC_peak_E4F1$V5 <-E4F1[i,1]
    ATAC_peak_E4F1$V6 <- E4F1[i,2]
    ATAC_peak_E4F1$V7 <- E4F1[i,3]
    ATAC_peak_E4F1$V8 <- E4F1[i,6]
    ATAC_peak_E4F1$V9 <- E4F1[i,7]
    ATAC_peak_E4F1_list = rbind(ATAC_peak_E4F1_list,ATAC_peak_E4F1)
  }
}
write.table(ATAC_peak_E4F1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_E4F1_list.txt",sep="\t")

ATAC_peak_E4F1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_E4F1_list.txt",sep="\t")

E4F1_features=list("E4F1" = ATAC_peak_E4F1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = E4F1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##PHF8
PHF8= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF727AEY_PHF8.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PHF8,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF727AEY_PHF8.txt",sep = "\t")

PHF8= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF727AEY_PHF8.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PHF8)){
  x=paste(PHF8[i,1],PHF8[i,2],PHF8[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PHF8_features = list(list$V1)

ATAC_peak_PHF8_list = data.frame()
for(i in 1:nrow(PHF8)){
  ATAC_peak_PHF8<- ATAC_peak[ATAC_peak$x1 ==PHF8[i,1] & ((PHF8[i,2] > ATAC_peak$start &PHF8[i,2] < ATAC_peak$end) | (PHF8[i,3] > ATAC_peak$start &PHF8[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PHF8) != 0){
    ATAC_peak_PHF8$V5 <-PHF8[i,1]
    ATAC_peak_PHF8$V6 <- PHF8[i,2]
    ATAC_peak_PHF8$V7 <- PHF8[i,3]
    ATAC_peak_PHF8$V8 <- PHF8[i,6]
    ATAC_peak_PHF8$V9 <- PHF8[i,7]
    ATAC_peak_PHF8_list = rbind(ATAC_peak_PHF8_list,ATAC_peak_PHF8)
  }
}
write.table(ATAC_peak_PHF8_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PHF8_list.txt",sep="\t")

ATAC_peak_PHF8_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PHF8_list.txt",sep="\t")

PHF8_features=list("PHF8" = ATAC_peak_PHF8_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PHF8_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CTCF
CTCF= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF736NYC_CTCF.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CTCF,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF736NYC_CTCF.txt",sep = "\t")

CTCF= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF736NYC_CTCF.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CTCF)){
  x=paste(CTCF[i,1],CTCF[i,2],CTCF[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CTCF_features = list(list$V1)

ATAC_peak_CTCF_list = data.frame()
for(i in 1:nrow(CTCF)){
  ATAC_peak_CTCF<- ATAC_peak[ATAC_peak$x1 ==CTCF[i,1] & ((CTCF[i,2] > ATAC_peak$start &CTCF[i,2] < ATAC_peak$end) | (CTCF[i,3] > ATAC_peak$start &CTCF[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CTCF) != 0){
    ATAC_peak_CTCF$V5 <-CTCF[i,1]
    ATAC_peak_CTCF$V6 <- CTCF[i,2]
    ATAC_peak_CTCF$V7 <- CTCF[i,3]
    ATAC_peak_CTCF$V8 <- CTCF[i,6]
    ATAC_peak_CTCF$V9 <- CTCF[i,7]
    ATAC_peak_CTCF_list = rbind(ATAC_peak_CTCF_list,ATAC_peak_CTCF)
  }
}
write.table(ATAC_peak_CTCF_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CTCF_list.txt",sep="\t")

ATAC_peak_CTCF_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CTCF_list.txt",sep="\t")

CTCF_features=list("CTCF" = ATAC_peak_CTCF_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CTCF_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##DACH1
DACH1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF742IRF_DACH1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(DACH1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF742IRF_DACH1.txt",sep = "\t")

DACH1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF742IRF_DACH1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(DACH1)){
  x=paste(DACH1[i,1],DACH1[i,2],DACH1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

DACH1_features = list(list$V1)

ATAC_peak_DACH1_list = data.frame()
for(i in 1:nrow(DACH1)){
  ATAC_peak_DACH1<- ATAC_peak[ATAC_peak$x1 ==DACH1[i,1] & ((DACH1[i,2] > ATAC_peak$start &DACH1[i,2] < ATAC_peak$end) | (DACH1[i,3] > ATAC_peak$start &DACH1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_DACH1) != 0){
    ATAC_peak_DACH1$V5 <-DACH1[i,1]
    ATAC_peak_DACH1$V6 <- DACH1[i,2]
    ATAC_peak_DACH1$V7 <- DACH1[i,3]
    ATAC_peak_DACH1$V8 <- DACH1[i,6]
    ATAC_peak_DACH1$V9 <- DACH1[i,7]
    ATAC_peak_DACH1_list = rbind(ATAC_peak_DACH1_list,ATAC_peak_DACH1)
  }
}
write.table(ATAC_peak_DACH1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DACH1_list.txt",sep="\t")

ATAC_peak_DACH1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DACH1_list.txt",sep="\t")

DACH1_features=list("DACH1" = ATAC_peak_DACH1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = DACH1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZBTB40
ZBTB40= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF745YDH_ZBTB40.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZBTB40,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF745YDH_ZBTB40.txt",sep = "\t")

ZBTB40= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF745YDH_ZBTB40.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZBTB40)){
  x=paste(ZBTB40[i,1],ZBTB40[i,2],ZBTB40[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZBTB40_features = list(list$V1)

ATAC_peak_ZBTB40_list = data.frame()
for(i in 1:nrow(ZBTB40)){
  ATAC_peak_ZBTB40<- ATAC_peak[ATAC_peak$x1 ==ZBTB40[i,1] & ((ZBTB40[i,2] > ATAC_peak$start &ZBTB40[i,2] < ATAC_peak$end) | (ZBTB40[i,3] > ATAC_peak$start &ZBTB40[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZBTB40) != 0){
    ATAC_peak_ZBTB40$V5 <-ZBTB40[i,1]
    ATAC_peak_ZBTB40$V6 <- ZBTB40[i,2]
    ATAC_peak_ZBTB40$V7 <- ZBTB40[i,3]
    ATAC_peak_ZBTB40$V8 <- ZBTB40[i,6]
    ATAC_peak_ZBTB40$V9 <- ZBTB40[i,7]
    ATAC_peak_ZBTB40_list = rbind(ATAC_peak_ZBTB40_list,ATAC_peak_ZBTB40)
  }
}
write.table(ATAC_peak_ZBTB40_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB40_list.txt",sep="\t")

ATAC_peak_ZBTB40_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB40_list.txt",sep="\t")

ZBTB40_features=list("ZBTB40" = ATAC_peak_ZBTB40_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZBTB40_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##BCOR
BCOR= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF746QQX_BCOR.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BCOR,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF746QQX_BCOR.txt",sep = "\t")

BCOR= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF746QQX_BCOR.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BCOR)){
  x=paste(BCOR[i,1],BCOR[i,2],BCOR[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BCOR_features = list(list$V1)

ATAC_peak_BCOR_list = data.frame()
for(i in 1:nrow(BCOR)){
  ATAC_peak_BCOR<- ATAC_peak[ATAC_peak$x1 ==BCOR[i,1] & ((BCOR[i,2] > ATAC_peak$start &BCOR[i,2] < ATAC_peak$end) | (BCOR[i,3] > ATAC_peak$start &BCOR[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BCOR) != 0){
    ATAC_peak_BCOR$V5 <-BCOR[i,1]
    ATAC_peak_BCOR$V6 <- BCOR[i,2]
    ATAC_peak_BCOR$V7 <- BCOR[i,3]
    ATAC_peak_BCOR$V8 <- BCOR[i,6]
    ATAC_peak_BCOR$V9 <- BCOR[i,7]
    ATAC_peak_BCOR_list = rbind(ATAC_peak_BCOR_list,ATAC_peak_BCOR)
  }
}
write.table(ATAC_peak_BCOR_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BCOR_list.txt",sep="\t")

ATAC_peak_BCOR_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BCOR_list.txt",sep="\t")

BCOR_features=list("BCOR" = ATAC_peak_BCOR_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BCOR_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##BRF2
BRF2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF748HCN_BRF2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BRF2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF748HCN_BRF2.txt",sep = "\t")

BRF2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF748HCN_BRF2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BRF2)){
  x=paste(BRF2[i,1],BRF2[i,2],BRF2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BRF2_features = list(list$V1)

ATAC_peak_BRF2_list = data.frame()
for(i in 1:nrow(BRF2)){
  ATAC_peak_BRF2<- ATAC_peak[ATAC_peak$x1 ==BRF2[i,1] & ((BRF2[i,2] > ATAC_peak$start &BRF2[i,2] < ATAC_peak$end) | (BRF2[i,3] > ATAC_peak$start &BRF2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BRF2) != 0){
    ATAC_peak_BRF2$V5 <-BRF2[i,1]
    ATAC_peak_BRF2$V6 <- BRF2[i,2]
    ATAC_peak_BRF2$V7 <- BRF2[i,3]
    ATAC_peak_BRF2$V8 <- BRF2[i,6]
    ATAC_peak_BRF2$V9 <- BRF2[i,7]
    ATAC_peak_BRF2_list = rbind(ATAC_peak_BRF2_list,ATAC_peak_BRF2)
  }
}
write.table(ATAC_peak_BRF2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRF2_list.txt",sep="\t")

ATAC_peak_BRF2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRF2_list.txt",sep="\t")

BRF2_features=list("BRF2" = ATAC_peak_BRF2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BRF2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##U2AF1
U2AF1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF752DQV_U2AF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(U2AF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF752DQV_U2AF1.txt",sep = "\t")

U2AF1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF752DQV_U2AF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(U2AF1)){
  x=paste(U2AF1[i,1],U2AF1[i,2],U2AF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

U2AF1_features = list(list$V1)

ATAC_peak_U2AF1_list = data.frame()
for(i in 1:nrow(U2AF1)){
  ATAC_peak_U2AF1<- ATAC_peak[ATAC_peak$x1 ==U2AF1[i,1] & ((U2AF1[i,2] > ATAC_peak$start &U2AF1[i,2] < ATAC_peak$end) | (U2AF1[i,3] > ATAC_peak$start &U2AF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_U2AF1) != 0){
    ATAC_peak_U2AF1$V5 <-U2AF1[i,1]
    ATAC_peak_U2AF1$V6 <- U2AF1[i,2]
    ATAC_peak_U2AF1$V7 <- U2AF1[i,3]
    ATAC_peak_U2AF1$V8 <- U2AF1[i,6]
    ATAC_peak_U2AF1$V9 <- U2AF1[i,7]
    ATAC_peak_U2AF1_list = rbind(ATAC_peak_U2AF1_list,ATAC_peak_U2AF1)
  }
}
write.table(ATAC_peak_U2AF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_U2AF1_list.txt",sep="\t")

ATAC_peak_U2AF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_U2AF1_list.txt",sep="\t")

U2AF1_features=list("U2AF1" = ATAC_peak_U2AF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = U2AF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CCAR2
CCAR2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF757ZBP_CCAR2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CCAR2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF757ZBP_CCAR2.txt",sep = "\t")

CCAR2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF757ZBP_CCAR2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CCAR2)){
  x=paste(CCAR2[i,1],CCAR2[i,2],CCAR2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CCAR2_features = list(list$V1)

ATAC_peak_CCAR2_list = data.frame()
for(i in 1:nrow(CCAR2)){
  ATAC_peak_CCAR2<- ATAC_peak[ATAC_peak$x1 ==CCAR2[i,1] & ((CCAR2[i,2] > ATAC_peak$start &CCAR2[i,2] < ATAC_peak$end) | (CCAR2[i,3] > ATAC_peak$start &CCAR2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CCAR2) != 0){
    ATAC_peak_CCAR2$V5 <-CCAR2[i,1]
    ATAC_peak_CCAR2$V6 <- CCAR2[i,2]
    ATAC_peak_CCAR2$V7 <- CCAR2[i,3]
    ATAC_peak_CCAR2$V8 <- CCAR2[i,6]
    ATAC_peak_CCAR2$V9 <- CCAR2[i,7]
    ATAC_peak_CCAR2_list = rbind(ATAC_peak_CCAR2_list,ATAC_peak_CCAR2)
  }
}
write.table(ATAC_peak_CCAR2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CCAR2_list.txt",sep="\t")

ATAC_peak_CCAR2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CCAR2_list.txt",sep="\t")

CCAR2_features=list("CCAR2" = ATAC_peak_CCAR2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CCAR2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CSDE1
CSDE1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF760YHF_CSDE1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CSDE1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF760YHF_CSDE1.txt",sep = "\t")

CSDE1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF760YHF_CSDE1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CSDE1)){
  x=paste(CSDE1[i,1],CSDE1[i,2],CSDE1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CSDE1_features = list(list$V1)

ATAC_peak_CSDE1_list = data.frame()
for(i in 1:nrow(CSDE1)){
  ATAC_peak_CSDE1<- ATAC_peak[ATAC_peak$x1 ==CSDE1[i,1] & ((CSDE1[i,2] > ATAC_peak$start &CSDE1[i,2] < ATAC_peak$end) | (CSDE1[i,3] > ATAC_peak$start &CSDE1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CSDE1) != 0){
    ATAC_peak_CSDE1$V5 <-CSDE1[i,1]
    ATAC_peak_CSDE1$V6 <- CSDE1[i,2]
    ATAC_peak_CSDE1$V7 <- CSDE1[i,3]
    ATAC_peak_CSDE1$V8 <- CSDE1[i,6]
    ATAC_peak_CSDE1$V9 <- CSDE1[i,7]
    ATAC_peak_CSDE1_list = rbind(ATAC_peak_CSDE1_list,ATAC_peak_CSDE1)
  }
}
write.table(ATAC_peak_CSDE1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CSDE1_list.txt",sep="\t")

ATAC_peak_CSDE1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CSDE1_list.txt",sep="\t")

CSDE1_features=list("CSDE1" = ATAC_peak_CSDE1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CSDE1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CEBPB
CEBPB= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF761UZG_CEBPB_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CEBPB,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF761UZG_CEBPB_age.txt",sep = "\t")

CEBPB= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF761UZG_CEBPB_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CEBPB)){
  x=paste(CEBPB[i,1],CEBPB[i,2],CEBPB[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CEBPB_features = list(list$V1)

ATAC_peak_CEBPB_list = data.frame()
for(i in 1:nrow(CEBPB)){
  ATAC_peak_CEBPB<- ATAC_peak[ATAC_peak$x1 ==CEBPB[i,1] & ((CEBPB[i,2] > ATAC_peak$start &CEBPB[i,2] < ATAC_peak$end) | (CEBPB[i,3] > ATAC_peak$start &CEBPB[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CEBPB) != 0){
    ATAC_peak_CEBPB$V5 <-CEBPB[i,1]
    ATAC_peak_CEBPB$V6 <- CEBPB[i,2]
    ATAC_peak_CEBPB$V7 <- CEBPB[i,3]
    ATAC_peak_CEBPB$V8 <- CEBPB[i,6]
    ATAC_peak_CEBPB$V9 <- CEBPB[i,7]
    ATAC_peak_CEBPB_list = rbind(ATAC_peak_CEBPB_list,ATAC_peak_CEBPB)
  }
}
write.table(ATAC_peak_CEBPB_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CEBPB_list.txt",sep="\t")

ATAC_peak_CEBPB_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CEBPB_list.txt",sep="\t")

CEBPB_features=list("CEBPB" = ATAC_peak_CEBPB_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CEBPB_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##GABPB1
GABPB1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF764HIS_GABPB1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GABPB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF764HIS_GABPB1.txt",sep = "\t")

GABPB1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF764HIS_GABPB1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GABPB1)){
  x=paste(GABPB1[i,1],GABPB1[i,2],GABPB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GABPB1_features = list(list$V1)

ATAC_peak_GABPB1_list = data.frame()
for(i in 1:nrow(GABPB1)){
  ATAC_peak_GABPB1<- ATAC_peak[ATAC_peak$x1 ==GABPB1[i,1] & ((GABPB1[i,2] > ATAC_peak$start &GABPB1[i,2] < ATAC_peak$end) | (GABPB1[i,3] > ATAC_peak$start &GABPB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GABPB1) != 0){
    ATAC_peak_GABPB1$V5 <-GABPB1[i,1]
    ATAC_peak_GABPB1$V6 <- GABPB1[i,2]
    ATAC_peak_GABPB1$V7 <- GABPB1[i,3]
    ATAC_peak_GABPB1$V8 <- GABPB1[i,6]
    ATAC_peak_GABPB1$V9 <- GABPB1[i,7]
    ATAC_peak_GABPB1_list = rbind(ATAC_peak_GABPB1_list,ATAC_peak_GABPB1)
  }
}
write.table(ATAC_peak_GABPB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GABPB1_list.txt",sep="\t")

ATAC_peak_GABPB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GABPB1_list.txt",sep="\t")

GABPB1_features=list("GABPB1" = ATAC_peak_GABPB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GABPB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##NBN
NBN= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF779XHR_NBN_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NBN,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF779XHR_NBN_age.txt",sep = "\t")

NBN= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF779XHR_NBN_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NBN)){
  x=paste(NBN[i,1],NBN[i,2],NBN[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NBN_features = list(list$V1)

ATAC_peak_NBN_list = data.frame()
for(i in 1:nrow(NBN)){
  ATAC_peak_NBN<- ATAC_peak[ATAC_peak$x1 ==NBN[i,1] & ((NBN[i,2] > ATAC_peak$start &NBN[i,2] < ATAC_peak$end) | (NBN[i,3] > ATAC_peak$start &NBN[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NBN) != 0){
    ATAC_peak_NBN$V5 <-NBN[i,1]
    ATAC_peak_NBN$V6 <- NBN[i,2]
    ATAC_peak_NBN$V7 <- NBN[i,3]
    ATAC_peak_NBN$V8 <- NBN[i,6]
    ATAC_peak_NBN$V9 <- NBN[i,7]
    ATAC_peak_NBN_list = rbind(ATAC_peak_NBN_list,ATAC_peak_NBN)
  }
}
write.table(ATAC_peak_NBN_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NBN_list.txt",sep="\t")

ATAC_peak_NBN_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NBN_list.txt",sep="\t")

NBN_features=list("NBN" = ATAC_peak_NBN_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NBN_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CC2D1A
CC2D1A= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF782WMI_CC2D1A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CC2D1A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF782WMI_CC2D1A.txt",sep = "\t")

CC2D1A= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF782WMI_CC2D1A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CC2D1A)){
  x=paste(CC2D1A[i,1],CC2D1A[i,2],CC2D1A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CC2D1A_features = list(list$V1)

ATAC_peak_CC2D1A_list = data.frame()
for(i in 1:nrow(CC2D1A)){
  ATAC_peak_CC2D1A<- ATAC_peak[ATAC_peak$x1 ==CC2D1A[i,1] & ((CC2D1A[i,2] > ATAC_peak$start &CC2D1A[i,2] < ATAC_peak$end) | (CC2D1A[i,3] > ATAC_peak$start &CC2D1A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CC2D1A) != 0){
    ATAC_peak_CC2D1A$V5 <-CC2D1A[i,1]
    ATAC_peak_CC2D1A$V6 <- CC2D1A[i,2]
    ATAC_peak_CC2D1A$V7 <- CC2D1A[i,3]
    ATAC_peak_CC2D1A$V8 <- CC2D1A[i,6]
    ATAC_peak_CC2D1A$V9 <- CC2D1A[i,7]
    ATAC_peak_CC2D1A_list = rbind(ATAC_peak_CC2D1A_list,ATAC_peak_CC2D1A)
  }
}
write.table(ATAC_peak_CC2D1A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CC2D1A_list.txt",sep="\t")

ATAC_peak_CC2D1A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CC2D1A_list.txt",sep="\t")

CC2D1A_features=list("CC2D1A" = ATAC_peak_CC2D1A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CC2D1A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##STAT5A
STAT5A= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF784XFN_STAT5A_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(STAT5A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF784XFN_STAT5A_age.txt",sep = "\t")

STAT5A= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF784XFN_STAT5A_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(STAT5A)){
  x=paste(STAT5A[i,1],STAT5A[i,2],STAT5A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

STAT5A_features = list(list$V1)

ATAC_peak_STAT5A_list = data.frame()
for(i in 1:nrow(STAT5A)){
  ATAC_peak_STAT5A<- ATAC_peak[ATAC_peak$x1 ==STAT5A[i,1] & ((STAT5A[i,2] > ATAC_peak$start &STAT5A[i,2] < ATAC_peak$end) | (STAT5A[i,3] > ATAC_peak$start &STAT5A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_STAT5A) != 0){
    ATAC_peak_STAT5A$V5 <-STAT5A[i,1]
    ATAC_peak_STAT5A$V6 <- STAT5A[i,2]
    ATAC_peak_STAT5A$V7 <- STAT5A[i,3]
    ATAC_peak_STAT5A$V8 <- STAT5A[i,6]
    ATAC_peak_STAT5A$V9 <- STAT5A[i,7]
    ATAC_peak_STAT5A_list = rbind(ATAC_peak_STAT5A_list,ATAC_peak_STAT5A)
  }
}
write.table(ATAC_peak_STAT5A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_STAT5A_list.txt",sep="\t")

ATAC_peak_STAT5A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_STAT5A_list.txt",sep="\t")

STAT5A_features=list("STAT5A" = ATAC_peak_STAT5A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = STAT5A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##TRIM25
TRIM25= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF792LVS_TRIM25.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TRIM25,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF792LVS_TRIM25.txt",sep = "\t")

TRIM25= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF792LVS_TRIM25.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TRIM25)){
  x=paste(TRIM25[i,1],TRIM25[i,2],TRIM25[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TRIM25_features = list(list$V1)

ATAC_peak_TRIM25_list = data.frame()
for(i in 1:nrow(TRIM25)){
  ATAC_peak_TRIM25<- ATAC_peak[ATAC_peak$x1 ==TRIM25[i,1] & ((TRIM25[i,2] > ATAC_peak$start &TRIM25[i,2] < ATAC_peak$end) | (TRIM25[i,3] > ATAC_peak$start &TRIM25[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TRIM25) != 0){
    ATAC_peak_TRIM25$V5 <-TRIM25[i,1]
    ATAC_peak_TRIM25$V6 <- TRIM25[i,2]
    ATAC_peak_TRIM25$V7 <- TRIM25[i,3]
    ATAC_peak_TRIM25$V8 <- TRIM25[i,6]
    ATAC_peak_TRIM25$V9 <- TRIM25[i,7]
    ATAC_peak_TRIM25_list = rbind(ATAC_peak_TRIM25_list,ATAC_peak_TRIM25)
  }
}
write.table(ATAC_peak_TRIM25_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TRIM25_list.txt",sep="\t")

ATAC_peak_TRIM25_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TRIM25_list.txt",sep="\t")

TRIM25_features=list("TRIM25" = ATAC_peak_TRIM25_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TRIM25_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZNF639
ZNF639= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF792AUO_ZNF639.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF639,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF792AUO_ZNF639.txt",sep = "\t")

ZNF639= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF792AUO_ZNF639.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF639)){
  x=paste(ZNF639[i,1],ZNF639[i,2],ZNF639[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF639_features = list(list$V1)

ATAC_peak_ZNF639_list = data.frame()
for(i in 1:nrow(ZNF639)){
  ATAC_peak_ZNF639<- ATAC_peak[ATAC_peak$x1 ==ZNF639[i,1] & ((ZNF639[i,2] > ATAC_peak$start &ZNF639[i,2] < ATAC_peak$end) | (ZNF639[i,3] > ATAC_peak$start &ZNF639[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF639) != 0){
    ATAC_peak_ZNF639$V5 <-ZNF639[i,1]
    ATAC_peak_ZNF639$V6 <- ZNF639[i,2]
    ATAC_peak_ZNF639$V7 <- ZNF639[i,3]
    ATAC_peak_ZNF639$V8 <- ZNF639[i,6]
    ATAC_peak_ZNF639$V9 <- ZNF639[i,7]
    ATAC_peak_ZNF639_list = rbind(ATAC_peak_ZNF639_list,ATAC_peak_ZNF639)
  }
}
write.table(ATAC_peak_ZNF639_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF639_list.txt",sep="\t")

ATAC_peak_ZNF639_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF639_list.txt",sep="\t")

ZNF639_features=list("ZNF639" = ATAC_peak_ZNF639_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF639_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##HNRNPH1
HNRNPH1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF797TVO_HNRNPH1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HNRNPH1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF797TVO_HNRNPH1.txt",sep = "\t")

HNRNPH1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF797TVO_HNRNPH1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HNRNPH1)){
  x=paste(HNRNPH1[i,1],HNRNPH1[i,2],HNRNPH1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HNRNPH1_features = list(list$V1)

ATAC_peak_HNRNPH1_list = data.frame()
for(i in 1:nrow(HNRNPH1)){
  ATAC_peak_HNRNPH1<- ATAC_peak[ATAC_peak$x1 ==HNRNPH1[i,1] & ((HNRNPH1[i,2] > ATAC_peak$start &HNRNPH1[i,2] < ATAC_peak$end) | (HNRNPH1[i,3] > ATAC_peak$start &HNRNPH1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HNRNPH1) != 0){
    ATAC_peak_HNRNPH1$V5 <-HNRNPH1[i,1]
    ATAC_peak_HNRNPH1$V6 <- HNRNPH1[i,2]
    ATAC_peak_HNRNPH1$V7 <- HNRNPH1[i,3]
    ATAC_peak_HNRNPH1$V8 <- HNRNPH1[i,6]
    ATAC_peak_HNRNPH1$V9 <- HNRNPH1[i,7]
    ATAC_peak_HNRNPH1_list = rbind(ATAC_peak_HNRNPH1_list,ATAC_peak_HNRNPH1)
  }
}
write.table(ATAC_peak_HNRNPH1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HNRNPH1_list.txt",sep="\t")

ATAC_peak_HNRNPH1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HNRNPH1_list.txt",sep="\t")

HNRNPH1_features=list("HNRNPH1" = ATAC_peak_HNRNPH1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HNRNPH1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##RB1
RB1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF805SGG_RB1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF805SGG_RB1_age.txt",sep = "\t")

RB1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF805SGG_RB1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RB1)){
  x=paste(RB1[i,1],RB1[i,2],RB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RB1_features = list(list$V1)

ATAC_peak_RB1_list = data.frame()
for(i in 1:nrow(RB1)){
  ATAC_peak_RB1<- ATAC_peak[ATAC_peak$x1 ==RB1[i,1] & ((RB1[i,2] > ATAC_peak$start &RB1[i,2] < ATAC_peak$end) | (RB1[i,3] > ATAC_peak$start &RB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RB1) != 0){
    ATAC_peak_RB1$V5 <-RB1[i,1]
    ATAC_peak_RB1$V6 <- RB1[i,2]
    ATAC_peak_RB1$V7 <- RB1[i,3]
    ATAC_peak_RB1$V8 <- RB1[i,6]
    ATAC_peak_RB1$V9 <- RB1[i,7]
    ATAC_peak_RB1_list = rbind(ATAC_peak_RB1_list,ATAC_peak_RB1)
  }
}
write.table(ATAC_peak_RB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RB1_list.txt",sep="\t")

ATAC_peak_RB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RB1_list.txt",sep="\t")

RB1_features=list("RB1" = ATAC_peak_RB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##FOXJ2
FOXJ2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF812WVH_FOXJ2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(FOXJ2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF812WVH_FOXJ2.txt",sep = "\t")

FOXJ2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF812WVH_FOXJ2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(FOXJ2)){
  x=paste(FOXJ2[i,1],FOXJ2[i,2],FOXJ2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

FOXJ2_features = list(list$V1)

ATAC_peak_FOXJ2_list = data.frame()
for(i in 1:nrow(FOXJ2)){
  ATAC_peak_FOXJ2<- ATAC_peak[ATAC_peak$x1 ==FOXJ2[i,1] & ((FOXJ2[i,2] > ATAC_peak$start &FOXJ2[i,2] < ATAC_peak$end) | (FOXJ2[i,3] > ATAC_peak$start &FOXJ2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_FOXJ2) != 0){
    ATAC_peak_FOXJ2$V5 <-FOXJ2[i,1]
    ATAC_peak_FOXJ2$V6 <- FOXJ2[i,2]
    ATAC_peak_FOXJ2$V7 <- FOXJ2[i,3]
    ATAC_peak_FOXJ2$V8 <- FOXJ2[i,6]
    ATAC_peak_FOXJ2$V9 <- FOXJ2[i,7]
    ATAC_peak_FOXJ2_list = rbind(ATAC_peak_FOXJ2_list,ATAC_peak_FOXJ2)
  }
}
write.table(ATAC_peak_FOXJ2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXJ2_list.txt",sep="\t")

ATAC_peak_FOXJ2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXJ2_list.txt",sep="\t")

FOXJ2_features=list("FOXJ2" = ATAC_peak_FOXJ2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = FOXJ2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##CBX5
CBX5= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF813SKN_CBX5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CBX5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF813SKN_CBX5.txt",sep = "\t")

CBX5= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF813SKN_CBX5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CBX5)){
  x=paste(CBX5[i,1],CBX5[i,2],CBX5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CBX5_features = list(list$V1)

ATAC_peak_CBX5_list = data.frame()
for(i in 1:nrow(CBX5)){
  ATAC_peak_CBX5<- ATAC_peak[ATAC_peak$x1 ==CBX5[i,1] & ((CBX5[i,2] > ATAC_peak$start &CBX5[i,2] < ATAC_peak$end) | (CBX5[i,3] > ATAC_peak$start &CBX5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CBX5) != 0){
    ATAC_peak_CBX5$V5 <-CBX5[i,1]
    ATAC_peak_CBX5$V6 <- CBX5[i,2]
    ATAC_peak_CBX5$V7 <- CBX5[i,3]
    ATAC_peak_CBX5$V8 <- CBX5[i,6]
    ATAC_peak_CBX5$V9 <- CBX5[i,7]
    ATAC_peak_CBX5_list = rbind(ATAC_peak_CBX5_list,ATAC_peak_CBX5)
  }
}
write.table(ATAC_peak_CBX5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX5_list.txt",sep="\t")

ATAC_peak_CBX5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CBX5_list.txt",sep="\t")

CBX5_features=list("CBX5" = ATAC_peak_CBX5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CBX5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##MYNN
MYNN= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF813XIL_MYNN.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MYNN,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF813XIL_MYNN.txt",sep = "\t")

MYNN= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF813XIL_MYNN.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MYNN)){
  x=paste(MYNN[i,1],MYNN[i,2],MYNN[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MYNN_features = list(list$V1)

ATAC_peak_MYNN_list = data.frame()
for(i in 1:nrow(MYNN)){
  ATAC_peak_MYNN<- ATAC_peak[ATAC_peak$x1 ==MYNN[i,1] & ((MYNN[i,2] > ATAC_peak$start &MYNN[i,2] < ATAC_peak$end) | (MYNN[i,3] > ATAC_peak$start &MYNN[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MYNN) != 0){
    ATAC_peak_MYNN$V5 <-MYNN[i,1]
    ATAC_peak_MYNN$V6 <- MYNN[i,2]
    ATAC_peak_MYNN$V7 <- MYNN[i,3]
    ATAC_peak_MYNN$V8 <- MYNN[i,6]
    ATAC_peak_MYNN$V9 <- MYNN[i,7]
    ATAC_peak_MYNN_list = rbind(ATAC_peak_MYNN_list,ATAC_peak_MYNN)
  }
}
write.table(ATAC_peak_MYNN_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MYNN_list.txt",sep="\t")

ATAC_peak_MYNN_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MYNN_list.txt",sep="\t")

MYNN_features=list("MYNN" = ATAC_peak_MYNN_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MYNN_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

#save(scmulti2,file = '/Users/ramzipit/Desktop/lnc_project/scmulti2_0530.Rdata')



##NR3C1
NR3C1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF818BIH_NR3C1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NR3C1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF818BIH_NR3C1_age.txt",sep = "\t")

NR3C1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF818BIH_NR3C1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NR3C1)){
  x=paste(NR3C1[i,1],NR3C1[i,2],NR3C1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NR3C1_features = list(list$V1)

ATAC_peak_NR3C1_list = data.frame()
for(i in 1:nrow(NR3C1)){
  ATAC_peak_NR3C1<- ATAC_peak[ATAC_peak$x1 ==NR3C1[i,1] & ((NR3C1[i,2] > ATAC_peak$start &NR3C1[i,2] < ATAC_peak$end) | (NR3C1[i,3] > ATAC_peak$start &NR3C1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NR3C1) != 0){
    ATAC_peak_NR3C1$V5 <-NR3C1[i,1]
    ATAC_peak_NR3C1$V6 <- NR3C1[i,2]
    ATAC_peak_NR3C1$V7 <- NR3C1[i,3]
    ATAC_peak_NR3C1$V8 <- NR3C1[i,6]
    ATAC_peak_NR3C1$V9 <- NR3C1[i,7]
    ATAC_peak_NR3C1_list = rbind(ATAC_peak_NR3C1_list,ATAC_peak_NR3C1)
  }
}
write.table(ATAC_peak_NR3C1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NR3C1_list.txt",sep="\t")

ATAC_peak_NR3C1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NR3C1_list.txt",sep="\t")

NR3C1_features=list("NR3C1" = ATAC_peak_NR3C1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NR3C1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SIRT6
SIRT6= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF821XJU_SIRT6_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SIRT6,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF821XJU_SIRT6_age.txt",sep = "\t")

SIRT6= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF821XJU_SIRT6_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SIRT6)){
  x=paste(SIRT6[i,1],SIRT6[i,2],SIRT6[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SIRT6_features = list(list$V1)

ATAC_peak_SIRT6_list = data.frame()
for(i in 1:nrow(SIRT6)){
  ATAC_peak_SIRT6<- ATAC_peak[ATAC_peak$x1 ==SIRT6[i,1] & ((SIRT6[i,2] > ATAC_peak$start &SIRT6[i,2] < ATAC_peak$end) | (SIRT6[i,3] > ATAC_peak$start &SIRT6[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SIRT6) != 0){
    ATAC_peak_SIRT6$V5 <-SIRT6[i,1]
    ATAC_peak_SIRT6$V6 <- SIRT6[i,2]
    ATAC_peak_SIRT6$V7 <- SIRT6[i,3]
    ATAC_peak_SIRT6$V8 <- SIRT6[i,6]
    ATAC_peak_SIRT6$V9 <- SIRT6[i,7]
    ATAC_peak_SIRT6_list = rbind(ATAC_peak_SIRT6_list,ATAC_peak_SIRT6)
  }
}
write.table(ATAC_peak_SIRT6_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIRT6_list.txt",sep="\t")

ATAC_peak_SIRT6_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIRT6_list.txt",sep="\t")

SIRT6_features=list("SIRT6" = ATAC_peak_SIRT6_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SIRT6_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NFXL1
NFXL1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF823RGT_NFXL1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NFXL1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF823RGT_NFXL1.txt",sep = "\t")

NFXL1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF823RGT_NFXL1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NFXL1)){
  x=paste(NFXL1[i,1],NFXL1[i,2],NFXL1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NFXL1_features = list(list$V1)

ATAC_peak_NFXL1_list = data.frame()
for(i in 1:nrow(NFXL1)){
  ATAC_peak_NFXL1<- ATAC_peak[ATAC_peak$x1 ==NFXL1[i,1] & ((NFXL1[i,2] > ATAC_peak$start &NFXL1[i,2] < ATAC_peak$end) | (NFXL1[i,3] > ATAC_peak$start &NFXL1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NFXL1) != 0){
    ATAC_peak_NFXL1$V5 <-NFXL1[i,1]
    ATAC_peak_NFXL1$V6 <- NFXL1[i,2]
    ATAC_peak_NFXL1$V7 <- NFXL1[i,3]
    ATAC_peak_NFXL1$V8 <- NFXL1[i,6]
    ATAC_peak_NFXL1$V9 <- NFXL1[i,7]
    ATAC_peak_NFXL1_list = rbind(ATAC_peak_NFXL1_list,ATAC_peak_NFXL1)
  }
}
write.table(ATAC_peak_NFXL1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFXL1_list.txt",sep="\t")

ATAC_peak_NFXL1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFXL1_list.txt",sep="\t")

NFXL1_features=list("NFXL1" = ATAC_peak_NFXL1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NFXL1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##CREM
CREM= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF832RYO_CREM.bed.gz",extraCols=extraCols_narrowPeak)
write.table(CREM,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF832RYO_CREM.txt",sep = "\t")

CREM= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF832RYO_CREM.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CREM)){
  x=paste(CREM[i,1],CREM[i,2],CREM[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CREM_features = list(list$V1)

ATAC_peak_CREM_list = data.frame()
for(i in 1:nrow(CREM)){
  ATAC_peak_CREM<- ATAC_peak[ATAC_peak$x1 ==CREM[i,1] & ((CREM[i,2] > ATAC_peak$start &CREM[i,2] < ATAC_peak$end) | (CREM[i,3] > ATAC_peak$start &CREM[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CREM) != 0){
    ATAC_peak_CREM$V5 <-CREM[i,1]
    ATAC_peak_CREM$V6 <- CREM[i,2]
    ATAC_peak_CREM$V7 <- CREM[i,3]
    ATAC_peak_CREM$V8 <- CREM[i,6]
    ATAC_peak_CREM$V9 <- CREM[i,7]
    ATAC_peak_CREM_list = rbind(ATAC_peak_CREM_list,ATAC_peak_CREM)
  }
}
write.table(ATAC_peak_CREM_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CREM_list.txt",sep="\t")

ATAC_peak_CREM_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CREM_list.txt",sep="\t")

CREM_features=list("CREM" = ATAC_peak_CREM_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CREM_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ARID4B
ARID4B= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF832YJP_ARID4B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ARID4B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF832YJP_ARID4B.txt",sep = "\t")

ARID4B= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF832YJP_ARID4B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ARID4B)){
  x=paste(ARID4B[i,1],ARID4B[i,2],ARID4B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ARID4B_features = list(list$V1)

ATAC_peak_ARID4B_list = data.frame()
for(i in 1:nrow(ARID4B)){
  ATAC_peak_ARID4B<- ATAC_peak[ATAC_peak$x1 ==ARID4B[i,1] & ((ARID4B[i,2] > ATAC_peak$start &ARID4B[i,2] < ATAC_peak$end) | (ARID4B[i,3] > ATAC_peak$start &ARID4B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ARID4B) != 0){
    ATAC_peak_ARID4B$V5 <-ARID4B[i,1]
    ATAC_peak_ARID4B$V6 <- ARID4B[i,2]
    ATAC_peak_ARID4B$V7 <- ARID4B[i,3]
    ATAC_peak_ARID4B$V8 <- ARID4B[i,6]
    ATAC_peak_ARID4B$V9 <- ARID4B[i,7]
    ATAC_peak_ARID4B_list = rbind(ATAC_peak_ARID4B_list,ATAC_peak_ARID4B)
  }
}
write.table(ATAC_peak_ARID4B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARID4B_list.txt",sep="\t")

ATAC_peak_ARID4B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARID4B_list.txt",sep="\t")

ARID4B_features=list("ARID4B" = ATAC_peak_ARID4B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ARID4B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##SNIP1
SNIP1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF846QZA_SNIP1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SNIP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF846QZA_SNIP1.txt",sep = "\t")

SNIP1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF846QZA_SNIP1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SNIP1)){
  x=paste(SNIP1[i,1],SNIP1[i,2],SNIP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SNIP1_features = list(list$V1)

ATAC_peak_SNIP1_list = data.frame()
for(i in 1:nrow(SNIP1)){
  ATAC_peak_SNIP1<- ATAC_peak[ATAC_peak$x1 ==SNIP1[i,1] & ((SNIP1[i,2] > ATAC_peak$start &SNIP1[i,2] < ATAC_peak$end) | (SNIP1[i,3] > ATAC_peak$start &SNIP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SNIP1) != 0){
    ATAC_peak_SNIP1$V5 <-SNIP1[i,1]
    ATAC_peak_SNIP1$V6 <- SNIP1[i,2]
    ATAC_peak_SNIP1$V7 <- SNIP1[i,3]
    ATAC_peak_SNIP1$V8 <- SNIP1[i,6]
    ATAC_peak_SNIP1$V9 <- SNIP1[i,7]
    ATAC_peak_SNIP1_list = rbind(ATAC_peak_SNIP1_list,ATAC_peak_SNIP1)
  }
}
write.table(ATAC_peak_SNIP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SNIP1_list.txt",sep="\t")

ATAC_peak_SNIP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SNIP1_list.txt",sep="\t")

SNIP1_features=list("SNIP1" = ATAC_peak_SNIP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SNIP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##IKZF1
IKZF1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF847SXQ_IKZF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(IKZF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF847SXQ_IKZF1.txt",sep = "\t")

IKZF1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF847SXQ_IKZF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(IKZF1)){
  x=paste(IKZF1[i,1],IKZF1[i,2],IKZF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

IKZF1_features = list(list$V1)

ATAC_peak_IKZF1_list = data.frame()
for(i in 1:nrow(IKZF1)){
  ATAC_peak_IKZF1<- ATAC_peak[ATAC_peak$x1 ==IKZF1[i,1] & ((IKZF1[i,2] > ATAC_peak$start &IKZF1[i,2] < ATAC_peak$end) | (IKZF1[i,3] > ATAC_peak$start &IKZF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_IKZF1) != 0){
    ATAC_peak_IKZF1$V5 <-IKZF1[i,1]
    ATAC_peak_IKZF1$V6 <- IKZF1[i,2]
    ATAC_peak_IKZF1$V7 <- IKZF1[i,3]
    ATAC_peak_IKZF1$V8 <- IKZF1[i,6]
    ATAC_peak_IKZF1$V9 <- IKZF1[i,7]
    ATAC_peak_IKZF1_list = rbind(ATAC_peak_IKZF1_list,ATAC_peak_IKZF1)
  }
}
write.table(ATAC_peak_IKZF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_IKZF1_list.txt",sep="\t")

ATAC_peak_IKZF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_IKZF1_list.txt",sep="\t")

IKZF1_features=list("IKZF1" = ATAC_peak_IKZF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = IKZF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZNF217
ZNF217= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF849FMB_ZNF217.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF217,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF849FMB_ZNF217.txt",sep = "\t")

ZNF217= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF849FMB_ZNF217.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF217)){
  x=paste(ZNF217[i,1],ZNF217[i,2],ZNF217[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF217_features = list(list$V1)

ATAC_peak_ZNF217_list = data.frame()
for(i in 1:nrow(ZNF217)){
  ATAC_peak_ZNF217<- ATAC_peak[ATAC_peak$x1 ==ZNF217[i,1] & ((ZNF217[i,2] > ATAC_peak$start &ZNF217[i,2] < ATAC_peak$end) | (ZNF217[i,3] > ATAC_peak$start &ZNF217[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF217) != 0){
    ATAC_peak_ZNF217$V5 <-ZNF217[i,1]
    ATAC_peak_ZNF217$V6 <- ZNF217[i,2]
    ATAC_peak_ZNF217$V7 <- ZNF217[i,3]
    ATAC_peak_ZNF217$V8 <- ZNF217[i,6]
    ATAC_peak_ZNF217$V9 <- ZNF217[i,7]
    ATAC_peak_ZNF217_list = rbind(ATAC_peak_ZNF217_list,ATAC_peak_ZNF217)
  }
}
write.table(ATAC_peak_ZNF217_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF217_list.txt",sep="\t")

ATAC_peak_ZNF217_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF217_list.txt",sep="\t")

ZNF217_features=list("ZNF217" = ATAC_peak_ZNF217_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF217_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##ZBTB11
ZBTB11= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF860DDA_ZBTB11.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZBTB11,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF860DDA_ZBTB11.txt",sep = "\t")

ZBTB11= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF860DDA_ZBTB11.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZBTB11)){
  x=paste(ZBTB11[i,1],ZBTB11[i,2],ZBTB11[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZBTB11_features = list(list$V1)

ATAC_peak_ZBTB11_list = data.frame()
for(i in 1:nrow(ZBTB11)){
  ATAC_peak_ZBTB11<- ATAC_peak[ATAC_peak$x1 ==ZBTB11[i,1] & ((ZBTB11[i,2] > ATAC_peak$start &ZBTB11[i,2] < ATAC_peak$end) | (ZBTB11[i,3] > ATAC_peak$start &ZBTB11[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZBTB11) != 0){
    ATAC_peak_ZBTB11$V5 <-ZBTB11[i,1]
    ATAC_peak_ZBTB11$V6 <- ZBTB11[i,2]
    ATAC_peak_ZBTB11$V7 <- ZBTB11[i,3]
    ATAC_peak_ZBTB11$V8 <- ZBTB11[i,6]
    ATAC_peak_ZBTB11$V9 <- ZBTB11[i,7]
    ATAC_peak_ZBTB11_list = rbind(ATAC_peak_ZBTB11_list,ATAC_peak_ZBTB11)
  }
}
write.table(ATAC_peak_ZBTB11_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB11_list.txt",sep="\t")

ATAC_peak_ZBTB11_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZBTB11_list.txt",sep="\t")

ZBTB11_features=list("ZBTB11" = ATAC_peak_ZBTB11_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZBTB11_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##MTA2
MTA2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF862EPJ_MTA2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MTA2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF862EPJ_MTA2.txt",sep = "\t")

MTA2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF862EPJ_MTA2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MTA2)){
  x=paste(MTA2[i,1],MTA2[i,2],MTA2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MTA2_features = list(list$V1)

ATAC_peak_MTA2_list = data.frame()
for(i in 1:nrow(MTA2)){
  ATAC_peak_MTA2<- ATAC_peak[ATAC_peak$x1 ==MTA2[i,1] & ((MTA2[i,2] > ATAC_peak$start &MTA2[i,2] < ATAC_peak$end) | (MTA2[i,3] > ATAC_peak$start &MTA2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MTA2) != 0){
    ATAC_peak_MTA2$V5 <-MTA2[i,1]
    ATAC_peak_MTA2$V6 <- MTA2[i,2]
    ATAC_peak_MTA2$V7 <- MTA2[i,3]
    ATAC_peak_MTA2$V8 <- MTA2[i,6]
    ATAC_peak_MTA2$V9 <- MTA2[i,7]
    ATAC_peak_MTA2_list = rbind(ATAC_peak_MTA2_list,ATAC_peak_MTA2)
  }
}
write.table(ATAC_peak_MTA2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MTA2_list.txt",sep="\t")

ATAC_peak_MTA2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MTA2_list.txt",sep="\t")

MTA2_features=list("MTA2" = ATAC_peak_MTA2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MTA2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##WHSC1
WHSC1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF862UUR_WHSC1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(WHSC1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF862UUR_WHSC1.txt",sep = "\t")

WHSC1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF862UUR_WHSC1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(WHSC1)){
  x=paste(WHSC1[i,1],WHSC1[i,2],WHSC1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

WHSC1_features = list(list$V1)

ATAC_peak_WHSC1_list = data.frame()
for(i in 1:nrow(WHSC1)){
  ATAC_peak_WHSC1<- ATAC_peak[ATAC_peak$x1 ==WHSC1[i,1] & ((WHSC1[i,2] > ATAC_peak$start &WHSC1[i,2] < ATAC_peak$end) | (WHSC1[i,3] > ATAC_peak$start &WHSC1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_WHSC1) != 0){
    ATAC_peak_WHSC1$V5 <-WHSC1[i,1]
    ATAC_peak_WHSC1$V6 <- WHSC1[i,2]
    ATAC_peak_WHSC1$V7 <- WHSC1[i,3]
    ATAC_peak_WHSC1$V8 <- WHSC1[i,6]
    ATAC_peak_WHSC1$V9 <- WHSC1[i,7]
    ATAC_peak_WHSC1_list = rbind(ATAC_peak_WHSC1_list,ATAC_peak_WHSC1)
  }
}
write.table(ATAC_peak_WHSC1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_WHSC1_list.txt",sep="\t")

ATAC_peak_WHSC1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_WHSC1_list.txt",sep="\t")

WHSC1_features=list("WHSC1" = ATAC_peak_WHSC1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = WHSC1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##GATAD2A
GATAD2A= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF863NVQ_GATAD2A.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GATAD2A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF863NVQ_GATAD2A.txt",sep = "\t")

GATAD2A= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF863NVQ_GATAD2A.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GATAD2A)){
  x=paste(GATAD2A[i,1],GATAD2A[i,2],GATAD2A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GATAD2A_features = list(list$V1)

ATAC_peak_GATAD2A_list = data.frame()
for(i in 1:nrow(GATAD2A)){
  ATAC_peak_GATAD2A<- ATAC_peak[ATAC_peak$x1 ==GATAD2A[i,1] & ((GATAD2A[i,2] > ATAC_peak$start &GATAD2A[i,2] < ATAC_peak$end) | (GATAD2A[i,3] > ATAC_peak$start &GATAD2A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GATAD2A) != 0){
    ATAC_peak_GATAD2A$V5 <-GATAD2A[i,1]
    ATAC_peak_GATAD2A$V6 <- GATAD2A[i,2]
    ATAC_peak_GATAD2A$V7 <- GATAD2A[i,3]
    ATAC_peak_GATAD2A$V8 <- GATAD2A[i,6]
    ATAC_peak_GATAD2A$V9 <- GATAD2A[i,7]
    ATAC_peak_GATAD2A_list = rbind(ATAC_peak_GATAD2A_list,ATAC_peak_GATAD2A)
  }
}
write.table(ATAC_peak_GATAD2A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GATAD2A_list.txt",sep="\t")

ATAC_peak_GATAD2A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GATAD2A_list.txt",sep="\t")

GATAD2A_features=list("GATAD2A" = ATAC_peak_GATAD2A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GATAD2A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##TRIM24
TRIM24= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF864JGI_TRIM24.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TRIM24,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF864JGI_TRIM24.txt",sep = "\t")

TRIM24= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF864JGI_TRIM24.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TRIM24)){
  x=paste(TRIM24[i,1],TRIM24[i,2],TRIM24[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TRIM24_features = list(list$V1)

ATAC_peak_TRIM24_list = data.frame()
for(i in 1:nrow(TRIM24)){
  ATAC_peak_TRIM24<- ATAC_peak[ATAC_peak$x1 ==TRIM24[i,1] & ((TRIM24[i,2] > ATAC_peak$start &TRIM24[i,2] < ATAC_peak$end) | (TRIM24[i,3] > ATAC_peak$start &TRIM24[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TRIM24) != 0){
    ATAC_peak_TRIM24$V5 <-TRIM24[i,1]
    ATAC_peak_TRIM24$V6 <- TRIM24[i,2]
    ATAC_peak_TRIM24$V7 <- TRIM24[i,3]
    ATAC_peak_TRIM24$V8 <- TRIM24[i,6]
    ATAC_peak_TRIM24$V9 <- TRIM24[i,7]
    ATAC_peak_TRIM24_list = rbind(ATAC_peak_TRIM24_list,ATAC_peak_TRIM24)
  }
}
write.table(ATAC_peak_TRIM24_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TRIM24_list.txt",sep="\t")

ATAC_peak_TRIM24_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TRIM24_list.txt",sep="\t")

TRIM24_features=list("TRIM24" = ATAC_peak_TRIM24_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TRIM24_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##MCM2
MCM2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF869KPH_MCM2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MCM2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF869KPH_MCM2.txt",sep = "\t")

MCM2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF869KPH_MCM2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MCM2)){
  x=paste(MCM2[i,1],MCM2[i,2],MCM2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MCM2_features = list(list$V1)

ATAC_peak_MCM2_list = data.frame()
for(i in 1:nrow(MCM2)){
  ATAC_peak_MCM2<- ATAC_peak[ATAC_peak$x1 ==MCM2[i,1] & ((MCM2[i,2] > ATAC_peak$start &MCM2[i,2] < ATAC_peak$end) | (MCM2[i,3] > ATAC_peak$start &MCM2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MCM2) != 0){
    ATAC_peak_MCM2$V5 <-MCM2[i,1]
    ATAC_peak_MCM2$V6 <- MCM2[i,2]
    ATAC_peak_MCM2$V7 <- MCM2[i,3]
    ATAC_peak_MCM2$V8 <- MCM2[i,6]
    ATAC_peak_MCM2$V9 <- MCM2[i,7]
    ATAC_peak_MCM2_list = rbind(ATAC_peak_MCM2_list,ATAC_peak_MCM2)
  }
}
write.table(ATAC_peak_MCM2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MCM2_list.txt",sep="\t")

ATAC_peak_MCM2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MCM2_list.txt",sep="\t")

MCM2_features=list("MCM2" = ATAC_peak_MCM2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MCM2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##ARID2
ARID2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF879ZMI_ARID2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ARID2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF879ZMI_ARID2.txt",sep = "\t")

ARID2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF879ZMI_ARID2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ARID2)){
  x=paste(ARID2[i,1],ARID2[i,2],ARID2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ARID2_features = list(list$V1)

ATAC_peak_ARID2_list = data.frame()
for(i in 1:nrow(ARID2)){
  ATAC_peak_ARID2<- ATAC_peak[ATAC_peak$x1 ==ARID2[i,1] & ((ARID2[i,2] > ATAC_peak$start &ARID2[i,2] < ATAC_peak$end) | (ARID2[i,3] > ATAC_peak$start &ARID2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ARID2) != 0){
    ATAC_peak_ARID2$V5 <-ARID2[i,1]
    ATAC_peak_ARID2$V6 <- ARID2[i,2]
    ATAC_peak_ARID2$V7 <- ARID2[i,3]
    ATAC_peak_ARID2$V8 <- ARID2[i,6]
    ATAC_peak_ARID2$V9 <- ARID2[i,7]
    ATAC_peak_ARID2_list = rbind(ATAC_peak_ARID2_list,ATAC_peak_ARID2)
  }
}
write.table(ATAC_peak_ARID2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARID2_list.txt",sep="\t")

ATAC_peak_ARID2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARID2_list.txt",sep="\t")

ARID2_features=list("ARID2" = ATAC_peak_ARID2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ARID2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##PRPF4
PRPF4= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF886UMM_PRPF4.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PRPF4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF886UMM_PRPF4.txt",sep = "\t")

PRPF4= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF886UMM_PRPF4.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PRPF4)){
  x=paste(PRPF4[i,1],PRPF4[i,2],PRPF4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PRPF4_features = list(list$V1)

ATAC_peak_PRPF4_list = data.frame()
for(i in 1:nrow(PRPF4)){
  ATAC_peak_PRPF4<- ATAC_peak[ATAC_peak$x1 ==PRPF4[i,1] & ((PRPF4[i,2] > ATAC_peak$start &PRPF4[i,2] < ATAC_peak$end) | (PRPF4[i,3] > ATAC_peak$start &PRPF4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PRPF4) != 0){
    ATAC_peak_PRPF4$V5 <-PRPF4[i,1]
    ATAC_peak_PRPF4$V6 <- PRPF4[i,2]
    ATAC_peak_PRPF4$V7 <- PRPF4[i,3]
    ATAC_peak_PRPF4$V8 <- PRPF4[i,6]
    ATAC_peak_PRPF4$V9 <- PRPF4[i,7]
    ATAC_peak_PRPF4_list = rbind(ATAC_peak_PRPF4_list,ATAC_peak_PRPF4)
  }
}
write.table(ATAC_peak_PRPF4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PRPF4_list.txt",sep="\t")

ATAC_peak_PRPF4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PRPF4_list.txt",sep="\t")

PRPF4_features=list("PRPF4" = ATAC_peak_PRPF4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PRPF4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##BACH1
BACH1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF887JWQ_BACH1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BACH1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF887JWQ_BACH1.txt",sep = "\t")

BACH1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF887JWQ_BACH1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BACH1)){
  x=paste(BACH1[i,1],BACH1[i,2],BACH1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BACH1_features = list(list$V1)

ATAC_peak_BACH1_list = data.frame()
for(i in 1:nrow(BACH1)){
  ATAC_peak_BACH1<- ATAC_peak[ATAC_peak$x1 ==BACH1[i,1] & ((BACH1[i,2] > ATAC_peak$start &BACH1[i,2] < ATAC_peak$end) | (BACH1[i,3] > ATAC_peak$start &BACH1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BACH1) != 0){
    ATAC_peak_BACH1$V5 <-BACH1[i,1]
    ATAC_peak_BACH1$V6 <- BACH1[i,2]
    ATAC_peak_BACH1$V7 <- BACH1[i,3]
    ATAC_peak_BACH1$V8 <- BACH1[i,6]
    ATAC_peak_BACH1$V9 <- BACH1[i,7]
    ATAC_peak_BACH1_list = rbind(ATAC_peak_BACH1_list,ATAC_peak_BACH1)
  }
}
write.table(ATAC_peak_BACH1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BACH1_list.txt",sep="\t")

ATAC_peak_BACH1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BACH1_list.txt",sep="\t")

BACH1_features=list("BACH1" = ATAC_peak_BACH1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BACH1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##GABPA
GABPA= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF889XXJ_GABPA.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GABPA,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF889XXJ_GABPA.txt",sep = "\t")

GABPA= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF889XXJ_GABPA.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GABPA)){
  x=paste(GABPA[i,1],GABPA[i,2],GABPA[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GABPA_features = list(list$V1)

ATAC_peak_GABPA_list = data.frame()
for(i in 1:nrow(GABPA)){
  ATAC_peak_GABPA<- ATAC_peak[ATAC_peak$x1 ==GABPA[i,1] & ((GABPA[i,2] > ATAC_peak$start &GABPA[i,2] < ATAC_peak$end) | (GABPA[i,3] > ATAC_peak$start &GABPA[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GABPA) != 0){
    ATAC_peak_GABPA$V5 <-GABPA[i,1]
    ATAC_peak_GABPA$V6 <- GABPA[i,2]
    ATAC_peak_GABPA$V7 <- GABPA[i,3]
    ATAC_peak_GABPA$V8 <- GABPA[i,6]
    ATAC_peak_GABPA$V9 <- GABPA[i,7]
    ATAC_peak_GABPA_list = rbind(ATAC_peak_GABPA_list,ATAC_peak_GABPA)
  }
}
write.table(ATAC_peak_GABPA_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GABPA_list.txt",sep="\t")

ATAC_peak_GABPA_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GABPA_list.txt",sep="\t")

GABPA_features=list("GABPA" = ATAC_peak_GABPA_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GABPA_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##GTF2F1
GTF2F1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF895ZCW_GTF2F1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(GTF2F1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF895ZCW_GTF2F1.txt",sep = "\t")

GTF2F1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF895ZCW_GTF2F1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(GTF2F1)){
  x=paste(GTF2F1[i,1],GTF2F1[i,2],GTF2F1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

GTF2F1_features = list(list$V1)

ATAC_peak_GTF2F1_list = data.frame()
for(i in 1:nrow(GTF2F1)){
  ATAC_peak_GTF2F1<- ATAC_peak[ATAC_peak$x1 ==GTF2F1[i,1] & ((GTF2F1[i,2] > ATAC_peak$start &GTF2F1[i,2] < ATAC_peak$end) | (GTF2F1[i,3] > ATAC_peak$start &GTF2F1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_GTF2F1) != 0){
    ATAC_peak_GTF2F1$V5 <-GTF2F1[i,1]
    ATAC_peak_GTF2F1$V6 <- GTF2F1[i,2]
    ATAC_peak_GTF2F1$V7 <- GTF2F1[i,3]
    ATAC_peak_GTF2F1$V8 <- GTF2F1[i,6]
    ATAC_peak_GTF2F1$V9 <- GTF2F1[i,7]
    ATAC_peak_GTF2F1_list = rbind(ATAC_peak_GTF2F1_list,ATAC_peak_GTF2F1)
  }
}
write.table(ATAC_peak_GTF2F1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF2F1_list.txt",sep="\t")

ATAC_peak_GTF2F1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_GTF2F1_list.txt",sep="\t")

GTF2F1_features=list("GTF2F1" = ATAC_peak_GTF2F1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = GTF2F1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZNF316
ZNF316= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF907ICG_ZNF316.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZNF316,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF907ICG_ZNF316.txt",sep = "\t")

ZNF316= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF907ICG_ZNF316.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZNF316)){
  x=paste(ZNF316[i,1],ZNF316[i,2],ZNF316[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZNF316_features = list(list$V1)

ATAC_peak_ZNF316_list = data.frame()
for(i in 1:nrow(ZNF316)){
  ATAC_peak_ZNF316<- ATAC_peak[ATAC_peak$x1 ==ZNF316[i,1] & ((ZNF316[i,2] > ATAC_peak$start &ZNF316[i,2] < ATAC_peak$end) | (ZNF316[i,3] > ATAC_peak$start &ZNF316[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZNF316) != 0){
    ATAC_peak_ZNF316$V5 <-ZNF316[i,1]
    ATAC_peak_ZNF316$V6 <- ZNF316[i,2]
    ATAC_peak_ZNF316$V7 <- ZNF316[i,3]
    ATAC_peak_ZNF316$V8 <- ZNF316[i,6]
    ATAC_peak_ZNF316$V9 <- ZNF316[i,7]
    ATAC_peak_ZNF316_list = rbind(ATAC_peak_ZNF316_list,ATAC_peak_ZNF316)
  }
}
write.table(ATAC_peak_ZNF316_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF316_list.txt",sep="\t")

ATAC_peak_ZNF316_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZNF316_list.txt",sep="\t")

ZNF316_features=list("ZNF316" = ATAC_peak_ZNF316_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZNF316_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RLF
RLF= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF910PKV_RLF.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RLF,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF910PKV_RLF.txt",sep = "\t")

RLF= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF910PKV_RLF.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RLF)){
  x=paste(RLF[i,1],RLF[i,2],RLF[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RLF_features = list(list$V1)

ATAC_peak_RLF_list = data.frame()
for(i in 1:nrow(RLF)){
  ATAC_peak_RLF<- ATAC_peak[ATAC_peak$x1 ==RLF[i,1] & ((RLF[i,2] > ATAC_peak$start &RLF[i,2] < ATAC_peak$end) | (RLF[i,3] > ATAC_peak$start &RLF[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RLF) != 0){
    ATAC_peak_RLF$V5 <-RLF[i,1]
    ATAC_peak_RLF$V6 <- RLF[i,2]
    ATAC_peak_RLF$V7 <- RLF[i,3]
    ATAC_peak_RLF$V8 <- RLF[i,6]
    ATAC_peak_RLF$V9 <- RLF[i,7]
    ATAC_peak_RLF_list = rbind(ATAC_peak_RLF_list,ATAC_peak_RLF)
  }
}
write.table(ATAC_peak_RLF_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RLF_list.txt",sep="\t")

ATAC_peak_RLF_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RLF_list.txt",sep="\t")

RLF_features=list("RLF" = ATAC_peak_RLF_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RLF_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##BRCA1
BRCA1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF914PHM_BRCA1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BRCA1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF914PHM_BRCA1_age.txt",sep = "\t")

BRCA1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF914PHM_BRCA1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BRCA1)){
  x=paste(BRCA1[i,1],BRCA1[i,2],BRCA1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BRCA1_features = list(list$V1)

ATAC_peak_BRCA1_list = data.frame()
for(i in 1:nrow(BRCA1)){
  ATAC_peak_BRCA1<- ATAC_peak[ATAC_peak$x1 ==BRCA1[i,1] & ((BRCA1[i,2] > ATAC_peak$start &BRCA1[i,2] < ATAC_peak$end) | (BRCA1[i,3] > ATAC_peak$start &BRCA1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BRCA1) != 0){
    ATAC_peak_BRCA1$V5 <-BRCA1[i,1]
    ATAC_peak_BRCA1$V6 <- BRCA1[i,2]
    ATAC_peak_BRCA1$V7 <- BRCA1[i,3]
    ATAC_peak_BRCA1$V8 <- BRCA1[i,6]
    ATAC_peak_BRCA1$V9 <- BRCA1[i,7]
    ATAC_peak_BRCA1_list = rbind(ATAC_peak_BRCA1_list,ATAC_peak_BRCA1)
  }
}
write.table(ATAC_peak_BRCA1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRCA1_list.txt",sep="\t")

ATAC_peak_BRCA1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BRCA1_list.txt",sep="\t")

BRCA1_features=list("BRCA1" = ATAC_peak_BRCA1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BRCA1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##DPF2
DPF2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF917LAX_DPF2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(DPF2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF917LAX_DPF2.txt",sep = "\t")

DPF2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF917LAX_DPF2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(DPF2)){
  x=paste(DPF2[i,1],DPF2[i,2],DPF2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

DPF2_features = list(list$V1)

ATAC_peak_DPF2_list = data.frame()
for(i in 1:nrow(DPF2)){
  ATAC_peak_DPF2<- ATAC_peak[ATAC_peak$x1 ==DPF2[i,1] & ((DPF2[i,2] > ATAC_peak$start &DPF2[i,2] < ATAC_peak$end) | (DPF2[i,3] > ATAC_peak$start &DPF2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_DPF2) != 0){
    ATAC_peak_DPF2$V5 <-DPF2[i,1]
    ATAC_peak_DPF2$V6 <- DPF2[i,2]
    ATAC_peak_DPF2$V7 <- DPF2[i,3]
    ATAC_peak_DPF2$V8 <- DPF2[i,6]
    ATAC_peak_DPF2$V9 <- DPF2[i,7]
    ATAC_peak_DPF2_list = rbind(ATAC_peak_DPF2_list,ATAC_peak_DPF2)
  }
}
write.table(ATAC_peak_DPF2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DPF2_list.txt",sep="\t")

ATAC_peak_DPF2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DPF2_list.txt",sep="\t")

DPF2_features=list("DPF2" = ATAC_peak_DPF2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = DPF2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##NR0B1
NR0B1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF931LUC_NR0B1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(NR0B1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF931LUC_NR0B1.txt",sep = "\t")

NR0B1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF931LUC_NR0B1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NR0B1)){
  x=paste(NR0B1[i,1],NR0B1[i,2],NR0B1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NR0B1_features = list(list$V1)

ATAC_peak_NR0B1_list = data.frame()
for(i in 1:nrow(NR0B1)){
  ATAC_peak_NR0B1<- ATAC_peak[ATAC_peak$x1 ==NR0B1[i,1] & ((NR0B1[i,2] > ATAC_peak$start &NR0B1[i,2] < ATAC_peak$end) | (NR0B1[i,3] > ATAC_peak$start &NR0B1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NR0B1) != 0){
    ATAC_peak_NR0B1$V5 <-NR0B1[i,1]
    ATAC_peak_NR0B1$V6 <- NR0B1[i,2]
    ATAC_peak_NR0B1$V7 <- NR0B1[i,3]
    ATAC_peak_NR0B1$V8 <- NR0B1[i,6]
    ATAC_peak_NR0B1$V9 <- NR0B1[i,7]
    ATAC_peak_NR0B1_list = rbind(ATAC_peak_NR0B1_list,ATAC_peak_NR0B1)
  }
}
write.table(ATAC_peak_NR0B1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NR0B1_list.txt",sep="\t")

ATAC_peak_NR0B1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NR0B1_list.txt",sep="\t")

NR0B1_features=list("NR0B1" = ATAC_peak_NR0B1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NR0B1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##KDM5B
KDM5B= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF933VQA_KDM5B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(KDM5B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF933VQA_KDM5B.txt",sep = "\t")

KDM5B= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF933VQA_KDM5B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(KDM5B)){
  x=paste(KDM5B[i,1],KDM5B[i,2],KDM5B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

KDM5B_features = list(list$V1)

ATAC_peak_KDM5B_list = data.frame()
for(i in 1:nrow(KDM5B)){
  ATAC_peak_KDM5B<- ATAC_peak[ATAC_peak$x1 ==KDM5B[i,1] & ((KDM5B[i,2] > ATAC_peak$start &KDM5B[i,2] < ATAC_peak$end) | (KDM5B[i,3] > ATAC_peak$start &KDM5B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_KDM5B) != 0){
    ATAC_peak_KDM5B$V5 <-KDM5B[i,1]
    ATAC_peak_KDM5B$V6 <- KDM5B[i,2]
    ATAC_peak_KDM5B$V7 <- KDM5B[i,3]
    ATAC_peak_KDM5B$V8 <- KDM5B[i,6]
    ATAC_peak_KDM5B$V9 <- KDM5B[i,7]
    ATAC_peak_KDM5B_list = rbind(ATAC_peak_KDM5B_list,ATAC_peak_KDM5B)
  }
}
write.table(ATAC_peak_KDM5B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KDM5B_list.txt",sep="\t")

ATAC_peak_KDM5B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KDM5B_list.txt",sep="\t")

KDM5B_features=list("KDM5B" = ATAC_peak_KDM5B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = KDM5B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##POLR2B
POLR2B= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF935YOF_POLR2B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(POLR2B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF935YOF_POLR2B.txt",sep = "\t")

POLR2B= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF935YOF_POLR2B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(POLR2B)){
  x=paste(POLR2B[i,1],POLR2B[i,2],POLR2B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

POLR2B_features = list(list$V1)

ATAC_peak_POLR2B_list = data.frame()
for(i in 1:nrow(POLR2B)){
  ATAC_peak_POLR2B<- ATAC_peak[ATAC_peak$x1 ==POLR2B[i,1] & ((POLR2B[i,2] > ATAC_peak$start &POLR2B[i,2] < ATAC_peak$end) | (POLR2B[i,3] > ATAC_peak$start &POLR2B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_POLR2B) != 0){
    ATAC_peak_POLR2B$V5 <-POLR2B[i,1]
    ATAC_peak_POLR2B$V6 <- POLR2B[i,2]
    ATAC_peak_POLR2B$V7 <- POLR2B[i,3]
    ATAC_peak_POLR2B$V8 <- POLR2B[i,6]
    ATAC_peak_POLR2B$V9 <- POLR2B[i,7]
    ATAC_peak_POLR2B_list = rbind(ATAC_peak_POLR2B_list,ATAC_peak_POLR2B)
  }
}
write.table(ATAC_peak_POLR2B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLR2B_list.txt",sep="\t")

ATAC_peak_POLR2B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLR2B_list.txt",sep="\t")

POLR2B_features=list("POLR2B" = ATAC_peak_POLR2B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = POLR2B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##COPS2
COPS2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF942JZR_COPS2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(COPS2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF942JZR_COPS2.txt",sep = "\t")

COPS2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF942JZR_COPS2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(COPS2)){
  x=paste(COPS2[i,1],COPS2[i,2],COPS2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

COPS2_features = list(list$V1)

ATAC_peak_COPS2_list = data.frame()
for(i in 1:nrow(COPS2)){
  ATAC_peak_COPS2<- ATAC_peak[ATAC_peak$x1 ==COPS2[i,1] & ((COPS2[i,2] > ATAC_peak$start &COPS2[i,2] < ATAC_peak$end) | (COPS2[i,3] > ATAC_peak$start &COPS2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_COPS2) != 0){
    ATAC_peak_COPS2$V5 <-COPS2[i,1]
    ATAC_peak_COPS2$V6 <- COPS2[i,2]
    ATAC_peak_COPS2$V7 <- COPS2[i,3]
    ATAC_peak_COPS2$V8 <- COPS2[i,6]
    ATAC_peak_COPS2$V9 <- COPS2[i,7]
    ATAC_peak_COPS2_list = rbind(ATAC_peak_COPS2_list,ATAC_peak_COPS2)
  }
}
write.table(ATAC_peak_COPS2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_COPS2_list.txt",sep="\t")

ATAC_peak_COPS2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_COPS2_list.txt",sep="\t")

COPS2_features=list("COPS2" = ATAC_peak_COPS2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = COPS2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##SUZ12
SUZ12= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF943SUL_SUZ12.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SUZ12,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF943SUL_SUZ12.txt",sep = "\t")

SUZ12= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF943SUL_SUZ12.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SUZ12)){
  x=paste(SUZ12[i,1],SUZ12[i,2],SUZ12[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SUZ12_features = list(list$V1)

ATAC_peak_SUZ12_list = data.frame()
for(i in 1:nrow(SUZ12)){
  ATAC_peak_SUZ12<- ATAC_peak[ATAC_peak$x1 ==SUZ12[i,1] & ((SUZ12[i,2] > ATAC_peak$start &SUZ12[i,2] < ATAC_peak$end) | (SUZ12[i,3] > ATAC_peak$start &SUZ12[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SUZ12) != 0){
    ATAC_peak_SUZ12$V5 <-SUZ12[i,1]
    ATAC_peak_SUZ12$V6 <- SUZ12[i,2]
    ATAC_peak_SUZ12$V7 <- SUZ12[i,3]
    ATAC_peak_SUZ12$V8 <- SUZ12[i,6]
    ATAC_peak_SUZ12$V9 <- SUZ12[i,7]
    ATAC_peak_SUZ12_list = rbind(ATAC_peak_SUZ12_list,ATAC_peak_SUZ12)
  }
}
write.table(ATAC_peak_SUZ12_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SUZ12_list.txt",sep="\t")

ATAC_peak_SUZ12_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SUZ12_list.txt",sep="\t")

SUZ12_features=list("SUZ12" = ATAC_peak_SUZ12_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SUZ12_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZKSCAN1
ZKSCAN1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF948WSE_ZKSCAN1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZKSCAN1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF948WSE_ZKSCAN1.txt",sep = "\t")

ZKSCAN1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF948WSE_ZKSCAN1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZKSCAN1)){
  x=paste(ZKSCAN1[i,1],ZKSCAN1[i,2],ZKSCAN1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZKSCAN1_features = list(list$V1)

ATAC_peak_ZKSCAN1_list = data.frame()
for(i in 1:nrow(ZKSCAN1)){
  ATAC_peak_ZKSCAN1<- ATAC_peak[ATAC_peak$x1 ==ZKSCAN1[i,1] & ((ZKSCAN1[i,2] > ATAC_peak$start &ZKSCAN1[i,2] < ATAC_peak$end) | (ZKSCAN1[i,3] > ATAC_peak$start &ZKSCAN1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZKSCAN1) != 0){
    ATAC_peak_ZKSCAN1$V5 <-ZKSCAN1[i,1]
    ATAC_peak_ZKSCAN1$V6 <- ZKSCAN1[i,2]
    ATAC_peak_ZKSCAN1$V7 <- ZKSCAN1[i,3]
    ATAC_peak_ZKSCAN1$V8 <- ZKSCAN1[i,6]
    ATAC_peak_ZKSCAN1$V9 <- ZKSCAN1[i,7]
    ATAC_peak_ZKSCAN1_list = rbind(ATAC_peak_ZKSCAN1_list,ATAC_peak_ZKSCAN1)
  }
}
write.table(ATAC_peak_ZKSCAN1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZKSCAN1_list.txt",sep="\t")

ATAC_peak_ZKSCAN1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZKSCAN1_list.txt",sep="\t")

ZKSCAN1_features=list("ZKSCAN1" = ATAC_peak_ZKSCAN1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZKSCAN1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##MCM5
MCM5= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF949IHL_MCM5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(MCM5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF949IHL_MCM5.txt",sep = "\t")

MCM5= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF949IHL_MCM5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MCM5)){
  x=paste(MCM5[i,1],MCM5[i,2],MCM5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MCM5_features = list(list$V1)

ATAC_peak_MCM5_list = data.frame()
for(i in 1:nrow(MCM5)){
  ATAC_peak_MCM5<- ATAC_peak[ATAC_peak$x1 ==MCM5[i,1] & ((MCM5[i,2] > ATAC_peak$start &MCM5[i,2] < ATAC_peak$end) | (MCM5[i,3] > ATAC_peak$start &MCM5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MCM5) != 0){
    ATAC_peak_MCM5$V5 <-MCM5[i,1]
    ATAC_peak_MCM5$V6 <- MCM5[i,2]
    ATAC_peak_MCM5$V7 <- MCM5[i,3]
    ATAC_peak_MCM5$V8 <- MCM5[i,6]
    ATAC_peak_MCM5$V9 <- MCM5[i,7]
    ATAC_peak_MCM5_list = rbind(ATAC_peak_MCM5_list,ATAC_peak_MCM5)
  }
}
write.table(ATAC_peak_MCM5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MCM5_list.txt",sep="\t")

ATAC_peak_MCM5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MCM5_list.txt",sep="\t")

MCM5_features=list("MCM5" = ATAC_peak_MCM5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MCM5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##SRSF3
SRSF3= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF959YOR_SRSF3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SRSF3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF959YOR_SRSF3.txt",sep = "\t")

SRSF3= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF959YOR_SRSF3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SRSF3)){
  x=paste(SRSF3[i,1],SRSF3[i,2],SRSF3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SRSF3_features = list(list$V1)

ATAC_peak_SRSF3_list = data.frame()
for(i in 1:nrow(SRSF3)){
  ATAC_peak_SRSF3<- ATAC_peak[ATAC_peak$x1 ==SRSF3[i,1] & ((SRSF3[i,2] > ATAC_peak$start &SRSF3[i,2] < ATAC_peak$end) | (SRSF3[i,3] > ATAC_peak$start &SRSF3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SRSF3) != 0){
    ATAC_peak_SRSF3$V5 <-SRSF3[i,1]
    ATAC_peak_SRSF3$V6 <- SRSF3[i,2]
    ATAC_peak_SRSF3$V7 <- SRSF3[i,3]
    ATAC_peak_SRSF3$V8 <- SRSF3[i,6]
    ATAC_peak_SRSF3$V9 <- SRSF3[i,7]
    ATAC_peak_SRSF3_list = rbind(ATAC_peak_SRSF3_list,ATAC_peak_SRSF3)
  }
}
write.table(ATAC_peak_SRSF3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SRSF3_list.txt",sep="\t")

ATAC_peak_SRSF3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SRSF3_list.txt",sep="\t")

SRSF3_features=list("SRSF3" = ATAC_peak_SRSF3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SRSF3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ILF3
ILF3= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF965MIZ_ILF3.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ILF3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF965MIZ_ILF3.txt",sep = "\t")

ILF3= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF965MIZ_ILF3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ILF3)){
  x=paste(ILF3[i,1],ILF3[i,2],ILF3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ILF3_features = list(list$V1)

ATAC_peak_ILF3_list = data.frame()
for(i in 1:nrow(ILF3)){
  ATAC_peak_ILF3<- ATAC_peak[ATAC_peak$x1 ==ILF3[i,1] & ((ILF3[i,2] > ATAC_peak$start &ILF3[i,2] < ATAC_peak$end) | (ILF3[i,3] > ATAC_peak$start &ILF3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ILF3) != 0){
    ATAC_peak_ILF3$V5 <-ILF3[i,1]
    ATAC_peak_ILF3$V6 <- ILF3[i,2]
    ATAC_peak_ILF3$V7 <- ILF3[i,3]
    ATAC_peak_ILF3$V8 <- ILF3[i,6]
    ATAC_peak_ILF3$V9 <- ILF3[i,7]
    ATAC_peak_ILF3_list = rbind(ATAC_peak_ILF3_list,ATAC_peak_ILF3)
  }
}
write.table(ATAC_peak_ILF3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ILF3_list.txt",sep="\t")

ATAC_peak_ILF3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ILF3_list.txt",sep="\t")

ILF3_features=list("ILF3" = ATAC_peak_ILF3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ILF3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##BCLAF1
BCLAF1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF966JHD_BCLAF1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(BCLAF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF966JHD_BCLAF1.txt",sep = "\t")

BCLAF1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF966JHD_BCLAF1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BCLAF1)){
  x=paste(BCLAF1[i,1],BCLAF1[i,2],BCLAF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BCLAF1_features = list(list$V1)

ATAC_peak_BCLAF1_list = data.frame()
for(i in 1:nrow(BCLAF1)){
  ATAC_peak_BCLAF1<- ATAC_peak[ATAC_peak$x1 ==BCLAF1[i,1] & ((BCLAF1[i,2] > ATAC_peak$start &BCLAF1[i,2] < ATAC_peak$end) | (BCLAF1[i,3] > ATAC_peak$start &BCLAF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BCLAF1) != 0){
    ATAC_peak_BCLAF1$V5 <-BCLAF1[i,1]
    ATAC_peak_BCLAF1$V6 <- BCLAF1[i,2]
    ATAC_peak_BCLAF1$V7 <- BCLAF1[i,3]
    ATAC_peak_BCLAF1$V8 <- BCLAF1[i,6]
    ATAC_peak_BCLAF1$V9 <- BCLAF1[i,7]
    ATAC_peak_BCLAF1_list = rbind(ATAC_peak_BCLAF1_list,ATAC_peak_BCLAF1)
  }
}
write.table(ATAC_peak_BCLAF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BCLAF1_list.txt",sep="\t")

ATAC_peak_BCLAF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BCLAF1_list.txt",sep="\t")

BCLAF1_features=list("BCLAF1" = ATAC_peak_BCLAF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BCLAF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ZKSCAN8
ZKSCAN8= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF967XKC_ZKSCAN8.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZKSCAN8,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF967XKC_ZKSCAN8.txt",sep = "\t")

ZKSCAN8= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF967XKC_ZKSCAN8.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZKSCAN8)){
  x=paste(ZKSCAN8[i,1],ZKSCAN8[i,2],ZKSCAN8[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZKSCAN8_features = list(list$V1)

ATAC_peak_ZKSCAN8_list = data.frame()
for(i in 1:nrow(ZKSCAN8)){
  ATAC_peak_ZKSCAN8<- ATAC_peak[ATAC_peak$x1 ==ZKSCAN8[i,1] & ((ZKSCAN8[i,2] > ATAC_peak$start &ZKSCAN8[i,2] < ATAC_peak$end) | (ZKSCAN8[i,3] > ATAC_peak$start &ZKSCAN8[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZKSCAN8) != 0){
    ATAC_peak_ZKSCAN8$V5 <-ZKSCAN8[i,1]
    ATAC_peak_ZKSCAN8$V6 <- ZKSCAN8[i,2]
    ATAC_peak_ZKSCAN8$V7 <- ZKSCAN8[i,3]
    ATAC_peak_ZKSCAN8$V8 <- ZKSCAN8[i,6]
    ATAC_peak_ZKSCAN8$V9 <- ZKSCAN8[i,7]
    ATAC_peak_ZKSCAN8_list = rbind(ATAC_peak_ZKSCAN8_list,ATAC_peak_ZKSCAN8)
  }
}
write.table(ATAC_peak_ZKSCAN8_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZKSCAN8_list.txt",sep="\t")

ATAC_peak_ZKSCAN8_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZKSCAN8_list.txt",sep="\t")

ZKSCAN8_features=list("ZKSCAN8" = ATAC_peak_ZKSCAN8_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZKSCAN8_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ARID1B
ARID1B= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF967XTM_ARID1B.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ARID1B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF967XTM_ARID1B.txt",sep = "\t")

ARID1B= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF967XTM_ARID1B.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ARID1B)){
  x=paste(ARID1B[i,1],ARID1B[i,2],ARID1B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ARID1B_features = list(list$V1)

ATAC_peak_ARID1B_list = data.frame()
for(i in 1:nrow(ARID1B)){
  ATAC_peak_ARID1B<- ATAC_peak[ATAC_peak$x1 ==ARID1B[i,1] & ((ARID1B[i,2] > ATAC_peak$start &ARID1B[i,2] < ATAC_peak$end) | (ARID1B[i,3] > ATAC_peak$start &ARID1B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ARID1B) != 0){
    ATAC_peak_ARID1B$V5 <-ARID1B[i,1]
    ATAC_peak_ARID1B$V6 <- ARID1B[i,2]
    ATAC_peak_ARID1B$V7 <- ARID1B[i,3]
    ATAC_peak_ARID1B$V8 <- ARID1B[i,6]
    ATAC_peak_ARID1B$V9 <- ARID1B[i,7]
    ATAC_peak_ARID1B_list = rbind(ATAC_peak_ARID1B_list,ATAC_peak_ARID1B)
  }
}
write.table(ATAC_peak_ARID1B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARID1B_list.txt",sep="\t")

ATAC_peak_ARID1B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARID1B_list.txt",sep="\t")

ARID1B_features=list("ARID1B" = ATAC_peak_ARID1B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ARID1B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##TCF12
TCF12= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF970YIG_TCF12.bed.gz",extraCols=extraCols_narrowPeak)
write.table(TCF12,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF970YIG_TCF12.txt",sep = "\t")

TCF12= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF970YIG_TCF12.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TCF12)){
  x=paste(TCF12[i,1],TCF12[i,2],TCF12[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TCF12_features = list(list$V1)

ATAC_peak_TCF12_list = data.frame()
for(i in 1:nrow(TCF12)){
  ATAC_peak_TCF12<- ATAC_peak[ATAC_peak$x1 ==TCF12[i,1] & ((TCF12[i,2] > ATAC_peak$start &TCF12[i,2] < ATAC_peak$end) | (TCF12[i,3] > ATAC_peak$start &TCF12[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TCF12) != 0){
    ATAC_peak_TCF12$V5 <-TCF12[i,1]
    ATAC_peak_TCF12$V6 <- TCF12[i,2]
    ATAC_peak_TCF12$V7 <- TCF12[i,3]
    ATAC_peak_TCF12$V8 <- TCF12[i,6]
    ATAC_peak_TCF12$V9 <- TCF12[i,7]
    ATAC_peak_TCF12_list = rbind(ATAC_peak_TCF12_list,ATAC_peak_TCF12)
  }
}
write.table(ATAC_peak_TCF12_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TCF12_list.txt",sep="\t")

ATAC_peak_TCF12_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TCF12_list.txt",sep="\t")

TCF12_features=list("TCF12" = ATAC_peak_TCF12_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TCF12_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##KHSRP
KHSRP= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF973UWZ_KHSRP.bed.gz",extraCols=extraCols_narrowPeak)
write.table(KHSRP,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF970YIG_KHSRP.txt",sep = "\t")

KHSRP= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF970YIG_KHSRP.txt",sep = "\t")

list = vector()
for (i in  1:nrow(KHSRP)){
  x=paste(KHSRP[i,1],KHSRP[i,2],KHSRP[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

KHSRP_features = list(list$V1)

ATAC_peak_KHSRP_list = data.frame()
for(i in 1:nrow(KHSRP)){
  ATAC_peak_KHSRP<- ATAC_peak[ATAC_peak$x1 ==KHSRP[i,1] & ((KHSRP[i,2] > ATAC_peak$start &KHSRP[i,2] < ATAC_peak$end) | (KHSRP[i,3] > ATAC_peak$start &KHSRP[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_KHSRP) != 0){
    ATAC_peak_KHSRP$V5 <-KHSRP[i,1]
    ATAC_peak_KHSRP$V6 <- KHSRP[i,2]
    ATAC_peak_KHSRP$V7 <- KHSRP[i,3]
    ATAC_peak_KHSRP$V8 <- KHSRP[i,6]
    ATAC_peak_KHSRP$V9 <- KHSRP[i,7]
    ATAC_peak_KHSRP_list = rbind(ATAC_peak_KHSRP_list,ATAC_peak_KHSRP)
  }
}
write.table(ATAC_peak_KHSRP_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KHSRP_list.txt",sep="\t")

ATAC_peak_KHSRP_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KHSRP_list.txt",sep="\t")

KHSRP_features=list("KHSRP" = ATAC_peak_KHSRP_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = KHSRP_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##SMAD5
SMAD5= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF973ZVM_SMAD5.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMAD5,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF973ZVM_SMAD5.txt",sep = "\t")

SMAD5= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF973ZVM_SMAD5.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMAD5)){
  x=paste(SMAD5[i,1],SMAD5[i,2],SMAD5[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMAD5_features = list(list$V1)

ATAC_peak_SMAD5_list = data.frame()
for(i in 1:nrow(SMAD5)){
  ATAC_peak_SMAD5<- ATAC_peak[ATAC_peak$x1 ==SMAD5[i,1] & ((SMAD5[i,2] > ATAC_peak$start &SMAD5[i,2] < ATAC_peak$end) | (SMAD5[i,3] > ATAC_peak$start &SMAD5[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMAD5) != 0){
    ATAC_peak_SMAD5$V5 <-SMAD5[i,1]
    ATAC_peak_SMAD5$V6 <- SMAD5[i,2]
    ATAC_peak_SMAD5$V7 <- SMAD5[i,3]
    ATAC_peak_SMAD5$V8 <- SMAD5[i,6]
    ATAC_peak_SMAD5$V9 <- SMAD5[i,7]
    ATAC_peak_SMAD5_list = rbind(ATAC_peak_SMAD5_list,ATAC_peak_SMAD5)
  }
}
write.table(ATAC_peak_SMAD5_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMAD5_list.txt",sep="\t")

ATAC_peak_SMAD5_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMAD5_list.txt",sep="\t")

SMAD5_features=list("SMAD5" = ATAC_peak_SMAD5_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMAD5_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##ZMIZ1
ZMIZ1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF976WVR_ZMIZ1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(ZMIZ1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF976WVR_ZMIZ1.txt",sep = "\t")

ZMIZ1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF976WVR_ZMIZ1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ZMIZ1)){
  x=paste(ZMIZ1[i,1],ZMIZ1[i,2],ZMIZ1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ZMIZ1_features = list(list$V1)

ATAC_peak_ZMIZ1_list = data.frame()
for(i in 1:nrow(ZMIZ1)){
  ATAC_peak_ZMIZ1<- ATAC_peak[ATAC_peak$x1 ==ZMIZ1[i,1] & ((ZMIZ1[i,2] > ATAC_peak$start &ZMIZ1[i,2] < ATAC_peak$end) | (ZMIZ1[i,3] > ATAC_peak$start &ZMIZ1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ZMIZ1) != 0){
    ATAC_peak_ZMIZ1$V5 <-ZMIZ1[i,1]
    ATAC_peak_ZMIZ1$V6 <- ZMIZ1[i,2]
    ATAC_peak_ZMIZ1$V7 <- ZMIZ1[i,3]
    ATAC_peak_ZMIZ1$V8 <- ZMIZ1[i,6]
    ATAC_peak_ZMIZ1$V9 <- ZMIZ1[i,7]
    ATAC_peak_ZMIZ1_list = rbind(ATAC_peak_ZMIZ1_list,ATAC_peak_ZMIZ1)
  }
}
write.table(ATAC_peak_ZMIZ1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZMIZ1_list.txt",sep="\t")

ATAC_peak_ZMIZ1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ZMIZ1_list.txt",sep="\t")

ZMIZ1_features=list("ZMIZ1" = ATAC_peak_ZMIZ1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ZMIZ1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##L3MBTL2
L3MBTL2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF979SYP_L3MBTL2.bed.gz",extraCols=extraCols_narrowPeak)
write.table(L3MBTL2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF979SYP_L3MBTL2.txt",sep = "\t")

L3MBTL2= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF979SYP_L3MBTL2.txt",sep = "\t")

list = vector()
for (i in  1:nrow(L3MBTL2)){
  x=paste(L3MBTL2[i,1],L3MBTL2[i,2],L3MBTL2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

L3MBTL2_features = list(list$V1)

ATAC_peak_L3MBTL2_list = data.frame()
for(i in 1:nrow(L3MBTL2)){
  ATAC_peak_L3MBTL2<- ATAC_peak[ATAC_peak$x1 ==L3MBTL2[i,1] & ((L3MBTL2[i,2] > ATAC_peak$start &L3MBTL2[i,2] < ATAC_peak$end) | (L3MBTL2[i,3] > ATAC_peak$start &L3MBTL2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_L3MBTL2) != 0){
    ATAC_peak_L3MBTL2$V5 <-L3MBTL2[i,1]
    ATAC_peak_L3MBTL2$V6 <- L3MBTL2[i,2]
    ATAC_peak_L3MBTL2$V7 <- L3MBTL2[i,3]
    ATAC_peak_L3MBTL2$V8 <- L3MBTL2[i,6]
    ATAC_peak_L3MBTL2$V9 <- L3MBTL2[i,7]
    ATAC_peak_L3MBTL2_list = rbind(ATAC_peak_L3MBTL2_list,ATAC_peak_L3MBTL2)
  }
}
write.table(ATAC_peak_L3MBTL2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_L3MBTL2_list.txt",sep="\t")

ATAC_peak_L3MBTL2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_L3MBTL2_list.txt",sep="\t")

L3MBTL2_features=list("L3MBTL2" = ATAC_peak_L3MBTL2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = L3MBTL2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##SMAD1
SMAD1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF980CWC_SMAD1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SMAD1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF980CWC_SMAD1.txt",sep = "\t")

SMAD1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF980CWC_SMAD1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SMAD1)){
  x=paste(SMAD1[i,1],SMAD1[i,2],SMAD1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SMAD1_features = list(list$V1)

ATAC_peak_SMAD1_list = data.frame()
for(i in 1:nrow(SMAD1)){
  ATAC_peak_SMAD1<- ATAC_peak[ATAC_peak$x1 ==SMAD1[i,1] & ((SMAD1[i,2] > ATAC_peak$start &SMAD1[i,2] < ATAC_peak$end) | (SMAD1[i,3] > ATAC_peak$start &SMAD1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SMAD1) != 0){
    ATAC_peak_SMAD1$V5 <-SMAD1[i,1]
    ATAC_peak_SMAD1$V6 <- SMAD1[i,2]
    ATAC_peak_SMAD1$V7 <- SMAD1[i,3]
    ATAC_peak_SMAD1$V8 <- SMAD1[i,6]
    ATAC_peak_SMAD1$V9 <- SMAD1[i,7]
    ATAC_peak_SMAD1_list = rbind(ATAC_peak_SMAD1_list,ATAC_peak_SMAD1)
  }
}
write.table(ATAC_peak_SMAD1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMAD1_list.txt",sep="\t")

ATAC_peak_SMAD1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SMAD1_list.txt",sep="\t")

SMAD1_features=list("SMAD1" = ATAC_peak_SMAD1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SMAD1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##DNMT1
DNMT1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF981HQX_DNMT1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(DNMT1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF981HQX_DNMT1.txt",sep = "\t")

DNMT1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF981HQX_DNMT1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(DNMT1)){
  x=paste(DNMT1[i,1],DNMT1[i,2],DNMT1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

DNMT1_features = list(list$V1)

ATAC_peak_DNMT1_list = data.frame()
for(i in 1:nrow(DNMT1)){
  ATAC_peak_DNMT1<- ATAC_peak[ATAC_peak$x1 ==DNMT1[i,1] & ((DNMT1[i,2] > ATAC_peak$start &DNMT1[i,2] < ATAC_peak$end) | (DNMT1[i,3] > ATAC_peak$start &DNMT1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_DNMT1) != 0){
    ATAC_peak_DNMT1$V5 <-DNMT1[i,1]
    ATAC_peak_DNMT1$V6 <- DNMT1[i,2]
    ATAC_peak_DNMT1$V7 <- DNMT1[i,3]
    ATAC_peak_DNMT1$V8 <- DNMT1[i,6]
    ATAC_peak_DNMT1$V9 <- DNMT1[i,7]
    ATAC_peak_DNMT1_list = rbind(ATAC_peak_DNMT1_list,ATAC_peak_DNMT1)
  }
}
write.table(ATAC_peak_DNMT1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DNMT1_list.txt",sep="\t")

ATAC_peak_DNMT1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DNMT1_list.txt",sep="\t")

DNMT1_features=list("DNMT1" = ATAC_peak_DNMT1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = DNMT1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SP1
SP1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF984FGT_SP1_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF984FGT_SP1_age.txt",sep = "\t")

SP1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF984FGT_SP1_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SP1)){
  x=paste(SP1[i,1],SP1[i,2],SP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SP1_features = list(list$V1)

ATAC_peak_SP1_list = data.frame()
for(i in 1:nrow(SP1)){
  ATAC_peak_SP1<- ATAC_peak[ATAC_peak$x1 ==SP1[i,1] & ((SP1[i,2] > ATAC_peak$start &SP1[i,2] < ATAC_peak$end) | (SP1[i,3] > ATAC_peak$start &SP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SP1) != 0){
    ATAC_peak_SP1$V5 <-SP1[i,1]
    ATAC_peak_SP1$V6 <- SP1[i,2]
    ATAC_peak_SP1$V7 <- SP1[i,3]
    ATAC_peak_SP1$V8 <- SP1[i,6]
    ATAC_peak_SP1$V9 <- SP1[i,7]
    ATAC_peak_SP1_list = rbind(ATAC_peak_SP1_list,ATAC_peak_SP1)
  }
}
write.table(ATAC_peak_SP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SP1_list.txt",sep="\t")

ATAC_peak_SP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SP1_list.txt",sep="\t")

SP1_features=list("SP1" = ATAC_peak_SP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##DDIT3
DDIT3= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF989RVM_DDIT3_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(DDIT3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF989RVM_DDIT3.txt",sep = "\t")

DDIT3= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF989RVM_DDIT3.txt",sep = "\t")

list = vector()
for (i in  1:nrow(DDIT3)){
  x=paste(DDIT3[i,1],DDIT3[i,2],DDIT3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

DDIT3_features = list(list$V1)

ATAC_peak_DDIT3_list = data.frame()
for(i in 1:nrow(DDIT3)){
  ATAC_peak_DDIT3<- ATAC_peak[ATAC_peak$x1 ==DDIT3[i,1] & ((DDIT3[i,2] > ATAC_peak$start &DDIT3[i,2] < ATAC_peak$end) | (DDIT3[i,3] > ATAC_peak$start &DDIT3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_DDIT3) != 0){
    ATAC_peak_DDIT3$V5 <-DDIT3[i,1]
    ATAC_peak_DDIT3$V6 <- DDIT3[i,2]
    ATAC_peak_DDIT3$V7 <- DDIT3[i,3]
    ATAC_peak_DDIT3$V8 <- DDIT3[i,6]
    ATAC_peak_DDIT3$V9 <- DDIT3[i,7]
    ATAC_peak_DDIT3_list = rbind(ATAC_peak_DDIT3_list,ATAC_peak_DDIT3)
  }
}
write.table(ATAC_peak_DDIT3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DDIT3_list.txt",sep="\t")

ATAC_peak_DDIT3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_DDIT3_list.txt",sep="\t")

DDIT3_features=list("DDIT3" = ATAC_peak_DDIT3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = DDIT3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##PML
PML= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF989RWS_PML_age.bed.gz",extraCols=extraCols_narrowPeak)
write.table(PML,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF989RWS_PML_age.txt",sep = "\t")

PML= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF989RWS_PML_age.txt",sep = "\t")

list = vector()
for (i in  1:nrow(PML)){
  x=paste(PML[i,1],PML[i,2],PML[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

PML_features = list(list$V1)

ATAC_peak_PML_list = data.frame()
for(i in 1:nrow(PML)){
  ATAC_peak_PML<- ATAC_peak[ATAC_peak$x1 ==PML[i,1] & ((PML[i,2] > ATAC_peak$start &PML[i,2] < ATAC_peak$end) | (PML[i,3] > ATAC_peak$start &PML[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_PML) != 0){
    ATAC_peak_PML$V5 <-PML[i,1]
    ATAC_peak_PML$V6 <- PML[i,2]
    ATAC_peak_PML$V7 <- PML[i,3]
    ATAC_peak_PML$V8 <- PML[i,6]
    ATAC_peak_PML$V9 <- PML[i,7]
    ATAC_peak_PML_list = rbind(ATAC_peak_PML_list,ATAC_peak_PML)
  }
}
write.table(ATAC_peak_PML_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PML_list.txt",sep="\t")

ATAC_peak_PML_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_PML_list.txt",sep="\t")

PML_features=list("PML" = ATAC_peak_PML_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = PML_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##RBM34
RBM34= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF990RDF_RBM34.bed.gz",extraCols=extraCols_narrowPeak)
write.table(RBM34,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF990RDF_RBM34.txt",sep = "\t")

RBM34= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF990RDF_RBM34.txt",sep = "\t")

list = vector()
for (i in  1:nrow(RBM34)){
  x=paste(RBM34[i,1],RBM34[i,2],RBM34[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

RBM34_features = list(list$V1)

ATAC_peak_RBM34_list = data.frame()
for(i in 1:nrow(RBM34)){
  ATAC_peak_RBM34<- ATAC_peak[ATAC_peak$x1 ==RBM34[i,1] & ((RBM34[i,2] > ATAC_peak$start &RBM34[i,2] < ATAC_peak$end) | (RBM34[i,3] > ATAC_peak$start &RBM34[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_RBM34) != 0){
    ATAC_peak_RBM34$V5 <-RBM34[i,1]
    ATAC_peak_RBM34$V6 <- RBM34[i,2]
    ATAC_peak_RBM34$V7 <- RBM34[i,3]
    ATAC_peak_RBM34$V8 <- RBM34[i,6]
    ATAC_peak_RBM34$V9 <- RBM34[i,7]
    ATAC_peak_RBM34_list = rbind(ATAC_peak_RBM34_list,ATAC_peak_RBM34)
  }
}
write.table(ATAC_peak_RBM34_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RBM34_list.txt",sep="\t")

ATAC_peak_RBM34_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_RBM34_list.txt",sep="\t")

RBM34_features=list("RBM34" = ATAC_peak_RBM34_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = RBM34_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HES1
HES1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF991IRG_HES1.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HES1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF991IRG_HES1.txt",sep = "\t")

HES1= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF991IRG_HES1.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HES1)){
  x=paste(HES1[i,1],HES1[i,2],HES1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HES1_features = list(list$V1)

ATAC_peak_HES1_list = data.frame()
for(i in 1:nrow(HES1)){
  ATAC_peak_HES1<- ATAC_peak[ATAC_peak$x1 ==HES1[i,1] & ((HES1[i,2] > ATAC_peak$start &HES1[i,2] < ATAC_peak$end) | (HES1[i,3] > ATAC_peak$start &HES1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HES1) != 0){
    ATAC_peak_HES1$V5 <-HES1[i,1]
    ATAC_peak_HES1$V6 <- HES1[i,2]
    ATAC_peak_HES1$V7 <- HES1[i,3]
    ATAC_peak_HES1$V8 <- HES1[i,6]
    ATAC_peak_HES1$V9 <- HES1[i,7]
    ATAC_peak_HES1_list = rbind(ATAC_peak_HES1_list,ATAC_peak_HES1)
  }
}
write.table(ATAC_peak_HES1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HES1_list.txt",sep="\t")

ATAC_peak_HES1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HES1_list.txt",sep="\t")

HES1_features=list("HES1" = ATAC_peak_HES1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HES1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HLTF
HLTF= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF992VUZ_HLTF.bed.gz",extraCols=extraCols_narrowPeak)
write.table(HLTF,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF992VUZ_HLTF.txt",sep = "\t")

HLTF= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF992VUZ_HLTF.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HLTF)){
  x=paste(HLTF[i,1],HLTF[i,2],HLTF[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HLTF_features = list(list$V1)

ATAC_peak_HLTF_list = data.frame()
for(i in 1:nrow(HLTF)){
  ATAC_peak_HLTF<- ATAC_peak[ATAC_peak$x1 ==HLTF[i,1] & ((HLTF[i,2] > ATAC_peak$start &HLTF[i,2] < ATAC_peak$end) | (HLTF[i,3] > ATAC_peak$start &HLTF[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HLTF) != 0){
    ATAC_peak_HLTF$V5 <-HLTF[i,1]
    ATAC_peak_HLTF$V6 <- HLTF[i,2]
    ATAC_peak_HLTF$V7 <- HLTF[i,3]
    ATAC_peak_HLTF$V8 <- HLTF[i,6]
    ATAC_peak_HLTF$V9 <- HLTF[i,7]
    ATAC_peak_HLTF_list = rbind(ATAC_peak_HLTF_list,ATAC_peak_HLTF)
  }
}
write.table(ATAC_peak_HLTF_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HLTF_list.txt",sep="\t")

ATAC_peak_HLTF_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HLTF_list.txt",sep="\t")

HLTF_features=list("HLTF" = ATAC_peak_HLTF_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HLTF_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##SFPQ
SFPQ= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF995MDH_SFPQ.bed.gz",extraCols=extraCols_narrowPeak)
write.table(SFPQ,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF995MDH_SFPQ.txt",sep = "\t")

SFPQ= read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ENCFF995MDH_SFPQ.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SFPQ)){
  x=paste(SFPQ[i,1],SFPQ[i,2],SFPQ[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SFPQ_features = list(list$V1)

ATAC_peak_SFPQ_list = data.frame()
for(i in 1:nrow(SFPQ)){
  ATAC_peak_SFPQ<- ATAC_peak[ATAC_peak$x1 ==SFPQ[i,1] & ((SFPQ[i,2] > ATAC_peak$start &SFPQ[i,2] < ATAC_peak$end) | (SFPQ[i,3] > ATAC_peak$start &SFPQ[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SFPQ) != 0){
    ATAC_peak_SFPQ$V5 <-SFPQ[i,1]
    ATAC_peak_SFPQ$V6 <- SFPQ[i,2]
    ATAC_peak_SFPQ$V7 <- SFPQ[i,3]
    ATAC_peak_SFPQ$V8 <- SFPQ[i,6]
    ATAC_peak_SFPQ$V9 <- SFPQ[i,7]
    ATAC_peak_SFPQ_list = rbind(ATAC_peak_SFPQ_list,ATAC_peak_SFPQ)
  }
}
write.table(ATAC_peak_SFPQ_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SFPQ_list.txt",sep="\t")

ATAC_peak_SFPQ_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SFPQ_list.txt",sep="\t")

SFPQ_features=list("SFPQ" = ATAC_peak_SFPQ_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SFPQ_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


#########
##KLF4

KLF4= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.KLF4.AllCell.bed")
write.table(KLF4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.KLF4.AllCell.txt",sep = "\t")

KLF4 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.KLF4.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(KLF4)){
  x=paste(KLF4[i,1],KLF4[i,2],KLF4[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

KLF4_features = list(list$V1)

ATAC_peak_KLF4_list = data.frame()
for(i in 1:nrow(KLF4)){
  ATAC_peak_KLF4<- ATAC_peak[ATAC_peak$x1 ==KLF4[i,1] & ((KLF4[i,2] > ATAC_peak$start &KLF4[i,2] < ATAC_peak$end) | (KLF4[i,3] > ATAC_peak$start &KLF4[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_KLF4) != 0){
    ATAC_peak_KLF4$V5 <-KLF4[i,1]
    ATAC_peak_KLF4$V6 <- KLF4[i,2]
    ATAC_peak_KLF4$V7 <- KLF4[i,3]
    ATAC_peak_KLF4$V8 <- KLF4[i,6]
    ATAC_peak_KLF4$V9 <- KLF4[i,7]
    ATAC_peak_KLF4_list = rbind(ATAC_peak_KLF4_list,ATAC_peak_KLF4)
  }
}
write.table(ATAC_peak_KLF4_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KLF4_list.txt",sep="\t")

ATAC_peak_KLF4_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_KLF4_list.txt",sep="\t")

KLF4_features=list("KLF4" = ATAC_peak_KLF4_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = KLF4_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##LMNA

LMNA= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.LMNA_age.AllCell.bed")
write.table(LMNA,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.LMNA_age.AllCell.txt",sep = "\t")

LMNA = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.LMNA_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(LMNA)){
  x=paste(LMNA[i,1],LMNA[i,2],LMNA[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

LMNA_features = list(list$V1)

ATAC_peak_LMNA_list = data.frame()
for(i in 1:nrow(LMNA)){
  ATAC_peak_LMNA<- ATAC_peak[ATAC_peak$x1 ==LMNA[i,1] & ((LMNA[i,2] > ATAC_peak$start &LMNA[i,2] < ATAC_peak$end) | (LMNA[i,3] > ATAC_peak$start &LMNA[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_LMNA) != 0){
    ATAC_peak_LMNA$V5 <-LMNA[i,1]
    ATAC_peak_LMNA$V6 <- LMNA[i,2]
    ATAC_peak_LMNA$V7 <- LMNA[i,3]
    ATAC_peak_LMNA$V8 <- LMNA[i,6]
    ATAC_peak_LMNA$V9 <- LMNA[i,7]
    ATAC_peak_LMNA_list = rbind(ATAC_peak_LMNA_list,ATAC_peak_LMNA)
  }
}
write.table(ATAC_peak_LMNA_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_LMNA_list.txt",sep="\t")

ATAC_peak_LMNA_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_LMNA_list.txt",sep="\t")

LMNA_features=list("LMNA" = ATAC_peak_LMNA_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = LMNA_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NFKB2

NFKB2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.NFKB2.AllCell.bed")
write.table(NFKB2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.NFKB2.AllCell.txt",sep = "\t")

NFKB2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.NFKB2.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NFKB2)){
  x=paste(NFKB2[i,1],NFKB2[i,2],NFKB2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NFKB2_features = list(list$V1)

ATAC_peak_NFKB2_list = data.frame()
for(i in 1:nrow(NFKB2)){
  ATAC_peak_NFKB2<- ATAC_peak[ATAC_peak$x1 ==NFKB2[i,1] & ((NFKB2[i,2] > ATAC_peak$start &NFKB2[i,2] < ATAC_peak$end) | (NFKB2[i,3] > ATAC_peak$start &NFKB2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NFKB2) != 0){
    ATAC_peak_NFKB2$V5 <-NFKB2[i,1]
    ATAC_peak_NFKB2$V6 <- NFKB2[i,2]
    ATAC_peak_NFKB2$V7 <- NFKB2[i,3]
    ATAC_peak_NFKB2$V8 <- NFKB2[i,6]
    ATAC_peak_NFKB2$V9 <- NFKB2[i,7]
    ATAC_peak_NFKB2_list = rbind(ATAC_peak_NFKB2_list,ATAC_peak_NFKB2)
  }
}
write.table(ATAC_peak_NFKB2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFKB2_list.txt",sep="\t")

ATAC_peak_NFKB2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFKB2_list.txt",sep="\t")

NFKB2_features=list("NFKB2" = ATAC_peak_NFKB2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NFKB2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##STAT3

STAT3= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.STAT3.AllCell.bed")
write.table(STAT3,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.STAT3.AllCell.txt",sep = "\t")

STAT3 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.STAT3.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(STAT3)){
  x=paste(STAT3[i,1],STAT3[i,2],STAT3[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

STAT3_features = list(list$V1)

ATAC_peak_STAT3_list = data.frame()
for(i in 1:nrow(STAT3)){
  ATAC_peak_STAT3<- ATAC_peak[ATAC_peak$x1 ==STAT3[i,1] & ((STAT3[i,2] > ATAC_peak$start &STAT3[i,2] < ATAC_peak$end) | (STAT3[i,3] > ATAC_peak$start &STAT3[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_STAT3) != 0){
    ATAC_peak_STAT3$V5 <-STAT3[i,1]
    ATAC_peak_STAT3$V6 <- STAT3[i,2]
    ATAC_peak_STAT3$V7 <- STAT3[i,3]
    ATAC_peak_STAT3$V8 <- STAT3[i,6]
    ATAC_peak_STAT3$V9 <- STAT3[i,7]
    ATAC_peak_STAT3_list = rbind(ATAC_peak_STAT3_list,ATAC_peak_STAT3)
  }
}
write.table(ATAC_peak_STAT3_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_STAT3_list.txt",sep="\t")

ATAC_peak_STAT3_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_STAT3_list.txt",sep="\t")

STAT3_features=list("STAT3" = ATAC_peak_STAT3_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = STAT3_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TP63

TP63= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.TP63.AllCell.bed")
write.table(TP63,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.TP63.AllCell.txt",sep = "\t")

TP63 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.20.TP63.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TP63)){
  x=paste(TP63[i,1],TP63[i,2],TP63[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TP63_features = list(list$V1)

ATAC_peak_TP63_list = data.frame()
for(i in 1:nrow(TP63)){
  ATAC_peak_TP63<- ATAC_peak[ATAC_peak$x1 ==TP63[i,1] & ((TP63[i,2] > ATAC_peak$start &TP63[i,2] < ATAC_peak$end) | (TP63[i,3] > ATAC_peak$start &TP63[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TP63) != 0){
    ATAC_peak_TP63$V5 <-TP63[i,1]
    ATAC_peak_TP63$V6 <- TP63[i,2]
    ATAC_peak_TP63$V7 <- TP63[i,3]
    ATAC_peak_TP63$V8 <- TP63[i,6]
    ATAC_peak_TP63$V9 <- TP63[i,7]
    ATAC_peak_TP63_list = rbind(ATAC_peak_TP63_list,ATAC_peak_TP63)
  }
}
write.table(ATAC_peak_TP63_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TP63_list.txt",sep="\t")

ATAC_peak_TP63_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TP63_list.txt",sep="\t")

TP63_features=list("TP63" = ATAC_peak_TP63_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TP63_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CDKN1B

CDKN1B= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.50.CDKN1B.AllCell.bed")
write.table(CDKN1B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.50.CDKN1B.AllCell.txt",sep = "\t")

CDKN1B = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.50.CDKN1B.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CDKN1B)){
  x=paste(CDKN1B[i,1],CDKN1B[i,2],CDKN1B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CDKN1B_features = list(list$V1)

ATAC_peak_CDKN1B_list = data.frame()
for(i in 1:nrow(CDKN1B)){
  ATAC_peak_CDKN1B<- ATAC_peak[ATAC_peak$x1 ==CDKN1B[i,1] & ((CDKN1B[i,2] > ATAC_peak$start &CDKN1B[i,2] < ATAC_peak$end) | (CDKN1B[i,3] > ATAC_peak$start &CDKN1B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CDKN1B) != 0){
    ATAC_peak_CDKN1B$V5 <-CDKN1B[i,1]
    ATAC_peak_CDKN1B$V6 <- CDKN1B[i,2]
    ATAC_peak_CDKN1B$V7 <- CDKN1B[i,3]
    ATAC_peak_CDKN1B$V8 <- CDKN1B[i,6]
    ATAC_peak_CDKN1B$V9 <- CDKN1B[i,7]
    ATAC_peak_CDKN1B_list = rbind(ATAC_peak_CDKN1B_list,ATAC_peak_CDKN1B)
  }
}
write.table(ATAC_peak_CDKN1B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CDKN1B_list.txt",sep="\t")

ATAC_peak_CDKN1B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CDKN1B_list.txt",sep="\t")

CDKN1B_features=list("CDKN1B" = ATAC_peak_CDKN1B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CDKN1B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##AR

AR= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.05.AR.AllCell.bed")
write.table(AR,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.AR_age.AllCell.txt",sep = "\t")

AR = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.AR_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(AR)){
  x=paste(AR[i,1],AR[i,2],AR[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

AR_features = list(list$V1)

ATAC_peak_AR_list = data.frame()
for(i in 1:nrow(AR)){
  ATAC_peak_AR<- ATAC_peak[ATAC_peak$x1 ==AR[i,1] & ((AR[i,2] > ATAC_peak$start &AR[i,2] < ATAC_peak$end) | (AR[i,3] > ATAC_peak$start &AR[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_AR) != 0){
    ATAC_peak_AR$V5 <-AR[i,1]
    ATAC_peak_AR$V6 <- AR[i,2]
    ATAC_peak_AR$V7 <- AR[i,3]
    ATAC_peak_AR$V8 <- AR[i,6]
    ATAC_peak_AR$V9 <- AR[i,7]
    ATAC_peak_AR_list = rbind(ATAC_peak_AR_list,ATAC_peak_AR)
  }
}
write.table(ATAC_peak_AR_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_AR_list.txt",sep="\t")

ATAC_peak_AR_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_AR_list.txt",sep="\t")

AR_features=list("AR" = ATAC_peak_AR_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = AR_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##ARNTL

ARNTL= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.ARNTL_age.AllCell.bed")
write.table(ARNTL,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.ARNTL_age.AllCell.txt",sep = "\t")

ARNTL = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.ARNTL_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(ARNTL)){
  x=paste(ARNTL[i,1],ARNTL[i,2],ARNTL[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

ARNTL_features = list(list$V1)

ATAC_peak_ARNTL_list = data.frame()
for(i in 1:nrow(ARNTL)){
  ATAC_peak_ARNTL<- ATAC_peak[ATAC_peak$x1 ==ARNTL[i,1] & ((ARNTL[i,2] > ATAC_peak$start &ARNTL[i,2] < ATAC_peak$end) | (ARNTL[i,3] > ATAC_peak$start &ARNTL[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_ARNTL) != 0){
    ATAC_peak_ARNTL$V5 <-ARNTL[i,1]
    ATAC_peak_ARNTL$V6 <- ARNTL[i,2]
    ATAC_peak_ARNTL$V7 <- ARNTL[i,3]
    ATAC_peak_ARNTL$V8 <- ARNTL[i,6]
    ATAC_peak_ARNTL$V9 <- ARNTL[i,7]
    ATAC_peak_ARNTL_list = rbind(ATAC_peak_ARNTL_list,ATAC_peak_ARNTL)
  }
}
write.table(ATAC_peak_ARNTL_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARNTL_list.txt",sep="\t")

ATAC_peak_ARNTL_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_ARNTL_list.txt",sep="\t")

ARNTL_features=list("ARNTL" = ATAC_peak_ARNTL_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = ARNTL_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CDK7

CDK7= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.CDK7_age.AllCell.bed")
write.table(CDK7,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.CDK7_age.AllCell.txt",sep = "\t")

CDK7 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.CDK7_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CDK7)){
  x=paste(CDK7[i,1],CDK7[i,2],CDK7[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CDK7_features = list(list$V1)

ATAC_peak_CDK7_list = data.frame()
for(i in 1:nrow(CDK7)){
  ATAC_peak_CDK7<- ATAC_peak[ATAC_peak$x1 ==CDK7[i,1] & ((CDK7[i,2] > ATAC_peak$start &CDK7[i,2] < ATAC_peak$end) | (CDK7[i,3] > ATAC_peak$start &CDK7[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CDK7) != 0){
    ATAC_peak_CDK7$V5 <-CDK7[i,1]
    ATAC_peak_CDK7$V6 <- CDK7[i,2]
    ATAC_peak_CDK7$V7 <- CDK7[i,3]
    ATAC_peak_CDK7$V8 <- CDK7[i,6]
    ATAC_peak_CDK7$V9 <- CDK7[i,7]
    ATAC_peak_CDK7_list = rbind(ATAC_peak_CDK7_list,ATAC_peak_CDK7)
  }
}
write.table(ATAC_peak_CDK7_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CDK7_list.txt",sep="\t")

ATAC_peak_CDK7_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CDK7_list.txt",sep="\t")

CDK7_features=list("CDK7" = ATAC_peak_CDK7_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CDK7_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##CTNNB1

CTNNB1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.CTNNB1_age.AllCell.bed")
write.table(CTNNB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.CTNNB1_age.AllCell.txt",sep = "\t")

CTNNB1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.CTNNB1_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(CTNNB1)){
  x=paste(CTNNB1[i,1],CTNNB1[i,2],CTNNB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

CTNNB1_features = list(list$V1)

ATAC_peak_CTNNB1_list = data.frame()
for(i in 1:nrow(CTNNB1)){
  ATAC_peak_CTNNB1<- ATAC_peak[ATAC_peak$x1 ==CTNNB1[i,1] & ((CTNNB1[i,2] > ATAC_peak$start &CTNNB1[i,2] < ATAC_peak$end) | (CTNNB1[i,3] > ATAC_peak$start &CTNNB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_CTNNB1) != 0){
    ATAC_peak_CTNNB1$V5 <-CTNNB1[i,1]
    ATAC_peak_CTNNB1$V6 <- CTNNB1[i,2]
    ATAC_peak_CTNNB1$V7 <- CTNNB1[i,3]
    ATAC_peak_CTNNB1$V8 <- CTNNB1[i,6]
    ATAC_peak_CTNNB1$V9 <- CTNNB1[i,7]
    ATAC_peak_CTNNB1_list = rbind(ATAC_peak_CTNNB1_list,ATAC_peak_CTNNB1)
  }
}
write.table(ATAC_peak_CTNNB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CTNNB1_list.txt",sep="\t")

ATAC_peak_CTNNB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_CTNNB1_list.txt",sep="\t")

CTNNB1_features=list("CTNNB1" = ATAC_peak_CTNNB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = CTNNB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##FOXO1

FOXO1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.FOXO1_age.AllCell.bed")
write.table(FOXO1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.FOXO1_age.AllCell.txt",sep = "\t")

FOXO1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.FOXO1_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(FOXO1)){
  x=paste(FOXO1[i,1],FOXO1[i,2],FOXO1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

FOXO1_features = list(list$V1)

ATAC_peak_FOXO1_list = data.frame()
for(i in 1:nrow(FOXO1)){
  ATAC_peak_FOXO1<- ATAC_peak[ATAC_peak$x1 ==FOXO1[i,1] & ((FOXO1[i,2] > ATAC_peak$start &FOXO1[i,2] < ATAC_peak$end) | (FOXO1[i,3] > ATAC_peak$start &FOXO1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_FOXO1) != 0){
    ATAC_peak_FOXO1$V5 <-FOXO1[i,1]
    ATAC_peak_FOXO1$V6 <- FOXO1[i,2]
    ATAC_peak_FOXO1$V7 <- FOXO1[i,3]
    ATAC_peak_FOXO1$V8 <- FOXO1[i,6]
    ATAC_peak_FOXO1$V9 <- FOXO1[i,7]
    ATAC_peak_FOXO1_list = rbind(ATAC_peak_FOXO1_list,ATAC_peak_FOXO1)
  }
}
write.table(ATAC_peak_FOXO1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXO1_list.txt",sep="\t")

ATAC_peak_FOXO1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_FOXO1_list.txt",sep="\t")

FOXO1_features=list("FOXO1" = ATAC_peak_FOXO1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = FOXO1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##HIC1

HIC1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.HIC1_age.AllCell.bed")
write.table(HIC1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.HIC1_age.AllCell.txt",sep = "\t")

HIC1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.HIC1_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(HIC1)){
  x=paste(HIC1[i,1],HIC1[i,2],HIC1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HIC1_features = list(list$V1)

ATAC_peak_HIC1_list = data.frame()
for(i in 1:nrow(HIC1)){
  ATAC_peak_HIC1<- ATAC_peak[ATAC_peak$x1 ==HIC1[i,1] & ((HIC1[i,2] > ATAC_peak$start &HIC1[i,2] < ATAC_peak$end) | (HIC1[i,3] > ATAC_peak$start &HIC1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HIC1) != 0){
    ATAC_peak_HIC1$V5 <-HIC1[i,1]
    ATAC_peak_HIC1$V6 <- HIC1[i,2]
    ATAC_peak_HIC1$V7 <- HIC1[i,3]
    ATAC_peak_HIC1$V8 <- HIC1[i,6]
    ATAC_peak_HIC1$V9 <- HIC1[i,7]
    ATAC_peak_HIC1_list = rbind(ATAC_peak_HIC1_list,ATAC_peak_HIC1)
  }
}
write.table(ATAC_peak_HIC1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HIC1_list.txt",sep="\t")

ATAC_peak_HIC1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HIC1_list.txt",sep="\t")

HIC1_features=list("HIC1" = ATAC_peak_HIC1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HIC1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##LMNB1

LMNB1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.LMNB1_age.AllCell.bed")
write.table(LMNB1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.LMNB1_age.AllCell.txt",sep = "\t")

LMNB1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.LMNB1_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(LMNB1)){
  x=paste(LMNB1[i,1],LMNB1[i,2],LMNB1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

LMNB1_features = list(list$V1)

ATAC_peak_LMNB1_list = data.frame()
for(i in 1:nrow(LMNB1)){
  ATAC_peak_LMNB1<- ATAC_peak[ATAC_peak$x1 ==LMNB1[i,1] & ((LMNB1[i,2] > ATAC_peak$start &LMNB1[i,2] < ATAC_peak$end) | (LMNB1[i,3] > ATAC_peak$start &LMNB1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_LMNB1) != 0){
    ATAC_peak_LMNB1$V5 <-LMNB1[i,1]
    ATAC_peak_LMNB1$V6 <- LMNB1[i,2]
    ATAC_peak_LMNB1$V7 <- LMNB1[i,3]
    ATAC_peak_LMNB1$V8 <- LMNB1[i,6]
    ATAC_peak_LMNB1$V9 <- LMNB1[i,7]
    ATAC_peak_LMNB1_list = rbind(ATAC_peak_LMNB1_list,ATAC_peak_LMNB1)
  }
}
write.table(ATAC_peak_LMNB1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_LMNB1_list.txt",sep="\t")

ATAC_peak_LMNB1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_LMNB1_list.txt",sep="\t")

LMNB1_features=list("LMNB1" = ATAC_peak_LMNB1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = LMNB1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##MAPK14

MAPK14= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.MAPK14_age.AllCell.bed")
write.table(MAPK14,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.MAPK14_age.AllCell.txt",sep = "\t")

MAPK14 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.MAPK14_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(MAPK14)){
  x=paste(MAPK14[i,1],MAPK14[i,2],MAPK14[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

MAPK14_features = list(list$V1)

ATAC_peak_MAPK14_list = data.frame()
for(i in 1:nrow(MAPK14)){
  ATAC_peak_MAPK14<- ATAC_peak[ATAC_peak$x1 ==MAPK14[i,1] & ((MAPK14[i,2] > ATAC_peak$start &MAPK14[i,2] < ATAC_peak$end) | (MAPK14[i,3] > ATAC_peak$start &MAPK14[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_MAPK14) != 0){
    ATAC_peak_MAPK14$V5 <-MAPK14[i,1]
    ATAC_peak_MAPK14$V6 <- MAPK14[i,2]
    ATAC_peak_MAPK14$V7 <- MAPK14[i,3]
    ATAC_peak_MAPK14$V8 <- MAPK14[i,6]
    ATAC_peak_MAPK14$V9 <- MAPK14[i,7]
    ATAC_peak_MAPK14_list = rbind(ATAC_peak_MAPK14_list,ATAC_peak_MAPK14)
  }
}
write.table(ATAC_peak_MAPK14_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MAPK14_list.txt",sep="\t")

ATAC_peak_MAPK14_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_MAPK14_list.txt",sep="\t")

MAPK14_features=list("MAPK14" = ATAC_peak_MAPK14_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = MAPK14_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NCOR2

NCOR2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.NCOR2_age.AllCell.bed")
write.table(NCOR2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.NCOR2_age.AllCell.txt",sep = "\t")

NCOR2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.NCOR2_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NCOR2)){
  x=paste(NCOR2[i,1],NCOR2[i,2],NCOR2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NCOR2_features = list(list$V1)

ATAC_peak_NCOR2_list = data.frame()
for(i in 1:nrow(NCOR2)){
  ATAC_peak_NCOR2<- ATAC_peak[ATAC_peak$x1 ==NCOR2[i,1] & ((NCOR2[i,2] > ATAC_peak$start &NCOR2[i,2] < ATAC_peak$end) | (NCOR2[i,3] > ATAC_peak$start &NCOR2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NCOR2) != 0){
    ATAC_peak_NCOR2$V5 <-NCOR2[i,1]
    ATAC_peak_NCOR2$V6 <- NCOR2[i,2]
    ATAC_peak_NCOR2$V7 <- NCOR2[i,3]
    ATAC_peak_NCOR2$V8 <- NCOR2[i,6]
    ATAC_peak_NCOR2$V9 <- NCOR2[i,7]
    ATAC_peak_NCOR2_list = rbind(ATAC_peak_NCOR2_list,ATAC_peak_NCOR2)
  }
}
write.table(ATAC_peak_NCOR2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOR2_list.txt",sep="\t")

ATAC_peak_NCOR2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NCOR2_list.txt",sep="\t")

NCOR2_features=list("NCOR2" = ATAC_peak_NCOR2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NCOR2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##NFE2L2

NFE2L2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.NFE2L2_age.AllCell.bed")
write.table(NFE2L2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.NFE2L2_age.AllCell.txt",sep = "\t")

NFE2L2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.NFE2L2_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(NFE2L2)){
  x=paste(NFE2L2[i,1],NFE2L2[i,2],NFE2L2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

NFE2L2_features = list(list$V1)

ATAC_peak_NFE2L2_list = data.frame()
for(i in 1:nrow(NFE2L2)){
  ATAC_peak_NFE2L2<- ATAC_peak[ATAC_peak$x1 ==NFE2L2[i,1] & ((NFE2L2[i,2] > ATAC_peak$start&NFE2L2[i,2] < ATAC_peak$end) | (NFE2L2[i,3] > ATAC_peak$start &NFE2L2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_NFE2L2) != 0){
    ATAC_peak_NFE2L2$V5 <-NFE2L2[i,1]
    ATAC_peak_NFE2L2$V6 <- NFE2L2[i,2]
    ATAC_peak_NFE2L2$V7 <- NFE2L2[i,3]
    ATAC_peak_NFE2L2$V8 <- NFE2L2[i,6]
    ATAC_peak_NFE2L2$V9 <- NFE2L2[i,7]
    ATAC_peak_NFE2L2_list = rbind(ATAC_peak_NFE2L2_list,ATAC_peak_NFE2L2)
  }
}
write.table(ATAC_peak_NFE2L2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFE2L2_list.txt",sep="\t")

ATAC_peak_NFE2L2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_NFE2L2_list.txt",sep="\t")

NFE2L2_features=list("NFE2L2" = ATAC_peak_NFE2L2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = NFE2L2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##POLA1

POLA1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.POLA1_age.AllCell.bed")
write.table(POLA1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.POLA1_age.AllCell.txt",sep = "\t")

POLA1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.POLA1_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(POLA1)){
  x=paste(POLA1[i,1],POLA1[i,2],POLA1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

POLA1_features = list(list$V1)

ATAC_peak_POLA1_list = data.frame()
for(i in 1:nrow(POLA1)){
  ATAC_peak_POLA1<- ATAC_peak[ATAC_peak$x1 ==POLA1[i,1] & ((POLA1[i,2] > ATAC_peak$start &POLA1[i,2] < ATAC_peak$end) | (POLA1[i,3] > ATAC_peak$start &POLA1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_POLA1) != 0){
    ATAC_peak_POLA1$V5 <-POLA1[i,1]
    ATAC_peak_POLA1$V6 <- POLA1[i,2]
    ATAC_peak_POLA1$V7 <- POLA1[i,3]
    ATAC_peak_POLA1$V8 <- POLA1[i,6]
    ATAC_peak_POLA1$V9 <- POLA1[i,7]
    ATAC_peak_POLA1_list = rbind(ATAC_peak_POLA1_list,ATAC_peak_POLA1)
  }
}
write.table(ATAC_peak_POLA1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLA1_list.txt",sep="\t")

ATAC_peak_POLA1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_POLA1_list.txt",sep="\t")

POLA1_features=list("POLA1" = ATAC_peak_POLA1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = POLA1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SIRT1

SIRT1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.SIRT1_age.AllCell.bed")
write.table(SIRT1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.SIRT1_age.AllCell.txt",sep = "\t")

SIRT1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.SIRT1_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SIRT1)){
  x=paste(SIRT1[i,1],SIRT1[i,2],SIRT1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SIRT1_features = list(list$V1)

ATAC_peak_SIRT1_list = data.frame()
for(i in 1:nrow(SIRT1)){
  ATAC_peak_SIRT1<- ATAC_peak[ATAC_peak$x1 ==SIRT1[i,1] & ((SIRT1[i,2] > ATAC_peak$start &SIRT1[i,2] < ATAC_peak$end) | (SIRT1[i,3] > ATAC_peak$start &SIRT1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SIRT1) != 0){
    ATAC_peak_SIRT1$V5 <-SIRT1[i,1]
    ATAC_peak_SIRT1$V6 <- SIRT1[i,2]
    ATAC_peak_SIRT1$V7 <- SIRT1[i,3]
    ATAC_peak_SIRT1$V8 <- SIRT1[i,6]
    ATAC_peak_SIRT1$V9 <- SIRT1[i,7]
    ATAC_peak_SIRT1_list = rbind(ATAC_peak_SIRT1_list,ATAC_peak_SIRT1)
  }
}
write.table(ATAC_peak_SIRT1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIRT1_list.txt",sep="\t")

ATAC_peak_SIRT1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIRT1_list.txt",sep="\t")

SIRT1_features=list("SIRT1" = ATAC_peak_SIRT1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SIRT1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##SIRT7

SIRT7= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.SIRT7_age.AllCell.bed")
write.table(SIRT7,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.SIRT7_age.AllCell.txt",sep = "\t")

SIRT7 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.SIRT7_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(SIRT7)){
  x=paste(SIRT7[i,1],SIRT7[i,2],SIRT7[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

SIRT7_features = list(list$V1)

ATAC_peak_SIRT7_list = data.frame()
for(i in 1:nrow(SIRT7)){
  ATAC_peak_SIRT7<- ATAC_peak[ATAC_peak$x1 ==SIRT7[i,1] & ((SIRT7[i,2] > ATAC_peak$start &SIRT7[i,2] < ATAC_peak$end) | (SIRT7[i,3] > ATAC_peak$start &SIRT7[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_SIRT7) != 0){
    ATAC_peak_SIRT7$V5 <-SIRT7[i,1]
    ATAC_peak_SIRT7$V6 <- SIRT7[i,2]
    ATAC_peak_SIRT7$V7 <- SIRT7[i,3]
    ATAC_peak_SIRT7$V8 <- SIRT7[i,6]
    ATAC_peak_SIRT7$V9 <- SIRT7[i,7]
    ATAC_peak_SIRT7_list = rbind(ATAC_peak_SIRT7_list,ATAC_peak_SIRT7)
  }
}
write.table(ATAC_peak_SIRT7_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIRT7_list.txt",sep="\t")

ATAC_peak_SIRT7_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_SIRT7_list.txt",sep="\t")

SIRT7_features=list("SIRT7" = ATAC_peak_SIRT7_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = SIRT7_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##TERF1

TERF1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TERF1_age.AllCell.bed")
write.table(TERF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TERF1_age.AllCell.txt",sep = "\t")

TERF1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TERF1_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TERF1)){
  x=paste(TERF1[i,1],TERF1[i,2],TERF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TERF1_features = list(list$V1)

ATAC_peak_TERF1_list = data.frame()
for(i in 1:nrow(TERF1)){
  ATAC_peak_TERF1<- ATAC_peak[ATAC_peak$x1 ==TERF1[i,1] & ((TERF1[i,2] > ATAC_peak$start &TERF1[i,2] < ATAC_peak$end) | (TERF1[i,3] > ATAC_peak$start &TERF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TERF1) != 0){
    ATAC_peak_TERF1$V5 <-TERF1[i,1]
    ATAC_peak_TERF1$V6 <- TERF1[i,2]
    ATAC_peak_TERF1$V7 <- TERF1[i,3]
    ATAC_peak_TERF1$V8 <- TERF1[i,6]
    ATAC_peak_TERF1$V9 <- TERF1[i,7]
    ATAC_peak_TERF1_list = rbind(ATAC_peak_TERF1_list,ATAC_peak_TERF1)
  }
}
write.table(ATAC_peak_TERF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TERF1_list.txt",sep="\t")

ATAC_peak_TERF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TERF1_list.txt",sep="\t")

TERF1_features=list("TERF1" = ATAC_peak_TERF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TERF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



##TERF2

TERF2= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TERF2_age.AllCell.bed")
write.table(TERF2,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TERF2_age.AllCell.txt",sep = "\t")

TERF2 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TERF2_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TERF2)){
  x=paste(TERF2[i,1],TERF2[i,2],TERF2[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TERF2_features = list(list$V1)

ATAC_peak_TERF2_list = data.frame()
for(i in 1:nrow(TERF2)){
  ATAC_peak_TERF2<- ATAC_peak[ATAC_peak$x1 ==TERF2[i,1] & ((TERF2[i,2] > ATAC_peak$start &TERF2[i,2] < ATAC_peak$end) | (TERF2[i,3] > ATAC_peak$start &TERF2[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TERF2) != 0){
    ATAC_peak_TERF2$V5 <-TERF2[i,1]
    ATAC_peak_TERF2$V6 <- TERF2[i,2]
    ATAC_peak_TERF2$V7 <- TERF2[i,3]
    ATAC_peak_TERF2$V8 <- TERF2[i,6]
    ATAC_peak_TERF2$V9 <- TERF2[i,7]
    ATAC_peak_TERF2_list = rbind(ATAC_peak_TERF2_list,ATAC_peak_TERF2)
  }
}
write.table(ATAC_peak_TERF2_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TERF2_list.txt",sep="\t")

ATAC_peak_TERF2_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TERF2_list.txt",sep="\t")

TERF2_features=list("TERF2" = ATAC_peak_TERF2_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TERF2_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)


##TOP1

TOP1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP1_age.AllCell.bed")
write.table(TOP1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP1_age.AllCell.txt",sep = "\t")

TOP1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP1_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TOP1)){
  x=paste(TOP1[i,1],TOP1[i,2],TOP1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TOP1_features = list(list$V1)

ATAC_peak_TOP1_list = data.frame()
for(i in 1:nrow(TOP1)){
  ATAC_peak_TOP1<- ATAC_peak[ATAC_peak$x1 ==TOP1[i,1] & ((TOP1[i,2] > ATAC_peak$start &TOP1[i,2] < ATAC_peak$end) | (TOP1[i,3] > ATAC_peak$start &TOP1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TOP1) != 0){
    ATAC_peak_TOP1$V5 <-TOP1[i,1]
    ATAC_peak_TOP1$V6 <- TOP1[i,2]
    ATAC_peak_TOP1$V7 <- TOP1[i,3]
    ATAC_peak_TOP1$V8 <- TOP1[i,6]
    ATAC_peak_TOP1$V9 <- TOP1[i,7]
    ATAC_peak_TOP1_list = rbind(ATAC_peak_TOP1_list,ATAC_peak_TOP1)
  }
}
write.table(ATAC_peak_TOP1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TOP1_list.txt",sep="\t")

ATAC_peak_TOP1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TOP1_list.txt",sep="\t")

TOP1_features=list("TOP1" = ATAC_peak_TOP1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TOP1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TOP2A

TOP2A= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP2A_age.AllCell.bed")
write.table(TOP2A,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP2A_age.AllCell.txt",sep = "\t")

TOP2A = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP2A_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TOP2A)){
  x=paste(TOP2A[i,1],TOP2A[i,2],TOP2A[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TOP2A_features = list(list$V1)

ATAC_peak_TOP2A_list = data.frame()
for(i in 1:nrow(TOP2A)){
  ATAC_peak_TOP2A<- ATAC_peak[ATAC_peak$x1 ==TOP2A[i,1] & ((TOP2A[i,2] > ATAC_peak$start &TOP2A[i,2] < ATAC_peak$end) | (TOP2A[i,3] > ATAC_peak$start &TOP2A[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TOP2A) != 0){
    ATAC_peak_TOP2A$V5 <-TOP2A[i,1]
    ATAC_peak_TOP2A$V6 <- TOP2A[i,2]
    ATAC_peak_TOP2A$V7 <- TOP2A[i,3]
    ATAC_peak_TOP2A$V8 <- TOP2A[i,6]
    ATAC_peak_TOP2A$V9 <- TOP2A[i,7]
    ATAC_peak_TOP2A_list = rbind(ATAC_peak_TOP2A_list,ATAC_peak_TOP2A)
  }
}
write.table(ATAC_peak_TOP2A_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TOP2A_list.txt",sep="\t")

ATAC_peak_TOP2A_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TOP2A_list.txt",sep="\t")

TOP2A_features=list("TOP2A" = ATAC_peak_TOP2A_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TOP2A_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

##TOP2B

TOP2B= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP2B_age.AllCell.bed")
write.table(TOP2B,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP2B_age.AllCell.txt",sep = "\t")

TOP2B = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.Bld.500.TOP2B_age.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(TOP2B)){
  x=paste(TOP2B[i,1],TOP2B[i,2],TOP2B[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

TOP2B_features = list(list$V1)

ATAC_peak_TOP2B_list = data.frame()
for(i in 1:nrow(TOP2B)){
  ATAC_peak_TOP2B<- ATAC_peak[ATAC_peak$x1 ==TOP2B[i,1] & ((TOP2B[i,2] > ATAC_peak$start &TOP2B[i,2] < ATAC_peak$end) | (TOP2B[i,3] > ATAC_peak$start &TOP2B[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_TOP2B) != 0){
    ATAC_peak_TOP2B$V5 <-TOP2B[i,1]
    ATAC_peak_TOP2B$V6 <- TOP2B[i,2]
    ATAC_peak_TOP2B$V7 <- TOP2B[i,3]
    ATAC_peak_TOP2B$V8 <- TOP2B[i,6]
    ATAC_peak_TOP2B$V9 <- TOP2B[i,7]
    ATAC_peak_TOP2B_list = rbind(ATAC_peak_TOP2B_list,ATAC_peak_TOP2B)
  }
}
write.table(ATAC_peak_TOP2B_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TOP2B_list.txt",sep="\t")

ATAC_peak_TOP2B_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_TOP2B_list.txt",sep="\t")

TOP2B_features=list("TOP2B" = ATAC_peak_TOP2B_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = TOP2B_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

#####BANF1
BANF1= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.ALL.05.BANF1.AllCell.bed")
write.table(BANF1,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.ALL.05.BANF1.AllCell.txt",sep = "\t")

BANF1 = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/Oth.ALL.05.BANF1.AllCell.txt",sep = "\t")

list = vector()
for (i in  1:nrow(BANF1)){
  x=paste(BANF1[i,1],BANF1[i,2],BANF1[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

BANF1_features = list(list$V1)

ATAC_peak_BANF1_list = data.frame()
for(i in 1:nrow(BANF1)){
  ATAC_peak_BANF1<- ATAC_peak[ATAC_peak$x1 ==BANF1[i,1] & ((BANF1[i,2] > ATAC_peak$start &BANF1[i,2] < ATAC_peak$end) | (BANF1[i,3] > ATAC_peak$start &BANF1[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_BANF1) != 0){
    ATAC_peak_BANF1$V5 <-BANF1[i,1]
    ATAC_peak_BANF1$V6 <- BANF1[i,2]
    ATAC_peak_BANF1$V7 <- BANF1[i,3]
    ATAC_peak_BANF1$V8 <- BANF1[i,6]
    ATAC_peak_BANF1$V9 <- BANF1[i,7]
    ATAC_peak_BANF1_list = rbind(ATAC_peak_BANF1_list,ATAC_peak_BANF1)
  }
}
write.table(ATAC_peak_BANF1_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BANF1_list.txt",sep="\t")

ATAC_peak_BANF1_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_BANF1_list.txt",sep="\t")

BANF1_features=list("BANF1" = ATAC_peak_BANF1_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = BANF1_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



#####HOTCHIRP
HOTCHIRP= import("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/macs2_sort_e3_5_50_dup3_summits_final_formerge.bed")
write.table(HOTCHIRP,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/macs2_sort_e3_5_50_dup3_summits_final_formerge.txt",sep = "\t")

HOTCHIRP = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/macs2_sort_e3_5_50_dup3_summits_final_formerge.txt",sep = "\t")



list = vector()
for (i in  1:nrow(HOTCHIRP)){
  x=paste(HOTCHIRP[i,1],HOTCHIRP[i,2],HOTCHIRP[i,3],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

HOTCHIRP_features = list(list$V1)

ATAC_peak_HOTCHIRP_list = data.frame()
for(i in 1:nrow(HOTCHIRP)){
  ATAC_peak_HOTCHIRP<- ATAC_peak[ATAC_peak$x1 ==HOTCHIRP[i,1] & ((HOTCHIRP[i,2] > ATAC_peak$start &HOTCHIRP[i,2] < ATAC_peak$end) | (HOTCHIRP[i,3] > ATAC_peak$start &HOTCHIRP[i,3] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_HOTCHIRP) != 0){
    ATAC_peak_HOTCHIRP$V5 <-HOTCHIRP[i,1]
    ATAC_peak_HOTCHIRP$V6 <- HOTCHIRP[i,2]
    ATAC_peak_HOTCHIRP$V7 <- HOTCHIRP[i,3]
    ATAC_peak_HOTCHIRP$V8 <- HOTCHIRP[i,6]
    ATAC_peak_HOTCHIRP$V9 <- HOTCHIRP[i,7]
    ATAC_peak_HOTCHIRP_list = rbind(ATAC_peak_HOTCHIRP_list,ATAC_peak_HOTCHIRP)
  }
}
write.table(ATAC_peak_HOTCHIRP_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HOTCHIRP_list.txt",sep="\t")

ATAC_peak_HOTCHIRP_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_HOTCHIRP_list.txt",sep="\t")
HOTCHIRP_features=list("HOTCHIRP" = ATAC_peak_HOTCHIRP_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = HOTCHIRP_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)

save(scmulti2,file = "/Users/ramzipit/Desktop/scmulti2_temp.Rdata")


#######triplex
#LINC01578tri
LINC01578tri = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/triplex/LINC01578/TE/LINC01578_TTS_2.txt",sep ="\t")

list = vector()
for (i in  1:nrow(LINC01578tri)){
  x=paste(LINC01578tri[i,4],LINC01578tri[i,5],LINC01578tri[i,6],sep='-')
  list = cbind(list,x)
}
list = as.data.frame(t(list))

LINC01578tri_features = list(list$V1)

ATAC_peak_LINC01578tri_list = data.frame()
for(i in 1:nrow(LINC01578tri)){
  ATAC_peak_LINC01578tri<- ATAC_peak[ATAC_peak$x1 ==LINC01578tri[i,4] & ((LINC01578tri[i,5] > ATAC_peak$start & LINC01578tri[i,5] < ATAC_peak$end) | (LINC01578tri[i,6] > ATAC_peak$start &LINC01578tri[i,6] < ATAC_peak$end)),]
  if (nrow(ATAC_peak_LINC01578tri) != 0){
    ATAC_peak_LINC01578tri$V5 <-LINC01578tri[i,4]
    ATAC_peak_LINC01578tri$V6 <- LINC01578tri[i,5]
    ATAC_peak_LINC01578tri$V7 <- LINC01578tri[i,6]
    ATAC_peak_LINC01578tri_list = rbind(ATAC_peak_LINC01578tri_list,ATAC_peak_LINC01578tri)
  }
}
write.table(ATAC_peak_LINC01578tri_list,"/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_LINC01578tri_list.txt",sep="\t")

ATAC_peak_LINC01578tri_list = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/CHIP-bed/ATAC_peak_LINC01578tri_list.txt",sep="\t")

LINC01578tri_features=list("LINC01578tri" = ATAC_peak_LINC01578tri_list$ATAC_peak)

DefaultAssay(scmulti2) <- 'peaks'
scmulti2 = AddChromatinModule(
  object = scmulti2,
  features = LINC01578tri_features,
  genome = BSgenome.Hsapiens.UCSC.hg38,
)



