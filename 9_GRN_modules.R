library(ComplexHeatmap)
library(circlize)

#####
NT_0.05ATAC_p = read.table("NT_0.05ATACmatrix_p.txt",sep = "\t", header = T)
colnames(NT_0.05ATAC_rho) = rownames(NT_0.05ATAC_rho)
NT_0.05ATAC_rho = read.table("NT_0.05ATACmatrix_rho.txt",sep = "\t", header = T)
colnames(NT_0.05ATAC_rho) = rownames(NT_0.05ATAC_rho)


#motif.list 640
motif.list2 = as.data.frame(unique(motif.list))
colnames(motif.list2)[1] = "motif.list"
NT_0.05ATAC_p = cbind(rownames(NT_0.05ATAC_p),NT_0.05ATAC_p)
colnames(NT_0.05ATAC_p)[1] = "motif.list"

library(dplyr)
NT_0.05ATAC_p = left_join(motif.list2,NT_0.05ATAC_p,by='motif.list')
rownames(NT_0.05ATAC_p) = NT_0.05ATAC_p$motif.list
NT_0.05ATAC_p = NT_0.05ATAC_p[,-1]

NT_0.05ATAC_p = as.data.frame(t(NT_0.05ATAC_p))
NT_0.05ATAC_p = cbind(rownames(NT_0.05ATAC_p),NT_0.05ATAC_p)
colnames(NT_0.05ATAC_p)[1] = "motif.list"
NT_0.05ATAC_p = left_join(motif.list2,NT_0.05ATAC_p,by='motif.list')
rownames(NT_0.05ATAC_p) = NT_0.05ATAC_p$motif.list
NT_0.05ATAC_p = NT_0.05ATAC_p[,-1]
NT_0.05ATAC_p2 = NT_0.05ATAC_p
NT_0.05ATAC_p2 = as.data.frame(t(NT_0.05ATAC_p2))

#motif.list 640
motif.list2 = as.data.frame(unique(motif.list))
colnames(motif.list2)[1] = "motif.list"
NT_0.05ATAC_rho = cbind(rownames(NT_0.05ATAC_rho),NT_0.05ATAC_rho)
colnames(NT_0.05ATAC_rho)[1] = "motif.list"

library(dplyr)
NT_0.05ATAC_rho = left_join(motif.list2,NT_0.05ATAC_rho,by='motif.list')
rownames(NT_0.05ATAC_rho) = NT_0.05ATAC_rho$motif.list
NT_0.05ATAC_rho = NT_0.05ATAC_rho[,-1]

NT_0.05ATAC_rho = as.data.frame(t(NT_0.05ATAC_rho))
NT_0.05ATAC_rho = cbind(rownames(NT_0.05ATAC_rho),NT_0.05ATAC_rho)
colnames(NT_0.05ATAC_rho)[1] = "motif.list"
NT_0.05ATAC_rho = left_join(motif.list2,NT_0.05ATAC_rho,by='motif.list')
rownames(NT_0.05ATAC_rho) = NT_0.05ATAC_rho$motif.list
NT_0.05ATAC_rho = NT_0.05ATAC_rho[,-1]
NT_0.05ATAC_rho2 = NT_0.05ATAC_rho
NT_0.05ATAC_rho2 = as.data.frame(t(NT_0.05ATAC_rho2))

NT_ATAC05ind = read.table(file = "NT_ATAC05_DAmotifs_ind.txt",sep =",")
NT_ATAC05ind[1,1] = "254"
NT_ATAC05ind[1,ncol(NT_ATAC05ind)] = "289"
NT_ATAC05ind = as.numeric(NT_ATAC05ind)
NT_ATAC05ind = as.data.frame(NT_ATAC05ind+1)

colnames(NT_ATAC05ind) = 'order'
order = c(1:430)
NT_0.05ATACDA_rho = as.data.frame(cbind(order,rownames(NT_0.05ATAC_rho2),NT_0.05ATAC_rho2))

library(dplyr)
NT_0.05ATACDA_rho2 = left_join(NT_ATAC05ind,NT_0.05ATACDA_rho,by='order')
NT_0.05ATACDA_rho3 = NT_0.05ATACDA_rho2
rownames(NT_0.05ATACDA_rho3) = NT_0.05ATACDA_rho3$`rownames(NT_0.05ATAC_rho2)`
NT_0.05ATACDA_rho3 = NT_0.05ATACDA_rho3[,-1]
NT_0.05ATACDA_rho3 = NT_0.05ATACDA_rho3[,-1]

NT_0.05ATACDA_rho3 = as.data.frame(t(rbind(order,NT_0.05ATACDA_rho3)))
colnames(NT_0.05ATACDA_rho3)[1] = 'order'
NT_0.05ATACDA_rho3 = left_join(NT_ATAC05ind,NT_0.05ATACDA_rho3,by='order')
rownames(NT_0.05ATACDA_rho3) = NT_0.05ATACDA_rho2$`rownames(NT_0.05ATAC_rho2)`
NT_0.05ATACDA_rho3 = NT_0.05ATACDA_rho3[,-1]
library(ComplexHeatmap)
library(circlize)

NT_0.05ATACDA_rho4 = NT_0.05ATACDA_rho3[rownames(NT_0.05ATACDA_rho3) != "LINC01638tri",colnames(NT_0.05ATACDA_rho3) != "LINC01638tri"]
NT_0.05ATACDA_rho4 = NT_0.05ATACDA_rho4[rownames(NT_0.05ATACDA_rho4) != "ZNF274",colnames(NT_0.05ATACDA_rho4) != "ZNF274"]


write.table(NT_0.05ATACDA_rho4,file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/GRN/ATAC/NT_0.05ATACDA_rho_cluster_without0.txt",sep = "\t")


cluster = c(rep('module1',51),rep('module2',24),rep('module3',52),rep('module4',51),rep('module5',53),
            rep('module6',45),rep('module7',87),rep('module8',23),rep('module9',42))#428
#1:51 52:75 76:127 128:178 179:231 232:276 277:363 364:386 387:428

col_fun = colorRamp2(c(-0.85,-0.5,-0.3,0,0.3,0.5,0.85), c("#3B2061",
                                                        "#694B9A",
                                                        "#9A90BE",
                                                        "#F7F3ED",
                                                        "#FBC075",
                                                        "#ED9B3D","#B75E25"))
#col_fun = colorRamp2(c(-0.8,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.8), c("#312584","#3C40C3","#2E6DD2","#1D90C0","#0DAA9E","#47B967","#A9AD21","#CFD61D","#CFE017"))

col_fun(seq(-0.85, 0.85,0.05),)

#27,75,140,202,297,425,470,573,613
Heatmap(as.matrix(NT_0.05ATACDA_rho4), 
        row_split = cluster, 
        #column_split = c(27,48,65,62,95,128,45,103,40),
        #column_split =c(letters[1:27],letters[28:75],letters[76:140],letters[141:202],letters[203:297],letters[298:425],letters[426:470],letters[471:573],letters[574:613]),
        column_split = cluster, 
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        border = TRUE,
        col = col_fun,
        use_raster = F,
        heatmap_legend_param = list(col_fun = col_fun, title = "Correlation",title_position = "leftcenter-rot",at = c(-0.85,0.85),
                                    labels = c("-0.85", "0.85"),border="black",legend_height = unit(20, "mm"),legend_width = unit(1, "mm"),
                                    labels_gp = gpar(col = "black", font = 2),#direction = "horizontal",
                                    title_position = "topcenter")
        
)

#
ATAC_1module = NT_0.05ATACDA_rho_6module[1:24,]
age = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/genage_human.txt",sep = '\t',header = T)
sene = read.table(file = "/Users/ramzipit/Desktop/lnc_project/bioinfo_results/seurat/ATAC/annotation/sene_senesig_merge.txt",sep ="\t",header = T,quote = "")

age_1module = intersect(age$symbol,rownames(ATAC_1module))
sene_1module = intersect(sene$Gene.symbol,rownames(ATAC_1module))

#
ATAC_2module = NT_0.05ATACDA_rho_6module[25:77,]
age_2module = intersect(age$symbol,rownames(ATAC_2module))
sene_2module = intersect(sene$Gene.symbol,rownames(ATAC_2module))

#
ATAC_3module = NT_0.05ATACDA_rho_6module[78:122,]
age_3module = intersect(age$symbol,rownames(ATAC_3module))
sene_3module = intersect(sene$Gene.symbol,rownames(ATAC_3module))

#
ATAC_4module = NT_0.05ATACDA_rho_6module[123:209,]
age_4module = intersect(age$symbol,rownames(ATAC_4module))
sene_4module = intersect(sene$Gene.symbol,rownames(ATAC_4module))

#
ATAC_5module = NT_0.05ATACDA_rho_6module[210:232,]
age_5module = intersect(age$symbol,rownames(ATAC_5module))
sene_5module = intersect(sene$Gene.symbol,rownames(ATAC_5module))

#
ATAC_6module = NT_0.05ATACDA_rho_6module[233:274,]
age_6module = intersect(age$symbol,rownames(ATAC_6module))
sene_6module = intersect(sene$Gene.symbol,rownames(ATAC_6module))

###ATAC&RNA integration
library(Seurat)

NT_RNAmatrix_mid05_rho_7module = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/GRN/RNA_important/NT_RNAmatrix_mid05_rho_7module.txt",sep = "\t")

module1 = rownames(NT_RNAmatrix_mid05_rho_7module)[1:247]
module2 = rownames(NT_RNAmatrix_mid05_rho_7module)[248:456]
module3 = rownames(NT_RNAmatrix_mid05_rho_7module)[457:646]
module4 = rownames(NT_RNAmatrix_mid05_rho_7module)[647:723]
module5 = rownames(NT_RNAmatrix_mid05_rho_7module)[724:1036]
module6 = rownames(NT_RNAmatrix_mid05_rho_7module)[1037:1145]
module7 = rownames(NT_RNAmatrix_mid05_rho_7module)[1146:1249]

scmulti4 = scmulti2
DefaultAssay(scmulti4) <- "RNA"
RNA_module1 <- list(module1)
scmulti4 <- AddModuleScore(scmulti4,
                          features = RNA_module1,
                          ctrl = 100,
                          name = "RNA_module_one")

RNA_module2 <- list(module2)
scmulti4 <- AddModuleScore(scmulti4,
                           features = RNA_module2,
                           ctrl = 100,
                           name = "RNA_module_two")

RNA_module3 <- list(module3)
scmulti4 <- AddModuleScore(scmulti4,
                           features = RNA_module3,
                           ctrl = 100,
                           name = "RNA_module_three")

RNA_module4 <- list(module4)
scmulti4 <- AddModuleScore(scmulti4,
                           features = RNA_module4,
                           ctrl = 100,
                           name = "RNA_module_four")

RNA_module5 <- list(module5)
scmulti4 <- AddModuleScore(scmulti4,
                           features = RNA_module5,
                           ctrl = 100,
                           name = "RNA_module_five")

RNA_module6 <- list(module6)
scmulti4 <- AddModuleScore(scmulti4,
                           features = RNA_module6,
                           ctrl = 100,
                           name = "RNA_module_six")

RNA_module7 <- list(module7)
scmulti4 <- AddModuleScore(scmulti4,
                           features = RNA_module7,
                           ctrl = 100,
                           name = "RNA_module_seven")

####
NT_0.05ATACDA_rho_6module  = read.table("/Users/ramzipit/Desktop/lnc_project/bioinfo_results/GRN/ATAC/NT_0.05ATACDA_rho_6module.txt",sep = "\t")

ATACmodule1 = rownames(NT_0.05ATACDA_rho_6module)[1:24]
ATACmodule2 = rownames(NT_0.05ATACDA_rho_6module)[25:77]
ATACmodule3 = rownames(NT_0.05ATACDA_rho_6module)[78:122]
ATACmodule4 = rownames(NT_0.05ATACDA_rho_6module)[123:209]
ATACmodule5 = rownames(NT_0.05ATACDA_rho_6module)[210:232]
ATACmodule6 = rownames(NT_0.05ATACDA_rho_6module)[233:274]


DefaultAssay(scmulti4) <- "newcharmVAR"
ATAC_module1 <- list(ATACmodule1)
scmulti4 <- AddModuleScore(scmulti4,
                           features = ATAC_module1,
                           ctrl = 20,
                           name = "ATAC_module_one")

ATAC_module2 <- list(ATACmodule2)
scmulti4 <- AddModuleScore(scmulti4,
                           features = ATAC_module2,
                           ctrl = 20,
                           name = "ATAC_module_two")

ATAC_module3 <- list(ATACmodule3)
scmulti4 <- AddModuleScore(scmulti4,
                           features = ATAC_module3,
                           ctrl = 20,
                           name = "ATAC_module_three")

ATAC_module4 <- list(ATACmodule4)
scmulti4 <- AddModuleScore(scmulti4,
                           features = ATAC_module4,
                           ctrl = 20,
                           name = "ATAC_module_four")

ATAC_module5 <- list(ATACmodule5)
scmulti4 <- AddModuleScore(scmulti4,
                           features = ATAC_module5,
                           ctrl = 20,
                           name = "ATAC_module_five")

ATAC_module6 <- list(ATACmodule6)
scmulti4 <- AddModuleScore(scmulti4,
                           features = ATAC_module6,
                           ctrl = 20,
                           name = "ATAC_module_six")
####
library(magrittr)

scmulti4_meta3 <- scmulti4@meta.data %>%
  t() %>%
  as.data.frame() %>%
  { .[332:nrow(.), , drop = FALSE] } %>%
  { replace(., is.na(.), 0) } %>%
  as.matrix() %>%
  apply(2, as.numeric)

rownames(scmulti4_meta3) <- rownames(scmulti4@meta.data)[332:nrow(scmulti4@meta.data)]
scmulti4[["module"]] <- CreateAssayObject(data = scmulti4_meta3)




##
all_control_0.05
NT05_module = as.data.frame(all_control_0.05@assays$module@data)

library(psych)

t_NT05_module  = t(NT05_module)
y1 = corr.test(t_NT05_module)

t_NT05_module_rho = y1$r
t_NT05_module_rho[is.na(t_NT05_module_rho)] = 0

t_NT05_module_rho2 = t_NT05_module_rho[1:7,8:13]


t_NT05_module_p = y1$p
t_NT05_module_p2 = t_NT05_module_p[1:7,8:13]

t_NT05_module_padj = y1$p.adj
t_NT05_module_padj2 = t_NT05_module_padj[1:7,8:13]

##
library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(-1,-0.7,-0.4,0,0.4,0.7,1), c("#124985","#256DAF","#5BA2CB","#F3F4F5","#F5AE8D","#D35A4A","#670020"))
#col_fun = colorRamp2(c(-0.8,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.8), c("#312584","#3C40C3","#2E6DD2","#1D90C0","#0DAA9E","#47B967","#A9AD21","#CFD61D","#CFE017"))

col_fun(seq(-1, 1,0.05),)

Heatmap(as.matrix(t_NT05_module_rho2), 
        #row_split = cluster, 
        #column_split = cluster, 
        cluster_rows = T,
        cluster_columns = T,
        show_column_names = T,
        show_row_names = T,
        border = TRUE,
        col = col_fun,
        use_raster = F,
        heatmap_legend_param = list(col_fun = col_fun, title = "Correlation",title_position = "leftcenter-rot",at = c(-1,1),
                                    labels = c("-1", "1"),border="black",legend_height = unit(20, "mm"),legend_width = unit(1, "mm"),
                                    labels_gp = gpar(col = "black", font = 2),#direction = "horizontal",
                                    title_position = "topcenter")
)


####
gene <- c("POU1F1", "HDAC3","BRCA1","SP1","TFAP2A","CTNNB1","NFE2L1","NEAT1tri","H3K27ac",
          "BRCA1",  "CTNNB1", "EHF", "ETV6" , "HDAC3" ,"MEF2A","OTX2" ,
          "STAT5A", "HDAC3" ,"NCOR1", "NBN", "EGR1","AR","FOXM1","TBP","H3K27me3", "LINC00578tri",
          "AR","ETS1" ,"FOXM1" ,"HDAC3" ,"NBN", "PHB2" ,
          "NBN" ,"JUND","FOS" ,"TOP2A","TOP2B" ,"TAF1" ,"ATF2","TBP", "HBP1", "TFDP1" ,"DDIT3", "SIRT6" , "ARNTL" , "NCOR2" , "TP73","NFE2L2","H3K79me2" ,"H3K27me3",
          "ARID1B","CBX8","FOS" ,"HBP1","JUNB","MYCN" ,"NBN","NFE2L2",
          "NFKB2" , "STAT5A", "HDAC3",  "NCOR1" , "HIF1A",  "NFKB1" , "RELA"  , "JUN" ,   "HESX1" , "CREB1" , "DDIT3", 
          "ARID1B", "ARID4B", "EZH2", "HDAC3" ,"HIF1A" ,"IRF1","JUN" ,"KDM1A" ,"MNX1" ,"NFKB2","NRF1" ,
          "TP53", "JUND" ,"FOS" ,
          "E2F7","FOS","FOXA1" ,"JUNB" ,
          "STAT5B", "ESR1","NFE2L2",
          "ESR1", "ESRRB" ,"KDM1A","NFE2L2","H3K4me3",
          "NBN",   "HIF1A", "EGR1","SP1"  ,"H3K4me1", "H3K9ac",
          "BHLHE40","EHF","HIF1A","NBN" ,
          "STAT3",  "STAT5A", "NBN" ,"EP300",  "AR","HDAC1" ,"POLA1" ,
          "AR"  ,  "EP300" ,"EZH2" , "HDAC1", "NBN"  ,
          "NFKB2",  "MYC"  ,  "MXI1" ,  "TFAP2A" ,"HBP1"  ,
          "ESRRB", "HBP1" , "MYC" ,  "NFKB2",
          "ARNTL",
          "BCLAF1",
          "FOXO3",  "SP1","FOXO4" ,"HIC1" ,"NFE2L2","H3K9me3",
          "ARID4B", "EHMT2",  "ETS2"  ,"FOXD1" , "FOXO3" , "FOXO4" , "FOXP1" , "FOXP3" , "MEF2A" , "NFE2L2",
          "EGR1",  "FOXO1" ,"HDAC1" ,"TAF1",
          "FOXO1", "HDAC1", "MEOX1" ,"MNX1" ,
          "BRD4",  "EZH2" , "FOXA1","H3K79me2",  "H3K9ac", "BRD4",
          "FOXM1", "HDAC1" ,"HDAC2", "ARNTL","H3K9ac",  "H3K4me3",
          "BCLAF1", "E2F7", "FOXA1", "FOXM1","HDAC1","HDAC2","MITF","MNX1","CLOCK1",
          "TP53",   "HDAC3" , "RB1"  ,  "MAPK14" , "PPARG" ,
          "ARID1B", "E2F7"  , "HDAC3" , "MAPK14" ,"PPARG" ,
          "CTNNB1",
          "CTNNB1", "E2F3"  ,"LINC00578tri", "LINC02693tri",
          "CTCF","MAX::MYC","FOS::JUNB"
          
          
)
gene = unique(gene)
gene = as.data.frame(gene)

chromvar = as.data.frame(all_control_0.05@assays$newcharmVAR@data)
chromvar = cbind(rownames(chromvar),chromvar)
colnames(chromvar)[1] = 'gene'

age_chromvar = left_join(gene,chromvar,by = 'gene')
rownames(age_chromvar) = age_chromvar[,1]
age_chromvar = age_chromvar[,-1]

row = vector()
inte = data.frame()
for(i in 1:7){
  for(j in 1:nrow(age_chromvar)){
    x = cor.test(as.numeric(NT05_module[i,]),as.numeric(age_chromvar[j,]),method = "pearson")
    xc = x$estimate
    xp = x$p.value
    row = c(rownames(NT05_module[i,]),rownames(age_chromvar[j,]),xc,xp)
    inte = rbind(inte,row)
  }
}


RNAmodule1_inte = inte[inte$X.RNA.module.one1. == "RNA-module-one1",]
rownames(RNAmodule1_inte) = RNAmodule1_inte$X.POU1F1.
RNAmodule1_inte2 = as.data.frame(as.numeric(unlist(RNAmodule1_inte[,3])))
rownames(RNAmodule1_inte2) = rownames(RNAmodule1_inte)
colnames(RNAmodule1_inte2) = 'RNAmodule1'

RNAmodule2_inte = inte[inte$X.RNA.module.one1. == "RNA-module-two1",]
rownames(RNAmodule2_inte) = RNAmodule2_inte$X.POU1F1.
RNAmodule2_inte2 = as.data.frame(as.numeric(unlist(RNAmodule2_inte[,3])))
rownames(RNAmodule2_inte2) = rownames(RNAmodule2_inte)
colnames(RNAmodule2_inte2) = 'RNAmodule2'

RNAmodule3_inte = inte[inte$X.RNA.module.one1. == "RNA-module-three1",]
rownames(RNAmodule3_inte) = RNAmodule3_inte$X.POU1F1.
RNAmodule3_inte2 = as.data.frame(as.numeric(unlist(RNAmodule3_inte[,3])))
rownames(RNAmodule3_inte2) = rownames(RNAmodule3_inte)
colnames(RNAmodule3_inte2) = 'RNAmodule3'

RNAmodule4_inte = inte[inte$X.RNA.module.one1. == "RNA-module-four1",]
rownames(RNAmodule4_inte) = RNAmodule4_inte$X.POU1F1.
RNAmodule4_inte2 = as.data.frame(as.numeric(unlist(RNAmodule4_inte[,3])))
rownames(RNAmodule4_inte2) = rownames(RNAmodule4_inte)
colnames(RNAmodule4_inte2) = 'RNAmodule4'

RNAmodule5_inte = inte[inte$X.RNA.module.one1. == "RNA-module-five1",]
rownames(RNAmodule5_inte) = RNAmodule5_inte$X.POU1F1.
RNAmodule5_inte2 = as.data.frame(as.numeric(unlist(RNAmodule5_inte[,3])))
rownames(RNAmodule5_inte2) = rownames(RNAmodule5_inte)
colnames(RNAmodule5_inte2) = 'RNAmodule5'

RNAmodule6_inte = inte[inte$X.RNA.module.one1. == "RNA-module-six1",]
rownames(RNAmodule6_inte) = RNAmodule6_inte$X.POU1F1.
RNAmodule6_inte2 = as.data.frame(as.numeric(unlist(RNAmodule6_inte[,3])))
rownames(RNAmodule6_inte2) = rownames(RNAmodule6_inte)
colnames(RNAmodule6_inte2) = 'RNAmodule6'

RNAmodule7_inte = inte[inte$X.RNA.module.one1. == "RNA-module-seven1",]
rownames(RNAmodule7_inte) = RNAmodule7_inte$X.POU1F1.
RNAmodule7_inte2 = as.data.frame(as.numeric(unlist(RNAmodule7_inte[,3])))
rownames(RNAmodule7_inte2) = rownames(RNAmodule7_inte)
colnames(RNAmodule7_inte2) = 'RNAmodule7'

RNAmodule_inte = as.data.frame(cbind(RNAmodule1_inte2,RNAmodule2_inte2,RNAmodule3_inte2,
                                     RNAmodule4_inte2,RNAmodule5_inte2,RNAmodule6_inte2,
                                     RNAmodule7_inte2))
RNAmodule_inte2 = apply(RNAmodule_inte,2,as.numeric)
rownames(RNAmodule_inte2) = rownames(RNAmodule_inte)

col_fun = colorRamp2(c(-1,-0.5,-0.3,0,0.3,0.5,1), c("#124985","#256DAF","#5BA2CB","#F3F4F5","#F5AE8D","#D35A4A","#670020"))
#col_fun = colorRamp2(c(-0.8,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.8), c("#312584","#3C40C3","#2E6DD2","#1D90C0","#0DAA9E","#47B967","#A9AD21","#CFD61D","#CFE017"))

col_fun(seq(-1, 1,0.05),)

Heatmap(as.matrix(RNAmodule_inte2), 
        #row_split = cluster, 
        #column_split = cluster, 
        cluster_rows = T,
        cluster_columns = T,
        show_column_names = T,
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 8),
        show_row_names = T,
        border = TRUE,
        col = col_fun,
        use_raster = F,
        heatmap_legend_param = list(col_fun = col_fun, title = "Correlation",title_position = "leftcenter-rot",at = c(-1,1),
                                    labels = c("-1", "1"),border="black",legend_height = unit(20, "mm"),legend_width = unit(1, "mm"),
                                    labels_gp = gpar(col = "black", font = 2),#direction = "horizontal",
                                    title_position = "topcenter")
)
