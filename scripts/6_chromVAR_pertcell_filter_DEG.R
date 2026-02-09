library(ggplot2)
library(Seurat)


library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(Seurat)

DefaultAssay(object = T2K2) <- "RNA"
T2K2$replicate = 'rep1'
T2K2 <- CalcPerturbSig(object = T2K2,
                           assay = "RNA",
                           slot = "data",
                           gd.class ="perturb",
                           nt.cell.class = "NT",
                           reduction = "pca",
                           split.by = "replicate",
                           ndims = 15,
                           num.neighbors = 10,
                           new.assay.name = "PRTB_RNA",
)

T2K2 <- RunMixscape(
  object = T2K2, 
  assay = "PRTB", 
  slot = "scale.data", 
  labels = "perturb", 
  nt.class.name = "NT", 
  min.de.genes = 2, 
  min.cells = 5,
  iter.num = 10, 
  de.assay = "RNA", 
  verbose = F,
  prtb.type = "KD",
  logfc.threshold = 0.1)

library(Matrix)
library(matrixStats)

T2K2_test = T2K2
object = T2K2_test
assay = "PRTB"
slot = "scale.data"
labels = "perturb"
nt.class.name = "NT"
new.class.name = "mixscape_class"
min.de.genes = 2
min.cells = 5
de.assay = "RNA"
logfc.threshold = 0.1
iter.num = 10
verbose = FALSE
split.by = NULL
fine.mode = FALSE
fine.mode.labels = "guide_ID"
prtb.type = "KD"
#DefaultAssay(object) = 'RNA'
{
  mixtools.installed <- PackageCheck("mixtools", error = FALSE)
  if (!mixtools.installed[1]) {
    stop("Please install the mixtools package to use RunMixscape", 
         "\nThis can be accomplished with the following command: ", 
         "\n----------------------------------------", "\ninstall.packages('mixtools')", 
         "\n----------------------------------------", call. = FALSE)
  }
  assay <- assay %||% DefaultAssay(object = object)
  if (is.null(x = labels)) {
    stop("Please specify target gene class metadata name")
  }
  prtb_markers <- list()
  object[[new.class.name]] <- object[[labels]]
  object[[new.class.name]][, 1] <- as.character(x = object[[new.class.name]][, 
                                                                             1])
  object[[paste0(new.class.name, "_p_", tolower(x = prtb.type))]] <- 0
  gv.list <- list()
  if (is.null(x = split.by)) {
    split.by <- splits <- "con1"
  }else {splits <- as.character(x = unique(x = object[[split.by]][, 
                                                                  1]))}
  cells.s.list <- list()
  for (s in splits) {
    Idents(object = object) <- split.by
    cells.s <- WhichCells(object = object, idents = s)
    cells.s.list[[s]] <- cells.s
    genes <- setdiff(x = unique(x = object[[labels]][cells.s, 
                                                     1]), y = nt.class.name)
    Idents(object = object) <- labels
    for (gene in genes) {
      if (isTRUE(x = verbose)) {
        message("Processing ", gene)
      }
      orig.guide.cells <- intersect(x = WhichCells(object = object, 
                                                   idents = gene), y = cells.s)
      nt.cells <- intersect(x = WhichCells(object = object, 
                                           idents = nt.class.name), y = cells.s)
      if (isTRUE(x = fine.mode)) {
        guides <- setdiff(x = unique(x = object[[fine.mode.labels]][orig.guide.cells, 
                                                                    1]), y = nt.class.name)
        all.de.genes <- c()
        for (gd in guides) {
          gd.cells <- rownames(x = object[[]][orig.guide.cells, 
          ])[which(x = object[[]][orig.guide.cells, 
                                  fine.mode.labels] == gd)]
          de.genes <- TopDEGenesMixscape(object = object, 
                                         ident.1 = gd.cells, ident.2 = nt.cells, de.assay = de.assay, 
                                         logfc.threshold = logfc.threshold, labels = fine.mode.labels, 
                                         verbose = verbose)
          all.de.genes <- c(all.de.genes, de.genes)
        }
        all.de.genes <- unique(all.de.genes)
      }
      else {
        all.de.genes <- TopDEGenesMixscape(object = object, 
                                           ident.1 = orig.guide.cells, ident.2 = nt.cells, 
                                           de.assay = de.assay, logfc.threshold = logfc.threshold, 
                                           labels = labels, verbose = verbose)
      }
      prtb_markers[[s]][[gene]] <- all.de.genes
      if (length(x = all.de.genes) < min.de.genes) {
        prtb_markers[[s]][[gene]] <- character()
      }
    }
  }
  all_markers <- unique(x = unlist(x = prtb_markers))
  missing_genes <- all_markers[!all_markers %in% rownames(x = object[[assay]])]
  object <- GetMissingPerturb(object = object, assay = assay, 
                              features = missing_genes, verbose = verbose)
  for (s in splits) {
    cells.s <- cells.s.list[[s]]
    genes <- setdiff(x = unique(x = object[[labels]][cells.s, 
                                                     1]), y = nt.class.name)
    if (verbose) {
      message("Classifying cells for: ")
    }
    for (gene in genes) {
      Idents(object = object) <- labels
      post.prob <- 0
      orig.guide.cells <- intersect(x = WhichCells(object = object, 
                                                   idents = gene), y = cells.s)
      nt.cells <- intersect(x = WhichCells(object = object, 
                                           idents = nt.class.name), y = cells.s)
      all.cells <- c(orig.guide.cells, nt.cells)
      if (length(x = prtb_markers[[s]][[gene]]) == 0) {
        if (verbose) {
          message("  Fewer than ", min.de.genes, " DE genes for ", 
                  gene, ". Assigning cells as NP.")
        }
        object[[new.class.name]][orig.guide.cells, 1] <- paste0(gene, 
                                                                " NP")
      }
      else {
        if (verbose) {
          message("  ", gene)
        }
        de.genes <- prtb_markers[[s]][[gene]]
        dat <- GetAssayData(object = object[[assay]], 
                            slot = "data")[de.genes, all.cells, drop = FALSE]
        if (slot == "scale.data") {
          dat <- ScaleData(object = dat, features = de.genes, 
                           verbose = FALSE)
        }
        converged <- FALSE
        n.iter <- 0
        old.classes <- object[[new.class.name]][all.cells, 
        ]
        while (!converged && n.iter < iter.num) {
          Idents(object = object) <- new.class.name
          guide.cells <- intersect(x = WhichCells(object = object, 
                                                  idents = gene), y = cells.s)
          vec <- rowMeans2(x = dat[, guide.cells, drop = FALSE]) - 
            rowMeans2(x = dat[, nt.cells, drop = FALSE])
          pvec <- apply(X = dat, MARGIN = 2, FUN = ProjectVec, 
                        v2 = vec)
          if (n.iter == 0) {
            gv <- as.data.frame(x = pvec)
            gv[, labels] <- nt.class.name
            gv[intersect(x = rownames(x = gv), y = guide.cells), 
               labels] <- gene
            gv.list[[gene]][[s]] <- gv
          }
          guide.norm <- DefineNormalMixscape(pvec[guide.cells])
          nt.norm <- DefineNormalMixscape(pvec[nt.cells])
          mm <- mixtools::normalmixEM(x = pvec, mu = c(nt.norm$mu, 
                                                       guide.norm$mu), sigma = c(nt.norm$sd, guide.norm$sd), 
                                      k = 2, mean.constr = c(nt.norm$mu, NA), sd.constr = c(nt.norm$sd, 
                                                                                            NA), verb = FALSE, maxit = 5000, maxrestarts = 100)
          lik.ratio <- dnorm(x = pvec[orig.guide.cells], 
                             mean = mm$mu[1], sd = mm$sigma[1])/dnorm(x = pvec[orig.guide.cells], 
                                                                      mean = mm$mu[2], sd = mm$sigma[2])
          post.prob <- 1/(1 + lik.ratio)
          object[[new.class.name]][names(x = which(post.prob > 
                                                     0.5)), 1] <- gene
          object[[new.class.name]][names(x = which(post.prob < 
                                                     0.5)), 1] <- paste(gene, " NP", sep = "")
          if (length(x = which(x = object[[new.class.name]] == 
                               gene & Cells(x = object) %in% cells.s)) < 
              min.de.genes) {
            if (verbose) {
              message("Fewer than ", min.cells, " cells assigned as ", 
                      gene, "Assigning all to NP.")
            }
            object[[new.class.name]][guide.cells, 1] <- "NP"
            converged <- TRUE
          }
          if (all(object[[new.class.name]][all.cells, 
          ] == old.classes)) {
            converged <- TRUE
          }
          old.classes <- object[[new.class.name]][all.cells, 
          ]
          n.iter <- n.iter + 1
        }
        object[[new.class.name]][which(x = object[[new.class.name]] == 
                                         gene & Cells(x = object) %in% cells.s), 1] <- paste(gene, 
                                                                                             prtb.type, sep = " ")
      }
      object[[paste0(new.class.name, ".global")]] <- as.character(x = sapply(X = as.character(x = object[[new.class.name]][, 
                                                                                                                           1]), FUN = function(x) {
                                                                                                                             strsplit(x = x, split = " (?=[^ ]+$)", perl = TRUE)[[1]][2]
                                                                                                                           }))
      object[[paste0(new.class.name, ".global")]][which(x = is.na(x = object[[paste0(new.class.name, 
                                                                                     ".global")]])), 1] <- nt.class.name
      object[[paste0(new.class.name, "_p_", tolower(prtb.type))]][names(x = post.prob), 
                                                                  1] <- post.prob
    }
  }
  
  ProcessTool <- function(object, gv.list) {
    Tool(object = object) <- gv.list
    return(object)
  }
  
  object <- ProcessTool(object = object, gv.list = gv.list)
  
  
  #Tool(object = object) <- gv.list
  Idents(object = object) <- new.class.name
  return(object)
}

T2K2_test = object

p_kd <- as.vector(T2K2_test$mixscape_class_p_kd)
mixscape_class.global <- as.data.frame(T2K2_test$mixscape_class.global)


mixscape_class.global2 = vector()
j = 1
for (i in p_kd){
  if (mixscape_class.global[j,1] == "NT"){
    x = "NT"
    j = j+1
  }
  
  else if (i > 0.9){
    x = "KD"
    j = j+1
  }
  
  else {
    x = "NP"
    j = j+1
  }
  mixscape_class.global2 = c(mixscape_class.global2,x)
}

T2K2_test$mixscape_class.global2 = mixscape_class.global2

#
mix_class <- as.data.frame(as.vector(T2K2_test$mixscape_class.global2))
mixscape_class2 = vector()
gene_class <- as.data.frame(as.vector(T2K2_test$perturb))

for (i in 1:nrow(mix_class)){
  if (gene_class[i,1] == "NT"){
    x = "NT"}
  
  else {x = paste(gene_class[i,1],mix_class[i,1], sep=" ")
  }
  mixscape_class2 = c(mixscape_class2,x)
}

T2K2_test$mixscape_class2 = mixscape_class2

# Calculate percentage of KO cells for all target gene classes.
df <- prop.table(table(T2K2_test$mixscape_class.global2, T2K2_test$perturb),2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "KD"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c("NT", "NP", "KD"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "g")[[1]][1])
df2$guide_number <- sapply(as.character(df2$Var2), 
                           function(x) strsplit(x, split = "g")[[1]][2])
df3 <- df2[-c(which(df2$gene == "NT")),]

p1 <- ggplot(df3, aes(x = "", y = value*100, fill= Var1)) +
  geom_bar(stat= "identity") +
  theme_classic()+
  scale_fill_manual(values = c("grey49", "grey79","coral1")) +
  #scale_fill_manual(values = c("grey79","coral1")) + 
  ylab("% of cells") +
  xlab("sgRNA")

p1 + theme(axis.text.x = element_text(size = 18, hjust = 1), 
           axis.text.y = element_text(size = 18), 
           axis.title = element_text(size = 16), 
           strip.text = element_text(size=8, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 8, scales = "free") +
  labs(fill = "mixscape class") +theme(legend.title = element_text(size = 14),
                                       legend.text = element_text(size = 12)) 





##
T2K2 = T2K2_test
T2K2 <- SetIdent(T2K2, value = "Default")
DefaultAssay(object = T2K2) <- "peaks"
table(Idents(T2K2))



LINC00511 <- WhichCells(object = T2K2, cells = c(
  "T2_AAACCGAAGGACCGCT-1","T2_ACTAACTCATTATGCG-1","T2_ATATAGGCAGGACACA-1","T2_CAAACGCGTCAGGCCA-1","T2_CGGCAATGTTGGTTCT-1",
  "T2_GTAGTTATCAGGATGA-1","T2_TCAGGTTAGTATCGCG-1","T2_TCCTTAGTCGCTAGTG-1","T2_TGATCAGGTGCATTTC-1"
))

T2K2 <- SetIdent(object = T2K2, cells = LINC00511, value = "LINC00511")

markers_NEAT1 <- FindMarkers(object = T2K2, ident.1 = "NEAT1", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                             test.use = "wilcox")
markers_LINC00578 <- FindMarkers(object = T2K2, ident.1 = "0578", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                                 test.use = "wilcox")
markers_LINC01578 <- FindMarkers(object = T2K2, ident.1 = "1578", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                                 test.use = "wilcox")
markers_FAM66C <- FindMarkers(object = T2K2, ident.1 = "FAM66C", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                              test.use = "wilcox")
markers_MIR155HG <- FindMarkers(object = T2K2, ident.1 = "MIR155HG", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                                test.use = "wilcox")
markers_MIR222HG <- FindMarkers(object = T2K2, ident.1 = "MIR222HG", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                                test.use = "wilcox")
markers_SNHG1 <- FindMarkers(object = T2K2, ident.1 = "SNHG1", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                             test.use = "wilcox")
markers_PURPL <- FindMarkers(object = T2K2, ident.1 = "PURPL", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                             test.use = "wilcox")
#markers_BRD4 <- FindMarkers(object = T2K2, ident.1 = "BRD4", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                            test.use = "wilcox")

markers_LINC01638 <- FindMarkers(object = T2K2, ident.1 = "1638", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                                 test.use = "wilcox")
markers_XIST <- FindMarkers(object = T2K2, assay = "RNA", ident.1 = "XIST",ident.2 = "NT", test.use = "wilcox",
                            min.pct = 0.1,logfc.threshold = 0.25)

markers_LINC00511 <- FindMarkers(object = T2K2, assay = "RNA", ident.1 = "LINC00511",ident.2 = "NT", test.use = "wilcox",
                                 min.pct = 0.1,logfc.threshold = 0.25)


####
markers_MEG3 <- FindMarkers(object = T2K2, ident.1 = "MEG3", ident.2 = "NT", only.pos = FALSE, logfc.threshold = 0.25, 
                            test.use = "wilcox")

markers_GAS5 <- FindMarkers(object = T2K2, assay = "peaks", ident.1 = "GAS5",ident.2 = "NT", test.use = "wilcox",
                            min.pct = 0.1,logfc.threshold = 0.25)

markers_HOTAIRM1 <- FindMarkers(object = T2K2, assay = "peaks", ident.1 = "HOTAIRM1",ident.2 = "NT", test.use = "wilcox",
                                min.pct = 0.1,logfc.threshold = 0.25)


markers_LINC02693<- FindMarkers(object = T2K2, assay = "peaks", ident.1 = "2693",ident.2 = "NT", test.use = "wilcox",
                                min.pct = 0.1,logfc.threshold = 0.25)

######
fc.name = "avg_log2FC"
pos.markers_NEAT1 <- markers_NEAT1[which(x = markers_NEAT1[,fc.name] > 0.25), ]
neg.markers_NEAT1 <- markers_NEAT1[which(x = markers_NEAT1[,fc.name] < -0.25), ]

pos.markers_LINC00578 <- markers_LINC00578[which(x = markers_LINC00578[,fc.name] > 0.25), ]
neg.markers_LINC00578 <- markers_LINC00578[which(x = markers_LINC00578[,fc.name] < -0.25), ]

pos.markers_LINC01578 <- markers_LINC01578[which(x = markers_LINC01578[,fc.name] > 0.25), ]
neg.markers_LINC01578 <- markers_LINC01578[which(x = markers_LINC01578[,fc.name] < -0.25), ]

pos.markers_FAM66C  <- markers_FAM66C[which(x = markers_FAM66C[,fc.name] > 0.25), ]
neg.markers_FAM66C  <- markers_FAM66C[which(x = markers_FAM66C[,fc.name] < -0.25), ]

pos.markers_MIR155HG  <- markers_MIR155HG[which(x = markers_MIR155HG[,fc.name] > 0.25), ]
neg.markers_MIR155HG  <- markers_MIR155HG[which(x = markers_MIR155HG[,fc.name] < -0.25), ]

pos.markers_MIR222HG  <- markers_MIR222HG[which(x = markers_MIR222HG[,fc.name] > 0.25), ]
neg.markers_MIR222HG  <- markers_MIR222HG[which(x = markers_MIR222HG[,fc.name] < -0.25), ]

pos.markers_SNHG1 <- markers_SNHG1[which(x = markers_SNHG1[,fc.name] > 0.25), ]
neg.markers_SNHG1  <- markers_SNHG1[which(x = markers_SNHG1[,fc.name] < -0.25), ]

pos.markers_PURPL <- markers_PURPL[which(x = markers_PURPL[,fc.name] > 0.25), ]
neg.markers_PURPL <- markers_PURPL[which(x = markers_PURPL[,fc.name] < -0.25), ]

#pos.markers_BRD4  <- markers_BRD4[which(x = markers_BRD4[,fc.name] > 0.25), ]
#neg.markers_BRD4  <- markers_BRD4[which(x = markers_BRD4[,fc.name] < -0.25), ]

pos.markers_LINC01638  <- markers_LINC01638[which(x = markers_LINC01638[,fc.name] > 0.25), ]
neg.markers_LINC01638  <- markers_LINC01638[which(x = markers_LINC01638[,fc.name] < -0.25), ]

pos.markers_MEG3  <- markers_MEG3[which(x = markers_MEG3[,fc.name] > 0.25), ]
neg.markers_MEG3  <- markers_MEG3[which(x = markers_MEG3[,fc.name] < -0.25), ]

pos.markers_GAS5  <- markers_GAS5[which(x = markers_GAS5[,fc.name] > 0.25), ]
neg.markers_GAS5  <- markers_GAS5[which(x = markers_GAS5[,fc.name] < -0.25), ]

pos.markers_HOTAIRM1 <- markers_HOTAIRM1[which(x = markers_HOTAIRM1[,fc.name] > 0.25), ]
neg.markers_HOTAIRM1 <- markers_HOTAIRM1[which(x = markers_HOTAIRM1[,fc.name] < -0.25), ]

pos.markers_LINC02693  <- markers_LINC02693[which(x = markers_LINC02693[,fc.name] > 0.25), ]
neg.markers_LINC02693  <- markers_LINC02693[which(x = markers_LINC02693[,fc.name] < -0.25), ]

pos.markers_XIST  <- markers_XIST[which(x = markers_XIST[,fc.name] > 0.25), ]
neg.markers_XIST  <- markers_XIST[which(x = markers_XIST[,fc.name] < -0.25), ]

pos.markers_LINC00511  <- markers_LINC00511[which(x = markers_LINC00511[,fc.name] > 0.25), ]
neg.markers_LINC00511  <- markers_LINC00511[which(x = markers_LINC00511[,fc.name] < -0.25), ]

pval.cutoff = 0.05
max.genes = 30
marker.list <- c(rownames(x = subset(x = pos.markers_NEAT1,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_NEAT1,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_PURPL,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_PURPL,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_MIR222HG,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_MIR222HG,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_MIR155HG,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_MIR155HG,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_FAM66C,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_FAM66C,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_LINC01578,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_LINC01578,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_LINC00578,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_LINC00578,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_SNHG1,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_SNHG1,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_LINC01638,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_LINC01638,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_MEG3,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_MEG3,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_GAS5,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_GAS5,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_HOTAIRM1,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_HOTAIRM1,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_LINC02693,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_LINC02693,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_XIST,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_XIST,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = pos.markers_LINC00511,p_val < pval.cutoff))[1:max.genes],
                 rownames(x = subset(x = neg.markers_LINC00511,p_val < pval.cutoff))[1:max.genes]
                 #rownames(x = subset(x = pos.markers_BRD4,p_val < pval.cutoff))[1:max.genes],
                 #rownames(x = subset(x = neg.markers_BRD4,p_val < pval.cutoff))[1:max.genes]
)



##
max.cells.group = 20
sub <- subset(x = T2K2, idents = c("NT03","NEAT1","0578","1578","FAM66C",
                                        "MIR155HG","MIR222HG","SNHG1","PURPL",
                                        #"BRD4_KD",
                                        "1638","MEG3","GAS5","HOTAIRM1","2693",
                                        "XIST","0511"), downsample = max.cells.group)

peak_features <- rownames(sub[["peaks"]])
# 确保marker.list中的特征与peaks assay中的顺序一致
marker.list <- marker.list[marker.list %in% peak_features]
marker.list <- marker.list[order(match(marker.list, peak_features))]
sub <- ScaleData(object = sub, features = marker.list,assay = "peaks")

###peakheatmap
#isTRUE(x = order.by.prob)) 
prtb.type = "KD"
p_ko <- sub[[paste0("mixscape_class", "_p_", tolower(x = prtb.type))]][, 
                                                                        1, drop = FALSE]
ordered.cells <- rownames(x = p_ko)[order(p_ko[,1], decreasing = TRUE)]

#####
pp = as.data.frame(Idents(sub))
pp2 =vector()
pp2 = pp$`Idents(sub)`
sub$level = pp2

cluster_info <- sub$level
cluster_info <- factor(cluster_info, levels = c("NEAT1", "PURPL",  "MIR222HG", "MIR155HG", "FAM66C", "1578",
                                                "0578", "SNHG1",
                                                "1638","MEG3","GAS5","HOTAIRM1","2693",
                                                "XIST","0511",
                                                #"BRD4_KD",
                                                "NT"))
gene_features <- marker.list

sub2 <- GetAssayData(sub, slot = "scale.data",assay = 'peaks')
#sub3 = cbind(rownames(sub3),sub3)
sub2 <- as.matrix(sub2[marker.list, names(cluster_info)])
#gene_pos <- which(rownames(sub2)%in%gene)
#row_anno <- rowAnnotation(gene = anno_mark(at=gene_pos,labels = gene),annotation_name_gp= gpar(fontsize = 4))
#
library(RColorBrewer)

col <- c("#9E0142","#8C510A","#5E4FA2","#01665E","#542788","#4D9221","#8E0152","#225EA8",
         "#C2A5CF","#D0D1E6","#67A9CF","#DFC27D","#74A9CF","#BF812D","#8C96C6",
         "#D4B9DA","#BDBDBD")
names(col) <- levels(cluster_info)

library(circlize)
col_fun = colorRamp2(c(-4,-2,0, 2,4), c("#053061","#4393C3","#F7F7F7","#D6604D","#67001F"))
col_fun(seq(-3, 3,0.5),)
lgd = Legend(col_fun = col_fun, title = "scaled mRNA",at = c(-4,4),
             labels = c("low", "high"),border="black",legend_height = unit(3, "cm"),legend_width = unit(2, "mm") ,
             labels_gp = gpar(col = "black", font = 2),direction = "horizontal",title_position = "topcenter")

df = data.frame(cluster_info)

top_anno <- HeatmapAnnotation(df = df,
                              annotation_height = unit(5, "mm"),show_legend = T, #annotation_label = "",
                              col = list(cluster_info = c("NEAT1" = "#9E0142","PURPL" ="#8C510A","MIR222HG" ="#5E4FA2",
                                                          "MIR155HG" = "#01665E",  "FAM66C" = "#542788","1578" ="#4D9221",
                                                          "0578" = "#8E0152", "SNHG1" = "#225EA8",
                                                          
                                                          "1638" = "#C2A5CF","MEG3" = "#D0D1E6","GAS5" = "#67A9CF","HOTAIRM1" = "#DFC27D",
                                                          "2693" = "#74A9CF","XIST" = "#BF812D","0511" = "#8C96C6",
                                                          "NT" = "#BDBDBD"))#,height = 1, width = 1, gap = unit(1, "mm")
)                             
p <- Heatmap(sub2,
             cluster_rows = F,
             cluster_columns = FALSE,
             show_column_names = FALSE,
             show_row_names = F,
             column_split = cluster_info,
             top_annotation = top_anno,
             #right_annotation  = row_anno,
             column_title = NULL,
             column_order = ordered.cells,
             col = col_fun,
             border = TRUE,
             heatmap_legend_param = list(col_fun = col_fun, title = "Scaled mRNA",title_position = "leftcenter-rot",at = c(-4,4),
                                         labels = c("low", "high"),border="black",legend_height = unit(20, "mm"),legend_width = unit(1, "mm"),
                                         labels_gp = gpar(col = "black", font = 2),#direction = "horizontal",
                                         title_position = "topcenter")
)

p






