if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(Seurat)

library(tidyverse)
library(Signac)
library(hdf5r)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biovizBase)
library(patchwork)
library(GenomeInfoDb)
library(dplyr)


dir = "filtered_feature_bc_matrix"
list.files(dir)

counts <- Read10X(data.dir = dir)
class(counts)

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

fragment.path <- 'atac_fragments.tsv.gz'

scmulti <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

scmulti[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragment.path,
  annotation = annotation
)

DefaultAssay(scmulti) <- "ATAC"

scmulti <- NucleosomeSignal(scmulti)
scmulti <- TSSEnrichment(scmulti)

# call peaks using MACS2
peaks <- CallPeaks(scmulti, macs2.path = "/Users/ramzipit/opt/anaconda3/bin/MACS2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)


# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(scmulti),
  features = peaks,
  cells = colnames(scmulti)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
scmulti[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragment.path,
  annotation = annotation
)

##
DefaultAssay(scmulti) = 'RNA'

scmulti <- NormalizeData(scmulti, normalization.method = "LogNormalize", scale.factor = 10000)


# cellcycle
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match = rownames(scmulti))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search=s_genes, match = rownames(scmulti))
scmulti <- CellCycleScoring(scmulti,g2m.features=g2m_genes,s.features=s_genes)

table(scmulti$Phase)

scmulti <- ScaleData(scmulti) 
scmulti <- FindVariableFeatures(scmulti, selection.method = "vst", nfeatures = 2000) 
scmulti<- RunPCA(object = scmulti, pc.genes = VariableFeatures(scmulti)) 
scmulti <- FindNeighbors(scmulti, dims = 1:15)

#scmulti <- FindClusters(scmulti, resolution = 0.2)
#table(scmulti@meta.data$RNA_snn_res.0.2) 
#set.seed(123)
scmulti <- RunTSNE(object = scmulti, dims = 1:15, do.fast = TRUE)
DimPlot(scmulti,reduction = "tsne",label=T)

scmulti  <- RunUMAP(scmulti, dims = 1:15)
DimPlot(scmulti, reduction = "umap")+ggtitle("UMAP")


########scATAC 
DefaultAssay(scmulti) <- "peaks"
scmulti <- FindTopFeatures(scmulti, min.cutoff = 5)
scmulti <- RunTFIDF(scmulti)
scmulti <- RunSVD(scmulti)


#使用DepthCor函数评估每个LSI成分与测序深度之间的相关性：
DepthCor(scmulti) #我们可以看到第一个LSI成分与细胞中counts的总数之间有非常强的相关性，因此我们将在后续的分析中删除此成分。

# 使用RunUMAP函数进行UMAP非线性降维
scmulti <- RunUMAP(object = scmulti, reduction = 'lsi', dims = 2:30)
# 对细胞执行基于图的聚类
scmulti <- FindNeighbors(object = scmulti, reduction = 'lsi', dims = 2:30)
scmulti <- FindClusters(object = scmulti, verbose = FALSE, algorithm = 3)
# 使用DimPlot函数进行数据可视化
DimPlot(object = scmulti, label = TRUE) + NoLegend()


####
GBC_sig
GBC_sig = GBC_sig[GBC_sig$identity != 'blank',]
x = as.data.frame(rownames(GBC_sig))

##
DefaultAssay(scmulti) <- "RNA"

ID = vector()

for (i in 1:nrow(x)){
  k = paste(x[i,1],"1",sep="-")
  ID = c(ID,k)
}

#
scmulti <- scmulti[,WhichCells(object = scmulti, cells = ID, invert = F)]

IDnew = as.data.frame(rownames(as.data.frame(scmulti@active.ident)))
colnames(IDnew)[1] = 'ID'

GBC_sig = cbind(GBC_sig,ID)

GBC_sig2 = left_join(IDnew,GBC_sig,by = 'ID')

#####add GBC info
GBC_signature = vector()
GBC_signature = paste(GBC_sig2$Call_1,GBC_sig2$Call_2,GBC_sig2$Call_3,sep='+')

#scmulti_K22$GBC_signature = GBC_signature 
GBC_sig2 = cbind(GBC_sig2[,2:5],GBC_sig2[1])

NT_freq = table(GBC_signature)
NT_freq[names(NT_freq) == "nontarget-1+NA+NA"]
NT_freq[names(NT_freq) == "nontarget-2+NA+NA"]

GBC_signature2 = as.matrix(GBC_signature)
xx = vector()
for (i in 1:nrow(GBC_sig2)){
  if ((GBC_sig2[i,1] == "nontarget-1" & GBC_sig2[i,4] == 'single') | 
      (GBC_sig2[i,1] == "nontarget-2" & GBC_sig2[i,4] == 'single') | 
      (GBC_sig2[i,1] == "nontarget-1" & (GBC_sig2[i,2] == "nontarget-1" | GBC_sig2[i,2] == "nontarget-2") & GBC_sig2[i,4] == 'double') |
      (GBC_sig2[i,1] == "nontarget-2" & (GBC_sig2[i,2] == "nontarget-1" | GBC_sig2[i,2] == "nontarget-2") & GBC_sig2[i,4] == 'double')){
    x = "NT"}
  else if(GBC_sig2[i,4] == "single"){
    x = GBC_sig2[i,1]}
  else if(GBC_sig2[i,4] == "double"){
    x = paste(GBC_sig2[i,1],GBC_sig2[i,2],sep = '+')
  }
  else if(GBC_sig2[i,4] == "triple"){
    x = paste(GBC_sig2[i,1],GBC_sig2[i,2],GBC_sig2[i,3],sep = '+')
  }
  else {
    x = 'NA'
  }
  
  xx = c(xx,x)
}
GBC_sig2 =cbind(GBC_sig2,xx)

scmulti$perturb = GBC_sig2$xx


oneGBClist_0.05 <- as.data.frame(scmulti$perturb)
oneGBClist_0.05 <- cbind(rownames(oneGBClist_0.05),oneGBClist_0.05)
colnames(oneGBClist_0.05)[1] = 'ID'

oneGBClist_0.05 = left_join(oneGBClist_0.05,GBC_sig2,by = 'ID') 

scmulti$perturb_number = oneGBClist_0.05$identity

x <- as.data.frame(scmulti$perturb)
x <- cbind(rownames(x),x)
colnames(x)[1] = 'ID'
x = left_join(x,GBC_sig2,by = 'ID') 

y = x$xx
z<- ifelse(y == 'NT', y, 'perturb')
oneGBClist_0.05$cripsr = z

scmulti$crispr = oneGBClist_0.05$cripsr

#GBC_BJ_raw = cbind(rownames(GBC_BJ_raw),GBC_BJ_raw)
oneGBClist_0.05_single =oneGBClist_0.05[oneGBClist_0.05$identity == 'single',]
table(oneGBClist_0.05_single$Call_2)

#scmulti <- subset(x = scmulti,subset = nCount_ATAC < 700000 & nCount_RNA < 130000 & nCount_ATAC > 80 & nCount_RNA > 600 & nucleosome_signal < 18 &TSS.enrichment > 1)


scmulti_single <- subset(scmulti, perturb== "NT")

scmulti_single <- subset(scmulti, perturb == "7SK" | perturb == "BRD4" | perturb == "LINC00189"| perturb == "LINC00511"| perturb == "LINC00578" | perturb == "LINC00657" | perturb == "LINC01011" 
                            | perturb == "LINC01012" | perturb == "LINC01133" | perturb == "LINC01578" | perturb == "LINC02693" | perturb == "CASC15" | perturb == "FAM66c" | perturb == "FENDR" | perturb == "GAS5"| perturb == "NEAT1"
                            | perturb == "HOTAIRM1"| perturb == "JPX"| perturb == "LRRC75a"| perturb == "MIAT"| perturb == "MALAT1"| perturb == "MEG3"| perturb == "MIR222HG"| perturb == "MIR155HG"| perturb == "PURPL"
                            | perturb == "SNHG3"| perturb == "RMRP"| perturb == "SNHG1"| perturb == "XIST" | perturb == "SNHG16"| perturb == "LINC01615"| perturb == "LINC01638" |
                              perturb == "NT")

## CalcPerturbSig-Best:ndims:15; KNN:10
scmulti_single <- CalcPerturbSig(object = scmulti_single,
                                    assay = "RNA",
                                    slot = "data",
                                    gd.class ="perturb",
                                    nt.cell.class = "NT",
                                    reduction = "pca",
                                    #split.by = "replicate",
                                    ndims = 15,
                                    num.neighbors = 10,
                                    new.assay.name = "PRTB",
)

# Prepare PRTB assay for dimensionality reduction: 
# Normalize data, find variable features and center data.


# Use variable features from RNA assay.
VariableFeatures(object = scmulti_single) <- VariableFeatures(object = scmulti_single[["RNA"]])
#scmulti3 <- ScaleData(object = scmulti3, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
scmulti_single <- RunPCA(object = scmulti_single, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
scmulti_single <- RunUMAP(
  object = scmulti_single, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')


# Generate plots to check if clustering is driven by biological replicate ID, 
# cell cycle phase or target gene class.
q2 <- DimPlot(
  object = scmulti_single, 
  group.by = 'Phase', 
  reduction = 'prtbumap', 
  pt.size = 0.2, label = F, repel = T) +
  ggtitle("Cell Cycle Phase") +
  ylab("UMAP 2") +
  xlab("UMAP 1") 

q2

q3 <- DimPlot(
  object = scmulti_single,
  group.by = 'crispr',
  reduction = 'prtbumap', 
  split.by = "crispr", 
  ncol = 1, 
  pt.size = 0.2, 
  cols = c("grey39","goldenrod3")) +
  ggtitle("Perturbation Status") +
  ylab("UMAP 2") +
  xlab("UMAP 1") 

q3

# Run mixscape 

PerturbDiff <- function(object, assay, slot, all_cells, nt_cells, features, neighbors, verbose) {
  nt_data <- as.matrix(x = expm1(x = GetAssayData(object = object, assay = assay, slot = slot)[features, nt_cells, drop = FALSE]))
  if (verbose) {
    mysapply <- pbsapply
  } else {
    mysapply <- sapply
  }
  
  # new_expr <- mysapply(X = all_cells, FUN = function(i) {
  #   index <- Indices(object = neighbors)[i, ]
  #   nt_cells20 <- nt_cells[index]
  #   avg_nt <- rowMeans2(x = nt_data[, nt_cells20, drop = FALSE])
  #   avg_nt <- as.matrix(x = avg_nt)
  #   colnames(x = avg_nt) <- i
  #   return(avg_nt)
  # })
  idx <- Indices(object = neighbors)[all_cells,]
  model.matrix <- sparseMatrix(i = as.vector(idx), j = rep(1:nrow(x = idx), times = ncol(x = idx)), x = 1, dims = c(length(x = nt_cells), nrow(x = idx)))
  model.matrix <- model.matrix/rep(colSums(model.matrix), each = nrow(x = model.matrix))
  new_expr <- nt_data %*% model.matrix
  
  new_expr <- matrix(data = new_expr, nrow = length(x = features))
  new_expr <- log1p(x = new_expr)
  rownames(x = new_expr) <- rownames(x = nt_data)
  colnames(x = new_expr) <- all_cells
  diff <- new_expr - as.matrix(GetAssayData(object = object, slot = slot, assay = assay)[features, colnames(x = new_expr), drop = FALSE])
  return(diff)
}

assay = pc.assay
features = missing_genes
new.assay.name = "PRTB_inte"
GetMissingPerturb <- function(object, assay, features, verbose = TRUE) {
  if (length(x = features) == 0) {
    return(object)
  }
  if (verbose) {
    message("Computing perturbation signature for missing features.")
  }
  command <- grep(pattern = "CalcPerturbSig", x = Command(object = object), value = TRUE)
  command.match <- sapply(X = command, FUN = function(x) {
    Command(object = object, command = x, value = "new.assay.name") == assay
  })
  if (length(x = which(x = command.match)) > 1) {
    stop("Ambiguous command log.")
  }
  if(length(x = which(x = command.match)) == 0) {
    stop("Cannot find previously run CalcPertubSig command. Please make sure you've run CalcPerturbSig to create the provided assay.")
  }
  command <- names(x = command.match)
  command <- "CalcPerturbSig.RNA.pca"  
  if ("split.by" %in% names(x = slot(object = Command(object = object, command = command), name ="params"))) {
    split.by <- Command(object = object, command = command, value = "split.by")
  } else {
    split.by <- NULL
  }
  gd.class <- Command(object = object, command = command, value = "gd.class")
  nt.cell.class <- Command(object = object, command = command, value = "nt.cell.class")
  slot <- Command(object = object, command = command, value = "slot")
  assay.orig <- Command(object = object, command = command, value = "assay")
  old.idents <- Idents(object = object)
  if (! is.null(x = split.by)) {
    Idents(object = object) <-  split.by
  } else {
    Idents(object = object) <- "rep1"
  }
  replicate <- unique(x = Idents(object = object))
  all_diff <- list()
  all_nt_cells <- Cells(x = object)[which(x = object[[]][gd.class] == nt.cell.class)]
  features <- setdiff(x = features, y = rownames(x = object[[assay]]))
  for (r in replicate) {
    # isolate nt cells
    all_cells <- WhichCells(object = object, idents = r)
    nt_cells <- intersect(x = all_nt_cells, all_cells)
    # pull previously computed neighbors
    neighbors <- Tool(object = object, slot = command)[[make.names(names = paste0(assay, "_", r))]]
    diff <- PerturbDiff(
      object = object,
      assay = assay.orig,
      slot = slot,
      all_cells = all_cells,
      nt_cells = nt_cells,
      features = features,
      neighbors = neighbors,
      verbose = verbose
    )
    all_diff[[r]] <- diff
  }
  all_diff <- do.call(what = cbind, args = all_diff)
  all_diff <- all_diff[, colnames(x = object[[assay]]), drop = FALSE]
  new.assay <- CreateAssayObject(
    data = rbind(
      GetAssayData(object = object[[assay]], slot = "data"),
      all_diff
    ),
    min.cells = 0,
    min.features = 0,
    check.matrix = FALSE
  )
  new.assay <- SetAssayData(
    object = new.assay,
    slot = "scale.data",
    new.data = GetAssayData(object = object[[assay]], slot = "scale.data")
  )
  object[[assay]] <- new.assay
  Idents(object = object) <- old.idents
  return(object)
}



#####p_val_TopDEGenesMixscape
TopDEGenesMixscape <- function(
    object,
    ident.1,
    ident.2 = NULL,
    labels = 'gene',
    de.assay = "RNA",
    test.use = "wilcox",
    pval.cutoff = 5e-2,
    logfc.threshold = 0.25,
    verbose = TRUE
) {
  if (verbose) {
    message("Finding new perturbation gene set")
  }
  de.genes <- data.frame()
  tryCatch(
    expr = {
      de.genes <- FindMarkers(
        object = object,
        ident.1 = ident.1,
        ident.2 = ident.2,
        group.by = labels,
        assay = de.assay,
        test.use = test.use,
        logfc.threshold = logfc.threshold,
        verbose = verbose,
        min.pct = 0.1
      )
      de.genes <- de.genes[de.genes$p_val < pval.cutoff, ]
    },
    error = function(e) {}
  )
  return(rownames(x = de.genes))
}

ProjectVec <- function(v1, v2) {
  return(as.vector(x = (v1 %*% v2) / (v2 %*% v2)))
}

DefineNormalMixscape <- function(x) {
  mu <- mean(x)
  sd <- sd(x)
  return(list(mu = mu, sd = sd))
}


library(tidyverse)
library(Signac)
library(Seurat)
library(hdf5r)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biovizBase)
library(patchwork)
library(GenomeInfoDb)
library(dplyr)
library(matrixStats)
library(Matrix)

scmulti_single_test = scmulti_single

#p_val_adj_TopDEGenesMixscape
TopDEGenesMixscape <- function(
    object,
    ident.1,
    ident.2 = NULL,
    labels = 'gene',
    de.assay = "RNA",
    test.use = "wilcox",
    pval.cutoff = 5e-2,
    logfc.threshold = 0.25,
    verbose = TRUE
) {
  if (verbose) {
    message("Finding new perturbation gene set")
  }
  de.genes <- data.frame()
  tryCatch(
    expr = {
      de.genes <- FindMarkers(
        object = object,
        ident.1 = ident.1,
        ident.2 = ident.2,
        group.by = labels,
        assay = de.assay,
        test.use = test.use,
        logfc.threshold = logfc.threshold,
        verbose = verbose,
        min.pct = 0.1
      )
      de.genes <- de.genes[de.genes$p_val_adj < pval.cutoff, ]
    },
    error = function(e) {}
  )
  return(rownames(x = de.genes))
}


#####
scmulti_single_test = scmulti_single
object = scmulti_single_test
assay = "PRTB"
slot = "scale.data"
labels = "perturb"
nt.class.name = "NT"
new.class.name = "mixscape_class"
min.de.genes = 40
min.cells = 5
de.assay = "RNA"
logfc.threshold = 0.15
iter.num = 10
verbose = FALSE
split.by = NULL
fine.mode = FALSE
fine.mode.labels = "guide_ID"
prtb.type = "KD"

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

scmulti_single_test = object

p_kd <- as.vector(scmulti_single_test$mixscape_class_p_kd)
mixscape_class.global <- as.data.frame(scmulti_single_test$mixscape_class.global)


mixscape_class.global2 = vector()
j = 1
for (i in p_kd){
  if (mixscape_class.global[j,1] == "NT"){
    x = "NT"
    j = j+1
  }
  
  else if (i > 0.6){
    x = "KD"
    j = j+1
  }
  
  else {
    x = "NP"
    j = j+1
  }
  mixscape_class.global2 = c(mixscape_class.global2,x)
}

scmulti_single_test$mixscape_class.global2 = mixscape_class.global2

#
mix_class <- as.data.frame(as.vector(scmulti_single_test$mixscape_class.global2))
mixscape_class2 = vector()
gene_class <- as.data.frame(as.vector(scmulti_single_test$perturb))

for (i in 1:nrow(mix_class)){
  if (gene_class[i,1] == "NT"){
    x = "NT"}
  
  else {x = paste(gene_class[i,1],mix_class[i,1], sep=" ")
  }
  mixscape_class2 = c(mixscape_class2,x)
}

scmulti_single_test$mixscape_class2 = mixscape_class2




# Calculate percentage of KD cells for all target gene classes.
df <- prop.table(table(scmulti_single_test$mixscape_class.global2, scmulti_single_test$perturb),2)

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
T2_K2_integrated
object = T2_K2_integrated 

assay = 'RNA'
de.assay = "RNA"
reduction.key = "LDA_"
seed = 42
pc.assay = "PRTB_inte"
labels = "perturb"
nt.label = "NT"
npcs = 30
verbose = TRUE
logfc.threshold = 0.25

library(pbapply)


      projected_pcs <- list()
      gene_list <- setdiff(x = unique(x = object[[labels]][, 1]), y = nt.label)
      Idents(object = object) <- labels
      DefaultAssay(object = object) <- pc.assay
      all_genes <- list()
      nt.cells <- WhichCells(object = object, idents = nt.label)
      for (g in gene_list) {
        if (verbose) {
          message(g)
        }
        gd.cells <- WhichCells(object = object, idents = g)
        gene_set <- TopDEGenesMixscape(
          object = object,
          ident.1 = gd.cells,
          ident.2 = nt.cells,
          de.assay = de.assay,
          logfc.threshold = logfc.threshold,
          labels = labels,
          verbose = verbose
        )
        if (length(x = gene_set) < (npcs + 1)) {
          all_genes[[g]] <- character()
          next
        }
        all_genes[[g]] <- gene_set
      }
      all_markers <- unique(x = unlist(x = all_genes))
      missing_genes <- all_markers[!all_markers %in% rownames(x = object[[pc.assay]])]
      object <- GetMissingPerturb(object = object, assay = pc.assay, features = missing_genes, verbose = verbose)
      for (g in gene_list) {
        if (verbose) {
          message(g)
        }
        gene_subset <- subset(x = object, idents = c(g, nt.label))
        gene_set <- all_genes[[g]]
        if (length(x = gene_set) == 0) {
          next
        }
        gene_subset <- ScaleData(
          object = gene_subset,
          features = gene_set,
          verbose = FALSE
        )
        gene_subset <- RunPCA(
          object = gene_subset,
          features = gene_set,
          npcs = npcs,
          verbose = FALSE
        )
        project_pca <- ProjectCellEmbeddings(
          reference = gene_subset,
          query = object,
          dims = 1:npcs,
          verbose = FALSE
        )
        colnames(x = project_pca) <- paste(g, colnames(x = project_pca), sep = "_")
        projected_pcs[[g]] <- project_pca
      }
      return(projected_pcs)
    

projected_pcs

lda.lables <- object[[labels]][,]

object = projected_pcs
labels = lda.lables
assay = "RNA"

verbose = TRUE
ndims.print = 1:5
nfeatures.print = 30
reduction.key = "LDA_"
seed = 42

library(MASS)


      if (!is.null(x = seed)) {
        set.seed(seed = seed)
      }
      object <- data.frame(object)
      var_names <- colnames(x = object)
      object$lda_cluster_label <- labels
      lda_results <- lda(formula = lda_cluster_label ~ ., data = object)
      lda_predictions <- predict(object = lda_results, newdata = object)
      
      lda_cv <-lda(
        formula = lda_cluster_label ~ .,
        data = object,
        CV = TRUE
      )$posterior
      feature.loadings <- lda_results$scaling
      cell.embeddings <- lda_predictions$x
      lda.assignments <- lda_predictions$class
      lda.posterior <- lda_predictions$posterior
      colnames(x = lda.posterior) <- paste0("LDAP_", colnames(x = lda.posterior))
      rownames(x = feature.loadings) <- var_names
      colnames(x = feature.loadings) <- paste0(reduction.key, 1:ncol(x = cell.embeddings))
      rownames(x = cell.embeddings) <- rownames(x = object)
      colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
      reduction.data <- CreateDimReducObject(
        embeddings = cell.embeddings,
        loadings = feature.loadings,
        assay = assay,
        key = reduction.key,
        misc = list(
          assignments = lda.assignments,
          posterior = lda.posterior,
          model = lda_results,
          cv = lda_cv
        )
      )
      if (verbose) {
        print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
      }
      #return(reduction.data)

object_lda = reduction.data

T2_K2_integrated[["lda"]] <- object_lda


T2_K2_integrated <- RunUMAP(
  object = T2_K2_integrated,
  dims = 1:15,
  reduction = 'lda',
  reduction.key = 'ldaumap',
  reduction.name = 'ldaumap')



# 
all_cells <- Cells(T2_K2_integrated)

# 找出需要保留的细胞
cells_to_keep <- setdiff(all_cells, cells_to_remove)

T2_K2_integrated <- subset(T2_K2_integrated, cells = cells_to_keep)

# 
perturb_values <- as.data.frame(T2_K2_integrated@meta.data$perturb)
ident = as.data.frame(T2_K2_integrated@active.ident)
rownames(perturb_values) = rownames(ident)
perturb_values = cbind(rownames(perturb_values),perturb_values)
colnames(perturb_values)[1] = 'id'

x = perturb_values[!( perturb_values$`T2_K2_integrated3@meta.data$perturb` %in% c(
  "NT" ,"LINC00511" , "NEAT1" ,"XIST","FAM66c" ,"MIR155HG",
  "MIR222HG","LINC01638","GAS5" ,"LINC00578" ,"HOTAIRM1",
  "LINC02693","LINC01578",'PURPL',"MEG3", "SNHG1")),]

cells_to_remove = x$id
all_cells <- Cells(T2_K2_integrated)
cells_to_keep <- setdiff(all_cells, cells_to_remove)
T2_K2_integrated <- subset(T2_K2_integrated, cells = cells_to_keep)

cell_type_cols <- c("NT" = "#AFB1B7","LINC00511" = "#BA4A94", 
                    "NEAT1" =  "#4AA123","XIST" = "#CB1D1B",
                    "FAM66c" = "#7573AD","MIR155HG" = "#1927aa",
                    "MIR222HG" =  "#DDE9B0","LINC01638" = "#72bcd9",
                    "GAS5" = "#F4C685",
                    
                    "LINC00578" = "#F7BEAE","HOTAIRM1" = "#e56e24",
                    "LINC02693" = "#7BAE92",
                    "LINC01578" = '#2B543A',
                    'PURPL'="#812F33","MEG3" =  "#C9A47D", "SNHG1" = '#EAAFCD'
                   # "LINC01012" = "#7A8EBA", "LRRC75a" = "#C7BCB1","MIAT" = "#6C4C73",
                    #"RMRP" = "#b4446c","FAM66C" ="#EAD9E1" ,"SNHG3" = "#4D4C72"
                    )

Idents(T2_K2_integrated) <- "perturb"
plot7 <- DimPlot(T2_K2_integrated, label = F, 
                 reduction = "ldaumap", repel = T, 
                 pt.size = 0.2,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot7


##
T2_K2_integrated3 <- T2_K2_integrated3[c("K2", "T2"), ]

cell_type_cols <- c("K2" = "#3081BC",
                    "T2" =  "#AC4A4B"
                    )


Idents(T2_K2_integrated3) <- "batch"
plot8 <- DimPlot(T2_K2_integrated3, label = F, 
                 reduction = "ldaumap", repel = T, 
                 pt.size = 0.06,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot8 


####ATAC

# Calculate perturbation signature (PRTB).
T2_K2_integrated<- CalcPerturbSig(
  object = T2_K2_integrated, 
  assay = "peaks", 
  slot = "data", 
  gd.class ="perturb", 
  nt.cell.class = "NT", 
  reduction = "pca", 
  ndims = 15, 
  num.neighbors = 10, 
  split.by = "batch", 
  new.assay.name = "PRTB_inte_peaks")

# Prepare PRTB assay for dimensionality reduction: 
# Normalize data, find variable features and center data.
DefaultAssay(object = T2_K2_integrated) <- 'PRTB_inte_peaks'

# Use variable features from RNA assay.
#T2_K2_integrated <- FindVariableFeatures(T2_K2_integrated, selection.method = "vst", nfeatures = 2000) 


VariableFeatures(object = T2_K2_integrated) <- VariableFeatures(object = T2_K2_integrated[["peaks"]])
T2_K2_integrated <- ScaleData(object = T2_K2_integrated, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
T2_K2_integrated <- RunPCA(object = T2_K2_integrated, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
T2_K2_integrated <- RunUMAP(
  object = T2_K2_integrated, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

# Generate plots to check if clustering is driven by biological replicate ID, 
# cell cycle phase or target gene class.
q1 <- DimPlot(
  object = T2_K2_integrated, 
  group.by = 'batch', 
  reduction = 'prtbumap', 
  pt.size = 0.2, cols = "Dark2", label = F, repel = T) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("Biological Replicate") +
  ylab("UMAP 2") +
  xlab("UMAP 1") 
q1

###
table(T2_K2_integrated$perturb)
object = T2_K2_integrated

assay = 'peaks'
de.assay = "peaks"
reduction.key = "LDA_"
seed = 42
pc.assay = "PRTB_inte_peaks"
labels = "perturb"
nt.label = "NT"
npcs = 30
verbose = TRUE
logfc.threshold = 0.25

library(pbapply)

      projected_pcs <- list()
      gene_list <- setdiff(x = unique(x = object[[labels]][, 1]), y = nt.label)
      Idents(object = object) <- labels
      DefaultAssay(object = object) <- pc.assay
      all_genes <- list()
      nt.cells <- WhichCells(object = object, idents = nt.label)
      for (g in gene_list) {
        if (verbose) {
          message(g)
        }
        gd.cells <- WhichCells(object = object, idents = g)
        gene_set <- TopDEGenesMixscape(
          object = object,
          ident.1 = gd.cells,
          ident.2 = nt.cells,
          de.assay = de.assay,
          logfc.threshold = logfc.threshold,
          labels = labels,
          verbose = verbose
        )
        if (length(x = gene_set) < (npcs + 1)) {
          all_genes[[g]] <- character()
          next
        }
        all_genes[[g]] <- gene_set
      }
      all_markers <- unique(x = unlist(x = all_genes))
      missing_genes <- all_markers[!all_markers %in% rownames(x = object[[pc.assay]])]
      object <- GetMissingPerturb(object = object, assay = pc.assay, features = missing_genes, verbose = verbose)
      Idents(object) <- object$perturb
      for (g in gene_list) {
        if (verbose) {
          message(g)
        }
        gene_subset <- subset(x = object, idents = c(g, nt.label))
        gene_set <- all_genes[[g]]
        if (length(x = gene_set) == 0) {
          next
        }
        gene_subset <- ScaleData(
          object = gene_subset,
          features = gene_set,
          verbose = FALSE
        )
        gene_subset <- RunPCA(
          object = gene_subset,
          features = gene_set,
          npcs = npcs,
          verbose = FALSE
        )
        project_pca <- ProjectCellEmbeddings(
          reference = gene_subset,
          query = object,
          dims = 1:npcs,
          verbose = FALSE
        )
        colnames(x = project_pca) <- paste(g, colnames(x = project_pca), sep = "_")
        projected_pcs[[g]] <- project_pca
      }
      return(projected_pcs)


projected_pcs

lda.lables <- object[[labels]][,]

object = projected_pcs
labels = lda.lables
assay = "peaks"

verbose = TRUE
ndims.print = 1:5
nfeatures.print = 30
reduction.key = "LDA_"
seed = 42

library(MASS)

      if (!is.null(x = seed)) {
        set.seed(seed = seed)
      }
      object <- data.frame(object)
      var_names <- colnames(x = object)
      object$lda_cluster_label <- labels
      lda_results <- lda(formula = lda_cluster_label ~ ., data = object)
      lda_predictions <- predict(object = lda_results, newdata = object)
      
      lda_cv <-lda(
        formula = lda_cluster_label ~ .,
        data = object,
        CV = TRUE
      )$posterior
      feature.loadings <- lda_results$scaling
      cell.embeddings <- lda_predictions$x
      lda.assignments <- lda_predictions$class
      lda.posterior <- lda_predictions$posterior
      colnames(x = lda.posterior) <- paste0("LDAP_", colnames(x = lda.posterior))
      rownames(x = feature.loadings) <- var_names
      colnames(x = feature.loadings) <- paste0(reduction.key, 1:ncol(x = cell.embeddings))
      rownames(x = cell.embeddings) <- rownames(x = object)
      colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
      reduction.data <- CreateDimReducObject(
        embeddings = cell.embeddings,
        loadings = feature.loadings,
        assay = assay,
        key = reduction.key,
        misc = list(
          assignments = lda.assignments,
          posterior = lda.posterior,
          model = lda_results,
          cv = lda_cv
        )
      )
      if (verbose) {
        print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
      }

object_lda = reduction.data
T2_K2_integrated2 =T2_K2_integrated
T2_K2_integrated2[["lda_atac"]] <- object_lda

T2_K2_integrated2 <- RunUMAP(
  object = T2_K2_integrated2,
  dims = 1:15,
  reduction = 'lda_atac',
  reduction.key = 'ldaumap_atac_',
  reduction.name = 'ldaumap_atac')


all_cells <- Cells(T2_K2_integrated2)

perturb_values <- as.data.frame(T2_K2_integrated2@meta.data$perturb)
ident = as.data.frame(T2_K2_integrated2@active.ident)
rownames(perturb_values) = rownames(ident)
perturb_values = cbind(rownames(perturb_values),perturb_values)
colnames(perturb_values)[1] = 'id'

x = perturb_values[!( perturb_values$`T2_K2_integrated2@meta.data$perturb` %in% c(
  "NT" ,"LINC00511" , "NEAT1" ,"XIST","FAM66c" ,"MIR155HG",
  "MIR222HG","LINC01638","GAS5" ,"LINC00578" ,"HOTAIRM1",
  "LINC02693","LINC01578",'PURPL',"MEG3", "SNHG1")),]

cells_to_remove = x$id
all_cells <- Cells(T2_K2_integrated2)

cells_to_keep <- setdiff(all_cells, cells_to_remove)

T2_K2_integrated2 <- subset(T2_K2_integrated2, cells = cells_to_keep)


cell_type_cols <- c("NT" = "#AFB1B7","LINC00511" = "#BA4A94", 
                    "NEAT1" =  "#4AA123","XIST" = "#CB1D1B",
                    "FAM66c" = "#7573AD","MIR155HG" = "#1927aa",
                    "MIR222HG" =  "#DDE9B0","LINC01638" = "#72bcd9",
                    "GAS5" = "#F4C685",
                    
                    "LINC00578" = "#F7BEAE","HOTAIRM1" = "#e56e24",
                    "LINC02693" = "#7BAE92",
                    "LINC01578" = '#2B543A',
                    'PURPL'="#812F33","MEG3" =  "#C9A47D", "SNHG1" = '#EAAFCD'
                    # "LINC01012" = "#7A8EBA", "LRRC75a" = "#C7BCB1","MIAT" = "#6C4C73",
                    #"RMRP" = "#b4446c","FAM66C" ="#EAD9E1" ,"SNHG3" = "#4D4C72"
)

Idents(T2_K2_integrated2) <- "perturb"
plot7 <- DimPlot(T2_K2_integrated2, label = F, 
                 reduction = "ldaumap_atac", repel = T, 
                 pt.size = 0.2,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot7

##

cell_type_cols <- c("K2" = "#3081BC",
                    "T2" =  "#AC4A4B"
                    
)


Idents(T2_K2_integrated2) <- "batch"
plot8 <- DimPlot(T2_K2_integrated2, label = F, 
                 reduction = "ldaumap_atac", repel = T, 
                 pt.size = 0.2,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot8 

#####integrated_atac_rna_embedding

T2_K2_integrated3 = T2_K2_integrated2


object = T2_K2_integrated3

assay = 'RNA'
de.assay = "RNA"
reduction.key = "LDA_"
seed = 42
pc.assay = "PRTB_inte"
labels = "perturb"
nt.label = "NT"
npcs = 30
verbose = TRUE
logfc.threshold = 0.25

library(pbapply)

      projected_pcs <- list()
      gene_list <- setdiff(x = unique(x = object[[labels]][, 1]), y = nt.label)
      Idents(object = object) <- labels
      DefaultAssay(object = object) <- pc.assay
      all_genes <- list()
      nt.cells <- WhichCells(object = object, idents = nt.label)
      for (g in gene_list) {
        if (verbose) {
          message(g)
        }
        gd.cells <- WhichCells(object = object, idents = g)
        gene_set <- TopDEGenesMixscape(
          object = object,
          ident.1 = gd.cells,
          ident.2 = nt.cells,
          de.assay = de.assay,
          logfc.threshold = logfc.threshold,
          labels = labels,
          verbose = verbose
        )
        if (length(x = gene_set) < (npcs + 1)) {
          all_genes[[g]] <- character()
          next
        }
        all_genes[[g]] <- gene_set
      }
      all_markers <- unique(x = unlist(x = all_genes))
      missing_genes <- all_markers[!all_markers %in% rownames(x = object[[pc.assay]])]
      
      
      object <- GetMissingPerturb(object = object, assay = pc.assay, features = missing_genes, verbose = verbose)
      for (g in gene_list) {
        if (verbose) {
          message(g)
        }
        gene_subset <- subset(x = object, idents = c(g, nt.label))
        gene_set <- all_genes[[g]]
        if (length(x = gene_set) == 0) {
          next
        }
        gene_subset <- ScaleData(
          object = gene_subset,
          features = gene_set,
          verbose = FALSE
        )
        gene_subset <- RunPCA(
          object = gene_subset,
          features = gene_set,
          npcs = npcs,
          verbose = FALSE
        )
        project_pca <- ProjectCellEmbeddings(
          reference = gene_subset,
          query = object,
          dims = 1:npcs,
          verbose = FALSE
        )
        colnames(x = project_pca) <- paste(g, colnames(x = project_pca), sep = "_")
        projected_pcs[[g]] <- project_pca
      }
      return(projected_pcs)

projected_pcs

lda.lables <- object[[labels]][,]

object = projected_pcs
labels = lda.lables
assay = "RNA"

verbose = TRUE
ndims.print = 1:5
nfeatures.print = 30
reduction.key = "LDA_"
seed = 42

library(MASS)

      if (!is.null(x = seed)) {
        set.seed(seed = seed)
      }
      object <- data.frame(object)
      var_names <- colnames(x = object)
      object$lda_cluster_label <- labels
      lda_results <- lda(formula = lda_cluster_label ~ ., data = object)
      lda_predictions <- predict(object = lda_results, newdata = object)
      
      lda_cv <-lda(
        formula = lda_cluster_label ~ .,
        data = object,
        CV = TRUE
      )$posterior
      feature.loadings <- lda_results$scaling
      cell.embeddings <- lda_predictions$x
      lda.assignments <- lda_predictions$class
      lda.posterior <- lda_predictions$posterior
      colnames(x = lda.posterior) <- paste0("LDAP_", colnames(x = lda.posterior))
      rownames(x = feature.loadings) <- var_names
      colnames(x = feature.loadings) <- paste0(reduction.key, 1:ncol(x = cell.embeddings))
      rownames(x = cell.embeddings) <- rownames(x = object)
      colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
      reduction.data <- CreateDimReducObject(
        embeddings = cell.embeddings,
        loadings = feature.loadings,
        assay = assay,
        key = reduction.key,
        misc = list(
          assignments = lda.assignments,
          posterior = lda.posterior,
          model = lda_results,
          cv = lda_cv
        )
      )
      if (verbose) {
        print(x = reduction.data, dims = ndims.print, nfeatures = nfeatures.print)
      }

object_lda = reduction.data

T2_K2_integrated3[["lda_rna"]] <- object_lda


T2_K2_integrated3 <- RunUMAP(
  object = T2_K2_integrated3,
  dims = 1:15,
  reduction = 'lda_rna',
  reduction.key = 'ldaumap_rna_',
  reduction.name = 'ldaumap_rna')


perturb_values <- as.data.frame(T2_K2_integrated3@meta.data$perturb)
ident = as.data.frame(T2_K2_integrated3@active.ident)
rownames(perturb_values) = rownames(ident)
perturb_values = cbind(rownames(perturb_values),perturb_values)
colnames(perturb_values)[1] = 'id'

x = perturb_values[!( perturb_values$`T2_K2_integrated3@meta.data$perturb` %in% c(
  "NT" ,"LINC00511" , "NEAT1" ,"XIST","FAM66c" ,"MIR155HG",
  "MIR222HG","LINC01638","GAS5" ,"LINC00578" ,"HOTAIRM1",
  "LINC02693","LINC01578",'PURPL',"MEG3", "SNHG1")),]

cells_to_remove = x$id
all_cells <- Cells(T2_K2_integrated3)

cells_to_keep <- setdiff(all_cells, cells_to_remove)

T2_K2_integrated3 <- subset(T2_K2_integrated3, cells = cells_to_keep)

cell_type_cols <- c("NT" = "#AFB1B7","LINC00511" = "#BA4A94", 
                    "NEAT1" =  "#4AA123","XIST" = "#CB1D1B",
                    "FAM66c" = "#7573AD","MIR155HG" = "#1927aa",
                    "MIR222HG" =  "#DDE9B0","LINC01638" = "#72bcd9",
                    "GAS5" = "#F4C685",
                    
                    "LINC00578" = "#F7BEAE","HOTAIRM1" = "#e56e24",
                    "LINC02693" = "#7BAE92",
                    "LINC01578" = '#2B543A',
                    'PURPL'="#812F33","MEG3" =  "#C9A47D", "SNHG1" = '#EAAFCD'
                    # "LINC01012" = "#7A8EBA", "LRRC75a" = "#C7BCB1","MIAT" = "#6C4C73",
                    #"RMRP" = "#b4446c","FAM66C" ="#EAD9E1" ,"SNHG3" = "#4D4C72"
)

Idents(T2_K2_integrated3) <- "perturb"
plot7 <- DimPlot(T2_K2_integrated3, label = F, 
                 reduction = "ldaumap_rna", repel = T, 
                 pt.size = 0.2,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot7

cell_type_cols <- c("K2" = "#3081BC",
                    "T2" =  "#AC4A4B"
                    
)


Idents(T2_K2_integrated3) <- "batch"
plot8 <- DimPlot(T2_K2_integrated3, label = F, 
                 reduction = "ldaumap_rna", repel = T, 
                 pt.size = 0.2,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot8 

# WNN ####inte_atac_rna
T2_K2_integrated3 <- FindMultiModalNeighbors(
  T2_K2_integrated3,
  k.nn = 10,
  reduction.list = list("lda_rna", "lda_atac"),
  dims.list = list(1:15, 2:20),  # 设置使用的维度数
  weighted.nn.name = "weighted.nn",
  knn.graph.name = "wknn",
  snn.graph.name = "wsnn",
  knn.range = 200,
  prune.SNN = 1/150,
  sd.scale = 1,
  smooth = F
  
  
)

T2_K2_integrated3  <- RunUMAP(
  object = T2_K2_integrated3,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(T2_K2_integrated3, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()


library(scales)
library(ggplot2)
# Set colors for each perturbation.
col = setNames(object = hue_pal()(22),nm = levels(scmulti_bj_single_test$perturb))
names(col) <- c("all_control_0.05",names(col)[2:16])
col[1] <- "grey39"


p <- DimPlot(object = scmulti_bj_single_test, 
             reduction = "umap", 
             repel = T, 
             label.size = 5, 
             label = T, 
             cols = col) + NoLegend()


p2 <- p+ 
  scale_color_manual(values=col, drop=FALSE) + 
  ylab("UMAP 2") +
  xlab("UMAP 1") 

p2

cell_type_cols <- c("NT" = "#AFB1B7","LINC00511" = "#BA4A94", 
                    "NEAT1" =  "#4AA123","XIST" = "#CB1D1B",
                    "FAM66c" = "#7573AD","MIR155HG" = "#1927aa",
                    "MIR222HG" =  "#DDE9B0","LINC01638" = "#72bcd9",
                    "GAS5" = "#F4C685",
                    
                    "LINC00578" = "#F7BEAE","HOTAIRM1" = "#e56e24",
                    "LINC02693" = "#7BAE92",
                    "LINC01578" = '#2B543A',
                    'PURPL'="#812F33","MEG3" =  "#C9A47D", "SNHG1" = '#EAAFCD'
                    # "LINC01012" = "#7A8EBA", "LRRC75a" = "#C7BCB1","MIAT" = "#6C4C73",
                    #"RMRP" = "#b4446c","FAM66C" ="#EAD9E1" ,"SNHG3" = "#4D4C72"
)

Idents(T2_K2_integrated3) <- "perturb"
plot7 <- DimPlot(T2_K2_integrated3, label = F, 
                 reduction = "umap", repel = T, 
                 pt.size = 0.15,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot7

cell_type_cols <- c("K2" = "#3081BC",
                    "T2" =  "#AC4A4B"
                    
)


Idents(T2_K2_integrated3) <- "batch"
plot8 <- DimPlot(T2_K2_integrated3, label = F, 
                 reduction = "umap", repel = T, 
                 pt.size = 0.2,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
plot8 



##quantify
emb=Embeddings(T2_K2_integrated4, reduction = "umap")
perturb = T2_K2_integrated4@meta.data$perturb
batch = T2_K2_integrated4@meta.data$batch

out = cbind(emb,perturb,batch)

##欧氏距离
#rna
obj <- T2_K2_integrated4
red_name <- "umap"        
cluster_col <- "perturb"  #

# UMAP 坐标与 cluster 标签
emb <- Embeddings(obj, red_name)            
clu <- Idents(obj)                           

# 2) cluster质心
min_cells <- 1
clu_lvls <- levels(clu)
clu_lvls <- clu_lvls[sapply(clu_lvls, function(cl) sum(clu == cl) >= min_cells)]

centroids <- t(sapply(clu_lvls, function(cl) {
  colMeans(emb[clu == cl, , drop = FALSE])
}))

# 质心欧距矩阵
dist_mat <- as.matrix(dist(centroids, method = "euclidean"))


##per-cluster silhouette
library(cluster)
library(dplyr)
library(ggplot2)
library(readr)

obj <- T2_K2_integrated4        
reduction  <- "umap"            
Idents(obj) <- "perturb"               
# 或者：cluster_col <- "perturb"       #   如果想用 meta.data 里的某一列

prefix     <- "RNA_sil"                
min_cells  <- 2                        

## ==== 1. 取嵌入矩阵 ====
emb <- Embeddings(obj, reduction)
stopifnot(!is.null(emb), nrow(emb) > 0)

## ==== 2. 取簇标签 ====
# 用 Idents：
cl <- Idents(obj)

# 对齐顺序（以嵌入的行名为准）
cl <- cl[rownames(emb)]

## ==== 3. 过滤过小簇 ====
tab <- table(cl)
keep_levels <- names(tab[tab >= min_cells])
keep_idx <- cl %in% keep_levels
emb_f <- emb[keep_idx, , drop = FALSE]
cl_f  <- droplevels(cl[keep_idx])
stopifnot(nlevels(cl_f) >= 2)

## ==== 4. 计算 silhouette ====
d <- dist(emb_f, method = "euclidean")
sil <- silhouette(as.integer(cl_f), d)
sil_df <- as.data.frame(sil)
sil_df$cell    <- rownames(emb_f)
sil_df$cluster <- levels(cl_f)[sil_df$cluster]
names(sil_df)[names(sil_df) == "sil_width"] <- "silhouette"

## ==== 5. 每簇汇总 ====
sil_sum <- sil_df |>
  group_by(cluster) |>
  summarise(n_cells   = n(),
            mean_sil  = mean(silhouette),
            median_sil= median(silhouette),
            sd_sil    = sd(silhouette),
            .groups = "drop") |>
  arrange(desc(mean_sil))

####
library(ComplexHeatmap)
library(circlize)
dist_mat <- as.matrix(dist)

# 定义颜色梯度（根据欧式距离范围调整）
max_val <- max(dist_mat)
col_fun = colorRamp2(c(0, max_val/2, max_val),
                     c("#ffffff", "#cbc9e2", "#6a51a3"))

# 画热图
Heatmap(
  dist_mat,
  name = " ",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  border = TRUE,
  col = col_fun,
  column_title = "Euclidean Distance"
)


