library(ggplot2)
library(Seurat)

DefaultAssay( scmulti4) = 'newcharmVAR'
scmulti4 <- CalcPerturbSig(object = scmulti4,
                           assay = "newcharmVAR",
                           slot = "data",
                           gd.class ="NT_sig_1",
                           nt.cell.class = "NT",
                           reduction = "pca",
                           split.by = "replicate",
                           ndims = 15,
                           num.neighbors = 10,
                           new.assay.name = "PRTB_newcharmVAR",
)

VariableFeatures(object = scmulti4) <- VariableFeatures(object = scmulti4[["newcharmVAR"]])
#scmulti4 <- ScaleData(object = scmulti4, do.scale = F, do.center = T)

# Run PCA to reduce the dimensionality of the data.
scmulti4 <- RunPCA(object = scmulti4, reduction.key = 'prtbpca', reduction.name = 'prtbpca')

# Run UMAP to visualize clustering in 2-D.
scmulti4 <- RunUMAP(
  object = scmulti4, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')

#
sub <- subset(scmulti4, idents = c("NT","NEAT1","LINC00578","LINC01578","FAM66C", "MIR155HG",
                                   "MIR222HG","SNHG1","PURPL",
                                   #"BRD4_KD",
                                   "LINC01638","XIST","LINC00511",
                                   "MEG3","GAS5","HOTAIRM1","LINC02693"))
sub <- AddMetaData(object = sub,     #seurat对象
                   metadata = Idents(sub),    #需要添加的metadata
                   col.name = "ID")
# Run LDA.


#####

object = sub
assay = "newcharmVAR"
pc.assay = "PRTB_newcharmVAR"
labels = "ID" 
nt.label = "all_control_0.05" 
npcs = 10
logfc.threshold = 0.25
verbose = T
pval.cutoff = 0.05
de.assay = "PRTB_newcharmVAR"
test.use = "LR"

##
TopDEGenesMixscape <- function(
    object,
    ident.1,
    ident.2 = NULL,
    labels = "ID",
    de.assay = "newcharmVAR",
    test.use = "LR",
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
        verbose = verbose
      )
      de.genes <- de.genes[de.genes$p_val < pval.cutoff, ]
    },
    error = function(e) {}
  )
  return(rownames(x = de.genes))
}
#
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
#
ProjectCellEmbeddings <- function(
    reference,
    query,
    reduction = "pca",
    reference.assay = NULL,
    query.assay = NULL,
    dims = 1:50,
    scale = TRUE,
    verbose = TRUE,
    feature.mean = NULL,
    feature.sd = NULL
) {
  if (verbose) {
    message("Projecting cell embeddings")
  }
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  features <- rownames(x = Loadings(object = reference[[reduction]]))
  features <- intersect(x = features, y = rownames(x = query[[query.assay]]))
  proj.data <- RescaleQuery(
    reference = reference,
    query = query,
    features = features,
    scale = scale,
    feature.mean = feature.mean,
    feature.sd = feature.sd
  )
  ref.feature.loadings <- Loadings(object = reference[[reduction]])[features, dims]
  proj.pca <- t(crossprod(x = ref.feature.loadings, y = proj.data))
  return(proj.pca)
}
#
RescaleQuery <- function(
    reference,
    query,
    reference.assay = NULL,
    query.assay = NULL,
    features = NULL,
    feature.mean = NULL,
    feature.sd = NULL,
    scale = TRUE
) {
  reference.assay <- reference.assay %||% DefaultAssay(object = reference)
  query.assay <- query.assay %||% DefaultAssay(object = query)
  features <- features %||% intersect(
    rownames(x = reference[[reference.assay]]),
    rownames(x = query[[query.assay]])
  )
  reference.data <-  GetAssayData(
    object = reference,
    assay = reference.assay,
    slot = "data")[features, ]
  query.data <- GetAssayData(
    object = query,
    assay = query.assay,
    slot = "data")[features, ]
  if (is.null(x = feature.mean)) {
    feature.mean <- rowMeans(x = reference.data)
    if (scale) {
      feature.sd <- sqrt(
        x = SparseRowVar2(
          mat = as.sparse(x = reference.data),
          mu = feature.mean,
          display_progress = FALSE
        )
      )
      feature.sd[is.na(x = feature.sd)] <- 1
    } else {
      feature.sd <- rep(x = 1, nrow( reference.data))
    }
    feature.mean[is.na(x = feature.mean)] <- 1
  }
  proj.data <- GetAssayData(
    object = query,
    assay = query.assay,
    slot = "data"
  )[features, ]
  store.names <- dimnames(x = proj.data)
  if (is.numeric(x = feature.mean) && feature.mean[[1]] != "SCT") {
    proj.data <- FastSparseRowScaleWithKnownStats(
      mat = as.sparse(x = proj.data),
      mu = feature.mean,
      sigma = feature.sd,
      display_progress = FALSE
    )
  }
  dimnames(x = proj.data) <- store.names
  return(proj.data)
}
#
SparseRowVar2 <- function(mat, mu, display_progress) {
  .Call('_Seurat_SparseRowVar2', PACKAGE = 'Seurat', mat, mu, display_progress)
}
FastSparseRowScaleWithKnownStats <- function(mat, mu, sigma, scale = TRUE, center = TRUE, scale_max = 10, display_progress = TRUE) {
  .Call('_Seurat_FastSparseRowScaleWithKnownStats', PACKAGE = 'Seurat', mat, mu, sigma, scale, center, scale_max, display_progress)
}

#
projected_pcs <- list()
gene_list <- setdiff(x = unique(x = object[[labels]][, 1]), 
                     y = nt.label)
Idents(object = object) <- labels
DefaultAssay(object = object) <- pc.assay
#DefaultAssay(object = object) <- de.assay
all_genes <- list()
nt.cells <- WhichCells(object = object, idents = nt.label)
#

for (g in gene_list) {
  if (verbose) {
    message(g)
  }
  gd.cells <- WhichCells(object = object, idents = g)
  gene_set <- TopDEGenesMixscape(object = object, ident.1 = gd.cells, 
                                 ident.2 = nt.cells, de.assay = de.assay, logfc.threshold = logfc.threshold, 
                                 labels = labels, verbose = verbose)
  if (length(x = gene_set) < (npcs + 1)) {
    all_genes[[g]] <- character()
    next
  }
  all_genes[[g]] <- gene_set
}

all_markers <- unique(x = unlist(x = all_genes))
missing_genes <- all_markers[!all_markers %in% rownames(x = object[[pc.assay]])]
object <- GetMissingPerturb(object = object, assay = pc.assay, 
                            features = missing_genes, verbose = verbose)
for (g in gene_list) {
  if (verbose) {
    message(g)
  }
  gene_subset <- subset(x = object, idents = c(g, nt.label))
  gene_set <- all_genes[[g]]
  if (length(x = gene_set) == 0) {
    next
  }
  gene_subset <- ScaleData(object = gene_subset, features = gene_set, 
                           verbose = FALSE)
  gene_subset <- RunPCA(object = gene_subset, features = gene_set, 
                        npcs = npcs, verbose = FALSE)
  project_pca <- ProjectCellEmbeddings(reference = gene_subset, 
                                       query = object, dims = 1:npcs, verbose = FALSE)
  colnames(x = project_pca) <- paste(g, colnames(x = project_pca), 
                                     sep = "_")
  projected_pcs[[g]] <- project_pca
}


##
lda.lables <- object[[labels]][, ]
object_lda <- RunLDA(object = projected_pcs, labels = lda.lables, 
                     assay = assay, verbose = verbose)
#object[["lda"]] <- object_lda
sub[["lda"]] <- object_lda

DefaultAssay(sub) =  'newcharmVAR'
sub <- RunUMAP(
  object = sub,
  dims = 1:11,
  reduction = 'lda',
  reduction.key = 'ldaumap',
  reduction.name = 'ldaumap')

Idents(sub) <- "ID"
sub$mixscape_class3 <- as.factor(sub$ID)


library(scales)
library(ggplot2)

cell_type_cols <- c("NT" = "#AFB1B7","NEAT1" = "#9c52f2", 
                    "LINC00578" =  "#4AA123","LINC01578" = "#CC071C",
                    "FAM66C" = "#7573AD","MIR155HG" = "#210aee",
                    "MIR222HG" =  "#DDE9B0","SNHG1" = "#00ddff",
                    "PURPL" = "#EF5376","LINC01638" = "#C9A47D",
                    "XIST" = "#F7BEAE","LINC00511" = "#e56e24",
                    "MEG3" = "#ffc800","GAS5" = "#7BAE92",
                    "HOTAIRM1" = "#74507B","LINC02693" = "#7A8EBA")


plot7 <- DimPlot(sub, label = F, 
                 reduction = "ldaumap", repel = T, 
                 pt.size = 0.8,cols = cell_type_cols)+
  labs(x = "UMAP1", y = "UMAP2",title = "Celltype") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

plot7


