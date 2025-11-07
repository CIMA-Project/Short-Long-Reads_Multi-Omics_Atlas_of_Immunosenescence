library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(R.utils)
.cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                   "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                   "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                   "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                   "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                   "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

#  Function for sub-clustering
subcluster <- function(seurat_obj, dims_use, resolution, nfeatures = 2000, n_neigh = 30, min_dist = 0.3) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident", max_iter = 30) # nolint
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = dims_use)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = dims_use, reduction = "harmony",
                        n.neighbors = n_neigh, min.dist = min_dist)
  return(seurat_obj)
}
# NK annotation
all_cell_data <- readRDS("/data/RDS/Annotation/4_PBMC_annotation_L1.rds")

NK <- subset(all_cell_data, subset = celltype == "NK")
NK <- subcluster(NK, dims_use = 1:10, resolution = 0.8, nfeatures = 1500)

NK_marker <- c("GATA3", "KIT", "FCGR3A", "NCAM1", "CD160", "KLRC1", "XCL1", "XCL2",
               "GZMK", "B3GAT1", "CX3CR1", "MKI67", "TYMS", "STMN1")
DimPlot(NK, label = TRUE, cols = .cluster_cols, pt.size = 0.5)
FeaturePlot(NK, features = NK_marker, ncol = 4)
# NKdim-1/NKdim-2/NKbright
NK <- RenameIdents(NK,
                       '0' = "NKdim-2", '1' = "NKdim-2", '2' = "NKdim-2", '3' = "NKdim-2", '4' = "NKdim-2",
                       '5' = "NKdim-1", '6' = "NKdim-2", '7' = "NKbright", '8' = "NKdim-2")
NK$L3_celltype <- Idents(NK)
#NKdim
NKdim_clusters <- c("0", "1", "2", "3", "4", "5", "6", "8", "9")
NKdim <- subset(NK, subset = seurat_clusters %in% NKdim_clusters)
NKdim <- subcluster(NKdim, dims_use = 1:10, resolution = 0.8, nfeatures = 1500)
DimPlot(NKdim, label = TRUE, cols = .cluster_cols, pt.size = 0.5)
FeaturePlot(NKdim, features = NK_marker, ncol = 4)
NKdim <- RenameIdents(NKdim,
                      '0' = "Terminal NK", '1' = "XCL-NKdim", '2' = "XCL-NKdim", '3' = "Terminal NK",
                      '4' = "XCL-NKdim", '5' = "XCL+NKdim", '6' = "Terminal NK", '7' = "Terminal NK")
NKdim$L3_celltype <- Idents(NKdim)
NK$L3_celltype[colnames(NKdim)] <- NKdim$L3_celltype
saveRDS(NK, file = "/data/RDS/Annotation/9_NK_annotation_end.rds")
