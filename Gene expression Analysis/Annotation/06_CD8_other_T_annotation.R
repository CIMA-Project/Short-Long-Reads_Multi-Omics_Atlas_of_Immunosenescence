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

# Function
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

# L2
all_cell_data <- readRDS("/data/RDS/Annotation/4_PBMC_annotation_L1.rds")
CD8_T_otherT <- subset(all_cell_data, subset = celltype == "CD8 T & other T")
CD8_T_otherT <- subcluster(CD8_T_otherT, dims_use = 1:8, resolution = 1)
CD8T_otherT_marker <- c("TRABD2A", "MAL", "CCR7", "SELL", "LEF1", "TCF7", "GPR183",
                        "CD44", "S100A4", "GZMA", "GZMH", "GZMK", "NKG7", "GNLY", "MKI67",
                        "CABP4", "NOTCH1", "TRDC", "SLC4A10", "KLRB1", "KLRF1", "KLRC2", "TYROBP", "FCGR3A")
FeaturePlot(CD8_T_otherT, features = CD8T_otherT_marker, ncol = 4)
CD8_T_otherT <- RenameIdents(CD8_T_otherT,
                             '0' = "NKT", '1' = "NKT",
                             '2' = "CD8 CTL", '3' = "CD8 Naïve", '4' = "CD8 CTL", '5' = "CD8 CTL",
                             '6' = "CD8 Memory", '7' = "CD8 Naïve", '8' = "CD8 Memory",
                             '9' = "MAIT", '10' = "CD8 CTL", '11' = "gdT", '12' = "Cycling T")
CD8_T_otherT$L3_celltype <- Idents(CD8_T_otherT)
#L3
CD8_Memory <- subset(CD8_T_otherT, subset = L3_celltype == "CD8 Memory")
CD8_Memory <- subcluster(CD8_Memory, dims_use = 1:4, resolution = 0.5)
FeaturePlot(CD8_Memory, features = "CCR7", pt.size = 0.5)
CD8_Memory <- RenameIdents(CD8_Memory,
                           '0' = "CD8 Tem", '1' = "CD8 Tcm",
                           '2' = "CD8 Tem", '3' = "CD8 Tcm",
                           '4' = "CD8 Tem", '5' = "CD8 Tem")
CD8_Memory$L3_celltype <- Idents(CD8_Memory)
CD8_T_otherT$L3_celltype[colnames(CD8_Memory)] <- CD8_Memory$L3_celltype
saveRDS(CD8_T_otherT, "/data/RDS/Annotation/7_CD8T_otherT_annotation.rds")
