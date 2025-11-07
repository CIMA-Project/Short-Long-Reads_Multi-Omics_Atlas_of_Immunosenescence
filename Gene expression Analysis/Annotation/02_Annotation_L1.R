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
# =======================================
pbmc <- readRDS("/data/RDS/pbmc_remove_doublets.rds")
#genes to exclude
exclude_patterns <- c(
  "^RPL", "^RPS", "\\d+\\.\\d+", "^RNU.*P$", "^RN7SL.*P$", "^LINC", "^MT"
)
exclude_pattern <- paste(exclude_patterns, collapse = "|")
genes_to_exclude <- rownames(pbmc)[grepl(exclude_pattern, rownames(pbmc))]
#TCR/BCR/HIST genes from files
bcrgene  <- read.table("/data/Files/Annotation/BCR_gene.txt", stringsAsFactors = FALSE)
histgene <- read.table("/data/Files/Annotation/HIST_gene.txt", stringsAsFactors = FALSE)
tcr_gene <- read.table("/data/Files/Annotation/TCR_gene.txt", stringsAsFactors = FALSE)
file_genes <- c(bcrgene[[1]], histgene[[1]], tcr_gene[[1]])
#combine all excluded genes and save
all_excluded_genes <- union(genes_to_exclude, file_genes)
write.csv(all_excluded_genes, file = "/data/Files/Annotation/all_excluded_genes.csv", row.names = FALSE)

# Normalize and identify variable features while excluding unwanted genes
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10, n.neighbors = 30, min.dist = 0.5)
# Harmony
pbmc <- RunHarmony(pbmc, group.by.vars = "orig.ident", max_iter = 30)
DimPlot(pbmc, reduction = "harmony", group.by = "group", pt.size = 0.5)
ggsave("/data/Figure/Annotation/pbmc_after_harmony_by_group.png", width = 9, height = 8, dpi = 300)

# UMAP & Clustering after Harmony
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:16)
pbmc <- FindClusters(pbmc, resolution = 1.5)
pbmc <- RunUMAP(pbmc, dims = 1:16, reduction = "harmony", n.neighbors = 30, min.dist = 0.5)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
L1_marker <- c("CD3D","CD4","CD8A","CD8B","NCAM1","FCGR3A","NCR1","CD79A","CD34","CD68","AIF1","IL7R","KLRC1","KLRD1")
FeaturePlot(pbmc, features = L1_marker, ncol = 4)
# define platelet clusters and remove
Idents(pbmc) <- "seurat_clusters"
pbmc <- RenameIdents(pbmc,
  '15' = "platelet", '20' = "platelet", '25' = "platelet", '27' = "platelet",
  .default = "PBMC"
)
pbmc$celltype <- Idents(pbmc)
DimPlot(pbmc, reduction = "umap", pt.size = 0.5)
ggsave("/data/Figure/Annotation/PBMC_with_platelet.pdf", height = 8, width = 8)

pbmc <- subset(pbmc, subset = celltype == "PBMC")
saveRDS(pbmc, "/data/RDS/Annotation/1_PBMC_without_platelet.rds")
# annotation step 1: major cell types (B, T&NK, Myeloid)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunHarmony(pbmc, group.by.vars = "orig.ident", max_iter = 30)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:14)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- RunUMAP(pbmc, dims = 1:14, reduction = "harmony", n.neighbors = 30, min.dist = 0.5)
pbmc <- RenameIdents(pbmc,
  '4' = "Myeloid", '5' = "B", '8' = "B", '12' = "Myeloid", '13' = "Myeloid",
  '14' = "Myeloid", '15' = "B", '17' = "Myeloid", '18' = "B", '19' = "Myeloid",
  .default = "T&NK"
)
pbmc$celltype <- Idents(pbmc)
saveRDS(pbmc, "/data/RDS/Annotation/2_PBMC_annotation_L1.rds")

# Below is the sub-clustering and annotation of major cell types
# B subclustering
B_sub <- subset(pbmc, subset = celltype == "B") 
B_sub <- NormalizeData(B_sub)
B_sub <- FindVariableFeatures(B_sub, selection.method = "vst", nfeatures = 2000)
B_sub <- ScaleData(B_sub)
B_sub <- RunPCA(B_sub, features = VariableFeatures(object = B_sub))
B_sub <- RunHarmony(B_sub, group.by.vars = "orig.ident", max_iter = 30)
B_sub <- FindNeighbors(B_sub, dims = 1:7, reduction = "harmony")
B_sub <- FindClusters(B_sub, resolution = 0.8)
B_sub <- RunUMAP(B_sub, dims = 1:7, reduction = "harmony", n.neighbors = 30, min.dist = 0.5)
DimPlot(B_sub, reduction = "umap", label=T, cols = .cluster_cols,pt.size=0.1)
L1_marker <- c("CD3D", "CD4", "CD8A", "CD8B", "NCAM1",
               "FCGR3A", "NCR1", "CD79A", "CD34", "CD68",
               "AIF1", "IL7R", "KLRC1", "KLRD1")
FeaturePlot(B_sub, feature = L1_marker ,reduction = "umap",ncol =4, raster=FALSE)               
B_sub <- RenameIdents(B_sub, '7' = "B_contaminate_CD3D", .default = "B")
pbmc$celltype[colnames(B_sub)] <- Idents(B_sub)

# Myeloid subclustering
Myeloid_sub <- subset(pbmc, subset = celltype == "Myeloid")
Myeloid_sub <- NormalizeData(Myeloid_sub)
Myeloid_sub <- FindVariableFeatures(Myeloid_sub, selection.method = "vst", nfeatures = 2000)
Myeloid_sub <- ScaleData(Myeloid_sub)
Myeloid_sub <- RunPCA(Myeloid_sub, features = VariableFeatures(object = Myeloid_sub))
Myeloid_sub <- RunHarmony(Myeloid_sub, group.by.vars = "orig.ident", max_iter = 30)
Myeloid_sub <- FindNeighbors(Myeloid_sub, dims = 1:17, reduction = "harmony")
Myeloid_sub <- FindClusters(Myeloid_sub, resolution = 0.8)
Myeloid_sub <- RunUMAP(Myeloid_sub, dims = 1:17, reduction = "harmony", n.neighbors = 30, min.dist = 0.5)
DimPlot(Myeloid_sub, reduction = "umap", label=T, cols = .cluster_cols,pt.size=0.1) 
FeaturePlot(Myeloid_sub, feature = L1_marker ,reduction = "umap",ncol =4, raster=FALSE)
Myeloid_sub <- RenameIdents(Myeloid_sub,
  '8' = "Myeloid_contaminate_CD3D", '13' = "Myeloid_contaminate_CD3D",
  '14' = "Myeloid_contaminate_CD79A", .default = "Myeloid"
)
pbmc$celltype[colnames(Myeloid_sub)] <- Idents(Myeloid_sub)

# Eextract major cell types for final annotation
pbmc <- subset(pbmc, subset = celltype %in% c("B", "T&NK", "Myeloid")) 
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunHarmony(pbmc, group.by.vars = "orig.ident", max_iter = 30)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:14)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- RunUMAP(pbmc, dims = 1:14, reduction = "harmony", n.neighbors = 30, min.dist = 0.5)
DimPlot(pbmc, reduction = "umap", label = T, cols = .cluster_cols, pt.size = 0.5)
FeaturePlot(pbmc, feature = L1_marker, reduction = "umap", ncol = 4, raster = FALSE)
pbmc <- RenameIdents(pbmc,
  '0' = "CD4 T", '1' = "CD4 T", '4' = "Myeloid", '5' = "B", '7' = "B", 
  '8' = "CD4 T", '10' = "Myeloid", '13' = "Myeloid", '14' = "B",
  '16' = "Myeloid", '17' = "HSPC", .default = "T&NK"
)
pbmc$celltype <- Idents(pbmc)
saveRDS(pbmc, "/data/RDS/Annotation/3_PBMC_annotation_L1.rds")
