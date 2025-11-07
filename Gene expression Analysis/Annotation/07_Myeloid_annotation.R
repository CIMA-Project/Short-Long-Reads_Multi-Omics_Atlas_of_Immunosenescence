library(Seurat)
library(harmony)

# Function
subcluster <- function(seurat_obj, dims_use, resolution, nfeatures = 2000, n_neigh = 30, min_dist = 0.5) {
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

# L2 Annotation
all_cell_data <- readRDS("/data/RDS/Annotation/4_PBMC_annotation_L1.rds")
Myeloid <- subset(all_cell_data, subset = celltype == "Myeloid")
Myeloid <- subcluster(Myeloid, dims_use = 1:13, resolution = 0.8)
Myeloid_marker <- c("GP9", "TUBB1", "PF4", "PPBP", "ITGAX", "CD74", "HLA-DPB1",
                    "FCER1A", "CD1C", "BATF3", "IDO1", "LILRA4", "IRF4", "SPIB", "IRF8",
                    "IL17RA", "FCN1", "LYZ", "CTSD")
FeaturePlot(Myeloid, features = Myeloid_marker, ncol = 4)
Myeloid <- RenameIdents(Myeloid,
                        '0' = "Monocyte",'1' = "Monocyte",'2' = "Monocyte",'3' = "Monocyte",'4' = "Monocyte",
                        '5' = "cDC",'6' = "cDC",'7' = "Monocyte",'8' = "pDC",
                        '9' = "Monocyte",'10' = "Monocyte",'11' = "MK",'12' = "cDC")
Myeloid$L3_celltype <- Idents(Myeloid)

# cDC Annotation
cDC <- subset(Myeloid, subset = L3_celltype == "cDC")
cDC <- subcluster(cDC, dims_use = 1:7, resolution = 0.5)
FeaturePlot(cDC, features = Myeloid_marker, ncol = 4)
cDC <- RenameIdents(cDC,
                    '0' = "cDC2",'1' = "cDC2",'2' = "cDC2",'3' = "cDC2",'4' = "cDC1",'5' = "cDC2")
cDC$L3_celltype <- Idents(cDC)
Myeloid$L3_celltype[colnames(cDC)] <- cDC$L3_celltype
# Monocyte Annotation
Mono <- subset(Myeloid, subset = L3_celltype == "Monocyte")
Mono <- subcluster(Mono, dims_use = 1:12, resolution = 0.5)
L3_Myeloid_marker <- c("FCN1", "LYZ", "CTSD", "CD14", "FCGR3A", "VCAN",
                       "S100A8", "S100A9", "S100A12", "CDKN1C", "SIGLEC10")
FeaturePlot(Mono, features = L3_Myeloid_marker, ncol = 4)
# L3 Annotation
Mono <- RenameIdents(Mono,
                     '0' = "cMono",'1' = "cMono",'2' = "ncMono",'3' = "ncMono",
                     '4' = "cMono",'5' = "intMono",'6' = "ncMono")
Mono$L3_celltype <- Idents(Mono)
Myeloid$L3_celltype[colnames(Mono)] <- Mono$L3_celltype
saveRDS(Myeloid, "/data/RDS/Annotation/8_Myeloid_annotation.rds")
