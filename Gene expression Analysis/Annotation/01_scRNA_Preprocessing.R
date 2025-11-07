library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(R.utils)
.cluster_cols <- c(
    "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", # nolint
    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

load_seurat_group <- function(base_path, gender_vector, group_name = NULL,
                              min_cells = 3, min_features = 200, gene_column = 1) {
  dir_names <- list.files(base_path)
  scRNA_list <- vector("list", length(dir_names))
  for (i in seq_along(dir_names)) {
    # read
    sample_path <- file.path(base_path, dir_names[i])
    counts <- Read10X(data.dir = sample_path, gene.column = gene_column)
    project_name <- strsplit(dir_names[i], "_")[[1]][3]
    # create Seurat object
    scRNA_list[[i]] <- CreateSeuratObject(
      counts = counts,
      project = project_name,
      min.cells = min_cells,
      min.features = min_features
    )
    scRNA_list[[i]]@meta.data$gender <- gender_vector[[i]]
  }
  if (length(scRNA_list) == 1) {
    merged <- scRNA_list[[1]]
  } else {
    merged <- merge(
      x = scRNA_list[[1]],
      y = scRNA_list[2:length(scRNA_list)]
    )
  }
  # add group info
  if (!is.null(group_name)) {
    merged@meta.data$group <- group_name
  }
  return(merged)
}
# 30-40 samples
gender_young <- factor(c(
  "Male", "Female", "Female",
  "Male", "Female", "Male",
  "Female", "Female", "Male"
))
YOUNG_path <- "/data/matrix/YOUNG/"
YOUNG <- load_seurat_group(
  base_path = YOUNG_path,
  gender_vector = gender_young,
  group_name = "YOUNG"
)
head(YOUNG@meta.data)

# 60-70 samples
gender_old <- factor(c(
  "Female", "Male", "Female",
  "Male", "Female", "Male",
  "Female", "Male", "Female", "Male"
))
OLD_path <- "/data/matrix/OLD/"
OLD <- load_seurat_group(
  base_path = OLD_path,
  gender_vector = gender_old,
  group_name = "OLD"
)
head(OLD@meta.data)
# merge young and old
pbmc <- merge(
  x = YOUNG,
  y = OLD
)

#QC and filtering
# mitochondrial percentage
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
# ERCC percentage
pbmc[["percent.ercc"]] <- PercentageFeatureSet(pbmc, pattern = "^ERCC-")
# ribosomal percentage
rp_genes <- grep("^RP[SL][[:digit:]]", x = rownames(pbmc@assays$RNA), value = TRUE)
pbmc <- PercentageFeatureSet(pbmc, features = rp_genes, col.name = "percent.rp")
# non-coding RNA percentage
ncRNA_genes <- grep("^[A-Z][A-Z][0-9]*\\.[0-9]", x = rownames(pbmc@assays$RNA), value = TRUE)
pbmc <- PercentageFeatureSet(pbmc, features = ncRNA_genes, col.name = "percent.ncRNA")
# lincRNA percentage
lincRNA_genes <- grep("(^LOC|LINC)[1-9]*", x = rownames(pbmc@assays$RNA), value = TRUE)
pbmc <- PercentageFeatureSet(pbmc, features = lincRNA_genes, col.name = "percent.LOC")
# hemoglobin percentage
HB_gene_list <- c(
  "HBA1","HBA2","HBB","HBD","HBE1",
  "HBG1","HBG2","HBM","HBQ1","HBZ"
)
HB_matches <- match(HB_gene_list, rownames(pbmc@assays$RNA))
HB_genes <- rownames(pbmc@assays$RNA)[HB_matches]
HB_genes <- HB_genes[!is.na(HB_genes)] # remove NA values
pbmc[["percent.HB"]] <- PercentageFeatureSet(
  pbmc,
  features = HB_genes
)

# Violin plot function for QC metrics
qc_features <- c(
  "nFeature_RNA", "nCount_RNA", 
  "percent.mt", "percent.rp", 
  "percent.ncRNA", "percent.HB"
)
plot_qc_violin <- function(seurat_obj, file_name, width = 20, height = 10) {
  p <- VlnPlot(
    seurat_obj,
    features = qc_features,
    ncol = 3,
    pt.size = 0,
    group.by = "orig.ident"
  )
  ggsave(file_name, plot = p, width = width, height = height)
}
# before filtering
plot_qc_violin(pbmc, "/data/Figure/QC/Raw.pdf")
saveRDS(pbmc, file = "/data/RDS/scRNA_merged_raw.rds")
# filtering
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 400 & nFeature_RNA < 5000 &
    percent.mt < 10 &
    percent.HB < 1.5
)
# after filtering
pbmc
#Filtered out 5579 cells. Original number of cells: 121830, after filtering:116251
plot_qc_violin(pbmc, "/data/Figure/QC/QC_Filter.pdf")
saveRDS(pbmc, file = '/data/RDS/FilterData_YOUNG_OLD.rds')

# Standard preprocessing workflow
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, vars.to.regress = c("percent.mt"))
pbmc <- RunPCA(pbmc)
ElbowPlot(pbmc, ndims = 50)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:12)
pbmc <- FindClusters(pbmc, resolution = 0.7)
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:12, seed.use = 123)
p_umap <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
# UMAP plot after QC
ggsave(filename = "/data/Figure/QC/umap_after_QC.pdf", plot = p_umap, width = 8, height = 8)

# doublet deletion using DoubletFinder
library(DoubletFinder)
pbmc_list <- SplitObject(pbmc, split.by = "orig.ident")
pkchoose_list <- list()
for (i in seq_along(pbmc_list)) {
  seu_temp <- pbmc_list[[i]]
  sweep_res_list <- paramSweep_v3(seu_temp, PCs = 1:12, sct = TRUE)
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  # PK selection
  sweep_stats_pk <- find.pK(sweep_stats)
  pK <- as.numeric(as.character(sweep_stats_pk$pK))
  BCmetric <- sweep_stats_pk$BCmetric
  pK_choose <- pK[which(BCmetric %in% max(BCmetric))]
  # anticipated number of doublets
  homotypic_prop <- modelHomotypic(seu_temp$seurat_clusters)
  doublet_rate <- ncol(seu_temp) * 8 * 1e-6
  nExp_poi <- round(doublet_rate * nrow(seu_temp@meta.data))
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))
  # DoubletFinder
  seu_temp <- doubletFinder_v3(
    seu_temp,
    PCs = 1:12,
    sct = TRUE,
    pN = 0.25,
    pK = as.numeric(pK_choose),
    nExp = nExp_poi_adj,
    reuse.pANN = FALSE
  )
  pbmc_list[[i]] <- seu_temp
  pkchoose_list[[i]] <- pK_choose
  message(sprintf("Doublet detection completed for sample %d", i))
}
# colname name standardization
for (i in seq_along(pbmc_list)) {
  sc_obj <- pbmc_list[[i]]
  doublet_col <- grep("DF.classifications", colnames(sc_obj@meta.data), value = TRUE)
  if (length(doublet_col) > 0) {
    colnames(sc_obj@meta.data)[colnames(sc_obj@meta.data) == doublet_col] <- "DF.classifications"
  }
  pbmc_list[[i]] <- sc_obj
}
# delete pANN columns
for (i in seq_along(pbmc_list)) {
  sc_obj <- pbmc_list[[i]]
  pann_cols <- grep("pANN_", colnames(sc_obj@meta.data), value = TRUE)
  if (length(pann_cols) > 0) {
    sc_obj@meta.data <- sc_obj@meta.data[, !colnames(sc_obj@meta.data) %in% pann_cols]
  } else {
    message(sprintf("No pANN columns found in sample %d", i))
  }
  pbmc_list[[i]] <- sc_obj
}

# merge back
pbmc <- merge(pbmc_list[[1]], pbmc_list[2:length(pbmc_list)])
saveRDS(pbmc, file = "/data/RDS/pbmc_with_doublets.rds")
# subset singlets
pbmc_singlet <- subset(pbmc, subset = DF.classifications == "Singlet")
saveRDS(pbmc_singlet, file = "/data/RDS/pbmc_remove_doublets.rds")
#Doublets Filtered out 4932 cells. Original number of cells: 116254, after filtering: 111319