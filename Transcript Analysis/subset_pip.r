# ------------------------------------------------------------------------------
# This script subsets a Seurat object to a specific cell type and performs isoform-level differential usage / isoform switching analysis
# using IsoformSwitchAnalyzeR with the satuRn statistical test. It performs basic per-sample integration and clustering of the subset,
# aggregates isoform counts per sample/cluster, constructs a design matrix, runs isoform-level filtering and tests, and writes out a table
# of significant isoform switches.
# ------------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(ggplot2)
library(IsoformSwitchAnalyzeR)
library(Matrix)
library(rtracklayer)

# Load full Seurat object and subset to target cell type
lr_raw <- readRDS('./data_lr.rds')
lr_subset <- subset(lr_raw, subset = L3_celltype == 'CD4 Tem')

# Split RNA assay by sample for integration
lr_subset[["RNA"]] <- split(lr_subset[["RNA"]], f = lr_subset$sample)

# Basic preprocessing
lr_subset <- NormalizeData(lr_subset)
lr_subset <- FindVariableFeatures(lr_subset)
lr_subset <- ScaleData(lr_subset)
lr_subset <- RunPCA(lr_subset)

# Integrate samples using CCA
lr_subset <- IntegrateLayers(
  object = lr_subset, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = TRUE, k.weight = 50
)
lr_subset[["RNA"]] <- JoinLayers(lr_subset[["RNA"]])

# Clustering and visualization
lr_subset <- FindNeighbors(lr_subset, reduction = "integrated.cca", dims = 1:30)
lr_subset <- FindClusters(lr_subset, resolution = 0.1, algorithm = 1)
lr_subset <- RunUMAP(lr_subset, dims = 1:6, n.neighbors = 30L, min.dist = 0.1, reduction = "integrated.cca")
saveRDS(lr_subset, file = "./CD4Tem.rds")

# Aggregate isoform counts per sample-cluster
lr_subset$sample_ct <- paste0(lr_subset$sample, "_", lr_subset$RNA_snn_res.0.1)
pb_list <- AggregateExpression(
  object = lr_subset, assays = "RNA", group.by = "sample_ct",
  slot = "counts", return.seurat = FALSE
)
iso_pb <- pb_list$RNA

# Build design and comparison matrices
cols <- colnames(iso_pb)
group <- ifelse(grepl("-0", cols), "0", "1")
design <- data.frame(sampleID = cols, condition = group, stringsAsFactors = FALSE)
comparisons <- data.frame(condition_1 = '0', condition_2 = '1', stringsAsFactors = FALSE)

# Import data into IsoformSwitchAnalyzeR
sar1 <- importRdata(
  isoformCountMatrix = as.matrix(iso_pb),
  designMatrix = design,
  comparisonsToMake = comparisons,
  isoformExonAnnoation = "./CD4Tem.gtf",
  isoformNtFasta = "./OUT.novel.transcript_models_corrected.fasta",
  detectUnwantedEffects = FALSE
)

# Predict ORFs and apply expression filtering
sar1 <- analyzeORF(sar1)
sar1 <- preFilter(
  switchAnalyzeRlist = sar1,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE,
  reduceToSwitchingGenes = FALSE
)

# Perform differential isoform usage test using satuRn
sar2 <- isoformSwitchTestSatuRn(
  switchAnalyzeRlist = sar1,
  reduceToSwitchingGenes = FALSE
)

# Extract and filter significant isoform switches
sw_tbl <- as.data.frame(sar2$isoformFeatures)
sig_tbl <- sw_tbl %>%
  filter(isoform_switch_q_value < 1, abs(dIF) >= 0.1) %>%
  arrange(isoform_switch_q_value, desc(abs(dIF)))

# Save result table
write.table(sig_tbl, file = "CD4Tem.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
