# ------------------------------------------------------------------------------
# This script performs Differential Transcript Usage (DTU) analysis 
# on single-cell RNA-seq data using Seurat and IsoformSwitchAnalyzeR.
# It aggregates isoform-level counts per cell type and sample, 
# constructs a design matrix, and identifies isoform switches between groups.
# ------------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(IsoformSwitchAnalyzeR)
library(Matrix)
library(rtracklayer)

# Load single-cell Seurat object containing isoform-level RNA data
lr <- readRDS("./data_lr.rds")

# Create a combined label for each sample and cell type
lr$sample_ct <- paste0(lr$sample, "_", lr$L3_celltype)

# Aggregate expression counts by sample-celltype combinations
# Returns a list containing summed expression matrices for each assay
pb_list <- AggregateExpression(
  object         = lr,
  assays         = "RNA",
  group.by       = "sample_ct",
  slot           = "counts",
  return.seurat  = FALSE
)

# Extract RNA assay matrix (genes/isoforms as rows, aggregated groups as columns)
iso_pb <- pb_list$RNA

# Filter out isoforms with zero total expression across all samples
total_expr <- rowSums(iso_pb)
keep <- total_expr > 0
iso_pb <- iso_pb[keep, ]

# Extract column names representing each aggregated group
cols <- colnames(iso_pb)

# Define experimental groups
group <- ifelse(grepl("^Y", cols), "Y", "O")

# Extract cell type labels
celltype <- sub("^[^-]+-", "", cols)

# Build design matrix describing sample IDs and corresponding conditions
design <- data.frame(
  sampleID  = cols,
  condition = paste0(group, "-", celltype),
  stringsAsFactors = FALSE
)

# Generate condition and cell type lists
conds <- unique(design$condition)
celltypes <- unique(sub("^[^-]+-", "", conds))

# Create comparison pairs between old and young samples for each cell type
comparisons <- data.frame(
  condition_1 = paste0("O-", celltypes),
  condition_2 = paste0("Y-", celltypes),
  stringsAsFactors = FALSE
)

# Import data into IsoformSwitchAnalyzeR
sar1 <- importRdata(
  isoformCountMatrix    = as.matrix(iso_pb),
  designMatrix          = design,
  comparisonsToMake     = comparisons,
  isoformExonAnnoation  = "./OUT.novel.transcript_models_common.gtf",
  isoformNtFasta        = "./OUT.novel.transcript_models_corrected.fasta",
  detectUnwantedEffects = FALSE
)

# Predict open reading frames (ORFs) for imported isoforms
sar1 <- analyzeORF(sar1)
saveRDS(sar1, file = "./sar1.rds")

sar2 <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist     = sar1,
  pathToOutput           = "./IsoSA_part1",
  prepareForWebServers   = FALSE,
  outputSequences        = TRUE
)

# Assess biological consequences of identified isoform switches
# such as changes in UTRs, ORFs, exon structure, etc.
sar2 <- analyzeSwitchConsequences(
  sar2, 
  onlySigIsoforms       = TRUE, 
  dIFcutoff             = 0, 
  consequencesToAnalyze = c(
    'tss','tts','last_exon','isoform_length','exon_number',
    'intron_structure','ORF_length','5_utr_seq_similarity','5_utr_length',
    '3_utr_seq_similarity','3_utr_length','ORF_seq_similarity','NMD_status'
  )
)
saveRDS(sar2, file = "./sar2.rds")

sar3 <- isoformSwitchAnalysisPart2(
  switchAnalyzeRlist           = sar2, 
  n                            = 10,  # For plotting, limits to top 10 switches
  removeNoncodinORFs           = FALSE,
  pathToOutput                 = "./IsoSA_part2",
  pathToPFAMresultFile         = "./IsoSA_part2/pfam_results.txt",
  pathToCPC2resultFile         = "./IsoSA_part2/cpc2_results.txt",
  pathToIUPred2AresultFile     = "./IsoSA_part2/iupred2a_results.txt",
  pathToSignalPresultFile      = "./IsoSA_part2/signalP_results.txt", 
  pathToDeepTMHMMresultFile    = "./IsoSA_part2/DeepTMHMM.gff3",
  pathToDeepLoc2resultFile     = "./IsoSA_part2/deeploc2_results.csv"
)

# Analyze additional isoform-level consequences incorporating new annotations
sar3 <- analyzeSwitchConsequences(
  sar3, 
  onlySigIsoforms       = TRUE, 
  dIFcutoff             = 0, 
  consequencesToAnalyze = c(
    'tss','tts','last_exon','isoform_length','exon_number','intron_structure',
    'ORF_length','5_utr_seq_similarity','5_utr_length','3_utr_seq_similarity',
    '3_utr_length','coding_potential','ORF_seq_similarity','NMD_status',
    'domains_identified','signal_peptide_identified'
  )
)

# Save final IsoformSwitchAnalyzeR object
saveRDS(sar3, file = "./sar3.rds")
