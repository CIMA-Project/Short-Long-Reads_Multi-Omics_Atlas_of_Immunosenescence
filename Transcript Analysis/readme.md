# Single-cell Full-length Transcriptome Analysis Pipeline

## Overview

This pipeline performs comprehensive analysis of single-cell full-length transcriptome data, including preprocessing, metacell assignment, differential transcript usage, and differential splicing analysis. The workflow consists of five main components:

1. **Data Preprocessing** (`preprocess.py`)
2. **Metacell Assignment** (`metacell.py`)
3. **Differential Transcript Usage Analysis** (`Differential_Transcript_Usage.r`)
4. **Subset Analysis** (`subset_pip.r`)
5. **Differential Splicing Analysis** (`SUPPA.sh`)

---

## Data Preprocessing

### Script: `preprocess.py`


```bash
python preprocess.py
```

### Input Files Required in Current Directory:

* `./OUT.discovered_transcript_grouped_counts.barcodes.tsv` - Cell barcodes file
* `./OUT.discovered_transcript_grouped_counts.features.tsv` - Feature (isoform) names file
* `./OUT.discovered_transcript_grouped_counts.matrix.mtx` - Count matrix in Matrix Market format
* `./OUT.novel.transcript_models_classification.txt` - SQANTI3 classification results
* `./OUT.novel.transcript_models_corrected.gtf` - GTF file with corrected transcript models
* `./data_sr.h5ad` - Short-read single-cell data for cell type annotation transfer

### Output:

* `./data_lr.h5ad` - Preprocessed AnnData object containing long-read single-cell data with cell type annotations

### Key Processing Steps:

1. Reads the count matrix and associated metadata
2. Adds SQANTI3 annotations for transcript classification
3. Performs basic quality control filtering (minimum 3 cells per gene, 5 genes per cell)
4. Maps gene names from GTF file to transcript IDs
5. Filters transcripts based on structural categories
6. Calculates QC metrics (mitochondrial, ribosomal, hemoglobin, ncRNA content)
7. Transfers cell type annotations from short-read data
8. Filters genes expressed in at least 3 samples

---

## Metacell Assignment

### Script: `metacell.py`

```bash
python metacell.py
```

### Input:

* `./data_lr.h5ad` - Preprocessed AnnData object from Data Preprocessing

### Output:

* `./metacell_assignment.csv` - CSV file mapping cells to metacells
* `./metacell.h5ad` - AnnData object with aggregated metacell expression
* `./metacell/` - Directory containing intermediate files

### Key Processing Steps:

1. Creates metacells using SEACells algorithm separately for each sample and cell type
2. Applies TF-IDF transformation and normalization to expression data
3. Performs dimensionality reduction using Truncated SVD
4. Assigns cells to metacells based on transcriptional similarity
5. Aggregates expression profiles at metacell level
6. Handles edge cases (small cell populations, convergence issues)

---

## Differential Transcript Usage Analysis

### Script: `Differential_Transcript_Usage.r`


```r
Rscript Differential_Transcript_Usage.r
```

### Input:

* `./data_lr.rds` - Seurat object containing isoform-level RNA data
* `./OUT.novel.transcript_models_common.gtf` - GTF annotation file
* `./OUT.novel.transcript_models_corrected.fasta` - Transcript nucleotide sequences
* External annotation files for IsoformSwitchAnalyzeR (PFAM, CPC2, IUPred2A, SignalP, DeepTMHMM, DeepLoc2 results)

### Output:

* `./sar1.rds` - IsoformSwitchAnalyzeR object after ORF prediction
* `./sar2.rds` - IsoformSwitchAnalyzeR object after consequence analysis
* `./sar3.rds` - Final IsoformSwitchAnalyzeR object with all annotations
* `./IsoSA_part1/` - Directory containing initial analysis results
* `./IsoSA_part2/` - Directory containing detailed analysis results

### Key Processing Steps:

1. Aggregates isoform expression by sample and cell type
2. Builds design matrix comparing 30-40 age group vs 60-70 age group samples for each cell type
3. Imports data into IsoformSwitchAnalyzeR framework
4. Predicts open reading frames (ORFs)
5. Analyzes isoform switch consequences (TSS, TTS, UTR changes, ORF length, NMD status, etc.)
6. Incorporates external functional annotations

---

## Subset Analysis

### Key Processing Steps:

1. Subsets data to CD4 Tem cell type
2. Performs sample integration using CCA
3. Clusters cells and runs UMAP visualization
4. Aggregates isoform counts by sample and cluster
5. Runs differential usage analysis using satuRn statistical test
6. Filters significant isoform switches (|dIF| ≥ 0.1)

---

## Differential Splicing Analysis

### Script: `SUPPA.sh`


```bash
bash SUPPA.sh
```

### Input Files Required in Current Directory:

* `./model.ioe` - Event annotation file in ioe format
* `./suppa_expression.txt` - Expression matrix (TPM values) with samples as columns

### Output:

* `./events.psi` - PSI (Percent Spliced In) values per event
* `./60_70.psi` - PSI values for the 60-70 age group (columns 1-11)
* `./60_70.tpm` - TPM values for the 60-70 age group (columns 1-11)
* `./30_40.psi` - PSI values for the 30-40 age group (columns 1,12-20)
* `./30_40.tpm` - TPM values for the 30-40 age group (columns 1,12-20)
* `./dpsi.dpsi` - Differential splicing results
* `./output2.txt` - Filtered differential splicing results (non-nan values)

### Key Processing Steps:

1. Calculates PSI values for each splicing event
2. Creates subset files for two age groups (30-40 and 60-70)
3. Runs differential splicing analysis using empirical method
4. Filters out results with 'nan' values while preserving header

---

## Required External Tools and Dependencies

### Python Packages:

* anndata, scanpy, scCyclone, SEACells, numpy, pandas, scipy, sklearn, tqdm, psutil

### R Packages:

* Seurat, ggplot2, IsoformSwitchAnalyzeR, Matrix, rtracklayer, dplyr

### External Bioinformatics Tools:

* SUPPA2 (for splicing analysis)
* IsoformSwitchAnalyzeR dependencies:
  * PFAM for domain annotation
  * CPC2 for coding potential prediction
  * IUPred2A for disordered regions
  * SignalP for signal peptide prediction
  * DeepTMHMM for transmembrane helices
  * DeepLoc2 for subcellular localization
