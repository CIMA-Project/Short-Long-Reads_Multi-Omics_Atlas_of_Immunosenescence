# Full-length Transcript TSS and PAS Site Analysis Pipeline

## Overview

This pipeline processes long-read sequencing data to identify Transcription Start Sites (TSS) and Polyadenylation Sites (PAS), followed by UTR length analysis. The workflow consists of five main steps:

1. **Extract TSS and PAS information from BAM files** (`01.Get_read_TSS_PAS.py`)
2. **Convert read information to BED format** (`02.From_readsinfo_GetBed.py`)
3. **Filter TSS and PAS sites** (`03.Filter_TSS_site.py` and `03.Filter_PAS_site.py`)
4. **Filter annotated transcripts** (`04.Filter_annot_transcript.py`)
5. **Calculate UTR5 and UTR3 lengths and generate matrices** (`05.Get_read_utr5_utr3_length.py`)

---

## Step 1: Extract TSS and PAS Information from BAM Files

### Script: `01.Get_read_TSS_PAS.py`


```bash
python 01.Get_read_TSS_PAS.py -i <input_bam> -o <output_prefix> -r <reference_genome>
```

### Input:

* `$input_bam` - Input BAM file (aligned to reference genome)
* `$reference_genome` - Reference genome FASTA file

### Output:

* `$output_prefix_TSS_PAS.txt` - Text file containing TSS and PAS site information for each read

### Required Files:

* BAM file (must contain alignment information)
* Reference genome FASTA file
* BAM file must contain NB and GN tags for proper processing

---

## Step 2: Convert Read Information to BED Format

### Script: `02.From_readsinfo_GetBed.py`


```bash
python 02.From_readsinfo_GetBed.py -r <read_assignments> -i <info_files> -o <output_dir>
```

### Input:

* `$read_assignments` - Gzipped read assignments file (e.g., `read_assignments.txt.gz`)
* `$info_files` - One or more gzipped info files (output from Step 1)
* `$output_dir` - Output directory (optional, defaults to current directory)

### Output:

* `$prefix.TSS.bed` - BED file for Transcription Start Sites
* `$prefix.TES.bed` - BED file for Transcription End Sites (Polyadenylation Sites)

### Required Files:

* `read_assignments.txt.gz` (contains read-to-transcript mapping)
* Info files generated in Step 1

---

## Step 3: Filter TSS and PAS Sites

### Scripts: `03.Filter_TSS_site.py` and `03.Filter_PAS_site.py`


```bash
# Filter TSS sites
python 03.Filter_TSS_site.py <TSS_FilterInfo_file> <output_file> [num_processes]

# Filter PAS sites
python 03.Filter_PAS_site.py <TES_FilterInfo_file> <output_file> [num_processes]
```

### Input:

* `$TSS_FilterInfo_file` - TSS information file for filtering (BED format from Step 2)
* `$TES_FilterInfo_file` - TES/PAS information file for filtering (BED format from Step 2)
* `$num_processes` - Optional parameter for number of parallel processes (default: CPU core count)

### Output:

* `$output_file` - Filtered result file

### Required Files:

* TSS.bed or TES.bed format files from Step 2
* Scripts support multi-processing for improved performance

---

## Step 4: Filter Annotated Transcripts

### Script: `04.Filter_annot_transcript.py`


```bash
python 04.Filter_annot_transcript.py <long_file> <annot_file> <output_file> [num_processes]
```

### Input:

* `$long_file` - File containing transcript-to-gene mapping (format: transcript_name\tgene_name\tlength)
* `$annot_file` - Annotation information file to be filtered
* `$output_file` - Output result file
* `$num_processes` - Optional parameter for number of parallel processes (default: CPU core count)

### Output:

* `$output_file` - Filtered annotation information file

### Required Files:

* Transcript-gene mapping file (containing transcript name, gene name, and transcript length)
* Annotation file to be filtered

---

## Step 5: Calculate UTR5 and UTR3 Lengths and Generate Matrices

### Script: `05.Get_read_utr5_utr3_length.py`


```bash
python 05.Get_read_utr5_utr3_length.py -utr5 <utr5_file> -utr3 <utr3_file> -o <output_dir>
```

### Input:

* `$utr5_file` - UTR5 distance measurements file
* `$utr3_file` - UTR3 distance measurements file
* `$output_dir` - Output directory (optional, defaults to "utr_results")

### Output:

* `gene_utr5_matrix.txt` - Gene-level UTR5 length and relative value matrix
* `transcript_utr5_matrix.txt` - Transcript-level UTR5 length and relative value matrix
* `cell_utr5_relative.txt` - Cell-level UTR5 relative value statistics
* `gene_utr3_matrix.txt` - Gene-level UTR3 length and relative value matrix
* `transcript_utr3_matrix.txt` - Transcript-level UTR3 length and relative value matrix
* `cell_utr3_relative.txt` - Cell-level UTR3 relative value statistics
* `reads_utr5_utr3_length.txt` - UTR5 and UTR3 information for each read

### Required Files:

* UTR5 distance measurements file
* UTR3 distance measurements file

---

## Complete Pipeline Example


```bash
# Step 1: Extract TSS and PAS information
python 01.Get_read_TSS_PAS.py -i sample.bam -o sample -r reference.fa

# Step 2: Convert to BED format
python 02.From_readsinfo_GetBed.py -r read_assignments.txt.gz -i sample.txt.gz -o bed_files/

# Step 3: Filter sites
python 03.Filter_TSS_site.py bed_files/sample.TSS.bed filtered_TSS.bed 8
python 03.Filter_PAS_site.py bed_files/sample.TES.bed filtered_TES.bed 8

# Step 4: Filter based on annotation
python 04.Filter_annot_transcript.py gene_transcript_map.txt filtered_TSS.bed final_TSS.bed 8
python 04.Filter_annot_transcript.py gene_transcript_map.txt filtered_TES.bed final_TES.bed 8

# Step 5: Calculate UTR lengths
python 05.Get_read_utr5_utr3_length.py -utr5 utr5_distances.txt -utr3 utr3_distances.txt -o utr_results/
```
