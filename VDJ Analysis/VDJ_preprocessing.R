## ------------------------------
library(scRepertoire)
library(data.table)
library(Seurat)
library(ggplot2) 
library(ggpubr) 
library(dplyr)
library(purrr)
library(stringr)
library(mascarade)
library(tidydr)
library(ggbreak)
library(reshape2)
library(rstatix)
library(circlize)
library(RColorBrewer)
## VDJ all packages
## ------------------------------
## ------------------------------
## BGI_VDJseq_preprocessing
## Preprocess TCR data and combine paired receptors
## ------------------------------
## ------------------------------
metadata =read.csv('/ALSM/10vs10_metadata.csv')
combined.TCR <- list.files(
  path = "/ALSM/TCR/data/",
  pattern = "\\.csv$", 
  full.names = TRUE
) %>%
  set_names(~ gsub("^(.*?)_TCR_.*$", "\\1", basename(.))) %>%
  map(~ fread(.x, header = TRUE, stringsAsFactors = FALSE)) %>%
  map(~ mutate(.x, reads = umis)) %>%
  combineTCR(
    samples = names(.),
    removeNA = TRUE,
    removeMulti = FALSE,
    filterMulti = TRUE
  )           
##merge metadata
for (i in seq_along(combined.TCR)) {
  combined.TCR[[i]] <- merge(combined.TCR[[i]], metadata, by = "sample", all.x = TRUE)
}

alsm_object_t<- readRDS('/ALSM/RNA/ALL_T_L3_annotation.rds')
alsm_object_t@meta.data$type <- "pbmc"

sce_combined_tcr <- combineExpression(combined.TCR, 
                         alsm_object_t, 
                         cloneCall = "aa",
                         chain='both', 
                         group.by = "sample",         
                         proportion = FALSE, 
                        cloneSize=c("1"=1, "2"=2, "3"=3, "4"=4, "5"=5,"6-10"=10,">10"=100))
sce_combined_tcr$sample <- sub("_.*$", "", colnames(sce_combined_tcr))
tcr_cells_info <- sce_combined_tcr@meta.data %>%
  dplyr::select(L3_celltype)%>%
  mutate(barcode=rownames(.))
combined_tcr_df <- bind_rows(combined.TCR, .id = "sample")
combined.TCR <- combined_tcr_df %>%
  merge(tcr_cells_info, by ="barcode")%>%
  filter(!is.na(L3_celltype))%>%
 split(.,.$sample)
saveRDS(sce_combined_tcr,"ASLM/TCR/GEX_TCR.rds",compress='bzip2')
saveRDS(combined.TCR,"ASLM/TCR/TCR.rds",compress='bzip2')



## ------------------------------
## bcr
## ------------------------------
combined.BCR <- list.files(
  path = "ALSM/BCR/data/",
  pattern = "\\.csv$", 
  full.names = TRUE
) %>%
  set_names(~ gsub("^(.*?)_BCR_.*$", "\\1", basename(.))) %>%
  map(~ fread(.x, header = TRUE, stringsAsFactors = FALSE)) %>%
  map(~ mutate(.x, reads = umis)) %>%
  combineBCR(
    samples = names(.),
    removeNA = TRUE,
    removeMulti = FALSE,
    filterMulti = TRUE
  )
#merge metadata
for (i in seq_along(combined.BCR)) {
  combined.BCR[[i]] <- merge(combined.BCR[[i]], metadata, by = "sample", all.x = TRUE)
}   
alsm_object_b<- readRDS('ALSM/RNA/PBMC_L1_L2_L3_end_annotation_Bcells.rds')
alsm_object_b@meta.data$type <- "pbmc"
sce_combined_bcr <- combineExpression(combined.BCR, 
                         alsm_object_b, 
                         cloneCall = "strict",
                         chain='both', 
                         group.by = "sample",         
                         proportion = FALSE, 
                        cloneSize=c("1"=1, "2"=2, "3"=3, "4"=4, "5"=5,"6-10"=10,">10"=100))
sce_combined_bcr$sample <- sub("_.*$", "", colnames(sce_combined_bcr))
##fliter non BCR CELLs
bcr_cells_info <- sce_combined_bcr@meta.data %>%
  dplyr::select(L3_celltype)%>%
  mutate(barcode=rownames(.))
combined_bcr_df <- bind_rows(combined.BCR, .id = "sample")
combined.BCR <- combined_bcr_df %>%
  merge(bcr_cells_info, by ="barcode")%>%
  filter(!is.na(L3_celltype))%>%
 split(.,.$sample)
saveRDS(sce_combined_bcr,"ASLM/BCR/GEX_BCR.rds",compress='bzip2')
saveRDS(combined.BCR,"ASLM/TCR/BCR.rds",compress='bzip2')
## ------------------------------
## cd4 tem 0_1
## ------------------------------

cd4tem <- fread('/long_read/CD4Tem_subcelltype.csv')
combined.TCR_CD4_Tem_01 <- lapply(combined.TCR, function(sample_df) {
  result <- merge(sample_df, cd4tem, by = "barcode") %>%
    filter(!is.na(cell_type)) 
  return(result)
})
alsm_object_t_cd4tem <- subset(alsm_object_t,L3_celltype=="CD4 Tem")
sce_combined_tcr_cd4tem <- combineExpression(combined.TCR_CD4_Tem_01, 
                         alsm_object_t_cd4tem, 
                         cloneCall = "aa",
                         chain='TRB', 
                         group.by = "sample",         
                         proportion = FALSE, 
                        cloneSize=c("1"=1, "2"=2, "3"=3, "4"=4, "5"=5,"6-10"=10,">10"=100))
sce_combined_tcr_cd4tem$barcode <- colnames(sce_combined_tcr_cd4tem)
sce_combined_tcr_cd4tem@metadata <- merge(sce_combined_tcr_cd4tem@metadata,cd4tem,by="barcode")
saveRDS(combined.TCR_CD4_Tem_01,"ASLM/TCR/long_read_tcr.rds",compress="bzip2")
saveRDS(sce_combined_tcr_cd4tem,"ASLM/TCR/long_read_tcr_gex.rds",compress="bzip2")












