#gene and transcript expression Dotplot
#short_read
library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize) 
library(ggplotify)
cell <- readRDS("/data/RDS/Annotation/10_PBMC_annotation_L3.rds")
options(repr.plot.height = 7, repr.plot.width = 12)
all_marker <- c("MS4A1", "IGHD", "CD79A", "MZB1", "MKI67", "CD34",
  "GP9", "CD14", "FCGR3A", "LYZ", "IDO1", "CD1C", "LILRA4",
  "CD3D", "CD4", "FOXP3", "CCR7", "IL7R",
  "CD8A", "SLC4A10",
  "NCR1", "GZMA", "GZMH", "NKG7", "CD160", "GZMK", "XCL1"
)
genes_to_check <- lapply(all_marker, function(x) str_to_upper(x))
#Figure1G: Dotplot of gene expression by cell type
p <- DotPlot(pbmc, features = all_marker,
  assay = "RNA", cols = "RdBu", col.min = -1, col.max = 1
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, hjust = 1)
  ) +
  scale_y_discrete(limits = order) +
  guides(size = guide_legend("Percent Expression")) +
  labs(x = NULL, y = NULL)
print(p)

ggsave("/data/Figure/PBMC_L3_gene_dotplot.pdf", 
       plot = p, width = 12, height = 7)

#long_read
pbmc <- readRDS ("/data/RDS/1010umap.rds")

gene_order <- all_marker
L3_features <- pbmc@assays$RNA@meta.features %>%
  filter(gene_name %in% gene_order) %>%
  mutate(gene_name = factor(gene_name, levels = gene_order)) %>%
  arrange(gene_name) %>% 
  rownames()
L3_gene_names <- pbmc@assays$RNA@meta.features[L3_features, "gene_name"]
# gene and transcript mapping table
L3_gene_transcript_map <- data.frame(
  Gene_Name = L3_gene_names,
  Transcript_ID = L3_features
)

DefaultAssay(pbmc) <- "RNA"
gene_cell_exp <- AverageExpression(pbmc,
                                   features = L3_features,
                                   group.by = 'L3_celltype',
                                   slot = "data"
)
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
marker_exp <- scale(t(gene_cell_exp), scale = T, center = T)
marker_exp

data <- L3_gene_transcript_map
transcript_to_gene <- setNames(data$Gene_Name, data$Transcript_ID) 
Transcript_groups <- factor(transcript_to_gene[colnames(marker_exp)], 
                      levels = gene_order)
Transcript_groups
# Column annotation for heatmap
col_anno <- HeatmapAnnotation(
  GeneGroup = anno_block(
    labels = levels(Transcript_groups),
    labels_gp = gpar(col = "black", fontsize = 9, fontface = "bold"),
    labels_rot = 90,
    labels_just = "center",
    gp = gpar(
      col = "gray",
      lwd = 1
    ),
    width = unit(1.2, "cm")
  ),
  show_annotation_name = FALSE
)
# Color function for heatmap
col_fun <- colorRamp2(
  breaks = c(-1, 0, 2),
  colors = c('#0B75B3', 'white', '#C31E1F')
)

#Figure1H: Heatmap of Transcript expression by cell type
p3 <- as.ggplot(Heatmap(marker_exp,
        name = "Expression",
        col = col_fun,  
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_column_slices = FALSE,
        column_split = Transcript_groups,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_title = "Cell Type",
        row_order = rev(order), 
        column_title = NULL,
        top_annotation = col_anno,
        column_names_gp = gpar(fontsize = 8),
        column_names_rot = 90,
        column_names_centered = F,
        column_gap = unit(1, "mm"),
        row_names_gp = gpar(fontsize = 10),
        row_names_side = "left", 
        heatmap_legend_param = list(
          title = "Expression Level",
          at = c(-1, 0, 2),  # 确保包含0值
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8)
        ),
        #rect_gp = gpar(col = "gray", lwd = 0.5),
        rect_gp = gpar(col = NA, lwd = 0),
        border = TRUE,
        border_gp = gpar(col = "darkgray", lwd = 1))
)
p3
ggsave("/data/Figure/PBMC_L3_heatmap_order_celltype_gene.pdf", plot= p3,width = 22, height = 10)