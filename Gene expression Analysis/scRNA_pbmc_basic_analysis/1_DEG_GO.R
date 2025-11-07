library(Seurat)
library(tibble)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

pbmc <- readRDS("/data/RDS/Annotation/10_PBMC_annotation_L3.rds")
exclude_patterns <- c(
  "^RPL",
  "^RPS",
  "\\d+\\.\\d+",
  "^RNU.*P$",
  "^RN7SL.*P$",
  "^LINC",
  "^MT"
)
exclude_pattern <- paste(exclude_patterns, collapse = "|")
genes_to_exclude <- rownames(pbmc)[grepl(exclude_pattern, rownames(pbmc))]
pbmc <- pbmc[!rownames(pbmc) %in% genes_to_exclude, ]
write_csv(as_tibble(genes_to_exclude,
                    .name_repair = "Gene"), "/data/Files/DEG/DEG_gene_to_exclude.csv")

# Differential Expression Analysis by Cell Type
outpath <- "/data/Files/DEG/"
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

all_celltypes <- as.character(unique(pbmc@meta.data$L3_celltype))
summary_table <- tibble(CellType = character(), 
                        UP_Count = integer(), DOWN_Count = integer())
for (celltype in all_celltypes) {
  cell_indices <- which(pbmc@meta.data$L3_celltype == celltype)
  if (length(cell_indices) == 0) {
    cat("No cells found for celltype:", celltype, "Skipping...\n")
    next
  }
  pbmc_subset <- pbmc[, cell_indices]
  Idents(pbmc_subset) <- pbmc_subset$group
  Markers <- FindAllMarkers(pbmc_subset, 
                            only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0)
  Markers$pvalue_L <- -log10(Markers$p_val_adj + 10^-200)
  Markers$Class <- ifelse(
    Markers$avg_log2FC > 0.15 & Markers$p_val_adj < 0.05 &
      Markers$pct.1 > 0.1 & Markers$pct.2 > 0.1 & 
      Markers$cluster == 'YOUNG', "Down",
    ifelse(
           Markers$avg_log2FC > 0.15 & Markers$p_val_adj < 0.05 &
           Markers$pct.1 > 0.1 & Markers$pct.2 > 0.1 & 
           Markers$cluster == 'OLD', "Up", "None")
    )
  Markers$avg_log2FC <- ifelse(Markers$cluster == 'YOUNG',
                               -Markers$avg_log2FC, Markers$avg_log2FC)
  up_count <- sum(Markers$Class == "Up", na.rm = TRUE)
  down_count <- sum(Markers$Class == "Down", na.rm = TRUE)
  summary_table <- bind_rows(summary_table,tibble(CellType = celltype, UP_Count = up_count, DOWN_Count = down_count))

  write_csv(Markers, file.path(outpath, paste0(celltype, "_diff_genes.csv")))
}
write_csv(summary_table, file.path(outpath, "DEG_summary_by_L3_celltype.csv"))

# Combine and Filter DEG Results
directory <- "/data/Files/DEG/"
csv_files <- list.files(directory, pattern = "\\.csv$", full.names = TRUE)
all_df <- lapply(csv_files, function(file) {
  celltype <- sub("_.*", "", basename(file))
  df <- read_csv(file, show_col_types = FALSE)
  df$Celltype <- celltype
  df
})
merged_df <- bind_rows(all_df)
write_csv(merged_df, file.path(directory, "combined_data.csv"))
filtered_df <- merged_df %>%
  filter(Class %in% c("Up", "Down"))
write_csv(filtered_df, file.path(directory, "combined_data_filtered.csv"))

#GO Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(ggplot2)
library(readr)
df_sig <- read_csv("/data/Files/DEG/combined_data_filtered.csv") %>%
  filter(Class == "Down")
group_df <- df_sig %>% select(gene, Celltype) %>% rename(group = Celltype)
gene_id_df <- bitr(group_df$gene,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

data <- left_join(gene_id_df, group_df, by = c("SYMBOL" = "gene"))
data_GO <- compareCluster(ENTREZID ~ group,
                          data          = data,
                          fun           = "enrichGO",
                          OrgDb         = org.Hs.eg.db,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
data_GO_sim <- simplify(data_GO, cutoff = 0.7, by = "p.adjust", select_fun = min)
selected_terms <- readLines("/data/Files/GO/selected_GO_term - down.txt") %>% str_trim()
df_go <- as.data.frame(data_GO_sim)
plot_data <- df_go %>% filter(Description %in% selected_terms)
write_csv(plot_data, "/data/Files/GO/plot_filtered_GO_terms_DOWN.csv")

#Figure 2A summary of DEG numbers
deg_summary_path <- "/data/Files/DEG/DEG_summary_by_L3_celltype.csv"
data <- read_csv(deg_summary_path)
cell_colors <- c(
  `HSPC` = "#FAD7D7", `pDC` = "#FBC5AB", `MK` = "#F6BABB",
  `cDC1` = "#EF8787", `intMono` = "#DE6A69", `cDC2` = "#F29667",
  `cMono` = "#CF4457", `ncMono` = "#EA5958", `Plasmablast` = "#E6E6F0",
  `Plasma` = "#DBC6D9", `Swit. Bm` = "#D48CB1", `Unswit. Bm` = "#8085C0",
  `Naïve B` = "#D4A7C4", `NKbright` = "#F07D50", `XCL+NKdim` = "#FEC928",
  `Terminal NK` = "#F4EB1B", `NKdim` = "#FBAA34", `CD4 Regulatory` = "#89CAEA",
  `CD4 CTL` = "#BBE6FA", `CD4 Tcm` = "#0B75B3", `CD4 Tem` = "#4596CD",
  `CD4 Naïve` = "#015696", `Cycling T` = "#77C3BC", `gdT` = "#24858E",
  `MAIT` = "#47A183", `CD8 Tem` = "#2BB17E", `CD8 Tcm` = "#50C56A",
  `CD8 CTL` = "#C3E023", `CD8 Naïve` = "#8DCA9D", `NKT` = "#87D649"
)
cell_order <- names(cell_colors)
data_long <- data %>%
  pivot_longer(
    cols = c(UP_Count, DOWN_Count),
    names_to = "DEG_type",
    values_to = "Count"
  ) %>%
  mutate(
    Cell_type = factor(CellType, levels = rev(cell_order)),
    DEG_type = factor(DEG_type, levels = c("UP_Count", "DOWN_Count"))
  )
# plotting
p <- ggplot(data_long, aes(x = Cell_type, y = Count)) +
  geom_bar(
    aes(alpha = DEG_type, fill = Cell_type),
    stat = "identity", 
    position = position_stack(reverse = TRUE), 
    width = 0.9,
    show.legend = FALSE
  ) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = NA, color = "black", linewidth = 0.1) +
  annotate("segment",
           y = 0, yend = -max(data_long$Count) * 0.02,
           x = unique(data_long$Cell_type),
           xend = unique(data_long$Cell_type),
           color = "black", linewidth = 0.5) +
  scale_alpha_manual(values = c("UP_Count" = 1, "DOWN_Count" = 0.5)) +
  scale_fill_manual(values = cell_colors) +
  labs(
    x = NULL,
    y = "Number of DEGs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(color = "black",
                               size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 20, 15, 15),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  coord_cartesian(clip = "off")

output_path <- "/data/work/Figure/F2A_cell_type_l3_DEG.pdf"
ggsave(output_path, plot = p, width = 12, height = 6)

#Figure 2B GO enrichment terms
plot_data <- read.csv("/data/Files/GO/plot_filtered_GO_terms_DOWN.csv")
plot_data_sorted <- plot_data %>%
  arrange(Cluster, p.adjust) %>%
  mutate(
    Cluster      = factor(Cluster, levels = unique(Cluster)),  
    Description  = str_wrap(Description, width = 30),
    Description  = factor(Description, levels = rev(unique(Description)))
  )
p2 <- ggplot(plot_data_sorted, aes(x = group, y = Description)) + 
  geom_point(aes(size = GeneRatio, fill = p.adjust),
             alpha = 0.8, shape = 21, stroke = 0.5) +
  scale_fill_gradient(
    low   = "#DC050C", high  = "#1965B0",
    name  = "Adjusted p-value",
    guide = guide_colorbar(barwidth = 15, barheight = 1,
                           title.position = "top", title.hjust = 0.5)) +
  scale_size_continuous(
    name  = "Gene Ratio",
    guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  ggtitle("Downregulated GO Terms: Age 30–40 Cohort") +
  theme_bw(base_size = 12) +
  theme(
    axis.title       = element_blank(),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 1, size = 10, face = "bold"),
    axis.text.y      = element_text(size = 10, color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.position  = "top",
    legend.box       = "horizontal",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.spacing.x = unit(0.5, "cm"),
    legend.box.margin = margin(0, 0, 10, 0),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5)
  ) +
  xlab("Celltype") +
  ylab("GO PATHWAY")
ggsave("/data/work/Figure/F2B_L3_celltype_GO.pdf",
       plot = p2, width = 7, height = 8, dpi = 300)
