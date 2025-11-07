#gene_set_score-CD4_Tem
sce <- readRDS("/data/work/Tem_01_02.rds")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)

input_dir <- "/data/Files/GENE_SET/PAN_cancer_T_cell_atlas_links_a_cellular_CD4T"
output_dir <- "/data/Figure/CD4_Tem_pathway/"
gene_set_files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)

#dataframe to store results
results <- data.frame(
  pathway = character(),
  n_genes = integer(),
  n_common = integer(),
  missing_rate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)
pathway_gene_list <- list()
all_genes <- rownames(sce)
for (file in gene_set_files) {
  pathway_name <- tools::file_path_sans_ext(basename(file))
  gene_list <- readLines(file)
  pathway_gene_list[[pathway_name]] <- gene_list
  common_genes <- intersect(gene_list, all_genes)
  n_total <- length(gene_list)
  n_common <- length(common_genes)
  if (n_common < 5) {
    message(sprintf("step skipped for pathway %s: not enough common genes (%d)", pathway_name, n_common))
    next
  }
  #analyze module scores
  sce_temp <- AddModuleScore(
    object = sce,
    features = list(common_genes),
    ctrl = 50,
    name = "pathway_score_",
    slot = 'data'
  )
  scores <- sce_temp@meta.data$pathway_score_1
  group1_scores <- scores[sce_temp$cell_type == "CD4 Tem_0"]
  group2_scores <- scores[sce_temp$cell_type == "CD4 Tem_1"]
  test_result <- wilcox.test(group1_scores, group2_scores)
  results <- rbind(results, data.frame(
    pathway = pathway_name,
    n_genes = n_total,
    n_common = n_common,
    missing_rate = 1 - n_common/n_total,
    p_value = test_result$p.value
  ))
}
# Adjust p-values for multiple testing
results <- results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)
sig_results <- filter(results, p_adj < 0.05)
sig_results <- sig_results %>%
  filter(pathway %in% c("Naïve",
                        "TCR signaling",
                        "Stress response",
                        "Cytokine",
                        "Pro-apoptotic",
                        "cytotoxicity"))
# Plotting significant pathways
if (nrow(sig_results) > 0) {
  for (i in 1:nrow(sig_results)) { # nolint
    pathway_name <- sig_results$pathway[i]
    gene_list <- pathway_gene_list[[pathway_name]]
    common_genes <- intersect(gene_list, all_genes)
    sce_temp <- AddModuleScore(
      object = sce,
      features = list(common_genes),
      ctrl = 50,
      name = "pathway_score_",
      slot = 'data'
    )
    # Prepare data for plotting
    plot_data <- sce_temp@meta.data %>%
      select(cell_type, pathway_score_1) %>%
      rename(Score = pathway_score_1)
    # Generate violin plot
    p <- ggplot(plot_data, aes(x = cell_type, y = Score, fill = cell_type)) +
      geom_violin(
        position = position_dodge(0.9),
        alpha = 0.7,
        width = 0.7,
        trim = TRUE,
        scale = "width"
      ) +
      geom_boxplot(
        width = 0.15,
        position = position_dodge(0.9),
        alpha = 0.3,
        outlier.shape = NA
      ) +
      stat_compare_means(
        aes(group = cell_type),
        method = "wilcox.test",
        label = "p.signif",
        hide.ns = FALSE,
        show.legend = FALSE,
         comparisons = list(c("CD4 Tem_1", "CD4 Tem_0")),
        step.increase = 0.1
      ) +
      labs(
           title = pathway_name,
           x = "celltype", y = "Pathway Score", fill = "celltype") +
        scale_fill_manual(values = c("CD4 Tem_1" = "#054168ff", "CD4 Tem_0" = "#84e2f0ff")) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA, size = 0.5),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      scale_y_continuous(breaks = seq(-1, 1.5, by = 0.05))
    file_name <- file.path(output_dir, paste0(gsub("[^[:alnum:]]", "_", pathway_name), ".pdf"))
    ggsave(file_name, p, width = 6, height = 12)
  }
} else {
  message("no significant pathways found.")
}