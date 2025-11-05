## ------------------------------
## define plotting 
## ------------------------------
index_plot_boxplot <- function(df,
                                 yvar = "shannon",
                                 group_var = "age_group",
                                 group_colors = c("old" = "#D8383A", "young" = "#2F7FC1"),
                                 ylab_text = NULL,
                                 title_text = NULL,
                                 p_method = "wilcox.test") {
  df[[group_var]] <- factor(df[[group_var]], levels = c("old", "young"))
  
  if (is.null(ylab_text)) {
    ylab_text <- yvar
  }

  max_y <- max(df[[yvar]], na.rm = TRUE) * 1.05
  line_y <- max_y - (max(df[[yvar]]) - min(df[[yvar]])) * 0.02
  # 构建 ggplot
  p <- ggplot(df, aes_string(x = group_var, y = yvar)) +
    geom_boxplot(
      width = 0.6,
      color = "black",
      fill = "white",
      outlier.shape = NA
    ) +
       geom_jitter(
      aes(fill = .data[[group_var]]),  
      width = 0.15,
      size = 3.5,
      shape = 21,       
      stroke = 0.6,    
      color = "black",  
      alpha = 0.9       
    ) +
    scale_color_manual(values = group_colors) +
    stat_compare_means(
      aes(label = ..p.signif..),
      method = p_method,
      label.x = 1.5,
      label.y = max_y,
      size = 8)+
       annotate("segment",
             x = 1, xend = 2,
             y = line_y, yend = line_y,
             size = 0.8, color = "black")+
       labs(
      x = "Age Group", 
      y = ylab_text,
      title = title_text    
    )  +
    theme_pubr(base_size = 16) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 14, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  return(p)
}

## ------------------------------
## CTL diversity 
## ------------------------------

path <-"/ASLM/TCR/pdf/"
index_list <- c("shannon","inv.simpson","chao1","gini")
metadata = read.csv('/ALSM/10vs10_metadata.csv')
##tcr diversity plot
sce_combined_tcr.TCR <- readRDS("/ASLM/TCR/GEX_TCR.rds")
sce_combined_tcr_ctl <- subset(sce_combined_tcr, L3_celltype=="CTL")

for (index in index_list) { 
    ctl_diversity <- clonalDiversity(sce_combined_tcr_ctl, 
                cloneCall = "aa",
                chain = "both",
                metric = index,
                exportTable = TRUE,
                skip.boots = TRUE)
    ctl_diversity <- merge(ctl_diversity, metadata, by = 'sample')
    pdf(paste0(path,"CTLs_", index, "_ageGroup.pdf"), width = 3.5, height = 5.5)
    index_plot_boxplot(ctl_diversity,
                                 yvar = index,
                                 group_var = "age_group",
                                 group_colors = c("old" = "#D8383A", "young" = "#2F7FC1"),
                                 ylab_text = paste0(index, " index"),  
                                 p_method = "wilcox.test")
    dev.off()
}
## ------------------------------
##cd4tem_0 and cd4tem_1 diversity 
## ------------------------------
sce_combined_tcr_cd4tem  <- readRDS("/ASLM/TCR/long_read_tcr_gex.rds")
for (index in index_list) { 
    cd4tem_diversity <- clonalDiversity(sce_combined_tcr_cd4tem, 
                cloneCall = "aa",
                chain = "TRB",
                metric = index,
                exportTable = TRUE,
                skip.boots = TRUE)
    cd4tem_diversity <- merge(cd4tem_diversity, metadata, by = 'sample')
    pdf(paste0(path,"CTLs_", index, "_ageGroup.pdf"), width = 3.5, height = 5.5)
    index_plot_boxplot(cd4tem_diversity,
                                 yvar = index,
                                 group_var = "age_group",
                                 group_colors = c("CD4 Tem_0" = "#b6dae8", "CD4 Tem_1" = "#51a3ca"),
                                 ylab_text = paste0(index, " index"),  
                                 p_method = "wilcox.test")
    dev.off()
}