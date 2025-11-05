## ------------------------------
## tcr L3_celltype  STARTRAC 
## ------------------------------
sce_combined_tcr <- readRDS("ASLM/TCR/GEX_TCR.rds")
sce_combined_tcr@meta.data <- sce_combined_tcr@meta.data %>%
  mutate(L3_celltype = str_replace(L3_celltype, ".*CTL.*", "CTL"),
         L3_celltype = str_replace(L3_celltype, "γδT", "gdT"))
Idents(sce_combined_tcr) <- sce_combined_tcr$L3_celltype
data_Startrac <- StartracDiversity(sce_combined_tcr,
                       type   = "type",
                       cloneCall = "aa",
                       chain = "both",
                       index = c("expa"),
                       exportTable = TRUE,
                       group.by = "orig.ident")
fwrite(data_Startrac,"ASLM/TCR/tcr_STARTRAC.csv")

data_Startrac <- data_Startrac %>%
  mutate(
    majorCluster = ifelse(str_detect(majorCluster, "Regulatory"), "CD4 Treg", majorCluster),
    expa = as.numeric(expa)  
  )

cell_type_order <- c("Naïve_T", "CD4 Treg", "Tcm", "Tem", "CTL", "Cycling T", "MAIT", "NKT", "gdT")
data_Startrac <- data_Startrac %>%
  mutate(majorCluster = factor(majorCluster, levels = cell_type_order))
kw_res <- data_Startrac %>%
  kruskal_test(expa ~ majorCluster)
dunn_res <- data_Startrac %>%
  dunn_test(expa ~ majorCluster, p.adjust.method = "BH")


naive_vs_other_expansion <- dunn_res %>%
  filter((group1 == "Naïve_T" | group2 == "Naïve_T") & p.adj.signif != "ns")
comparisons_list <- map2(naive_vs_other_expansion$group1, 
                          naive_vs_other_expansion$group2, 
                          ~c(.x, .y))
y_max <- max(data_Startrac$expa, na.rm = TRUE)
y_positions <- y_max * (1.05 + 0.03 * (0:(nrow(naive_vs_other_expansion)-1)))
pdf("ALSM/TCR/pdf/STARTRAC_celltype_TCR.pdf", width = 7.5, height = 5.5)
ggplot(data_Startrac, aes(x = majorCluster, y = expa)) +
  geom_boxplot(fill = "#E64B35", color = "black", width = 0.6, outlier.shape = NA) +
  geom_jitter(aes(fill = "#E64B35"), width = 0.15, size = 3.5, shape = 21, 
              stroke = 0.6, color = "black", alpha = 0.9) +
  labs(x = "", y = "STARTRAC expansion score", title = "") +
  theme_pubr() +               
  theme(
    plot.title = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),    
    axis.text.y = element_text(size = 10),                             
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  )+ 
    geom_signif(
      comparisons = comparisons_list,
      annotations = naive_vs_other_expansion$p.adj.signif,
      y_position = y_positions,
      tip_length = 0.01,
      textsize = 4,
      vjust = 0.5
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))
dev.off()

## ------------------------------
## tcr CTL STARTRAC expansion score(age_group)  
## ------------------------------
metadata =read.csv('ALSM/10vs10_metadata.csv')
data_Startrac_CTL <- data_Startrac[data_Startrac$majorCluster=="CTL",] %>% merge(.,metadata,by.x='group',by.y='sample')
  max_y <- max(data_Startrac_CTL$expa, na.rm = TRUE) * 1.05
  line_y <- max_y - (max(data_Startrac_CTL$expa) - min(data_Startrac_CTL$expa)) * 0.02
pdf("ALSM/TCR/pdf/STARTRAC_expansion_CTL_age_group.pdf",width  = 5.5, height = 5.5)  
plot(ggplot(data_Startrac_CTL, aes(x = age_group, y = expa, fill = age_group)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", fill = "white") +
       geom_jitter(
      aes(fill = age_group),  
      width = 0.15,
      size = 3.5,
      shape = 21,       
      stroke = 0.6,    
      color = "black",  
      alpha = 0.9       
    ) +
    stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test", label.x = 1.5,label.y = max_y, size = 5) +
    scale_color_manual(values = c("old" = "#D8383A", "young" = "#2F7FC1")) +
    labs(
      y = "STARTRAC expansion score",
      title ="CTL"
    ) +
       annotate("segment",
             x = 1, xend = 2,
             y = line_y, yend = line_y,
             size = 0.8, color = "black")+
    theme_pubr(base_size = 16) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    ))
dev.off()

## ------------------------------
## tcr cd4 tem STARTRAC expansion score(celltype)  
## ------------------------------
metadata =read.csv('ALSM/10vs10_metadata.csv')
sce_combined_tcr_cd4tem <- readRDS("ASLM/TCR/long_read_tcr_gex.rds")
startrac_cd4tem <- StartracDiversity(sce_combined_tcr_cd4tem, 
                  type = "type", 
                  exportTable = TRUE) %>% merge(.,metadata,by.x='group',by.y='sample')
max_y <- max(startrac_cd4tem$expa, na.rm = TRUE) * 1.05
line_y <- max_y - (max(startrac_cd4tem$expa) - min(startrac_cd4tem$expa)) * 0.02
pdf("ASLM/TCR/pdf/STARTRAC_expansion_cd4tem_age_group.pdf",width  = 5.5, height = 5.5)

plot(ggplot(startrac_cd4tem, aes(x = majorCluster, y = expa, fill = majorCluster)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, color = "black", fill = "white") +
       geom_jitter(
      aes(fill = age_group),  
      width = 0.15,
      size = 3.5,
      shape = 21,       
      stroke = 0.6,    
      color = "black",  
      alpha = 0.9       
    ) +
    stat_compare_means(aes(label = ..p.signif..),method = "wilcox.test", label.x = 1.5,label.y = max_y, size = 5) +
    scale_color_manual(values = c("CD4 Tem_0" = "#D8383A", "CD4 Tem_1" = "#2F7FC1")) +
    labs(
      y = "STARTRAC expansion score",
      title ="CD4 Tem_0 and CD4 Tem_1"
    ) +
       annotate("segment",
             x = 1, xend = 2,
             y = line_y, yend = line_y,
             size = 0.8, color = "black")+
    theme_pubr(base_size = 16) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    ))
dev.off()
