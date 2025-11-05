## ------------------------------
## CTL cdr3 length  
## ------------------------------
combined.TCR <- readRDS("ASLM/TCR/TCR.rds")
combined.TCR_clonalLength_TRB<- clonalLength(combined.TCR, 
             cloneCall="aa", 
             group.by = c("age_group","celltype_L3","barcode"),
              exportTable = TRUE,
             chain = "TRB")
combined.TCR_clonalLength_TRB<- combined.TCR_clonalLength_TRB %>%
  mutate(celltype_L3 = str_replace(celltype_L3, ".*CTL.*", "CTL"))
df_subset_ctl <- subset(combined.TCR_clonalLength_TRB,celltype_L3=="CTL")

age_colors <- c("old" = "#D8383A", "young" = "#2F7FC1")
pdf("./pdf/TRB_cdr3length_aa_ageGroup_CTL.pdf", width  = 3, height = 5.5)
plot(ggplot(df_subset_ctl, aes(x = .data[["age_group"]], 
                                             y = .data[["length"]],
                                             fill = .data[["age_group"]])) +
  geom_violin(
    width = 0.8,
    color = "black",
    alpha = 0.8,
    trim = FALSE,      
    scale = "width",
    linewidth = 0.5,
      adjust=1.5
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    alpha = 0.7,
    outlier.shape = NA 
  ) +
  geom_signif(
    comparisons = list(c("old", "young")), 
    map_signif_level = TRUE,               
    textsize = 4,                         
    tip_length = 0,                       
    y_position = 25,                      
    color = "black",                      
    vjust = 0.5                          
  ) +
  labs(
    x = "", 
    y = "CDR3 Length(aa)",
    fill = "
"
  ) +
  scale_fill_manual(values = age_colors) +
  theme_pubr(base_size = 16) +
  coord_cartesian(ylim = c(5,25)) + 
  scale_y_continuous(breaks = seq(10, 20, by = 10))+
  theme(
    legend.title = element_blank(),
    legend.position = c(0.9,0.8),
    legend.background = element_blank(),
    text = element_text(size = 12, family = "Arial", color = "black"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.text.x = element_text(size = 12, color = "black", vjust = 0.2),
    axis.title.x = element_text(size = 15, color = "black", vjust = -2),
    axis.text.y = element_text(size = 12, color = "black", hjust = 1),
    axis.title.y = element_text(size = 15, color = "black", vjust = 2.5),
    axis.ticks = element_line(size = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15)), 
    labels = scales::comma
  ))
dev.off()
## ------------------------------
## Switch.B cdr3 length  
## ------------------------------
combined.BCR <- readRDS("ASLM/TCR/BCR.rds")
combined.BCR_clonalLength_IGH<- clonalLength(combined.BCR, 
             cloneCall="aa", 
             group.by = c("age_group","celltype_L3","barcode"),
              exportTable = TRUE,
             chain = "IGH")

df_subset_Swit_Bm <- subset(combined.BCR_clonalLength_IGH,celltype_L3=="Swit. Bm")
##CTL cdr3 length
age_colors <- c("old" = "#D8383A", "young" = "#2F7FC1")
pdf("ASLM/BCR/pdf/IGH_cdr3length_aa_ageGroup_Swit_Bm.pdf", width  = 3, height = 5.5)
plot(ggplot(df_subset_Swit_Bm, aes(x = .data[["age_group"]], 
                                             y = .data[["length"]],
                                             fill = .data[["age_group"]])) +
  geom_violin(
    width = 0.8,
    color = "black",
    alpha = 0.8,
    trim = FALSE,      
    scale = "width",
    linewidth = 0.5,
      adjust=1.5
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    alpha = 0.7,
    outlier.shape = NA 
  ) +
  geom_signif(
    comparisons = list(c("old", "young")), 
    map_signif_level = TRUE,               
    textsize = 4,                         
    tip_length = 0,                       
    y_position = 25,                      
    color = "black",                      
    vjust = 0.5                          
  ) +
  labs(
    x = "", 
    y = "CDR3 Length(aa)",
    fill = "
"
  ) +
  scale_fill_manual(values = age_colors) +
  theme_pubr(base_size = 16) +
  coord_cartesian(ylim = c(5,25)) + 
  scale_y_continuous(breaks = seq(10, 20, by = 10))+
  theme(
    legend.title = element_blank(),
    legend.position = c(0.9,0.8),
    legend.background = element_blank(),
    text = element_text(size = 12, family = "Arial", color = "black"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.text.x = element_text(size = 12, color = "black", vjust = 0.2),
    axis.title.x = element_text(size = 15, color = "black", vjust = -2),
    axis.text.y = element_text(size = 12, color = "black", hjust = 1),
    axis.title.y = element_text(size = 15, color = "black", vjust = 2.5),
    axis.ticks = element_line(size = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15)), 
    labels = scales::comma
  ))
dev.off()
## ------------------------------
## CD4 tem 0_1 CDR3 length  
## ------------------------------
age_colors <- c("CD4 Tem_0" = "#b8c7de", "CD4 Tem_1" = "#a2a6c8")
combined.TCR_CD4_Tem_01<- readRDS("ASLM/TCR/long_read_tcr.rds")
combined.TCR_CD4_Tem_01_clonalLength_TRB<- clonalLength(combined.TCR_CD4_Tem_01, 
             cloneCall='aa', 
             group.by = c("barcode","cell_type"),
              exportTable = TRUE,
             chain = "TRB") 
pdf("ASLM/TCR/pdf/TRB_cdr3length_aa_ageGroup_CD4_Tem_01.pdf", width  = 3, height = 5.5)
plot(ggplot(combined.TCR_CD4_Tem_01_clonalLength_TRB, aes(x = .data[["cell_type"]], 
                                             y = .data[["length"]],
                                             fill = .data[["cell_type"]])) +
  geom_violin(
    width = 0.8,
    color = "black",
    alpha = 0.8,
    trim = FALSE,      
    scale = "width",
    linewidth = 0.5,
      adjust=1.5
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    alpha = 0.7,
    outlier.shape = NA 
  ) +
  geom_signif(
    comparisons = list(c("old", "young")), 
    map_signif_level = TRUE,               
    textsize = 4,                         
    tip_length = 0,                       
    y_position = 25,                      
    color = "black",                      
    vjust = 0.5                          
  ) +
  labs(
    x = "", 
    y = "CDR3 Length(aa)",
    fill = "
"
  ) +
  scale_fill_manual(values = age_colors) +
  theme_pubr(base_size = 16) +
  coord_cartesian(ylim = c(5,25)) + 
  scale_y_continuous(breaks = seq(10, 20, by = 10))+
  theme(
    legend.title = element_blank(),
    legend.position = c(0.9,0.8),
    legend.background = element_blank(),
    text = element_text(size = 12, family = "Arial", color = "black"),
    axis.line = element_line(size = 0.5, colour = "black"),
    axis.text.x = element_text(size = 12, color = "black", vjust = 0.2),
    axis.title.x = element_text(size = 15, color = "black", vjust = -2),
    axis.text.y = element_text(size = 12, color = "black", hjust = 1),
    axis.title.y = element_text(size = 15, color = "black", vjust = 2.5),
    axis.ticks = element_line(size = 0.5),
    axis.ticks.length = unit(0.3, "cm"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15)), 
    labels = scales::comma
  ))
dev.off()

