## ------------------------------
## tcr umap
## ------------------------------
sce_combined_tcr <- readRDS("ASLM/TCR/GEX_TCR.rds")
# 定义克隆大小的顺序（包括 NA）
clone_size_order <- c(
  "NA",  
  "1 (0 < X <= 1)",
  "2 (1 < X <= 2)",
  "3 (2 < X <= 3)",
  "4 (3 < X <= 4)",
  "5 (4 < X <= 5)",
  "6-10 (5 < X <= 10)",
  ">10 (10 < X <= 100)"
)
# clonesize factor
colorblind_vector_tcr <- c(
  `>10 (10 < X <= 100)` = "#660B24",
    `6-10 (5 < X <= 10)` = "#9F182C",
    `5 (4 < X <= 5)` = "#C73242",
    `4 (3 < X <= 4)` = "#E26158",
    `3 (2 < X <= 3)` = "#F88E74",
    `2 (1 < X <= 2)` = "#FFB8A0",
    `1 (0 < X <= 1)` = "#FED9C9",
    `NA` = "#d7dae0ff")
##提取坐标信息
umap_tcr <- data.frame(
  barcode = colnames(sce_combined_tcr),
  cloneSize = sce_combined_tcr$cloneSize,
  Embeddings(sce_combined_tcr, reduction = "umap"),
  L3_celltype = sce_combined_tcr$L3_celltype,
  age_group  = sce_combined_tcr$group
) %>%
  mutate(
    L3_celltype_new = str_replace(L3_celltype, ".*CTL.*", "CTL"),
    cloneSize_plot = factor(
      case_when(
        is.na(cloneSize) ~ "NA",
        TRUE ~ as.character(cloneSize)
      ),
      levels = clone_size_order  
    )
  )
## 圈细胞
maskTable_TCR <- generateMask(dims=umap_tcr[,c("UMAP_1","UMAP_2")],clusters=umap_tcr$L3_celltype_new,
                        minDensity=0.1,gridSize=200,smoothSigma=0.015,type = 'partition')
# 标签定位
ctl_label_pos <- maskTable_TCR %>%
  filter(cluster == "CTL") %>%    
  group_by(group) %>%             
  summarise(
    x = max(UMAP_1)*1.05,            
    y = max(UMAP_2)*1.05          
  ) %>%
  slice(1)                       

pdf("/ALSM/TCR/pdf/clonesize_TCR.pdf", width = 7.5, height = 5.5)
plot(ggplot() +
  geom_point(
    data = subset(umap_tcr, cloneSize_plot == "NA"),
    aes(x = UMAP_1, y = UMAP_2, color = "NA"),
    size = 2,
    alpha = 1  ) +
  geom_point(
    data = subset(umap_tcr, cloneSize_plot != "NA"),
    aes(x = UMAP_1, y = UMAP_2, color = cloneSize_plot),
    size = 2,
    alpha = 1
  ) +
scale_color_manual(
    name = "Clone Size",
    values = colorblind_vector_tcr,
    breaks = clone_size_order,  
    drop = FALSE  
  ) +
  geom_path(
    data = maskTable_TCR %>% filter(cluster == "CTL"),
    aes(x = UMAP_1, y = UMAP_2, group = group),
    color = "#d42c24",
    linewidth = 0.8,
    linetype =5,
    alpha = 100
  ) +
  annotate(
    "text",
    x = ctl_label_pos$x, 
    y = ctl_label_pos$y,  
    label = "CTL",
    color = "#d42c24",
    size = 5,
    fontface = "bold"
  )+
  labs(x = "UMAP 1", y = "UMAP 2") +
 theme_classic() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.line.x.top = element_blank(),    
    axis.line.y.right = element_blank(), 
    axis.ticks = element_blank(),         
    axis.text = element_blank(),  
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white")
  )
)
dev.off()

## ------------------------------
##bcr umap
## ------------------------------
sce_combined_bcr <- readRDS("ASLM/BCR/GEX_BCR.rds")
colorblind_vector_bcr <- c(
  `>10 (10 < X <= 100)` =  "#061c36ff",
  `6-10 (5 < X <= 10)` =  "#163b65ff",
  `5 (4 < X <= 5)` = "#2e5786ff",
  `4 (3 < X <= 4)` = "#4872a1ff",
  `3 (2 < X <= 3)` = "#7397c0ff",
  `2 (1 < X <= 2)` = "#9db8d7ff",
  `1 (0 < X <= 1)` = "#ccdef3ff",
  "NA" = "#B7BBC4"
)
###umap information 
umap_bcr <- data.frame(
  barcode = colnames(sce_combined_bcr),
  cloneSize = sce_combined_bcr$cloneSize,
  Embeddings(sce_combined_bcr, reduction = "umap"),
  L3_celltype = sce_combined_bcr$L3_celltype,
  age_group  = sce_combined_bcr$group
) %>%
  mutate(
    cloneSize_plot = factor(
      case_when(
        is.na(cloneSize) ~ "NA",
        TRUE ~ as.character(cloneSize)
      ),
      levels = clone_size_order  
    )
  )
maskTable_BCR <- generateMask(dims=umap_bcr[,c("UMAP_1","UMAP_2")],clusters=umap_bcr$L3_celltype,
                        minDensity=0.1,gridSize=200,smoothSigma=0.015,type = 'partition')

plasma_label_pos <- maskTable_BCR %>%
  filter(cluster == "Plasma") %>%    
  group_by(group) %>%             
  summarise(
    x = max(UMAP_1)*1.05,            
    y = max(UMAP_2)*1.05            
  ) %>%
  slice(1) 

pdf("ALSM/BCR/pdf/clonesize_BCR.pdf", width = 7.5, height = 5.5)
plot(ggplot() +
  geom_point(
    data = subset(umap_bcr, cloneSize_plot == "NA"),
    aes(x = UMAP_1, y = UMAP_2, color = "NA"),
    size = 2,
    alpha = 1  ) +
  geom_point(
    data = subset(umap_bcr, cloneSize_plot != "NA"),
    aes(x = UMAP_1, y = UMAP_2, color = cloneSize_plot),
    size = 2,
    alpha = 1
  ) +
scale_color_manual(
    name = "Clone Size",
    values = colorblind_vector_bcr,
    breaks = clone_size_order,
    drop = FALSE  
  ) +
  geom_path(
    data = maskTable_BCR %>% filter(cluster == "Plasma"),
    aes(x = UMAP_1, y = UMAP_2, group = group),
    color = "#061c36ff",
    linewidth = 0.8,
    linetype =5,
    alpha = 100
  ) +
  annotate(
    "text",
    x = plasma_label_pos$x, 
    y = plasma_label_pos$y,  
    label = "Plasma",
    color = "#061c36ff",
    size = 5,
    fontface = "bold"
  )+
  labs(x = "UMAP 1", y = "UMAP 2") +
 theme_classic() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.line.x.top = element_blank(),    
    axis.line.y.right = element_blank(), 
    axis.ticks = element_blank(),         
    axis.text = element_blank(),  
    panel.background = element_rect(fill = 'white'), # 白色背景
    plot.background = element_rect(fill = "white"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white")
  )
)
dev.off()


## ------------------------------
## cd4 tem0_1 umap
## ------------------------------
sce_combined_tcr_cd4tem <- readRDS("ASLM/TCR/long_read_tcr_gex.rds")
cd4tem_umap <- fread('long_read/CD4Tem_umap_coords.csv')
umap_cd4tem <- data.frame(
  barcode = colnames(sce_combined_tcr_cd4tem),
  cloneSize = sce_combined_tcr_cd4tem$cloneSize,
  celltype = sce_combined_tcr_cd4tem$celltype
) %>% merge(.,cd4tem_umap,by='barcode')%>%
  mutate(
    cloneSize_plot = factor(
      case_when(
        is.na(cloneSize) ~ "NA",
        TRUE ~ as.character(cloneSize)
      ),
      levels = clone_size_order  
    )
  )
colorblind_vector_cd4tem <- c(
  `>10 (10 < X <= 100)` =  "#061c36ff",
  `6-10 (5 < X <= 10)` =  "#163b65ff",
  `5 (4 < X <= 5)` = "#2e5786ff",
  `4 (3 < X <= 4)` = "#4872a1ff",
  `3 (2 < X <= 3)` = "#7397c0ff",
  `2 (1 < X <= 2)` = "#9db8d7ff",
  `1 (0 < X <= 1)` = "#ccdef3ff",
  "NA" = "#B7BBC4"
)
pdf("ALSM/TCR/pdf/clonesize_cd4tem_umap.pdf", width = 7.5, height = 5.5)
plot(ggplot() +
  geom_point(
    data = subset(umap_cd4tem, cloneSize_plot == "NA"),
    aes(x = UMAP_1, y = UMAP_2, color = "NA"),
    size = 2,
    alpha = 1  ) +
  geom_point(
    data = subset(umap_cd4tem, cloneSize_plot != "NA"),
    aes(x = UMAP_1, y = UMAP_2, color = cloneSize_plot),
    size = 2,
    alpha = 1
  ) +
scale_color_manual(
    name = "Clone Size",
    values = colorblind_vector_cd4tem,
    breaks = clone_size_order,  
    drop = FALSE  
  ) +
  labs(x = "UMAP 1", y = "UMAP 2") +
 theme_classic() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.line.x.top = element_blank(),    
    axis.line.y.right = element_blank(), 
    axis.ticks = element_blank(),         
    axis.text = element_blank(),  
    panel.background = element_rect(fill = 'white'), 
    plot.background = element_rect(fill = "white"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white")
  )
)
dev.off()

## ------------------------------
## tcr age_group barplot 
## ------------------------------
df_age_tcr <- umap_tcr
freq_table_age_tcr <- table(df_age_tcr$age_group, df_age_tcr$cloneSize, useNA = "no")
freq_percent_age_tcr <- as.data.frame.matrix(prop.table(freq_table_age_tcr, margin = 1) )
freq_percent_age_tcr$age_group <- rownames(freq_percent_age_tcr)
freq_percent_age_tcr <- freq_percent_age_tcr[ , !names(freq_percent_age_tcr) %in% c("NA")]
freq_long_age_tcr <- freq_percent_age_tcr %>%
  melt(id.vars = "age_group", 
       variable.name = "cloneSize", 
       value.name = "abundance") %>%
  filter(!is.na(cloneSize) & cloneSize != "None ( < X <= 0)")
my_colors_age_tcr <- c("#FED9C9","#FFB8A0","#F88E74","#E26158","#C73242","#9F182C","#660B24")
freq_long_age_tcr <- freq_long_age_tcr %>%
  mutate(cloneSize = factor(cloneSize, levels = clone_size_order))
pdf("ALSM/TCR/pdf/clonesize_barplot_TCR.pdf", width = 3.5, height = 5.5)
ggplot(freq_long_age_tcr, aes(x = age_group, y = abundance, fill = cloneSize)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.25) + 
  scale_fill_manual(name = "Clone Size", values = my_colors_age_tcr) + 
  labs(
    title = "",
    x = "",
    y = "Relative Percentage of Cells" 
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_break(c(0.2,0.8), scales=0.2, ticklabels=c(0.8,0.1), space=0.2)
dev.off()

## ------------------------------
##bcr age_group barplot
## ------------------------------
df_age_bcr <- umap_bcr
freq_table_age_bcr <- table(df_age_bcr$age_group, df_age_bcr$cloneSize, useNA = "no")
freq_percent_age_bcr <- as.data.frame.matrix(prop.table(freq_table_age_bcr, margin = 1) )
freq_percent_age_bcr$age_group <- rownames(freq_percent_age_bcr)
freq_percent_age_bcr <- freq_percent_age_bcr[ , !names(freq_percent_age_bcr) %in% c("NA")]
freq_long_age_bcr <- freq_percent_age_bcr %>%
  melt(id.vars = "age_group", 
       variable.name = "cloneSize", 
       value.name = "abundance") %>%
  filter(!is.na(cloneSize) & cloneSize != "None ( < X <= 0)")
my_colors_age_bcr <- c("#ccdef3ff","#9db8d7ff","#7397c0ff","#4872a1ff","#2e5786ff","#163b65ff","#061c36ff")
freq_long_age_bcr <- freq_long_age_bcr %>%
  mutate(cloneSize = factor(cloneSize, levels = clone_size_order))
pdf("ALSM/BCR/pdf/clonesize_barplot_BCR.pdf", width = 3.5, height = 5.5)
ggplot(freq_long_age_bcr, aes(x = age_group, y = abundance, fill = cloneSize)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.25) + 
  scale_fill_manual(name = "Clone Size", values = my_colors_age_bcr) + 
  labs(
    title = "",
    x = "",
    y = "Relative Percentage of Cells" 
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_break(c(0.2,0.8), scales=0.2, ticklabels=c(0.8,0.1), space=0.2)
dev.off()

## ------------------------------
##tcr age_group+l3_celltype  barplot
## ------------------------------
df_age_L3_tcr <- umap_tcr %>%
  mutate(
    L3_celltype = case_when(
      L3_celltype == "γδT" ~ "gdT",
      TRUE ~ as.character(L3_celltype)
    ),
    group_label = paste(L3_celltype, age_group, sep = ".")
  )

freq_table_L3_age_tcr <- table(df_age_L3_tcr$group_label, df_age_L3_tcr$cloneSize, useNA = "no")
freq_percent_L3_age_tcr <- as.data.frame.matrix(prop.table(freq_table_L3_age_tcr, margin = 1) )
freq_percent_L3_age_tcr$group_label <- rownames(freq_percent_L3_age_tcr)
freq_percent_L3_age_tcr <- freq_percent_L3_age_tcr[ , !names(freq_percent_L3_age_tcr) %in% c("NA")]
# 转换为长格式
freq_long_L3_age_tcr <- freq_percent_L3_age_tcr %>%
  melt(id.vars = "group_label", 
       variable.name = "cloneSize", 
       value.name = "abundance") %>%
  filter(!is.na(cloneSize) & cloneSize != "None ( < X <= 0)")
celltype_order <- c("CD4 CTL.YOUNG","CD4 CTL.OLD","CD8 CTL.YOUNG","CD8 CTL.OLD", "NKT.YOUNG", "NKT.OLD","CD8 Tem.YOUNG","CD8 Tem.OLD", 
                    "Cycling T.YOUNG","Cycling T.OLD",
                    "CD8 Tcm.YOUNG","CD8 Tcm.OLD", "MAIT.YOUNG", "MAIT.OLD","gdT.YOUNG","gdT.OLD", "CD4 Tcm.YOUNG","CD4 Tcm.OLD", "CD4 Regulatory.YOUNG","CD4 Regulatory.OLD",
                    "CD4 Tem.YOUNG", "CD4 Tem.OLD",
                    "CD8 Naïve.YOUNG","CD8 Naïve.OLD", "CD4 Naïve.YOUNG","CD4 Naïve.OLD")

my_colors_L3_age_tcr <- c("white","#FFE5E5","#F9B7B7","#F38181","#E74C3C","#D32F2F","#8B0000")
freq_long_L3_age_tcr <- freq_long_L3_age_tcr %>%
  mutate(cloneSize = factor(cloneSize, levels = clone_size_order),
         group_label=factor(group_label,levels = celltype_order))
pdf("ALSM/TCR/pdf/clonesize_barplot_L3+age_TCR.pdf", width = 7.5, height = 5.5)
ggplot(freq_long_L3_age_tcr, aes(x = group_label, y = abundance, fill = cloneSize)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.25) + 
  scale_fill_manual(name = "Clone Size", values = my_colors_L3_age_tcr) + 
  labs(
    title = "",
    x = "",
    y = "Relative Percentage of Cells" 
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_break(c(0.75,0.95), scales=0.05, ticklabels=c(0.75,0.95,1), space=0.01)



## ------------------------------
## cd4 tem0_1 barplot
## ------------------------------
sce_combined_tcr_cd4tem <- readRDS("ASLM/TCR/long_read_tcr_gex.rds")
df_cd4tem <- umap_cd4tem
freq_table_cd4tem <- table(df_cd4tem$age_group, df_cd4tem$cloneSize, useNA = "no")
freq_percent_cd4tem <- as.data.frame.matrix(prop.table(freq_table_cd4tem, margin = 1) )
freq_percent_cd4tem$age_group <- rownames(freq_percent_cd4tem)
freq_percent_cd4tem <- freq_percent_cd4tem[ , !names(freq_percent_cd4tem) %in% c("NA")]
freq_long_cd4tem <- freq_percent_cd4tem %>%
  melt(id.vars = "celltype", 
       variable.name = "cloneSize", 
       value.name = "abundance") %>%
  filter(!is.na(cloneSize) & cloneSize != "None ( < X <= 0)")
my_colors_cd4tem <- c("#ccdef3ff","#9db8d7ff","#7397c0ff","#4872a1ff","#2e5786ff","#163b65ff","#061c36ff")
freq_long_cd4tem <- freq_long_cd4tem %>%
  mutate(cloneSize = factor(cloneSize, levels = clone_size_order))
pdf("ALSM/TCR/pdf/clonesize_barplot_CD4Tem.pdf", width = 3.5, height = 5.5)
ggplot(freq_long_cd4tem, aes(x = celltype, y = abundance, fill = cloneSize)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.25) + 
  scale_fill_manual(name = "Clone Size", values = my_colors_cd4tem) + 
  labs(
    title = "",
    x = "",
    y = "Relative Percentage of Cells" 
  ) +
  theme_minimal() +
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_break(c(0.1,0.8), scales=0.2, ticklabels=c(0.8,0.1), space=0.2)
dev.off()
