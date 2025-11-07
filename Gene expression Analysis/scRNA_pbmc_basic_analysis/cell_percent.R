
library(Seurat) 
library(tidyverse)
library(ggpubr) 
library(ggalluvial)
library(viridis)  
library(RColorBrewer) 
library(reshape2) 

L3_color_map <- c(`HSPC` = "#FAD7D7", `cMono` = "#CF4457", `intMono` = "#DE6A69", `ncMono` = "#EA5958", `cDC2` = "#F29667", `cDC1` = "#EF8787", `pDC` = "#FBC5AB", `MK` = "#F6BABB",
                  `Swit. Bm` = "#D48CB1", `Unswit. Bm` = "#8085C0", `Naïve B` = "#D4A7C4", `Plasma` = "#DBC6D9", `Plasmablast` = "#E6E6F0",
                  `NKbright` = "#F07D50", `NKdim` = "#FBAA34", `XCL+NKdim` = "#FEC928", `Terminal NK` = "#F4EB1B",
                  `CD4 Naïve` = "#015696", `CD4 Regulatory` = "#89CAEA", `CD4 CTL` = "#BBE6FA", `CD4 Tcm` = "#0B75B3", `CD4 Tem` = "#4596CD",
                  `NKT` = "#87D649", `MAIT` = "#47A183", `γδT` = "#24858E", `Cycling T` = "#77C3BC",
                  `CD8 CTL` = "#C3E023", `CD8 Naïve` = "#8DCA9D", `CD8 Tem` = "#2BB17E", `CD8 Tcm` = "#50C56A"
)
#All PBMCs cell type color map

cell <- readRDS("/data/RDS/Annotation/10_PBMC_annotation_L3.rds")
sce_cell <- data.frame(
  persentage = numeric(),
  celltype = character(),
  group = character()
)
for (g in levels(factor(pbmc$group))) {
  celltype_counts <- table(pbmc$L3_celltype[pbmc$group == g])
  relative_abundance <- prop.table(celltype_counts)
  for (ct in names(relative_abundance)) {
    sce_cell <- rbind(sce_cell, data.frame(
      persentage = as.numeric(relative_abundance[ct]),
      celltype = ct,
      group = g
    ))
  }
}
is_alluvia_form(sce_cell, silent = TRUE)
sce_cell$group <- factor(sce_cell$group, levels = c("YOUNG", "OLD"))
sce_cell$celltype <- factor(sce_cell$celltype,
                            levels = c("HSPC", "cMono", "ncMono", "intMono","cDC1", "cDC2", "pDC","MK",
                                       "Swit. Bm","Unswit. Bm", "Naïve B", "Plasma","Plasmablast",
                                       "NKbright", "NKdim", "XCL+NKdim","Terminal NK",
                                       "CD4 Naïve","CD4 Tcm","CD4 Tem" , "CD4 Treg", "CD4 CTL",
                                       "NKT", "MAIT", "gdT", "Cycling T", "CD8 CTL", "CD8 Naïve", "CD8 Tcm", "CD8 Tem"))
options(repr.plot.height = 8,repr.plot.width = 8)
#Figure1C: L3 cell type relative abundance alluvial plot
p <- ggplot(as.data.frame(sce_cell),
            aes(x=group,
                y=persentage*100,
                fill=celltype,
                stratum = celltype,
                alluvium = celltype)) +
  geom_bar(stat = "identity", width = 0.45) +
  geom_alluvium() +
  geom_stratum(width = 0.45, size = 0.1) +
  scale_fill_manual(values = L3_color_map) +
  labs(x = "Age", y = "Relative Abundance (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p
ggsave("/data/Figure/F1C_L3_Relative_Abundance_celltype.pdf", plot = p, width = 8, height = 8)


# CD8 naive of PBMCs cell : percent boxplot
cell <-  readRDS("/data/RDS/Annotation/10_PBMC_annotation_L3.rds")
Idents(cell) <- "L3_celltype"
Cellratio <- prop.table(table(Idents(cell), cell$orig.ident), margin = 2)
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio, Var2 ~ Var1, value.var = "Freq")
rownames(cellper) <- cellper[, 1]
cellper <- cellper[, -1]
sample <- c("Y01" ,"Y02","Y03","Y04","Y05","Y06","Y07","Y08","Y09","O01","O02","O03","O04","O05","O06","O07","O08","O09","O10") 
group <- c("YOUNG", "YOUNG", "YOUNG","YOUNG", "YOUNG", "YOUNG","YOUNG", "YOUNG", "YOUNG", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD")
samples <- data.frame(sample, group)
rownames(samples) <- samples$sample
cellper$sample <- samples[rownames(cellper), 'sample']
cellper$group <- samples[rownames(cellper), 'group']
#CD8 Naïve cell percent
cellper_ <- cellper %>% select(sample, group, "CD8 Naïve")
colnames(cellper_) <- c('sample', 'group', 'percent')
cellper_$percent <- as.numeric(cellper_$percent)
cellper_ <- cellper_ %>% group_by(group) %>% mutate(
      upper = quantile(percent, 0.75),
      lower = quantile(percent, 0.25),
      mean = mean(percent),
      median = median(percent))
max_y <- max(cellper_$percent, na.rm = TRUE) * 1.05
line_y <- max_y - (max(cellper_$percent) - min(cellper_$percent)) * 0.02
#Figure
pp1 <- ggplot(cellper_, aes(x = group, y = percent, fill = group)) +
  geom_boxplot(position = position_dodge(0.9),
               alpha = 0.7,
               outlier.shape = NA,
               width = 0.7) +
  geom_jitter(shape = 21, width = 0.15, stroke = 0.6, color = "black",
              alpha = 0.9, aes(fill = group)) +
  labs(title = paste0("CD8 Naïve of PBMCs "), y = "Percentage", x = "Group") +
  scale_fill_manual(values = c("OLD" = "#E41A1C", "YOUNG" = "#377EB8")) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", label.x = 1.5,
                     label.y = max_y, size = 5) +
  annotate("segment",
           x = 1, xend = 2,
           y = line_y, yend = line_y,
           size = 0.8, color = "black") +
  theme_pubr(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
pp1
ggsave("/data/Figure/sF_CD8_naive_cellpercent.pdf", pp1, width = 6, height = 6)


#CD4 Treg cell of CD4 T cell : percent boxplot
Idents(cell) <- "L3_celltype"
Cellratio <- prop.table(table(Idents(cell), cell$orig.ident), margin = 2) 
Cellratio <- data.frame(Cellratio)
cellper <- dcast(Cellratio, Var2 ~ Var1, value.var = "Freq") 
rownames(cellper) <- cellper[, 1]
cellper <- cellper[, -1]
sample <- c("Y01" ,"Y02","Y03","Y04","Y05","Y06","Y07","Y08","Y09","O01","O02","O03","O04","O05","O06","O07","O08","O09","O10") 
group <- c("YOUNG", "YOUNG", "YOUNG","YOUNG", "YOUNG", "YOUNG","YOUNG", "YOUNG", "YOUNG", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD", "OLD")
samples <- data.frame(sample, group)
rownames(samples) <- samples$sample
cellper$sample <- samples[rownames(cellper), "sample"]
cellper$group <- samples[rownames(cellper), "group"]
#CD4 Treg
cellper_ <- cellper %>% select(sample, group, "CD4 Treg")
colnames(cellper_) <- c("sample", "group", "percent")
cellper_$percent <- as.numeric(cellper_$percent)
cellper_ <- cellper_ %>% group_by(group) %>% mutate(
      upper = quantile(percent, 0.75),
      lower = quantile(percent, 0.25),
      mean = mean(percent),
      median = median(percent))
max_y <- max(cellper_$percent, na.rm = TRUE) * 1.05
line_y <- max_y - (max(cellper_$percent) - min(cellper_$percent)) * 0.02

#Figure
pp1 <- ggplot(cellper_, aes(x = group, y = percent, fill = group)) +
  geom_boxplot(position = position_dodge(0.9),
               alpha = 0.7,
               outlier.shape = NA,
               width = 0.7) +
  geom_jitter(shape = 21, width = 0.15,stroke = 0.6, color = "black",
              alpha = 0.9, aes(fill = group)) +
  labs(title = paste0(" CD4 Regulatory of total ", type),
       y = "Percentage", x = "Group") +
  scale_fill_manual(values = c("OLD" = "#E41A1C", "YOUNG" = "#377EB8")) +
  stat_compare_means(aes(label = ..p.signif..),
                     method = "wilcox.test", label.x = 1.5,
                     label.y = max_y, size = 5) +
  annotate("segment",
           x = 1, xend = 2,
           y = line_y, yend = line_y, size = 0.8, color = "black") +
  theme_pubr(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )
pp1
ggsave("/data/Figure/sF_CD4_T_reg_cellpercent.pdf", pp1, width = 6, height = 6)


# density plot
library(tidyverse)
library(Seurat)
library(ggplot2)
library(viridis)
set.seed(1)
# plot Function
plot_umap_density <- function(data, type, base_size = 12, base_family = "") {
  ggplot(data = data, aes(x = UMAP_1, y = UMAP_2)) + # nolint
    theme_black(base_size, base_family) +
    xlim(c(min(data$UMAP_1) - 2, max(data$UMAP_1) + 2)) +
    ylim(c(min(data$UMAP_2) - 2, max(data$UMAP_2) + 2)) +
    stat_density_2d(aes(fill=..density..), geom = "raster", contour = FALSE, alpha = 1) + # nolint
    geom_point(color = "grey90", size = 0.05, alpha = 0.4) +
    scale_fill_viridis(option = "magma", alpha = 1) +
    #labs(title = paste0(celltype, " UMAP Density Plot for ", type)) +
    theme(legend.position = "none")
}

# theme Function
theme_black <- function(base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white", fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size * 0.8, color = "white"),
      legend.title = element_text(size = base_size * 0.8,
                                  face = "bold", hjust = 0, color = "white"),
      legend.position = "none",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      panel.background = element_rect(fill = "black", color = NA),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, "lines"),
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size * 0.8, color = "white"),
      strip.text.y = element_text(size = base_size * 0.8, 
                                  color = "white", angle = -90),
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size * 1.2, color = "white"),
      plot.margin = unit(rep(0, 4), "lines")
    )
}
pbmc <- readRDS("/data/RDS/Annotation/10_PBMC_annotation_L3.rds")
pbmc
table(pbmc$celltype)

options(repr.plot.height = 8, repr.plot.width = 8)
pbmc <- SetIdent(pbmc, value = "group")
pbmc <- subset(pbmc, downsample = 40000)
pbmc <- SetIdent(pbmc, value = "celltype")
pbmc$type <- pbmc$group
meta <- as.data.frame(pbmc@meta.data)
meta <- data.frame(ID = rownames(meta), meta, stringsAsFactors = FALSE)

coord <- Embeddings(object = pbmc, reduction = "umap")
coord <- coord[, c(1, 2)]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)
meta <- left_join(meta, coord, by = "ID")
unique_types <- unique(meta$type)
lapply(unique_types, function(tissue_type) {
  tissue_data <- meta %>% filter(type == tissue_type)

  p <- plot_umap_density(tissue_data, tissue_type)
  ggsave(filename = paste0("/data/Figure/", 
         tissue_type, "_umap_density.png"),
         width = 10, height = 10, plot = p, dpi = 600)
  p
})

# CD4T
CD4T <- readRDS("/data/RDS/Annotation/6_CD4T_annotation.rds")
CD4T
table(CD4T$group)
table(CD4T$L3_celltype)
options(repr.plot.height = 10, repr.plot.width = 10)
cell <- CD4T
cell$type <- cell$group 
minNum <- min(table(cell$type))
cell <- SetIdent(cell, value = "type")
cell <- subset(cell, downsample = minNum)
cell <- SetIdent(cell, value = "celltype")
meta <- as.data.frame(cell@meta.data)
meta <- data.frame(ID = rownames(meta), meta, stringsAsFactors = FALSE)
coord <- Embeddings(object = cell, reduction = "umap")
coord <- coord[, c(1, 2)]
colnames(coord) <- c("UMAP_1", "UMAP_2")
coord <- data.frame(ID = rownames(coord), coord)
meta <- left_join(meta, coord, by = "ID") 

unique_types <- unique(meta$type)
lapply(unique_types, function(tissue_type) {
  tissue_data <- meta %>% filter(type == tissue_type)
  p <- plot_umap_density(tissue_data, tissue_type)
  ggsave(filename = paste0("/data/Figure/",
         tissue_type, "_CD4T_umap_density.png"),
         width = 10, height = 10, plot = p)
  p
})