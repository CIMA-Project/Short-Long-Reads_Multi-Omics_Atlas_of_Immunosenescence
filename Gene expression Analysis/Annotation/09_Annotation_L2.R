pbmc <- readRDS("/data/RDS/Annotation/4_PBMC_annotation_L1.rds")
B <- readRDS("/data/RDS/Annotation/5_B_annotation.rds")
Myeloid <- readRDS("/data/RDS/Annotation/8_Myeloid_annotation.rds")
NK <- readRDS('/data/work/Annotation/RDS_Files/T&NK/7_NK_L2_annotation.rds')
CD8T <- readRDS('/data/work/Annotation/RDS_Files/T&NK/9_CD8T_otherT_L3_annotation.rds')
CD4T <- readRDS('/data/work/Annotation/RDS_Files/T&NK/9_CD4T_L3_annotation.rds')
HSPC <- subset(pbmc, subset = celltype == "HSPC")

Idents(HSPC) = "celltype"
Idents(pbmc) = "L3_celltype"
Idents(B) = "L3_celltype"
Idents(Myeloid) = "L3_celltype"
Idents(NK) = "L3_celltype"
Idents(CD8T) = "L3_celltype"
Idents(CD4T) = "L3_celltype"

Idents(pbmc, cells = colnames(B)) <- Idents(B)
Idents(pbmc, cells = colnames(Myeloid)) <- Idents(Myeloid)
Idents(pbmc, cells = colnames(NK)) <- Idents(NK)
Idents(pbmc, cells = colnames(CD8T)) <- Idents(CD8T)
Idents(pbmc, cells = colnames(CD4T)) <- Idents(CD4T)
Idents(pbmc, cells = colnames(HSPC)) <- Idents(HSPC)
pbmc$L3_celltype <- Idents(pbmc)
saveRDS(pbmc, "/data/RDS/Annotation/10_PBMC_annotation_L3.rds")

# Figure1A short-read scRNA UMAP
L3_color_map <- c(
  `HSPC` = "#FAD7D7",`cMono` = "#CF4457",`intMono` = "#DE6A69",`ncMono` = "#EA5958",`cDC2` = "#F29667",`cDC1` = "#EF8787",`pDC` = "#FBC5AB",`MK` = "#F6BABB",
    `Swit. Bm` = "#D48CB1",`Unswit. Bm` = "#8085C0",`Naïve B` = "#D4A7C4",`Plasma` = "#DBC6D9",`Plsamablast` = "#E6E6F0",
    `NKbright` = "#F07D50",`NKdim` = "#FBAA34",`XCL+NKdim` = "#FEC928",`Terminal NK` = "#F4EB1B",
    `CD4 Naïve` = "#015696",`CD4 Regulatory` = "#89CAEA",`CD4 CTL` = "#BBE6FA",`CD4 Tcm` = "#0B75B3", `CD4 Tem` = "#4596CD",   
    `NKT` = "#87D649",`MAIT` = "#47A183",`γδT` = "#24858E",`Cycling T` = "#77C3BC",
    `CD8 CTL` = "#C3E023",`CD8 Naïve` ="#8DCA9D" ,`CD8 Tem` = "#2BB17E",`CD8 Tcm` = "#50C56A"
)
DimPlot(pbmc, reduction = "umap", order = FALSE,
        shuffle = TRUE, pt.size = 0.5, label = TRUE, cols = L3_color_map)

umap = pbmc@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(cell_type = pbmc@meta.data$L3_celltype)
head(umap)
p2 <- ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  geom_point(size = 0.1, alpha = 1) +
  scale_color_manual(values = L3_color_map) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    legend.key = element_blank()
)
p2
ggsave(
  filename = "/data/Figure/Annotation/1_PBMC_UMAP_F1D.png",
  plot = p2, width = 6, height = 6, dpi = 600
)
