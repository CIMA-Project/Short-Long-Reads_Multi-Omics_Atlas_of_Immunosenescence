## ------------------------------
## MAIT VJ PAIR Usage
## ------------------------------
combined.TCR <-readRDS("ASLM/TCR/TCR.rds")
df<- bind_rows(combined.TCR,.id="sample")%>%
  mutate(
    TRAV = str_split_fixed(TCR1, "\\.", 3)[,1],
    TRAJ = str_split_fixed(TCR1, "\\.", 3)[,2]
  )

vj_summary <- df %>%
  group_by(celltype_L3, TRAV, TRAJ) %>%
  summarise(count = n(), .groups = "drop")
sub_data <- vj_summary %>% filter(celltype_L3 == "MAIT")

all_genes <- unique(c(sub_data$TRAV, sub_data$TRAJ))
grid_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(all_genes))
names(grid_colors) <- all_genes
circos.clear()
par(mar = c(1, 1, 3, 1))
pdf("./pdf/supplemental/MAIT_VJ_Usage.pdf", width = 10, height = 10)
title("V–J Pair Usage in MAIT Cells", cex.main = 1.5)
chordDiagram(
  x = sub_data[, c("TRAV", "TRAJ", "count")],
  grid.col = grid_colors,
  transparency = 0.4,
  annotationTrack = "grid", 
  preAllocateTracks = list(
    list(track.height = 0.05)  
  )
)

circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim), 
      y = ylim[2]*1.1 , 
      labels = sector.index,
      facing = "clockwise",  
      niceFacing = TRUE,    
      adj = c(0.5, 0),       
      cex = 0.6,            
      col = "black"          
    )
  }, 
  bg.border = NA
)
dev.off()
