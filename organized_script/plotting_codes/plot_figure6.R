library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(readr)
library(enrichR)
library(cowplot)

set.seed(1) 
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


####################################################################################

# For figure 6: 
## generate scatter plot and marker gene expression for the MuSCs between tail and quad


# Step 1: read in the object and define the customized order
obj <-readRDS("rds/annotated_obj_anocarv3_nfeat500_mt5.rds")

# define order
celltype.order <- c(
  "FAPs 1",
  "FAPs 2",
  "FAPs 3",
  "MuSCs",
  "Smooth muscle cells",
  "Endothelial cells",
  "LECs",
  "Neural/Glial cells",
  "Oligodendrocytes",
  "Macrophages",
  "T cells",
  "B cells",
  "Natural killer cells",
  "Bone marrow stem cells",
  "Megakaryocytes"
)
obj$Annotation_1 <- factor(
  obj$Annotation_1,
  levels = celltype.order
)

####################################################################################
# Pearson correlation of average expression profile between tail and quad

## Make combined cell type + condition identity
obj$celltype_condition <- paste(
  obj$Annotation_1,
  obj$orig.ident,
  sep = "_"
)

avg.exp <- AverageExpression(
  obj,
  assays = "RNA",
  group.by = "celltype_condition",
  slot = "data"
)$RNA

colnames(avg.exp)

## separate tail and quad
tail.cols <- paste0(celltype.order, "-tail")
quad.cols <- paste0(celltype.order, "-quad")
tail.avg <- as.matrix(avg.exp[, tail.cols])
quad.avg <- as.matrix(avg.exp[, quad.cols])

## calculate the Pearson correlation
cor.mat <- cor(
  tail.avg,
  quad.avg,
  method = "pearson"
)

# Clean names
rownames(cor.mat) <- celltype.order
colnames(cor.mat) <- celltype.order

library(ComplexHeatmap)
library(circlize)

col_fun <- colorRamp2(
  c(0.3, 0.65, 1),
  c("#2166AC", "white", "#B2182B")
)

p1 <- Heatmap(
  cor.mat,
  name = "Pearson r",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_title = "Tail",
  column_title = "Quad",
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%.2f", cor.mat[i, j]),
      x, y,
      gp = gpar(fontsize = 8)
    )
  }
)

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_6/heatmap_Pearson_tail_quad_average_exp.pdf", 
    width = 8, height = 7)
draw(p1)
dev.off()




####################################################################################
# Plot MuSC Pearson correlation

# Extract MuSC average expression
df.musc <- data.frame(
  quad = avg.exp[, "MuSCs-quad"],
  tail = avg.exp[, "MuSCs-tail"]
)

# Pearson correlation
r <- cor(
  df.musc$quad,
  df.musc$tail,
  method = "pearson"
)

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_6/musc_tail_quad_scatter_Pearson.pdf", 
    width = 6, height = 5)

p2 <- ggplot(df.musc, aes(x = quad, y = tail)) +
  geom_point(size = 0.8, alpha = 0.4) +
  theme_classic() +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "#D55E5E"
  ) +
  labs(
    title = "MuSCs",
    subtitle = paste0("Pearson r = ", round(r, 3)),
    x = "Quad average expression",
    y = "Tail average expression"
  )
print(p2)

dev.off()


####################################################################################
# Violin plot on MuSC marker genes
pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_6/musc_markers_across_tissues_vlnplot.pdf", width = 8, height = 4)

print(VlnPlot(obj, cols = c(tail = "#4C78A8",quad = "#F58518"), features = "pax7", group.by = "Annotation_1", split.by = "orig.ident", pt.size=0))
print(VlnPlot(obj, cols = c(tail = "#4C78A8",quad = "#F58518"), features = "calcr", group.by = "Annotation_1", split.by = "orig.ident", pt.size=0))
print(VlnPlot(obj, cols = c(tail = "#4C78A8",quad = "#F58518"), features = "myf5", group.by = "Annotation_1", split.by = "orig.ident", pt.size=0))
print(VlnPlot(obj, cols = c(tail = "#4C78A8",quad = "#F58518"), features = "myod1", group.by = "Annotation_1", split.by = "orig.ident", pt.size=0))

dev.off()