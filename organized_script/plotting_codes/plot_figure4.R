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

# Codes to generate plot in figure 4

####################################################################################

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



# plot heatmap of the top 3 surface markers sorted by log2fc
## top findallmarker heatmap: ranked by both log2fc and pct.1 > 0.4

genes <- c(
  "cmklr2", "wnt11", "nt5e",    # FAP1
  "postn", "rcn3", "col6a1",  # FAP2
  "itga8", "col12a1", "col11a1",    # FAP3
  "fgfr4", "calcr", "megf10",  # MuSC
  "hspb1", "pnck", "stac",  #Smooth muscle
  "plvap", "ehd3", "ptprb", # Endothelial
  "mmrn1", "piezo2", "hapln3", # LEC
  "ptch2", "bgn", "ptch1", #Neural
  "ugt8", "mpz", "septin3", # Oligo
  "lamp3", "csf1r", "ptafr", # Macrophages
  "il2rb","cd2", "il7r", # T
  "cd79b", "cd19", "pld4", # B
  "flt3", "arhgap28", "itgb7", # NK cells
  "trhr", "kit", "il2ra", # Bone marrow
  "atp8b1", "mpl", "cd226" # Megakaryocytes
)


# compute average accessibility per subtype
avg_mat <- AverageExpression(
  obj,
  assays = "RNA",
  features = genes,
  group.by = "Annotation_1",
  slot = "data"
)$RNA


# keep only genes that exist and reorder to match your vector
genes.use <- genes[genes %in% rownames(avg_mat)]
avg_mat <- avg_mat[genes.use, ]

# Z-scale across genes
avg_mat_scaled <- t(scale(t(avg_mat)))


# plot the heatmap
library(ComplexHeatmap)
library(circlize)

col_fun = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

p3 <- Heatmap(
  avg_mat_scaled,
  name = "Z-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 10),
  col = col_fun,
  rect_gp = grid::gpar(col = "grey70", lwd = 0.3)
)

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_4/Top_surface_markers_RBH_heatmap_080526.pdf", width = 8, height = 8)
draw(p3)
dev.off()




#################################################################################
# Plot the literature-supported marker genes featureplot

surface.genes <- c(
  "pdgfra",
  "postn",
  "itga8",
  "calcr",
  "csf1r",
  "flt1",
  "trpc6",
  "cd8a",
  "cd226",
  "ngfr",
  "cd180",
  "kit",
  "cldn11",
  "mmrn1",
  "itga2b"
)

p2 <- FeaturePlot(obj, features=surface.genes, reduction="tsne", ncol = 3, order=TRUE)

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_4/surface_markers_tsne_073026.pdf", width = 12, height = 14)
print(p2)
dev.off()

