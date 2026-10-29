sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


# Codes to generate plot in figure 3


#################################################################################
# tSNE plot with cell type annotated

obj <- readRDS("rds/annotated_obj_anocarv3_nfeat500_mt5.rds")

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


# define color codes
col <- c(
  "FAPs 1"           = "#5CC8D7",
  "FAPs 2"           = "#AADB84",
  "FAPs 3"           = "#D69C6C",
  "MuSCs"            = "#B294C7",
  "Macrophages"      = "#F88A89",
  "T cells"        = "#F06C45",
  "B cells"         = "#4992C2",
  "Natural killer cells" = "#FF66CC",
  "Bone marrow stem cells" = "#52AF43",
  "Megakaryocytes"      = "#D2042D",
  "Endothelial cells" = "#FE870D",
  "LECs"   = "#9966CC",
  "Smooth muscle cells" = "#569EA4",
  "Neural/Glial cells"  = "#E6CB75",
  "Oligodendrocytes"   = "#C7B699"
)



p1 <- DimPlot(obj, reduction="tsne", cols = col, group.by="Annotation_1", label=T, raster=F, order=T) + NoLegend()
p1

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_3/celltype_tsne_072626.pdf", width = 7, height = 6)
print(p1)
dev.off()


#################################################################################
# dotplot marker genes (literature supported)

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


Idents(obj) <- "Annotation_1"

genes <- c(
  "pdgfra", "mfap5",
  "pcolce2", "postn",
  "wnt5a", "lox",
  "pax7", "myf5",
  "lmod1", "acta2",
  "flt1", "pecam1", 
  "prox1", "mmrn1",
  "prickle1", "ptn",
  "plp1", "bcas1",
  "marco", "csf1r",
  "skap1","batf",
  "cd19", "pax5",
  "cd244", "flt3",
  "il2ra", "kit",
  "itga2b", "tbxas1"
)


p2 <- DotPlot(
  obj,
  features = genes,
  group.by = "Annotation_1"
)+ RotatedAxis() +
  scale_y_discrete(limits = rev(celltype.order))


p2



pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_3/marker_dotplot_072626.pdf", width = 10, height = 6)
print(p2)
dev.off()



#################################################################################
# top findallmarker heatmap: ranked by both log2fc and pct.1 > 0.4

genes <- c(
  "edn3", "mfap5", "wnt11",    # FAP1
  "cryaa", "lag3", "pcolce2",  # FAp1
  "fmod", "comp", "ptx4",    # FAP3
  "c1qtnf3", "myf5", "msc",  # MuSC
  "hspb1", "jph2", "colec10",  #Smooth muscle
  "sox7", "plvap", "dipk2b", # Endothelial
  "stab2", "tbx1", "myo7a", # LEC
  "krt18", "cldn1", "sfn", #Neural
  "ugt8", "mpz", "septin3", # Oligo
  "marco", "acod1", "acp3", # Macrophages
  "sh2d1b","il2rb", "ctsw", # T
  "cd79b", "pax5", "cd19", # B
  "flt3", "ubxn11", "arap2", # NK cells
  "rgs13", "trhr", "fgf21", # Bone marrow
  "gp1bb", "tbxas1", "tubb1" # Megakaryocytes
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

col_fun = colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "white", "#B2182B")
)

p3 <- Heatmap(
  t(avg_mat_scaled),  # transpose matrix
  name = "Z-score",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  row_names_side = "left",
  row_names_gp = grid::gpar(fontsize = 10),
  column_names_gp = grid::gpar(fontsize = 10),
  col = col_fun,
  rect_gp = grid::gpar(col = "grey70", lwd = 0.3)
)

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_3/Top_FindAllMarkers_heatmap_072626.pdf", width = 12, height = 5)
draw(p3)
dev.off()
