sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")

# Codes to generate plot in supplementary figure

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

#################################################################################
# cell type level RNA QC metric
pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/supp_figure/celltype_QC_082026.pdf", width = 8, height = 4)
print(VlnPlot(obj, features = "nFeature_RNA", cols = col, group.by = "Annotation_1", pt.size = 0) + NoLegend())
print(VlnPlot(obj, features = "nCount_RNA", cols = col, group.by = "Annotation_1", pt.size = 0) + NoLegend())
print(VlnPlot(obj, features = "percent.mt",cols = col, group.by = "Annotation_1", pt.size = 0) + NoLegend())
dev.off()


#################################################################################
# number of cells per tissue/library


