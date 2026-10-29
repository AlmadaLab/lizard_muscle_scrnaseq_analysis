sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(harmony)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


print("Run FindAllMarkers to get DEG per cell type")

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


marker_genes <- function(data) {
  # set to seurat clusters
  Idents(data) <- "Annotation_1"
  markers <- FindAllMarkers(data, only.pos = TRUE)
  filtered_markers <- markers %>%
    dplyr::filter(p_val_adj < 0.05,
                  avg_log2FC > 1) %>%
    arrange(desc(avg_log2FC))
  return (filtered_markers)
}

cluster_markers <- marker_genes(obj)
write.csv(cluster_markers, "tables/annotated_obj_anocarv3_FindAllMarkers_log2fc1.csv")
