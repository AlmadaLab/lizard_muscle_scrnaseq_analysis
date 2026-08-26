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
## conduct FindMarkers between tail and quads within the same cell type


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
# DEG Tail vs Quad within each matched cell type

celltypes <- unique(obj$Annotation_1)
DEG.list <- list()

for (ct in celltypes) {
  
  # Subset one cell type
  obj.sub <- subset(obj, subset = Annotation_1 == ct)
  
  # Set condition as identity
  Idents(obj.sub) <- "orig.ident"
  
  # Make sure both tail and quad are present
  if (all(c("tail", "quad") %in% levels(Idents(obj.sub)))) {
    
    DEG <- FindMarkers(
      obj.sub,
      ident.1 = "tail",
      ident.2 = "quad",
      test.use = "wilcox",
      logfc.threshold = 0.25,
      min.pct = 0.1
    )
    
    DEG$gene <- rownames(DEG)
    DEG$celltype <- ct
    
    DEG.list[[ct]] <- DEG
  }
}

DEG.all <- bind_rows(DEG.list)

# Add direction
DEG.filtered <- DEG.all %>%
  mutate(
    direction = ifelse(
      avg_log2FC > 0,
      "Tail",
      "Quad"
    )
  )

# Filter: apply very stringent threshold
DEG.filtered <- DEG.filtered %>%
  filter(
    p_val_adj < 0.01,
    abs(avg_log2FC) > 2,
    celltype == "FAPs 2"
  )


musc_subset <- subset(obj, idents = "MuSCs")

