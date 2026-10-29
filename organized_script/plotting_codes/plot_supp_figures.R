sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")

# Codes to generate plot in supplementary figure



#################################################################################
# Supp fig1a: reference genome comparison
ref_v2 <- read.delim(
  "/project2/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/Quads/outs/filtered_feature_bc_matrix/features.tsv.gz",
  header = FALSE
)

ref_v3 <- read.delim(
  "/project2/aealmada_561/judyz/single_cell/anocar_v3/cell_ranger/Quads/outs/filtered_feature_bc_matrix/features.tsv.gz",
  header = FALSE
)

# Add column names
colnames(ref_v2) <- c("gene_id", "gene_name", "feature_type")
colnames(ref_v3) <- c("gene_id", "gene_name", "feature_type")

head(ref_v2)
head(ref_v3)


# Count how many genes have gene name in AnoCar v2
v2_named <- ref_v2 %>%
  dplyr::filter(gene_name != gene_id)

nrow(v2_named)


# Count how many genes have gene name in AnoCar v3
v3_named <- ref_v3 %>%
  dplyr::filter(!grepl("^LOC", gene_name))

nrow(v3_named)



# Count overlap: unique annotated genes
# Convert gene names to uppercase
v2_genes <- unique(toupper(v2_named$gene_name))
v3_genes <- unique(toupper(v3_named$gene_name))

# Shared genes
shared_genes <- intersect(v2_genes, v3_genes)

# Number of shared genes
length(shared_genes)

# V3 only
v3_only <- setdiff(v3_genes, v2_genes)
length(v3_only)

# V2 only
v2_only <- setdiff(v2_genes, v3_genes)
length(v2_only)





#################################################################################
# Supp fig1b: Cell level QC

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
# Supp fig4: cellchat outgoing signaling

obj.cellChat <- readRDS("rds/processed_cellchat_obj_reordered.rds")


mat <- obj.cellChat@net$weight
groupSize <- as.numeric(table(obj.cellChat@idents))



# define color codes
col <- c(
  "FAPs 1"           = "#5CC8D7",
  "FAPs 2"           = "#AADB84",
  "FAPs 3"           = "#D69C6C",
  "MuSCs"            = "#B294C7",
  "Smooth muscle cells" = "#569EA4",
  "Endothelial cells" = "#FE870D",
  "LECs"   = "#9966CC",
  "Neural/Glial cells"  = "#E6CB75",
  "Oligodendrocytes"   = "#C7B699",
  "Macrophages"      = "#F88A89",
  "T cells"        = "#F06C45",
  "B cells"         = "#4992C2",
  "Natural killer cells" = "#FF66CC",
  "Bone marrow stem cells" = "#52AF43",
  "Megakaryocytes"      = "#D2042D"
)



# the dot size correlate with the number of cells in that cell type
# all panels share the same edge-width scale, so arrow thickness is visually comparable across cell types.
pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/supp_figure/CellChat_outgoing_networks.pdf", width = 15, height = 10)

par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, color.use=col, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()



# investigate the Top 3 receiver cell types for each sender
# Top 3 receiver cell types for each sender
top_receivers <- do.call(
  rbind,
  lapply(rownames(mat), function(sender) {
    
    receiver_weights <- mat[sender, ]
    
    # Optional: exclude autocrine/self-signaling
    receiver_weights[sender] <- NA
    
    top3 <- sort(receiver_weights, decreasing = TRUE, na.last = NA)[1:3]
    
    data.frame(
      sender = sender,
      receiver = names(top3),
      communication_weight = as.numeric(top3),
      rank = seq_along(top3),
      row.names = NULL
    )
  })
)

top_receivers


write.csv(top_receivers, "/project2/aealmada_561/judyz/20260621_scRNA_lizard/tables/cellchat_top3_receivers.csv")







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


#################################################################################
# Supp fig6: cell type proportion

obj <- readRDS("rds/annotated_obj_anocarv3_nfeat500_mt5.rds")

prop.original <- obj@meta.data %>%
  count(orig.ident, Annotation_1) %>%
  group_by(orig.ident) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

head(prop.original)


pdf("figures/figure_6/tail_quad_actual_composition.pdf", width = 8, height = 6)

p1 <- ggplot(
  prop.original,
  aes(x = Annotation_1, y = proportion, fill = orig.ident)
) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  scale_fill_manual(values = c(
    tail = "#4C78A8",
    quad = "#F58518"
  )) +
  labs(
    x = NULL,
    y = "Cell proportion"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p1)
dev.off()
