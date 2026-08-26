sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/figure6_v3_output.log", append = TRUE)
pdf("/project/aealmada_561/judyz/single_cell/anocar_v3/plot/general_difference_plots.pdf", width = 10, height = 8)

# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(svglite)
library(scales)
library(DESeq2)
library(enrichR)

set.seed(1234) 


# load RDS
ident_joint_combined <- readRDS("RDS/ident_joint_combined.rds")
tail_markers <- readRDS("RDS/tail_markers.rds")
quads_markers <- readRDS("RDS/quads_markers.rds")


#######################################################################################
# calculate cell types proportion of each sample
# add cell type labels
ident_joint_combined$cell_type <- Idents(ident_joint_combined)

proportions <- ident_joint_combined@meta.data %>%
  group_by(orig.ident, cell_type) %>%
  summarize(count = n()) %>%  # For each group (sample + cluster), calculates the number of rows (cells)
  group_by(orig.ident) %>%
  mutate(proportion = count / sum(count)) %>% #Creates a new column, proportion, that divides the count for each cluster by the total count of cells for the corresponding sample (sum(count))
  ungroup()

print(proportions)

# extract seurat color codes
cell_types <- levels(Idents(ident_joint_combined))
seurat_colors <- hue_pal()(length(cell_types))  # Generate the default Seurat color palette
names(seurat_colors) <- cell_types  # assign cell types to their Seurat colors


# plot the cell proportion

print(ggplot(proportions, aes(x = orig.ident, y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", color = "black") + # Use "fill" to stack bars proportionally
  labs(
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Cluster",
    title = "Proportion of Cells in Each Cluster by Sample"
  ) +
  scale_fill_manual(values = seurat_colors) +
  theme_classic()
)


#######################################################################################
## visualize tSNE distribution
print(DimPlot(ident_joint_combined, reduction = "tsne_integrated", group.by = "orig.ident", pt.size = 0.5))



#######################################################################################
## compare pseudobulk profile of tail and quads

ident_joint_combined$cell_type <- Idents(ident_joint_combined)
Idents(ident_joint_combined) <- "orig.ident" 

aggr <- AggregateExpression(ident_joint_combined, assays = "RNA", slot = "counts", return.seurat = TRUE, group.by = c("orig.ident", "cell_type"))


print(CellScatter(aggr, "lizard-quads_FAPs 1", "lizard-tail_FAPs 1"))
print(CellScatter(aggr, "lizard-quads_MuSCs", "lizard-tail_MuSCs"))
print(CellScatter(aggr, "lizard-quads_Macrophages", "lizard-tail_Macrophages"))
print(CellScatter(aggr, "lizard-quads_Endothelial cells", "lizard-tail_Endothelial cells"))
print(CellScatter(aggr, "lizard-quads_FAPs 2", "lizard-tail_FAPs 2"))
print(CellScatter(aggr, "lizard-quads_FAPs 3", "lizard-tail_FAPs 3"))
print(CellScatter(aggr, "lizard-quads_Smooth muscle cells", "lizard-tail_Smooth muscle cells"))
print(CellScatter(aggr, "lizard-quads_T-cells", "lizard-tail_T-cells"))
print(CellScatter(aggr, "lizard-quads_Natural killer cells", "lizard-tail_Natural killer cells"))
colnames(aggr) <- gsub("/", ".", colnames(aggr))
print(CellScatter(aggr, "lizard-quads_Neural.Glial cells", "lizard-tail_Neural.Glial cells"))
print(CellScatter(aggr, "lizard-quads_B-cells", "lizard-tail_B-cells"))
print(CellScatter(aggr, "lizard-quads_Bone marrow stem cells", "lizard-tail_Bone marrow stem cells"))
print(CellScatter(aggr, "lizard-quads_LECs", "lizard-tail_LECs"))
print(CellScatter(aggr, "lizard-quads_Megakaryocytes", "lizard-tail_Megakaryocytes"))


#######################################################################################
# subset to MuSCs
# change cluster labels to name

Idents(ident_joint_combined) <- ident_joint_combined$cell_type

musc_subset <- subset(ident_joint_combined, idents = "MuSCs")
print(RidgePlot(musc_subset, features = "pax7", group.by = "orig.ident"))
print(RidgePlot(musc_subset, features = "myf5", group.by = "orig.ident"))
print(RidgePlot(musc_subset, features = "myod1", group.by = "orig.ident"))
print(RidgePlot(musc_subset, features = "calcr", group.by = "orig.ident"))
print(RidgePlot(musc_subset, features = "notch3", group.by = "orig.ident"))
print(RidgePlot(musc_subset, features = "bmp5", group.by = "orig.ident"))

# close pdf
dev.off()

sink()