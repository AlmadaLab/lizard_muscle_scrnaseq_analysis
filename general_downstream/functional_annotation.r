sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/figure3_v3_output.log", append = TRUE)
pdf("/project/aealmada_561/judyz/single_cell/anocar_v3/plot/figure3_plots_v3.pdf", width = 10, height = 8)

# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(svglite)


set.seed(1234) 





# load RDS
quad_O_qc <- readRDS("RDS/quad_O_qc.rds")
tail_O_qc <- readRDS("RDS/tail_O_qc.rds")
joint_combined <- readRDS("RDS/joint_combined_after_integration.rds")
combined.markers <- readRDS("RDS/combined.markers.rds")



#######################################################################################

# assign cell names to cluster (default cluster order)
# cluster 0: FAPs 1
# cluster 1: MuSCs
# cluster 2: Macrophages
# cluster 3: Endothelial cells
# cluster 4: FAPs 2
# cluster 5: FAPs 3
# cluster 6: Smooth muscle cells
# cluster 7: T-cells
# cluster 8: Natural killer cells
# cluster 9: Neural/Glial cells
# cluster 10: B-cells
# cluster 11: Bone marrow stem cells
# cluster 12: Oligodendrocytes
# cluster 13: LECs
# cluster 14: Megakaryocytes

annotate.ids <- c("FAPs 1", "MuSCs", "Macrophages", "Endothelial cells", "FAPs 2", "FAPs 3",
    "Smooth muscle cells", "T-cells", "Natural killer cells", "Neural/Glial cells",
    "B-cells", "Bone marrow stem cells", "Oligodendrocytes", "LECs", "Megakaryocytes")
names(annotate.ids) <- levels(joint_combined)
joint_combined <- RenameIdents(joint_combined, annotate.ids)



# set custom order to make three FAPs close to each other
Idents(joint_combined) <- factor(Idents(joint_combined),
  levels = c("FAPs 1", "FAPs 2", "FAPs 3", "MuSCs", "Macrophages", "Endothelial cells", 
             "Smooth muscle cells", "T-cells", "Natural killer cells", 
             "Neural/Glial cells", "B-cells", "Bone marrow stem cells", 
             "Oligodendrocytes", "LECs", "Megakaryocytes"))


print(DimPlot(joint_combined, reduction = "tsne_integrated", label = TRUE, pt.size = 0.5) + NoLegend())


#######################################################################################
# add custom label to DE table

## pair up biological label with Seurat default label
cluster_to_label <- c(
  "0" = "FAPs 1",
  "1" = "MuSCs",
  "2" = "Macrophages",
  "3" = "Endothelial cells",
  "4" = "FAPs 2",
  "5" = "FAPs 3",
  "6" = "Smooth muscle cells",
  "7" = "T-cells",
  "8" = "Natural killer cells",
  "9" = "Neural/Glial cells",
  "10" = "B-cells",
  "11" = "Bone marrow stem cells",
  "12" = "Oligodendrocytes",
  "13" = "LECs",
  "14" = "Megakaryocytes"
)
## add to DE table
combined.markers$custom_cluster <- cluster_to_label[as.character(combined.markers$cluster)]

custom_order <- c("FAPs 1", "FAPs 2", "FAPs 3", "MuSCs", "Macrophages", "Endothelial cells", 
                  "Smooth muscle cells", "T-cells", "Natural killer cells", 
                  "Neural/Glial cells", "B-cells", "Bone marrow stem cells", 
                  "Oligodendrocytes", "LECs", "Megakaryocytes")
## change the order
combined.markers$custom_cluster <- factor(combined.markers$custom_cluster, levels = custom_order)


#######################################################################################
# visualize each marker using heatmap

# top 3 genes for each cluster (exclude unannotated genes)
top_a <- combined.markers %>%
  filter(avg_log2FC > 1, !grepl("^LOC", gene)) %>%
  group_by(custom_cluster) %>%
  slice_head(n = 3) %>%
  ungroup()

print(DoHeatmap(object = joint_combined, features = top_a$gene, size = 4, draw.lines = TRUE, raster = FALSE))




#######################################################################################
# visualize cell surface markers using heatmap

# get two surface markers for each cluster
surface_markers <- c("PDGFRA", "CD248",
         "ITGA11", "POSTN", 
         "ITGA8", "TNC", 
          "CALCR", "ITGA7",
          "CD68", "CSF1R",
          "FLT1", "PECAM1",
          "ANO1", "TRPC6", 
          "CD8A", "IL2RB", 
          "CD244", "CD226",
          "NGFR", "JAM3", 
          "CD79b", "CD79a",
          "KIT", "SELPLG", 
          "CLDN11", "MBP", 
          "MMRN1", "FLT4", 
          "ITGA2B", "MPL")

# plot tSNE with the top 10 cell surface markers
for (surface in tolower(surface_markers)) {
  print(FeaturePlot(joint_combined, features = surface, reduction = "tsne_integrated", slot = "scale.data", pt.size = 0.5))
}

# Plot heatmap with the top 10 cell surface markers
print(DoHeatmap(
  object = joint_combined,
  features = tolower(surface_markers),
  group.by = "ident",
  size = 4,
  draw.lines = TRUE,
  raster=FALSE
))

## downsample to make equal cell numbers for heatmap visualization
downsampled_cells <- colnames(subset(joint_combined, downsample = 50))
print(DoHeatmap(joint_combined, features = tolower(surface_markers), cells = downsampled_cells, raster=FALSE))



#######################################################################################
# visualize each marker using dotplot

# visualize each marker genes
marker_genes <- c("MFAP5", "PDGFRA",
                "PCOLCE2", "POSTN",
               "lox", "WNT5A",
               "MYF5", "PAX7",
               "CSF1R", "marco",
               "flt1", "pecam1",
               "ACTA2", "Lmod1",
               "SKAP1", "BATF",
               "FLT3", "cd244",
               "PTN", "prickle1",
               "PAX5", "cd19",
               "KIT", "IL2RA",
               "BCAS1", "plp1",
               "PROX1", "mmrn1",
               "ITGA2B", "Tbxas1")

## need to resort clusters creating the organized dotplot
print(DotPlot(joint_combined, features = tolower(marker_genes)) + RotatedAxis())

## switch the axis
print(DotPlot(joint_combined, features = tolower(marker_genes)) + RotatedAxis() + coord_flip())




# save the labeled data
ident_joint_combined <- joint_combined

if (!dir.exists("RDS")) {
  dir.create("RDS")
}

saveRDS(combined.markers, file = "RDS/combined.markers.rds")
saveRDS(ident_joint_combined, file = "RDS/ident_joint_combined.rds")


# close pdf
dev.off()

sink()