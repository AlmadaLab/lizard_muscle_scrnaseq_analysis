sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/integration_v3_output.log", append = TRUE)
pdf("/project/aealmada_561/judyz/single_cell/anocar_v3/plot/integration_plots_v3.pdf", width = 10, height = 8)

# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)

# load RDS
quad_O_qc <- readRDS("RDS/quad_O_qc.rds")
tail_O_qc <- readRDS("RDS/tail_O_qc.rds")

set.seed(1234) 

# combine datasets without integration
combined <- merge(tail_O_qc, y = quad_O_qc, add.cell.ids = c("Tail", "Quads"), project = "Lizard")
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000, slot = "data")
combined <- ScaleData(combined, features = rownames(combined))
combined <- RunPCA(combined, verbose = FALSE, features = VariableFeatures(combined))


combined <- FindNeighbors(combined, dims = 1:30, reduction = "pca")
combined <- FindClusters(combined, resolution = 0.2)
combined <- RunUMAP(combined, reduction = 'pca', dims = 1:30, reduction.name = "umap_unintegrated")
combined <- RunTSNE(combined, reduction = 'pca', dims = 1:30, reduction.name = "tsne_unintegrated")


#print(DimPlot(combined, reduction = "pca", group.by = "orig.ident"))
#print(DimPlot(combined, reduction = "umap_unintegrated", group.by = "orig.ident"))
print(DimPlot(combined, reduction = "tsne_unintegrated", group.by = "orig.ident", pt.size = 0.5)+ ggtitle("tSNE before integration colored by samples"))

print(combined)




# add Harmony integration 
combined <- RunHarmony(combined, group.by.vars = "orig.ident")
#print(DimPlot(combined, reduction = "harmony", group.by = "orig.ident"))


combined <- FindNeighbors(combined, dims = 1:30, reduction = "harmony")
combined <- FindClusters(combined, resolution = 0.2)
combined <- RunUMAP(combined, reduction = 'harmony', dims = 1:30, reduction.name = "umap_integrated")
combined <- RunTSNE(combined, reduction = 'harmony', dims = 1:30, reduction.name = "tsne_integrated")

#print(DimPlot(combined, reduction = "umap.integrated", group.by = "orig.ident"))
print(DimPlot(combined, reduction = "tsne_integrated", group.by = "orig.ident", pt.size = 0.5) + ggtitle("tSNE after integration colored by samples"))
print(DimPlot(combined, reduction = "tsne_integrated", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) + ggtitle("tSNE after integration colored by cell types"))




# plot quads and tails on separate plots to compare cell cluster composition
tail_subset <- subset(combined, subset = orig.ident == "lizard_tail")
#print(DimPlot(tail_subset, reduction = "tsne_integrated", group.by = "orig.ident"))
print(DimPlot(tail_subset, reduction = "tsne_integrated", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) + ggtitle("tail cluster after integration"))


quad_subset <- subset(combined, subset = orig.ident == "lizard_quads")
#print(DimPlot(quad_subset, reduction = "tsne.integrated", group.by = "orig.ident"))
print(DimPlot(quad_subset, reduction = "tsne_integrated", group.by = "seurat_clusters", label = TRUE, pt.size = 0.5) + ggtitle("quads cluster after integration"))






if (!dir.exists("RDS")) {
  dir.create("RDS")
}

saveRDS(combined, file = "RDS/combined_after_integration.rds")

# close pdf
dev.off()
sink()