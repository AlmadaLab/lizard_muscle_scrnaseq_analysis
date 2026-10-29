sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


print("Calculate %MT from AnoCar v2")


# Step 1: read in data and create seurat object
tail.matrix <- Read10X_h5("/project2/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/Tail/outs/filtered_feature_bc_matrix.h5")
quad.matrix <- Read10X_h5("/project2/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/Quads/outs/filtered_feature_bc_matrix.h5")

tail_obj <- CreateSeuratObject(counts = tail.matrix, project = "tail", min.cells = 3, min.features = 200)
quad_obj <- CreateSeuratObject(counts = quad.matrix, project = "quad", min.cells = 3, min.features = 200)

## merge the two data object together
merged_obj <- merge(
  x = tail_obj, 
  y = quad_obj, 
  add.cell.ids = c("tail", "quad"),
  project = "lizard")

print(merged_obj)


# Step 2: save the barcodes with high mt%
print("calculate mt%")
merged_obj[["percent.mt"]] <- PercentageFeatureSet(merged_obj, pattern = "^MT-")

print(summary(merged_obj$percent.mt))
print(summary(merged_obj$nCount_RNA))
print(summary(merged_obj$nFeature_RNA))
table(merged_obj$orig.ident)


# Step 3: QC
obj <- subset(merged_obj, subset = nFeature_RNA > 500 & percent.mt < 5)
table(obj$orig.ident)


seurat_mat <- NormalizeData(obj, normalization.method = "LogNormalize")
seurat_mat <- FindVariableFeatures(seurat_mat, selection.method='vst', nfeatures=2000)
seurat_mat <- ScaleData(seurat_mat)
seurat_mat <- RunPCA(seurat_mat, features = VariableFeatures(seurat_mat))
seurat_mat <- RunHarmony(seurat_mat, group.by.vars = "orig.ident")
seurat_mat <- FindNeighbors(seurat_mat, dims = 1:30, reduction = "harmony")
seurat_mat <- FindClusters(seurat_mat, resolution = 0.2)
seurat_mat <- RunTSNE(seurat_mat, reduction = 'harmony', dims = 1:30)

saveRDS(seurat_mat, "rds/integrated_obj_nfeat500_mt5_anocarV2.rds")