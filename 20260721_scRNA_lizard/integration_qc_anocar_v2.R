sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


print("Integrate the tail and quad libraires mapped to AnoCar reference genome V2")


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


# Step 3: save the merged object (anocar v2 unfiltered)
saveRDS(merged_obj, "rds/merged_obj_unfiltered_anocarV2.rds") # will be used for plotting the mitochondrial distribution

# save the list of barcodes with percent.mt > 10, these would be removed
mt_5_barcodes <- Cells(merged_obj)[merged_obj$percent.mt > 5]
mt_10_barcodes <- Cells(merged_obj)[merged_obj$percent.mt > 10]

# Check how many
length(mt_5_barcodes)
length(mt_10_barcodes)

saveRDS(mt_5_barcodes, "rds/mt_barcodes/barcodes_greater_than_5.rds")
saveRDS(mt_10_barcodes, "rds/mt_barcodes/barcodes_greater_than_10.rds")


# Step 3: QC
obj <- subset(merged_obj, subset = nFeature_RNA > 500 & percent.mt < 5)
table(obj$orig.ident)


# Step 4: normalize and select variable features before integraion
print("start normalizing and scaling")
seurat_mat <- NormalizeData(obj, normalization.method = "LogNormalize")
seurat_mat <- FindVariableFeatures(seurat_mat, selection.method='vst', nfeatures=2000)


# Step 5: scale data and run PCA before integration
seurat_mat <- ScaleData(seurat_mat) # scale only the variable features
seurat_mat <- RunPCA(seurat_mat, features = VariableFeatures(seurat_mat))


# Step 6: apply harmony integration
seurat_mat <- RunHarmony(seurat_mat, group.by.vars = "orig.ident")


# Step 7: low-dimensional embeddings and clustering analysis, on the integrated reduction
seurat_mat <- FindNeighbors(seurat_mat, dims = 1:30, reduction = "harmony")
seurat_mat <- FindClusters(seurat_mat, resolution = 0.2)

seurat_mat <- RunUMAP(seurat_mat, reduction = 'harmony', dims = 1:30)
seurat_mat <- RunTSNE(seurat_mat, reduction = 'harmony', dims = 1:30)


# Step 8: for seurat v5, join the layers
seurat_mat[["RNA"]] <- JoinLayers(seurat_mat[["RNA"]])


# save the object

saveRDS(seurat_mat, "rds/integrated_obj_nfeat500_mt5_anocarV2.rds")
