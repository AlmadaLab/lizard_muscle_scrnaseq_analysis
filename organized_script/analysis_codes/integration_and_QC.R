sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(harmony)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


print("Integrate the tail and quad libraires mapped to AnoCar reference genome V3")


# Step 1: read in data and create seurat object
tail.matrix <- Read10X_h5("/project2/aealmada_561/judyz/single_cell/anocar_v3/cell_ranger/Tail/outs/filtered_feature_bc_matrix.h5")
quad.matrix <- Read10X_h5("/project2/aealmada_561/judyz/single_cell/anocar_v3/cell_ranger/Quads/outs/filtered_feature_bc_matrix.h5")

tail_obj <- CreateSeuratObject(counts = tail.matrix, project = "tail", min.cells = 3, min.features = 200)
quad_obj <- CreateSeuratObject(counts = quad.matrix, project = "quad", min.cells = 3, min.features = 200)

## merge the two data object together
merged_obj <- merge(
  x = tail_obj, 
  y = quad_obj, 
  add.cell.ids = c("tail", "quad"),
  project = "lizard")

print(merged_obj)





# Step 2: get the calculated percent.mt from anocar v2
anocar.v2 <- readRDS("rds/merged_obj_unfiltered_anocarV2.rds")

shared.barcodes <- intersect(
  colnames(merged_obj),
  colnames(anocar.v2)
)

length(shared.barcodes)

# Subset merged_obj to shared barcodes
merged_obj.shared <- subset(
  merged_obj,
  cells = shared.barcodes
)

# Add the percent.mt metadata column to the merged_obj from anocar.v2
merged_obj.shared$percent.mt <- anocar.v2$percent.mt[
  match(
    colnames(merged_obj.shared),
    colnames(anocar.v2)
  )
]


print(summary(merged_obj.shared$percent.mt))
print(summary(merged_obj.shared$nCount_RNA))
print(summary(merged_obj.shared$nFeature_RNA))
table(merged_obj.shared$orig.ident)


# Step 3: save the merged object (anocar v2 unfiltered)
saveRDS(merged_obj.shared, "rds/merged_obj_unfiltered_anocarV3.rds") # will be used for plotting the mitochondrial distribution


# Step 3: QC
obj <- subset(merged_obj.shared, subset = nFeature_RNA > 500 & percent.mt < 5)
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


# Step 9: add annotation to seurat clusters
annotation <- read.csv("/project2/aealmada_561/judyz/20260621_scRNA_lizard/rscripts/annotation.csv")
add_annotation<- function(data, ann, levels) {
  Idents(data) <- data$seurat_clusters
  # change numeric cluster to level1 labels
  data@active.ident <- plyr::mapvalues(
    x = data@active.ident,
    from = ann$Cluster,
    to = ann$Annotation
  )
  
  # add column of level1_round1 annotation
  data[[levels]] <- Idents(data)
  
  return(data)
}

seurat_mat <- add_annotation(seurat_mat, annotation, "Annotation_1")

# save the object

saveRDS(seurat_mat, "rds/integrated_obj_nfeat500_mt5_anocarV3.rds")




