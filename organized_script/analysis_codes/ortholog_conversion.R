library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(CellChat)
library(Matrix)
set.seed(1) 
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


obj <-readRDS("rds/annotated_obj_anocarv3_nfeat500_mt5.rds")

# Step 1: extract counts and normalized data
counts <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "counts"
)

data <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)

dim(counts)
dim(data)

identical(rownames(counts), rownames(data))



# Step 2: match lizard genes to human orthologs
ortholog <- read.csv("database/human_anole_union_RBH_filtered.csv")

idx <- match(
  tolower(rownames(counts)),
  tolower(ortholog$anole)
)

# Genes that have a human ortholog
keep <- !is.na(idx)

# Get corresponding human gene symbols
human_genes <- ortholog$human[idx[keep]]

# Remove missing/empty human symbols
valid <- !is.na(human_genes) & human_genes != ""

keep_idx <- which(keep)[valid]
human_genes <- human_genes[valid]




# Step 3: convert both layers
counts_human <- counts[keep_idx, , drop = FALSE]
data_human   <- data[keep_idx, , drop = FALSE]

rownames(counts_human) <- human_genes
rownames(data_human)   <- human_genes


# Check
identical(
  rownames(counts_human),
  rownames(data_human)
)



# Step 5: create the new seurat object in human ortholog
human_ortho_obj <- CreateSeuratObject(
  counts = counts_human,
  meta.data = obj@meta.data
)

# Add the original normalized data
human_ortho_obj[["RNA"]]$data <- data_human


# save the object: 
saveRDS(human_ortho_obj, "rds/annotated_obj_in_human_orthologs.rds")
