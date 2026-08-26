library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(CellChat)
library(Matrix)
set.seed(1) 
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


# read in seurat object in human orthologs
obj <- readRDS("rds/annotated_obj_in_human_orthologs.rds")

# Step 1: create cellchat object from seurat object
data.input <- LayerData(obj, assay = "RNA", layer = "data")

labels <- Idents(obj)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels


obj.cellChat <- createCellChat(object = obj, group.by = "Annotation_1", assay = "RNA")

## change the order
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

obj.cellChat@idents <- factor(
  obj.cellChat@idents,
  levels = celltype.order
)

levels(obj.cellChat@idents)


# Step 2: set to human database
CellChatDB <- CellChatDB.human
unique(CellChatDB$interaction$annotation)

# Step 3: subset to everything except non-protein signaling
CellChatDB.use <- subsetDB(CellChatDB)
unique(CellChatDB.use$interaction$annotation)

# Set database
obj.cellChat@DB <- CellChatDB.use



# Step 4: process the single-cell expression data
obj.cellChat <- subsetData(obj.cellChat) # subset the expression data of signaling genes for saving computation cost

obj.cellChat <- identifyOverExpressedGenes(
  obj.cellChat,
  do.fast = FALSE
)

# The number of highly variable ligand-receptor pairs used for signaling inference is 1226 
obj.cellChat <- identifyOverExpressedInteractions(obj.cellChat)

# Project gene expression onto human PPI network: Apply the smoothing to reduce the dropout effects
obj.cellChat <- smoothData(obj.cellChat,adj = PPI.human)

## drop any unknown levels
obj.cellChat@idents <- droplevels(obj.cellChat@idents)


# Step 5: infer cell-cell communication network,  individual ligand–receptor interactions
obj.cellChat <- computeCommunProb(obj.cellChat, type = "triMean", raw.use = FALSE) # use the projected data
obj.cellChat@net



# Step 5: infer cell-cell communication network,  signaling pathway level
obj.cellChat <- computeCommunProbPathway(obj.cellChat, pairLR.use=NULL)
obj.cellChat@netP


# Step 6: aggregate cell-cell communication networks by summarizing the communication weights/probabilities between cell types
obj.cellChat <- aggregateNet(obj.cellChat)
# head(obj.cellChat@net$count) # show the number of interactions between any two cell groups

# head(obj.cellChat@net$weight) # show the matrix containing wieghts between any two cell groups


# save the cellchat object
saveRDS(obj.cellChat, "rds/processed_cellchat_obj_reordered.rds")

