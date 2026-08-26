sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/CellChat_v3_output.log", append = TRUE)
pdf("/project/aealmada_561/judyz/single_cell/anocar_v3/plot/CellChat_object_v3.pdf", width = 10, height = 8)

# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(svglite)
library(CellChat)


set.seed(1234) 



# load RDS
quad_O_qc <- readRDS("RDS/quad_O_qc.rds")
tail_O_qc <- readRDS("RDS/tail_O_qc.rds")
joint_combined <- readRDS("RDS/joint_combined_after_integration.rds")
combined.markers <- readRDS("RDS/combined.markers.rds")
ident_joint_combined <- readRDS("RDS/ident_joint_combined.rds")



#######################################################################################
# load from Seurat
seurat.input <- ident_joint_combined[["RNA"]]$data # normalized data matrix
labels <- Idents(ident_joint_combined)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels


#######################################################################################
# create cellchat object
cellChat <- createCellChat(object = ident_joint_combined, group.by = "ident", assay = "RNA")

## extra step: convert to uppercase
rownames(cellChat@data) <- toupper(rownames(cellChat@data))

#######################################################################################
# select cellchat database
CellChatDB <- CellChatDB.human

# select subset of database for cell chat: only protein signaling
CellChatDB.use <- subsetDB(CellChatDB)
cellChat@DB <- CellChatDB.use

#######################################################################################
# Preprocessing the expression data
cellChat <- subsetData(cellChat)
cellChat <- identifyOverExpressedGenes(cellChat, do.fast = FALSE)
cellChat <- identifyOverExpressedInteractions(cellChat)
cellChat <- smoothData(cellChat, adj = PPI.human)


#######################################################################################
# Compute communication probability and infer network
cellChat <- computeCommunProb(cellChat, type = "triMean", raw.use = FALSE)  # default is triMean (25% truncated)

# compute communication probability on signaling pathway level
cellChat <- computeCommunProbPathway(cellChat)

# compute aggregated cell-cell communication network
cellChat <- aggregateNet(cellChat)
groupSize <- as.numeric(table(cellChat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
#print(netVisual_circle(cellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions"))
print(netVisual_circle(cellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength"))




saveRDS(cellChat, file = "RDS/cellChat.rds")

# close pdf
dev.off()

sink()