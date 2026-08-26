library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(CellChat)
library(Matrix)
library(tidyr)
set.seed(1) 
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")


# read the cellchat object
obj.cellChat <- readRDS("rds/processed_cellchat_obj_reordered.rds")


levels(obj.cellChat@idents)

#######################################################################################
# for tables: extract all significant LR pairs and pathways
LR.sig <- subsetCommunication(obj.cellChat)
write.csv(LR.sig, "tables/cellChat_all_pathways.csv")


#######################################################################################
# Summarize the overall signaling pathways on cell level
communication <- data.frame(
  celltype = rownames(obj.cellChat@net$weight),
  Outgoing = rowSums(obj.cellChat@net$weight),
  Incoming = colSums(obj.cellChat@net$weight)
)


communication.long <- communication %>%
  pivot_longer(
    cols = c(Outgoing, Incoming),
    names_to = "Direction",
    values_to = "Strength"
  )

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

communication.long$celltype <- factor(
  communication.long$celltype,
  levels = rev(celltype.order)
)

p1 <-ggplot(
      communication.long,
      aes(x = celltype, y = Strength, fill = Direction)
    ) +
      geom_col(position = "dodge") +
      coord_flip() +
      scale_fill_manual(
        values = c(
          "Outgoing" = "#E67F7D",
          "Incoming" = "#78A6C8"
        )
      ) +
      theme_classic() +
      labs(
        x = NULL,
        y = "Aggregated communication strength",
        title = "Incoming and outgoing communication strength"
      )

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_5/sum_overall_signalings.pdf", width = 7, height = 6)
print(p1)
dev.off()


#######################################################################################
# relative contribution of each LR pair in each signaling pathway

# select three pathways for visualization
pathways <- c("CALCR", "NOTCH")

p2 <- netAnalysis_contribution(obj.cellChat,signaling = "CALCR")+
  geom_text(
    aes(label = paste0(round(contribution * 100, 1), "%")),
    hjust = -0.1,
    size = 4
  )


p3 <- netAnalysis_contribution(obj.cellChat, signaling = "NOTCH")+
  geom_text(
    aes(label = paste0(round(contribution * 100, 1), "%")),
    hjust = -0.1,
    size = 4
  )


pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_5/CALCR_NOTCH_relative_contribution.pdf", width = 6, height = 4)
print(p2)
print(p3)
dev.off()



#######################################################################################
# CALCR signaling pathway: chord diagram show the pathway communication between selected cell types

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_5/CALCR_chord_diagram.pdf", width = 10, height = 6)

# select cell types based on the rank of sender and receiver strength
calcr.df <- subsetCommunication(
  obj.cellChat,
  signaling = "CALCR"
)


# Major CALCR senders: FAPs 1, FAPs 2, Neural/Glial cells, FAPs 3
calcr.df %>%
  group_by(source) %>%
  summarise(strength = sum(prob)) %>%
  arrange(desc(strength))

# Major CALCR receivers: MuSCs, LECs, Endothelial cells, FAPs 3, Natural killer cells, 
calcr.df %>%
  group_by(target) %>%
  summarise(strength = sum(prob)) %>%
  arrange(desc(strength))


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

netVisual_chord_gene(
  obj.cellChat,
  signaling = "CALCR",
  color.use = col,
  sources.use = c("FAPs 1", "FAPs 2", "Neural/Glial cells", "FAPs 3"),
  targets.use = c("MuSCs", "LECs", "Endothelial cells", "FAPs 3")
)
dev.off()


#######################################################################################
# NOTCH signaling pathway: chord diagram show the pathway communication between selected cell types
pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_5/NOTCH_chord_diagram.pdf", width = 10, height = 6)

# select cell types based on the rank of sender and receiver strength
notch.df <- subsetCommunication(
  obj.cellChat,
  signaling = "NOTCH"
)


# Major CALCR senders: Endothelial cells, LECs, Smooth muscle cells
notch.df %>%
  group_by(source) %>%
  summarise(strength = sum(prob)) %>%
  arrange(desc(strength))

# Major CALCR receivers: Smooth muscle cells, MuSCs, LECs
notch.df %>%
  group_by(target) %>%
  summarise(strength = sum(prob)) %>%
  arrange(desc(strength))


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

netVisual_chord_gene(
  obj.cellChat,
  signaling = "NOTCH",
  color.use = col,
  sources.use = c("Endothelial cells", "LECs", "Smooth muscle cells"),
  targets.use = c("Smooth muscle cells", "MuSCs", "LECs")
)

dev.off()


#######################################################################################
# identify signaling roles of cell groups
# 2026.8.10 need to fix the bug

obj.cellChat <- netAnalysis_computeCentrality(
  obj.cellChat,
  slot.name = "netP"
)


pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_5/CALCR_NOTCH_signaling_roles_heatmap.pdf", width = 10, height = 6)

outgoing <- netAnalysis_signalingRole_heatmap(obj.cellChat, signaling = c("CALCR", "NOTCH"), color.use = col, pattern = "outgoing")
incoming <- netAnalysis_signalingRole_heatmap(obj.cellChat, signaling = c("CALCR", "NOTCH"), color.use = col, pattern = "incoming")
print(outgoing + incoming)

dev.off()







#######################################################################################
# Plot the signaling gene expression distribution
obj <- readRDS("rds/annotated_obj_in_human_orthologs.rds")

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

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_5/signaling_gene_vlnplot_calcr_notch.pdf", width = 8, height = 4)

## CALCR: plot ADM as signaling molecule and CALCR as the receptor
print(VlnPlot(obj, features = "ADM", group.by = "Annotation_1", pt.size=0, cols = col) + NoLegend())
print(VlnPlot(obj, features = "CALCR", group.by = "Annotation_1", pt.size=0, cols = col) + NoLegend())


## NOTCH: plot DLK1 as signaling molecule and NOTCH3 as the receptor
print(VlnPlot(obj, features = "DLK1", group.by = "Annotation_1", pt.size=0, cols = col) + NoLegend())
print(VlnPlot(obj, features = "NOTCH3", group.by = "Annotation_1", pt.size=0, cols = col) + NoLegend())

dev.off()







#######################################################################################
#Maybe supplementary figures: signaling sent from each cell type
mat <- obj.cellChat@net$weight
groupSize <- as.numeric(table(obj.cellChat@idents))


# the dot size correlate with the number of cells in that cell type
# all panels share the same edge-width scale, so arrow thickness is visually comparable across cell types.
pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/supp_figure/CellChat_outgoing_networks.pdf", width = 15, height = 10)

par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()

