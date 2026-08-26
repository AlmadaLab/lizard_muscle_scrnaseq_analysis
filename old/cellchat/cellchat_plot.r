sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/CellChat_v3_output.log", append = TRUE)
pdf("/project/aealmada_561/judyz/single_cell/anocar_v3/plot/CellChat_plots_v3.pdf", width = 10, height = 8)

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
cellChat <- readRDS("RDS/cellChat.rds")



#######################################################################################
# visualizations: Hierarchy plot
# change the receiver/sender numeric order to match the new customized order
# cellChat@netP$pathways, check all pathways
pathways <- c("CALCR", "NOTCH")


#######################################################################################
# extract all LR pairs and pathways
pathways.all <- cellChat@netP$pathways
write.csv(pathways.all, "/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/cellchat/signaling_pathways.csv")

all_LR_pairs <- extractEnrichedLR(cellChat, signaling = pathways.all, geneLR.return = FALSE)
write.csv(all_LR_pairs, "/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/cellchat/all_LR_pairs.csv")



#######################################################################################
# contribution of each ligand-receptor pair to the overall signaling pathway: selected three 
print(netAnalysis_contribution(cellChat, signaling = "CALCR"))
print(netAnalysis_contribution(cellChat, signaling = "NOTCH"))

#print(netAnalysis_contribution(cellChat, signaling =pathways.all))



#######################################################################################

# show all the interactions received by MuSCs, FAPs1
#print(netVisual_chord_gene(cellChat, sources.use = c(1,5,6,10), targets.use = c(2,4,14), signaling = "CALCR", legend.pos.x = 15))
#print(netVisual_chord_gene(cellChat, sources.use = c(4,14), targets.use = c(2, 4, 7, 14), signaling = "NOTCH", legend.pos.x = 15))
#print(netVisual_chord_gene(cellChat, sources.use = c(2,6), targets.use = c(1, 5, 6), signaling = "BMP", legend.pos.x = 15))

print(netVisual_chord_gene(cellChat, sources.use = c(1,2,3,10), targets.use = c(4,6,14), signaling = "CALCR", legend.pos.x = 15))
print(netVisual_chord_gene(cellChat, sources.use = c(6,14), targets.use = c(4, 6, 7, 14), signaling = "NOTCH", legend.pos.x = 15))

#######################################################################################
# identify signaling roles of cell groups
## compute network centrality scores
cellChat <- netAnalysis_computeCentrality(cellChat, slot.name = "netP")

## visualize scores using heatmap for each pathway and scatter plot for all pathways
print(netAnalysis_signalingRole_network(cellChat, signaling = pathways, width = 8, height = 2.5, font.size = 10))

#outgoing <- netAnalysis_signalingRole_heatmap(cellChat, signaling = c("CALCR", "BMP", "NOTCH"), pattern = "outgoing")
#incoming <- netAnalysis_signalingRole_heatmap(cellChat, signaling = c("CALCR", "BMP", "NOTCH"), pattern = "incoming")
#print(outgoing + incoming)


# poster, delete for manuscript
outgoing <- netAnalysis_signalingRole_heatmap(cellChat, signaling = c("CALCR", "NOTCH"), pattern = "outgoing")
incoming <- netAnalysis_signalingRole_heatmap(cellChat, signaling = c("CALCR", "NOTCH"), pattern = "incoming")
print(outgoing + incoming)

# close pdf
dev.off()

sink()