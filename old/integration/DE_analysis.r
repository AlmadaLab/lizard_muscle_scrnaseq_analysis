sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/DE_analysis_v3_output.log", append = TRUE)



# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(presto)

# load RDS
quad_O_qc <- readRDS("RDS/quad_O_qc.rds")
tail_O_qc <- readRDS("RDS/tail_O_qc.rds")
combined <- readRDS("RDS/combined_after_integration.rds")

set.seed(1234) 

print(combined)

# combine layers for integration data
DefaultAssay(combined) <- "RNA"
joint_combined <- JoinLayers(combined)

# only reports upregulated genes
combined.markers <- FindAllMarkers(joint_combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# sort based on avg_log2FC for each cluster 
combined.markers <- combined.markers %>%
  arrange(cluster, desc(avg_log2FC), p_val) %>%
  filter(avg_log2FC > 0.58, p_val_adj < 0.05)

print(head(combined.markers))


####################################################################################################################################################################
# check each cluster DE genes, total of 14 clusters

num_cluster <- unique(combined.markers$cluster)

for (cluster in num_cluster) {
    cluster_markers <- combined.markers[combined.markers$cluster == cluster, ]
    file_name <- paste0("/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/gene_list/cluster_", cluster, "_markers.csv")
    write.csv(cluster_markers, file = file_name, row.names = TRUE)
}




if (!dir.exists("RDS")) {
  dir.create("RDS")
}

saveRDS(combined.markers, file = "RDS/combined.markers.rds")
saveRDS(joint_combined, file = "RDS/joint_combined_after_integration.rds")


sink()
