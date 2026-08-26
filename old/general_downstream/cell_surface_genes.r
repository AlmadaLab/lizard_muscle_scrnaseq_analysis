sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/cell_surface_markers_v3_output.log", append = TRUE)


# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(readr)

set.seed(1234) 

# load RDS
quad_O_qc <- readRDS("RDS/quad_O_qc.rds")
tail_O_qc <- readRDS("RDS/tail_O_qc.rds")
joint_combined <- readRDS("RDS/joint_combined_after_integration.rds") # seurat object 
combined.markers <- readRDS("RDS/combined.markers.rds") # DE gene list


num_cluster <- unique(combined.markers$cluster)


###################################
####FIND CELL SURFACE MARKERS######
###################################

#######################################################################################

# use only cell surface marker datasets: human protein atlas subcellular location
subcellular_db <- read_tsv("/project/aealmada_561/judyz/single_cell/datasets/subcellular_location.tsv")
head(subcellular_db) #subcellular_db$`Gene name`

# since the file contains all subcellular data, need process file making it only includes plasma membrane
plasma_membrane <- grepl("Plasma membrane|Cell Junctions", subcellular_db$`Main location`)

# View the result (give index): TRUE indicates rows contain plasma membrane
#print(plasma_membrane)

subset_df <- subcellular_db[plasma_membrane, ]
print(subset_df)


for (c in num_cluster) {

  overlap_surface_1 <- intersect(
    tolower(subset_df$`Gene name`),
    tolower((combined.markers[combined.markers$cluster == c, ])$gene)
  )

  overlap_surface_df_1 <- data.frame(gene = overlap_surface_1)

  fp_1 <- paste0("/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/cell_surface_markers/subcellular_db/intersect_", c, ".csv")
  
  # Write the data frame to CSV
  write.csv(overlap_surface_df_1, fp_1, row.names = FALSE)
}


#######################################################################################
# use Cell_surface_protein_atlas (CSPA)
CSPA_db <- read.csv("/project/aealmada_561/judyz/single_cell/datasets/Cell_surface_protein_atlas.csv")
head(CSPA_db) #CSPA_db$ENTREZ.gene.symbol


for (c in num_cluster) {

  # Convert both vectors to lowercase for a case-insensitive intersection
  overlap_surface_2 <- intersect(
    tolower(CSPA_db$ENTREZ.gene.symbol),
    tolower((combined.markers[combined.markers$cluster == c, ])$gene)
  )

  overlap_surface_df_2 <- data.frame(gene = overlap_surface_2)


  fp_2 <- paste0("/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/cell_surface_markers/CSPA_db/intersect_", c, ".csv")

  # Write the data frame to CSV
  write.csv(overlap_surface_df_2, fp_2, row.names = FALSE)
}



#######################################################################################
# combine two lists and remove duplicates, getting the new intersected list

list_surface_markers <- list()

for (c in num_cluster) {
  f1 <- paste0("/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/cell_surface_markers/subcellular_db/intersect_", c, ".csv")
  list_f1 <- read.csv(f1, header = TRUE)

  f2 <- paste0("/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/cell_surface_markers/CSPA_db/intersect_", c, ".csv")
  list_f2 <- read.csv(f2, header = TRUE)

  f12 <- unique(tolower(c(list_f1$gene, list_f2$gene)))

  ### Sort genes based on log2fold change ###
  # Read the DE gene file for the current cluster
  de_gene_file <- paste0("/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/gene_list/cluster_", c, "_markers.csv")
  de_gene <- read.csv(de_gene_file, header = TRUE)

  # Create a data frame with gene names and their log2fold change
  de_gene_df <- data.frame(gene = de_gene$gene, fold = de_gene$avg_log2FC)

  # Match genes from f12 to de_gene_df
  matched_data <- de_gene_df[de_gene_df$gene %in% f12, ]

  # Sort the matched genes by avg_log2FC in ascending order
  sorted_genes <- matched_data[order(-matched_data$fold), ]

  # View the sorted gene names
  print(sorted_genes$gene)


  # Add sorted gene list to list_surface_markers
  list_surface_markers[[paste0("Cluster_", c)]] <- sorted_genes$gene
}

# Find the maximum length among all the gene lists
max_length_1 <- max(sapply(list_surface_markers, length))

# Pad shorter lists with NA to make all lists the same length
list_surface_markers <- lapply(list_surface_markers, function(y) {
  length(y) <- max_length_1
  return(y)
})

# Convert the list to a matrix first, then to a data frame
final_surface_table <- do.call(cbind, list_surface_markers)

# Convert matrix to data frame
final_surface_table <- as.data.frame(final_surface_table)

# Write the data frame to CSV
write.csv(final_surface_table, "/project/aealmada_561/judyz/single_cell/anocar_v3/result_data/cell_surface_markers/surface_markers.csv", row.names = FALSE)



sink()


