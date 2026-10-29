library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(readr)
library(enrichR)
library(cowplot)

set.seed(1) 
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")

####################################################################################
print("Get cell surface markers from database and ortholog conversion (using reciprocal match)")


# Step 1: read in two databases with surface genes 
subcellular_db <- read_tsv(
  "/project2/aealmada_561/judyz/single_cell/datasets/subcellular_location.tsv"
)

CSPA_db <- read.csv(
  "/project2/aealmada_561/judyz/single_cell/datasets/Cell_surface_protein_atlas.csv"
)


## get union set of cell surface gene sets from two databases
surface1 <- subcellular_db %>%
  filter(grepl("Plasma membrane|Cell Junctions", `Main location`)) %>%
  pull(`Gene name`) %>%
  unique()

surface2 <- CSPA_db %>%
  pull(ENTREZ.gene.symbol) %>%
  unique()

surface.genes <- union(surface1, surface2)

writeLines(
  surface.genes,
  "database/surface_genes.txt"
)

# Step 2: read in lizard-human ortholog table
ortholog <- read.csv("database/human_anole_union_RBH_filtered.csv")

# Find human surface genes with an anole ortholog
surface.gene.ortholog <- ortholog[
  ortholog$human %in% surface.genes,
  c("human", "anole")
]

# Get anole genes corresponding to human surface genes
surface.anole.genes <- unique(
  ortholog$anole[ortholog$human %in% surface.genes]
)


writeLines(
  surface.anole.genes,
  "database/surface_genes_anole_RBH_filtered_08052026.txt"
)



####################################################################################

# Step 1: read in the object and define the customized order
obj <-readRDS("rds/annotated_obj_anocarv3_nfeat500_mt5.rds")

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


# Step 2: read in FindAllMarkers and surface genes
markers <- read.csv("tables/annotated_obj_anocarv3_FindAllMarkers_log2fc1.csv")

# Read the surface gene list
surface.genes <- readLines("database/surface_genes_anole_RBH_filtered_08052026.txt")
surface.genes <- tolower(surface.genes)

length(surface.genes)
head(surface.genes)


# Step 3: filter to only surface markers
markers.surface <- markers %>%
  filter(gene %in% surface.genes)

write.csv(markers.surface, "tables/surface_markers_by_celltype_RBH_log2fc1_10092026.csv")

