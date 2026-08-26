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
print("Get cell surface markers from database and ortholog conversion")


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

# Step 2: upload those surface genes in human ortholog to ensembl and get the human to lizard ortholog conversion table
# https://useast.ensembl.org/info/data/biomart/index.html? 

# Step 3: filter the surface genes to only one to one ortholog and confidence = 1
ortholog <- read.delim(
  "database/surface_genes_human_to_lizard_biomart_export.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ortholog.filtered <- ortholog %>%
  filter(
    `Green anole homology type` == "ortholog_one2one",
    `Green anole orthology confidence [0 low, 1 high]` == 1
  )
length(unique(ortholog.filtered$`Green anole gene name`)) # number of unique genes with confident orthologs

# save those genes as the surface marker database (after ortholog converion with high confident of one to one)
length(surface.genes) # length of total surface genes: 2845
surface.genes.ortho <- unique(ortholog.filtered$`Green anole gene name`)

writeLines(
  surface.genes.ortho,
  "database/surface_genes_one2one_ortholog.txt"
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
markers <- read.csv("tables/annotated_obj_anocarv3_FindAllMarkers.csv")

# Read the surface gene list
surface.genes <- readLines("database/surface_genes_one2one_ortholog.txt")
surface.genes <- tolower(surface.genes)

length(surface.genes)
head(surface.genes)


# Step 3: filter to only surface markers
markers.surface <- markers %>%
  filter(gene %in% surface.genes)

write.csv(markers.surface, "tables/surface_markers_by_celltype_07302026.csv")




