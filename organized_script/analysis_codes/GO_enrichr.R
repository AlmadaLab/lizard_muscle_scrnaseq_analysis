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


print("GO enrichment analysis on cell type")


# Step 1: read in the object in human orthologs and define the customized order
obj <-readRDS("rds/annotated_obj_in_human_orthologs.rds")

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



# Step 2: use human ortholog, run FindAllMarkers

marker_genes <- function(data) {
  # set to seurat clusters
  Idents(data) <- "Annotation_1"
  markers <- FindAllMarkers(data, only.pos = TRUE)
  filtered_markers <- markers %>%
    dplyr::filter(p_val_adj < 0.05,
                  avg_log2FC > 0.5) %>%
    arrange(desc(avg_log2FC))
  return (filtered_markers)
}

cluster_markers <- marker_genes(obj)
write.csv(cluster_markers, "rds/anocarv3_human_orthologs_FindAllMarkers_log2fc0.5.csv")




# Step 3: search desired enrichr GO databases
dbs <- listEnrichrDbs()
# GO
dbs[grep("GO", dbs$libraryName), "libraryName"]

# select enrichr database
dbs <- c("GO_Biological_Process_2026")


# Step 3: define the function to perform GO analysis on each cell type (plotting)
go_analysis <- function(cluster_gene, name) {
  # enrichment analysis
  
  dbs <- c("GO_Biological_Process_2026")
  enrichment_results <- enrichr(cluster_gene, dbs)
  BP <- enrichment_results[[dbs]]
  
  # Calculate Gene Ratio
  BP <- BP %>%
    mutate(
      GeneCount = as.numeric(sub("/.*", "", Overlap)),
      GeneRatio = GeneCount / length(cluster_gene)
    )
  
  # Select top 10 terms by adjusted p-value
  BP_top <- BP %>%
    arrange(Adjusted.P.value) %>%
    slice_head(n = 10)
  
  # BP bar plot
  BP_plot <- ggplot(
    BP_top,
    aes(
      x = GeneRatio,
      y = reorder(Term, GeneRatio),
      fill = -log10(Adjusted.P.value)
    )
  ) +
    geom_col(width = 0.8) +
    labs(
      title = paste("Top 10 GO Biological Processes:", name),
      x = "Gene Ratio",
      y = NULL,
      fill = expression(-log[10]("Adjusted p-value"))
    ) +
    scale_fill_gradient(
      low = "#4C78A8",
      high = "#D73027"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      legend.title = element_text(size = 12)
    )
  
  return(BP_plot)
}

# define input genes for GO analysis: human ortholog DEG tables
# select more stringent markers based on log2fc: keep log2fc >=1 & p_adj < 0.05


DEG_markers <- read.csv("rds/anocarv3_human_orthologs_FindAllMarkers_log2fc0.5.csv")

filtered.markers <- DEG_markers %>%
  filter(avg_log2FC > 1)



# loop over all cell types to generate BP GO enrichment plot
celltypes <- unique(filtered.markers$cluster)

pdf(
  "/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_3/GO_BP_all_celltypes09152026.pdf",
  width = 12,
  height = 6,
  onefile = TRUE
)

for (ct in celltypes) {
  message("Processing: ", ct)
  cluster_gene <- filtered.markers %>%
    filter(cluster == ct) %>%
    pull(gene)
  p <- go_analysis(cluster_gene, ct)
  print(p)   # each print() creates a new page
}

dev.off()




