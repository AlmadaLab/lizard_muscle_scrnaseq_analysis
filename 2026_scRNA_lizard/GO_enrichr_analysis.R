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




# Step 2: search desired enrichr GO databases
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

# define input genes for GO analysis
# select more stringent markers based on log2fc: keep log2fc >=1 & p_adj < 0.05


DEG_markers <- read.csv("rds/annotated_obj_anocarv3_FindAllMarkers.csv")

filtered.markers <- DEG_markers %>%
  filter(avg_log2FC > 1)



# loop over all cell types to generate BP GO enrichment plot
celltypes <- unique(filtered.markers$cluster)

pdf(
  "/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_3/GO_BP_all_celltypes07262026.pdf",
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





# Step 4: define the function to perform GO analysis on each cell type (save the full GO terms)

save_go_table <- function(cluster_gene, name, outdir = "/project2/aealmada_561/judyz/20260621_scRNA_lizard/tables/GO_tables") {
  
  dbs <- "GO_Biological_Process_2026"
  enrichment_results <- enrichr(cluster_gene, dbs)
  
  # calculate gene ratio and sorted by Adjusted.P.value
  BP <- enrichment_results[[dbs]] %>%
    mutate(
      GeneCount = as.numeric(sub("/.*", "", Overlap)),
      GeneRatio = GeneCount / length(cluster_gene)
    ) %>%
    arrange(Adjusted.P.value)
  
  write.csv(
    BP,
    file = file.path(
      outdir,
      paste0(gsub("[ /]", "_", name), "_GO_BP.csv")
    ),
    row.names = FALSE
  )
}


celltypes <- unique(filtered.markers$cluster)

for (ct in celltypes) {
  message("Saving GO table: ", ct)
  cluster_gene <- filtered.markers %>%
    filter(cluster == ct) %>%
    pull(gene)
  save_go_table(cluster_gene, ct)
}

