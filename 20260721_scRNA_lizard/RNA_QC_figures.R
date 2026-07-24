sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")

# Codes to generate plot in figure 2 (general QC plots)


#################################################################################
# Plot the mitochondrial distribution using the unfiltered data

obj <- readRDS("rds/merged_obj_unfiltered_anocarV3.rds")


# Extract metadata
df <- obj@meta.data

# Plot
ggplot(df, aes(x = percent.mt, fill = orig.ident)) +
  geom_histogram(
    position = "identity",
    bins = 50,
    color = "black",
    alpha = 0.45
  ) +
  geom_vline(
    xintercept = 5,
    color = "red",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c(
    tail = "#4C78A8",
    quad = "#F58518"
  )) +
  labs(
    x = "Mitochondrial RNA (%)",
    y = "Number of Cells",
    fill = "Sample"
  ) +
  coord_cartesian(xlim = c(0, 28)) +
  theme_classic(base_size = 18) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )


#################################################################################
# Plot nfeature and ncount RNA on the QCed object
obj <- readRDS("rds/integrated_obj_nfeat500_mt5_anocarV3.rds")

VlnPlot(obj, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0,
  cols = c(
    tail = "#4C78A8",
    quad = "#F58518"
  )
) + NoLegend()


VlnPlot(obj, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0,
        cols = c(
          tail = "#4C78A8",
          quad = "#F58518"
        )
) + NoLegend()


VlnPlot(obj, features = "percent.mt", group.by = "orig.ident", pt.size = 0,
        cols = c(
          tail = "#4C78A8",
          quad = "#F58518"
        )
) + NoLegend()
