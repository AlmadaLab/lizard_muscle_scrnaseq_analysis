sink(stderr(), type = "output")

# load packages
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(1)
setwd("/project2/aealmada_561/judyz/20260621_scRNA_lizard")

# Codes to generate QC plot
## Plots for figure 2 and supplementary figure 1


#################################################################################
# Plot the mitochondrial distribution using the unfiltered data

obj <- readRDS("rds/merged_obj_unfiltered_anocarV3.rds")


# Extract metadata
df <- obj@meta.data

# Plot
p1 <- ggplot(df, aes(x = percent.mt, fill = orig.ident)) +
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

pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_2/percent_mt_histogram.pdf", width = 7, height = 6)
print(p1)
dev.off()



#################################################################################
# Plot nfeature and ncount RNA on the QCed object
obj <- readRDS("rds/annotated_obj_anocarv3_nfeat500_mt5.rds")

p2 <- VlnPlot(obj, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0,
              cols = c(
                tail = "#4C78A8",
                quad = "#F58518"
              )
) + NoLegend()


p3 <- VlnPlot(obj, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0,
              cols = c(
                tail = "#4C78A8",
                quad = "#F58518"
              )
) + NoLegend()


p4 <- VlnPlot(obj, features = "percent.mt", group.by = "orig.ident", pt.size = 0,
              cols = c(
                tail = "#4C78A8",
                quad = "#F58518"
              )
) + NoLegend()


p5 <- FeatureScatter(obj,
                     cols=c(tail = "#4C78A8", quad = "#F58518"), 
                     feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                     group.by="orig.ident")


pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_2/RNA_QC_by_samples_vlnplot.pdf", width = 7, height = 6)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()


# obtain the median nFeature, nCount
summary(obj$nCount_RNA)
summary(obj$nFeature_RNA)

obj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(
    median_nCount_RNA = median(nCount_RNA),
    median_nFeature_RNA = median(nFeature_RNA)
  )

#################################################################################
# Plot the tsne before running harmony integration

obj <- readRDS("rds/Notintegrated_obj_nfeat500_mt5_anocarV3.rds")
p6 <- DimPlot(obj, reduction="tsne", group.by="orig.ident", raster=F, order=T, cols = c(tail = "#4C78A8", quad = "#F58518"))


obj.2 <- readRDS("rds/annotated_obj_anocarv3_nfeat500_mt5.rds")
p7 <- DimPlot(obj.2, reduction="tsne", group.by="orig.ident", raster=F, order=T, cols = c(tail = "#4C78A8", quad = "#F58518"))


pdf("/project2/aealmada_561/judyz/20260621_scRNA_lizard/figures/figure_2/batch_effect_tsne.pdf", width = 7, height = 6)
print(p6)
print(p7)
dev.off()

