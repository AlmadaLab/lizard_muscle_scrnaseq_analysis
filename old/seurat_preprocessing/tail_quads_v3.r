sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/tail_quads_v3_output.log", append = TRUE)
pdf("/project/aealmada_561/judyz/single_cell/anocar_v3/plot/individual_qc_plots_v3.pdf", width = 10, height = 8)

# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(svglite)
library(cowplot)

# save plots

#svglite("/project/aealmada_561/judyz/single_cell/anocar_v3/plot/individual_qc_plots_v3.svg")


# set seed
set.seed(1234)






# read in both datasets (quads and tails)
quad_O.data <- Read10X(data.dir = "/project/aealmada_561/judyz/single_cell/anocar_v3/cell_ranger/Quads/outs/filtered_feature_bc_matrix") 

tail_O.data <- Read10X(data.dir = "/project/aealmada_561/judyz/single_cell/anocar_v3/cell_ranger/Tail/outs/filtered_feature_bc_matrix") 


####################################################################################################################################################################
# raw data
quad_raw <- CreateSeuratObject(counts = quad_O.data, project = "lizard_quads")
tail_raw <- CreateSeuratObject(counts = tail_O.data, project = "lizard_tail")


# create seurat objects
## remove low quality cells & features
quad_O <- CreateSeuratObject(counts = quad_O.data, project = "lizard_quads", min.cells = 3, min.features = 200)
tail_O <- CreateSeuratObject(counts = tail_O.data, project = "lizard_tail", min.cells = 3, min.features = 200)


## Statistics before 1st step QC
print("for tail before 1st step qc: ")
print(tail_raw)

print("for quads before 1st step qc: ")
print(quad_raw)


## Statistics after 1st step QC
print("for tail after 1st step qc: ")
print(tail_O)

print("for quads after 1st step qc: ")
print(quad_O)


## check original nFeatures (genes) in tails and quads
print("for tail original data, nFeature: ")
print(summary(tail_O@meta.data$nFeature_RNA))

print("for quads original data, nFeature: ")
print(summary(quad_O@meta.data$nFeature_RNA))


## check original nCount (UMI) in tails and quads
print("for tail original data, nCount: ")
print(summary(tail_O@meta.data$nCount_RNA))

print("for quads original data, nCount: ")
print(summary(quad_O@meta.data$nCount_RNA))


####################################################################################################################################################################
# violin plot showing the cutoffs for tail

# feature
med_feature_tail <- median(tail_O@meta.data$nFeature_RNA)
MAD_feature_tail <- mad(tail_O@meta.data$nFeature_RNA)

## set lower and upper bound of MAD: +- 5*MAD
lower_feature_tail <- med_feature_tail - 5 * MAD_feature_tail
upper_feature_tail <- med_feature_tail + 5 * MAD_feature_tail

# count
med_count_tail <- median(tail_O@meta.data$nCount_RNA)
MAD_count_tail <- mad(tail_O@meta.data$nCount_RNA)

## set lower and upper bound of MAD: +- 5*MAD
lower_count_tail <- med_count_tail - 5 * MAD_count_tail
upper_count_tail <- med_count_tail + 5 * MAD_count_tail

# Tail plots
p1 <- VlnPlot(tail_O, features = "nFeature_RNA") +
  geom_hline(yintercept = c(lower_feature_tail, upper_feature_tail), linetype = "dashed", color = "black") +
  ggtitle("Tail: nFeature_RNA with Cutoff") +
  annotate("text", x=1, y = upper_feature_tail + 10, label = "5 MAD threshold", color = "black", angle = 0) +
  theme(legend.position = "none") + theme(axis.title.x = element_blank())

p2 <- VlnPlot(tail_O, features = "nCount_RNA") +
  geom_hline(yintercept = c(lower_count_tail, upper_count_tail), linetype = "dashed", color = "black") +
  ggtitle("Tail: nCount_RNA with Cutoff") +
  annotate("text", x=1, y = upper_count_tail + 10, label = "5 MAD threshold", color = "black", angle = 0) +
  theme(legend.position = "none") + theme(axis.title.x = element_blank())


print(VlnPlot(tail_O, features = "nFeature_RNA"))
print(VlnPlot(tail_O, features = "nCount_RNA"))
################################################################################################################
# violin plot showing the cutoffs for quads

# feature
med_feature_quad <- median(quad_O@meta.data$nFeature_RNA)
MAD_feature_quad <- mad(quad_O@meta.data$nFeature_RNA)

## set lower and upper bound of MAD: +- 5*MAD
lower_feature_quad <- med_feature_quad - 5 * MAD_feature_quad
upper_feature_quad <- med_feature_quad + 5 * MAD_feature_quad

# count
med_count_quad <- median(quad_O@meta.data$nCount_RNA)
MAD_count_quad <- mad(quad_O@meta.data$nCount_RNA)

## set lower and upper bound of MAD: +- 5*MAD
lower_count_quad <- med_count_quad - 5 * MAD_count_quad
upper_count_quad <- med_count_quad + 5 * MAD_count_quad


# Quad plots
p3 <- VlnPlot(quad_O, features = "nFeature_RNA") +
  geom_hline(yintercept = c(lower_feature_quad, upper_feature_quad), linetype = "dashed", color = "black") +
  ggtitle("Quad: nFeature_RNA with Cutoff") +
  annotate("text", x = 1, y = upper_feature_quad + 10, label = "5 MAD threshold", color = "black", angle = 0) +
  theme(legend.position = "none") + theme(axis.title.x = element_blank())

p4 <- VlnPlot(quad_O, features = "nCount_RNA") +
  geom_hline(yintercept = c(lower_count_quad, upper_count_quad), linetype = "dashed", color = "black") +
  ggtitle("Quad: nCount_RNA with Cutoff") +
  annotate("text", x = 1, y = upper_count_quad + 10, label = "5 MAD threshold", color = "black", angle = 0) +
  theme(legend.position = "none") + theme(axis.title.x = element_blank())

print(VlnPlot(quad_O, features = "nFeature_RNA"))
print(VlnPlot(quad_O, features = "nCount_RNA"))

################################################################################################################
# combine plots
combine_feature <- plot_grid(p1, p3, ncol = 2)
combine_count <- plot_grid(p2, p4, ncol = 2)
print(combine_feature)
print(combine_count)



####################################################################################################################################################################
combined <- merge(tail_O, y = quad_O, add.cell.ids = c("Tail", "Quads"), project = "Lizard")

med_feature_tail <- median(tail_O@meta.data$nFeature_RNA)
MAD_feature_tail <- mad(tail_O@meta.data$nFeature_RNA)
upper_feature_tail <- med_feature_tail + 5 * MAD_feature_tail

med_feature_quad <- median(quad_O@meta.data$nFeature_RNA)
MAD_feature_quad <- mad(quad_O@meta.data$nFeature_RNA)
upper_feature_quad <- med_feature_quad + 5 * MAD_feature_quad

p_combined_feature <- VlnPlot(combined, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0.3) +
  geom_segment(aes(x = 0.6, xend = 1.4, y = upper_feature_quad, yend = upper_feature_quad), 
               linetype = "dashed", color = "red") +
  geom_segment(aes(x = 1.6, xend = 2.4, y = upper_feature_tail, yend = upper_feature_tail), 
               linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = upper_feature_quad + 100, label = "Quads: 5 MAD threshold", size = 3) +
  annotate("text", x = 2, y = upper_feature_tail + 100, label = "Tail: 5 MAD threshold", size = 3) +
  theme_classic() +
  labs(title = "nFeature_RNA with 5 MAD Thresholds", x = NULL, y = "nFeature_RNA")

print(p_combined_feature)





med_count_tail <- median(tail_O@meta.data$nCount_RNA)
MAD_count_tail <- mad(tail_O@meta.data$nCount_RNA)
upper_count_tail <- med_count_tail + 5 * MAD_count_tail

med_count_quad <- median(quad_O@meta.data$nCount_RNA)
MAD_count_quad <- mad(quad_O@meta.data$nCount_RNA)
upper_count_quad <- med_count_quad + 5 * MAD_count_quad

# Create violin plot for nCount_RNA with thresholds
p_combined_count <- VlnPlot(combined, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0.3) +
  geom_segment(aes(x = 0.6, xend = 1.4, y = upper_count_quad, yend = upper_count_quad), 
               linetype = "dashed", color = "red") +
  geom_segment(aes(x = 1.6, xend = 2.4, y = upper_count_tail, yend = upper_count_tail), 
               linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = upper_count_quad + 100, label = "Quads: 5 MAD threshold", size = 3) +
  annotate("text", x = 2, y = upper_count_tail + 100, label = "Tail: 5 MAD threshold", size = 3) +
  theme_classic() +
  labs(title = "nCount_RNA with 5 MAD Thresholds", x = NULL, y = "nCount_RNA")


print(p_combined_count)

####################################################################################################################################################################
## use MAD to process data
## funcion taking each dataset
individual_dataset <- function(dataset) {
  set.seed(1234)
  
  # feature
  med_feature <- median(dataset@meta.data$nFeature_RNA)
  MAD_feature <- mad(dataset@meta.data$nFeature_RNA)

  ## set lower and upper bound of MAD: +- 5*MAD
  lower_feature <- med_feature - 5*MAD_feature
  upper_feature <- med_feature + 5*MAD_feature
  
  # count
  med_count <- median(dataset@meta.data$nCount_RNA)
  MAD_count <- mad(dataset@meta.data$nCount_RNA)

  ## set lower and upper bound of MAD: +- 5*MAD
  lower_count <- med_count - 5*MAD_count
  upper_count <- med_count + 5*MAD_count


  # 5MAD, 15% as threshold
  # for v3, exclude mt threshold
  dataset <- subset(dataset, subset = nFeature_RNA > lower_feature & nFeature_RNA < upper_feature & nCount_RNA > lower_count & nCount_RNA < upper_count)
  
  #print(VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  print(FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + xlim(0, 40000) + ylim(0, 5000))
  
  
  # normalization
  dataset <- NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # highly variable feature
  dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 2000, slot = "data")
  
  # scale data
  #dataset <- ScaleData(dataset, features = rownames(dataset))
  
  return(dataset)
}


quad_O_qc <- individual_dataset(quad_O)
tail_O_qc <- individual_dataset(tail_O)


## Statistics after 2nd step QC
print("for tail after 2nd step qc: ")
print(tail_O_qc)

print("for quads after 2nd step qc: ")
print(quad_O_qc)


## check qc nFeatures (genes) in tails and quads
print("for tail qc data, nFeature: ")
print(summary(tail_O_qc@meta.data$nFeature_RNA))

print("for quads qc data, nFeature: ")
print(summary(quad_O_qc@meta.data$nFeature_RNA))


## check qc nCount (UMI) in tails and quads
print("for tail qc data, nCount: ")
print(summary(tail_O_qc@meta.data$nCount_RNA))

print("for quads qc data, nCount: ")
print(summary(quad_O_qc@meta.data$nCount_RNA))




####################################################################################################################################################################
# combine processed tail and quads to visualize nCount and nFeature scatter plot
combined <- merge(tail_O_qc, y = quad_O_qc, add.cell.ids = c("Tail", "Quads"), project = "Lizard")
print(FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="orig.ident") + xlim(0, 40000) + ylim(0, 5000))


# save data

if (!dir.exists("RDS")) {
  dir.create("RDS")
}

saveRDS(quad_O_qc, file = "RDS/quad_O_qc.rds")
saveRDS(tail_O_qc, file = "RDS/tail_O_qc.rds")
saveRDS(tail_O, file = "RDS/tail_O.rds")
saveRDS(quad_O, file = "RDS/quad_O.rds")
saveRDS(tail_raw, file = "RDS/tail_raw.rds")
saveRDS(quad_raw, file = "RDS/quad_raw.rds")


# close svg and pdf
dev.off()


sink()