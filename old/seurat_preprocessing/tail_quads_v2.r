sink("/project/aealmada_561/judyz/single_cell/anocar_v2/r_script_output/tail_quads_v2_output.log", append = TRUE)
pdf("/project/aealmada_561/judyz/single_cell/anocar_v2/plot/individual_qc_plots_v2.pdf", width = 10, height = 8)

# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)


# set seed
set.seed(1234)

# read in both datasets (quads and tails)
quad_O.data <- Read10X(data.dir = "/project/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/Quads/outs/filtered_feature_bc_matrix") 

tail_O.data <- Read10X(data.dir = "/project/aealmada_561/judyz/single_cell/anocar_v2/cell_ranger/Tail/outs/filtered_feature_bc_matrix") 


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



# quality check
## calculate percent.mt for both datastes
quad_O[["percent.mt"]] <- PercentageFeatureSet(quad_O, pattern = "^MT-")
tail_O[["percent.mt"]] <- PercentageFeatureSet(tail_O, pattern = "^MT-")


## check original mitochondrial genes in tails and quads
print("for tail original data, percent.mt: ")
print(summary(tail_O@meta.data$percent.mt))

print("for quads original data, percent.mt: ")
print(summary(quad_O@meta.data$percent.mt))


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
# plot mitochondrial gene distribution and cutoffs together
combined <- merge(tail_O, y = quad_O, add.cell.ids = c("Tail", "Quads"), project = "Lizard")
combined_mt_density_plot <- ggplot(combined@meta.data, aes(x = percent.mt, color = orig.ident, fill = orig.ident)) +
                            geom_density(alpha = 0.3, size = 1) +  # Add transparency for fill
                            #geom_vline(xintercept = 15, color = "black", linetype = "dashed") +
                            #annotate("text", x = 20, y=0.2, label="15% threshold", angle=0, color="black") +
                            labs(
                              title = "Mitochondrial Gene Percentage Density Curve",
                              x = "Mitochondrial Gene Percentage (%)",
                              y = "Density",
                              fill = "Identity",    # Add custom label for the fill legend
                              color = "Identity" 
                            ) +
                            scale_x_continuous(n.breaks=5)+
                            theme_classic() +
                            scale_color_manual(values = c("lizard_tail" = "steelblue", "lizard_quads" = "coral")) +
                            scale_fill_manual(values = c("lizard_tail" = "steelblue", "lizard_quads" = "coral")) +
                            theme(legend.title = element_text(), legend.position = "right")
print(combined_mt_density_plot)


####################################################################################################################################################################
ggplot(combined@meta.data, aes(x = percent.mt, fill = orig.ident)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50, color = "black") +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed", size = 1) +
  labs(
    title = "Distribution of Mitochondrial Gene Percentage",
    x = "Mitochondrial Gene Percentage (%)",
    y = "Number of Cells",
    fill = "Sample"
  ) +
  scale_fill_manual(values = c("lizard_tail" = "steelblue", "lizard_quads" = "coral")) +
  theme_classic() +
  theme(legend.title = element_text(), legend.position = "right")

####################################################################################################################################################################
hist_mt_plot <- function(dataset, name) {
  ggplot(dataset@meta.data, aes(x = percent.mt)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +
  geom_vline(xintercept = 15, color = "black", linetype = "dashed", size = 1)+
  labs(title = paste("Mitochondrial Gene Percentage Distribution ", name),
       x = "Mitochondrial Gene Percentage (%)",
       y = "Cell Count") +
  theme_minimal()
}

print(hist_mt_plot(quad_O, "lizard quads"))
print(hist_mt_plot(tail_O, "lizard tail"))

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
  dataset <- subset(dataset, subset = nFeature_RNA > lower_feature & nFeature_RNA < upper_feature & nCount_RNA > lower_count & nCount_RNA < upper_count & percent.mt < 15)
  
  print(VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  print(FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + xlim(0, 40000) + ylim(0, 5000))
  
  
  # normalization
  dataset <- NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # highly variable feature
  dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 2000, slot = "data")
  
  # scale data
  dataset <- ScaleData(dataset, features = rownames(dataset))
  
  return(dataset)
}


quad_O_qc <- individual_dataset(quad_O)
tail_O_qc <- individual_dataset(tail_O)


## Statistics after 2nd step QC
print("for tail after 2nd step qc: ")
print(tail_O_qc)

print("for quads after 2nd step qc: ")
print(quad_O_qc)


# check qc mt gene in tails and quads
print("for tail qc data, percent.mt: ")
print(summary(tail_O_qc@meta.data$percent.mt))

print("for quads qc data, percent.mt: ")
print(summary(quad_O_qc@meta.data$percent.mt))


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
# visualization of the distribution of cells with high mt genes
visual_mt_thresholds <- function(raw, mt_threshold, sample_name) {

    dataset <- raw

    # label cells with mt% greater and under than threshold
    dataset$high_mt <- ifelse(dataset$percent.mt > mt_threshold, paste("Over ",mt_threshold), paste("Under ",mt_threshold))

    # count how many cells have mt% over and under threshold
    cell_counts <- table(dataset$high_mt)
    print(cell_counts)
    
    # normalization
    dataset <- NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  
    # highly variable feature
    dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 2000)
  
    # scale data
    all.genes <- rownames(dataset)
	dataset <- ScaleData(dataset, features = all.genes)
  
    # dimensional analysis
    dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset))
  
    # cluster analysis
	dataset <- FindNeighbors(dataset, dims = 1:30)
	dataset <- FindClusters(dataset, resolution = 0.2)

    # tsne for visualization
    dataset <- RunTSNE(dataset, dims = 1:30)
    cells_over_threshold <- WhichCells(dataset, expression = percent.mt > mt_threshold)
    print(
      DimPlot(
        dataset, 
        reduction = "tsne", 
        cells.highlight = cells_over_threshold,  
        cols.highlight = "red", 
        cols = "gray"
      ) +
      scale_color_manual(
        values = c("gray", "red"), 
        labels = c(paste("Under ", mt_threshold), paste("Over ", mt_threshold))
      ) +
      ggtitle(paste("Mitochondrial Threshold ", mt_threshold, "%", "for ", sample_name))
    )
}

# test 5%, 10%, 15%, 20% for quads
visual_mt_thresholds(quad_O, 5, "Quads")
visual_mt_thresholds(quad_O, 10, "Quads")
visual_mt_thresholds(quad_O, 15, "Quads")
visual_mt_thresholds(quad_O, 20, "Quads")


# test 5%, 10%, 15%, 20% for tail
visual_mt_thresholds(tail_O, 5, "Tail")
visual_mt_thresholds(tail_O, 10, "Tail")
visual_mt_thresholds(tail_O, 15, "Tail")
visual_mt_thresholds(tail_O, 20, "Tail")





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


# close pdf
dev.off()

sink()