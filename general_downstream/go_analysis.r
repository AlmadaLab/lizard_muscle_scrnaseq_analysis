sink("/project/aealmada_561/judyz/single_cell/anocar_v3/r_script_output/go_analysis_v3_output.log", append = TRUE)

pdf("/project/aealmada_561/judyz/single_cell/anocar_v3/plot/go_analysis_plot_v3.pdf", width = 16, height = 10)

# load packages (already installed on HPC)
.libPaths("/home1/difeizhu/R/x86_64-pc-linux-gnu-library/4.4")

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(harmony)
library(readr)
library(enrichR)
library(cowplot)

set.seed(1234) 

# load RDS
quad_O_qc <- readRDS("RDS/quad_O_qc.rds")
tail_O_qc <- readRDS("RDS/tail_O_qc.rds")
joint_combined <- readRDS("RDS/joint_combined_after_integration.rds") # seurat object 
combined.markers <- readRDS("RDS/combined.markers.rds") # DE gene list


annotate.ids <- c("FAPs 1", "MuSCs", "Macrophages", "Endothelial cells", "FAPs 2", "FAPs 3",
    "Smooth muscle cells", "T-cells", "Natural killer cells", "Neural/Glial cells",
    "B-cells", "Bone marrow stem cells", "Oligodendrocytes", "LECs", "Megakaryocytes")
names(annotate.ids) <- levels(joint_combined)
joint_combined <- RenameIdents(joint_combined, annotate.ids)

# set custom order to make three FAPs close to each other
Idents(joint_combined) <- factor(Idents(joint_combined),
  levels = c("FAPs 1", "FAPs 2", "FAPs 3", "MuSCs", "Macrophages", "Endothelial cells", 
             "Smooth muscle cells", "T-cells", "Natural killer cells", 
             "Neural/Glial cells", "B-cells", "Bone marrow stem cells", 
             "Oligodendrocytes", "LECs", "Megakaryocytes"))


# select enrichr database
dbs <- c("GO_Biological_Process_2023", "GO_Cellular_Component_2023", "GO_Molecular_Function_2023")



# function of GO analysis
go_analysis <- function(cluster_gene, name) {
    # enrichment analysis
    enrichment_results <- enrichr(cluster_gene, dbs)


    # BP
    BP <- enrichment_results[["GO_Biological_Process_2023"]]
    BP_top <- BP[1:10, ]
  
  	# Calculate Gene Ratio
	BP_top$GeneCount <- as.numeric(sapply(strsplit(BP_top$Overlap, "/"), `[`, 1))
	BP_top$GeneRatio <- BP_top$GeneCount / length(cluster_gene)
  
    # BP bar plot
    BP_plot <- ggplot(BP_top, aes(x = GeneRatio, y = reorder(Term, GeneRatio), fill = -log10(P.value))) + 
                geom_bar(stat = "identity") +
                labs(title = paste("Top 10 GO Biological Processes", name), x = "Gene Ratio", y = "GO Term", fill = "-log10(p-value)") +
                scale_fill_gradient(low = "blue", high = "red") +
                theme_classic() +
                    theme(axis.text.y = element_text(size = 14),
                    axis.text.x = element_text(size = 14),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14))
  
  
    print(BP_plot)




    # CC
    CC <- enrichment_results[["GO_Cellular_Component_2023"]]
    CC_top <- CC[1:10, ]

    # Calculate gene ratio
    CC_top$GeneCount <- as.numeric(sapply(strsplit(CC_top$Overlap, "/"), `[`, 1))
	CC_top$GeneRatio <- CC_top$GeneCount / length(cluster_gene)

    # CC bar plot
    CC_plot <- ggplot(CC_top, aes(x = GeneRatio, y = reorder(Term, GeneRatio), fill = -log10(P.value))) + 
                geom_bar(stat = "identity") +
                labs(title = paste("Top 10 GO Cellular Component", name), x = "Gene Ratio", y = "GO Term", fill = "-log10(p-value)") +
                scale_fill_gradient(low = "blue", high = "red") +
                theme_classic() +
                    theme(axis.text.y = element_text(size = 14),
                    axis.text.x = element_text(size = 14),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14))

    print(CC_plot)



    # MF
    MF <- enrichment_results[["GO_Molecular_Function_2023"]]
    MF_top <- MF[1:10, ]

    # Calculate gene ratio
    MF_top$GeneCount <- as.numeric(sapply(strsplit(MF_top$Overlap, "/"), `[`, 1))
	MF_top$GeneRatio <- MF_top$GeneCount / length(cluster_gene)

    # MF bar plot
    MF_plot <- ggplot(MF_top, aes(x = GeneRatio, y = reorder(Term, GeneRatio), fill = -log10(P.value))) + 
                geom_bar(stat = "identity") +
                labs(title = paste("Top 10 GO Molecular Function", name), x = "Gene Ratio", y = "GO Term", fill = "-log10(p-value)") +
                scale_fill_gradient(low = "blue", high = "red") +
                theme_classic() +
                    theme(axis.text.y = element_text(size = 14),
                    axis.text.x = element_text(size = 14),
                    axis.title.x = element_text(size = 14),
                    axis.title.y = element_text(size = 14))

    print(MF_plot)
}

# select more stringent markers based on log2fc: keep log2fc >=1
filtered.markers <- combined.markers %>%
  filter(avg_log2FC > 1)

# call function
go_analysis((filtered.markers[filtered.markers$cluster == 0, ])$gene, "FAPs 1")
go_analysis((filtered.markers[filtered.markers$cluster == 1, ])$gene, "MuSCs")
go_analysis((filtered.markers[filtered.markers$cluster == 2, ])$gene, "Macrophages")
go_analysis((filtered.markers[filtered.markers$cluster == 3, ])$gene, "Endothelial cells")
go_analysis((filtered.markers[filtered.markers$cluster == 4, ])$gene, "FAPs 2")
go_analysis((filtered.markers[filtered.markers$cluster == 5, ])$gene, "FAPs 3")
go_analysis((filtered.markers[filtered.markers$cluster == 6, ])$gene, "Smooth muscle cells")
go_analysis((filtered.markers[filtered.markers$cluster == 7, ])$gene, "T-cells")
go_analysis((filtered.markers[filtered.markers$cluster == 8, ])$gene, "Natural killer cells")
go_analysis((filtered.markers[filtered.markers$cluster == 9, ])$gene, "Neural/Glial cells")
go_analysis((filtered.markers[filtered.markers$cluster == 10, ])$gene, "B-cells")
go_analysis((filtered.markers[filtered.markers$cluster == 11, ])$gene, "Bone marrow stem cells")
go_analysis((filtered.markers[filtered.markers$cluster == 12, ])$gene, "Oligodendrocytes")
go_analysis((filtered.markers[filtered.markers$cluster == 13, ])$gene, "LECs")
go_analysis((filtered.markers[filtered.markers$cluster == 14, ])$gene, "Megakaryocytes")








# close pdf
dev.off()

sink()