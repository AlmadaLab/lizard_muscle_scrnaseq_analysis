# Singe-Cell RNA Sequencing of Muscle-Associated Cells in the Tail and Quadriceps in the Highly Regenerative Green Anole Lizard

Difei Zhu, Kartik Soni, Jisun Chin, Thomas P. Lozito, Albert E. Almada



## Overview

This folder contains the R scripts used to characterize muscle-associated cells in green anole (*Anolis carolinensis*) tail and quadriceps. Start with [](analysis_codes/) to process the datasets. Use [](plotting_codes/) to find the script for a specific manuscript figure.





## Analysis codes

### Analysis workflow

These R scripts analyze single-cell RNA sequencing data from green anole (*Anolis carolinensis*) tail and quadriceps skeletal muscle. The workflow integrates the two tissues, assigns cell type labels, and identifies markers for each annotated population. Two downstream analyses use human orthologs to interpret biological functions and infer cell–cell signaling.



| Script                                                       | Summary                                                      |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [calculate_percent.mt.R](analysis_codes/calculate_percent.mt.R) | Calculate mitochondrial RNA percentages from tail and quadriceps matrices aligned to AnoCar 2.0. |
| [integration_and_QC.R](analysis_codes/integration_and_QC.R)  | Transfer mitochondrial percentages to matching rAnoCar3.1.pri barcodes. Retain cells with more than 500 detected genes and less than 5% mitochondrial RNA, integrate tissues with Harmony, run clustering analysis, annotate cell types with marker genes. |
| [DEG_wilcox_findallmarkers.R](analysis_codes/DEG_wilcox_findallmarkers.R) | Identify markers for each annotated cell type using a Wilcoxon rank-sum test. |
| [cell_surface_markers.R](analysis_codes/cell_surface_markers.R) | Filter the cell type marker table to anole orthologs of human surface proteins from the Human Protein Atlas and Cell Surface Protein Atlas. |
| [ortholog_conversion.R](analysis_codes/ortholog_conversion.R) | Use a reciprocal best hit table to express anole counts and normalized expression in human gene symbols for the two analyses below. |
| [GO_enrichr.R](analysis_codes/GO_enrichr.R)                  | Identify markers in human orthologs and test each cell type’s marker set for enrichment of Gene Ontology biological processes. |
| [cellchat.R](analysis_codes/cellchat.R)                      | Infer ligand–receptor interactions and signaling pathways between cell types using the human CellChat database. |



### Manuscript figures

| Figure        | Script                                                    | Panels covered                                               |
| ------------- | --------------------------------------------------------- | ------------------------------------------------------------ |
| 2             | [plot_figure2.R](plotting_codes/plot_figure2.R)           | Mitochondrial RNA distribution, per-sample QC metrics, and t-SNE comparison before and after Harmony integration. |
| 3             | [plot_figure3.R](plotting_codes/plot_figure3.R)           | Annotated cell type t-SNE, literature-supported marker dot plot, and top marker heatmap. |
| 4             | [plot_figure4.R](plotting_codes/plot_figure4.R)           | Heatmap of candidate cell-surface markers and t-SNE expression plots for selected markers. |
| 5             | [plot_figure5.R](plotting_codes/plot_figure5.R)           | Inferred cell–cell signaling, including ligand–receptor contributions and CALCR and NOTCH pathway networks. |
| 6             | [plot_figure6.R](plotting_codes/plot_figure6.R)           | Tail–quadriceps correlations in average gene expression and comparisons of muscle stem cell (MuSC) marker expression. |
| Supplementary | [plot_supp_figures.R](plotting_codes/plot_supp_figures.R) | Reference genome annotation overlap, QC by cell type, outgoing CellChat networks, and cell type proportions by tissue. |





## Software Version

| Software    | Version |
| ----------- | ------- |
| R           | 4.4.2   |
| Cell Ranger | 8.0.1   |
| Seurat      | 5.2.1   |
| Harmony     | 1.2.3   |
| enrichR     | 3.4     |
| CellChat    | 2.1.2   |
