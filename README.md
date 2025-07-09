# **Singe-Cell RNA Sequencing of Muscle-Associated Cells in the tail and quadriceps of the Super-Regenerating Green Anole Lizard** 

Difei Zhu, Kartik Soni, Jisun Chin, Thomas P. Lozito, Albert E. Almada



## Overview

This repository contains the code for processing and analyzing single-cell RNA sequencing (scRNA-seq) data from the tail and quadriceps muscles of the green anole lizard (*Anolis carolinensis*). Analyses include data preprocessing and quality check, integration, clustering, differential expression, functional annotation, cell-cell communication, and dataset comparison between tail and quadriceps.





## Repository Content

### cellchat/

Code for analyzing and visualizing cell-cell communication networks in the combined scRNA-seq dataset.

- Related to Figure 5

### cellranger_processing/

Code for processing raw sequencing data using AnoCar 2.0 and rAnoCar3.1.pri reference genomes. Also includes comparison of the number of annonated genes between the two reference genome versions.

- Related to Table 1, 2, and Figure 2a

### comparison/

Code for comparing tail and quadriceps datasets, including cell type proportions, bulk-level gene expression, and myogenic gene expression.

- Related to Figure 6

### general_downstream/

Code for downstream analyses such as functional annotation, extraction of cell surface markers, and gene ontology enrichment.

- Related to Figures 3, 4, and Supplementary Figure 3
- `cell_surface_genes.r`: Extracts surface-expressed genes for each cluster using a publicly available database

### integration/

Code for data integration (tail and quadriceps), clustering, and differential gene expression analysis across clusters

- Related to Supplementary Figure 1

### seurat_preprocessing/

Code for preprocessing scRNA-seq matrices using AnoCar 2.0 and rAnoCar3.1.pri genomes.

- AnoCar2.0: used for Figure 2b and Supplementary Figure 2 (mitochondrial gene percentages)
- rAnoCar3.2: used for figures 2c-e



## Software and Package Versions

**R**: 4.4.2

**CellRanger**: 8.0.1

**Seurat**: 5.2.1

**Harmony**: 1.2.3

**enrichR**: 3.4

**CellChat**: 2.1.2

