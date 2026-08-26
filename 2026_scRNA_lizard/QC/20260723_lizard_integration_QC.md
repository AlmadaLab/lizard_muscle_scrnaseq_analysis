20260723_lizard_integration_QC
================
Judy Zhu
2026-07-23

# Introduction

QC the integrated object (AnoCar V3)

# QC metrics

nFeature = 500, percent.mt = 5

``` r
knitr::opts_chunk$set(
  echo = FALSE,      # hide code (so md shows only text + plots)
  message = FALSE,   # suppress messages
  warning = FALSE    # suppress warnings
)
```

# read in integrated object and number of cells

    ## An object of class Seurat 
    ## 22670 features across 11688 samples within 1 assay 
    ## Active assay: RNA (22670 features, 2000 variable features)
    ##  3 layers present: data, counts, scale.data
    ##  4 dimensional reductions calculated: pca, harmony, umap, tsne

# check RNA QC metrics

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

by seurat clusters

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

cluster 15 is a low-qual cluster? if lack of a clear marker genes we can
just drop

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

# Check integration and potential batch effect

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## cells per library passed initial qc

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

# Library composition per seurat clusters (barplots)

no batch effect, every seurat cluster is well represented by both
samples

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

# Visualize cell type specific marker genes

## FAP 1 and FAP 2

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

## FAP 3

pdgfra-

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-13-5.png)<!-- -->

## MuSC and Smooth muscle cells

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->

## Macrophages, B and T cells

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-17-3.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-17-4.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-17-5.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-17-6.png)<!-- -->

## endothelial cells and NK cells, neural

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-19-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-19-3.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-19-4.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-19-5.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-19-6.png)<!-- -->

## Bone marrow stem cells, oligo

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-21-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-21-3.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-21-4.png)<!-- -->

## LECs and megakaryocytes

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-23-2.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-23-3.png)<!-- -->![](20260723_lizard_integration_QC_files/figure-gfm/unnamed-chunk-23-4.png)<!-- -->
