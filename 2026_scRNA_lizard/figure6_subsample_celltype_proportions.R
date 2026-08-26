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


####################################################################################

# For figure 6: 
## randomly retain 80% of all cells within tail and quad separately
## repeat 1000 times
## calculate the cell-type proportions in each iteration
## plot the proportion



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


####################################################################################
# plot the cell type proportion
prop.original <- obj@meta.data %>%
  count(orig.ident, Annotation_1) %>%
  group_by(orig.ident) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

head(prop.original)


pdf("figures/figure_6/tail_quad_actual_composition.pdf", width = 8, height = 6)

p1 <- ggplot(
      prop.original,
      aes(x = Annotation_1, y = proportion, fill = orig.ident)
    ) +
      geom_col(position = "dodge") +
      scale_y_continuous(labels = scales::percent) +
      theme_classic() +
      scale_fill_manual(values = c(
        tail = "#4C78A8",
        quad = "#F58518"
      )) +
      labs(
        x = NULL,
        y = "Cell proportion"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
print(p1)
dev.off()

####################################################################################
# set threshold for iteration
set.seed(1)
n_iter <- 1000
sample_frac <- 0.8

# Get metadata
meta <- obj@meta.data
# Check
table(obj$orig.ident, obj$Annotation_1)


# subsampling
prop.results <- vector("list", n_iter) # store the final results
for (i in seq_len(n_iter)) {
  
  # Cells from each tissue
  tail.cells <- rownames(meta)[meta$orig.ident == "tail"]
  quad.cells <- rownames(meta)[meta$orig.ident == "quad"]
  
  # Randomly retain 80% without replacement
  tail.sample <- sample(
    tail.cells,
    size = floor(length(tail.cells) * sample_frac),
    replace = FALSE
  )
  
  quad.sample <- sample(
    quad.cells,
    size = floor(length(quad.cells) * sample_frac),
    replace = FALSE
  )
  
  # Combine sampled cells
  sampled.cells <- c(tail.sample, quad.sample)
  
  # Get metadata
  meta.sub <- meta[sampled.cells, , drop = FALSE]
  
  # Calculate cell-type proportions separately within Tail and Quad
  tmp <- meta.sub %>%
    count(orig.ident, Annotation_1, name = "n") %>%
    group_by(orig.ident) %>%
    mutate(
      total_cells = sum(n),
      proportion = n / total_cells
    ) %>%
    ungroup() %>%
    mutate(iteration = i)
  
  prop.results[[i]] <- tmp
}

# Combine all 1000 iterations
prop.results <- bind_rows(prop.results)



####################################################################################
# # plot the box plot of cell type proportion
# celltypes <- unique(prop.results$Annotation_1)
# 
# pdf("figures/celltype_proportion_subsampling.pdf", width = 4, height = 4)
# 
# for (ct in celltypes) {
#   
#   df <- prop.results %>%
#     filter(Annotation_1 == ct)
#   
#   p <- ggplot(df, aes(x = orig.ident, y = proportion, fill = orig.ident)) +
#     geom_boxplot(width = 0.5, outlier.shape = NA) +
#     theme_classic() +
#     scale_y_continuous(labels = scales::percent) +
#     labs(
#       title = ct,
#       x = NULL,
#       y = "Cell proportion"
#     ) +
#     theme(legend.position = "none")
#   
#   print(p)
# }
# 
# dev.off()



####################################################################################
# plot the delta distribution of tail and quad

pdf("figures/figure_6/tail_quad_delta_subsampling_iterate1000.pdf", width = 6, height = 6)

library(tidyr)
prop.diff <- prop.results %>%
  select(iteration, orig.ident, Annotation_1, proportion) %>%
  pivot_wider(
    names_from = orig.ident,
    values_from = proportion
  ) %>%
  mutate(diff = tail - quad)


p1 <- ggplot(prop.diff, aes(x = Annotation_1, y = diff)) +
        geom_boxplot(outlier.shape = NA) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_y_continuous(labels = scales::percent) +
        theme_classic() +
        labs(
          x = NULL,
          y = "Tail - Quad cell proportion"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
print(p1)
dev.off()

