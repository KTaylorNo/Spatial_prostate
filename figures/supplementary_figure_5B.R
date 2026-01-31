# Supplementary Figure 5B
# Estimated tumour proportion across all samples using cluster association with marker genes

library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(purrr)

# Define marker gene classes

gene_classes <- list(
  TUMOUR      = c("PRAC1","HPN","PCAT14","AMACR","PCA3"),
  LUMINAL     = c("MSMB","ACPP","RDH11","NKX3-1","AZGP1"),
  BASAL       = c("KRT5","KRT15","TRIM29","TP63","SLC14A1"),
  CLUB        = c("MMP7","PIGR","LTF","CP","KRT7"),
  IMMUNE      = c("CXCR4","TRBC1","CD3D","LYZ","CD79A"),
  ENDOTHELIUM = c("EPAS1","EMP1","VWF","PECAM1","IFI27"),
  FIBROBLAST  = c("DCN","LUM","FBLN1","SFRP2"),
  MUSCLE      = c("TAGLN","ACTA2","ACTG2","MYH11","MYL9")
)

# Locate Visium Seurat objects

rds_files <- list.files(
  path       = "</Path/To/Seurat_Objects>/",
  pattern    = "visium\\.rds$",
  full.names = TRUE
)

# Compute tumour proportion per sample

tumour_props <- map_dfr(rds_files, function(f) {
  
  prostate <- readRDS(f)
  sample_id <- sub("^(P[12]_.*)_visium\\.rds$", "\\1", basename(f))
  
  # Normalized expression (SCT)
  expr   <- GetAssayData(prostate, assay = "SCT", layer = "data")
  means  <- Matrix::rowMeans(expr)
  expr_n <- expr / means[rownames(expr)]
  
  # Cluster IDs
  spot_clusts <- as.character(prostate$seurat_clusters)
  cluster_ids <- sort(unique(spot_clusts))
  
  # Score clusters by marker class
  cluster_scores <- sapply(gene_classes, function(genes) {
    genes_ok <- intersect(genes, rownames(expr_n))
    sapply(cluster_ids, function(cl) {
      cells <- WhichCells(prostate, idents = cl)
      mean(Matrix::rowMeans(expr_n[genes_ok, cells, drop = FALSE]))
    })
  })
  
  rownames(cluster_scores) <- cluster_ids
  
  # Assign strongest marker class per cluster
  cluster_labels <- apply(cluster_scores, 1, function(x) {
    names(gene_classes)[which.max(x)]
  })
  
  # Map cluster labels to spots
  idx <- match(spot_clusts, names(cluster_labels))
  spot_labels <- cluster_labels[idx]
  
  # Tumour proportion
  tumour_prop <- mean(spot_labels == "TUMOUR", na.rm = TRUE)
  
  tibble(Sample = sample_id, TumourProp = tumour_prop)
})

# Plot tumour proportions

p <- ggplot(tumour_props, aes(x = Sample, y = TumourProp)) +
  geom_col(fill = "red") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Proportion of spots annotated as tumour",
    x = "Sample",
    y = "Tumour Spot Percentage"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save figure

ggsave(
  filename = "</Path/To/Output>/TumourProportion_all_samples.png",
  plot     = p,
  width    = 10,
  height   = 2,
  dpi      = 600
)
