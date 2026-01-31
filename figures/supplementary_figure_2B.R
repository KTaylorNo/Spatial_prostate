# Supplementary Figure 2B
# Spatial molecular clusters annotated by strongest marker gene class association

library(Seurat)
library(Matrix)
library(ggplot2)

# Load annotated Seurat object

prostate <- readRDS(
  file = "</Path/To/Annotated_Seurat_Object>.rds"
)

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

# Extract SCT-normalized expression and gene-normalize

expr   <- GetAssayData(prostate, assay = "SCT", layer = "data")
means  <- Matrix::rowMeans(expr)
expr_n <- expr / means[rownames(expr)]

# Retrieve cluster IDs for each spot

spot_clusts <- as.character(prostate$seurat_clusters)

# Compute average expression per marker class per cluster

cluster_ids <- sort(unique(spot_clusts))

cluster_scores <- sapply(gene_classes, function(genes) {
  genes_ok <- intersect(genes, rownames(expr_n))
  sapply(cluster_ids, function(cl) {
    cells <- WhichCells(prostate, idents = cl)
    mean(Matrix::rowMeans(expr_n[genes_ok, cells, drop = FALSE]))
  })
})

rownames(cluster_scores) <- cluster_ids

# Assign each cluster to its strongest marker class

cluster_labels <- apply(cluster_scores, 1, function(x) {
  names(gene_classes)[which.max(x)]
})

print(cluster_labels)

# Map cluster labels back to each spot

idx <- match(spot_clusts, names(cluster_labels))
spot_labels <- cluster_labels[idx]

# Add annotation to Seurat metadata

prostate$ClusterAnnotation <- spot_labels

# Quick check
head(prostate@meta.data[, c("seurat_clusters", "ClusterAnnotation")])

# Define colour palette

palette <- c(
  TUMOUR      = "#E41A1C",
  LUMINAL     = "#377EB8",
  BASAL       = "#F781BF",
  CLUB        = "#984EA3",
  IMMUNE      = "#FFD92F",
  ENDOTHELIUM = "#A65628",
  FIBROBLAST  = "#FF7F00",
  MUSCLE      = "#4DAF4A"
)

# Build plotting dataframe

coords  <- GetTissueCoordinates(prostate)
plot_df <- data.frame(coords, Annotation = prostate$ClusterAnnotation)
plot_df$Annotation <- factor(plot_df$Annotation, levels = names(palette))

# Plot spatial cluster annotation

p <- ggplot(plot_df, aes(x = x, y = y, colour = Annotation)) +
  geom_point(size = 2.5, alpha = 1) +
  coord_fixed() +
  theme_void() +
  scale_color_manual(values = palette) +
  labs(title = "Spatial Plot: Clusters Coloured by Marker Class")

# Save figure

ggsave(
  filename = "</Path/To/Output>/ClusterAnnotation_SpatialPlot.png",
  plot     = p,
  width    = 7,
  height   = 6,
  dpi      = 300
)
