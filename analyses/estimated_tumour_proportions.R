# Estimated proportion of tumour-associated spots (Supplementary 5B).

# Define cell type marker gene classes
gene_classes <- list(
  Tumour = c("PRAC1", "HPN", "PCAT14", "AMACR", "PCA3"),
  Luminal = c("MSMB", "ACPP", "RDH11", "NKX3-1", "AZGP1"),
  Basal = c("KRT5", "KRT15", "TRIM29", "TP63", "SLC14A1"),
  Club = c("MMP7", "PIGR", "LTF", "CP", "KRT7"),
  Immune = c("CXCR4", "TRBC1", "CD3D", "LYZ", "CD79A"),
  Endothelium = c("EPAS1", "EMP1", "VWF", "PECAM1", "IFI27"),
  Fibroblast = c("DCN", "LUM", "FBLN1", "SFRP2"),
  Muscle = c("TAGLN", "ACTA2", "ACTG2", "MYH11", "MYL9")
)

prostate <- readRDS(file = </Path/To/tissue_section_Seurat_object>)
expr   <- GetAssayData(prostate, assay = "SCT", slot = "data")
means  <- Matrix::rowMeans(expr)
expr_n <- expr / means[rownames(expr)]
spot_clusts <- as.character(prostate$seurat_clusters)
cluster_ids <- sort(unique(spot_clusts))
cluster_scores <- sapply(gene_classes, function(genes) {
  genes_ok <- intersect(genes, rownames(expr_n))
  sapply(cluster_ids, function(cl) {
    cells <- WhichCells(prostate, idents = cl)
    mean(Matrix::rowMeans(expr_n[genes_ok, cells, drop = FALSE]))
  })
})
rownames(cluster_scores) <- cluster_ids
cluster_labels <- apply(cluster_scores, 1, function(x) {
  names(gene_classes)[which.max(x)]
})
idx <- match(spot_clusts, names(cluster_labels))
spot_labels <- cluster_labels[idx]
tumour_prop <- mean(spot_labels == "TUMOUR", na.rm = TRUE)
