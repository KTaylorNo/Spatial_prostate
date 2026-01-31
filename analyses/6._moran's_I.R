#To calculate Global Moran’s I (a measure of spatial autocorrelation)
#Used in Figure 2 and Supplementary Figures 4 and 5A.

# Downstream analysis steps following SCTransform normalisation

prostate <- RunPCA(
  prostate,
  assay   = "SCT",
  verbose = FALSE
)

prostate <- FindNeighbors(
  prostate,
  reduction = "pca",
  dims      = 1:30
)

prostate <- FindClusters(
  prostate,
  verbose = FALSE
)

prostate <- RunUMAP(
  prostate,
  reduction = "pca",
  dims      = 1:30,
  verbose   = FALSE
)


# Calculate Moran’s I and adjust p-values (Benjamini–Hochberg)

prostate <- FindSpatiallyVariableFeatures(
  prostate,
  assay            = "SCT",
  slot             = "scale.data",
  selection.method = "moransi"
)

prostate[["SCT"]]@meta.features$Morans_I_p_value_adj <- p.adjust(
  prostate[["SCT"]]@meta.features$Morans_I_p_value,
  method = "BH"
)
