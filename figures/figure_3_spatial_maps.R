# Figure 3 & Supplementary Figure 6
# Spatial expression maps of specific genes in a given section (here showing MSMB)

library(Seurat)
library(ggplot2)

# Load Seurat object

prostate_sample <- readRDS(
  file = "</Path/To/Visium_Seurat_Object>.rds"
)

# Use SCT assay for visualization
DefaultAssay(prostate_sample) <- "SCT"

# Spatial feature plot for MSMB

plot_msmb <- SpatialFeaturePlot(
  prostate_sample,
  features       = "MSMB",
  slot           = "data",
  pt.size.factor = 2.5,
  alpha          = c(1, 1)
) +
  guides(
    fill = guide_colorbar(
      title = expression(italic(MSMB))
    )
  )

# Save output

ggsave(
  filename = "</Path/To/Output_Directory>/MSMB_spatial_expression.png",
  plot     = plot_msmb,
  width    = 6,
  height   = 6,
  dpi      = 300
)

# Statistics for each visualised gene

library(dplyr)
seurat_obj <- readRDS(
  file = "</Path/To/Seurat_Object>/P1_H1_4_visium.rds"
)
DefaultAssay(seurat_obj) <- "SCT"

# Define genes of interest
genes <- c("MSMB", "TPM2", "KLF13", "UROD", "PSMD4", "NFIB")

meta_features <- seurat_obj[["SCT"]]@meta.features
expr_matrix   <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
summary_df <- data.frame(
  Gene            = genes,
  Spots_Detected  = rowSums(expr_matrix[genes, ] > 0),
  Avg_UMI         = round(rowMeans(expr_matrix[genes, ]), 2),
  Morans_I        = round(meta_features[genes, "MoransI_observed"], 3),
  Adj_P_Value     = signif(meta_features[genes, "MoransI_p.value.adj"], 3)
)
print(summary_df)
