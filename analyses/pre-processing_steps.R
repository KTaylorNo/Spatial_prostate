#To create a Seurat object and conduct filtering and normalisation.

library(Seurat)

# Load Visium data for a single tissue section
prostate <- Load10X_Spatial(
  data.dir = "</Path/To/Pipeline_output_directory>"
)

# Filter out spots with fewer than 500 UMIs
prostate <- subset(
  prostate,
  subset = nCount_Spatial >= 500
)

# SCTransform normalisation
prostate <- SCTransform(
  prostate,
  assay = "Spatial",
  verbose = TRUE,
  return.only.var.genes = FALSE
)
