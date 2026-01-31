# To integrate and group the available spot-level histopathological annotations
# Used in Figure 1A, Figure 1B and Supplementary Figure 1.

library(tibble)
library(dplyr)

# Load Seurat object and annotations

prostate <- readRDS(
  file = "</Path/To/tissue_section_Seurat_object>"
)

annotations <- read.csv(
  "</Path/To/tissue_section_consensus_annotations_file>"
)

# Integrate annotations into Seurat metadata

prostate@meta.data <- prostate@meta.data %>%
  mutate(Barcode = rownames(.)) %>%
  left_join(annotations, by = "Barcode") %>%
  column_to_rownames("Barcode")
