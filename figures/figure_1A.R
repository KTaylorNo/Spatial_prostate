# Visualisation of filtered histopathological annotations
# Here showing for sample P1_H1_4

library(Seurat)
library(ggplot2)

# Load annotated Seurat object

seurat_obj <- readRDS(
  file = "</Path/To/Annotated_Seurat_Object/P1_H1_4_annotated.rds>"
)

# Identify spots with relevant histopathological annotations

valid_cells <- WhichCells(
  seurat_obj,
  expression = Final_Annotations %in% c(
    "Benign", "Benign*", "Stroma",
    "GG1", "GG2", "GG4", "GG4 Cribriform"
  )
)

filtered_obj <- seurat_obj[, valid_cells]

# Create simplified annotation labels

filtered_obj$Filtered_Annotations <- NA

filtered_obj$Filtered_Annotations[
  filtered_obj$Final_Annotations %in% c("Benign", "Benign*")
] <- "Benign"

filtered_obj$Filtered_Annotations[
  filtered_obj$Final_Annotations == "Stroma"
] <- "Stroma"

filtered_obj$Filtered_Annotations[
  filtered_obj$Final_Annotations == "GG1"
] <- "GG1"

filtered_obj$Filtered_Annotations[
  filtered_obj$Final_Annotations == "GG2"
] <- "GG2"

filtered_obj$Filtered_Annotations[
  filtered_obj$Final_Annotations %in% c("GG4", "GG4 Cribriform")
] <- "GG4"

# Define custom colour palette

custom_colors <- c(
  "Benign" = "#1f77b4",   # blue
  "Stroma" = "#2ca02c",   # green
  "GG1"    = "#ff7f0e",   # orange
  "GG2"    = "#d62728",   # red
  "GG4"    = "#800000"    # dark red
)

# Spatial plot with enhanced legend and larger spot size

plot_filtered <- SpatialDimPlot(
  filtered_obj,
  group.by       = "Filtered_Annotations",
  label          = FALSE,
  pt.size.factor = 3.2,
  image.alpha    = 0.6
) +
  scale_fill_manual(
    values = custom_colors,
    name   = "Annotation type"
  ) +
  scale_color_manual(
    values = custom_colors,
    name   = "Annotation type"
  ) +
  guides(
    fill  = guide_legend(override.aes = list(size = 6)),
    color = guide_legend(override.aes = list(size = 6))
  ) +
  theme(
    legend.key.size = unit(1.5, "lines"),
    legend.text     = element_text(size = 12),
    legend.title    = element_text(size = 13, face = "bold")
  )

# Save output

ggsave(
  filename = "</Path/To/Output_Directory/P1_H1_4_filtered_annotations.png>",
  plot     = plot_filtered,
  width    = 8,
  height   = 6,
  dpi      = 300
)
