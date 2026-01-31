# Supplementary Figure 8 dot plot
# Marker gene expression across pooled Patient 1 samples

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define marker genes (Kiviaho et al)

marker_genes <- c(
  "PRAC1", "HPN", "PCAT14", "AMACR", "PCA3",
  "MSMB", "ACPP", "RDH11", "NKX3-1", "AZGP1",
  "KRT5", "KRT15", "TRIM29", "TP63", "SLC14A1",
  "MMP7", "PIGR", "LTF", "CP", "KRT7",
  "CXCR4", "TRBC1", "CD3D", "LYZ", "CD79A",
  "EPAS1", "EMP1", "VWF", "PECAM1", "IFI27",
  "DCN", "LUM", "FBLN1", "SFRP2", "COL1A1",
  "TAGLN", "ACTA2", "ACTG2", "MYH11", "MYL9"
)

# Patient 1 samples with spot-level annotations

samples <- c(
  "P1_H1_2", "P1_H1_4", "P1_H1_5",
  "P1_H2_1", "P1_H2_2", "P1_H2_5",
  "P1_V1_2"
)

# Load and pool samples

all_spots <- list()

for (sample in samples) {
  
  file_path <- paste0(
    "</Path/To/Annotated_Seurat_Objects>/",
    sample, "_annotated.rds"
  )
  
  if (!file.exists(file_path)) next
  
  obj <- readRDS(file_path)
  DefaultAssay(obj) <- "SCT"
  
  # Keep marker genes present in this sample
  genes_present <- intersect(marker_genes, rownames(obj))
  
  expr <- GetAssayData(obj, layer = "data")[genes_present, , drop = FALSE]
  ann_raw <- obj$Final_Annotations
  
  # Collapse annotations
  ann_clean <- case_when(
    ann_raw %in% c("NA", "Exclude", "", NA) ~ NA_character_,
    ann_raw %in% c("Benign", "Benign*") ~ "Benign",
    ann_raw %in% c("GG4", "GG4 Cribriform") ~ "GG4",
    ann_raw %in% c(
      "Vessel", "Transition_State", "PIN", "Nerve",
      "Inflammation", "Fat", "Chronic inflammation"
    ) ~ "Other",
    TRUE ~ ann_raw
  )
  
  df <- as.data.frame(t(expr))
  df$Annotation <- ann_clean
  df$Sample <- sample
  
  df <- df %>% filter(!is.na(Annotation))
  all_spots[[sample]] <- df
}

# Combine all samples

combined_df <- bind_rows(all_spots)

# Keep marker genes present in ANY sample

marker_genes <- marker_genes[marker_genes %in% colnames(combined_df)]

# Calculate % expressing and average expression

dot_df <- combined_df %>%
  pivot_longer(
    cols      = all_of(marker_genes),
    names_to  = "Gene",
    values_to = "Expression"
  ) %>%
  group_by(Annotation, Gene) %>%
  summarise(
    pct_expr = mean(Expression > 0) * 100,
    avg_expr = mean(Expression),
    .groups  = "drop"
  )

# Fill missing Annotation × Gene combinations with zeros

all_annotations <- unique(combined_df$Annotation)

complete_df <- expand.grid(
  Annotation       = all_annotations,
  Gene             = marker_genes,
  stringsAsFactors = FALSE
)

dot_df <- complete_df %>%
  left_join(dot_df, by = c("Annotation", "Gene")) %>%
  mutate(
    pct_expr = replace_na(pct_expr, 0),
    avg_expr = replace_na(avg_expr, 0)
  )

# Prepare plotting subset

genes_part1 <- marker_genes

annotation_order <- c("GG4", "GG2", "GG1", "Benign", "Stroma", "Other")

dot_df_part1 <- dot_df %>%
  filter(Gene %in% genes_part1) %>%
  mutate(
    Gene       = factor(Gene, levels = rev(genes_part1)),
    Annotation = factor(Annotation, levels = annotation_order)
  )

# Italic labels for y‑axis

italic_labels <- setNames(
  lapply(genes_part1, function(g) bquote(italic(.(g)))),
  genes_part1
)

# Build dot plot

p1 <- ggplot(dot_df_part1, aes(x = Annotation, y = Gene)) +
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_size(range = c(1, 8)) +
  scale_color_gradientn(
    colours = c("#2c7bb6", "#7b3294", "#d7191c"),
    limits  = c(0, 6)
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = italic_labels) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x.top    = element_text(angle = 45, hjust = 0, vjust = 0.5),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.title         = element_blank(),
    plot.margin        = margin(t = 30, r = 10, b = 10, l = 10)
  ) +
  labs(
    size  = "% spots with detection",
    color = "Avg expression"
  )

# Save plot

ggsave(
  filename = "</Path/To/Output_Directory>/supplementary_figure_8.png",
  plot     = p1,
  width    = 5,
  height   = 6,
  dpi      = 600
)
