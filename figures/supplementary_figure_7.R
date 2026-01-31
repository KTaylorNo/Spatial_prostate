# Supplementary figure 7
# Spatial detectability heatmaps for biomarker panels per annotation type

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define samples and genes

samples <- c(
  "P1_H1_2", "P1_H1_4", "P1_H1_5",
  "P1_H2_1", "P1_H2_2", "P1_H2_5",
  "P1_V1_2"
)

genes <- c(
  "KLK2", "AZGP1", "TPM2", "GSN", "BGN", "FLNC", "COL1A1", "SFRP4", "SRD5A2", "FAM13C",
  "GSTM2", "TPX2", "TK1", "TOP2A", "RRM2", "CENPM", "PCLAF", "CDK1", "NUSAP1", "BIRC5",
  "PTTG1", "CDC20", "CDCA3", "PBK", "CDKN3", "ASF1B", "FOXM1", "KIF20A", "CENPF", "PLK1",
  "PRC1", "ORC6", "ASPM", "RAD54L", "RAD51", "DLGAP5", "CDCA8", "BUB1B", "CEP55", "KIF11",
  "SKA1", "DTL", "MCM10", "ANO7", "MYBPC1", "NFIB", "TNFRSF19", "CAMK2N1", "PBX1", "PCDH7",
  "THBS2", "LASP1", "RABGAP1", "UBE2C", "EPPK1", "ZWILCH", "S1PR4", "IQGAP3", "THNSL2",
  "KLF13", "MICA", "C17orf97", "CCDC163", "ACOT1"
)

# Define biomarker panels

oncotype_dx <- c("KLK2","AZGP1","TPM2","GSN","BGN","FLNC","COL1A1","SFRP4",
                 "SRD5A2","FAM13C","GSTM2","TPX2")

prolaris <- c("TK1","TOP2A","RRM2","CENPM","PCLAF","CDK1","NUSAP1","BIRC5",
              "PTTG1","CDC20","CDCA3","PBK","CDKN3","ASF1B","FOXM1","KIF20A",
              "CENPF","PLK1","PRC1","ORC6","ASPM","RAD54L","RAD51","DLGAP5",
              "CDCA8","BUB1B","CEP55","KIF11","SKA1","DTL","MCM10")

decipher <- c("ANO7","MYBPC1","NFIB","TNFRSF19","CAMK2N1","PBX1","PCDH7",
              "THBS2","LASP1","RABGAP1","UBE2C","EPPK1","ZWILCH","S1PR4",
              "IQGAP3")

proclass <- c("THNSL2","KLF13","MICA","C17orf97","CCDC163","ACOT1")

panel_list <- list(
  "Oncotype Dx" = oncotype_dx,
  "Prolaris"    = prolaris,
  "Decipher"    = decipher,
  "ProClass"    = proclass
)

# Annotation order (top → bottom)
annotation_levels <- c("Stroma", "Other", "Benign", "GG1", "GG2", "GG4")

# Pool all spots across samples

all_spots_list <- list()

for (sample in samples) {
  
  path <- paste0(
    "</Path/To/Annotated_Seurat_Objects>/",
    sample, "_annotated.rds"
  )
  
  seurat_obj <- readRDS(path)
  
  counts <- GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  ann_raw <- seurat_obj$Final_Annotations
  
  ann_grouped <- case_when(
    ann_raw == "Stroma" ~ "Stroma",
    ann_raw %in% c("Benign", "Benign*") ~ "Benign",
    ann_raw %in% c("GG4", "GG4 Cribriform") ~ "GG4",
    ann_raw == "GG2" ~ "GG2",
    ann_raw == "GG1" ~ "GG1",
    TRUE ~ "Other"
  )
  
  genes_present <- genes[genes %in% rownames(counts)]
  counts_sub <- counts[genes_present, , drop = FALSE]
  
  df <- data.frame(
    Spot = colnames(counts_sub),
    Annotation = ann_grouped,
    Sample = sample,
    t(as.matrix(counts_sub > 0))
  )
  
  all_spots_list[[sample]] <- df
}

all_spots_df <- bind_rows(all_spots_list)

# Compute proportions per gene × annotation

detect_long <- all_spots_df %>%
  select(Annotation, all_of(genes[genes %in% colnames(all_spots_df)])) %>%
  pivot_longer(
    cols = -Annotation,
    names_to = "Gene",
    values_to = "Detected"
  )

detect_long$Annotation <- factor(detect_long$Annotation, levels = annotation_levels)

prop_df <- detect_long %>%
  group_by(Annotation, Gene) %>%
  summarise(Proportion = mean(Detected), .groups = "drop") %>%
  complete(
    Annotation = factor(annotation_levels, levels = annotation_levels),
    Gene = genes,
    fill = list(Proportion = 0)
  )

# Panel‑specific heatmap function

plot_panel_heatmap <- function(panel_name, gene_list) {
  
  df <- prop_df %>%
    filter(Gene %in% gene_list) %>%
    mutate(Gene = factor(Gene, levels = gene_list))
  
  ggplot(df, aes(x = Gene, y = Annotation, fill = Proportion)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Proportion\nof spots", limits = c(0, 1)) +
    scale_y_discrete(drop = FALSE, limits = rev(annotation_levels)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6, face = "italic"),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title   = element_text(size = 14, face = "bold")
    ) +
    labs(title = paste("Spatial detectability of", panel_name, "genes"))
}


# Save each panel heatmap

ggsave(
  filename = "</Path/To/Output>/OncotypeDx_heatmap.png",
  plot     = plot_panel_heatmap("Oncotype Dx", oncotype_dx),
  width    = 4, height = 2, dpi = 300
)

ggsave(
  filename = "</Path/To/Output>/Prolaris_heatmap.png",
  plot     = plot_panel_heatmap("Prolaris", prolaris),
  width    = 6, height = 2, dpi = 300
)

ggsave(
  filename = "</Path/To/Output>/Decipher_heatmap.png",
  plot     = plot_panel_heatmap("Decipher", decipher),
  width    = 4, height = 2, dpi = 300
)

ggsave(
  filename = "</Path/To/Output>/ProClass_heatmap.png",
  plot     = plot_panel_heatmap("ProClass", proclass),
  width    = 3, height = 2, dpi = 300
)
