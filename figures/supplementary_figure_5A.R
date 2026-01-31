# Supplementary Figure 5A
# Moran's I distribution of biomarker genes across all samples

library(Seurat)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

# Define biomarker gene panels

oncotype_dx <- c(
  "KLK2","AZGP1","TPM2","GSN","BGN","FLNC","COL1A1","SFRP4",
  "SRD5A2","FAM13C","GSTM2","TPX2"
)

prolaris <- c(
  "TK1","TOP2A","RRM2","CENPM","PCLAF","CDK1","NUSAP1","BIRC5",
  "PTTG1","CDC20","CDCA3","PBK","CDKN3","ASF1B","FOXM1","KIF20A",
  "CENPF","PLK1","PRC1","ORC6","ASPM","RAD54L","RAD51","DLGAP5",
  "CDCA8","BUB1B","CEP55","KIF11","SKA1","DTL","MCM10"
)

decipher <- c(
  "ANO7","MYBPC1","NFIB","TNFRSF19","CAMK2N1","PBX1","PCDH7",
  "THBS2","LASP1","RABGAP1","UBE2C","EPPK1","ZWILCH","S1PR4","IQGAP3"
)

proclass <- c("THNSL2","KLF13","MICA","C17orf97","CCDC163")

panel_colors <- c(
  "Oncotype Dx" = "purple",
  "Prolaris"    = "navyblue",
  "Decipher"    = "skyblue",
  "ProClass"    = "deeppink"
)

# Combine all genes and assign panel labels

genes_to_plot <- c(oncotype_dx, prolaris, decipher, proclass)

gene_panel_map <- tibble(
  Gene = genes_to_plot,
  Panel = factor(
    case_when(
      Gene %in% oncotype_dx ~ "Oncotype Dx",
      Gene %in% prolaris    ~ "Prolaris",
      Gene %in% decipher    ~ "Decipher",
      Gene %in% proclass    ~ "ProClass"
    ),
    levels = c("Oncotype Dx", "Prolaris", "Decipher", "ProClass")
  )
)

# Locate Visium Seurat objects

rds_files <- list.files(
  path       = "</Path/To/Seurat_Objects>/",
  pattern    = "visium.*\\.rds$",
  full.names = TRUE
)

# Extract Moran’s I values for biomarker genes

plot_df <- map_dfr(rds_files, function(f) {
  
  obj <- readRDS(f)
  sample_id <- sub("^(P[12]_.*)_visium.*\\.rds$", "\\1", basename(f))
  mf <- obj[["SCT"]]@meta.features
  
  tibble(
    Sample = sample_id,
    Gene   = genes_to_plot,
    MoranI = ifelse(
      genes_to_plot %in% rownames(mf),
      mf[genes_to_plot, "MoransI_observed"],
      NA_real_
    )
  )
}) %>%
  left_join(gene_panel_map, by = "Gene")


# Violin plot with panel‑coloured points

p <- ggplot(plot_df, aes(x = Sample, y = MoranI)) +
  geom_violin(
    fill  = "gold",
    alpha = 0.6,
    scale = "width",
    trim  = FALSE
  ) +
  geom_jitter(
    aes(color = Panel),
    width = 0.2,
    size  = 1.5,
    alpha = 0.9
  ) +
  scale_color_manual(values = panel_colors) +
  theme_minimal(base_size = 11) +
  labs(
    title    = "Moran’s I Distribution of Biomarker Genes per Sample",
    subtitle = "Points colored by gene panel",
    x        = "Sample",
    y        = "Moran's I",
    color    = "Panel"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save figure

ggsave(
  filename = "</Path/To/Output>/MoransI_violinplot_biomarkers.png",
  plot     = p,
  width    = 10,
  height   = 6,
  dpi      = 600
)
