# Figure 2A
# Distribution of Moran’s I per biomarker gene across all samples

library(Seurat)
library(dplyr)
library(purrr)
library(ggplot2)

# Define biomarker gene list

biomarkers <- c(
  "KLK2","AZGP1","TPM2","GSN","BGN","FLNC","COL1A1","SFRP4","SRD5A2","FAM13C",
  "GSTM2","TPX2","TK1","TOP2A","RRM2","CENPM","PCLAF","CDK1","NUSAP1","BIRC5",
  "PTTG1","CDC20","CDCA3","PBK","CDKN3","ASF1B","FOXM1","KIF20A","CENPF","PLK1",
  "PRC1","ORC6","ASPM","RAD54L","RAD51","DLGAP5","CDCA8","BUB1B","CEP55","KIF11",
  "SKA1","DTL","MCM10","ANO7","MYBPC1","NFIB","TNFRSF19","CAMK2N1","PBX1","PCDH7",
  "THBS2","LASP1","RABGAP1","UBE2C","EPPK1","ZWILCH","S1PR4","IQGAP3","TSBP1",
  "THNSL2","KLF13","MICA","C17orf97","CCDC163","ACOT1"
)

# Locate Visium Seurat objects (both patients)

rds_files <- list.files(
  path       = "</Path/To/Seurat_Objects>/",
  pattern    = "^P[12]_.*visium.*\\.rds$",
  full.names = TRUE
)

# Extract Moran’s I values for biomarker genes

moran_df <- map_dfr(rds_files, function(f) {
  
  obj <- readRDS(f)
  mf  <- as.data.frame(obj[["SCT"]]@meta.features)
  mf$Gene <- rownames(mf)
  
  sample_short <- sub("^P[12]_(.*)_visium.*$", "\\1", basename(f))
  
  mf %>%
    filter(Gene %in% biomarkers) %>%
    transmute(
      Sample = sample_short,
      Gene,
      MoranI = MoransI_observed
    )
})

# Lock gene order for plotting

moran_df <- moran_df %>%
  mutate(Gene = factor(Gene, levels = biomarkers))

# Boxplot of Moran’s I per gene

p_box <- ggplot(moran_df, aes(x = Gene, y = MoranI)) +
  geom_boxplot(
    outlier.shape = 16,
    outlier.size  = 1,
    fill          = "steelblue",
    color         = "black"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      face  = "italic"
    ),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Distribution of Moran’s I per Biomarker Gene Across All Samples",
    x     = "Biomarker Gene",
    y     = "Moran's I"
  ) +
  coord_cartesian(ylim = c(-0.02, 0.8))

# Save figure

ggsave(
  filename = "</Path/To/Output>/MoranI_boxplot_all_samples.png",
  plot     = p_box,
  width    = 8,
  height   = 4,
  dpi      = 800
)
