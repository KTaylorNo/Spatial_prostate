# Supplementary figure 4
# Moran's I heatmap for biomarker genes across patient samples (here showing patient 1)

library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(viridis)

# Define biomarker gene list

biomarkers <- c(
  "KLK2","AZGP1","TPM2","GSN","BGN","FLNC","COL1A1","SFRP4","SRD5A2","FAM13C",
  "GSTM2","TPX2","TK1","TOP2A","RRM2","CENPM","PCLAF","CDK1","NUSAP1","BIRC5",
  "PTTG1","CDC20","CDCA3","PBK","CDKN3","ASF1B","FOXM1","KIF20A","CENPF","PLK1",
  "PRC1","ORC6","ASPM","RAD54L","RAD51","DLGAP5","CDCA8","BUB1B","CEP55","KIF11",
  "SKA1","DTL","MCM10","ANO7","MYBPC1","NFIB","TNFRSF19","CAMK2N1","PBX1","PCDH7",
  "THBS2","LASP1","RABGAP1","UBE2C","EPPK1","ZWILCH","S1PR4","IQGAP3","THNSL2",
  "KLF13","MICA","C17orf97","CCDC163","ACOT1"
)

# Locate Patient 1 Visium Seurat objects

rds_files <- list.files(
  path   = "</Path/To/Seurat_Objects>/",
  pattern = "^P1_.*visium.*\\.rds$",
  full.names = TRUE
)

# Extract Moran's I values for each gene × sample

moran_df <- map_dfr(rds_files, function(f) {
  
  obj <- readRDS(f)
  mf  <- as.data.frame(obj[["SCT"]]@meta.features)
  mf$Gene <- rownames(mf)
  
  sample_short <- sub("^P1_(.*)_visium.*$", "\\1", basename(f))
  
  mf %>%
    filter(Gene %in% biomarkers) %>%
    transmute(
      Sample = sample_short,
      Gene,
      MoranI = MoransI_observed
    )
})

# Ensure all Sample × Gene combinations exist

moran_df <- expand_grid(
  Sample = unique(moran_df$Sample),
  Gene   = biomarkers
) %>%
  left_join(moran_df, by = c("Sample", "Gene"))

# Lock ordering for heatmap axes

moran_df <- moran_df %>%
  mutate(
    Gene   = factor(Gene,   levels = biomarkers),
    Sample = factor(Sample, levels = unique(Sample))
  )

# Build heatmap (missing values shown in grey)

p_heat <- ggplot(moran_df, aes(x = Gene, y = Sample, fill = MoranI)) +
  geom_tile(color = "white", size = 0.2, na.rm = FALSE) +
  scale_fill_viridis_c(
    option   = "C",
    name     = "Moran's I",
    limits   = c(-0.02, 0.8),
    na.value = "grey80"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5,
      face = "italic"
    ),
    axis.text.y = element_text(size = 8),
    panel.grid  = element_blank()
  ) +
  labs(
    title = "Moran’s I for biomarker genes across Patient 1 samples",
    x     = "Biomarker Gene",
    y     = "Tissue Section"
  )

# Save heatmap

ggsave(
  filename = "</Path/To/Output>/MoranI_heatmap_P1.png",
  plot     = p_heat,
  width    = 8,
  height   = 2.5,
  dpi      = 300
)
