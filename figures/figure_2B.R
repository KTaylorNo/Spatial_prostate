# Figure 2B
# Average Moran’s I vs. spot coverage for biomarker panels
# Here showing samples from patient 1

library(ggplot2)
library(dplyr)
library(ggrepel)
library(Seurat)

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
  "THBS2","LASP1","RABGAP1","UBE2C","EPPK1","ZWILCH","S1PR4",
  "IQGAP3"
)

proclass <- c("THNSL2","KLF13","MICA","C17orf97","CCDC163","ACOT1")

# Build unified panel map

panel_map <- c(
  setNames(rep("Oncotype Dx", length(oncotype_dx)), oncotype_dx),
  setNames(rep("Prolaris",    length(prolaris)),    prolaris),
  setNames(rep("Decipher",    length(decipher)),    decipher),
  setNames(rep("ProClass",    length(proclass)),    proclass)
)

all_genes <- names(panel_map)

panel_colors <- c(
  "Oncotype Dx" = "purple",
  "Prolaris"    = "navyblue",
  "Decipher"    = "skyblue",
  "ProClass"    = "deeppink"
)

# Patient 1 samples

samples <- c(
  "P1_H1_1","P1_H1_2","P1_H1_4","P1_H1_5",
  "P1_H2_1","P1_H2_2","P1_H2_5",
  "P1_V1_1","P1_V1_2"
)

# Extract per‑sample values

all_data <- data.frame()

for (sample in samples) {
  
  file_path <- paste0("</Path/To/Seurat_Objects>/", sample, "_visium.rds")
  if (!file.exists(file_path)) next
  
  prostate <- readRDS(file_path)
  DefaultAssay(prostate) <- "SCT"
  
  available_genes <- intersect(all_genes, rownames(prostate[["SCT"]]))
  if (length(available_genes) == 0) next
  
  counts <- GetAssayData(prostate, assay = "SCT", layer = "counts")
  prop_detected <- rowSums(counts[available_genes, ] > 0) / ncol(counts)
  
  meta <- prostate[["SCT"]]@meta.features
  moransI <- meta[available_genes, "MoransI_observed"]
  
  sample_df <- data.frame(
    Gene        = available_genes,
    Proportion  = prop_detected,
    MoransI     = moransI,
    Panel       = panel_map[available_genes],
    Sample      = sample
  )
  
  all_data <- rbind(all_data, sample_df)
}

# Average across samples

avg_df <- all_data %>%
  group_by(Gene, Panel) %>%
  summarise(
    Avg_Proportion = mean(Proportion, na.rm = TRUE),
    Avg_MoransI    = mean(MoransI, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Label = ifelse(Avg_MoransI > 0.1, Gene, NA))

# Scatter plot (colored by panel)

p <- ggplot(avg_df, aes(x = Avg_MoransI, y = Avg_Proportion, color = Panel)) +
  geom_point(size = 2) +
  geom_text_repel(
    aes(label = Label),
    fontface = "italic",
    size = 3,
    na.rm = TRUE
  ) +
  scale_color_manual(values = panel_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 0.5)) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Average Moran's I",
    y = "Average spot coverage",
    color = "Panel"
  )

# Save figure

ggsave(
  filename = "</Path/To/Output>/Figure2B_patient1.png",
  plot     = p,
  width    = 6,
  height   = 4,
  dpi      = 300
)
