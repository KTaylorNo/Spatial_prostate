# FIGURE 1C
# Average fraction of spots expressing biomarker panel genes

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


# Compute per-gene proportions for each sample

results_list <- list()

for (sample in samples) {
  
  path <- paste0(
    "</Path/To/Annotated_Seurat_Objects>/",
    sample, "_annotated.rds"
  )
  
  seurat_obj <- readRDS(path)
  counts <- GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  
  genes_present <- genes[genes %in% rownames(counts)]
  
  proportions <- sapply(genes_present, function(gene) {
    sum(counts[gene, ] > 0) / ncol(counts)
  })
  
  df <- data.frame(
    Gene       = genes_present,
    Proportion = proportions,
    Sample     = sample
  )
  
  results_list[[sample]] <- df
}

combined_df <- bind_rows(results_list)

# Compute average per gene across samples

avg_df <- combined_df %>%
  group_by(Gene) %>%
  summarise(AvgProportion = mean(Proportion)) %>%
  ungroup() %>%
  mutate(Panel = case_when(
    Gene %in% oncotype_dx ~ "Oncotype Dx",
    Gene %in% prolaris    ~ "Prolaris",
    Gene %in% decipher    ~ "Decipher",
    Gene %in% proclass    ~ "ProClass",
    TRUE                  ~ "Other"
  ))

# Panel colours

panel_colors <- c(
  "Oncotype Dx" = "purple",
  "Prolaris"    = "navyblue",
  "Decipher"    = "skyblue",
  "ProClass"    = "deeppink"
)

# Order genes within and across panels

panel_order <- c("Oncotype Dx", "Decipher", "Prolaris", "ProClass")

avg_df <- avg_df %>%
  filter(Panel != "Other") %>%
  mutate(Panel = factor(Panel, levels = panel_order))

avg_df <- avg_df %>%
  group_by(Panel) %>%
  arrange(desc(AvgProportion), .by_group = TRUE) %>%
  mutate(Gene = factor(Gene, levels = Gene)) %>%
  ungroup()

combined_gene_order <- avg_df %>%
  arrange(Panel, desc(AvgProportion)) %>%
  pull(Gene)

avg_df$Gene <- factor(avg_df$Gene, levels = combined_gene_order)


# Build combined bar plot

combined_plot <- ggplot(avg_df, aes(x = Gene, y = AvgProportion, fill = Panel)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = panel_colors) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = function(x) sub("\\.0+$", "", format(x, trim = TRUE))
  ) +
  labs(
    x = "Gene",
    y = "Average fraction of sample spots"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6, face = "italic"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "top"
  )

# Save plot

ggsave(
  filename = "/Path/To/Output_Directory>/Figure1C_AllPanels_Combined.png",
  plot     = combined_plot,
  width    = 6,
  height   = 3,
  dpi      = 600
)
