# Supplementary Figure 9
# Average spot coverage and expression of biomarker genes
# across annotation groups (Patient 1)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define biomarker gene list

marker_genes <- c(
  "KLK2","AZGP1","TPM2","GSN","BGN","FLNC","COL1A1","SFRP4","SRD5A2","FAM13C",
  "GSTM2","TPX2","TK1","TOP2A","RRM2","CENPM","PCLAF","CDK1","NUSAP1","BIRC5",
  "PTTG1","CDC20","CDCA3","PBK","CDKN3","ASF1B","FOXM1","KIF20A","CENPF","PLK1",
  "PRC1","ORC6","ASPM","RAD54L","RAD51","DLGAP5","CDCA8","BUB1B","CEP55","KIF11",
  "SKA1","DTL","MCM10","ANO7","MYBPC1","NFIB","TNFRSF19","CAMK2N1","PBX1","PCDH7",
  "THBS2","LASP1","RABGAP1","UBE2C","EPPK1","ZWILCH","S1PR4","IQGAP3","THNSL2",
  "KLF13","MICA","C17orf97","CCDC163","ACOT1"
)

# 2. Patient 1 samples with annotations

samples <- c(
  "P1_H1_2","P1_H1_4","P1_H1_5",
  "P1_H2_1","P1_H2_2","P1_H2_5",
  "P1_V1_2"
)

# Load and pool samples

all_spots <- list()

for (sample in samples) {
  
  file_path <- paste0("</Path/To/Annotated_Seurat_Objects>/", sample, "_annotated.rds")
  if (!file.exists(file_path)) next
  
  obj <- readRDS(file_path)
  DefaultAssay(obj) <- "SCT"
  
  genes_present <- intersect(marker_genes, rownames(obj))
  if (length(genes_present) == 0) next
  
  expr <- GetAssayData(obj, layer = "data")[genes_present, , drop = FALSE]
  ann_raw <- obj$Final_Annotations
  
  ann_clean <- case_when(
    ann_raw %in% c("NA","Exclude","",NA) ~ NA_character_,
    ann_raw %in% c("Benign","Benign*") ~ "Benign",
    ann_raw %in% c("GG4","GG4 Cribriform") ~ "GG4",
    ann_raw %in% c("Vessel","Transition_State","PIN","Nerve",
                   "Inflammation","Fat","Chronic inflammation") ~ "Other",
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

# Keep only biomarker genes present in the dataset

marker_genes <- intersect(marker_genes, colnames(combined_df))

# Compute % expressing and average expression

dot_df <- combined_df %>%
  pivot_longer(
    cols = all_of(marker_genes),
    names_to = "Gene",
    values_to = "Expression"
  ) %>%
  group_by(Annotation, Gene) %>%
  summarise(
    pct_expr = mean(Expression > 0) * 100,
    avg_expr = mean(Expression),
    .groups = "drop"
  )

# Fill missing annotation–gene combinations
all_annotations <- unique(combined_df$Annotation)

complete_df <- expand.grid(
  Annotation = all_annotations,
  Gene       = marker_genes,
  stringsAsFactors = FALSE
)

dot_df <- complete_df %>%
  left_join(dot_df, by = c("Annotation","Gene")) %>%
  mutate(
    pct_expr = replace_na(pct_expr, 0),
    avg_expr = replace_na(avg_expr, 0)
  )

# Factor ordering
desired_order <- c("Other","Stroma","Benign","GG1","GG2","GG4")

dot_df$Gene       <- factor(dot_df$Gene, levels = rev(marker_genes))
dot_df$Annotation <- factor(dot_df$Annotation, levels = rev(desired_order))

# Italic gene labels
italic_labels <- setNames(
  lapply(marker_genes, function(g) bquote(italic(.(g)))),
  marker_genes
)

# Full dot‑plot

p <- ggplot(dot_df, aes(x = Annotation, y = Gene)) +
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_size(range = c(1, 8)) +
  scale_color_viridis_c(option = "magma") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(labels = italic_labels) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x.top    = element_text(angle = 45, hjust = 0, vjust = 0.5),
    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.title         = element_blank()
  ) +
  labs(size = "% expressing", color = "Avg expression")

ggsave(
  filename = "</Path/To/Output>/marker_gene_dotplot_patient1.png",
  plot     = p,
  width    = 5,
  height   = 10,
  dpi      = 300
)

# Panel definitions

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

annotation_levels_flipped <- c("GG4","GG2","GG1","Benign","Other","Stroma")

# Panel dot‑plot function

make_dotplot <- function(panel_name, gene_list, dot_df) {
  
  df <- dot_df %>% filter(Gene %in% gene_list)
  
  italic_labels <- setNames(
    lapply(gene_list, function(g) bquote(italic(.(g)))),
    gene_list
  )
  
  ggplot(df, aes(x = Gene, y = Annotation)) +
    geom_point(aes(size = pct_expr, color = avg_expr)) +
    scale_size(range = c(1, 8), limits = c(0, 100),
               breaks = c(0, 25, 75, 100), name = "% expressing") +
    scale_color_viridis_c(option = "magma", limits = c(0, 5),
                          name = "Avg expression") +
    scale_x_discrete(limits = gene_list, labels = italic_labels,
                     position = "bottom") +
    scale_y_discrete(limits = annotation_levels_flipped) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.x.top    = element_blank(),
      axis.ticks.x.top   = element_blank(),
      axis.title         = element_blank()
    ) +
    labs(title = paste(panel_name, "marker gene expression"))
}

# Save panel plots

ggsave("</Path/To/Output>/DotPlot_OncotypeDx.png",
       make_dotplot("Oncotype Dx", oncotype_dx, dot_df),
       width = 6, height = 4, dpi = 300)

ggsave("</Path/To/Output>/DotPlot_Prolaris.png",
       make_dotplot("Prolaris", prolaris, dot_df),
       width = 10, height = 4, dpi = 300)

ggsave("</Path/To/Output>/DotPlot_Decipher.png",
       make_dotplot("Decipher", decipher, dot_df),
       width = 6, height = 4, dpi = 300)

ggsave("</Path/To/Output>/DotPlot_ProClass.png",
       make_dotplot("ProClass", proclass, dot_df),
       width = 4, height = 4, dpi = 300)
