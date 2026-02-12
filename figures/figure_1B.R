# Figure 1B
# Proportion of tissue annotation types per sample

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# Define sample names and annotation types

sample_names <- c(
  "P1_H1_2", "P1_H1_4", "P1_H1_5",
  "P1_H2_5", "P1_H2_1", "P1_H2_2",
  "P1_V1_2"
)

# Annotation order (top → bottom in legend)
all_annotations <- c("Stroma", "Other", "Benign", "GG1", "GG2", "GG4")

# Matching colors
custom_colors <- c(
  "Stroma" = "#2ca02c",
  "Other"  = "grey70",
  "Benign" = "#1f77b4",
  "GG1"    = "#ff7f0e",
  "GG2"    = "#d62728",
  "GG4"    = "#800000"
)

# Compute annotation proportions per sample

all_samples_df <- list()
cancer_fractions <- tibble(Sample = character(), CancerFraction = numeric())

for (sample in sample_names) {
  
  seurat_path <- paste0("</Path/To/Annotated_Seurat_Objects>/", sample, "_annotated.rds")
  seurat_obj  <- readRDS(seurat_path)
  
  total_spots <- ncol(seurat_obj)
  
  grouped_annotations <- case_when(
    seurat_obj$Final_Annotations == "Stroma" ~ "Stroma",
    seurat_obj$Final_Annotations %in% c("Benign", "Benign*") ~ "Benign",
    seurat_obj$Final_Annotations == "GG1" ~ "GG1",
    seurat_obj$Final_Annotations == "GG2" ~ "GG2",
    seurat_obj$Final_Annotations %in% c("GG4", "GG4 Cribriform") ~ "GG4",
    TRUE ~ "Other"
  )
  
  annotation_df <- as.data.frame(table(grouped_annotations)) %>%
    rename(Annotation = grouped_annotations, Count = Freq) %>%
    mutate(Proportion = Count / total_spots) %>%
    complete(Annotation = all_annotations, fill = list(Count = 0, Proportion = 0))
  
  annotation_df$Sample     <- sample
  annotation_df$Annotation <- factor(annotation_df$Annotation, levels = all_annotations)
  
  all_samples_df[[sample]] <- annotation_df
  
  cancer_fraction <- annotation_df %>%
    filter(Annotation %in% c("GG1", "GG2", "GG4")) %>%
    summarise(CancerFraction = sum(Proportion)) %>%
    mutate(Sample = sample)
  
  cancer_fractions <- bind_rows(cancer_fractions, cancer_fraction)
}

# Sort samples by cancer fraction (descending)

sorted_samples <- cancer_fractions %>%
  arrange(desc(CancerFraction)) %>%
  pull(Sample)

combined_df <- bind_rows(all_samples_df)

combined_df$Sample_clean <- gsub("^P1_", "", combined_df$Sample)
combined_df$Sample_clean <- factor(
  combined_df$Sample_clean,
  levels = gsub("^P1_", "", sorted_samples)
)

# Plot stacked barplot

p <- ggplot(combined_df, aes(x = Sample_clean, y = Proportion, fill = Annotation)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    labels = function(x) ifelse(x %% 1 == 0, as.integer(x), x)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  ) +
  ylab("Proportion of tissue annotations") +
  guides(fill = guide_legend(title = "Annotation"))

# Save figure

ggsave(
  filename = "/Path/To/Output>/Figure1B_stacked_barplot.png",
  plot     = p,
  width    = 7,
  height   = 4,
  dpi      = 600
)
