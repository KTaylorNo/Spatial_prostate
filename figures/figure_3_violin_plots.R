# # Figure 3 and Supplementary Figure 6
# Violin plots showing spot-level gene expression distribution across different annotation groups
# This code is an example for NFIB in sample H1_4 from patient 1

library(Seurat)
library(dplyr)
library(ggplot2)

# Load annotated Seurat object

seurat_obj <- readRDS(
  file = "</Path/To/Annotated_Seurat_Object>/P1_H1_4_annotated.rds"
)

# Extract NFIB expression from Spatial assay

counts <- GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")

if (!("NFIB" %in% rownames(counts))) {
  stop("NFIB not found in counts matrix.")
}

expr <- counts["NFIB", ]

# Collapse annotation categories

ann_raw <- seurat_obj$Final_Annotations

ann_grouped <- case_when(
  ann_raw == "Stroma" ~ "Stroma",
  ann_raw %in% c("Benign", "Benign*") ~ "Benign",
  ann_raw %in% c("GG4", "GG4 Cribriform") ~ "GG4",
  ann_raw == "GG2" ~ "GG2",
  ann_raw == "GG1" ~ "GG1",
  TRUE ~ "Other"
)

# Build dataframe and add whole-sample group

df <- data.frame(
  NFIB_Expression = expr,
  Annotation = ann_grouped
)

df_whole <- df %>% mutate(Annotation = "Whole\nsample")
df_combined <- bind_rows(df, df_whole)

# Annotation ordering and colors

annotation_levels <- c("Whole\nsample", "Stroma", "Other", "Benign", "GG1", "GG2", "GG4")

df_combined$Annotation <- factor(df_combined$Annotation, levels = annotation_levels)

custom_colors <- c(
  "Whole\nsample" = "gold",
  "Stroma"        = "#2ca02c",
  "Other"         = "grey70",
  "Benign"        = "#1f77b4",
  "GG1"           = "#ff7f0e",
  "GG2"           = "#d62728",
  "GG4"           = "#800000"
)

# Violin plot - showing density instead of individual points

p <- ggplot(df_combined, aes(x = Annotation, y = NFIB_Expression, fill = Annotation)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.7) +
  scale_fill_manual(values = custom_colors) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.y  = element_text(color = "black"),
    axis.text.x  = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "none",
    plot.title = element_blank()
  )

# Save figure

ggsave(
  filename = "</Path/To/Output>/NFIB_violin_P1_H1_4.png",
  plot     = p,
  width    = 5,
  height   = 3,
  dpi      = 600
)


# statistical tests for investigating annotation-specific expression

analyze_gene <- function(gene) {
  
  library(dplyr)
  library(effsize)
  library(Seurat)
  
  seurat_obj <- readRDS(
    file = "</Path/To/Annotated_Seurat_Object>/P1_H1_4_annotated.rds"
  )
  counts <- GetAssayData(seurat_obj, assay = "Spatial", layer = "counts")
  
  if (!(gene %in% rownames(counts))) {
    stop(paste("Gene", gene, "not found in counts matrix."))
  }
  
  expr <- counts[gene, ]
  
  ann_raw <- seurat_obj$Final_Annotations
  
  ann_grouped <- case_when(
    ann_raw == "Stroma" ~ "Stroma",
    ann_raw %in% c("Benign", "Benign*") ~ "Benign",
    ann_raw %in% c("GG4", "GG4 Cribriform") ~ "GG4",
    ann_raw == "GG2" ~ "GG2",
    ann_raw == "GG1" ~ "GG1",
    TRUE ~ "Other"
  )
  
  df <- data.frame(
    Expression = expr,
    Annotation = ann_grouped
  )
  
  # Remove annotation groups with zero total expression
  df <- df %>%
    group_by(Annotation) %>%
    filter(sum(Expression) > 0) %>%
    ungroup()

  cat("\nGene:", gene, "\n")
  cat("--------------------------------------------------\n")
  
  # Global Kruskal–Wallis test
  
  cat("Global test (Kruskal–Wallis):\n")
  kw <- kruskal.test(Expression ~ Annotation, data = df)
  print(kw)
  cat("\n")
  
  # Group medians
  
  group_medians <- df %>%
    group_by(Annotation) %>%
    summarise(median_expr = median(Expression))
  
  cat("Group medians:\n")
  print(group_medians)
  cat("\n")
  
  # Identify enriched group (highest median)
  top_group <- group_medians %>%
    filter(median_expr == max(median_expr)) %>%
    pull(Annotation)
  
  cat("Enriched group (highest median):", top_group, "\n\n")
  
  # Pairwise Wilcoxon tests vs enriched group

  cat("Pairwise Wilcoxon tests vs enriched group:\n")
  
  other_groups <- setdiff(unique(df$Annotation), top_group)
  
  results <- lapply(other_groups, function(g) {
    list(
      group = g,
      wilcox = wilcox.test(
        df$Expression[df$Annotation == top_group],
        df$Expression[df$Annotation == g]
      ),
      delta = effsize::cliff.delta(
        df$Expression[df$Annotation == top_group],
        df$Expression[df$Annotation == g]
      )
    )
  })
  
  # Print results
  for (res in results) {
    cat("\nComparison:", top_group, "vs", res$group, "\n")
    print(res$wilcox)
    print(res$delta)
  }
  
  invisible(results)
}
# enter gene name to show gene-specific results (here showing results for MSMB)
analyze_gene("MSMB")
