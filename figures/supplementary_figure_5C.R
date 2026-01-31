# Supplementary Figure 5C
# Number of Visium spots per sample remaining after filtering

library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)

# Locate Visium Seurat objects

rds_files <- list.files(
  path       = "</Path/To/Seurat_Objects>/",
  pattern    = "visium.*\\.rds$",
  full.names = TRUE
)

# Count number of spots per sample

spot_counts <- map_dfr(rds_files, function(f) {
  obj <- readRDS(f)
  sample_id <- sub("^(P[12]_.*)_visium.*\\.rds$", "\\1", basename(f))
  tibble(Sample = sample_id, NumSpots = ncol(obj))
})

# Plot spot counts

p <- ggplot(spot_counts, aes(x = Sample, y = NumSpots)) +
  geom_col(fill = "steelblue") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Number of Spots per Sample",
    x     = "Sample",
    y     = "Number of Spots after Filtering"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save figure

ggsave(
  filename = "</Path/To/Output>/SpotCounts_per_sample.png",
  plot     = p,
  width    = 10,
  height   = 2,
  dpi      = 600
)
