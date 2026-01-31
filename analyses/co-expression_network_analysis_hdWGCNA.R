# hdWGCNA pipeline for each Visium sample

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(dplyr)
library(tidyr)

# Use cowplot theme globally
theme_set(theme_cowplot())

# Load Visium data

prostate <- Load10X_Spatial(
  data.dir = "</Path/To/Visium_Output>/outs"
)

# Basic filtering

# Remove low-quality spots
prostate <- subset(prostate, subset = nCount_Spatial >= 500)

# Remove genes detected in <5 spots
spots_per_gene <- rowSums(GetAssayData(prostate, assay = "Spatial", layer = "counts") > 0)
low_spot_genes <- names(spots_per_gene[spots_per_gene < 5])

prostate <- subset(
  prostate,
  features = setdiff(rownames(prostate), low_spot_genes)
)

# Add barcode column
prostate$barcode <- colnames(prostate)

# Add spatial metadata

tissue_positions <- read.csv(
  "</Path/To/Visium_Output>/outs/spatial/tissue_positions.csv"
)

tissue_positions <- subset(tissue_positions, barcode %in% prostate$barcode)

new_meta <- left_join(
  prostate@meta.data,
  tissue_positions,
  by = "barcode"
)

prostate$row      <- new_meta$array_row
prostate$col      <- new_meta$array_col
prostate$imagerow <- new_meta$pxl_row_in_fullres
prostate$imagecol <- new_meta$pxl_col_in_fullres

# Normalization and dimensionality reduction

prostate <- prostate %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

# Clustering and UMAP

prostate <- FindNeighbors(prostate, dims = 1:30)
prostate <- FindClusters(prostate)
prostate <- RunUMAP(prostate, dims = 1:30)

# Prepare for hdWGCNA

prostate <- SetupForWGCNA(
  prostate,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "vis"
)

prostate <- MetaspotsByGroups(
  prostate,
  group.by = "seurat_clusters",
  ident.group = "seurat_clusters",
  assay = "Spatial"
)

prostate <- NormalizeMetacells(prostate)

# Inspect metacell object
GetMetacellObject(prostate)

prostate <- SetDatExpr(
  prostate,
  group.by = NULL,
  group_name = NULL
)

# Soft-threshold selection

prostate <- TestSoftPowers(prostate)
plot_list <- PlotSoftPowers(prostate)

# Save soft-threshold plots
for (i in seq_along(plot_list)) {
  ggsave(
    filename = paste0("</Path/To/Output>/soft_power_plot_", i, ".png"),
    plot = plot_list[[i]],
    width = 8,
    height = 6
  )
}

# Network construction and module detection

prostate <- ConstructNetwork(
  prostate,
  tom_name = "test",
  overwrite_tom = TRUE
)

prostate <- ModuleEigengenes(prostate)
prostate <- ModuleConnectivity(prostate)

prostate <- ResetModuleNames(
  prostate,
  new_name = "SM"
)

# Figure 4B and Supplementary Figure 3B
# Show spatial modules on sample's spatial map

MEs     <- GetMEs(prostate)
modules <- GetModules(prostate)
mods <- setdiff(levels(modules$module), "grey")
prostate@meta.data <- cbind(prostate@meta.data, MEs)
p <- SpatialFeaturePlot(
  prostate,
  features        = mods,
  alpha           = c(0.9, 1),
  pt.size.factor  = 3,
  ncol            = 8
)
ggsave(
  filename = "</Path/To/Output>/spatial_feature_plot_modules.png",
  plot     = p,
  width    = 25,
  height   = 16
)

# Identify top marker genes and biomarker genes per spatial module (kME > 0.1)
# See Supplementary Table 2

# Marker genes of interest
marker_genes <- c(
  "PRAC1","HPN","PCAT14","AMACR","PCA3","MSMB","ACPP","RDH11","NKX3-1","KRT5",
  "KRT15","TRIM29","TP63","SLC14A1","MMP7","PIGR","LTF","CP","KRT7","CXCR4",
  "TRBC1","CD3D","LYZ","CD79A","EPAS1","EMP1","VWF","PECAM1","IFI27","DCN",
  "LUM","FBLN1","SFRP2","TAGLN","ACTA2","ACTG2","MYH11","MYL9"
)

# Reshape module table and extract top marker genes
top_marker_genes_per_module <- modules %>%
  pivot_longer(
    cols      = starts_with("kME_SM"),
    names_to  = "module_kME",
    values_to = "kME"
  ) %>%
  mutate(module = gsub("kME_", "", module_kME)) %>%
  filter(gene_name %in% marker_genes) %>%
  group_by(module) %>%
  arrange(desc(kME), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  select(gene_name, module, kME)

# Display results
print(top_marker_genes_per_module)

# Biomarker gene list
selected_genes <- c(
  "KLK2","AZGP1","TPM2","GSN","BGN","FLNC","COL1A1","SFRP4","SRD5A2","FAM13C",
  "GSTM2","TPX2","TK1","TOP2A","RRM2","CENPM","PCLAF","CDK1","NUSAP1","BIRC5",
  "PTTG1","CDC20","CDCA3","PBK","CDKN3","ASF1B","FOXM1","KIF20A","CENPF","PLK1",
  "PRC1","ORC6","ASPM","RAD54L","RAD51","DLGAP5","CDCA8","BUB1B","CEP55","KIF11",
  "SKA1","DTL","MCM10","ANO7","MYBPC1","NFIB","TNFRSF19","CAMK2N1","PBX1","PCDH7",
  "THBS2","LASP1","RABGAP1","UBE2C","EPPK1","ZWILCH","S1PR4","IQGAP3","THNSL2",
  "KLF13","MICA","C17orf97","CCDC163","ACOT1"
)

# Extract biomarker genes with kME > 0.1 in any module
biomarker_kME_filtered <- modules %>%
  filter(gene_name %in% selected_genes) %>%
  pivot_longer(
    cols      = starts_with("kME_SM"),
    names_to  = "module_kME",
    values_to = "kME"
  ) %>%
  mutate(module = gsub("kME_", "", module_kME)) %>%
  filter(kME > 0.1) %>%
  arrange(module, desc(kME)) %>%
  select(gene_name, module, kME)

# Display full table
print(biomarker_kME_filtered, n = Inf)

# Figure 4C
# Biomarker gene co-expression with marker genes in spatial networks
# Here showing gene MYBPC1 in sample H1_4
# This analysis is repeated across all biomarker genes in each annotated sample from patient 1

library(Seurat)
library(hdWGCNA)
library(tidygraph)
library(dplyr)

# Define gene sets
custom_genes <- c(
  "KLK2","AZGP1","TPM2","GSN","BGN","FLNC","COL1A1","SFRP4","SRD5A2","FAM13C",
  "GSTM2","TPX2","TK1","TOP2A","RRM2","CENPM","PCLAF","CDK1","NUSAP1","BIRC5",
  "PTTG1","CDC20","CDCA3","PBK","CDKN3","ASF1B","FOXM1","KIF20A","CENPF","PLK1",
  "PRC1","ORC6","ASPM","RAD54L","RAD51","DLGAP5","CDCA8","BUB1B","CEP55","KIF11",
  "SKA1","DTL","MCM10","ANO7","MYBPC1","NFIB","TNFRSF19","CAMK2N1","PBX1","PCDH7",
  "THBS2","LASP1","RABGAP1","UBE2C","EPPK1","ZWILCH","S1PR4","IQGAP3","TSBP1",
  "THNSL2","KLF13","MICA","C17orf97","CCDC163","ACOT1",
  "PRAC1","HPN","PCAT14","AMACR","PCA3","MSMB","ACPP","RDH11","NKX3-1","KRT5",
  "KRT15","TRIM29","TP63","SLC14A1","MMP7","PIGR","LTF","CP","KRT7","CXCR4",
  "TRBC1","CD3D","LYZ","CD79A","EPAS1","EMP1","VWF","PECAM1","IFI27","DCN",
  "LUM","FBLN1","SFRP2","TAGLN","ACTA2","ACTG2","MYH11","MYL9"
)

marker_genes <- c(
  "PRAC1","HPN","PCAT14","AMACR","PCA3","MSMB","ACPP","RDH11","NKX3-1","KRT5",
  "KRT15","TRIM29","TP63","SLC14A1","MMP7","PIGR","LTF","CP","KRT7","CXCR4",
  "TRBC1","CD3D","LYZ","CD79A","EPAS1","EMP1","VWF","PECAM1","IFI27","DCN",
  "LUM","FBLN1","SFRP2","TAGLN","ACTA2","ACTG2","MYH11","MYL9"
)

# Extract modules and topoligcal overlap matrix (TOM)

modules <- GetModules(prostate) %>%
  filter(module != "grey") %>%
  mutate(module = droplevels(module))

TOM <- GetTOM(prostate)

# Subset TOM to custom genes

cur_genes <- intersect(custom_genes, modules$gene_name)
cur_TOM   <- TOM[cur_genes, cur_genes]

# Build graph object

graph <- cur_TOM %>%
  igraph::graph_from_adjacency_matrix(
    mode     = "undirected",
    weighted = TRUE
  ) %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(
    module = modules$module[match(name, modules$gene_name)],
    size   = 5
  )

# Compute marker‑gene connectivity

marker_genes_in_graph <- intersect(
  marker_genes,
  graph %>% activate(nodes) %>% pull(name)
)

edges_tbl <- graph %>%
  activate(edges) %>%
  as_tibble() %>%
  mutate(
    from_name = graph %>% activate(nodes) %>% pull(name) %>% .[from],
    to_name   = graph %>% activate(nodes) %>% pull(name) %>% .[to]
  )

marker_connectivity <- edges_tbl %>%
  filter(from_name %in% marker_genes_in_graph |
         to_name   %in% marker_genes_in_graph) %>%
  mutate(
    partner = ifelse(from_name %in% marker_genes_in_graph, from_name, to_name)
  ) %>%
  group_by(partner) %>%
  summarise(total_connectivity = sum(weight), .groups = "drop")

# Marker‑gene edges for a given biomarker gene (here showing MYBPC1)

mybpc1_edges <- edges_tbl %>%
  filter(
    (from_name == "MYBPC1" & to_name %in% marker_genes_in_graph) |
    (to_name   == "MYBPC1" & from_name %in% marker_genes_in_graph)
  ) %>%
  mutate(
    partner_gene = ifelse(from_name == "MYBPC1", to_name, from_name)
  )

# Normalize and rank MYBPC1 connections to identify strongest marker gene type connection (e.g. benign or tumour)

mybpc1_normalized <- mybpc1_edges %>%
  left_join(marker_connectivity, by = c("partner_gene" = "partner")) %>%
  mutate(
    normalized_weight = weight / total_connectivity
  ) %>%
  arrange(desc(normalized_weight)) %>%
  select(partner_gene, weight, total_connectivity, normalized_weight)

print(mybpc1_normalized)
