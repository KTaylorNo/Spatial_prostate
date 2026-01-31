#To calculate the baseline spot coverage metric, shown in Figure 1C.

tissue_section <- "P1_V1_2"
genes <- c("KLK2", "AZGP1", "TPM2", "GSN", "BGN", "FLNC", "COL1A1", "SFRP4", "SRD5A2", "FAM13C", 
           "GSTM2", "TPX2", "TK1", "TOP2A", "RRM2", "CENPM", "PCLAF", "CDK1", "NUSAP1", "BIRC5", 
           "PTTG1", "CDC20", "CDCA3", "PBK", "CDKN3", "ASF1B", "FOXM1", "KIF20A", "CENPF", "PLK1", 
           "PRC1", "ORC6", "ASPM", "RAD54L", "RAD51", "DLGAP5", "CDCA8", "BUB1B", "CEP55", "KIF11", 
           "SKA1", "DTL", "MCM10", "ANO7", "MYBPC1", "NFIB", "TNFRSF19", "CAMK2N1", "PBX1", "PCDH7", 
           "THBS2", "LASP1", "RABGAP1", "UBE2C", "EPPK1", "ZWILCH", "S1PR4", "IQGAP3", "THNSL2", 
           "KLF13", "MICA", "C17orf97", "CCDC163")

prostate <- readRDS(file = </Path/To/tissue_section_Seurat_object>)
counts <- GetAssayData(prostate, assay = "Spatial", layer = "counts")
genes_present <- genes[genes %in% rownames(counts)]
proportions <- sapply(genes_present, function(gene) {
  sum(counts[gene, ] > 0) / ncol(counts)
})
df <- data.frame(
  Sample = tissue_section,
  Gene = genes_present,
  Proportion = round(proportions, 4)
)
