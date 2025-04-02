# RNA TE Expression Aggregation from SoloTE Output (Subfamily Level)

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(enrichR)

setwd("~/projects/fmouse_multiome/SoloTE_out/")

# Get directories containing SoloTE results (subfamily level)
all_directories <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
data_directories <- all_directories[grep("_subfamily", all_directories)]

# Load metadata
meta <- read.csv("~/projects/combined_all/female_RNA/meta_final.csv")

# Initialize merged Seurat object
init_dir <- data_directories[1]
setwd(init_dir)
sample <- basename(init_dir)
solote_matrix <- ReadMtx("matrix.mtx", "barcodes.tsv", "features.tsv", feature.column = 1)
merged_seuratobj <- CreateSeuratObject(counts = solote_matrix, min.features = 100, project = sample)
merged_seuratobj$sample <- gsub("_SoloTE_output", "", merged_seuratobj$orig.ident)
merged_seuratobj$barcode <- paste(merged_seuratobj$sample, sapply(strsplit(rownames(merged_seuratobj@meta.data), "-"), "[[", 1), sep = ":")
merged_seuratobj$barcode <- gsub("02mo", "2mo", merged_seuratobj$barcode)
merged_seuratobj$barcode <- gsub("09mo", "9mo", merged_seuratobj$barcode)
mat <- match(merged_seuratobj$barcode, meta$barcode)
merged_seuratobj$celltype <- meta$celltype_final[mat]
merged_seuratobj <- merged_seuratobj[, !is.na(merged_seuratobj$celltype)]

# Loop through remaining directories
for (directory in data_directories[-1]) {
  setwd("~/projects/fmouse_multiome/SoloTE_out/")
  setwd(directory)
  sample <- basename(directory)
  cat("Processing:", sample, "\n")
  solote_matrix <- ReadMtx("matrix.mtx", "barcodes.tsv", "features.tsv", feature.column = 1)
  if (ncol(solote_matrix) < 25) next

  temp_seuratobj <- CreateSeuratObject(counts = solote_matrix, min.features = 100, project = sample)
  temp_seuratobj$sample <- gsub("_SoloTE_output", "", temp_seuratobj$orig.ident)
  temp_seuratobj$barcode <- paste(temp_seuratobj$sample, sapply(strsplit(rownames(temp_seuratobj@meta.data), "-"), "[[", 1), sep = ":")
  temp_seuratobj$barcode <- gsub("02mo", "2mo", temp_seuratobj$barcode)
  temp_seuratobj$barcode <- gsub("09mo", "9mo", temp_seuratobj$barcode)
  mat <- match(temp_seuratobj$barcode, meta$barcode)
  temp_seuratobj$celltype <- meta$celltype_final[mat]
  temp_seuratobj <- temp_seuratobj[, !is.na(temp_seuratobj$celltype)]

  merged_seuratobj <- merge(merged_seuratobj, temp_seuratobj)
}

# Assign to "seuratobj" and join layers
seuratobj <- merged_seuratobj
seuratobj <- JoinLayers(seuratobj)

# QC Plot
merged_seuratobj@meta.data %>% 
  ggplot(aes(color = sample, x = nFeature_RNA, fill = sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() + 
  scale_x_log10() + 
  geom_vline(xintercept = 100)

# Normalization and Dimensionality Reduction
seuratobj <- NormalizeData(seuratobj)
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
seuratobj <- ScaleData(seuratobj, features = rownames(seuratobj))
seuratobj <- RunPCA(seuratobj, features = VariableFeatures(seuratobj))
seuratobj <- FindNeighbors(seuratobj, dims = 1:45)
seuratobj <- FindClusters(seuratobj, resolution = 1)
seuratobj <- RunUMAP(seuratobj, dims = 1:45)
DimPlot(seuratobj, reduction = "umap")

# Save final object
saveRDS(seuratobj, '~/projects/combined_all/female_RNA/SoloTE/all_subfamily.RDS')
