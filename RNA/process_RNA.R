# RNA Processing Workflow using Seurat + SCT
# Author: Luisa Amaral (compiled and cleaned)
# Description: Consolidated pipeline for RNA-seq preprocessing, clustering,
# cell type annotation (using Allen Brain Institute reference), and differential expression.

# ==== Setup ====
library(Seurat)
library(tidyverse)
library(data.table)
library(Matrix)
library(patchwork)
library(sceasy)
library(reticulate)
library(ggplot2)

# ==== 1. Load Data (Post-DoubletFinder, from .RData files) ====
# This section loads pre-filtered Seurat objects where doublets were already removed.
data_dir <- "../Doublet_Finder/"
data_files <- list.files(data_dir, pattern = ".RData", full.names = TRUE)

# Load the first Seurat object and initialize the combined object
load(data_files[1])
combined_seurat <- obj

# Loop through remaining RData files and merge all objects
for (file in data_files[-1]) {
  load(file)
  seurat_obj <- obj
  sample_name <- gsub("../Doublet_Finder/", "", file)
  sample_name <- gsub(".RData", "", sample_name)
  seurat_obj$sample <- sample_name
  cat(sample_name, "\t")
  combined_seurat <- merge(combined_seurat, seurat_obj)
}

saveRDS(combined_seurat, file = "../female_RNA/combined_seurat.RDS")

# ==== 2. Quality Control ====
combined <- readRDS("../female_RNA/combined_seurat.RDS")
combined["percent.mt"] <- PercentageFeatureSet(combined, pattern = "^mt-")
VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
combined <- subset(combined, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10)

# ==== 3. Normalize with SCTransform ====
combined <- SCTransform(combined, verbose = FALSE)

# ==== 4. PCA, UMAP, Clustering ====
combined <- RunPCA(combined, verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:30)
combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
DimPlot(combined, reduction = "umap", label = TRUE)

# ==== 5. Save Pre-filtered ====
saveRDS(combined, file = "prefilter_RNA.rds")

# ==== 6. Manual Filtering (based on filtering_RNA.r) ====
# Read prefiltered object
combined <- readRDS("prefilter_RNA.rds")

# Remove predicted doublets
filtered_seurat <- subset(combined, subset = Combined_Classifications == "Singlet")

# Re-run normalization and clustering
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, nfeatures = 5000)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(filtered_seurat))
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:50)
filtered_seurat <- FindClusters(filtered_seurat, resolution = 3)
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:50)

# Manual filtering step:
# Subclusters with >50% low-quality cells (high doublet scores or low UMI) were manually identified and removed

# Save postfilter object
saveRDS(filtered_seurat, file = "postfilter3_RNA.RDS")

# ==== 7. Transfer Labels from Allen Brain Reference ====
# Allen data was preprocessed and converted to Seurat object using sceasy
allen <- readRDS("filtered_1000k_per_subclass.rds")

# Subsample Allen reference to max 250 cells per subclass
max_cells_per_cluster <- 250
indices_to_keep <- unlist(lapply(unique(allen$cl), function(cluster) {
  cluster_indices <- which(allen$cl == cluster)
  if (length(cluster_indices) > max_cells_per_cluster) sample(cluster_indices, max_cells_per_cluster)
  else cluster_indices
}))
allen <- subset(allen, cells = indices_to_keep)

# Process Allen reference
allen <- FindVariableFeatures(allen, nfeatures = 5000)
allen <- ScaleData(allen)
allen <- RunPCA(allen, features = VariableFeatures(allen))

# Label transfer
anchors <- FindTransferAnchors(reference = allen, query = filtered_seurat, dims = 1:40,
  reduction = "rpca", max.features = 250, k.anchor = 30, features = VariableFeatures(allen))
predictions <- TransferData(anchorset = anchors, refdata = allen$subclass_label, dims = 1:40, weight.reduction = "rpca.ref")
filtered_seurat <- AddMetaData(filtered_seurat, metadata = predictions)

# Save labeled Seurat object and metadata
saveRDS(filtered_seurat, file = "combined_final_SCT_labeled.rds")
write.csv(filtered_seurat@meta.data, "output/final_rna_metadata.csv")


# ==== 8. DEG with MAST ====

# Remove ribosomal genes (RPL and RPS)
all_genes <- rownames(filtered_seurat)
ribosomal_genes <- all_genes[grepl("^Rpl|^Rps", all_genes)]
cat("Ribosomal genes to be removed:
")
print(ribosomal_genes)
filtered_seurat <- subset(filtered_seurat, features = setdiff(all_genes, ribosomal_genes))

# === DEG by Cell Type ===
setwd("../results/DEG_celltype_filtered")
filtered_seurat$age_celltype <- paste(filtered_seurat$age, filtered_seurat$celltype_final, sep = "_")
Idents(filtered_seurat) <- "age_celltype"

for (celltype in unique(filtered_seurat$celltype_final)) {
  cat("Processing:", celltype, "
")
  ident1 <- gsub(" ", "_", paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_", paste("18mo", celltype, sep = "_"))

  if (sum(filtered_seurat$age_celltype == ident1) > 20 & sum(filtered_seurat$age_celltype == ident2) > 20) {
    de_results <- FindMarkers(
      filtered_seurat,
      ident.1 = ident1,
      ident.2 = ident2,
      test.use = "MAST",
      logfc.threshold = 0.1,
      min.pct = 0.01,
      latent.vars = c("region", "rep", "percent.mt", "percent.ribo")
    )
    file_name <- paste0(gsub("/", "-", gsub(" ", "_", celltype)), ".csv")
    write.csv(de_results, file = file_name, row.names = TRUE)
  }
}

# === DEG by Cell Type and Region ===
setwd("../results/DEG_celltype_region_filtered")
filtered_seurat$age_celltype_region <- gsub(" ", "_", paste(filtered_seurat$age, filtered_seurat$celltype_final, filtered_seurat$region, sep = "_"))
Idents(filtered_seurat) <- "age_celltype_region"

for (r in unique(filtered_seurat$region)) {
  for (celltype in unique(filtered_seurat$celltype_final)) {
    file_name <- paste0(gsub("/", "-", gsub(" ", "_", celltype)), "--", r, "_2vs18.csv")
    if (file.exists(file_name)) {
      cat(file_name, "file exists
")
      next
    }
    ident1 <- paste("2mo", celltype, r, sep = "_")
    ident2 <- paste("18mo", celltype, r, sep = "_")
    if (sum(filtered_seurat$age_celltype_region == ident1) > 20 & sum(filtered_seurat$age_celltype_region == ident2) > 20) {
      cat("Running", file_name, "
")
      de_results <- FindMarkers(
        filtered_seurat,
        ident.1 = ident1,
        ident.2 = ident2,
        test.use = "MAST",
        logfc.threshold = 0.1,
        min.pct = 0.01,
        latent.vars = c("rep", "percent.mt", "percent.ribo")
      )
      write.csv(de_results, file = file_name, row.names = TRUE)
    }
  }
}