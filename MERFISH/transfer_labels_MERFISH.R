# Load libraries
library(Seurat)
library(dplyr)

# Step 0: Read in Seurat objects
# Make sure to replace these paths with your actual saved RDS paths
fem <- readRDS("../female_RNA/combined_seurat.RDS")
mer <- readRDS("../MERFISH/merfish_seurat.RDS")

# Step 1: Subsample RNA cells to a max of 2000 per annotated cluster
max_cells_per_cluster <- 2000
indices_to_keep <- integer(0)
unique_clusters <- unique(fem$transfer_celltypes)

for (cluster in unique_clusters) {
  cluster_indices <- which(fem$transfer_celltypes == cluster)
  if (length(cluster_indices) > max_cells_per_cluster) {
    max_c = max(max_cells_per_cluster, ceiling(length(cluster_indices)/10))
    cat(cluster, max_c, "\n")
    sampled_indices <- sample(cluster_indices, max_cells_per_cluster)
    indices_to_keep <- c(indices_to_keep, sampled_indices)
  } else {
    indices_to_keep <- c(indices_to_keep, cluster_indices)
  }
}
fem <- subset(fem, cells = colnames(fem)[indices_to_keep])

# Step 2: Run label transfer using CCA
anchors <- FindTransferAnchors(
  reference = fem,
  query = mer,
  dims = 1:45,
  reduction = "cca",
  features = rownames(mer),
  normalization.method = "SCT"
)

predictions <- TransferData(
  anchorset = anchors,
  refdata = fem$transfer_celltypes,
  dims = 1:45
)
mer <- AddMetaData(mer, metadata = predictions)

# Step 3: Subcluster MERFISH data and assign predicted labels per subcluster
allmeta <- list()

for (cl in unique(mer$seurat_clusters)) {
  sub <- subset(mer, subset = seurat_clusters == cl)
  sub <- NormalizeData(sub)
  sub <- FindVariableFeatures(sub, nfeatures = 500)
  sub <- ScaleData(sub)
  sub <- RunPCA(sub, features = VariableFeatures(sub))
  sub <- FindNeighbors(sub, dims = 1:25)
  sub <- FindClusters(sub, resolution = 2)
  sub <- RunUMAP(sub, reduction = "X_pca", dims = 1:25)
  sub$sub_leiden <- paste(cl, sub$seurat_clusters)

  # Assign predicted ID based on majority vote in each subcluster
  metaf <- sub@meta.data
  metaf <- metaf[!is.na(sub$predicted.id) & metaf$prediction.score.max > 0.85, ]

  predictions_table <- table(metaf$seurat_clusters, metaf$predicted.id)
  predictions_table <- predictions_table / rowSums(predictions_table)
  predictions_df <- as.data.frame(predictions_table)

  new_df <- predictions_df %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq) %>%
    as.data.frame()

  mat <- match(sub$seurat_clusters, new_df$Var1)
  sub$predicted_id_ext <- as.character(new_df$Var2[mat])

  # Save UMAP and metadata for each subcluster
  pdf(paste0(cl, "_sub.pdf"))
  print(DimPlot(sub, group.by = "seurat_clusters", label = TRUE))
  print(DimPlot(sub, group.by = "subclass_label", label = TRUE))
  print(DimPlot(sub, group.by = "predicted_id_ext", label = TRUE))
  print(DimPlot(sub, group.by = "predicted.id", label = TRUE))
  print(DimPlot(sub, group.by = "age"))
  print(FeaturePlot(sub, "log1p_total_counts", max.cutoff = 2000))
  dev.off()

  meta <- sub@meta.data[, c("age", "predicted.id", "prediction.score.max",
                            "predicted_id_ext", "subclass_label", "sub_leiden")]
  write.table(meta, paste0(cl, "_sub_meta.txt"), sep = "\t", quote = FALSE)
  allmeta[[cl]] <- meta
}
