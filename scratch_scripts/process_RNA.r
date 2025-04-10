library(Seurat)
library(ggplot2)

# Directory containing your .RData files
data_dir <- "/oasis/tscc/scratch/smamde/Doublet_Corrected_Mouse_RNA/Doublet_Removal_Again"

# List all .RData files in the directory
data_files <- list.files(data_dir, pattern = ".RData", full.names = TRUE)


# Load the first Seurat object outside of the loop
load(data_files[1])
combined_seurat <- obj  # Initialize with the first Seurat object

# Loop through the files and load them into the combined Seurat object
for (file in data_files[-1]) {
  # Load the Seurat object from the .RData file
  load(file)
  seurat_obj <- obj
  # Add a unique sample name to each Seurat object (assuming you have a sample name)
  sample_name <- gsub("/oasis/tscc/scratch/smamde/Doublet_Corrected_Mouse_RNA/Doublet_Removal_Again/", "", file)
  sample_name <- gsub(".RData", "", sample_name)
  seurat_obj$sample <- sample_name
  cat(sample_name, "\t")
  # Add the Seurat object to the combined Seurat object
  combined_seurat <- merge(combined_seurat, seurat_obj)
}

saveRDS(combined_seurat, file = "/home/lamaral/oasis/female_RNA/combined_seurat.RDS")

combined_seurat <- NormalizeData(combined_seurat)

# Perform standard preprocessing steps (e.g., clustering, UMAP)
combined_seurat <- FindVariableFeatures(combined_seurat, selection.method = "vst", nfeatures = 5000)
combined_seurat <- ScaleData(combined_seurat, features = VariableFeatures(object = combined_seurat))
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(object = combined_seurat))
combined_seurat <- FindNeighbors(combined_seurat, dims = 1:45)
combined_seurat <- FindClusters(combined_seurat)
combined_seurat <- RunUMAP(combined_seurat,dims = 1:45)

# You can also perform additional analyses or visualizations here

# Save the combined Seurat object to a file for future use
saveRDS(combined_seurat, file = "/home/lamaral/oasis/female_RNA/combined_seurat_UMAP.RDS")


combined_seurat <- FindClusters(combined_seurat, resolution = 3)


combined_seurat = readRDS(file = "/home/lamaral/oasis/female_RNA/combined_seurat_UMAP.RDS")


head(combined_seurat@meta.data,5)

library(dplyr)
classification_columns <- grep("^DF.classifications_", colnames(combined_seurat@meta.data), value = TRUE)
metadata = combined_seurat@meta.data
metadata <- metadata %>%
  rowwise() %>%
  mutate(Combined_Classifications = coalesce(!!!syms(classification_columns)))

# Remove the individual classification columns
metadata <- metadata %>% select(-one_of(classification_columns))



combined_seurat@meta.data$combined_doublet_class = metadata$Combined_Classifications

table( metadata$Combined_Classifications)

options(repr.plot.width=14, repr.plot.height=8)

DimPlot(combined_seurat, group.by = "combined_doublet_class")


DimPlot(combined_seurat, group.by = "seurat_clusters", label = T)


table(combined_seurat@meta.data$seurat_clusters)

metadata = metadata[,-grep("pANN", colnames(metadata))]

metadata = metadata[,-grep("sample", colnames(metadata))]

combined_seurat@meta.data = metadata

head(metadata)

ggplot(metadata, aes(fill = orig.ident , x = factor(seurat_clusters))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(metadata, aes(fill = Combined_Classifications , x = factor(seurat_clusters))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

FeaturePlot(combined_seurat, "Ddit4", raster =F)



