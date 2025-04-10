library(ggplot2)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)

obj = readRDS("../../female_RNA/RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_final"


library(Seurat)
library(ggplot2)

# 1. Subset the Seurat object for the desired cell type
dg_obj_hca <- subset(obj, subset = (celltype_final %in% c("DG Glut") & region == "HCA"))

# 2. Preprocess the data
# If your object is not already normalized, then:
dg_obj_hca <- NormalizeData(dg_obj_hca)

# Identify variable features
dg_obj_hca <- FindVariableFeatures(dg_obj_hca)

# Optionally, you can visualize variable features
# VariableFeaturePlot(dg_obj)

# Scale the data (all genes or a subset)
dg_obj_hca <- ScaleData(dg_obj_hca)

# 3. Run PCA for dimensionality reduction
dg_obj_hca <- RunPCA(dg_obj_hca, features = VariableFeatures(object = dg_obj_hca))

# Optionally, check PCA results
# VizDimLoadings(gaba_obj, dims = 1:2, reduction = "pca")
# DimPlot(gaba_obj, reduction = "pca")

# 4. Determine the number of dimensions to use (here, we use the first 10)
dims_to_use <- 1:10

# 5. Find neighbors and clusters
dg_obj_hca <- FindNeighbors(dg_obj_hca, dims = dims_to_use)
dg_obj_hca <- FindClusters(dg_obj_hca, resolution = 0.5)  # Adjust resolution as needed

# 6. Run UMAP for visualization
dg_obj_hca <- RunUMAP(dg_obj_hca, dims = dims_to_use)

# 7. Visualize the clusters with UMAP
DimPlot(dg_obj_hca, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of STR_D12_Gaba Subsampled Cells")


library(Seurat)
library(ggplot2)

# 1. Subset the Seurat object for the desired cell type
dg_obj_hcp <- subset(obj, subset = (celltype_final %in% c("DG Glut") & region == "HCP"))

# 2. Preprocess the data
# If your object is not already normalized, then:
dg_obj_hcp <- NormalizeData(dg_obj_hcp,layer = 'counts')

# Identify variable features
dg_obj_hcp <- FindVariableFeatures(dg_obj_hcp)

# Optionally, you can visualize variable features
# VariableFeaturePlot(dg_obj)

# Scale the data (all genes or a subset)
dg_obj_hcp <- ScaleData(dg_obj_hcp)

# 3. Run PCA for dimensionality reduction
dg_obj_hcp <- RunPCA(dg_obj_hcp, features = VariableFeatures(object = dg_obj_hcp))

# Optionally, check PCA results
# VizDimLoadings(gaba_obj, dims = 1:2, reduction = "pca")
# DimPlot(gaba_obj, reduction = "pca")

# 4. Determine the number of dimensions to use (here, we use the first 10)
dims_to_use <- 1:10

# 5. Find neighbors and clusters
dg_obj_hcp <- FindNeighbors(dg_obj_hcp, dims = dims_to_use)
dg_obj_hcp <- FindClusters(dg_obj_hcp, resolution = 0.5)  # Adjust resolution as needed

# 6. Run UMAP for visualization
dg_obj_hcp <- RunUMAP(dg_obj_hcp, dims = dims_to_use)

# 7. Visualize the clusters with UMAP
DimPlot(dg_obj_hcp, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of STR_D12_Gaba Subsampled Cells")


DimPlot(dg_obj_hcp, reduction = "umap", label = TRUE, group.by = "batch")


dg_obj_hcp$batch = "no"
dg_obj_hcp$batch[which(dg_obj_hcp$sample %in% c("HCP_9mo_3" , "HCP_2mo_2"))] = "weird"

library(Seurat)
library(ggplot2)

# 1. Split the object by batch (assuming metadata column "sample")
obj_list <- SplitObject(dg_obj_hcp, split.by = "batch")

# 2. Normalize and find variable features for each batch
for (i in 1:length(obj_list)) {
  obj_list[[i]] <- NormalizeData(obj_list[[i]])
  obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], selection.method = "vst", nfeatures = 2000)
}

# 3. Find integration anchors
# You can adjust dims; here, we use the first 30 PCs.
anchors <- FindIntegrationAnchors(object.list = obj_list, dims = 1:30)

# 4. Integrate the data
obj_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Set the default assay to "integrated"
DefaultAssay(obj_integrated) <- "integrated"

# 5. Run the standard workflow on the integrated data
obj_integrated <- ScaleData(obj_integrated)
obj_integrated <- RunPCA(obj_integrated, npcs = 30)
obj_integrated <- RunUMAP(obj_integrated, dims = 1:30)
obj_integrated <- FindNeighbors(obj_integrated, dims = 1:30)
obj_integrated <- FindClusters(obj_integrated, resolution = 0.5)

# 6. Visualize the integrated data
DimPlot(obj_integrated, reduction = "umap", label = TRUE) +
  ggtitle("UMAP After Batch Correction")


obj_integrated
options(repr.plot.width=5, repr.plot.height=5)

DimPlot(obj_integrated, group.by = "age",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


obj_integrated
options(repr.plot.width=5, repr.plot.height=5)

DimPlot(obj_integrated, group.by = "sample",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=5, repr.plot.height=5)

DimPlot(dg_obj_hcp, group.by = "orig.ident",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=5, repr.plot.height=5)

DimPlot(dg_obj_hcp, group.by = "age",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=5, repr.plot.height=5)

DimPlot(dg_obj_hcp, group.by = "age",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


FeaturePlot(dg_obj_hca, "Meg3",reduction = "umap", label = TRUE) 


library(Seurat)
library(ggplot2)

# 1. Subset the Seurat object for the desired cell type
dg_obj_hcp <- subset(obj, subset = (celltype_final %in% c("DG Glut") & region == "HCP"))

# 2. Preprocess the data
# If your object is not already normalized, then:
dg_obj_hcp <- NormalizeData(dg_obj_hcp)

# Identify variable features
dg_obj_hcp <- FindVariableFeatures(dg_obj_hcp)

# Optionally, you can visualize variable features
# VariableFeaturePlot(dg_obj)

# Scale the data (all genes or a subset)
dg_obj_hcp <- ScaleData(dg_obj_hcp)

# 3. Run PCA for dimensionality reduction
dg_obj_hcp <- RunPCA(dg_obj_hcp, features = VariableFeatures(object = dg_obj_hcp))

# Optionally, check PCA results
# VizDimLoadings(gaba_obj, dims = 1:2, reduction = "pca")
# DimPlot(gaba_obj, reduction = "pca")

# 4. Determine the number of dimensions to use (here, we use the first 10)
dims_to_use <- 1:10

# 5. Find neighbors and clusters
dg_obj_hcp <- FindNeighbors(dg_obj_hcp, dims = dims_to_use)
dg_obj_hcp <- FindClusters(dg_obj_hcp, resolution = 0.5)  # Adjust resolution as needed

# 6. Run UMAP for visualization
dg_obj_hcp <- RunUMAP(dg_obj_hcp, dims = dims_to_use)

# 7. Visualize the clusters with UMAP
DimPlot(dg_obj_hcp, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of STR_D12_Gaba Subsampled Cells")


library(Seurat)
library(ggplot2)

# 1. Subset the Seurat object for the desired cell type
dg_obj <- subset(obj, subset = celltype_final %in% c("DG Glut", "DG-PIR_Ex IMN"))

# 2. Preprocess the data
# If your object is not already normalized, then:
dg_obj <- NormalizeData(dg_obj)

# Identify variable features
dg_obj <- FindVariableFeatures(dg_obj)

# Optionally, you can visualize variable features
# VariableFeaturePlot(dg_obj)

# Scale the data (all genes or a subset)
dg_obj <- ScaleData(dg_obj)

# 3. Run PCA for dimensionality reduction
dg_obj <- RunPCA(dg_obj, features = VariableFeatures(object = dg_obj))

# Optionally, check PCA results
# VizDimLoadings(gaba_obj, dims = 1:2, reduction = "pca")
# DimPlot(gaba_obj, reduction = "pca")

# 4. Determine the number of dimensions to use (here, we use the first 10)
dims_to_use <- 1:10

# 5. Find neighbors and clusters
dg_obj <- FindNeighbors(dg_obj, dims = dims_to_use)
dg_obj <- FindClusters(dg_obj, resolution = 0.5)  # Adjust resolution as needed

# 6. Run UMAP for visualization
dg_obj <- RunUMAP(dg_obj, dims = dims_to_use)

# 7. Visualize the clusters with UMAP
DimPlot(dg_obj, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of STR_D12_Gaba Subsampled Cells")


saveRDS(dg_obj , "dg_obj.RDS")

options(repr.plot.width=5, repr.plot.height=5)

DimPlot(dg_obj_hca, group.by = "age",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=5, repr.plot.height=5)

DimPlot(dg_obj_hca, group.by = "orig.ident",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=5, repr.plot.height=5)
dg_obj_hca <- FindClusters(dg_obj_hca, resolution = 0.25)  # Adjust resolution as needed

DimPlot(dg_obj_hca, group.by = "seurat_clusters",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=5, repr.plot.height=5)
dg_obj_hcp <- FindClusters(dg_obj_hcp, resolution = 0.05)  # Adjust resolution as needed

DimPlot(dg_obj_hcp, group.by = "seurat_clusters",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=19, repr.plot.height=4)

FeaturePlot(dg_obj_hcp, "Tsix",reduction = "umap", label = TRUE)


FeaturePlot(dg_obj_hcp, features = c("Rps8", "Rpl13", "Malat1"))


fam = FindAllMarkers(dg_obj_hcp)

head(fam[which(fam$cluster==0),],20)

options(repr.plot.width=5, repr.plot.height=5)

DimPlot(dg_obj_hcp, group.by = "orig.ident",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=5, repr.plot.height=5)

DimPlot(dg_obj_hca, group.by = "seurat_clusters",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=19, repr.plot.height=4)

FeaturePlot(dg_obj_hcp, "AC149090.1",split.by = "orig.ident",reduction = "umap", label = TRUE)


options(repr.plot.width=15, repr.plot.height=5)

FeaturePlot(dg_obj_hca, "Nrg3os",split.by = "orig.ident",reduction = "umap", label = TRUE) 

VlnPlot(dg_obj_hca, "Nrg3os",group.by = "orig.ident") 

VlnPlot(dg_obj_hca, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))


dg_obj_hca = AddModuleScore(dg_obj_hca , features = rownames(dg_obj_hca)[grep("Rp", rownames(dg_obj_hca))])

FeaturePlot(dg_obj_hca, features = "Cluster1")


options(repr.plot.width=5, repr.plot.height=5)

DimPlot(dg_obj_hcp, group.by = "orig.ident",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


DimPlot(dg_obj, group.by = "age",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


DimPlot(dg_obj, group.by = "region",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=20, repr.plot.height=5)

DimPlot(dg_obj, split.by = "region",group.by = "orig.ident",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=10, repr.plot.height=5)

FeaturePlot(dg_obj, "Xist",split.by = "age",reduction = "umap", label = TRUE) +
  ggtitle("UMAP of DG Glut Subsampled Cells")


options(repr.plot.width=10, repr.plot.height=5)

FeaturePlot(dg_obj, "Nrg3os",split.by = "age",reduction = "umap", label = TRUE)


options(repr.plot.width=10, repr.plot.height=5)

VlnPlot(obj_integrated, "Malat1",group.by = "sample")
VlnPlot(dg_obj_hcp, "Malat1",group.by = "sample")


options(repr.plot.width=10, repr.plot.height=5)

VlnPlot(obj_integrated, "Rpl13",group.by = "sample")
VlnPlot(dg_obj_hcp, "Rpl13",group.by = "sample")


options(repr.plot.width=10, repr.plot.height=5)

VlnPlot(obj_integrated, "Meg3",group.by = "sample")
VlnPlot(dg_obj_hcp, "Meg3",group.by = "sample")


options(repr.plot.width=10, repr.plot.height=5)

VlnPlot(obj_integrated, "Robo1",group.by = "sample")
VlnPlot(dg_obj_hcp, "Robo1",group.by = "sample")



