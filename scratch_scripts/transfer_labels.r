library(Seurat)
library(data.table)
library(sceasy)
library(reticulate)
use_condaenv("seurat")

sessionInfo()

sceasy::convertFormat("filtered_1000k_per_subclass.h5ad", from="anndata", to="seurat",outFile='filtered_1000k_per_subclass.rds')

allen <- readRDS('filtered_1000k_per_subclass.rds')

##Subsample Allen

# Set the maximum number of cells per cluster
max_cells_per_cluster <- 250

# Create a vector to store the cell indices you want to keep
indices_to_keep <- integer(0)

# Iterate through unique clusters and get the indices
unique_clusters <- unique(allen$cl)
for (cluster in unique_clusters) {
  # Get the indices of cells in the current cluster
  cluster_indices <- which(allen$cl == cluster)
  
  # If the cluster has more cells than the limit, sample the indices
  if (length(cluster_indices) > max_cells_per_cluster) {
    sampled_indices <- sample(cluster_indices, max_cells_per_cluster)
    indices_to_keep <- c(indices_to_keep, sampled_indices)
  } else {
    indices_to_keep <- c(indices_to_keep, cluster_indices)
  }
}

# Subset the Seurat object all at once using the collected indices
allen <- subset(allen, cells = indices_to_keep)

filtered_seurat = readRDS("../postfilter3.RDS")

filtered_seurat = subset(filtered_seurat, subset = nCount_RNA > 500)


###Subsample data
# Set the maximum number of cells per cluster
max_cells_per_cluster <- 2000

# Create a vector to store the cell indices you want to keep
indices_to_keep <- integer(0)

# Iterate through unique clusters and get the indices
unique_clusters <- unique(filtered_seurat$RNA_snn_res.3)
for (cluster in unique_clusters) {
  # Get the indices of cells in the current cluster
  cluster_indices <- which(filtered_seurat$RNA_snn_res.3 == cluster)
  
  # If the cluster has more cells than the limit, sample the indices
  if (length(cluster_indices) > max_cells_per_cluster) {
    sampled_indices <- sample(cluster_indices, max_cells_per_cluster)
    indices_to_keep <- c(indices_to_keep, sampled_indices)
  } else {
    indices_to_keep <- c(indices_to_keep, cluster_indices)
  }
}
filtered_seurat <- subset(filtered_seurat, cells = indices_to_keep)

allen <- FindVariableFeatures(allen, nfeatures = 5000)
allen <- ScaleData(allen)
allen <- RunPCA(allen, features = VariableFeatures(allen))

saveRDS(allen, file = "Female_10xV3_1000_per_cl_8k_gene_processed.rds")

m = fread("AIT21_annotation_freeze_081523.tsv")
m = as.data.frame(m)
rownames(m) = paste(m$cl)
mat = match(allen@meta.data$cl,m$cl)
#allen@meta.data$subclass_label = m[mat,"subclass_label"]
allen@meta.data$class_label = m[mat,"class_label"]

anchors <- FindTransferAnchors(reference = allen, query = filtered_seurat,
    dims = 1:40, reduction = "rpca", max.features = 250, k.anchor = 30, features = VariableFeatures(allen))

## An AnchorSet object containing 88172 anchors between the reference and query Seurat objects. 


predictions_cl <- TransferData(anchorset = anchors, refdata = allen$subclass_label,
    dims = 1:40, weight.reduction = "rpca.ref")
#filtered_seurat <- AddMetaData(filtered_seurat, metadata = predictions_cl)


filtered_seurat <- AddMetaData(filtered_seurat, metadata = predictions)


DimPlot(filtered_seurat, group.by = "predicted.id", label = T, raster=FALSE,legend = T)


predictions_class <- TransferData(anchorset = anchors, refdata = allen$class_label,
    dims = 1:40, weight.reduction = "rpca.ref")

anchors <- FindTransferAnchors(reference = al_glia, query = glia,
    dims = 1:30, reduction = "rpca", max.features = 250, k.anchor = 20, features = VariableFeatures(al_glia))

al_neu <- FindVariableFeatures(al_neu, nfeatures = 4000)
al_neu <- ScaleData(al_neu)
al_neu <- RunPCA(al_neu, features = VariableFeatures(al_neu))
al_neu <- RunUMAP(al_neu)

pdf("al_neu.pdf", height =6, width = 16)
DimPlot(al_neu, group.by = "subclass_label", label = T, raster=FALSE)
dev.off()

pdf("al_glia.pdf", height =6, width = 12)
DimPlot(al_glia, group.by = "subclass_label", label = T, raster=FALSE)
dev.off()

anchors <- FindTransferAnchors(reference = al_glia, query = glia,
    dims = 1:30, reduction = "cca", 
                               features = VariableFeatures(al_glia))

anchors <- FindTransferAnchors(reference = al_neu, query = neuron,
    dims = 1:40, reduction = "rpca", 
                               features = VariableFeatures(al_neu))



#rpca glia
#Finding anchors
#        Found 23636 anchors
#Filtering anchors
#        Retained 7996 anchors

#rpca neur
Finding anchors

        Found 40160 anchors
Filtering anchors
        Retained 13433 anchors


predictions_glcl <- TransferData(anchorset = anchors, refdata = al_glia$subclass_label,
    dims = 1:30, weight.reduction = "rpca.ref")


predictions_cca <- TransferData(anchorset = anchors, refdata = al_glia$subclass_label,
    dims = 1:30, weight.reduction = "cca")


glia <- AddMetaData(glia, metadata = predictions_glcl)


pdf("glia_pred.pdf", height = 6, width = 10)
DimPlot(glia, group.by = "predicted.id", label = T)
dev.off()

glia$cca_pred = predictions_cca$predicted.id
pdf("glia_cca.pdf", height = 6, width = 10)
DimPlot(glia, group.by = "cca_pred", label = T)
dev.off()



gaba = readRDS("../gaba.RDS")
anchors_glut = anchors
al_neu = allen[,which(!allen$subclass_label %in% unique(allen$subclass_label)[grep("NN", unique(allen$subclass_label))])]
al_neu <- FindVariableFeatures(al_neu, nfeatures = 5000)
al_neu <- ScaleData(al_neu)
al_neu <- RunPCA(al_neu, features = VariableFeatures(al_neu))
al_neu <- RunUMAP(al_neu,dims = 1:40)
pdf('al_neu_UMAP.pdf', height = 7, width = 35)
DimPlot(al_neu,group.by = "subclass_label", label = T)
dev.off()
anchors <- FindTransferAnchors(reference = al_neu, query = gaba,
    dims = 1:40, reduction = "rpca", max.features = 250, k.anchor = 25, features = VariableFeatures(al_neu))
predictions<- TransferData(anchorset = anchors, refdata = al_neu$subclass_label,
    dims = 1:40, weight.reduction = "rpca.ref")
gaba <- AddMetaData(gaba, metadata = predictions)
pdf('gaba_UMAP.pdf', height = 7, width = 30)
DimPlot(gaba,group.by = "predicted.id", label = T)
dev.off()
write.table(gaba@meta.data, file = "gaba_meta.txt", sep = "\t")
sink("gaba_anchors.txt")
cat(anchors)
sink()

glut = readRDS("../glut.RDS")
al_glut = allen[,which(allen$subclass_label %in% unique(allen$subclass_label)[grep("Glut", unique(allen$subclass_label))])]
al_glut <- FindVariableFeatures(al_glut, nfeatures = 5000)
al_glut <- ScaleData(al_glut)
al_glut <- RunPCA(al_glut, features = VariableFeatures(al_glut))
al_glut <- RunUMAP(al_glut,dims = 1:40)
pdf('al_glut_UMAP.pdf', height = 7, width = 35)
DimPlot(al_glut,group.by = "subclass_label", label = T)
dev.off()
anchors <- FindTransferAnchors(reference = al_glut, query = glut,
    dims = 1:40, reduction = "rpca", max.features = 250, k.anchor = 20, features = VariableFeatures(al_glut))
predictions<- TransferData(anchorset = anchors, refdata = al_glut$subclass_label,
    dims = 1:40, weight.reduction = "rpca.ref")
glut <- AddMetaData(glut, metadata = predictions)
pdf('glut_UMAP.pdf', height = 7, width = 30)
DimPlot(glut,group.by = "predicted.id", label = T)
dev.off()
write.table(glut@meta.data, file = "glut_meta.txt", sep = "\t")
sink("glut_anchors.txt")
cat(anchors)
sink()
