library(Seurat)


obj = readRDS("RNA_final_SCT.RDS")


obj
rna = obj@meta.data


rna=read.csv("meta_final.csv")

rna$match_barcode = paste("Female:",rna$barcode,"-1", sep= "")
head(rna$match_barcode )

meta = read.csv("../final_meta.csv")

head(meta)

mat = match( rna$match_barcode,meta$cell_id)

length(mat)

rna$tsse = meta[mat, "tsse"]
rna$n_fragment = meta[mat, "n_fragment"]

rna = rna[-which(is.na(rna$tsse)),]
nrow(rna)

quantile(rna$nCount_RNA,.99)
quantile(rna$n_fragment,.99)
quantile(rna$percent.mt,.99)
quantile(rna$percent.ribo,.99)
quantile(rna$nCount_RNA,c(.05,.99))



#rna =rna[which(rna$tsse>7&rna$nCount_RNA>1000&rna$nCount_RNA<30000&rna$n_fragment>1000&rna$n_fragment<70000&rna$percent.mt<.15&rna$percent.ribo<1),]
#rna =rna[which(rna$tsse>5.5&rna$nCount_RNA<35000&rna$n_fragment<75000&rna$percent.mt<.05&rna$percent.ribo<.5),]
rna=rna[which(rna$tsse>6&rna$percent.mt<.15&rna$percent.ribo<1&rna$nCount_RNA<35000&rna$n_fragment<75000),]
nrow(rna[which(rna$tsse>6&rna$percent.mt<.15&rna$percent.ribo<1&rna$nCount_RNA<35000&rna$n_fragment<75000),])

rna$celltype_age_region_rep = paste(rna$celltype_final,rna$age, rna$region,rna$rep)

max_cells_per_cluster <- 1000

# Create a score based on nCount_RNA and n_fragment
rna$score <- rna$nCount_RNA + rna$n_fragment

# Create a vector to store the cell indices you want to keep
indices_to_keep <- integer(0)

# Iterate through unique clusters and get the indices
unique_clusters <- unique(rna$celltype_age_region_rep)
for (cluster in unique_clusters) {
  # Get the indices of cells in the current cluster
  cluster_indices <- which(rna$celltype_age_region_rep == cluster)
  
  # If the cluster has more cells than the limit, select the top cells based on the score
  if (length(cluster_indices) > max_cells_per_cluster) {
    sorted_indices <- order(rna$score[cluster_indices], decreasing = TRUE)
    selected_indices <- cluster_indices[sorted_indices[1:max_cells_per_cluster]]
    indices_to_keep <- c(indices_to_keep, selected_indices)
  } else {
    indices_to_keep <- c(indices_to_keep, cluster_indices)
  }
}

# Subsetting the table based on the selected indices
selected_cells <- rna[indices_to_keep, ]


table(selected_cells$celltype_age_region_rep)

max_cells_per_cluster <- 500

# Create a vector to store the cell indices you want to keep
indices_to_keep <- integer(0)

# Iterate through unique clusters and get the indices
unique_clusters <- unique(rna$celltype_age_region_rep)
for (cluster in unique_clusters) {
  # Get the indices of cells in the current cluster
  cluster_indices <- which(rna$celltype_age_region_rep == cluster)
  
  # If the cluster has more cells than the limit, sample the indices
  if (length(cluster_indices) > max_cells_per_cluster) {
    sampled_indices <- sample(cluster_indices, max_cells_per_cluster)
    indices_to_keep <- c(indices_to_keep, sampled_indices)
  } else {
    indices_to_keep <- c(indices_to_keep, cluster_indices)
  }
}


max_cells_per_cluster <- 2000

# Create a vector to store the cell indices you want to keep
indices_to_keep <- integer(0)

# Iterate through unique clusters and get the indices
unique_clusters <- unique(rna$RNA_snn_res.3)
for (cluster in unique_clusters) {
  # Get the indices of cells in the current cluster
  cluster_indices <- which(rna$RNA_snn_res.3 == cluster)
  
  # If the cluster has more cells than the limit, sample the indices
  if (length(cluster_indices) > max_cells_per_cluster) {
    sampled_indices <- sample(cluster_indices, max_cells_per_cluster)
    indices_to_keep <- c(indices_to_keep, sampled_indices)
  } else {
    indices_to_keep <- c(indices_to_keep, cluster_indices)
  }
}


table(rna$celltype_final)
rna=rna[indices_to_keep,]
table(rna$celltype_final)

nrow(rna)

filtered_seurat = subset(obj, cells = which(obj$barcode%in% rna$barcode))


filtered_seurat

DimPlot(filtered_seurat)

head(rna)

write.csv(rna_ord[,c("match_barcode", "orig.ident","celltype_final", "nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo" ,"tsse", "n_fragment" )], file = "../subsampled_high_quality_ATAC_RNA_meta.csv")

all(filtered_seurat$barcode == rna_ord$barcode)

filtered_seurat$cell_id = rna_ord$match_barcode

head(filtered_seurat@meta.data$barcode)

head(rna$barcode)

rna_ord= rna[paste(rownames(filtered_seurat@meta.data)),]

head(rna_ord)

filtered_seurat

saveRDS(filtered_seurat, "subsampled_high_quality_RNA_seurat.RDS")

colnames(filtered_seurat) = filtered_seurat$cell_id

head(filtered_seurat@meta.data
    )

rm(filtered_seurat)

rm(obj)


