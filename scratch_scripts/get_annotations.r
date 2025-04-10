library(Seurat)
library(ggplot2)
library(dplyr)
setwd("~/ps-renlab2/projects/combined_all/female_RNA/")

obj = readRDS("postfilter4.RDS")

meta=read.table("subclustering/meta.txt")

obj@meta.data = meta

obj = subset(obj, subset = final_clusters != "doub")

obj$filt_doub = "no"
obj$filt_doub[which(obj$sub_leiden%in%c("1 3","13 5", "2 3", "21 1", "5 7", "8 6", "9 4", "17 4", "2 6", "5 9", "13 5","2 5"))]="doub"
table(obj$filt_doub)
DimPlot(obj,group.by='filt_doub', label = T)
#obj = subset(obj, subset = filt_doub == "no")


obj

obj = subset(obj, subset = filt_doub == "no")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 5000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, features = VariableFeatures(obj))
obj <- FindNeighbors(obj, dims = 1:45)  # You can adjust the number of dimensions
obj <- FindClusters(obj, resolution = 1)  # Adjust the resolution as needed
obj <- RunUMAP(obj, reduction = "pca", dims = 1:45)  # Adjust the number of dimensions as needed
saveRDS(obj, "postfilter4.RDS")

obj$clade = "Neuron"
obj$clade[which(obj$RNA_snn_res.0.1 %in% c(0,7,10,15,12,18))] = "Glial"
neuron = subset(obj, subset = clade == "Neuron")
neuron <- NormalizeData(neuron)
neuron <- FindVariableFeatures(neuron, nfeatures = 5000)
neuron <- ScaleData(neuron)
neuron <- RunPCA(neuron, features = VariableFeatures(neuron))
neuron <- FindNeighbors(neuron, dims = 1:40)  # You can adjust the number of dimensions
neuron <- FindClusters(neuron, resolution = 1)  # Adjust the resolution as needed
neuron <- RunUMAP(neuron, reduction = "pca", dims = 1:40)  # Adjust the number of dimensions as needed
saveRDS(neuron , file = "neuron_RNA.RDS")


neuron$filt_doub = "no"
neuron$filt_doub[which(neuron$sub_leiden%in%c("13 5", "2 3", "21 1", "5 7", "8 6", "9 4", "17 4", "2 6", "5 9", "13 5","2 5"))]="doub"
table(neuron$filt_doub)

neuron = subset(neuron, subset = filt_doub == "no")

neuron <- NormalizeData(neuron)
neuron <- FindVariableFeatures(neuron, nfeatures = 5000)
neuron <- ScaleData(neuron)
neuron <- RunPCA(neuron, features = VariableFeatures(neuron))
neuron <- FindNeighbors(neuron, dims = 1:40)  # You can adjust the number of dimensions
neuron <- FindClusters(neuron, resolution = 1)  # Adjust the resolution as needed
neuron <- RunUMAP(neuron, reduction = "pca", dims = 1:40)  # Adjust the number of dimensions as needed
saveRDS(neuron , file = "neuron_RNA.RDS")

DimPlot(obj,group.by='final_clusters', label = T)

options(repr.plot.width=15, repr.plot.height=7)

DimPlot(neuron,group.by='final_clusters', label = T)

DimPlot(neuron,group.by='seurat_clusters', label = T)


neuron$filt_doub = "no"
neuron$filt_doub[which(neuron$sub_leiden%in%c("13 5", "2 3", "21 1", "4 4","5 7", "8 6", "9 4", "17 4", "2 6", "5 9", "2 5"))]="doub"

table(neuron$filt_doub)

options(repr.plot.width=25, repr.plot.height=5)

VlnPlot(neuron, group.by = "sub_leiden", "Combined_pANN")


fm = FindMarkers(neuron, `ident.1` = 53, `ident.2` = 1, group.by = "seurat_clusters")

FeaturePlot(neuron,c("Erbb4"))
FeaturePlot(neuron,c("Cntn2"))
FeaturePlot(neuron,c("Mal"))




options(repr.plot.width=5, repr.plot.height=5)

DimPlot(neuron,group.by='filt_doub', label = T)

options(repr.plot.width=25, repr.plot.height=5)

VlnPlot(neuron, group.by = "sub_leiden", "Mog",pt.size = 0)
VlnPlot(neuron, group.by = "sub_leiden", "Mal",pt.size = 0)

VlnPlot(neuron, group.by = "sub_leiden", "C1qb",pt.size = 0)
VlnPlot(neuron, group.by = "sub_leiden", "Apoe",pt.size = 0)




DimPlot(neuron,group.by='sub_leiden', label = T)

DimPlot(neuron,group.by='best_celltype_fixed', label = T)

DimPlot(neuron,group.by='final_clusters', label = T)

DimPlot(neuron,group.by='age', label = T)

FeaturePlot(neuron,c("Apoe","Slc1a2"),ncol = 3)

FeaturePlot(neuron,"Combined_pANN")

neu_genes = c("Snap25", "Gad1", "Gad2", "Slc32a1", "Slc17a6", "Slc17a7", "Slc6a5", "Slc6a4", "Slc6a3", "Slc18a3", "Slc6a2")

FeaturePlot(neuron,neu_genes, ncol = 4)

FeaturePlot(glia,neu_genes, ncol = 4)

FeaturePlot(neuron,c("Mog", "Mal", "Ppp1r14a", "Snap25"), ncol = 2)

FeaturePlot(neuron,c("Mog", "Mal", "Ppp1r14a", "Snap25"), ncol = 2)

FeaturePlot(glia,c("Mog", "Mal", "Ppp1r14a", "Snap25"), ncol = 2)

neuron

DimPlot(glia,group.by='final_clusters', label = T)

DimPlot(glia,group.by='age', label = T)

DimPlot(glia,group.by='region', label = T)

atac = read.table("../meta_coclust.txt")
atac$barcode = gsub("Female:", "", atac$cell_id) 
atac$barcode = gsub("-1", "", atac$barcode)
key = atac[which(!is.na(atac$rachel_celltypes)),]
key = key[which(!duplicated(key$rachel_celltypes)),]
rownames(key)= key$rachel_celltypes
metaf = atac
metaf = metaf[which(!is.na(metaf$rachel_celltypes)),]

predictions <- table(metaf$leiden_subcluster,metaf$rachel_celltypes)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

new_df <- predictions %>%
group_by(Var1) %>%
filter(Freq == max(Freq)) %>%
select(Var1, Var2, Freq)

new_df = as.data.frame(new_df)
new_df = new_df[which(new_df$Freq>0.5),]
rownames(new_df) = new_df$Var1

atac$rachel_celltypes_ext = new_df[paste(atac$`leiden_subcluster`), "Var2"]
atac$rachel_celltypes_ext = as.character(atac$rachel_celltypes_ext)


mat = match( obj@meta.data$barcode, atac$barcode, nomatch = NA)
obj@meta.data$co_cluster_l1_ext = atac[mat, "co_cluster_l1_ext"]
obj@meta.data$co_cluster_l2_ext = atac[mat, "co_cluster_l2_ext"]
obj@meta.data$final_clusters_RNA = atac[mat, "final_clusters_RNA"]
obj@meta.data$atac_leiden_subcluster = atac[mat, "leiden_subcluster"]
obj@meta.data$rachel_celltypes_ext = atac[mat, "rachel_celltypes_ext"]


head(obj@meta.data)



DimPlot(obj,group.by='final_clusters', label = T)

DimPlot(obj,group.by='RNA_snn_res.0.1', label = T)

obj$clade = "Neuron"
obj$clade[which(obj$RNA_snn_res.0.15 %in% c(0,7,10,15,12,18))] = "Glial"
obj$clade[which(obj$RNA_snn_res.0.15 %in% c(0,9,11,13))] = "Glial"

DimPlot(obj,group.by='clade', label = T)

obj = FindClusters(obj, resolution = 0.15)

DimPlot(obj,group.by='RNA_snn_res.0.15', label = T)

obj$is_1_3 = "NO"
obj$is_1_3[which(obj$sub_leiden == "1 3" )] = "1 3"
DimPlot(obj,group.by='is_1_3', label = T)

table(obj$is_1_3, obj$age)

options(repr.plot.width=12, repr.plot.height=7)

DimPlot(obj,group.by='sub_leiden', label = T)

FeaturePlot(obj,'Mal')

#allmeta = list()
#for(cl in unique(obj$RNA_snn_res.0.15)) {
for(cl in 11) {
   # cl = "Astro"
    sub = subset(obj, subset = RNA_snn_res.0.15%in%c(cl))
    sub <- NormalizeData(sub)
    sub <- FindVariableFeatures(sub, nfeatures = 2000)
    sub <- ScaleData(sub)
    sub <- RunPCA(sub, features = VariableFeatures(sub))
    sub <- FindNeighbors(sub, dims = 1:25)  # You can adjust the number of dimensions
    sub <- FindClusters(sub, resolution = 0.5)  # Adjust the resolution as needed
    sub <- RunUMAP(sub, reduction = "pca", dims = 1:25)  # Adjust the number of dimensions as needed
    sub
    sub@meta.data$sub_leiden = paste(cl, sub$seurat_clusters)
    
    metaf = sub@meta.data
    metaf = metaf[which(!is.na(metaf$best_celltype_fixed)),]

    predictions <- table(metaf$seurat_clusters,metaf$best_celltype_fixed)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    rownames(new_df) = new_df$Var1

    sub@meta.data$celltype_from_ATAC = new_df[paste(sub@meta.data$`seurat_clusters`), "Var2"]
    sub@meta.data$celltype_from_ATAC = as.character(sub@meta.data$celltype_from_ATAC)

    metaf = sub@meta.data
    metaf = metaf[which(!is.na(sub$predicted.id)),]
    metaf = metaf[which(metaf$prediction.score.max>0.75),]

    predictions <- table(metaf$seurat_clusters,metaf$predicted.id)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    #new_df = new_df[which(new_df$Freq>0.8),]
    #rownames(new_df) = new_df$Var1

    mat = match(sub@meta.data$seurat_clusters, new_df$Var1, nomatch = NA)
    sub@meta.data$predicted_id_allen = new_df[mat, "Var2"]
    sub@meta.data$predicted_id_allen = as.character(sub@meta.data$predicted_id_allen)
    sub@meta.data$predicted_id_allen = new_df[mat, "Var2"]
    sub@meta.data$predicted_id_allen[which(is.na(sub@meta.data$predicted_id_allen))] = paste(cl, "unannot", sep = "_")

   # pdf(paste(cl,"_sub.pdf", sep = ""))
    print(DimPlot(sub, group.by = "seurat_clusters", label = T))
    print(DimPlot(sub, group.by = "celltype_from_ATAC", label = T))
    print(DimPlot(sub, group.by = "predicted_id_allen", label = T))
    print(DimPlot(sub, group.by = "region", label = T))
    print(DimPlot(sub, group.by = "age"))
    print(FeaturePlot(sub, "nCount_RNA", max.cutoff = 2000))
    print(DimPlot(sub, group.by = "best_celltype_fixed"))
    print(DimPlot(sub, group.by = "predicted.id"))
    print(DimPlot(sub, group.by = "RNA_snn_res.3"))
   # dev.off()

    meta = sub@meta.data[,c('orig.ident','nCount_RNA','nFeature_RNA','percent.mt',
                                 'percent.ribo','barcode', 'region', 'age', 'rep',
                                 'predicted.id','prediction.score.max','celltype_from_ATAC',
                                'predicted_id_allen','sub_leiden')]
 #   write.table(meta, paste(cl,"sub_meta.txt", sep = "_") , sep = "\t" )
}
   # allmeta[[cl]] = meta
#}

new_df

  metaf = sub@meta.data
    metaf = metaf[which(!is.na(sub$predicted.id)),]
    metaf = metaf[which(metaf$prediction.score.max>0.8),]
    metaf = metaf[which(metaf$nCount_RNA>1000),]

    predictions <- table(metaf$seurat_clusters,metaf$predicted.id)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.8),]
    #rownames(new_df) = new_df$Var1

    mat = match(sub@meta.data$seurat_clusters, new_df$Var1, nomatch = NA)
    sub@meta.data$predicted_id_allen = new_df[mat, "Var2"]
    sub@meta.data$predicted_id_allen = as.character(sub@meta.data$predicted_id_allen)
    sub@meta.data$predicted_id_allen = new_df[mat, "Var2"]
    sub@meta.data$predicted_id_allen = as.character(sub@meta.data$predicted_id_allen)
    sub@meta.data$predicted_id_allen[which(is.na(sub@meta.data$predicted_id_allen))] = paste(cl, "unannot", sep = "_")


print(DimPlot(sub, group.by = "predicted_id_allen", label = T))


print(DimPlot(sub, group.by = "predicted_id_allen", label = T))


sub$prediction.score.max

print(FeaturePlot(sub[,which(!is.na(sub$predicted.id))], "prediction.score.max"))

options(repr.plot.width=16, repr.plot.height=12)

print(DimPlot(sub[,which(!is.na(sub$predicted.id))], group.by = "predicted.id", label = T))


print(DimPlot(sub, group.by = "predicted_id_allen", label = T))



print(DimPlot(sub, group.by = "rachel_celltypes_ext", label = T))


options(repr.plot.width=10, repr.plot.height=7)

print(DimPlot(sub, group.by = "predicted_id_allen", label = T))


options(repr.plot.width=15, repr.plot.height=7)

print(DimPlot(sub, group.by = "predicted_id_allen", label = T))


options(repr.plot.width=35, repr.plot.height=7)

print(DimPlot(sub, group.by = "predicted.id", label = T))


new_df[which(new_df$Freq<0.7),]

print(FeaturePlot(sub, "Combined_pANN"))



