library(Seurat)
library(ggplot2)
library(dplyr)

al = read.table("../transfer_allen/AIT21_annotation_freeze_081523.tsv", header = T)



al = al[!duplicated(al$subclass_label),]
rownames(al) = al$subclass_label



al = al[which(al$subclass_label%in% unique(obj$final_clusters)),]

best_genes = c("Snap25",  "Gad1" ,   "Gad2"  ,  "Slc32a1" ,"Slc17a6" ,"Slc17a7", "Slc17a8",
           "Slc6a5" , "Slc6a4" , "Slc6a3",  "Slc18a3","Hdc", "Slc6a2"  )
length(best_genes

annot = al

#best_genes = c("Snap25",  "Gad1" ,   "Gad2"  ,  "Slc32a1" ,"Slc17a6" ,"Slc17a7", "Slc17a8",
#          "Slc6a5" , "Slc6a4" , "Slc6a3",  "Slc18a3","Hdc", "Slc6a2"  )
best_genes= c()

classes = al$subclass_label
for(i in classes){
    curr = annot[which(annot$subclass_label == i),]
    genes = paste(curr$cluster.markers, collapse = ",")
    genes = strsplit(genes, ",")[[1]]
    genes = genes[which(!(genes%in%best_genes))]
    best_gene = names(table(genes))[which(table(genes) == max(table(genes)))[1]]
    #best_gene = genes[1]
    best_genes = c(best_genes, best_gene)
}

best_genes

best_tab = cbind(classes,best_genes)

alcl = c(classes,unique(obj$final_clusters)[which(!unique(obj) %in% classes)])

alcl

obj = readRDS("postfilter3.RDS")

meta=read.table("subclustering/meta.txt")

nrow(meta)

nrow(obj@meta.data)

obj@meta.data = meta

obj = subset(obj, subset = final_clusters != "doub")

Idents(obj)="final_clusters"

levels(Idents(obj))

colnames(obj@meta.data)

atac = read.table("../meta_coclust.txt")
atac$barcode = gsub("Female:", "", atac$cell_id) 
atac$barcode = gsub("-1", "", atac$barcode) 

length(which(is.na(atac$leiden_subcluster)))

head(atac$rachel_celltypes)

key = atac[which(!is.na(atac$rachel_celltypes)),]
key = key[which(!duplicated(key$rachel_celltypes)),]
rownames(key)= key$rachel_celltypes

length(unique(key$co_cluster_l1))

meta_rach = read.csv("/home/lamaral/projects/brain_aging_mouse/analysis/integrtaion_adata/final_metadata.csv")

meta_rach$co2_ct = paste(meta_rach$co_cluster_l2 , meta_rach$new_celltypes)

key = meta_rach[which(!duplicated(meta_rach$co2_ct)),]
head(key)

length(unique(key$new_celltypes))

length(unique(key$co_cluster_l1))

atac$rachel_celltypes[1:5]
    metaf = atac


library(dplyr)
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


 library(dplyr)
    metaf = obj@meta.data
    metaf = metaf[which(!is.na(metaf$rachel_celltypes_ext)),]

    predictions <- table(metaf$sub_leiden,metaf$rachel_celltypes_ext)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.5),]
    rownames(new_df) = new_df$Var1

    obj@meta.data$rachel_celltypes_ext_ext = new_df[paste(obj@meta.data$`sub_leiden`), "Var2"]
    obj@meta.data$rachel_celltypes_ext_ext = as.character(obj@meta.data$rachel_celltypes_ext_ext)


mat = match( obj@meta.data$barcode, atac$barcode, nomatch = NA)
obj@meta.data$co_cluster_l1_ext = atac[mat, "co_cluster_l1_ext"]
obj@meta.data$co_cluster_l2_ext = atac[mat, "co_cluster_l2_ext"]
obj@meta.data$final_clusters_RNA = atac[mat, "final_clusters_RNA"]
obj@meta.data$atac_leiden_subcluster = atac[mat, "leiden_subcluster"]
obj@meta.data$rachel_celltypes_ext = atac[mat, "rachel_celltypes_ext"]


library(dplyr)
    metaf = obj@meta.data
    metaf = metaf[which(!is.na(metaf$)),]

    predictions <- table(metaf$sub_leiden,metaf$co_cluster_l2_ext)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.5),]
    rownames(new_df) = new_df$Var1

    obj@meta.data$co_cluster_l2_ext_ext = new_df[paste(obj@meta.data$`sub_leiden`), "Var2"]
    obj@meta.data$co_cluster_l2_ext_ext = as.character(obj@meta.data$co_cluster_l2_ext_ext)


library(dplyr)
    metaf = obj@meta.data
    metaf = metaf[which(!is.na(metaf$)),]

    predictions <- table(metaf$sub_leiden,metaf$co_cluster_l2_ext)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    new_df = new_df[which(new_df$Freq>0.5),]
    rownames(new_df) = new_df$Var1

    obj@meta.data$co_cluster_l2_ext_ext = new_df[paste(obj@meta.data$`sub_leiden`), "Var2"]
    obj@meta.data$co_cluster_l2_ext_ext = as.character(obj@meta.data$co_cluster_l2_ext_ext)


length(unique(obj$predicted.id))

options(repr.plot.width=42, repr.plot.height=8)

DimPlot(obj[,which(!is.na(obj$predicted.id))], group.by = "predicted.id", label = F, raster=FALSE)


options(repr.plot.width=12, repr.plot.height=8)

DimPlot(obj, group.by = "rachel_celltypes_ext_ext", label = T, raster=FALSE)


options(repr.plot.width=12, repr.plot.height=8)

DimPlot(obj, group.by = "rachel_celltypes_ext", label = T, raster=FALSE)


options(repr.plot.width=16, repr.plot.height=8)

DimPlot(obj, group.by = "co_cluster_l2_ext", label = T, raster=FALSE)


DimPlot(obj, group.by = "co_cluster_l2_ext_ext", label = T, raster=FALSE)


options(repr.plot.width=12, repr.plot.height=8)

DimPlot(obj, group.by = "co_cluster_l1_ext", label = T, raster=FALSE)


obj$final_clusters = factor(obj$final_clusters, levels = alcl)

options(repr.plot.width=30, repr.plot.height=30)

DotPlot(obj, group.by = "final_clusters",features=best_genes,scale = T,col.min = -2,
  col.max = 2)

obj = FindClusters(obj, resolution = 0.15)

DimPlot(obj,group.by='RNA_snn_res.0.1', label = T)

options(repr.plot.width=20, repr.plot.height=8)

DimPlot(obj,group.by='best_celltype_from_ATAC', label = T)

sub = subset(obj, subset = RNA_snn_res.0.1==10)

sub

setwd("subclustering")

#allmeta = list()
#for(cl in unique(obj$RNA_snn_res.0.1)) {
for(cl in 1) {
   # cl = "Astro"
    sub = subset(obj, subset = RNA_snn_res.0.1%in%c(cl))
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
    metaf = metaf[which(metaf$prediction.score.max>0.8),]

    predictions <- table(metaf$seurat_clusters,metaf$predicted.id)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    #rownames(new_df) = new_df$Var1

    mat = match(sub@meta.data$seurat_clusters, new_df$Var1, nomatch = NA)
    sub@meta.data$predicted_id_allen = new_df[mat, "Var2"]
    sub@meta.data$predicted_id_allen = as.character(sub@meta.data$predicted_id_allen)
    sub@meta.data$predicted_id_allen = new_df[mat, "Var2"]


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

    print(DimPlot(sub, group.by = "seurat_clusters", label = T))


options(repr.plot.width=18, repr.plot.height=8)

print(DimPlot(sub, group.by = "predicted_id_allen", label = T))


print(DimPlot(sub[,which(!is.na(sub$predicted.id))], group.by = "predicted.id"))


length(unique(sub$predicted.id))

 sub <- FindClusters(sub, resolution = 0.2)  # Adjust the resolution as needed
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
    metaf = metaf[which(metaf$prediction.score.max>0.8),]

    predictions <- table(metaf$seurat_clusters,metaf$predicted.id)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    #rownames(new_df) = new_df$Var1

    mat = match(sub@meta.data$seurat_clusters, new_df$Var1, nomatch = NA)
    sub@meta.data$predicted_id_allen = new_df[mat, "Var2"]
    sub@meta.data$predicted_id_allen = as.character(sub@meta.data$predicted_id_allen)
    sub@meta.data$predicted_id_allen = new_df[mat, "Var2"]


    pdf(paste(cl,"_sub.pdf", sep = ""))
    print(DimPlot(sub, group.by = "seurat_clusters", label = T))
    print(DimPlot(sub, group.by = "celltype_from_ATAC", label = T))
    print(DimPlot(sub, group.by = "predicted_id_allen", label = T))
    print(DimPlot(sub, group.by = "region", label = T))
    print(DimPlot(sub, group.by = "age"))
    print(FeaturePlot(sub, "nCount_RNA", max.cutoff = 2000))
    print(DimPlot(sub, group.by = "best_celltype_fixed"))
    print(DimPlot(sub, group.by = "predicted.id"))
    print(DimPlot(sub, group.by = "RNA_snn_res.3"))
    dev.off()

    meta = sub@meta.data[,c('orig.ident','nCount_RNA','nFeature_RNA','percent.mt',
                                 'percent.ribo','barcode', 'region', 'age', 'rep',
                                 'predicted.id','prediction.score.max','celltype_from_ATAC',
                                'predicted_id_allen','sub_leiden')]
    write.table(meta, paste(cl,"sub_meta.txt", sep = "_") , sep = "\t" )


table(meta$celltype_from_ATAC)

options(repr.plot.width=30, repr.plot.height=8)
print(DimPlot(sub, group.by = "predicted.id"))


print(FeaturePlot(sub, "Mog"))

dm = sub@meta.data[which(sub@meta.data$seurat_clusters %in% c(4,14,3,12)),]
nrow(dm)

table(dm$age)

mat = match(dm$barcode, obj@meta.data$barcode)
length(mat)
head(mat)

obj@meta.data[mat, "final_clusters"] = "doub"

options(repr.plot.width=40, repr.plot.height=8)
print(DimPlot(sub, group.by = "predicted.id"))

print(DimPlot(sub, group.by = "predicted.id"))

all = do.call(rbind,allmeta)

all$final_celltype = all$celltype_from_ATAC


use_allen = c("11 4", "10 5", "16 1", "18 0", "18 2", "1 6", "1 7","1 5", "1 11", 
              "16 1", "6 2", "6 6", "6 7", "6 9", "9 4")
for (us in use_allen){
    all[which(sub_leiden==us), "final_celltype"] = all[which(sub_leiden==us), "predicted_id_allen"]
    }

poss_doub = c("10 2", "10 6", "1 3")


poss_doub = c("6 9", "6 6", "6 7","2 7")
for (us in poss_doub){
    obj@meta.data$final_clusters[which(obj@meta.data$sub_leiden==us)] = "doub"
}


write.table(all, "all_meta.txt", sep = "\t") 





all =read.table("all_meta.txt")

head(all)

mat = match(sub$barcode, all$barcode )
head(mat)

table(sub$celltype_from_ATAC)

all[mat, "final_celltype"] = sub$celltype_from_ATAC

table(all$final_celltype)

mat = match( obj@meta.data$barcode,all$barcode, nomatch = NA)
length(mat)



obj@meta.data$sub_leiden = all[mat, "sub_leiden"]

obj@meta.data$final_clusters = all[mat, "final_celltype"]

write.table(obj@meta.data, "meta.txt", sep = "\t")

atac_meta = read.table("../../meta_subcluster_use.txt", sep = " ", header = T)

    head(atac_meta$cell_id[which(atac_meta$batch == "Female")])

length(which(unique(obj@meta.data$cell_id) %in% atac_meta$cell_id))

head(obj@meta.data$cell_id)

head(atac_meta)

    metaf = atac_meta
    metaf = metaf[which(!is.na(metaf$final_clusters_RNA)),]

    predictions <- table(metaf$leiden_subcluster,metaf$final_clusters_RNA)
    predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
    predictions <- as.data.frame(predictions)

    new_df <- predictions %>%
    group_by(Var1) %>%
    filter(Freq == max(Freq)) %>%
    select(Var1, Var2, Freq)

    new_df = as.data.frame(new_df)
    rownames(new_df) = new_df$Var1

    atac_meta$final_clusters_RNA_ext = new_df[paste(atac_meta$`leiden_subcluster`), "Var2"]
    atac_meta$final_clusters_RNA_ext = as.character(atac_meta$final_clusters_RNA_ext)


write.table(atac_meta, file = "../../meta_from_RNA.txt", sep = "\t")

length(unique(atac_meta$final_clusters_RNA_ext))



obj@meta.data$cell_id = paste("Female:", obj@meta.data$barcode,"-1", sep = "")

mat = match(obj@meta.data$cell_id , atac_meta$cell_id, nomatch = NA)

length(mat)

options(repr.plot.width=6, repr.plot.height=6)

DimPlot(obj, group.by = "RNA_snn_res.0.1", label = T, raster=FALSE)


DimPlot(obj, group.by = "RNA_snn_res.0.1", label = T, raster=FALSE)


options(repr.plot.width=36, repr.plot.height=36)

DimPlot(obj, split.by = "RNA_snn_res.0.1",group.by = "final_clusters" , label = T, raster=FALSE,ncol = 5)


options(repr.plot.width=22, repr.plot.height=12)

DimPlot(obj, group.by = "sub_leiden", label = T, raster=FALSE)


hist(obj$prediction.score.max)

options(repr.plot.width=12, repr.plot.height=12)

FeaturePlot(obj, "prediction.score.max", raster=FALSE)


length(unique(obj$final_clusters))

options(repr.plot.width=22, repr.plot.height=12)

DimPlot(obj, group.by = "final_clusters", label = T, raster=FALSE)


options(repr.plot.width=22, repr.plot.height=12)

DimPlot(obj, group.by = "final_clusters", label = T, raster=FALSE)




atac = read.table("../meta_subcluster_use.txt", sep = " ", header = T)
atac = atac[which(atac$batch == "Female"),]
atac$barcode = gsub("Female:", "", atac$cell_id) 
atac$barcode = gsub("-1", "", atac$barcode) 

mat = match( obj@meta.data$barcode, atac$barcode, nomatch = NA)
obj@meta.data$best_celltype_fixed = atac[mat, "best_celltype_fixed"]
obj@meta.data$leiden_subcluster_atac = atac[mat, "leiden_subcluster"]

metaf = obj@meta.data
metaf = metaf[which(!is.na(metaf$best_celltype_fixed)),]

predictions <- table(metaf$RNA_snn_res.3,metaf$best_celltype_fixed)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

new_df <- predictions %>%
  group_by(Var1) %>%
  filter(Freq == max(Freq)) %>%
  select(Var1, Var2, Freq)

new_df = as.data.frame(new_df)
rownames(new_df) = new_df$Var1

obj@meta.data$best_celltype_fixed_ext = new_df[paste(obj@meta.data$RNA_snn_res.3), "Var2"]
obj@meta.data$best_celltype_fixed_ext = as.character(obj@meta.data$best_celltype_fixed_ext)

allen = read.table("transfer_allen/meta_filt_pred.rpcaref_nov14.txt")
mat = match( rownames(obj@meta.data), rownames(allen), nomatch = NA)
obj@meta.data$predicted.id = allen[mat, "predicted.id"]
obj@meta.data$predicted_class = allen[mat, "predicted_class"]
obj@meta.data$prediction.score.max = allen[mat, "prediction.score.max"]

metaf = obj@meta.data
metaf = metaf[which(!is.na(metaf$predicted.id)),]
metaf = metaf[which(metaf$prediction.score.max>0.8),]

predictions <- table(metaf$RNA_snn_res.3,metaf$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

new_df <- predictions %>%
  group_by(Var1) %>%
  filter(Freq == max(Freq)) %>%
  select(Var1, Var2, Freq)

new_df = as.data.frame(new_df)
rownames(new_df) = new_df$Var1

obj@meta.data$predicted_id_all = new_df[paste(obj@meta.data$RNA_snn_res.3), "Var2"]
obj@meta.data$predicted_id_all = as.character(obj@meta.data$predicted_id_all)


length(unique(obj$predicted_id_all))

options(repr.plot.width=22, repr.plot.height=12)

DimPlot(obj, group.by = "predicted_id_all", label = T, raster=FALSE)


DimPlot(obj, group.by = "best_celltype_fixed_ext", label = T, raster=FALSE)


DimPlot(obj, group.by = "best_celltype_fixed", label = T, raster=FALSE)


write.table(obj@meta.data, "ext_meta.txt", sep = "\t", row.names = F)

DimPlot(obj, group.by = "RNA_snn_res.3", label = T, raster=FALSE)


min(obj$nFeature_RNA)

table(obj$predicted_id_all[which(obj$filt_lowQ=="Yes")])

options(repr.plot.width=20, repr.plot.height=10)

DimPlot(filtered_seurat, group.by = "predicted_id_all", label = T, raster=FALSE)


options(repr.plot.width=7, repr.plot.height=7)

DimPlot(filtered_seurat, group.by = "age", label = T, raster=FALSE)
FeaturePlot(filtered_seurat, features = "nCount_RNA",  raster=FALSE, max.cutoff = 3000)


mat = match( filtered_seurat@meta.data$barcode, atac$barcode, nomatch = NA)
filtered_seurat@meta.data$best_celltype_fixed = atac[mat, "best_celltype_fixed"]
filtered_seurat@meta.data$leiden_subcluster_atac = atac[mat, "leiden_subcluster"]

options(repr.plot.width=22, repr.plot.height=7)

DimPlot(filtered_seurat, group.by = "best_celltype_fixed", label = T, raster=FALSE)


saveRDS(filtered_seurat, "postfilter3.RDS")

head(obj@meta.data$orig.ident)

obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))

options(repr.plot.width=15, repr.plot.height=15)

obj@meta.data %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+ facet_wrap(~region+age, ncol = 3)


min(obj@meta.data$nCount_RNA)

obj@meta.data %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 300)+ facet_wrap(~region+age, ncol = 3)


length(which(obj$predicted_id_all == "Lymphoid NN"))

lym  = obj[,which(obj$predicted_id_all == "Lymphoid NN")]
lym

metaf = obj@meta.data
metaf = metaf[which(!is.na(metaf$best_celltype_fixed)),]

predictions <- table(metaf$RNA_snn_res.3,metaf$best_celltype_fixed)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

new_df <- predictions %>%
  group_by(Var1) %>%
  filter(Freq == max(Freq)) %>%
  select(Var1, Var2, Freq)

new_df = as.data.frame(new_df)
rownames(new_df) = new_df$Var1

obj@meta.data$best_celltype_fixed_ext = new_df[paste(obj@meta.data$RNA_snn_res.3), "Var2"]
obj@meta.data$best_celltype_fixed_ext = as.character(obj@meta.data$best_celltype_fixed_ext)


new_df[which(new_df$Freq<0.5),"Var1"]
obj@meta.data$unpred = new_df[paste(obj@meta.data$RNA_snn_res.3), "Var2"]
obj@meta.data$unpred[which(obj@meta.data$`RNA_snn_res.3` %in% paste(new_df[which(new_df$Freq<0.5),"Var1"]))] = "No pred"

new_df[which(new_df$Var1 ==132),]

DimPlot(obj, group.by = "unpred", label = T, raster=FALSE)


options(repr.plot.width=17, repr.plot.height=8)

obj$"celltype_from_ATAC" = obj$best_celltype_fixed_ext
DimPlot(obj, group.by = "celltype_from_ATAC", label = T, raster=FALSE, )


options(repr.plot.width=24, repr.plot.height=12)

DimPlot(obj, group.by = "best_celltype_fixed", label = T, raster=FALSE)


options(repr.plot.width=24, repr.plot.height=12)

DimPlot(obj, group.by = "leiden_subcluster_atac", label = T, raster=FALSE)


DimPlot(obj, group.by = "RNA_snn_res.3", label = T, raster=FALSE)


obj@meta.data$present_in_ATAC=FALSE
obj@meta.data$present_in_ATAC[which(!is.na(obj@meta.data$leiden_subcluster_atac))]=TRUE
table(obj@meta.data$present_in_ATAC)

ggplot(obj@meta.data, aes( x = present_in_ATAC, y = nCount_RNA)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +ylim(c(0,3000))

options(repr.plot.width=12, repr.plot.height=12)

FeaturePlot(obj, "nCount_RNA",max.cutoff = 2000)


DimPlot(obj, group.by = "age", label = T, raster=TRUE)



