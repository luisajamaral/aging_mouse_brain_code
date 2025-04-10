library(Seurat)
library(ggplot2)
#combined_seurat = readRDS(file = "/home/lamaral/oasis/female_RNA/combined_seurat_UMAP.RDS")
#head(combined_seurat@meta.data,5)




# Read data
m0 <- read.csv("subclustering3/9_combined_meta.csv")

# Assuming 'meta' is your metadata object
mat <- match(m0$cell_id, meta$cell_id)
m0$celltype <- meta$celltype[mat]
m0$leiden_subcluster <- meta$leiden_subcluster[mat]
m0$CellType_1127_ext <- meta$CellType_1127_ext[mat]
m0$CellType_1127 <- meta$CellType_1127[mat]
m0$predicted.id = meta$predicted.id[mat]
m0$predicted_id_sep <- meta$predicted_id_sep[mat]
m0$celltypes_ext = meta$celltypes_ext[mat]
m0$present = "No"
m0$present[which(!is.na(m0$celltype))]="Yes" 

#ext leiden

metaf = m0
metaf = metaf[which(!is.na(m0$celltype)),]
predictions <- table(metaf$leiden,metaf$celltype)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

new_df <- predictions %>%
group_by(Var1) %>%
filter(Freq == max(Freq)) %>%
select(Var1, Var2, Freq)

new_df = as.data.frame(new_df)
new_df = new_df[which(new_df$Freq>0.5),]
rownames(new_df) = new_df$Var1

meta$celltypes_ext = new_df[paste(meta$`leiden`), "Var2"]
meta$celltypes_ext = as.character(meta$celltypes_ext)
meta$celltypes_ext[which(is.na(meta$celltypes_ext))] = meta$leiden_subcluster[which(is.na(meta$celltypes_ext))]









# Create UMAP scatter plot with ggplot

gc = names(table(m0$predicted_id_sep)[which(table(m0$predicted_id_sep)>50)])


cluster_centers <- m0[which(m0$predicted_id_sep%in%gc),]%>%
  group_by(predicted_id_sep) %>%
  summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

umap_plot <- ggplot(m0[which(m0$predicted_id_sep%in%gc),], aes(x = umap_x, y = umap_y, color = predicted_id_sep)) +
  geom_point() + theme_bw()+
  # Add labels to the center of each cluster
  geom_text(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = predicted_id_sep),color = "black" ,
            position = position_nudge(y = 0.05),size = 2 )  # Adjust the 'y' v
print(umap_plot)


gc = names(table(m0$predicted.id)[which(table(m0$predicted.id)>50)])


cluster_centers <- m0[which(m0$predicted.id%in%gc),]%>%
  group_by(predicted.id) %>%
  summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

umap_plot <- ggplot(m0[which(m0$predicted.id%in%gc),], aes(x = umap_x, y = umap_y, color = predicted.id)) +
  geom_point() + theme_bw()+
  # Add labels to the center of each cluster
  geom_text(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = predicted.id),color = "black" ,
            position = position_nudge(y = 0.05),size = 2 )  # Adjust the 'y' v
print(umap_plot)


umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = celltype)) +
  geom_point() +
  theme_bw()
print(umap_plot)

gc = names(table(m0$celltype)[which(table(m0$celltype)>50)])

cluster_centers <- m0[which(m0$celltype%in%gc),] %>%
  group_by(celltype) %>%
  summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

umap_plot <- ggplot(m0[which(m0$celltype%in%gc),], aes(x = umap_x, y = umap_y, color = celltype)) +
  geom_point() +
  theme_bw() +
  # Add labels to the center of each cluster
  geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = celltype),
                  color = "black", position = position_nudge(y = 0.0), size = 2)
print(umap_plot)

cluster_centers <- m0[which(!is.na(m0$celltypes_ext)),] %>%
  group_by(celltypes_ext) %>%
  summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = celltypes_ext)) +
  geom_point() +
  theme_bw() +
  # Add labels to the center of each cluster
  geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = celltypes_ext),
                  color = "black", position = position_nudge(y = 0.0), size = 2)
print(umap_plot)


cluster_centers <- m0[which(!is.na(m0$CellType_1127_ext)),] %>%
  group_by(CellType_1127_ext) %>%
  summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

umap_plot <- ggplot(m0[which(!is.na(m0$CellType_1127_ext)),], aes(x = umap_x, y = umap_y, color = CellType_1127_ext)) +
  geom_point() +
  theme_bw() +
  # Add labels to the center of each cluster
  geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = CellType_1127_ext),
                  color = "black", position = position_nudge(y = 0.0), size = 2)
print(umap_plot)


cluster_centers <- m0%>%
  group_by(leiden_subcluster) %>%
  summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))


umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = leiden_subcluster)) +
  geom_point() +
  theme_bw() +
  # Add labels to the center of each cluster
  geom_text_repel(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = leiden_subcluster),
                  color = "black", position = position_nudge(y = 0.0), size = 2)
print(umap_plot)


umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = age)) +
  geom_point() +
  theme_bw()
print(umap_plot)


umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = region)) +
  geom_point() +
  theme_bw()
print(umap_plot)

umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = as.factor(present))) +
  geom_point() +
  theme_bw()
print(umap_plot)

umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, alpha = .5,color = doublet_probability)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "red") +  # Customize the color scale
  theme_bw()

print(umap_plot)

cluster_centers <- m0[which(!is.na(m0$predicted_id_sep)),]%>%
  group_by(predicted_id_sep) %>%
  summarize(umap_x_center = mean(umap_x), umap_y_center = mean(umap_y))

umap_plot <- ggplot(m0[which(!is.na(m0$predicted_id_sep)),], aes(x = umap_x, y = umap_y, color = predicted_id_sep)) +
  geom_point() + theme_bw()+
  # Add labels to the center of each cluster
  geom_text(data = cluster_centers, aes(x = umap_x_center, y = umap_y_center, label = predicted_id_sep),color = "black" ,
            position = position_nudge(y = 0.05),size = 2 )  # Adjust the 'y' v
print(umap_plot)

m0 = m0[which(m0$batch == "Female"), ]
pr_tab = table(m0$leiden, m0$present)/rowSums(table(m0$leiden, m0$present))
hist(pr_tab[,2])
pr_tab[which(pr_tab[,2] < .5),, drop = F]
m0$bad = "no"
m0$bad[which(m0$leiden %in% rownames(pr_tab[which(pr_tab[,2] < .12),,drop = F]))]= "yes"
bad_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = bad)) +
  geom_point() +
  theme_bw()
print(bad_plot)

rownames(pr_tab[which(pr_tab[,2] < .1),,drop = F])

umap_plot <- ggplot(m0, aes(x = umap_x, y = umap_y, color = as.factor(present), alpha = 0.5)) +
  geom_point() +
  theme_bw()
print(umap_plot)


library(dplyr)
classification_columns <- grep("^DF.classifications_", colnames(combined_seurat@meta.data), value = TRUE)
metadata = combined_seurat@meta.data
metadata <- metadata %>%
  rowwise() %>%
  mutate(Combined_Classifications = coalesce(!!!syms(classification_columns)))

# Remove the individual classification columns
metadata <- metadata %>% select(-one_of(classification_columns))
head(metadata)



classification_columns <- grep("^pANN", colnames(metadata), value = TRUE)
#metadata = combined_seurat@meta.data
metadata <- metadata %>%
  rowwise() %>%
  mutate(Combined_pANN = coalesce(!!!syms(classification_columns)))

# Remove the individual classification columns
metadata <- metadata %>% select(-one_of(classification_columns))


print("IN")

head(as.numeric(metadata$Combined_pANN))

FeaturePlot(combined_seurat, as.numeric(metadata$Combined_pANN))

atac = read.csv("../meta_bestcelltype.csv")

atac = atac[which(atac$batch == "Female"),]

sub(".*:", "", atac$cell_id)[1:5]


atac$barcode = gsub("Female:", "", atac$cell_id) 
atac$barcode = gsub("-1", "", atac$barcode) 

head(atac)


head(metadata)

metadata$barcode = paste(metadata$orig.ident, ":",sub("-.*", "",colnames(combined_seurat)), sep = "")

metadata$region = sapply(strsplit(metadata$orig.ident, "_"), "[[", 1)
metadata$age = sapply(strsplit(metadata$orig.ident, "_"), "[[", 2)
metadata$rep = sapply(strsplit(metadata$orig.ident, "_"), "[[", 3)

rownames(atac) = atac$barcode

metadata$best_celltype_from_ATAC= atac[paste(metadata$barcode),"best_celltype"]

head(combined_seurat@meta.data)

combined_seurat

write.table(metadata, "RNA_metadata.csv", sep = ",")

metadata = as.data.frame(metadata)
rownames(metadata) = rownames(combined_seurat@meta.data)

combined_seurat@meta.data = metadata

head(colnames(combined_seurat))

options(repr.plot.width=14, repr.plot.height=8)

DimPlot(combined_seurat, group.by = "seurat_clusters", label = T, raster=FALSE)


DimPlot(combined_seurat, group.by = "best_celltype_from_ATAC", label = T, raster=FALSE)


options(repr.plot.width=11, repr.plot.height=8)

DimPlot(combined_seurat, group.by = "age", label = T, raster=FALSE)


options(repr.plot.width=11, repr.plot.height=8)

DimPlot(combined_seurat, group.by = "region", label = T, raster=FALSE)


FeaturePlot(combined_seurat, "Combined_pANN", raster =F)


FeaturePlot(combined_seurat, "nCount_RNA", raster =F)


options(repr.plot.width=13, repr.plot.height=8)

ggplot(metadata, aes(fill = orig.ident , x = factor(seurat_clusters))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggplot(metadata, aes(fill = Combined_Classifications , x = factor(seurat_clusters))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


ggplot(combined_seurat@meta.data, aes(fill =  Combined_Classifications, x = factor(RNA_snn_res.3))) +
  geom_bar(position = "fill", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


doublet_percentage <- aggregate(combined_seurat@meta.data$Combined_Classifications, 
                                 by = list(combined_seurat@meta.data$RNA_snn_res.3), 
                                 FUN = function(x) sum(x == "Doublet") / length(x) * 100)


hist(doublet_percentage$x)

doublet_percentage[which(doublet_percentage$x > 60),'Group.1']

combined_seurat@meta.data$doub_cl = "No"
combined_seurat@meta.data$doub_cl[which(combined_seurat@meta.data$RNA_snn_res.3 %in% doublet_percentage[which(doublet_percentage$x > 60),'Group.1'])] = "Yes"
table(combined_seurat@meta.data$doub_cl)

DimPlot(combined_seurat, group.by = "RNA_snn_res.3", label = T, raster=FALSE)


DimPlot(combined_seurat, group.by = "Combined_Classifications", label = T, raster=FALSE)


DimPlot(combined_seurat, group.by = "doub_cl", label = T, raster=FALSE)


ggplot(combined_seurat@meta.data, aes(fill =  best_celltype_from_ATAC, x = factor(RNA_snn_res.3))) +
  geom_bar(position = "fill", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


length(which(is.na(filtered_seurat@meta.data$best_celltype_from_ATAC)))

ggplot(combined_seurat@meta.data, aes(fill =  age, x = factor(RNA_snn_res.3))) +
  geom_bar(position = "fill", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


table(combined_seurat@meta.data$notATAC, combined_seurat@meta.data$doub_cl)

length(which(!is.na(combined_seurat@meta.data$best_celltype_from_ATAC)))



notATAC_percentage <- aggregate(combined_seurat@meta.data$best_celltype_from_ATAC, 
                                 by = list(combined_seurat@meta.data$RNA_snn_res.3), 
                                 FUN = function(x) sum(is.na(x)) / length(x) * 100)
combined_seurat@meta.data$notATAC = "No"
combined_seurat@meta.data$notATAC[which(combined_seurat@meta.data$RNA_snn_res.3 %in% notATAC_percentage[which(notATAC_percentage$x > 85),'Group.1'])] = "Yes"
table(combined_seurat@meta.data$notATAC)

DimPlot(combined_seurat, group.by = "notATAC", label = T, raster=FALSE)


ggplot(combined_seurat@meta.data[which(combined_seurat@meta.data$notATAC=="Yes"),], aes(fill =  RNA_snn_res.3, x = orig.ident)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + facet_wrap(~notATAC)


ggplot(combined_seurat@meta.data[which(combined_seurat@meta.data$doub_cl=="Yes"),], aes(fill =  RNA_snn_res.3, x = orig.ident)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


lowCounts_percentage <- aggregate(combined_seurat@meta.data$nCount_RNA, 
                                 by = list(combined_seurat@meta.data$RNA_snn_res.3), 
                                 FUN = function(x) sum(x < 350) / length(x) * 100)

combined_seurat@meta.data$lowCounts = "No"
combined_seurat@meta.data$lowCounts[which(combined_seurat@meta.data$RNA_snn_res.3 %in% lowCounts_percentage[which(lowCounts_percentage$x > 50),'Group.1'])] = "Yes"
table(combined_seurat@meta.data$lowCounts)

DimPlot(combined_seurat, group.by = "lowCounts", label = T, raster=FALSE)


ggplot(combined_seurat@meta.data, aes(fill = best_celltype_from_ATAC , x = factor(RNA_snn_res.3))) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


ggplot(metadata, aes( x = factor(seurat_clusters), y = Combined_pANN)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))# + facet_wrap(~region, nrow = 8)

ggplot(metadata, aes( x = factor(orig.ident), y = nCount_RNA)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

ag_pann = aggregate(Combined_pANN ~ seurat_clusters, metadata, FUN = median)

hist(ag_pann$Combined_pANN)

ag_pann[which(ag_pann$Combined_pANN> 0.5),]

combined_seurat <- FindClusters(combined_seurat, resolution = 3)


saveRDS(combined_seurat , file = "prefilter_RNA.RDS")

print("in")

## filtering doublets and rerun umap

filtered_seurat <- subset(combined_seurat, subset = Combined_Classifications == "Singlet")
filtered_seurat
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, nfeatures = 5000)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(filtered_seurat))
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:50)  # You can adjust the number of dimensions
filtered_seurat <- FindClusters(filtered_seurat, resolution = 1)  # Adjust the resolution as needed
filtered_seurat <- RunUMAP(filtered_seurat, reduction = "pca", dims = 1:50)  # Adjust the number of dimensions as needed


filtered_seurat <- FindClusters(filtered_seurat, resolution = 3)  # Adjust the resolution as needed
saveRDS(filtered_seurat , file = "postfilter_RNA.RDS")

filtered_seurat

DimPlot(filtered_seurat, group.by = "region", label = T, raster=FALSE)


filtered_seurat$age = factor(filtered_seurat$age , levels = c("2mo", "9mo", "18mo"))
DimPlot(filtered_seurat, group.by = "age", label = T, raster=FALSE)


DimPlot(filtered_seurat, group.by = "best_celltype_from_ATAC", label = T, raster=FALSE)


DimPlot(filtered_seurat, group.by = "age", label = T, raster=FALSE)


DimPlot(filtered_seurat, group.by = "RNA_snn_res.3", label = T, raster=FALSE)


FeaturePlot(filtered_seurat, "nCount_RNA",raster=FALSE,max.cutoff = 5000)


notATAC_percentage[which(notATAC_percentage$x > 80),'Group.1']

notATAC_percentage <- aggregate(filtered_seurat@meta.data$best_celltype_from_ATAC, 
                                 by = list(filtered_seurat@meta.data$RNA_snn_res.3), 
                                 FUN = function(x) sum(is.na(x)) / length(x) * 100)
filtered_seurat@meta.data$notATAC = "No"
filtered_seurat@meta.data$notATAC[which(filtered_seurat@meta.data$RNA_snn_res.3 %in% notATAC_percentage[which(notATAC_percentage$x > 90),'Group.1'])] = "Yes"
table(filtered_seurat@meta.data$notATAC)

DimPlot(filtered_seurat, group.by = "notATAC", label = T, raster=FALSE)


ggplot(filtered_seurat@meta.data[which(filtered_seurat@meta.data$`RNA_snn_res.3`==53),], aes(fill =  RNA_snn_res.3, x = orig.ident)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + facet_wrap(~notATAC)


ggplot(filtered_seurat@meta.data[which(filtered_seurat@meta.data$notATAC=="Yes"),], aes(fill =  RNA_snn_res.3, x = orig.ident)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + facet_wrap(~notATAC)


ggplot(filtered_seurat@meta.data[which(filtered_seurat@meta.data$notATAC=="No"),], aes(fill =  RNA_snn_res.3, x = orig.ident)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + facet_wrap(~notATAC)


ggplot(filtered_seurat@meta.data, aes( x = RNA_snn_res.3, y = nCount_RNA)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +ylim(c(0,2000))

lowCounts_percentage <- aggregate(filtered_seurat@meta.data$nCount_RNA, 
                                 by = list(filtered_seurat@meta.data$RNA_snn_res.3), 
                                 FUN = function(x) sum(x < 1000) / length(x) * 100)

filtered_seurat@meta.data$lowCounts = "No"
filtered_seurat@meta.data$lowCounts[which(filtered_seurat@meta.data$RNA_snn_res.3 %in% lowCounts_percentage[which(lowCounts_percentage$x > 50),'Group.1'])] = "Yes"
table(filtered_seurat@meta.data$lowCounts)

lowCounts_percentage[which(lowCounts_percentage$x > 50),'Group.1']

DimPlot(filtered_seurat, group.by = "notATAC", label = T, raster=FALSE)


table(filtered_seurat@meta.data$notATAC)

filtered_seurat = readRDS("postfilter_RNA.RDS")

notATAC_percentage <- aggregate(filtered_seurat@meta.data$best_celltype_from_ATAC, 
                                 by = list(filtered_seurat@meta.data$RNA_snn_res.3), 
                                 FUN = function(x) sum(is.na(x)) / length(x) * 100)
filtered_seurat@meta.data$notATAC = "No"
filtered_seurat@meta.data$notATAC[which(filtered_seurat@meta.data$RNA_snn_res.3 %in% notATAC_percentage[which(notATAC_percentage$x > 90),'Group.1'])] = "Yes"
table(filtered_seurat@meta.data$notATAC)

## filtering doublets and rerun umap

filtered_seurat <- subset(filtered_seurat, subset = notATAC == "No")
filtered_seurat
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, nfeatures = 5000)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(filtered_seurat))
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:50)  # You can adjust the number of dimensions
filtered_seurat <- FindClusters(filtered_seurat, resolution = 1)  # Adjust the resolution as needed
filtered_seurat <- RunUMAP(filtered_seurat, reduction = "pca", dims = 1:50)  # Adjust the number of dimensions as needed
saveRDS(filtered_seurat , file = "postfilter_RNA.RDS")

DimPlot(filtered_seurat, group.by = "best_celltype_from_ATAC", label = T, raster=FALSE)


length(which(!is.na(filtered_seurat$best_celltype_from_ATAC)))

filtered_seurat = readRDS("postfilter_RNA.RDS")

filtered_seurat

options(repr.plot.width=12, repr.plot.height=8)

DimPlot(filtered_seurat, group.by = "RNA_snn_res.3", label = T, raster=FALSE)


ggplot(filtered_seurat@meta.data[which(filtered_seurat@meta.data$`RNA_snn_res.3`==93),], aes(fill =  predicted.id, x = orig.ident)) +
  geom_bar(position = "stack", stat = "count") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 


options(repr.plot.width=13, repr.plot.height=8)

DimPlot(filtered_seurat, group.by = "best_celltype_from_ATAC", label = T, raster=FALSE)


filtered_seurat$age = factor(filtered_seurat$age, levels = c("2mo" , "9mo", "18mo"))
DimPlot(filtered_seurat, group.by = "age", label = F, raster=T)


options(repr.plot.width=10, repr.plot.height=8)

FeaturePlot(filtered_seurat, "nCount_RNA",raster=FALSE,max.cutoff = 5000)


FeaturePlot(filtered_seurat, "Combined_pANN",raster=FALSE,max.cutoff = 5000)


filtered_seurat

filtered_seurat = readRDS("postfilter_RNA.RDS")

options(repr.plot.width=20, repr.plot.height=8)

DimPlot(filtered_seurat, group.by = "RNA_snn_res.3", label = T, raster=FALSE)


filtered_seurat
filtered_seurat <- subset(filtered_seurat, subset = RNA_snn_res.3 != "53")
filtered_seurat
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, nfeatures = 5000)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(filtered_seurat))
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:50)  # You can adjust the number of dimensions
filtered_seurat <- FindClusters(filtered_seurat, resolution = 1)  # Adjust the resolution as needed
filtered_seurat <- RunUMAP(filtered_seurat, reduction = "pca", dims = 1:50)  # Adjust the number of dimensions as needed
saveRDS(filtered_seurat , file = "postfilter2_RNA.RDS")

DimPlot(filtered_seurat, group.by = "age", label = T, raster=FALSE)


DimPlot(filtered_seurat, group.by = "age", label = T, raster=FALSE)


obj$filt_lowQ = "No"
obj$filt_lowQ[which(obj$nFeature_RNA<250 | obj$nCount_RNA < 300)] = "Yes"
obj$filt_lowQ[which(obj$`RNA_snn_res.3`==93)] = "Yes"
filtered_seurat <- subset(obj, subset = filt_lowQ == "No")
filtered_seurat
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, nfeatures = 5000)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(filtered_seurat))
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:50)  # You can adjust the number of dimensions
filtered_seurat <- FindClusters(filtered_seurat, resolution = 1)  # Adjust the resolution as needed
filtered_seurat <- FindClusters(filtered_seurat, resolution = 3)  # Adjust the resolution as needed
filtered_seurat <- RunUMAP(filtered_seurat, reduction = "pca", dims = 1:50)  # Adjust the number of dimensions as needed
saveRDS(filtered_seurat, "postfilter3.RDS")

### TRANSFER LABELS

allenSeu = readRDS("/projects/ps-renlab2/szu/projects/CEMBA2/19.snap2_integration/out/transferLabel_tscc/allen_ds1000_seurat.rds")

allenSeu

allenSeu <- NormalizeData(allenSeu)
allenSeu <- FindVariableFeatures(allenSeu, nfeatures = 5000)
allenSeu <- ScaleData(allenSeu)
allenSeu <- RunPCA(allenSeu, features = VariableFeatures(allenSeu))


anchors <- FindTransferAnchors(reference = allenSeu, query = filtered_seurat,
    dims = 1:40, reference.reduction = "pca", features = VariableFeatures(filtered_seurat))


predictions <- TransferData(anchorset = anchors, refdata = allenSeu$subclass,
    dims = 1:30)

filtered_seurat <- AddMetaData(filtered_seurat, metadata = predictions)


options(repr.plot.width=12, repr.plot.height=8)

DimPlot(filtered_seurat, group.by = "predicted.id", label = T, raster=FALSE,legend = T)


pdf("subclass_UMAP.pdf", height = 8, width = 50)
DimPlot(filtered_seurat, group.by = "predicted.id", label = T, raster=FALSE)
dev.off()

p = DimPlot(filtered_seurat, group.by = "predicted.id",label = T, label.size = .2)

p + theme(legend.position = "none")

nam= names(table(filtered_seurat$predicted.id))[which(table(filtered_seurat$predicted.id)<50)]

filtered_seurat$id_other = filtered_seurat$predicted.id
filtered_seurat$id_other[which(filtered_seurat$id_other %in%nam)] = "other"

DimPlot(filtered_seurat, group.by = "id_other", label = T, raster=FALSE)


head(allenSeu@meta.data)

library(data.table)
ann = fread("/projects/ps-renlab2/szu/projects/shared_data/allen/AIT21_annotation_freeze_081523.tsv")

length(table(ann$subclass_id_label))

length(table(filtered_seurat$predicted.id))


