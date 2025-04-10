library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(Seurat)

ameta = read.csv("../final_meta.csv")


age_color = read.csv("color_scheme/Age.csv")
region_color = read.csv("color_scheme/MajorRegion.csv")
celltype_color = read.csv("color_scheme/updated_celltype_palette.csv")
modality_color = read.csv("color_scheme/Modality.csv") 

rownames(age_color) = age_color$Age
rownames(region_color) = region_color$Region
rownames(celltype_color) = celltype_color$CellType
rownames(modality_color) = modality_color$Modality

obj = readRDS("../female_RNA/RNA_final_SCT.RDS")

obj$cell_id = paste("Female:", obj$barcode, "-1",sep = "")


umap_coordinates <- Embeddings(obj, reduction = "umap")

umap_coordinates = as.data.frame(umap_coordinates)
umap_coordinates$cell_id = obj$cell_id
#umap_coordinates = umap_coordinates[-which(is.na(umap_coordinates$cell_id)),]
rownames(umap_coordinates) = umap_coordinates$cell_id
nrow((umap_coordinates))

obj@meta.data["umap_1"] = umap_coordinates$umap_1
obj@meta.data["umap_2"] = umap_coordinates$umap_2

meta= obj@meta.data
meta$celltype_final[which(meta$final_clusters_final=="IOL")] = "IOL NN"
meta$region_name = ""
meta$region_name[which(meta$region=="HCA")] = "Anterior_Hippocampus"
meta$region_name[which(meta$region=="HCP")] = "Posterior_Hippocampus"
meta$region_name[which(meta$region=="ENT")] = "Entorhinal_Cortex"
meta$region_name[which(meta$region=="AMY")] = "Amygdala"
meta$region_name[which(meta$region=="FC")] = "Frontal_Cortex"
meta$region_name[which(meta$region=="NAC")] = "Nucleus_accumbens"
meta$region_name[which(meta$region=="CP")] = "Caudate_Putamen"
meta$region_name[which(meta$region=="RLP")] = "PAG-PCG"

sample = meta[sample(1:nrow(meta),100000),]


length(which(meta$final_clusters_final=="IOL"))

options(repr.plot.width=18, repr.plot.height=10)
color_vector <- setNames(celltype_color$Color, celltype_color$CellType)

ggplot(sample, aes(x = umap_1, y=umap_2, color = celltype_final)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 

options(repr.plot.width=14, repr.plot.height=12)
color_vector <- setNames(celltype_color$Color, celltype_color$CellType)
ord <- names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)

ggplot(meta, aes(x = celltype_final, fill = celltype_final)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + ylim(c(0,100000)) + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() +
  geom_text(stat='count', aes(label=after_stat(count)), position=position_stack(vjust=1), hjust=-0.1)

options(repr.plot.width=6, repr.plot.height=5)
color_vector <- setNames(age_color$Color, age_color$Age)
meta$age = factor(meta$age, levels = c("2mo", "9mo", "18mo"))

ggplot(sample, aes(x = umap_1, y=umap_2, color = age)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 



ggplot(meta, aes(x = region_name, fill = age)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + ylim(c(0,75000)) + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() 


options(repr.plot.width=6, repr.plot.height=15)


ggplot(meta, aes(x = celltype_final, fill = age)) +
  geom_bar(stat = "count",position="fill", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() 

options(repr.plot.width=7, repr.plot.height=5)
color_vector <- setNames(region_color$Color, region_color$Region)

ggplot(sample, aes(x = umap_1, y=umap_2, color = region_name)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() 


ggplot(meta, aes(x = region_name, fill = region_name)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + ylim(c(0,100000)) + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() +
  geom_text(stat='count', aes(label=after_stat(count)), position=position_stack(vjust=1), hjust=-0.1)



options(repr.plot.width=14, repr.plot.height=12)
color_vector <- setNames(celltype_color$Color, celltype_color$CellType)
ord <- names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)

ggplot(meta, aes(x = celltype_final, fill = celltype_final)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + ylim(c(0,100000)) + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() +
  geom_text(stat='count', aes(label=after_stat(count)), position=position_stack(vjust=1), hjust=-0.1)
