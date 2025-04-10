library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(Seurat)

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


write.table(file = "rna_meta.csv", obj@meta.data)

DimPlot(obj, group.by = "celltype_final", label =T)+NoLegend()

umap_coordinates <- Embeddings(obj, reduction = "umap")

umap_coordinates = as.data.frame(umap_coordinates)
umap_coordinates$cell_id = obj$cell_id
#umap_coordinates = umap_coordinates[-which(is.na(umap_coordinates$cell_id)),]
rownames(umap_coordinates) = umap_coordinates$cell_id
nrow((umap_coordinates))

obj@meta.data["umap_x"] = umap_coordinates$umap_1
obj@meta.data["umap_y"] = umap_coordinates$umap_2

length(which(meta$final_clusters_final=="IOL"))

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


options(repr.plot.width=18, repr.plot.height=10)
color_vector <- setNames(celltype_color$Color, celltype_color$CellType)

ggplot(sample, aes(x = umap_x, y=umap_y, color = celltype_final)) +
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

color_vector <- setNames(celltype_color$Color, celltype_color$CellType)

g = ggplot(sample, aes(x = umap_x, y=umap_y, color = celltype_final)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(legend.position = "none") 
g

pdf("celltype_umap_RNA.pdf", width=6, height =6)
print(g)
dev.off()
pdf("celltype_umap_RNA_labeled.pdf", width=10, height =10)
DimPlot(obj, group.by = "celltype_final", label =T)+NoLegend()
dev.off()

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


g = ggplot(sample, aes(x = umap_x, y=umap_y, color = age)) +
  geom_point(alpha=.8, size = .2) +
  scale_color_manual(values = color_vector) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_classic() +
  theme(legend.position = "none") 
g

pdf("age_umap_RNA.pdf", width=6, height =6)
print(g)
dev.off()

options(repr.plot.width=6, repr.plot.height=15)

color_vector <- setNames(age_color$Color, age_color$Age)

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

obj@meta.data = meta

fam = FindAllMarkers(obj, logfc.threshold = .5, only.pos = T,return.thresh = .00001)

top_markers <- fam %>%
  filter(!grepl("^Gm", gene) & !grepl("Rik$", gene)) %>%
  arrange(desc(avg_log2FC)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  group_by(cluster) %>%
  slice(1)

which(duplicated(top_markers$cluster))

#top_markers = top_markers[-which(duplicated(top_markers$cluster)),]
#top_markers = top_markers[-which(duplicated(top_markers$gene)),]

avg_exp <- AverageExpression(obj, return.seurat = TRUE)


dist_matrix <- dist(t(avg_exp@assays$RNA$data[unique(fam$gene),]))
clustering <- hclust(dist_matrix)
ordered_celltypes <- clustering$labels[clustering$order]


dendro_data <- as.dendrogram(clustering)
dendro_plot <- ggdendrogram(dendro_data, rotate = TRUE, theme_dendro = FALSE)


# Fetch the expression data for the top markers
expr_data <- FetchData(obj, vars = top_markers$gene)

# Calculate the percentage of cells expressing each gene in each cell type
expr_data$celltype <- obj$celltype_final

percent_exp_df <- expr_data %>%
  gather(gene, expression, -celltype) %>%
  group_by(celltype, gene) %>%
  summarize(percent_exp = mean(expression > 0) * 100,
            avg_exp = mean(expression))

# Ensure celltype is a factor and ordered correctly
percent_exp_df$celltype <- factor(percent_exp_df$celltype, levels = ordered_celltypes)


gene_order <- top_markers %>%
  arrange(match(cluster, ordered_celltypes)) %>%
  pull(gene)

percent_exp_df$gene <- factor(percent_exp_df$gene, levels = gene_order)

# Create the dot plot using ggplot2
dotplot <- ggplot(percent_exp_df, aes(x = gene, y = celltype, size = percent_exp, color = avg_exp)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(1, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cell Type Marker Gene Dot Plot")

# Plot the dot plot
print(dotplot)


options(repr.plot.width=25, repr.plot.height=13)

grid.arrange(dotplot, dendro_plot+theme_classic(), ncol = 2)

options(repr.plot.width=6, repr.plot.height=15)

color_vector <- setNames(age_color$Color, age_color$Age)
meta$celltype_final=factor(meta$celltype_final, levels = ordered_celltypes)


ggplot(meta, aes(x = celltype_final, fill = age)) +
  geom_bar(stat = "count",position="fill", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() 

Idents(obj) = obj$celltype_final
av = AggregateExpression(obj,normalization.method = "RC")
rn = av$RNA
rn = rn[VariableFeatures(obj),]
df = as.data.frame(rn)
#cor(df)
pheatmap(cor(df))

pheatmap(cor(df),clustering_method = "ward.D2")




