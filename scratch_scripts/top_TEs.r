library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)

cts = read.table("~/projects/combined_all/Figures/celltypes_figure6.txt")
cts = cts$x
cts

meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")
meta$celltype_final = gsub("Lymphoid", "Lymphoid NN", meta$celltype_final)

lkey=read.table("~/projects/combined_all/scTE/out.saf")
#lkey = lkey[-which(lkey$V5=="Retroposon"),]
ukey = lkey[-which(duplicated(lkey$V1)),]
rownames(ukey) = ukey$V1

bkey  = ukey[grep("_", rownames(ukey)), ]
rownames(bkey) = gsub("_", "-", rownames(bkey))
ukey = rbind(ukey,bkey)
ukey$V6 = gsub("[?]", "", ukey$V6)

setwd("~/projects/combined_all/female_RNA/SoloTE/subfamily_DESeq2/")

files = list.files(".", "only.txt", full.names = F)
cur = files[1]
tab = read.table(cur)
tab$gene = rownames(tab)
f = gsub("_DE_TEs_only.txt", "", files[1])

tab$celltype = f
#head(tab)
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_DE_TEs_only.txt", "", files[i])

    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}

tab$logFC = -tab$logFC
tab$gene = gsub("SoloTE-", "", tab$gene)

tab$type = "Gene"
tab$type[grep("SoloTE", rownames(tab))] = "SoloTE"
tab$dir = "up"
tab$dir[which(tab$logFC<0)] = "down"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("IMN", tab$celltype)]="NN"

#tab= tab[-grep("IT", tab$celltype),]

    
    
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","_", ord)


options(repr.plot.width=17, repr.plot.height=4)

tab = tab[which(tab$celltype%in% cts), ] 
rna= tab
d = tab
d = d[which(d$PValue<0.01 & abs(d$logFC) > 0.25),]
d = d[which(d$type!="Gene"),]
d$celltype = factor(d$celltype, levels = rev(ord))


options(repr.plot.width=7, repr.plot.height=8)

g1 = ggplot(d, aes(x = celltype, ,fill = dir)) +
  geom_bar(stat = "count",position = "stack", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_grid(clade~., scales = "free",space = "free")+ ggtitle(paste( "P < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))
print(g1)
    

setwd("~/projects/combined_all/scTE/celltype_DEseq2_results/")

files = list.files(".", "txt", full.names = F)
f = files[1]
tab = read.table(f)
tab$gene = rownames(tab)
f = gsub("_2vs18.txt", "", files[1])

tab$celltype = f
#head(tab)
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_2vs18.txt", "", files[i])

    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}

tab$logFC = -tab$logFC
tab$type = ukey[paste(tab$gene),"V6"]

tab$celltype = gsub("L2_3", "L2-3", tab$celltype)
tab$celltype = gsub("L2_3", "L2-3", tab$celltype)
tab$celltype = gsub("L5_6", "L5-6", tab$celltype)
tab$celltype = gsub("L6b_CT", "L6b-CT", tab$celltype)
tab$celltype = gsub("Lymphoid", "Lymphoid_NN", tab$celltype)


tab = tab[which(tab$celltype!="doublet"),]
unique(tab$celltype)[which(!unique(tab$celltype) %in% ord)]

tab$type = "Gene"
tab$type[grep("SoloTE", rownames(tab))] = "SoloTE"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("IMN", tab$celltype)]="NN"

ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","_", ord)


options(repr.plot.width=17, repr.plot.height=4)

tab = tab[which(tab$celltype%in% cts), ] 

d = tab
d = d[which(d$fdr<0.05 & abs(d$logFC) > 0.25),]
d$celltype = factor(d$celltype, levels = rev(ord))

options(repr.plot.width=17, repr.plot.height=7)

options(repr.plot.width=7, repr.plot.height=8)

g1 = ggplot(d, aes(x = celltype, ,fill = dir)) +
  geom_bar(stat = "count",position = "stack", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_grid(clade~., scales = "free",space = "free")+ ggtitle(paste( "P < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))
print(g1)
    

unique(tab$gene)%in%rna$gene

sig= tab[which(tab$fdr<.00001 & abs(tab$logFC) > 0.5 ),]
sig_rna= rna[which(rna$PValue<.001 & abs(rna$logFC) > 0.15),]
sig_genes = (unique(sig$gene))
sig_genes_rna = (unique(sig_rna$gene))
length(sig_genes_rna)
length(sig_genes)
sig_genes = sig_genes_rna[which(sig_genes_rna %in% sig_genes)]
length(sig_genes)
sig_genes

tab$type = "ATAC"
rna$type = "RNA"



combined = rbind(tab, rna)
combined$TE_type = ukey[paste(combined$gene),6]


combined = combined[which(combined$gene%in%sig_genes),]

head(combined)

options(repr.plot.width=8, repr.plot.height=9)

# Convert PValue into significance stars
data <- combined %>%
  mutate(significance_star = ifelse(PValue < 0.01, "*", ""))
data$logFC[which(data$logFC>1)] = 1
data$logFC[which(data$logFC< -1)] = -1

data$celltype = factor(as.character(data$celltype))
data$celltype = factor(data$celltype , levels = rev(levels(data$celltype)))
data$clade = factor(data$clade, levels = (c("NN", "Gaba", "Glut")))


# Create heatmap
g = ggplot(data, aes(x = gene, y = celltype, fill = logFC)) +
  geom_tile() +  # Heatmap tiles
  geom_text(aes(label = significance_star), color = "black", size = 5) +  # Add significance stars
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Diverging colors
  facet_grid(clade~type, scales = "free", space = "free") +  # Split by ATAC vs RNA
  theme_bw() +
  labs(
    title = "Top Age-differential TEs",
    x = "Gene",
    y = "Cell Type",
    fill = "Age-LogFC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), text = element_text(size = 14),
        legend.position = "right") 
g

g 
pdf("~/projects/combined_all/Figures/Figure6-Het-TEs/Top_TEs_heatmap.pdf", height =9, width = 7)
print(g)
dev.off()

library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)


heatmap_data <- combined %>%
  select(gene, celltype, type, logFC, TE_type, PValue) %>%
  pivot_wider(names_from = gene, values_from = logFC)

# Extract row annotation (TE_type)
row_annotation <- data %>%
  distinct(celltype, TE_type) %>%
  as.data.frame()
rownames(row_annotation) <- row_annotation$celltype
row_annotation$celltype <- NULL  # Remove redundant column

# Create significance matrix for star annotations
significance_matrix <- data %>%
  select(gene, celltype, PValue) %>%
  pivot_wider(names_from = gene, values_from = PValue) %>%
  mutate(across(where(is.numeric), ~ ifelse(. < 0.05, "*", "")))
rownames(significance_matrix) <- significance_matrix$celltype
significance_matrix$celltype <- NULL  # Remove redundant column

# Split data into ATAC and RNA
heatmap_ATAC <- heatmap_data %>% filter(type == "ATAC") %>% select(-type)
heatmap_RNA <- heatmap_data %>% filter(type == "RNA") %>% select(-type)

# Convert to matrix for pheatmap
heatmap_matrix_ATAC <- as.matrix(heatmap_ATAC[, -1])  # Remove celltype column
rownames(heatmap_matrix_ATAC) <- heatmap_ATAC$celltype

heatmap_matrix_RNA <- as.matrix(heatmap_RNA[, -1])  # Remove celltype column
rownames(heatmap_matrix_RNA) <- heatmap_RNA$celltype

# Define colors for annotation
ann_colors <- list(
  TE_type = c(LTR = "red", SINE = "blue", LINE = "green", DNA = "purple")
)

# Plot ATAC heatmap
pheatmap(
  heatmap_matrix_ATAC, 
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  cluster_rows = TRUE, cluster_cols = TRUE,
  display_numbers = significance_matrix[rownames(heatmap_matrix_ATAC), , drop = FALSE],  # Add stars for significance
  main = "ATAC"
)

# Plot RNA heatmap
pheatmap(
  heatmap_matrix_RNA, 
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  cluster_rows = TRUE, cluster_cols = TRUE,
  display_numbers = significance_matrix[rownames(heatmap_matrix_RNA), , drop = FALSE],  # Add stars for significance
  main = "RNA"
)



duplicates <- data %>%
  group_by(gene, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1L)
print(duplicates)

data = data[-which(is.na(data$fdr)),]
sig_matrix <- data %>%
  select(gene, celltype, fdr) %>%
  pivot_wider(names_from = celltype, values_from = fdr) 

options(repr.plot.width=12, repr.plot.height=15)

rna$TE_type = ukey[paste(rna$gene),6]

# Add a column for significance stars
data <- rna %>%
  mutate(star = ifelse(significance == "significant", "*", ""))
data$logFC[which(data$PValue>0.1)] = 0
head(data)
data = data[which(!is.na(data$gene)),]
sig_matrix <- data %>%
  select(gene, celltype, PValue) %>%
  pivot_wider(names_from = celltype, values_from = PValue) %>%
  as.data.frame()
head(sig_matrix)

rownames(sig_matrix) = sig_matrix$gene
sig_matrix = sig_matrix[,-c(1)] 
good = which(rowMins(as.matrix(sig_matrix))<0.1)



# Prepare the heatmap matrix (genes as rows, cell types as columns)
heatmap_matrix <- data %>%
  select(gene, celltype, logFC) %>%
  pivot_wider(names_from = celltype, values_from = logFC) %>%
  as.data.frame()
heatmap_matrix = heatmap_matrix[good,]
head(heatmap_matrix)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix <- heatmap_matrix[, -1]  # Remove gene column after setting rownames
heatmap_matrix <- as.matrix(heatmap_matrix)
head(heatmap_matrix)
# Create row annotations for TE_type
row_annotation <- data.frame(TE_type = data$TE_type, row.names = rownames(data))
head(row_annotation)
min_val <- -1 # Adjust based on your data range for increased intensity
max_val <- 1

breaks <- c(
  seq(min_val, 0, length.out = 26),  # From minimum to 0
  seq(0, max_val, length.out = 25)[-1]  # From 0 to maximum
)
heatmap_matrix = heatmap_matrix[rowSums(abs(heatmap_matrix))>0,]
heatmap_matrix = heatmap_matrix[,colSums(abs(heatmap_matrix))>0]


#col_annotation <- data.frame(clade = unique(data$clade), row.names = unique(data$celltype))
# Create the heatmap
pheatmap(clustering_method = "ward.D2",
  heatmap_matrix,
  cluster_rows = TRUE,  # Cluster rows
  cluster_cols = TRUE,  # Cluster columns
  #annotation_row = row_annotation,  # Add row annotation
 # annotation_col = col_anno, 
    color = colorRampPalette(c("blue", "white", "red"))(50),
      breaks = breaks,

  number_color = "black",
  main = "Heatmap of logFC with Clustering and Annotations"
)


options(repr.plot.width=7, repr.plot.height=8)

tab$TE_type = ukey[paste(tab$gene),6]

# Add a column for significance stars
data <- tab %>%
  mutate(star = ifelse(significance == "significant", "*", ""))
#data$logFC[which(data$fdr>.1)] = 0

sig_matrix <- data %>%
  select(gene, celltype, fdr) %>%
  pivot_wider(names_from = celltype, values_from = fdr) %>%
  replace(is.na(.), 1) %>%  # Replace NA with 1
  as.data.frame()
rownames(sig_matrix) = sig_matrix$gene
sig_matrix = sig_matrix[,-c(1)] 

sig_matrix_og = sig_matrix

good = which(rowMins(as.matrix(sig_matrix_og))<1e-12)
good_col = which(colMins(as.matrix(sig_matrix_og))<1e-5)

sig_matrix = as.data.frame(sig_matrix)
sig_matrix[which(sig_matrix<=0.001,arr.ind = T)] = '***'
sig_matrix[which(sig_matrix<=0.05,arr.ind = T)] = '*'
sig_matrix[which(sig_matrix>0.05,arr.ind = T)] = ''

# Prepare the heatmap matrix (genes as rows, cell types as columns)
heatmap_matrix <- data %>%
  select(gene, celltype, logFC) %>%
  pivot_wider(names_from = celltype, values_from = logFC) %>%
  replace(is.na(.), 0) %>%  # Replace NA with 0
  as.data.frame()
heatmap_matrix = heatmap_matrix[good,]


rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix <- heatmap_matrix[, -1]  # Remove gene column after setting rownames
heatmap_matrix <- as.matrix(heatmap_matrix)
head(heatmap_matrix)
# Create row annotations for TE_type
row_annotation <- data.frame(TE_type = data$TE_type, row.names = rownames(data))
head(row_annotation)
min_val <- -1 # Adjust based on your data range for increased intensity
max_val <- 1

breaks <- c(
  seq(min_val, 0, length.out = 26),  # From minimum to 0
  seq(0, max_val, length.out = 25)[-1]  # From 0 to maximum
)
#heatmap_matrix = heatmap_matrix[rowSums(abs(heatmap_matrix))>.75,]
#heatmap_matrix = heatmap_matrix[,colSums(abs(heatmap_matrix))>0.75]
heatmap_matrix= heatmap_matrix[,names(good_col)]

sig_matrix = sig_matrix[rownames(heatmap_matrix), colnames(heatmap_matrix)]
#col_annotation <- data.frame(clade = unique(data$clade), row.names = unique(data$celltype))

annotation_colors <- list(
  clade = c(
    "Glut" = "firebrick1",
    "Gaba" = "green2",
    "NN" = "blue3"  ), 
    TE_type = c("LINE"  = '#F8766D', "LTR"= '#00BA38')
)

# Create the heatmap
pheatmap(display_numbers = t(sig_matrix),
         clustering_method = "ward",
  t(heatmap_matrix), 
  cluster_rows = TRUE,  # Cluster rows
  cluster_cols = TRUE,  # Cluster columns
  annotation_row= col_anno, 
  annotation_col = row_annotation,  # Add row annotation
annotation_colors = annotation_colors,
    color = colorRampPalette(c("blue", "white", "red"))(50),
      breaks = breaks,

  number_color = "black",
  main = "TE accessibility age-logFC "
)



