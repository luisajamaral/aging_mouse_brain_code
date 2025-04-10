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
#tab$TE_type = ukey[paste(tab$gene), 6]


tab$type = "Gene"
tab$type[grep("SoloTE", rownames(tab))] = "SoloTE"
tab$direction = "Up in aging"
tab$direction[which(tab$logFC<0)] = "Down in aging"
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


d = tab
d = d[which(d$PValue<0.05),]
d = d[which(d$type!="Gene"),]
d$celltype = factor(d$celltype, levels = rev(ord))


g1 = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = position_stack(), color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "P < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))
print(g1)
    
#d = tab
##d = d[which(d$PValue<0.05),]
#d = d[which(d$type=="Gene"),]
#d$celltype = factor(d$celltype, levels = rev(ord))


#g = ggplot(d, aes(x = celltype, ,fill = direction)) +
#  geom_bar(stat = "count",position = "dodge", color = "black") +
#  theme_minimal() +
#  labs(title = "",
#       x = " Type",
#       y = "Number of DE") +
#  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "P < 0.05 genes"))+
#  theme(text = element_text(size = 18))
    
    #print(g)




options(repr.plot.width=18, repr.plot.height=8)

ggplot(d, aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", position = position_dodge(), color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "P < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))

rna = tab

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
#tab$direction = "Up in aging"
#tab$direction[which(tab$logFC>0)] = "Down in aging"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("IMN", tab$celltype)]="NN"
#tab$clade[grep("Inh-IMN", tab$celltype)]="Gaba"

#tab= tab[-grep("IT", tab$celltype),]

    
    
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","_", ord)


options(repr.plot.width=17, repr.plot.height=4)


d = tab
d = d[which(d$fdr<0.05),]
d$celltype = factor(d$celltype, levels = rev(ord))

options(repr.plot.width=17, repr.plot.height=7)





options(repr.plot.width=18, repr.plot.height=8)

ggplot(d, aes(x = celltype, fill = dir)) +
  geom_bar(stat = "count", position = position_dodge(), color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "fdr < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))

sig= tab[which(tab$fdr<1 & abs(tab$logFC) > 0.1 ),]
sig_rna= rna[which(rna$fdr<1 & abs(rna$logFC) > 0.1),]


sig_genes = (unique(sig$gene))
length(sig_genes)

sig_genes_rna = (unique(sig_rna$gene))
length(sig_genes_rna)

sig_genes = sig_genes_rna[which(sig_genes_rna %in% sig_genes)]
length(sig_genes)

tab$type = "ATAC"
rna$type = "RNA"


#rna = rna[,-c(11)]

head(rna)
head(tab)


rna$gene = ukey[paste(rna$gene ),1]
head(rna)

data = tab

col_anno = data[,c("celltype", "clade")]
head(col_anno)
col_anno = col_anno[-which(duplicated(col_anno$celltype)),]
col_anno$celltype = gsub(" ", "_", col_anno$celltype)
rownames(col_anno) = col_anno$celltype
col_anno = col_anno[,-c(1), drop = F]
head(col_anno)


head(col_anno)

head(rna)
head(tab)

head(data)
data$fdr = as.numeric(data$fdr)
head(data)

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


#head(heatmap_matrix)
pheatmap(clustering_method = "ward.D2",
  heatmap_matrix,
  cluster_rows = TRUE,  # Cluster rows
  cluster_cols = TRUE,  # Cluster columns
  annotation_row = row_annotation,  # Add row annotation
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


# Perform clustering via pheatmap to extract the clustering results
clustering <- pheatmap(
  t(heatmap_matrix),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "ward",
  silent = TRUE  # Suppress plot output for now
)

# Reverse the column order in the clustering result
clustering$tree_col$order <- rev(clustering$tree_col$order)

# Generate the heatmap with reversed column order
p = pheatmap(
  display_numbers = t(sig_matrix),
  t(heatmap_matrix),
  cluster_rows = clustering$tree_row,  # Use the precomputed row clustering
  cluster_cols = clustering$tree_col,  # Use the modified column clustering
  annotation_row = col_anno,
  annotation_col = row_annotation,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  breaks = breaks,
  number_color = "black",
  main = "TE accessibility age-logFC"
)

#ggsave(p , width = 17, height = 12.5,file = "../../Figures/Figure6-Het-TEs/TE_accessibility_LogFC_heatmap.pdf")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
gg_color_hue(3)

sig_matrix <- data %>%
  select(gene, celltype, fdr) %>%
  pivot_wider(names_from = celltype, values_from = fdr) %>%
  replace(is.na(.), 1) %>%  # Replace NA with 1
  as.data.frame()
rownames(sig_matrix) = sig_matrix$gene
sig_matrix = sig_matrix[,-c(1)] 
good = which(rowMins(as.matrix(sig_matrix))<1e-6)
head(sig_matrix)

sig_matrix = as.data.frame(sig_matrix)
sig_matrix[which(sig_matrix<=0.001,arr.ind = T)] = '***'
sig_matrix[which(sig_matrix<=0.05,arr.ind = T)] = '*'
sig_matrix[which(sig_matrix>0.05,arr.ind = T)] = ''

head(sig_matrix)

library(pheatmap)
library(dplyr)
library(tidyr)
library(matrixStats)

# Create a significance matrix with stars
stars_matrix <- data %>%
  select(gene, celltype, fdr) %>%
  pivot_wider(names_from = celltype, values_from = fdr) %>%
  replace(is.na(.), 1) %>%  # Replace NA with 1 (not significant)
  mutate(across(everything(), ~ ifelse(. < 0.05, "*", ""))) %>%  # Add stars for p < 0.05
  as.data.frame()

rownames(stars_matrix) <- stars_matrix$gene
stars_matrix <- stars_matrix[, -1]  # Remove gene column

# Filter the stars_matrix to match the heatmap_matrix
stars_matrix <- stars_matrix[rownames(heatmap_matrix), colnames(heatmap_matrix)]

# Generate heatmap with stars
pheatmap(
  t(heatmap_matrix),  # Transpose to match clustering expectations
  cluster_rows = TRUE,  # Cluster rows
  cluster_cols = TRUE,  # Cluster columns
  clustering_method = "ward.D2",  # Clustering method
  annotation_row = col_anno,  # Column annotation
  annotation_col = row_annotation,  # Row annotation
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Color scheme
  breaks = breaks,  # Custom breaks
  display_numbers = t(stars_matrix),  # Add stars for significance
  number_color = "black",  # Color of the stars
  main = "Heatmap of logFC with Significance"  # Heatmap title
)


head(sig_matrix)


library(MatrixGenerics)
head(sig_matrix)
rownames(sig_matrix) = sig_matrix$gene
sig_matrix = sig_matrix[,-c(1)] 
rowMins(as.matrix(sig_matrix))>0.01

max_val

heatmap_matrix = heatmap_matrix[rowSums(abs(heatmap_matrix))>0.5,]

# Add a column for significance stars
data <- tab %>%
  mutate(star = ifelse(significance == "significant", "*", ""))
data$TE_type = ukey[paste(data$gene),6]

# Prepare the heatmap matrix (genes as rows, cell types as columns)
heatmap_matrix <- data %>%
  select(gene, celltype, logFC) %>%
  pivot_wider(names_from = celltype, values_from = logFC) %>%
  replace(is.na(.), 0) %>%  # Replace NA with 0
  as.data.frame()

rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix <- heatmap_matrix[, -1]  # Remove gene column after setting rownames
heatmap_matrix <- as.matrix(heatmap_matrix)

# Create row annotations for TE_type
row_annotation <- data %>%
  distinct(gene, TE_type) %>%
  as.data.frame()

rownames(row_annotation) <- row_annotation$gene
row_annotation <- row_annotation[, -1]  # Remove the gene column

# Create the heatmap
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,  # Cluster rows
  cluster_cols = TRUE,  # Cluster columns
  annotation_row = row_annotation,  # Add row annotation
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Custom color scale
  display_numbers = matrix(
    ifelse(data$significance == "significant", "*", ""), 
    nrow = nrow(heatmap_matrix), 
    ncol = ncol(heatmap_matrix), 
    byrow = TRUE
  ),
  number_color = "black",
  main = "Heatmap of logFC with Clustering and Annotations"
)


library(tibble)
# Prepare the heatmap matrix
# Ensure there are no row names before using column_to_rownames
row_annotation <- data %>%
  distinct(gene, TE_type) %>%
  as.data.frame() %>%
  column_to_rownames("gene")


heatmap_matrix <- data %>%
  select(gene, celltype, logFC) %>%
  pivot_wider(names_from = celltype, values_from = logFC) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Row and column annotations
row_annotation <- data %>%
  distinct(gene, TE_type) %>%
  column_to_rownames("gene")

# Heatmap with clustering and annotations
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,  # Cluster rows
  cluster_cols = TRUE,  # Cluster columns
  annotation_row = row_annotation,  # Add row annotation
  color = colorRampPalette(c("blue", "white", "red"))(50),  # Custom color scale
  main = "Heatmap of logFC with Clustering and Annotations",
  fontsize_row = 10,
  fontsize_col = 10,
  display_numbers = matrix(ifelse(data$significance == "significant", "*", ""), nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix), byrow = TRUE),
  number_color = "black"
)

all = rbind(tab, rna)
#all$gene = gsub("-", "_", all$gene)

all$celltype_modality = paste(all$celltype, all$type)

all = all[which(!is.na(all$gene)),]

head(all[which(is.na(all$gene)),])

all[which(all$gene =="L1MA5A"& all$celltype=="CA3_Glut"),]

library(dplyr)
atac_data <- all %>% filter(type == "ATAC")
rna_data <- all %>% filter(type == "RNA")

# Merge the dataframes on the gene and celltype columns
merged_data <- merge(atac_data, rna_data, by = c("gene", "celltype"), suffixes = c("_ATAC", "_RNA"))


merged_data$TE_type = ukey[paste(merged_data$gene),6]


merged_data = merged_data[which(merged_data$fdr_ATAC<0.01 | merged_data$PValue_RNA<0.01),]


cor(merged_data$logFC_RNA,merged_data$logFC_ATAC)

options(repr.plot.width=12, repr.plot.height=5)

ggplot(merged_data, aes(x = logFC_ATAC, y = logFC_RNA, color =significance_RNA )) +
  geom_point() + xlim(c(-1,1
                       ))+ylim(c(-3,3))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +  
    geom_text_repel(max.overlaps = 100, aes(label = ifelse((logFC_ATAC > 0.25 & logFC_RNA >1.25) |  
                                                           (logFC_ATAC < -0.25 & logFC_RNA < -1.25) & 
        (fdr_ATAC<0.05 &fdr_RNA<0.05) , gene, "")),
            hjust = 1.1, vjust = 1.1, check_overlap = TRUE) +

  labs(x = "logFC (ATAC)", y = "logFC (RNA)", title = "Scatter plot of logFC for ATAC and RNA (TEs)") +
  theme_minimal()+facet_wrap(~clade_RNA)

head(all)


#sig= tab[which(tab$fdr<.01 & abs(tab$logFC) > .3 ),]
#sig_rna= rna[which(tab$fdr<.01 & abs(tab$logFC) > .3 ),]

#sig_genes = (unique(sig$gene))
#sig_genes_rna = (unique(sig_rna$gene))
#sig_genes = sig_genes_rna[which(sig_genes_rna %in% sig_genes)]
all$modality = all$type

filtered_df <- all %>%
  filter(PValue < 0.1 & abs(logFC)> 0.05 )

# Get the top 5 genes for each celltype_modality by fdr
top_genes <- filtered_df %>%
  group_by(celltype_modality, ) %>%
  arrange(fdr) %>%
  slice_head(n = 10) %>%
  ungroup()

# Print the result
nrow(top_genes)


sig_genes =unique(top_genes$gene)
sig_genes
sig_genes = c(sig_genes , "L1MA5A")
options(repr.plot.width=15, repr.plot.height=12)
# Reshape the data for the heatmap
sig_tab = all[which(all$gene %in% c(sig_genes)),]
df = sig_tab
df = df[which(df$PValue<0.01),]
#df$gene = gsub("-", "_", df$gene)
heatmap_data <- dcast(df, gene ~ celltype_modality, value.var = "logFC")
row.names(heatmap_data) <- heatmap_data$gene
heatmap_data <- heatmap_data[,-1]

df$type = ukey[paste(df$gene), "V6"]

# Annotations
gene_annotation <- df[, c("gene", "type")]
gene_annotation <- unique(gene_annotation)
row.names(gene_annotation) <- gene_annotation$gene
gene_annotation <- gene_annotation[,-1, drop = FALSE]

celltype_annotation <- df[, c("celltype_modality", "clade")]
celltype_annotation <- unique(celltype_annotation)
row.names(celltype_annotation) <- celltype_annotation$celltype
celltype_annotation <- celltype_annotation[,-1, drop = FALSE]

# Colors for annotations
ann_colors <- list(
  type = brewer.pal(8, "Set2")[1:length(unique(gene_annotation$type))],
  clade = brewer.pal(8, "Set1")[1:length(unique(celltype_annotation$clade))]
)
names(ann_colors$type) <- unique(gene_annotation$type)
names(ann_colors$clade) <- c("Glut", "NN","Gaba")



# Define the color palette for the heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
heatmap_data[which(heatmap_data>1, arr.ind = T)] = 1
heatmap_data[which(heatmap_data< -1, arr.ind = T)] = -1

# Define breaks for the heatmap so that 0 is the midpoint
max_abs_val <- max(abs(heatmap_data), na.rm = TRUE)
breaks <- seq(-max_abs_val, max_abs_val, length.out = 101)

heatmap_data = heatmap_data[,c(grep("Glut",colnames(heatmap_data)), 
                               grep("Gaba",colnames(heatmap_data)),
                                    grep("IMN",colnames(heatmap_data)),
                                         grep("NN",colnames(heatmap_data)))]

heatmap_data = heatmap_data[,c(grep("ATAC",colnames(heatmap_data)), grep("RNA",colnames(heatmap_data)))]


L1MA5A = heatmap_data[grep("L1MA5A" , rownames(heatmap_data)), ] 

heatmap_data <- heatmap_data[rowMeans(is.na(heatmap_data)) <= .9, ]

heatmap_data = heatmap_data[which(apply(abs(heatmap_data), 1, mean, na.rm = TRUE)>.1),]
heatmap_data = heatmap_data[which(apply(abs(heatmap_data), 1, sum, na.rm = TRUE)>1),]
heatmap_data = heatmap_data[which(apply(abs(heatmap_data), 1, max, na.rm = TRUE)>.1),]
heatmap_data = rbind(heatmap_data, L1MA5A) 

heatmap_data = heatmap_data[,which(apply(heatmap_data, 2, max, na.rm = TRUE)>0.1)]
heatmap_data = heatmap_data[,which(apply(heatmap_data, 2, sum, na.rm = TRUE)>.1)]
heatmap_data <- heatmap_data[,colMeans(is.na(heatmap_data)) <= .9]


cts = sapply(strsplit(as.character(colnames(heatmap_data)), " "), `[`, 1)
type = sapply(strsplit(as.character(colnames(heatmap_data)), " "), `[`, 2)


bad1 = cts[which(type == "ATAC")][which(!cts[which(type == "ATAC")]%in%cts[which(type == "RNA")])]
bad2 = cts[which(type == "RNA")][which(!cts[which(type == "RNA")]%in%cts[which(type == "ATAC")])]
bads = c(bad1,bad2)
heatmap_data = heatmap_data[,which(!cts%in%bads)]


celltype_annotation$clade = factor(celltype_annotation$clade, levels = rev(c("Gaba","Glut", "NN")))

#heatmap_data = heatmap_data[,which(colMaxs(heatmap_data)>0)]

heatmap_data = heatmap_data[order(rownames(heatmap_data)),]
#heatmap_data = heatmap_data[,order(celltype_annotation$clade)]


options(repr.plot.width=10, repr.plot.height=15)


# Create the heatmap
pheatmap(
  as.matrix(heatmap_data),gaps_col = c(grep("RNA",colnames(heatmap_data))[1]-1),
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = gene_annotation,
  annotation_col = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Age logFC - up in RNA"
)

options(repr.plot.width=12, repr.plot.height=12)

min_val <- -.75  # Adjust based on your data range for increased intensity
max_val <- .75

# Create custom breaks with a narrower range and more intervals for intensity
breaks <- seq(min_val, max_val, length.out = length(heatmap_colors) + 1)


p = pheatmap(
  t(as.matrix(heatmap_data)),gaps_row = c(grep("RNA",colnames(heatmap_data))[1]-1),
  cluster_rows = F,
  cluster_cols = F,
 # annotation_col = gene_annotation,
  annotation_row = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = T,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Age logFC"
)




