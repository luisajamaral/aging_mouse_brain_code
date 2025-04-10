library(ggplot2)
library(data.table)
library(UpSetR)
library(pheatmap)
library(Seurat)
library(enrichR)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(viridis)
library(ggrepel)
library(AnnotationDbi)


setwd("~/projects/combined_all/female_RNA/")

color_age = read.csv("../Figures/color_scheme/Age.csv")
color_region = read.csv("../Figures/color_scheme/MajorRegion-id.csv")


obj = readRDS("RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)

library(Seurat)
library(ggplot2)
options(repr.plot.width=5, repr.plot.height=5)

# 1. Subset the Seurat object for the desired cell type
dg_obj <- subset(obj, subset = ( celltype_final %in% c( "OB-STR-CTX Inh IMN") ))
dg_obj
# 2. Preprocess the data
# If your object is not already normalized, then:
#dg_obj <- NormalizeData(dg_obj,)

# Identify variable features
dg_obj <- FindVariableFeatures(dg_obj)

# Optionally, you can visualize variable features
# VariableFeaturePlot(dg_obj)

# Scale the data (all genes or a subset)
dg_obj <- ScaleData(dg_obj,vars.to.regress = c("Gapdh", "percent.mt", "percent.ribo"))

# 3. Run PCA for dimensionality reduction
dg_obj <- RunPCA(dg_obj, features = VariableFeatures(object = dg_obj))

# Optionally, check PCA results
VizDimLoadings(dg_obj, dims = 1:2, reduction = "pca")
DimPlot(dg_obj, reduction = "pca")

# 4. Determine the number of dimensions to use (here, we use the first 10)
dims_to_use <- 1:10

# 5. Find neighbors and clusters
dg_obj <- FindNeighbors(dg_obj, dims = dims_to_use)
dg_obj <- FindClusters(dg_obj, resolution = 0.5)  # Adjust resolution as needed

# 6. Run UMAP for visualization
dg_obj <- RunUMAP(dg_obj, dims = dims_to_use)

# 7. Visualize the clusters with UMAP
DimPlot(dg_obj, reduction = "umap", label = TRUE) +
  ggtitle("UMAP of dg_obj Subsampled Cells")


DimPlot(dg_obj,group.by = "celltype_final", reduction = "umap", label = TRUE) +
  ggtitle("UMAP of dg_obj Subsampled Cells")


DimPlot(dg_obj,group.by = "age", reduction = "umap", label = TRUE) +
  ggtitle("UMAP of dg_obj Subsampled Cells")


DimPlot(dg_obj,group.by = "region", reduction = "umap", label = TRUE) +
  ggtitle("UMAP of dg_obj Subsampled Cells")

DimPlot(dg_obj,group.by = "rep", reduction = "umap", label = TRUE) +
  ggtitle("UMAP of dg_obj Subsampled Cells")

FeaturePlot(dg_obj,"Pax6", reduction = "umap")
FeaturePlot(dg_obj,"Lhx2", reduction = "umap")
FeaturePlot(dg_obj,"Apoe", reduction = "umap")

VlnPlot(dg_obj,"Pax6", group.by = "age", split.by = "region")


qc_genes <- c("1700030F04Rik", "1810037I17Rik", "2010107E04Rik", "Acadl", "Ap2s1",
              "Atp5d", "Atp5e", "Atp5g1", "Atp5g3", "Atp5j", "Atp5j2", "Atp5k", 
              "Atp5l", "Atp6v1f", "Chchd10", "Coa3", "Cox5b", "Cox6a1", "Cox6b1", 
              "Cox6c", "Cox7a2", "Cox7a2l", "Cox7b", "Cycs", "Edf1", "Eef1b2", 
              "Eif5a", "Fau", "Fkbp3", "Ftl1", "Guk1", "Heph", "Hras", "Mif", 
              "Mrfap1", "Naca", "Ndufa1", "Ndufa2", "Ndufa4", "Ndufa5", "Ndufb7", 
              "Ndufc1", "Ndufs7", "Ndufv3", "Necap1", "Nlrc4", "Pdxp", "Pfn2", 
              "Polr2m", "Rab3a", "Rtl8a", "Slc16a2", "Snrpd2", "Snu13", "Taf1c", 
              "Timm8b", "Tpt1", "Ubb", "Uqcr11", "Uqcrb", "Uqcrq", "Usp50")

# Ensure genes exist in Seurat object
qc_genes <- qc_genes[qc_genes %in% rownames(obj)]

# Add module score to Seurat object
obj <- AddModuleScore(obj, features = list(qc_genes), name = "QC_Score")

# Plot module score across samples
ggplot(obj@meta.data, aes(x = sample, y = QC_Score1)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "QC Module Score Across Samples", x = "Sample", y = "QC Module Score")

options(repr.plot.width=8, repr.plot.height=15)

# Plot module score across samples
ggplot(obj@meta.data, aes(x = region, y = QC_Score1, fill = age)) +
  geom_boxplot() +
  theme_minimal() + coord_flip() +
  labs(title = "QC Module Score Across Samples", x = "Sample", y = "QC Module Score")

seurat_obj = obj
# Extract expression matrix and log-transform it
expr_matrix <- GetAssayData(seurat_obj, slot = "data")[qc_genes, , drop = FALSE]  # This retrieves log-normalized values

# Sum log-transformed expression for each cell to compute QC Score
seurat_obj$QC_Score <- colSums(expr_matrix)

# Fetch Malat1 expression (log-normalized)
if ("Malat1" %in% rownames(seurat_obj)) {
    seurat_obj$Malat1_Expression <- FetchData(seurat_obj, vars = "Malat1")
} else {
    warning("Malat1 is not in the Seurat object. Skipping Malat1 comparison.")
}

# Visualize QC Score across samples
ggplot(seurat_obj@meta.data, aes(x = sample, y = QC_Score)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "QC Score Across Samples", x = "Sample", y = "QC Score")

# Check correlation with Malat1 if available
if ("Malat1_Expression" %in% colnames(seurat_obj@meta.data)) {
    qc_correlation <- cor.test(seurat_obj$QC_Score, seurat_obj$Malat1_Expression, method = "pearson")
    print(qc_correlation)
}

library(ggplot2)

# Histogram of QC Scores
ggplot(seurat_obj@meta.data, aes(x = QC_Score)) +
  geom_histogram(bins = 100, fill = "blue", alpha = 0.5) +
  theme_minimal() + coord_flip()+
  labs(title = "Distribution of QC Scores", x = "QC Score", y = "Cell Count")

# Violin plot across samples
ggplot(seurat_obj@meta.data, aes(x = sample, y = QC_Score)) +
  geom_violin(fill = "lightblue", alpha = 0.5) +
  theme_minimal() + coord_flip()+
  labs(title = "QC Score Across Samples", x = "Sample", y = "QC Score")



qc_threshold <- quantile(seurat_obj$QC_Score, probs = 0.95)  # Bottom 5% cutoff
bad_cells <- rownames(seurat_obj@meta.data[seurat_obj$QC_Score > qc_threshold, ])

# Number of bad cells
length(bad_cells)

library(dplyr)
library(ggplot2)
library(tidyr)


# Count total cells per sample before filtering
total_cells_per_sample <- seurat_obj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(total_cells = n())

# Count removed cells per sample
removed_cells_per_sample <- seurat_obj@meta.data[bad_cells, ] %>%
  group_by(orig.ident) %>%
  summarise(removed_cells = n())

# Merge data and compute percentage lost
loss_summary <- left_join(total_cells_per_sample, removed_cells_per_sample, by = "orig.ident") %>%
  mutate(removed_cells = replace_na(removed_cells, 0),  # Replace NA with 0 for samples with no removed cells
         percent_lost = (removed_cells / total_cells) * 100)

# View data
print(loss_summary)


head(loss_summary,100)

ggplot(seurat_obj@meta.data, aes(x = orig.ident, y = QC_Score, fill = age)) +
  geom_boxplot() +
  theme_minimal() + coord_flip()+
  labs(title = "QC Score Across Samples", x = "Sample", y = "QC Score")


FeaturePlot(seurat_obj,"QC_Score", reduction = "umap", label = F) 

Idents(dg_obj) ="celltype_final"
fm= FindAllMarkers(dg_obj , only.pos = T)

head(fm[which(fm$cluster == "OB-STR-CTX Inh IMN"),],15)

FeaturePlot(dg_obj,"Aldh1l1", reduction = "umap", label = TRUE) 
DimPlot(dg_obj,group.by = "sample", reduction = "umap", label = TRUE) 
DimPlot(dg_obj,group.by = "age", reduction = "umap", label = TRUE) 
DimPlot(dg_obj,group.by = "celltype_final", reduction = "umap", label = TRUE) 

all_genes = rownames(obj)
head(all_genes)
write(all_genes , file = "../Figu

Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)
colnames(av$RNA)

obj$celltype_age = paste(obj$celltype_final,  obj$age)

Idents(obj) = "celltype_age"
av = AverageExpression(obj)

obj$celltype_final[which(obj$best_celltype_fixed == "IOL")] = "IOL NN"
obj$celltype_age = paste(obj$celltype_final,  obj$age)

Idents(obj) = "celltype_age"

obj$celltype_age = paste(obj$celltype_final,  obj$age)

Idents(obj) = "celltype_age"

fm = FindMarkers(test.use = "MAST",latent.vars = "orig.ident",obj, `ident.1` = "IOL NN 18mo", `ident.2` = "IOL NN 2mo", logfc.threshold = 0.01, min.pct = 0.01)

head(fm[order(fm$p_val_adj),],20)

head(fm[order(fm$p_val_adj),],20)

color_age = read.csv("../Figures/color_scheme/Age.csv")
color_region = read.csv("../Figures/color_scheme/MajorRegion-id.csv")

head(color_region)
head(color_age)

colnames(av)



deg_results = read.csv("DEG_results_.01_.01/IOL_NN.csv")

deg_results= fm



deg_results[1:5,]

head(df)
df[which(df$X=="Pax6"),]

options(repr.plot.width=7, repr.plot.height=4)
#deg_results = read.csv("progenitor_DEGs/OB-STR-CTX_Inh_IMN2vsOld.csv")
deg_results = read.csv("DEG_results_rm_rpl_rps/OB-STR-CTX_Inh_IMN.csv")

deg_results$avg_log2FC = -deg_results$avg_log2FC
# Add classification
df <- deg_results %>%
  mutate(
    significance = case_when(
      avg_log2FC > 0.25 & p_val_adj < 0.05 ~ "Upregulated",
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

df$pct_exp = (df$`pct.1`+  df$`pct.2`)/2
head(df)
# Select top 10 differential genes
top_up_genes <- df[which(df$avg_log2FC>.5 & df$p_val_adj<1e-3)[1:7], ]
top_down_genes <- df[which(df$avg_log2FC< (-.5) & df$p_val_adj<1e-3)[1:7], ]
top_down_genes <- rbind(top_down_genes,df[which(df$X=="Pax6"),])

top_genes <- rbind(top_down_genes,top_up_genes,df[which(df$X=="Pax6"),])

head(top_genes)
#top_genes = top_genes[-grep("Rik", top_genes$X),]
#top_genes = top_genes[-grep("Gm", top_genes$X),]

# MA Plot
ma_plot <- ggplot(df, aes(x = (pct_exp), y = avg_log2FC)) + # Log-transform base_mean
  geom_point(aes(color = significance), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "OB-STR-CTX Inh IMN Age-DEGs\nMA Plot",
    x = "% cells expression",
    y = "Log2 Fold Change (M)",
    color = "Significance"
  ) +
  theme(legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+
  geom_text_repel(nudge_y = 1,
    data = top_up_genes,
    aes(label = X),
    size = 5.25,
    max.overlaps = 100,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
geom_text_repel(nudge_y = -1,
    data = top_down_genes,
    aes(label = X),
    size = 5.5,
    max.overlaps = 100,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  )


print(ma_plot)

options(repr.plot.width=5, repr.plot.height=7.5)

ma_plot

pdf("../Figures/Figure2-IMN-IOL/OB_MA_Plot.pdf", width = 5, height = 7.5)
ma_plot
dev.off()

options(repr.plot.width=6, repr.plot.height=4)

#clust$Gene[which(clust$Cluster==1)]
c_up = enrichGO(universe = df$X,gene = df$X[which(df$significance=="Upregulated"
                                                 )],ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.1)
c_up@result$order= 1:nrow(c_up@result)


sc_up <- simplify(
                          c_up, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )
options(repr.plot.width=6, repr.plot.height=5)



d1 = dotplot(sc_up,showCategory = 5, title = "Age-Up genes\nGO BP Enrichment")
d1

options(repr.plot.width=6, repr.plot.height=3)
#clust$Gene[which(clust$Cluster==1)]
c_down = enrichGO(universe = df$X,gene = df$X[which(df$significance=="Downregulated"
                                                 )],ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.01)
c_down@result$order= 1:nrow(c_down@result)

sc_down <- simplify(
                          c_down, 
                          cutoff = 0.6,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )
d2 = dotplot(sc_down,showCategory = 7, title = "Age-Down genes\nGO BP Enrichment")
d2

options(repr.plot.width=6, repr.plot.height=5)
d1 = dotplot(sc_up,showCategory = 8, title = "Age-Up genes\nGO BP Enrichment")
d2 = dotplot(sc_down,showCategory = 8, title = "Age-Down genes\nGO BP Enrichment")

d1
d2


pdf("../Figures/Figure2-IMN-IOL/up_GO_BP.pdf", width = 6, height = 5)
d1
dev.off()


pdf("../Figures/Figure2-IMN-IOL/down_GO_BP.pdf", width = 6, height = 5)
d2
dev.off()


library(org.Mm.eg.db)
library(GO.db)
library(AnnotationDbi)

# Get all descendant GO terms of "learning or memory" (GO:0007611)
all_go_terms <- c("GO:0007611", GOBPOFFSPRING[["GO:0007611"]])

# Retrieve all genes mapped to the GO term and its children
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = all_go_terms, 
  keytype = "GO", 
  columns = c("SYMBOL")
)

# Get unique gene symbols
learning_or_memory <- unique(genes_in_go$SYMBOL)
length(learning_or_memory)  # Should now match enrichGO




library(org.Mm.eg.db)
library(GO.db)
library(AnnotationDbi)

# Get all descendant GO terms of "learning or memory" (GO:0007611)
all_go_terms <- c("GO:0050808", GOBPOFFSPRING[["GO:0050808"]])

# Retrieve all genes mapped to the GO term and its children
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = all_go_terms, 
  keytype = "GO", 
  columns = c("SYMBOL")
)

# Get unique gene symbols
synapse_organization <- unique(genes_in_go$SYMBOL)
length(synapse_organization)  # Should now match enrichGO




deg_results = read.csv("DEG_results_rm_rpl_rps/OB-STR-CTX_Inh_IMN.csv")

deg_results$avg_log2FC = -deg_results$avg_log2FC

deg_results[grep("Lhx", deg_results$X),]

deg_results = top_genes
deg_results = read.csv("DEG_results_rm_rpl_rps/OB-STR-CTX_Inh_IMN.csv")

deg_results$avg_log2FC = -deg_results$avg_log2FC
deg_results <- deg_results %>%
  filter(p_val_adj < 0.005
         & abs(avg_log2FC) > 0.5) # Adjust thresholds as needed
nrow(deg_results)
head(deg_results)
head(df)
deg_results = rbind(deg_results , df[which(df$X=="Pax6"),])
head(deg_results)

deg_results$pathway = NA
deg_results$pathway[which(deg_results$X %in% learning_or_memory)] = "learning or memory"
#deg_results$pathway[which(deg_results$X %in% synapse_organization)] = "synapse organization"


# Get all unique DEGs from your results

all_degs <- unique(deg_results$X)
deg_expression = av$RNA[all_degs,grep("OB-STR", colnames(av$RNA))]
deg_expression = deg_expression[,c(grep("NAC", colnames(deg_expression)),grep(" CP", colnames(deg_expression)))]

#deg_expression = deg_expression[,c(grep("NAC", colnames(deg_expression)),grep("AMY", colnames(deg_expression)), grep("CP", colnames(deg_expression)))]
#deg_expression = deg_expression[,c(grep("HCA", colnames(deg_expression)),grep("HCP", colnames(deg_expression)))]

# Normalize data for better visualization
deg_expression_scaled <- t(scale(t(as.matrix(deg_expression)))) # Scale by gene
deg_expression_scaled[is.na(deg_expression_scaled)] <- 0 # Handle NA values

(colSums(deg_expression)) 

dist_matrix <- dist(deg_expression_scaled)
hclust_result <- hclust(dist_matrix, method = "ward.D2")
num_clusters <- 4 # Choose the number of clusters manually or use a method like the elbow plot
clusters <- cutree(hclust_result, k = num_clusters)

set.seed(123)
kmeans_result <- kmeans(deg_expression_scaled, centers = 4, nstart = 25) # Adjust centers
clusters <- kmeans_result$cluster


#metadata <- data.frame(
#  original_order = colnames(deg_expression_scaled),
#  region = sub(".*Glut (.+?) \\d+mo", "\\1", colnames(deg_expression_scaled)),
#  age = as.numeric(sub(".* (\\d+)mo", "\\1", colnames(deg_expression_scaled)))
#)

metadata <- data.frame(
  original_order = colnames(deg_expression_scaled),
  region = sapply(strsplit(as.character(colnames(deg_expression_scaled)), " "), `[`, 4),
  age = sapply(strsplit(as.character(colnames(deg_expression_scaled)), " "), `[`, 5)
)

metadata$age = factor(metadata$age , levels = c("2mo", "9mo", "18mo"))

# Order by region alphabetically, then by age numerically
metadata <- metadata[order(metadata$region, metadata$age), ]
rownames(metadata)  = metadata$original_order
#metadata$age = paste(metadata$age , "mo" ,sep ="")


deg_expression_scaled <- deg_expression_scaled[, metadata$original_order]

# Confirm the new order
colnames(deg_expression_scaled)

metadata

options(repr.plot.width=5, repr.plot.height=6)
library(khroma)
PRGn <- color("PRGn")
BuRd <- color("BuRd")

num_clusters=2
# Convert to named color vectors
age_colors <- setNames(color_age$Color, color_age$Age)
region_colors <- setNames(color_region$Color, color_region$Region)

# Combine color annotations into a list
annotation_colors <- list(
  age = age_colors,
  region = region_colors
)

row_annotation = deg_results[,"pathway", drop = F]
rownames(row_annotation) = deg_results$X
metadata$age = as.factor(metadata$age)

#deg_expression_scaled = deg_expression_scaled[,c(1,4,2,5,3,6)]
# Generate the heatmap
p <- pheatmap(gaps_col = c(2,4),
  deg_expression_scaled[,],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,main  = "Top OB-STR-CTX Inh IMN DEGs",
  cutree_rows = num_clusters, color = viridis(100),
  #color = colorRampPalette(c("blue", "white", "red"))(100)  ,
  show_rownames = T,
  show_colnames = TRUE,
  #main = paste("DG Age-DEGs\n", table(deg_results$direction)[[2]], " up | ",table(deg_results$direction)[[1]]," down", sep = "") ,
  annotation_row = row_annotation,  # Row annotation
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)

head(deg_results[which(deg_results$p_val_adj<0.001),]

met = metadata[, c("region", "age")]
rownames(met) =  gsub("OB-STR-CTX Inh IMN " , "", rownames(met) )


annotation_colors

barplot(deg_expression[which(rownames(deg_expression_scaled)=="Pax6") ,])
deg_expression[which(rownames(deg_expression_scaled)=="Pax6") ,]

library(ggplot2)
library(dplyr)

# Create the dataset
pax6_expression <- data.frame(
  Region = rep(c("NAC", "CP"), each = 3),
  Age = rep(c("2mo", "9mo", "18mo"), times = 2),
  Expression = c(1.6328, 0.8357, 0.0370, 0.9717, 0.5464, 0.0585)  # Pax6 expression values
)


options(repr.plot.width=4, repr.plot.height=5)

# Ensure Age is an ordered factor
pax6_expression$Age <- factor(pax6_expression$Age, levels = c("2mo", "9mo", "18mo"))
# Plot the data as a grouped barplot
ggplot(pax6_expression, aes(x = Age, y = Expression, fill = Region)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.5) +  
  geom_point(position = position_dodge(width = 0.6), size = 2, color = "black") +  # Add points
  scale_fill_manual(values = c("NAC" = "#F3C8E8", "CP" = "#a56bd4")) +  # Custom colors
  theme_minimal(base_size = 16) +
  labs(title = "Pax6 Expression in OB-STR-CTX Inh IMN",
       x = "Age",
       y = "Average Normalized Expression",
       fill = "Brain Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


obj$Pax6 = obj@assays$RNA$scale.data["Pax6",]


obj$Mal = obj@assays$RNA$scale.data["Mal",]


library(purrr)  # Load the purrr package for map functions
library(broom)  # broom is useful for tidying up statistical test results
library(ggpubr) 
options(repr.plot.width=4.5, repr.plot.height=4)
subset_dg <- obj@meta.data %>%
  filter(celltype_final == "Oligo NN")



gene = "Mal" 



# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_dg, aes(x = age, y = Mal, fill = age
                           )) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "t.test", label = "p.format") +
  # Customize the theme
  theme_classic() +#ylim(c(0,3))+
  labs(
    x = "Age Group",
    y = paste("Expression of ",gene),
    fill = "Age Group",
    title = paste("Oligo NN",gene," (RNA)")
  ) +
 # scale_color_manual(values = c("Rep1" = "black", "Rep2" = "black","Rep3" = "black")) +

  scale_fill_manual(values = c("2mo" = '#480080', "9mo" = '#e23c5d', "18mo" = '#ffb42c')) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )#+
 # scale_fill_manual(values = c('NAC' = '#F3C8E8', 'CP'= '#a56bd4', 'HCA'= '#e8a57d','HCP'='#6f2632')) 
g

rownames(obj)

pdf(g, file = "../Figures/Figure3-oligos/Mal_RNA.pdf",height = 4, width = 5)
g
dev.off()

library(purrr)  # Load the purrr package for map functions
library(broom)  # broom is useful for tidying up statistical test results
library(ggpubr) 
options(repr.plot.width=4.5, repr.plot.height=4)
subset_dg <- obj@meta.data %>%
  filter(celltype_final == "OB-STR-CTX Inh IMN" & region %in% c('NAC', 'CP'))



gene = "Pax6" 

# Create a violin plot with boxplot-like statistics and statistical comparisons
g = ggplot(subset_dg, aes(x = age, y = Pax6, fill = age
                           )) +
  # Add violin plot
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  # Add boxplot overlay for summary statistics (median, IQR)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), alpha = 0.6, outlier.shape = NA) +
  # Add mean as a point for extra clarity
  stat_summary(fun = "mean", geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  # Add statistical comparisons (3 tests between age groups)
  stat_compare_means(comparisons = list(c("2mo", "9mo"), c("2mo", "18mo"), c("9mo", "18mo")), 
                     method = "t.test", label = "p.format") +
  # Customize the theme
  theme_classic() +#ylim(c(0,3))+
  labs(
    x = "Age Group",
    y = paste("Expression of ",gene),
    fill = "Age Group",
    title = paste("OB-STR-CTX Inh IMN ",gene," (RNA)")
  ) +
 # scale_color_manual(values = c("Rep1" = "black", "Rep2" = "black","Rep3" = "black")) +

  scale_fill_manual(values = c("2mo" = '#480080', "9mo" = '#e23c5d', "18mo" = '#ffb42c')) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16)
  )#+
 # scale_fill_manual(values = c('NAC' = '#F3C8E8', 'CP'= '#a56bd4', 'HCA'= '#e8a57d','HCP'='#6f2632')) 
g

pdf(g, file = "../Figures/Figure2-IMN-IOL/Pax6_RNA.pdf",height = 4, width = 5)
g
dev.off()

options(repr.plot.width=4, repr.plot.height=4)
colnames(deg_expression_scaled) = gsub("OB-STR-CTX Inh IMN " , "", colnames(deg_expression_scaled))
p <- pheatmap(gaps_row = c(2,4),cutree_cols = 2,
 # (deg_expression_scaled[which(row_annotation$pathway!="") ,]),
   (deg_expression_scaled[which(rownames(deg_expression_scaled)=="Pax6") ,]),

  clustering_method = "ward.D2",
  cluster_rows = T,
  cluster_cols = F,main  = "OB-STR-CTX Inh IMN DEGs\nlearning or memory",
  cutree_rows = num_clusters, color = viridis(100),
  #color = colorRampPalette(c("blue", "white", "red"))(100)  ,
  show_rownames = T,
  show_colnames = TRUE,
  #main = paste("DG Age-DEGs\n", table(deg_results$direction)[[2]], " up | ",table(deg_results$direction)[[1]]," down", sep = "") ,
 # annotation_col = row_annotation,  # Row annotation
  annotation_col = met,  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)

pdf(file = "../Figures/Figure2-IMN-IOL/OB_IMN_DEG_heat_horiz.pdf", height =4, width = 4)
print(p)
dev.off()

options(repr.plot.width=5, repr.plot.height=3)
colnames(deg_expression_scaled) = gsub("OB-STR-CTX Inh IMN " , "", colnames(deg_expression_scaled))
p <- pheatmap(gaps_row = c(2,4),cutree_cols = 2,
  t(deg_expression_scaled[which(row_annotation$pathway!="") ,]),
  clustering_method = "ward.D2",
  cluster_rows = F,
  cluster_cols = T,main  = "OB-STR-CTX Inh IMN DEGs\nlearning or memory",
  cutree_rows = num_clusters, color = viridis(100),
  #color = colorRampPalette(c("blue", "white", "red"))(100)  ,
  show_rownames = T,
  show_colnames = TRUE,
  #main = paste("DG Age-DEGs\n", table(deg_results$direction)[[2]], " up | ",table(deg_results$direction)[[1]]," down", sep = "") ,
 # annotation_col = row_annotation,  # Row annotation
  annotation_row = met,  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)

pdf(file = "../Figures/Figure2-IMN-IOL/OB_IMN_DEG_heat_horiz.pdf", height =3, width = 5)
print(p)
dev.off()

deg_expression_scaled["Pax6",]

deg_results[which(deg_results$X=="Pax6"),]

pdf(file = "../Figures/Figure2-IMN-IOL/OB_IMN_DEG_heat.pdf", height =8, width = 6)
print(p)
dev.off()

options(repr.plot.width=4, repr.plot.height=7)

# Perform hierarchical clustering
row_clustering <- hclust(dist(deg_expression_scaled), method = "ward.D2")

# Cut tree into the desired number of clusters
cluster_assignments <- cutree(row_clustering, k = num_clusters)

# Create a data frame mapping genes to clusters
gene_clusters <- data.frame(Gene = rownames(deg_expression_scaled), Cluster = cluster_assignments)

# Print or save the cluster assignments
#print(gene_clusters)
#write.csv(gene_clusters, "gene_clusters.csv", row.names = FALSE)

# Create a row annotation data frame
row_annotation <- data.frame(Cluster = factor(cluster_assignments))
rownames(row_annotation) <- rownames(deg_expression_scaled)

# Define colors for clusters
cluster_colors <- RColorBrewer::brewer.pal(num_clusters, "Set2")
names(cluster_colors) <- levels(row_annotation$Cluster)
annotation_colors$Cluster <- cluster_colors

# Replot the heatmap with cluster annotation
p <- pheatmap(gaps_col = c(3),
  deg_expression_scaled,
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  breaks = seq(-3, 3, length.out = 100),
  color = BuRd(100),
  show_rownames = F,
  show_colnames = TRUE,
   annotation_row = row_annotation,  # Row annotation with clusters
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)


which(deg_results$X=="Tnks")
nrow(deg_results)

clust1 = gene_clusters$Gene[which(gene_clusters$Cluster==1
                                 )] 
options(repr.plot.width=6, repr.plot.height=6)

#clust$Gene[which(clust$Cluster==1)]
c_1 = enrichGO(universe = df$X,gene = clust1,ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.7)
c_1@result$order= 1:nrow(c_1@result)
sc_1 <- simplify(
                          c_1, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_1,showCategory = 10, title = "Clust1 genes\nGO CC Enrichment")
d1

clust1 = gene_clusters$Gene[which(gene_clusters$Cluster==2
                                 )] 
options(repr.plot.width=6, repr.plot.height=6)

#clust$Gene[which(clust$Cluster==1)]
c_1 = enrichGO(universe = df$X,gene = clust1,ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.05)
c_1@result$order= 1:nrow(c_1@result)
sc_2 <- simplify(
                          c_1, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_2,showCategory = 10, title = "Clust2 genes\nGO CC Enrichment")
d1

clust1 = gene_clusters$Gene[which(gene_clusters$Cluster==3
                                 )] 
options(repr.plot.width=6, repr.plot.height=6)

#clust$Gene[which(clust$Cluster==1)]
c_1 = enrichGO(universe = df$X,gene = clust1,ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.05)
c_1@result$order= 1:nrow(c_1@result)
sc_3 <- simplify(
                          c_1, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_3,showCategory = 10, title = "Clust3 genes\nGO CC Enrichment")
d1

options(repr.plot.width=7, repr.plot.height=7)

# Perform hierarchical clustering
row_clustering <- hclust(dist(deg_expression_scaled), method = "ward.D2")

# Cut tree into the desired number of clusters
cluster_assignments <- cutree(row_clustering, k = num_clusters)

# Create a data frame mapping genes to clusters
gene_clusters <- data.frame(Gene = rownames(deg_expression_scaled), Cluster = cluster_assignments)

# Print or save the cluster assignments
#print(gene_clusters)
#write.csv(gene_clusters, "gene_clusters.csv", row.names = FALSE)

# Create a row annotation data frame
row_annotation <- data.frame(Cluster = factor(cluster_assignments))
rownames(row_annotation) <- rownames(deg_expression_scaled)

# Define colors for clusters
cluster_colors <- RColorBrewer::brewer.pal(num_clusters, "Set2")
names(cluster_colors) <- levels(row_annotation$Cluster)
annotation_colors$Cluster <- cluster_colors

# Replot the heatmap with cluster annotation
p <- pheatmap(gaps_col = c(6),
  deg_expression_scaled,
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  breaks = seq(-3, 3, length.out = 100),
  color = BuRd(100),
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = paste("DG Age-DEGs with Clusters\n",
               table(deg_results$direction)[[2]], " up | ",
               table(deg_results$direction)[[1]], " down", sep = ""),
  annotation_row = row_annotation,  # Row annotation with clusters
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)



