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

obj = readRDS("RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)

obj$age_binary = "old"
obj$age_binary[which(obj$age == "2mo")] = "young"
obj$celltype_age_binary = paste(obj$celltype_final,  obj$age_binary)
Idents(obj) = "celltype_age_binary"

fm = FindMarkers(test.use = "MAST",latent.vars = "region",obj, `ident.1` = "IOL NN old", `ident.2` = "IOL NN young", logfc.threshold = 0.01, min.pct = 0.01)

head(fm,20)

obj$celltype_final[which(obj$best_celltype_fixed == "IOL")] = "IOL NN"
obj$celltype_age = paste(obj$celltype_final,  obj$age)
Idents(obj) = "celltype_age"

obj$celltype_age = paste(obj$celltype_final,  obj$age)
Idents(obj) = "celltype_age"
sub = obj[which(obj$region %in% c("HCA", "HCP"))]
fm = FindMarkers(test.use = "MAST",latent.vars = "region",
                 sub, `ident.1` = "DG Glut 18mo", 
                 `ident.2` = "DG Glut 2mo", 
                 logfc.threshold = 0.01, min.pct = 0.01)

fm = FindMarkers(test.use = "MAST",latent.vars = "region",
                 sub, `ident.1` = "Oligo NN 18mo", 
                 `ident.2` = "Oligo NN 2mo", 
                 logfc.threshold = 0.01, min.pct = 0.01)

write.csv(fm, file = "../Figures/Figure3-oligos/Oligo_DEG_MAST_latent_region_.01_.01.csv")

sub = obj[which(obj$region %in% c("HCA", "HCP"))]

fm = FindMarkers(test.use = "MAST",latent.vars = "region",sub, `ident.1` = "DG Glut 18mo", `ident.2` = "DG Glut 2mo", logfc.threshold = 0.01, min.pct = 0.01)

write.csv(fm, file = "../Figures/Figure5-DG/DG_DEG_MAST_latent_region_.01_.01.csv")

color_age = read.csv("../Figures/color_scheme/Age.csv")
color_region = read.csv("../Figures/color_scheme/MajorRegion-id.csv")

head(color_region)
head(color_age)

deg_results = read.csv("DEG_results_rm_rpl_rps/L2-3_IT_CTX_Glut.csv")

head(deg_results)



deg_results$avg_log2FC = -deg_results$avg_log2FC

deg_results= fm
deg_results$X  =rownames(fm)

options(repr.plot.width=6, repr.plot.height=5)

# Add classification
df <- deg_results %>%
  mutate(
    significance = case_when(
      avg_log2FC > 0.25 & p_val_adj < 1e-15 ~ "Upregulated",
      avg_log2FC < -0.25 & p_val_adj < 1e-15 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )
cat(table(df$significance))
df$pct_exp = (df$`pct.1`+  df$`pct.2`)/2
head(df)
# Select top 10 differential genes
top_genes <- df[c(1:20), ]
#top_genes = top_genes[-grep("Rik", top_genes$X),]
#top_genes = top_genes[-grep("Gm", top_genes$X),]

# MA Plot
ma_plot <- ggplot(df, aes(x = (pct_exp), y = avg_log2FC)) + # Log-transform base_mean
  geom_point(aes(color = significance), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Oligodendrocytes MA Plot",
    x = "% cells expression",
    y = "Log2 Fold Change (M)",
    color = "Significance"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+
  geom_text_repel(
    data = top_genes,
    aes(label = X),
    size = 4,
    max.overlaps = 100,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  )


print(ma_plot)

head(df,30)

ma_plot +
  theme(legend.position = "top")+ ggtitle("Oligodendrocyte Age-DEGs MA Plot")

pdf("../Figures/Figure3-oligos/MA.pdf", height =7)
ma_plot+theme(legend.position = "top")
dev.off()

svg("../Figures/Figure3-oligos/MA.svg", height =6)
ma_plot
dev.off()

head(df[grep("Pcdh", df$X),],20)

head(df[which(df$significance=="Upregulated"),],10)

options(repr.plot.width=6, repr.plot.height=7)

#clust$Gene[which(clust$Cluster==1)]
c_up = enrichGO(universe = df$X,gene = df$X[which(df$significance=="Upregulated"
                                                 )],ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.5)
c_up@result$order= 1:nrow(c_up@result)
sc_up <- simplify(
                          c_up, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_up,showCategory = 10, title = "Age-Up genes\nGO CC Enrichment")
d1

svg("../Figures/Figure3-oligos/up_GO_CC.svg", height =3.5, width = 6)
d1
dev.off()

pdf("../Figures/Figure3-oligos/up_GO_CC.pdf", height =3.5, width = 6)
d1
dev.off()

options(repr.plot.width=6, repr.plot.height=6)

#clust$Gene[which(clust$Cluster==1)]
c_down = enrichGO(universe = df$X,gene = df$X[which(df$significance=="Downregulated"
                                                 )],ont = "CC",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.05)
c_down@result$order= 1:nrow(c_down@result)

sc_down <- simplify(
                          c_down, 
                          cutoff = 0.6,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )
d2 = dotplot(sc_down,showCategory = 10, title = "Age-Down genes\nGO CC Enrichment")
d2

svg("../Figures/Figure3-oligos/down_GO_CC.svg", height =3.5, width = 6)
d2
dev.off()

pdf("../Figures/Figure3-oligos/down_GO_CC.pdf", height =3.5, width = 6)
d2
dev.off()

d1 
d2

head(sc_down,10)

# synaptic membrane GO:0097060
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0097060", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

synaptic_membrane = unique(genes_in_go$SYMBOL)
head(synaptic_membrane)

# myelin sheath GO:0043209
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0043209", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

myelin_sheath = unique(genes_in_go$SYMBOL)
head(myelin_sheath)

# cytosolic ribosome GO:0022626
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0022626", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

cytosolic_ribosome = unique(genes_in_go$SYMBOL)
head(cytosolic_ribosome)

# mitochondrial respirasome GO:0005746
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0005746", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

mitochondrial_respirasome = unique(genes_in_go$SYMBOL)
head(mitochondrial_respirasome)


# presynapse GO:0098793
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0098793", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

presynapse = unique(genes_in_go$SYMBOL)
head(presynapse)

deg_results <- deg_results %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) # Adjust thresholds as needed

# Get all unique DEGs from your results

all_degs <- unique(deg_results$X)
deg_expression = av$RNA[all_degs,grep("DG Glut", colnames(av$RNA))]

#deg_expression = deg_expression[,c(grep("NAC", colnames(deg_expression)),grep("AMY", colnames(deg_expression)), grep("CP", colnames(deg_expression)))]
deg_expression = deg_expression[,c(grep("HCA", colnames(deg_expression)),grep("HCP", colnames(deg_expression)))]

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


metadata <- data.frame(
  original_order = colnames(deg_expression_scaled),
  region = sub(".*NN (.+?) \\d+mo", "\\1", colnames(deg_expression_scaled)),
  age = as.numeric(sub(".* (\\d+)mo", "\\1", colnames(deg_expression_scaled)))
)

# Order by region alphabetically, then by age numerically
metadata <- metadata[order(metadata$age, metadata$region), ]
rownames(metadata)  = metadata$original_order
metadata$age = paste(metadata$age , "mo" ,sep ="")


deg_expression_scaled <- deg_expression_scaled[, metadata$original_order]

# Confirm the new order
colnames(deg_expression_scaled)

row_annotation$pathway= ""
row_annotation$pathway[which(rownames(row_annotation)%in% myelin_vector)] = "myelin sheath"
row_annotation = row_annotation[,2, drop = F]

which(deg_results$X %in% myelin_sheath)

deg_results$pathway = NA
deg_results$pathway[which(deg_results$X %in% myelin_sheath)] = "myelin sheath"
deg_results$pathway[which(deg_results$X %in% cytosolic_ribosome)] = "cytosolic ribosome"
deg_results$pathway[which(deg_results$X %in% synaptic_membrane)] = "synaptic membrane"
deg_results$pathway[which(deg_results$X %in% mitochondrial_respirasome)] = "mitochondrial respirasome"
deg_results$pathway[which(deg_results$X %in% presynapse)] = "presynapse"

options(repr.plot.width=7, repr.plot.height=7)

# Create a data frame for the row annotation
row_annotation <- data.frame(
  Pathway = as.factor(deg_results$pathway)  # Convert to factor for better coloring
)
rownames(row_annotation) <- deg_results$X  # Add gene names as rownames
# Replot the heatmap with the row annotation
p <- pheatmap(
  deg_expression_scaled[ ,],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  show_rownames = T,
  show_colnames = TRUE,
  main = nrow(deg_expression_scaled),
  annotation_row = row_annotation,  # Add the cluster annotation as a row annotation
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  fontsize = 10
)


rownames(deg_results) = deg_results$X
deg_results$direction = "up"
deg_results$direction[which(deg_results$avg_log2FC>0)] = "down"
head(deg_results)


table(deg_results$direction)

row_annotation = deg_results[, "direction", drop = F]
head(row_annotation)

print(table(deg_results$direction))

library(khroma)
PRGn <- color("PRGn")
BuRd <- color("BuRd")


options(repr.plot.width=7, repr.plot.height=7)


# Convert to named color vectors
age_colors <- setNames(color_age$Color, color_age$Age)
region_colors <- setNames(color_region$Color, color_region$Region)

# Combine color annotations into a list
annotation_colors <- list(
  age = age_colors,
  region = region_colors
)
metadata$age = as.factor(metadata$age)
# Generate the heatmap
p <- pheatmap(
  deg_expression_scaled[,],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  breaks = seq(-3,3,length.out=100), color = PRGn(100),

  #color = colorRampPalette(c("blue", "white", "red"))(100)  ,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = paste("Oligodendrocyte Age-DEGs\n", table(deg_results$direction)[[2]], " up | ",table(deg_results$direction)[[1]]," down", sep = "") ,
#  annotation_row = row_annotation,  # Row annotation
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)

sub = deg_expression_scaled[-grep("Rik" , rownames(deg_expression_scaled)),]
sub = sub[-grep("Gm" , rownames(sub)),]
#sub = sub[grep("Mal", rownames(sub)),]

deg_results$pathway = NA
deg_results$pathway[which(deg_results$X %in% myelin_sheath)] = "myelin sheath"
#deg_results$pathway[which(deg_results$X %in% cytosolic_ribosome)] = "cytosolic ribosome"
#deg_results$pathway[which(deg_results$X %in% synaptic_membrane)] = "synaptic membrane"
#deg_results$pathway[which(deg_results$X %in% mitochondrial_respirasome)] = "mitochondrial respirasome"
#deg_results$pathway[which(deg_results$X %in% presynapse)] = "presynapse"

options(repr.plot.width=7, repr.plot.height=6)
# Create a data frame for the row annotation
row_annotation <- data.frame(
  Pathway = as.factor(deg_results$pathway)  # Convert to factor for better coloring
)
rownames(row_annotation) = deg_results$X
p <- pheatmap(gaps_col = c(8,16),
  sub[1:30,],
  clustering_method = "ward.D2",
  cluster_rows = T,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  breaks = seq(-3.5,3.5,length.out=25), color = viridis(25),

  #color = colorRampPalette(c("blue", "white", "red"))(100)  ,
  show_rownames = T,
  show_colnames = TRUE,
  main = paste("Oligodendrocyte top 50 Age-DEGs", sep = "") ,
  annotation_row = row_annotation,  # Row annotation
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)

svg("../../combined_all/Figures/Figure3-oligos/Oligo_DEG_top30_heatmap_myelin.svg", height = 6)
p
dev.off()

while(1)dev.off()

genes = read.csv("DEG_results_.01_.01/Oligo_NN.csv")
genes$avg_log2FC = -genes$avg_log2FC
genes = genes[order(genes$avg_log2FC,decreasing = T),]
genesList = genes$avg_log2FC
key = bitr( as.character(genes$X), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db,drop = F)
genes$ENTREZID = key$ENTREZID
names(genesList) = as.character(genes$ENTREZID)
head(genesList)


pdf("../Figures/Figure3-oligos/Oligo_DEG_heatmap2.pdf", width = 6)
p
dev.off()

svg("../../Figures/Figure3-oligos/Oligo_DEG_heatmap.svg", width = 6)
p
dev.off()

svg("Oligo_DEG_heatmap.svg", height = 7 ,width = 7)
print(p)
dev.off()


#clust$Gene[which(clust$Cluster==1)]
c1 = enrichGO(gene = deg_results$X[which(deg_results$direction=="up"
                                                 )],ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.01)

sc1 <- simplify(
                          c1, 
                          cutoff = 0.6,               # Similarity cutoff
                          by = "qvalue",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )


dotplot(sc1,showCategory = 15, title = "Age-Up genes\nGO BP Enrichment")


options(repr.plot.width=7, repr.plot.height=5)

d = dotplot(sc1,showCategory = 11, title = "Age-Up genes\nGO BP Enrichment")
d

svg("../../Figures/Figure3-oligos/up_go_bp.svg", height =5)
d
dev.off()

options(repr.plot.width=7, repr.plot.height=7)

#clust$Gene[which(clust$Cluster==1)]
c1 = enrichGO(gene = deg_results$X[which(deg_results$direction=="down"
                                                 )],ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.01)

sc1 <- simplify(
                          c1, 
                          cutoff = 0.65,               # Similarity cutoff
                          by = "qvalue",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )
grep("enshea", sc1$Description)
dotplot(sc1,showCategory = 20, title = "Age-Down genes\nGO BP Enrichment")




#head(sc1, 50)

sc_custom = sc1[which(sc1@result$ID %in% c("GO:0002181", "GO:0042063",
                                         "GO:0006119","GO:2001233",
                                         "GO:0050807","GO:0006644","GO:0016311","GO:0051402" , 
                                         "GO:0009060", "GO:0046034","GO:0006457", "GO:0007272",
                                         "GO:0008203","GO:0008366",   "GO:0034976"))
                    ,]
sc1@result = sc_custom                    
d=dotplot(sc1,showCategory = 20, title = "Age-Down genes\nGO BP Enrichment")
d

svg("../../Figures/Figure3-oligos/down_go_bp.svg", height =5)
d
dev.off()

options(repr.plot.width=7, repr.plot.height=5)

dotplot(sc1,showCategory = 11, title = "Age-Down genes\nGO BP Enrichment")



# Create a data frame with gene names and their cluster labels
cluster_result_df <- data.frame(
  Gene = rownames(deg_expression_scaled),
  Cluster = cluster_labels
)

# View or save the result
head(cluster_result_df)

options(repr.plot.width=5, repr.plot.height=6)

#clust$Gene[which(clust$Cluster==1)]
c1 = enrichGO(gene = cluster_result_df$Gene[which(cluster_result_df$Cluster==1
                                                 )],ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.01)

sc1 <- simplify(
                          c1, 
                          cutoff = 0.62,               # Similarity cutoff
                          by = "qvalue",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

dotplot(c1)
dotplot(sc1,showCategory = 22, title = "Cluster 1")


grep("myel" ,sc1$Description)


pdf("cluster_result_deg_D12MSN.pdf")
print(p)
dev.off()

write.table(cluster_result_df , "cluster_result_df_D12MSN.txt")

head(cluster_result_df[which(cluster_result_df$Cluster == 1),])

head(cluster_result_df[which(cluster_result_df$Cluster == 3),])


