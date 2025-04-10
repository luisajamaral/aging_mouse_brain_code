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
library(biomaRt)


setwd("~/projects/combined_all/female_RNA/")
obj = readRDS("RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)

color_age = read.csv("../Figures/color_scheme/Age.csv")
color_region = read.csv("../Figures/color_scheme/MajorRegion-id.csv")

head(color_region)
head(color_age)



color_region  = color_region[c(),]

dg_obj = readRDS("../Figures//Figure5-DG/dg_obj.RDS")

DimPlot(dg_obj,group.by = "age") 
DimPlot(dg_obj,group.by = "region") 
DimPlot(dg_obj,group.by = "sample",label=T) 

dg_obj$weird = "no"
dg_obj$weird[which(dg_obj$sample %in% c("HCP_9mo_3","HCP_9mo_2", "HCA_9mo_2", "HCA_2mo_2"))]= "weird"
DimPlot(dg_obj,group.by = "weird",label=T) 

DimPlot(dg_obj,group.by = "age",label=T) 

options(repr.plot.width=5, repr.plot.height=5)

FeaturePlot(dg_obj, "Rpl13")

options(repr.plot.width=18, repr.plot.height=5)
VlnPlot(dg_obj, "Rpl13", group.by = "sample")
VlnPlot(dg_obj, "Malat1", group.by = "sample")

norm_data <- GetAssayData(dg_obj, slot = "data")
# Check that MALAT1 is present
if (!"Malat1" %in% rownames(norm_data)) {
  stop("MALAT1 is not found in the normalized data. Please check the gene name.")
}

# Extract MALAT1 expression across all cells
malat1_exp <- norm_data["Malat1", ]
# Compute correlation for each gene with MALAT1
gene_correlations <- apply(norm_data, 1, function(gene_exp) {
  cor(gene_exp, malat1_exp, method = "pearson")
})

# Sort correlations in descending order (most positive correlation first)
cor_sorted <- sort(gene_correlations, decreasing = TRUE)

# Top 10 genes most positively correlated with MALAT1
top_positive <- head(cor_sorted, 10)
print("Top 10 genes positively correlated with MALAT1:")
print(top_positive)

# Top 10 genes most negatively correlated with MALAT1
top_negative <- tail(cor_sorted, 10)
print("Top 10 genes negatively correlated with MALAT1:")
print(top_negative)

# Top 10 genes most positively correlated with MALAT1
top_positive <- head(cor_sorted, 10)
print("Top 10 genes positively correlated with MALAT1:")
print(top_positive)

Idents(dg_obj) = "weird"
fm = FindAllMarkers(dg_obj[,which(dg_obj$age =="9mo")],only.pos = T)
head(fm)

head(fm[which(fm$cluster=="weird"),],50)

library(Seurat)
library(ggplot2)
options(repr.plot.width=5, repr.plot.height=5)

# 1. Subset the Seurat object for the desired cell type
#dg_obj <- subset(obj, subset = ( celltype_final =="DG Glut" & region %in% c("HCP", "HCA")))

# 2. Preprocess the data
# If your object is not already normalized, then:
dg_obj <- NormalizeData(dg_obj,)

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


options(repr.plot.width=10, repr.plot.height=5)

DimPlot(dg_obj , group.by = "age")
DimPlot(dg_obj , split.by = "region", group.by = "sample")

#DimPlot(dg_obj , group.by = "weird")
#FeaturePlot(dg_obj , "Gapdh")
FeaturePlot(dg_obj , "Malat1")

DimPlot(dg_obj , group.by = "sample")


obj$celltype_region_age_rep = paste(obj$celltype_final, obj$region, obj$age, obj$rep)
Idents(obj) = "celltype_region_age_rep"
av2 = AverageExpression(obj)


# Connect to Ensembl
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # Change to "hsapiens_gene_ensembl" for human

# Retrieve gene types for all genes in your dataset
gene_annotations <- getBM(attributes = c("external_gene_name", "gene_biotype","transcript_length"),
                          filters = "external_gene_name",
                          values = rownames(obj),
                          mart = mart)

# Filter for lncRNAs
lncRNA_genes <- gene_annotations$external_gene_name[gene_annotations$gene_biotype == "lncRNA"]

protein_coding_genes <- gene_annotations$external_gene_name[gene_annotations$gene_biotype == "protein_coding"]

ggplot(gene_annotations, aes(x = gene_biotype, y = transcript_length)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +  # Adds individual points for better visualization
  scale_y_log10() +  # Log scale to better visualize differences
  theme_minimal() +coord_flip()+
  labs(title = "Transcript Length by Gene Biotype",
       x = "Gene Biotype",
       y = "Transcript Length (log10)")

all_genes <- rownames(obj)

# Identify ribosomal genes (RPL and RPS)
ribosomal_genes <- all_genes[grepl("^Rpl|^Rps", all_genes)]

# Print the ribosomal genes for verification
cat("Ribosomal genes to be removed:\n")
print(ribosomal_genes)


genes_to_remove <- unique(c(ribosomal_genes))
obj
# Remove genes from the Seurat object
dg_obj <- subset(dg_obj, features = setdiff(rownames(dg_obj), genes_to_remove))
dg_obj

dg_obj = obj[,which(obj$celltype_final == "DG Glut")]
genes_to_remove <- unique(c(ribosomal_genes))
dg_obj
# Remove genes from the Seurat object
dg_obj <- subset(dg_obj, features = setdiff(rownames(dg_obj), genes_to_remove))
dg_obj
dg_obj$celltype_age = paste(dg_obj$celltype_final,  dg_obj$age)
Idents(dg_obj) = "celltype_age"
sub = dg_obj
fm = FindMarkers(latent.vars = c( "region", "rep"),test.use = "MAST",
                 obj, `ident.1` = "DG Glut 18mo", 
                 `ident.2` = "DG Glut 2mo", 
                 logfc.threshold = 0.1, min.pct = 0.01)
head(fm) 
fm$X = rownames(fm)
deg_results = fm

head(deg_results)

obj$celltype_age = paste(obj$celltype_final,  obj$age)
Idents(obj) = "celltype_age"
sub = obj[which(obj$celltype_final == "DG Glut" & obj$region %in% c("HCA", "HCP"))]
fm = FindMarkers(test.use = "MAST",latent.vars = "region",
                 sub, `ident.1` = "DG Glut 18mo", 
                 `ident.2` = "DG Glut 2mo", 
                 logfc.threshold = 0.25, min.pct = 0.01)
head(fm) 
fm$X = rownames(fm)
deg_results = fm



deg_results = read.csv("../Figures/Figure5-DG/DG_DEG_MAST_latent_region_.01_.01.csv")

deg_results = read.csv("DEG_results_.01_.01/DG_Glut.csv")
deg_results$avg_log2FC = -deg_results$avg_log2FC

deg_results = read.csv("DEG_results_latent_rep_mito_together//L2-3_IT_ENT_Glut.csv")
deg_results$avg_log2FC = -deg_results$avg_log2FC

deg_results = read.csv("DEG_results_region/DG_Glut--HCP.csv")
deg_results$avg_log2FC = -deg_results$avg_log2FC

options(repr.plot.width=5, repr.plot.height=5)

deg_results = read.csv("DEG_results_rm_rpl_rps_mrpl_mt/DG_Glut.csv")
#deg_results = read.csv("../Figures/Figure5-DG/DG_DEG_MAST_latent_region_.01_.01.csv")
deg_results$avg_log2FC = -deg_results$avg_log2FC
deg_results$annotation = "Other"
deg_results$annotation[which(deg_results$X%in%lncRNA_genes)] = "lncRNA"
deg_results$annotation[which(deg_results$X%in%protein_coding_genes)] = "protein_coding"
deg_results <- deg_results %>%
  mutate(
    significance = case_when(
      avg_log2FC > 0.25 & p_val_adj < 0.05 & annotation == "Other" ~ "Upregulated Other",
      avg_log2FC < -0.25 & p_val_adj < 0.05 & annotation == "Other" ~ "Downregulated Other",
      avg_log2FC > 0.25 & p_val_adj < 0.05 & annotation == "protein_coding" ~ "Upregulated protein_coding",
      avg_log2FC < -0.25 & p_val_adj < 0.05 & annotation == "protein_coding" ~ "Downregulated protein_coding",
      avg_log2FC > 0.25 & p_val_adj < 0.05 & annotation == "lncRNA" ~ "Upregulated lncRNA",
      avg_log2FC < -0.25 & p_val_adj < 0.05  & annotation == "lncRNA" ~ "Downregulated lncRNA",
      TRUE ~ "Not Significant"
    )
  )
# Example dataset (replace this with your actual DEG data)
deg_summary <- deg_results %>%
  arrange(desc(avg_log2FC)) %>%  # Sort genes by logFC (high to low)
  mutate(rank = row_number())  # Assign a rank to each gene

# Plot
ggplot(deg_summary, aes(x = rank, y = avg_log2FC, fill = annotation)) +
  geom_bar(stat = "identity", width = 1) +  # Thin bars to emphasize ranking
  #facet_wrap(~celltype, scales = "free_x") +  # Facet by cell type
  theme_minimal() +
  scale_fill_manual(values = c("lncRNA" = "red", "protein_coding" = "blue", "Other" = "gray")) +
  labs(
    title = "Ranked Gene Expression DG_Glut",
    x = "Ranked Genes (sorted by logFC)",
    y = "log2 Fold Change",
    fill = "Gene Annotation"
  ) +
  theme(
    axis.text.x = element_blank(),  # Hide x-axis labels for clarity
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )
#deg_results=deg_results[!which(deg_results$annotation == "lncRNA"),]

library(ggplot2)
library(dplyr)
library(tidyr)
# If available, consider readr::read_csv for faster file reading.

# Define directory containing DEG results
deg_dir <- "DEG_results_rm_rpl_rps_mrpl_mt/"

# Get all CSV files in the directory
deg_files <- list.files(deg_dir, full.names = TRUE, pattern = "*.csv")

# Placeholder for results
deg_summary <- data.frame()

# Loop through each DEG file
for (file in deg_files) {
  
  # Extract cell type name from file name
  celltype <- gsub(".csv", "", basename(file))
  
  # Read data
  deg_results <- read.csv(file)
  
  # Reverse the sign of avg_log2FC if needed
  deg_results$avg_log2FC <- -deg_results$avg_log2FC
  
  # Annotate gene types
  deg_results$annotation <- "Other"
  deg_results$annotation[deg_results$X %in% lncRNA_genes] <- "lncRNA"
  deg_results$annotation[deg_results$X %in% protein_coding_genes] <- "protein_coding"
  
  # Filter significant genes (p_val_adj < 0.05) to decide if there is enough signal
  sig_degs <- deg_results %>% filter(p_val_adj < 0.1)
  
  # Only proceed if there are at least 100 significant DEGs
  if (nrow(sig_degs) < 50) {
    next
  }
  
  # Select top 100 upregulated and top 100 downregulated genes from the significant ones
  top_up <- sig_degs %>% 
    arrange(desc(avg_log2FC)) %>% 
    head(100) %>% 
    mutate(category = "Upregulated")
  
  top_down <- sig_degs %>% 
    arrange(avg_log2FC) %>% 
    head(100) %>% 
    mutate(category = "Downregulated")
  
  # Instead of filtering only for non-significant genes, select all genes that are not in the top lists.
  # This will include genes that are significant (but not in the top 100) and non-significant genes.
  all_others <- deg_results %>% 
    filter(!(X %in% c(top_up$X, top_down$X))) %>% 
    mutate(category = "All Genes")
  
  # Combine the three groups: top up, top down, and all other genes.
  combined_genes <- bind_rows(top_up, top_down, all_others)
  
  # Count genes in each category per annotation type
  count_summary <- combined_genes %>%
    group_by(celltype = celltype, category, annotation) %>%
    summarise(count = n(), .groups = "drop")
  
  # Append to overall summary
  deg_summary <- bind_rows(deg_summary, count_summary)
}

# Assign clade labels based on celltype patterns
deg_summary$clade <- "Glut"
deg_summary$clade[grep("NN", deg_summary$celltype)] <- "NN"
deg_summary$clade[grep("IMN", deg_summary$celltype)] <- "NN"
deg_summary$clade[grep("Gaba", deg_summary$celltype)] <- "Gaba"

# Plot gene counts per cell type
ggplot(deg_summary, aes(x = celltype, y = count, fill = annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(clade ~ category, scales = "free", space = "free_y") +
  theme_minimal() +
  labs(title = "Top 100 Age-DEGs versus All Genes",
       x = "Cell Type",
       y = "Gene Proportion",
       fill = "Gene Annotation") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(ggplot2)
library(dplyr)
library(tidyr)
options(repr.plot.width=8, repr.plot.height=12)

# Define directory containing DEG results
deg_dir <- "DEG_results_rm_rpl_rps_mrpl_mt/DG-PIR_Ex_IMN.csv"

# Get all CSV files in the directory
deg_files <- list.files(deg_dir, full.names = TRUE, pattern = "*.csv")



# Placeholder for results
deg_counts <- data.frame()

# Loop through each DEG file
for (file in deg_files) {
  
  # Extract cell type name from file name
  celltype <- gsub(".csv", "", basename(file))

  # Read data
  deg_results <- read.csv(file)
  deg_results$avg_log2FC = -deg_results$avg_log2FC  # Reverse if needed
  
  # Count upregulated and downregulated genes
  up_count <- sum(deg_results$avg_log2FC > 0.25 & deg_results$p_val_adj < 0.05, na.rm = TRUE)
  down_count <- sum(deg_results$avg_log2FC < -0.25 & deg_results$p_val_adj < 0.05, na.rm = TRUE)

  # Store results
  temp_df <- data.frame(
    celltype = c(celltype, celltype),
    direction = c("Upregulated", "Downregulated"),
    count = c(up_count, -down_count)  # Negative count for mirror effect
  )

  deg_counts <- bind_rows(deg_counts, temp_df)
}

# Assign clades based on cell type naming pattern
deg_counts$clade <- "Glut"
deg_counts$clade[grep("NN", deg_counts$celltype)] <- "NN"
deg_counts$clade[grep("Gaba", deg_counts$celltype)] <- "Gaba"
deg_counts$clade[grep("IMN", deg_counts$celltype)] <- "NN"

# Plot mirror bar plot with facet by clade
ggplot(deg_counts, aes(x = celltype, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +  # Flip to make it horizontal
  facet_wrap(~clade, scales = "free_y",ncol = 1) +  # Facet by clade
  theme_minimal() +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(
    title = "Significantly Upregulated and Downregulated Genes per Cell Type",
    x = "Cell Type",
    y = "Number of DEGs",
    fill = "Direction"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


head(deg_counts)

library(ggplot2)
library(dplyr)
options(repr.plot.width=8, repr.plot.height=12)


# Define directory containing DEG results
deg_dir <- "DEG_results_latent_rep_mito///////"

# Get all CSV files in the directory
deg_files <- list.files(deg_dir, full.names = TRUE, pattern = "*.csv")

# Placeholder for results
logFC_summary <- data.frame()

# Gene of interest (Change this to any other gene)
gene_of_interest <- "Meg3"

# Loop through each DEG file
for (file in deg_files) {
  
  # Extract cell type name from file name
  celltype <- gsub(".csv", "", basename(file))

  # Read data
  deg_results <- read.csv(file)
  deg_results$avg_log2FC = -deg_results$avg_log2FC  # Reverse if needed
  deg_results$p_val_adj = as.numeric(deg_results$p_val_adj)
  # Extract logFC and p-value for the gene of interest
  gene_data <- deg_results %>%
    filter(X == gene_of_interest) %>%
    dplyr::select(X, avg_log2FC, p_val_adj) %>%
    mutate(celltype = celltype)

  # If the gene is found, append to summary
  if (nrow(gene_data) > 0) {
    logFC_summary <- bind_rows(logFC_summary, gene_data)
  }
}

# Convert adjusted p-values into significance categories
logFC_summary <- logFC_summary %>%
  mutate(significance = case_when(
    p_val_adj < 0.01 & avg_log2FC < 0 ~ "Downregulated",
    p_val_adj < 0.01 & avg_log2FC > 0 ~ "Upregulated",
    TRUE ~ "NS"
  ))

logFC_summary$clade =  "Glut"
logFC_summary$clade[grep("NN",logFC_summary$celltype)] =  "NN"
logFC_summary$clade[grep("Gaba",logFC_summary$celltype)] =  "Gaba"
logFC_summary$clade[grep("IMN",logFC_summary$celltype)] =  "NN"

logFC_summary$region = sapply(strsplit(as.character(logFC_summary$celltype), "--"), `[`, 2)


# Plot logFC across cell types
ggplot(logFC_summary, aes(x = celltype, y = avg_log2FC, color = significance, size = -log10(p_val_adj+1e-50))) +
  geom_point() +
  facet_grid( clade ~ region , scales = "free_y", space = "free") +  

  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Add a reference line at logFC = 0
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "NS" = "gray")) +
  theme_minimal() +coord_flip()+ylim(c(-3,3))+
  labs(
    title = paste("LogFC of", gene_of_interest, "Across Cell Types"),
    x = "Cell Type",
    y = "Log Fold Change",
    color = "Significance",
    size = "-log10(p_val_adj)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


grep("Pcd",deg_results$X)
deg_results$X[134]

library(ggplot2)
library(dplyr)
options(repr.plot.width=8, repr.plot.height=12)


# Define directory containing DEG results
deg_dir <- "../../merfish//"

# Get all CSV files in the directory
deg_files <- list.files(deg_dir, full.names = TRUE, pattern = "*.csv")

# Placeholder for results
logFC_summary <- data.frame()

# Gene of interest (Change this to any other gene)
gene_of_interest <- "Pcdhb7"

# Loop through each DEG file
for (file in deg_files) {
  
  # Extract cell type name from file name
  celltype <- gsub(".csv", "", basename(file))

  # Read data
  deg_results <- read.csv(file)
  deg_results$avg_log2FC = -deg_results$avg_log2FC  # Reverse if needed
  deg_results$p_val_adj = as.numeric(deg_results$p_val_adj)
  # Extract logFC and p-value for the gene of interest
  gene_data <- deg_results %>%
    filter(X == gene_of_interest) %>%
    dplyr::select(X, avg_log2FC, p_val_adj) %>%
    mutate(celltype = celltype)

  # If the gene is found, append to summary
  if (nrow(gene_data) > 0) {
    logFC_summary <- bind_rows(logFC_summary, gene_data)
  }
}

# Convert adjusted p-values into significance categories
logFC_summary <- logFC_summary %>%
  mutate(significance = case_when(
    p_val_adj < 0.01 & avg_log2FC < 0 ~ "Downregulated",
    p_val_adj < 0.01 & avg_log2FC > 0 ~ "Upregulated",
    TRUE ~ "NS"
  ))

logFC_summary$clade =  "Glut"
logFC_summary$clade[grep("NN",logFC_summary$celltype)] =  "NN"
logFC_summary$clade[grep("Gaba",logFC_summary$celltype)] =  "Gaba"
logFC_summary$clade[grep("IMN",logFC_summary$celltype)] =  "NN"

# Plot logFC across cell types
ggplot(logFC_summary, aes(x = celltype, y = avg_log2FC, color = significance, size = -log10(p_val_adj+1e-50))) +
  geom_point() +
  facet_grid( clade ~ 0 , scales = "free_y", space = "free") +  

  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +  # Add a reference line at logFC = 0
  scale_color_manual(values = c("Downregulated" = "blue", "Upregulated" = "red", "NS" = "gray")) +
  theme_minimal() +coord_flip()+ylim(c(-3,3))+
  labs(
    title = paste("LogFC of", gene_of_interest, "Across Cell Types (MERFISH)"),
    x = "Cell Type",
    y = "Log Fold Change",
    color = "Significance",
    size = "-log10(p_val_adj)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability


#deg_results = deg_results[which(is.na(deg_results$annotation)),]
head(deg_results)
table(deg_results$annotation)

options(repr.plot.width=6.5, repr.plot.height=5.5)
#deg_results = fm
deg_results = read.csv("///STR_D12_Gaba.csv")
deg_results$avg_log2FC = -deg_results$avg_log2FC 

# Add classification
df <- deg_results %>%
  mutate(
    significance = case_when(
      avg_log2FC > 0.25 & p_val_adj < 1e-15 ~ "Upregulated",
      avg_log2FC < -0.25 & p_val_adj < 1e-15 ~ "Downregulated",
       TRUE ~ "Not Significant"
    )
  )
table(df$significance)
df$pct_exp = (df$`pct.1`+  df$`pct.2`)/2
head(df)
# Select top 10 differential genes
top_genes <- rbind(df[which(df$avg_log2FC>0.5)[1:10], ],df[which(df$avg_log2FC< -0.5)[1:10], ])
#top_genes = top_genes[-grep("Rik", top_genes$X),]
#top_genes = top_genes[-grep("Gm", top_genes$X),]

# MA Plot
ma_plot <- ggplot(df, aes(x = (pct_exp), y = avg_log2FC)) + # Log-transform base_mean
  geom_point(aes(color = significance), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated" = "red","Upregulated lncRNA"="magenta","Upregulated Other" = "darkorange",
                                "Downregulated Other" = "deepskyblue",
                                "Downregulated lncRNA"="purple",
                                "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +

  labs(
    title = "STR D12 Gaba MA Plot",
    x = "% cells expression",
    y = "Log2 Fold Change (18mo/2mo)",
    color = "Significance"
  ) +ylim(c(-4,4))+
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


# Add classification
df <- deg_results %>%
  mutate(
    pathway = case_when(
      X %in% cell_adhesion ~ "cell_adhesion",
      TRUE ~ " "
    )
  )
table(df$significance)
df$pct_exp = (df$`pct.1`+  df$`pct.2`)/2
head(df)
# Select top 10 differential genes
top_genes <- df[which(df$avg_log2FC>0.5 & df$pathway == "DNA damage response" )[1:10], ]
#top_genes = top_genes[-grep("Rik", top_genes$X),]
#top_genes = top_genes[-grep("Gm", top_genes$X),]
df = df[order(df$pathway),]
# MA Plot
ma_plot2 <- ggplot(df, aes(x = (pct_exp), y = avg_log2FC)) + # Log-transform base_mean
  geom_point(aes(color = pathway), size = 2, alpha = 0.7) +
    theme_minimal() +
  labs(
    x = "% cells expression",
    y = "Log2 Fold Change (18mo/2mo)",
    color = "Pathway"
  ) +ylim(c(-4,4))+
  scale_color_manual(values = c(' ' = "gray", "DNA damage response" = "hotpink")) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
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

print(ma_plot2)

head(geneList)

geneList <- deg_results$avg_log2FC
names(geneList) <- deg_results$X
geneList <- na.omit(geneList)

# Sort the gene list in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# Convert gene symbols to Entrez IDs using bitr() from clusterProfiler.
# Change the OrgDb if you're working with mouse data.
gene.df <- bitr(names(geneList),
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)

# Filter geneList to include only genes with an Entrez ID
geneList_entrez <- geneList[names(geneList) %in% gene.df$SYMBOL]
# Replace the gene symbols with Entrez IDs:
geneList_entrez <- geneList_entrez[match(gene.df$SYMBOL, names(geneList_entrez))]
names(geneList_entrez) <- gene.df$ENTREZID

# Optionally, remove any NA values (if any conversion failed)
geneList_entrez <- na.omit(geneList_entrez)

# Now run GSEA with the KEGG database using gseKEGG()
# Change organism = "hsa" to "mmu" if you are analyzing mouse data.
kegg_gsea <- gseKEGG(geneList = geneList_entrez,
                     organism = "mmu",
                     nPerm = 1000,
                     minGSSize = 10,
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

# View the top pathways
head(kegg_gsea@result)

dotplot(kegg_gsea)


# Convert results to a data frame
kegg_df <- as.data.frame(kegg_result)

# Sort by enrichment score (-log10 p-value)
kegg_df$logp <- -log10(kegg_df$p.adjust)

# Identify upregulated and downregulated pathways
upregulated <- kegg_df[order(-kegg_df$logp), ][1:5, ]  # Top 5 upregulated
downregulated <- kegg_df[order(kegg_df$logp), ][1:5, ]  # Top 5 downregulated

# Combine both sets
top_kegg <- rbind(upregulated, downregulated)

# Dotplot of top pathways
dotplot_kegg <- ggplot(top_kegg, aes(x = logp, y = reorder(Description, logp), color = logp, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "Top 5 Up and Downregulated KEGG Pathways",
    x = "-log10(Adjusted P-value)",
    y = "KEGG Pathway"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Print plot
print(dotplot_kegg)


library(cowplot)
options(repr.plot.width=11, repr.plot.height=4)

# Combine both plots into one column
combined_plot <- plot_grid(ma_plot, ma_plot2, ncol = 1, align = "v" )

# Display the combined plot
print(combined_plot)




options(repr.plot.width=6, repr.plot.height=5)

#clust$Gene[which(clust$Cluster==1)]
c_up = enrichGO(universe = df$X, gene = df$X[which(df$p_val_adj<1e-15 & df$avg_log2FC>0.25)],ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.1)
c_up@result$order= 1:nrow(c_up@result)
sc_up <- simplify(
                          c_up, 
                          cutoff = 0.8,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_up,showCategory = 10, title = "Age-Up genes\nGO BP Enrichment")
d1

head(c_up,20)

d1 = dotplot(sc_up,showCategory = 10, title = "Age-Up genes\nGO BP Enrichment")
d1

pdf( file = "~/projects/combined_all/Figures/Figure4-Gaba/up_GO_BP.pdf", height = 5, width = 6.5)
d1
dev.off()

d1


#clust$Gene[which(clust$Cluster==1)]
c_down = enrichGO(universe = df$X, gene = df$X[which(df$p_val_adj<1e-5 & df$avg_log2FC< (-0.25))],ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.1)
c_down@result$order= 1:nrow(c_down@result)
sc_down <- simplify(
                          c_down, 
                          cutoff = 0.8,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_down,showCategory = 13, title = "Age-Down genes\nGO BP Enrichment")
pdf(d1, file = "~/projects/combined_all/Figures/Figure4-Gaba/down_GO_BP.pdf", height = 5, width = 6.5)
d1
dev.off()

d1

#DNA damage response GO:0006974
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0006974", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

DNA_damage_response = unique(genes_in_go$SYMBOL)
head(DNA_damage_response)

#cell adhesion GO:0098609
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0098609", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

cell_adhesion = unique(genes_in_go$SYMBOL)
head(cell_adhesion)

cell_adhesion

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

# myelin sheath GO:0043209
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0043209", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

myelin_sheath = unique(genes_in_go$SYMBOL)
head(myelin_sheath)


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
  filter(p_val_adj < 1e-25 & abs(avg_log2FC) > 0.25) # Adjust thresholds as needed

nrow(deg_results)

## Get all unique DEGs from your results

all_degs <- unique(deg_results$X)
deg_expression = av$RNA[all_degs,grep("DG Glut", colnames(av$RNA))]

#deg_expression = deg_expression[,c(grep("NAC", colnames(deg_expression)), grep("CP", colnames(deg_expression)))]
deg_expression = deg_expression[,c(grep("HCA", colnames(deg_expression)),grep("HCP", colnames(deg_expression)))]

#deg_expression = deg_expression[,c(grep("t FC ", colnames(deg_expression)),grep("t ENT ", colnames(deg_expression)))]
#deg_expression = deg_expression[,c(grep("t FC ", colnames(deg_expression)),grep("09d09 ", colnames(deg_expression)))]

# Normalize data for better visualization
deg_expression_scaled <- t(scale(t(as.matrix(deg_expression)))) # Scale by gene
deg_expression_scaled[is.na(deg_expression_scaled)] <- 0 # Handle NA values

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
  region = sapply(strsplit(as.character(colnames(deg_expression_scaled)), " "), `[`, 3),
  age = sapply(strsplit(as.character(colnames(deg_expression_scaled)), " "), `[`, 4)
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

options(repr.plot.width=7, repr.plot.height=10)

# Create a data frame for the row annotation
row_annotation <- data.frame(
  Pathway = as.factor(deg_results$pathway)  # Convert to factor for better coloring
)
#rownames(row_annotation) <- deg_results$X  # Add gene names as rownames
# Replot the heatmap with the row annotation
p <- pheatmap(
  deg_expression_scaled[1:1000
                        ,],
   # row_annotation = row_annotation,
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = T,
  cutree_rows = num_clusters,
  show_rownames = F,
  show_colnames = TRUE,
  main = nrow(deg_expression_scaled),
  fontsize = 10
)


deg_results$pathway = NA
#deg_results$pathway[which(deg_results$X %in% myelin_sheath)] = "myelin sheath"
#deg_results$pathway[which(deg_results$X %in% cytosolic_ribosome)] = "cytosolic ribosome"
#deg_results$pathway[which(deg_results$X %in% mitochondrial_respirasome)] = "mitochondrial respirasome"
#deg_results$pathway[which(deg_results$X %in% presynapse)] = "presynapse"
deg_results$pathway[which(deg_results$X %in% DNA_damage_response)] = "DNA damage response"

#deg_results$pathway = NA
#deg_results$pathway[which(deg_results$X == "Gapdh")] = "Gapdh"


which(deg_results$X == "Gapdh")

options(repr.plot.width=7, repr.plot.height=7)

# Create a data frame for the row annotation
row_annotation <- data.frame(
  Pathway = as.factor(deg_results$pathway)  # Convert to factor for better coloring
)
rownames(row_annotation) <- deg_results$X  # Add gene names as rownames
# Replot the heatmap with the row annotation
p <- pheatmap(gaps_col = c(3), 
  #deg_expression_scaled[which(rownames(deg_expression_scaled)%in%DNA_damage_response) ,-c(9)],
  deg_expression_scaled[c(1:150
                          
                         ) ,-c(9)],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  show_rownames = F,
  show_colnames = TRUE,
  main = "top 150 DEGs",
  annotation_row = row_annotation,  # Add the cluster annotation as a row annotation
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  fontsize = 10
)


##### rownames(deg_results) = deg_results$X
deg_results$direction = "up"
deg_results$direction[which(deg_results$avg_log2FC<0)] = "down"
head(deg_results,20)


row_annotation = deg_results[, "direction", drop = F]


library(khroma)
PRGn <- color("PRGn")
BuRd <- color("BuRd")


options(repr.plot.width=7, repr.plot.height=7)

num_clusters=4
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
p <- pheatmap(gaps_col = c(6),
  deg_expression_scaled[1:1000,-c(9)],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  breaks = seq(-3,3,length.out=100), color = BuRd(100),
  #color = colorRampPalette(c("blue", "white", "red"))(100)  ,
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = paste("DG Age-DEGs\n", table(deg_results$direction)[[2]], " up | ",table(deg_results$direction)[[1]]," down", sep = "") ,
  annotation_row = row_annotation,  # Row annotation
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)

table(obj$region[which(obj$celltype_final=="DG Glut")],obj$age[which(obj$celltype_final=="DG Glut")])

#deg_expression_scaled = deg_expression_scaled[,-c(9)]

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


anno = deg_results[,"annotation", drop = F]
rownames(anno) = deg_results$X

p <- pheatmap(gaps_col = c(6), 
  deg_expression_scaled[,],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  #breaks = seq(-3, 3, length.out = 100),
  color = BuRd(100),
  show_rownames = F,
  show_colnames = TRUE,
  annotation_row = anno,  # Row annotation with clusters
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  #annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)


p <- pheatmap(gaps_col = c(6), 
  deg_expression_scaled[,-c(9)],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  breaks = seq(-3, 3, length.out = 100),
  color = BuRd(100),
  show_rownames = FALSE,
  show_colnames = TRUE,
  main = paste("DG Glut Age-DEGs with Clusters\n",
               table(deg_results$direction)[[2]], " up | ",
               table(deg_results$direction)[[1]], " down", sep = ""),
  annotation_row = row_annotation,  # Row annotation with clusters
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)


table(row_annotation$Pathway)

options(repr.plot.width=6, repr.plot.height=5)

# Create a data frame for the row annotation
row_annotation <- data.frame(
  Pathway = as.factor(deg_results$pathway)  # Convert to factor for better coloring
)
rownames(row_annotation) <- deg_results$X  # Add gene names as rownames
# Replot the heatmap with the row annotation
p <- pheatmap(gaps_col = c(3),   color = viridis(100),

  #deg_expression_scaled[which(rownames(deg_expression_scaled)%in%DNA_damage_response) ,-c(9)],
  deg_expression_scaled[c(1:1000) ,],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = 2,  breaks = seq(-2, 2, length.out = 100),

  show_rownames = F,
  show_colnames = TRUE,
  main = "top 1000 DEGs",  annotation_colors = annotation_colors,  # Add the color schemes

  annotation_row = row_annotation,  # Add the cluster annotation as a row annotation
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  fontsize = 10
)


pdf( file = "~/projects/combined_all/Figures/Figure5-DG/heat.pdf", height = 5, width = 7)
p
dev.off()

while(1) dev.off()

options(repr.plot.width=10, repr.plot.height=3)

p <- pheatmap(gaps_col = c(3),   color = viridis(100),

  #deg_expression_scaled[which(rownames(deg_expression_scaled)%in%DNA_damage_response) ,-c(9)],
  t(deg_expression_scaled[c(1:50
                          
                         ) ,]),
  clustering_method = "ward.D2",
  cluster_rows = F,
  cluster_cols = T,
  cutree_rows = 2,  breaks = seq(-2, 2, length.out = 100),

  show_rownames = T,
  show_colnames = TRUE,
  main = "top 100 DEGs",  annotation_colors = annotation_colors,  # Add the color schemes

  annotation_col = row_annotation,  # Add the cluster annotation as a row annotation
  annotation_row = metadata[, c("region", "age")],  # Column annotations
  fontsize = 10
)

options(repr.plot.width=7, repr.plot.height=9)

p <- pheatmap(gaps_col = c(6), 
  deg_expression_scaled[1:50,],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = 2,
  breaks = seq(-2, 2, length.out = 100),
  color = viridis(100),
  show_rownames = T,
  show_colnames = TRUE,
  main = paste("DG Glut Top 50 Age-DEGs"),
  annotation_row = row_annotation,  # Row annotation with clusters
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)


options(repr.plot.width=7, repr.plot.height=9)

p <- pheatmap(gaps_col = c(3), 
  deg_expression_scaled[1:50,],
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cutree_rows = num_clusters,
  breaks = seq(-2, 2, length.out = 100),
  color = viridis(100),
  show_rownames = T,
  show_colnames = TRUE,
  main = paste("DG Glut Top 50 Age-DEGs"),
  annotation_row = row_annotation,  # Row annotation with clusters
  annotation_col = metadata[, c("region", "age")],  # Column annotations
  annotation_colors = annotation_colors,  # Add the color schemes
  fontsize = 10
)


bg = read.csv("DEG_results_latent_rep_mito_together/DG_Glut.csv")

clust1 = gene_clusters$Gene[which(gene_clusters$Cluster==1)] 
options(repr.plot.width=6, repr.plot.height=6)

#clust$Gene[which(clust$Cluster==1)]
c_1 = enrichGO(universe = df$X,gene = clust1,ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.5)
c_1@result$order= 1:nrow(c_1@result)
sc_1 <- simplify(
                          c_1, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_1,showCategory = 20, title = "Clust1 genes\nGO CC Enrichment")
d1

head(sc_1,20)

clust1 = gene_clusters$Gene[which(gene_clusters$Cluster==2)] 
options(repr.plot.width=6, repr.plot.height=4)

#clust$Gene[which(clust$Cluster==1)]
c_1 = enrichGO(universe = df$X,gene = clust1,ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.5)
c_1@result$order= 1:nrow(c_1@result)
sc_1 <- simplify(
                          c_1, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_1,showCategory = 10, title = "Clust2 genes\nGO CC Enrichment")
d1

clust1 = gene_clusters$Gene[which(gene_clusters$Cluster==3)] 
options(repr.plot.width=6, repr.plot.height=4)

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

d1 = dotplot(sc_1,showCategory = 10, title = "Clust3 genes\nGO CC Enrichment")
d1

clust1 = gene_clusters$Gene[which(gene_clusters$Cluster==4)] 
options(repr.plot.width=6, repr.plot.height=4)

#clust$Gene[which(clust$Cluster==1)]
c_1 = enrichGO(universe = bg$X,gene = clust1,ont = "BP",keyType = "SYMBOL",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.5)
c_1@result$order= 1:nrow(c_1@result)
sc_1 <- simplify(
                          c_1, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

d1 = dotplot(sc_1,showCategory = 10, title = "Clust4 genes\nGO CC Enrichment")
d1

options(repr.plot.width=6, repr.plot.height=6)

d1 = dotplot(sc_1,showCategory = 10, title = "Clust4 genes\nGO CC Enrichment")
d1


