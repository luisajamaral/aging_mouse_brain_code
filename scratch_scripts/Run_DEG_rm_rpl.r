library(Seurat)
library(ggplot2)

obj = readRDS("../female_RNA/RNA_final_SCT.RDS")

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


# Extract expression matrix and log-transform it
expr_matrix <- GetAssayData(obj, slot = "data")[qc_genes, , drop = FALSE]  # This retrieves log-normalized values

# Sum log-transformed expression for each cell to compute QC Score
obj$QC_Score <- colSums(expr_matrix)



# Visualize QC Score across samples
ggplot(obj@meta.data, aes(x = orig.ident, y = QC_Score, fill = age)) +
  geom_boxplot() +
  theme_minimal() +coord_flip()+
  labs(title = "QC Score Across Samples", x = "Sample", y = "QC Score")
# Visualize QC Score across samples
ggplot(obj@meta.data, aes(x = orig.ident, y = percent.ribo, fill = age)) +
  geom_boxplot() +
  theme_minimal() +coord_flip()+
  labs(title = "QC Score Across Samples", x = "Sample", y = "QC Score")




FeaturePlot(obj, features = "percent.ribo")

qc_threshold <- quantile(obj$percent.ribo, probs = 0.9)  # Bottom 5% cutoff
qc_threshold
bad_cells <- rownames(obj@meta.data[obj$percent.ribo > qc_threshold, ])

# Number of bad cells
length(bad_cells)


# Count total cells per sample before filtering
total_cells_per_sample <- obj@meta.data %>%
  group_by(orig.ident) %>%
  summarise(total_cells = n())

# Count removed cells per sample
removed_cells_per_sample <- obj@meta.data[bad_cells, ] %>%
  group_by(orig.ident) %>%
  summarise(removed_cells = n())

# Merge data and compute percentage lost
loss_summary <- left_join(total_cells_per_sample, removed_cells_per_sample, by = "orig.ident") %>%
  mutate(removed_cells = replace_na(removed_cells, 0),  # Replace NA with 0 for samples with no removed cells
         percent_lost = (removed_cells / total_cells) * 100)

# View data
print(loss_summary)


# Compute retained cells
loss_summary <- loss_summary %>%
  mutate(retained_cells = total_cells - removed_cells)

# Convert to long format for ggplot
loss_long <- loss_summary %>%
  dplyr::select(orig.ident, removed_cells, retained_cells) %>%
  tidyr::pivot_longer(cols = c(removed_cells, retained_cells), 
                       names_to = "status", values_to = "cell_count")

# Rename for better visualization
loss_long$status <- recode(loss_long$status, 
                           removed_cells = "Removed",
                           retained_cells = "Retained")

ggplot(loss_long, aes(x = orig.ident, y = cell_count, fill = status)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of Cells Retained vs. Removed",
       x = "Sample",
       y = "Number of Cells",
       fill = "Cell Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(loss_long, aes(x = orig.ident, y = cell_count, fill = status)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Proportion of Cells Retained vs. Removed",
       x = "Sample",
       y = "Number of Cells",
       fill = "Cell Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


 library(ggplot2)

ggplot(obj@meta.data[which(obj@meta.data$region == "HCA"),], aes(x = nCount_RNA, fill = age ,color = as.factor(rep))) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribution of nCount_RNA",
       x = "QC Score",
       y = "Density") + facet_wrap(~orig.ident,ncol = 1)


ggplot(obj@meta.data[which(obj@meta.data$region == "HCA"),], aes(x = QC_Score, fill = age ,color = as.factor(rep))) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribution of QC Scores Across Samples",
       x = "QC Score",
       y = "Density") + facet_wrap(~orig.ident,ncol = 1)


ggplot(obj@meta.data[which(obj@meta.data$region == "HCA"),], aes(x = percent.ribo, fill = age ,color = as.factor(rep))) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Distribution of QC Scores Across Samples",
       x = "QC Score",
       y = "Density") + facet_wrap(~orig.ident,ncol = 1)


Idents(obj) = "celltype_final"

setwd("DEG_results_latent_qc_region")

obj$celltype_final[which(obj$best_celltype_fixed == "IOL")] = "IOL NN"

obj$age_celltype = paste(obj$age , obj$celltype_final)
Idents(obj) = "age_celltype"
head(Idents(obj))

hca = obj[,which(obj$region == "HCA")]
hca
head(hca)

hca <- SCTransform(hca, vars.to.regress = "percent.ribo", verbose = FALSE)


 
de_results <- FindMarkers(
    hca,
    ident.1 = '18mo DG Glut',
    ident.2 = '2mo DG Glut',
    test.use = "MAST",
    logfc.threshold = 0.1, 
    min.pct = 0.05, latent.vars = c('percent.ribo')
  )

head(de_results)

head(de_results)
de_results$X= rownames(de_results)

library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)# Connect to Ensembl
library(ggrepel)
library(AnnotationDbi)
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # Change to "hsapiens_gene_ensembl" for human

# Retrieve gene types for all genes in your dataset
gene_annotations <- getBM(attributes = c("external_gene_name", "gene_biotype","transcript_length"),
                          filters = "external_gene_name",
                          values = rownames(obj),
                          mart = mart)

# Filter for lncRNAs
lncRNA_genes <- gene_annotations$external_gene_name[gene_annotations$gene_biotype == "lncRNA"]
protein_coding_genes <- gene_annotations$external_gene_name[gene_annotations$gene_biotype == "protein_coding"]

print(ma_plot)

deg_results = de_results
deg_results$annotation = "Other"

# Add classification
df <- deg_results %>%
  mutate(
    significance = case_when(
      avg_log2FC > 0.25 & p_val_adj < 1e-5 & annotation == "Other" ~ "Upregulated Other",
      avg_log2FC < -0.25 & p_val_adj < 1e-5 & annotation == "Other" ~ "Downregulated Other",
      avg_log2FC > 0.25 & p_val_adj < 1e-5 & annotation == "protein_coding" ~ "Upregulated protein_coding",
      avg_log2FC < -0.25 & p_val_adj < 1e-5 & annotation == "protein_coding" ~ "Downregulated protein_coding",
      avg_log2FC > 0.25 & p_val_adj < 1e-5 & annotation == "lncRNA" ~ "Upregulated lncRNA",
      avg_log2FC < -0.25 & p_val_adj < 1e-5  & annotation == "lncRNA" ~ "Downregulated lncRNA",
      TRUE ~ "Not Significant"
    )
  )
table(df$significance)
df$pct_exp = (df$`pct.1`+  df$`pct.2`)/2
head(df)
# Select top 10 differential genes
top_genes <- rbind(df[which(df$avg_log2FC>1)[1:10], ],df[which(df$avg_log2FC< -1)[1:10], ])
#top_genes = top_genes[-grep("Rik", top_genes$X),]
#top_genes = top_genes[-grep("Gm", top_genes$X),]

# MA Plot
ma_plot <- ggplot(df, aes(x = (pct_exp), y = avg_log2FC)) + # Log-transform base_mean
  geom_point(aes(color = significance), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Upregulated protein_coding" = "red","Upregulated lncRNA"="magenta","Upregulated Other" = "darkorange",
                                "Downregulated Other" = "deepskyblue",
                                "Downregulated lncRNA"="purple",
                                "Downregulated protein_coding" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +

  labs(
    title = "DG_Glut MA Plot",
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

#DNA damage response GO:0006974
genes_in_go <- AnnotationDbi::select(
  org.Mm.eg.db, 
  keys = "GO:0006974", 
  keytype = "GO", 
  columns = c("SYMBOL")
)

DNA_damage_response = unique(genes_in_go$SYMBOL)
head(DNA_damage_response)

# Add classification
df <- deg_results %>%
  mutate(
    pathway = case_when(
      X %in% DNA_damage_response ~ "DNA damage response",
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

print(ma_plot2)

gseaplot2(gsego, geneSetID = "GO:0006974",title = "Microglial activation")  # Replace with a specific pathway ID


genesList = de_results$avg_log2FC

key = bitr( as.character(de_results$X), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db,drop = F)
de_results$ENTREZID = key$ENTREZID
names(genesList) = as.character(de_results$ENTREZID)
genesList = genesList[order(genesList, decreasing = T)]
# Run GSEA analysis
gsego <- gseGO(geneList = genesList, minGSSize = 20,
               OrgDb = org.Mm.eg.db,
               pvalueCutoff = 0.3, 
               nPermSimple = 100000,
               verbose = FALSE)
gsego@result$order = 1:nrow(gsego@result)

head(gsea_df)

options(repr.plot.width=12, repr.plot.height=7)


# Simplify GO terms
simplified_gsea <- simplify(
                          gsego, 
                          cutoff = 0.7,               # Similarity cutoff
                          by = "order",                 # Use NES for filtering
                          select_fun = min,#function(x) x[which.max(abs(x))],  # Custom function for absolute NES
                          measure = "Wang"            # Semantic similarity measure
                        )

# Convert to a data frame
gsea_df <- as.data.frame(simplified_gsea@result)
# Clean up the pathway descriptions
gsea_df$Description <- gsub(" - Mus musculus", "", gsea_df$Description)
gsea_df$Description <- gsub("house mouse", "", gsea_df$Description)
gsea_df$Description <- gsub("[()]", "", gsea_df$Description)

# Add a column to classify pathways as Upregulated or Downregulated
gsea_df$Regulation <- ifelse(gsea_df$NES > 0, "Upregulated", "Downregulated")

#gsea_df = gsea_df[-grep("process", gsea_df$Description),]

# Select the top 5 upregulated and top 5 downregulated pathways
top_upregulated <- gsea_df %>%
  filter(Regulation == "Upregulated") %>%
  arrange(p.adjust) %>%
  head(25)

top_downregulated <- gsea_df %>%
  filter(Regulation == "Downregulated") %>%
  arrange(p.adjust) %>%
  head(25)

# Combine the top pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated)

# Plot using ggplot2
ggplot(top_pathways, aes(x = NES, 
                         y = reorder(Description, NES), 
                         size = -log10(p.adjust), 
                         color = Regulation)) +
  geom_point() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       size = "-log10(Adjusted P-value)",
       color = "Regulation") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


#clust$Gene[which(clust$Cluster==1)]
df = de_results
up = df$X[which(df$p_val_adj<1e-3 & df$avg_log2FC>0.1)]
length(up)
c_up = enrichGO(universe = df$X, gene = up,ont = "BP",
                keyType = "SYMBOL",
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

d1 = dotplot(sc_up,showCategory = 19, title = "Age-Up genes\nGO BP Enrichment")
d1


Idents(obj) = "age_celltype"
# Loop through each unique cell type
for (celltype in unique(obj$celltype_final)[1:length(unique(obj$celltype_final))]) {
#for (celltype in c("Astro-NT NN","DG Glut")) {
  # Define unique identity combining age and cell type
  cat(celltype , "\n")  
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("18mo", celltype, sep = "_"))

  de_results <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    logfc.threshold = 0, 
    min.pct = 0.01, latent.vars = c('region','QC_Score')
  )
  head(de_results)
  # Save the results to a CSV file for each cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"csv", sep = ".")
  write.csv(de_results, file = file_name, row.names = TRUE)
}




DefaultAssay(obj)="RNA"

DimPlot(obj, label = T)+NoLegend()

obj$age_celltype <- gsub(" ", "_", paste(obj$age, obj$celltype_final, sep = "_"))


setwd("DEG_results_rm_rpl_rps")

all_genes <- rownames(obj)

# Identify ribosomal genes (RPL and RPS)
ribosomal_genes <- all_genes[grepl("^Rpl|^Rps", all_genes)]

# Print the ribosomal genes for verification
cat("Ribosomal genes to be removed:\n")
print(ribosomal_genes)


mito_genes <- grep("^mt-", rownames(obj), value = TRUE)
mito_genes

mito_ribo <- grep("Mrp", rownames(obj), value = TRUE)
mito_ribo

genes_to_remove <- unique(c(ribosomal_genes, mito_genes,mito_ribo))
obj
# Remove genes from the Seurat object
obj <- subset(obj, features = setdiff(rownames(obj), genes_to_remove))
obj

colnames(obj@meta.data)

obj$age_celltype =  gsub(" ", "_",obj$age_celltype)
Idents(obj) = "age_celltype"

head(Idents(obj))


which(celltype == unique(obj$celltype_final))

setwd("../DEG_results_rm_rpl_rps_mrpl_mt_latent_rep_region_mito_ribo")


Idents(obj) = "age_celltype"
# Loop through each unique cell type
for (celltype in unique(obj$celltype_final)[11:length(unique(obj$celltype_final))]) {
#for (celltype in c("Astro-NT NN","DG Glut")) {
  # Define unique identity combining age and cell type
  cat(celltype , "\n")  
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("18mo", celltype, sep = "_"))
  if(length(which(obj$age_celltype ==ident1))>20 & length(which(obj$age_celltype ==ident2))>20){

  de_results <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    logfc.threshold = 0.1, 
    min.pct = 0.01, latent.vars = c('region','rep', 'percent.mt', 'percent.ribo')
  )
  head(de_results)
  # Save the results to a CSV file for each cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"csv", sep = ".")
  write.csv(de_results, file = file_name, row.names = TRUE)
 }
}


setwd("../DEG_results_region_rm_rpl_rps_mrpl_mt_latent_rep_region_mito_ribo/")

obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) = "age_celltype_region"
for(r in unique(obj$region)){
# Loop through each unique cell type
    #for (celltype in unique(obj$celltype_final)) {
    for (celltype in unique(obj$celltype_final)) {
      # Define unique identity combining age and cell type
      file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"--",r, "_2vs18.csv", sep = "")

      if(file.exists(file_name)){
          cat(file_name, "file exists")
          next
     }
      ident1 <- gsub(" ", "_",paste("2mo", celltype,r, sep = "_"))
      ident2 <- gsub(" ", "_",paste("18mo", celltype, r,sep = "_"))
      if(length(which(obj$age_celltype_region ==ident1))>20 & length(which(obj$age_celltype_region ==ident2))>20){
          cat("Running", file_name)
          de_results <- FindMarkers(
            obj,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            logfc.threshold = 0.1, 
            min.pct = 0.01, latent.vars = c('rep', 'percent.mt', 'percent.ribo')
          )
          head(de_results)
          # Save the results to a CSV file for each cell type
          write.csv(de_results, file = file_name, row.names = TRUE)
        }
    }
}


setwd("../DEG_results_rm_rpl_rps_mrpl_mt")


Idents(obj) = "age_celltype"
# Loop through each unique cell type
for (celltype in unique(obj$celltype_final)[1:length(unique(obj$celltype_final))]) {
#for (celltype in c("Astro-NT NN","DG Glut")) {
  # Define unique identity combining age and cell type
  cat(celltype , "\n")  
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("18mo", celltype, sep = "_"))

  de_results <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    logfc.threshold = 0, 
    min.pct = 0.01, latent.vars = c('region','rep')
  )
  head(de_results)
  # Save the results to a CSV file for each cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"csv", sep = ".")
  write.csv(de_results, file = file_name, row.names = TRUE)
}


length(unique(obj$celltype_final)==celltype)


Idents(obj) = "age_celltype"
# Loop through each unique cell type
for (celltype in unique(obj$celltype_final)[11:length(unique(obj$celltype_final))]) {
#for (celltype in c("Astro-NT NN","DG Glut")) {
  # Define unique identity combining age and cell type
  cat(celltype , "\n")  
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("18mo", celltype, sep = "_"))

  de_results <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    logfc.threshold = 0, 
    min.pct = 0.01, latent.vars = c('region','rep')
  )
  head(de_results)
  # Save the results to a CSV file for each cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"csv", sep = ".")
  write.csv(de_results, file = file_name, row.names = TRUE)
}



Idents(obj) = "age_celltype"
# Loop through each unique cell type
for (celltype in unique(obj$celltype_final)[11:length(unique(obj$celltype_final))]) {
#for (celltype in c("Astro-NT NN","DG Glut")) {
  # Define unique identity combining age and cell type
  cat(celltype , "\n")  
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("18mo", celltype, sep = "_"))

  de_results <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    logfc.threshold = 0, 
    min.pct = 0.01, latent.vars = c('region','rep')
  )
  head(de_results)
  # Save the results to a CSV file for each cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"csv", sep = ".")
  write.csv(de_results, file = file_name, row.names = TRUE)
}


unique(obj$celltype_final)

which(unique(obj$celltype_final) == celltype)

ident1

ident2

ident2%in%Idents(obj)

obj@meta.data[which(Idents(obj) == ident2),]


