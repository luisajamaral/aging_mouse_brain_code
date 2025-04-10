library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)


#genes = read.csv("../DEG_results_rm_rpl_rps/DG_Glut.csv")
genes = read.csv("../DEG_results_rm_rpl_rps//DG_Glut.csv")
genes$avg_log2FC = -genes$avg_log2FC
genes = genes[order(genes$avg_log2FC, decreasing = T),]
genesList = genes$avg_log2FC
key = bitr( as.character(genes$X), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db,drop = F)
genes$ENTREZID = key$ENTREZID
names(genesList) = as.character(genes$ENTREZID)
# Run GSEA analysis
gsego <- gseGO(geneList = genesList, minGSSize = 20,
               OrgDb = org.Mm.eg.db,
               pvalueCutoff = 0.1, 
               nPermSimple = 10000,
               verbose = FALSE)
gsego@result$order = 1:nrow(gsego@result)
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




#genes = genes[!grepl("^Rpl|^Rps|^Mrpl", genes$X),]

genes = genes[order(genes$avg_log2FC, decreasing = T),]
head(genes)

genes[which(genes$X=="Robo1"),]

genes[which(genes$X=="Zc3hav1"),]

genes[which(genes$X=="Mal"),]

genesList = genes$avg_log2FC

key = bitr( as.character(genes$X), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db,drop = F)
genes$ENTREZID = key$ENTREZID
names(genesList) = as.character(genes$ENTREZID)

# Run GSEA analysis
gsego <- gseGO(geneList = genesList, minGSSize = 20,
               OrgDb = org.Mm.eg.db,
               pvalueCutoff = 0.1, 
               nPermSimple = 10000,
               verbose = FALSE)
gsego@result$order = 1:nrow(gsego@result)

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



gsea_df <- gsea_df %>%
  mutate(core_enrichment_list = strsplit(core_enrichment, "/")) %>%
  unnest(core_enrichment_list)

# Convert Entrez IDs to Gene Symbols
converted_genes <- AnnotationDbi::select(org.Mm.eg.db,
                                         keys = unique(gsea_df$core_enrichment_list),
                                         column = "SYMBOL",
                                         keytype = "ENTREZID")

# Merge back the gene symbols into the original data
gsea_df <- gsea_df %>%
  left_join(converted_genes, by = c("core_enrichment_list" = "ENTREZID"))

# Group the gene symbols back into the core_enrichment column
gsea_df <- gsea_df %>%
  group_by(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge, Regulation) %>%
  summarise(core_enrichment_symbols = paste(SYMBOL, collapse = "/"), .groups = "drop")
gsea_df = gsea_df[order(gsea_df$p.adjust),]
# View the updated DataFrame
head(gsea_df)

head(gsea_df[which(gsea_df$NES>0),],20)

options(repr.plot.width=9, repr.plot.height=5)

library(ggplot2)

kk2 <- gseKEGG(geneList     = genesList,eps = 1e-50, minGSSize = 20,
               organism     = 'mmu',
               pvalueCutoff = 0.1,nPermSimple = 10000,
               verbose      = FALSE)

# Convert the GSEA result into a data frame
# Convert the GSEA result into a data frame
gsea_df <- as.data.frame(kk2@result)

# Clean up the pathway descriptions
gsea_df$Description <- gsub(" - Mus musculus", "", gsea_df$Description)
gsea_df$Description <- gsub("house mouse", "", gsea_df$Description)
gsea_df$Description <- gsub("[()]", "", gsea_df$Description)

# Add a column to classify pathways as Upregulated or Downregulated
gsea_df$Regulation <- ifelse(gsea_df$NES > 0, "Upregulated", "Downregulated")

# Select the top 5 upregulated and top 5 downregulated pathways
top_upregulated <- gsea_df %>%
  filter(Regulation == "Upregulated") %>%
  arrange(p.adjust) %>%
  head(10)

top_downregulated <- gsea_df %>%
  filter(Regulation == "Downregulated") %>%
  arrange(p.adjust) %>%
  head(10)

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


gsea_df <- gsea_df %>%
  mutate(core_enrichment_list = strsplit(core_enrichment, "/")) %>%
  unnest(core_enrichment_list)

# Convert Entrez IDs to Gene Symbols
converted_genes <- AnnotationDbi::select(org.Mm.eg.db,
                                         keys = unique(gsea_df$core_enrichment_list),
                                         column = "SYMBOL",
                                         keytype = "ENTREZID")

# Merge back the gene symbols into the original data
gsea_df <- gsea_df %>%
  left_join(converted_genes, by = c("core_enrichment_list" = "ENTREZID"))

# Group the gene symbols back into the core_enrichment column
gsea_df <- gsea_df %>%
  group_by(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue, rank, leading_edge, Regulation) %>%
  summarise(core_enrichment_symbols = paste(SYMBOL, collapse = "/"), .groups = "drop")
gsea_df = gsea_df[order(gsea_df$p.adjust),]
# View the updated DataFrame
head(gsea_df)

gsea_df[which(gsea_df$NES>0),]

gsea_df$Description

library("pathview")
hsa04110 <- pathview(gene.data  = genesList,
                     pathway.id = "mmu04080",
                     species    = "mmu",
                     limit      = list(gene=max(abs(genesList)), cpd=1))

gseaplot2(gsego, geneSetID = "GO:0001774",title = "Microglial activation")  # Replace with a specific pathway ID


gseaplot2(gsego, geneSetID = "GO:0051607",title = "Defense response to virus")  # Replace with a specific pathway ID


options(repr.plot.width=5, repr.plot.height=5)

gseaplot2(kk2, geneSetID = "mmu03010",title = "Ribosome - Mus musculus")  # Replace with a specific pathway ID
gseaplot2(kk2, geneSetID = "mmu00190",title = "Oxidative phosphorylation")  # Replace with a specific pathway ID
gseaplot2(kk2, geneSetID = "mmu04720",title = "Long Term Potentiation")  # Replace with a specific pathway ID
gseaplot2(kk2, geneSetID = "mmu05012",title = "Parkinson disease")  # Replace with a specific pathway ID
gseaplot2(kk2, geneSetID = "mmu04145",title = "Phagosome")  # Replace with a specific pathway ID




head(genes)

options(repr.plot.width=14, repr.plot.height=14)


go_bp_down <- enrichGO(gene = genes$ENTREZID[which(genes$avg_log2FC<(-.25) &genes$p_val_adj<0.05)],ont = "CC",
                       OrgDb = org.Mm.eg.db, 
                 pvalueCutoff = 0.05)
head(go_bp_down,5)

go_bp_up <- enrichGO(gene = genes$ENTREZID[which(genes$avg_log2FC>(.25) &genes$p_val_adj<0.05)],ont = "CC",OrgDb = org.Mm.eg.db,
                    ,  pvalueCutoff = 0.05)
head(go_bp_up,5)


options(repr.plot.width=7, repr.plot.height=10)

dotplot(go_bp_up, showCategory = 15, title  = "up")  # Show the top 20 pathways
dotplot(go_bp_down, showCategory = 15,title  = "down")  # Show the top 20 pathways


options(repr.plot.width=5, repr.plot.height=5)

gseaplot2(kk2, geneSetID = "mmu03010",title = "Ribosome - Mus musculus")  # Replace with a specific pathway ID
gseaplot2(kk2, geneSetID = "mmu00190",title = "Oxidative phosphorylation")  # Replace with a specific pathway ID
gseaplot2(kk2, geneSetID = "mmu04720",title = "Long Term Potentiation")  # Replace with a specific pathway ID



gseaplot2(gsego, geneSetID = "GO:0140236",title = "translation at presynapse")  # Replace with a specific pathway ID


gsego <- gseGO(geneList     = genesList,OrgDb = org.Mm.eg.db,
               pvalueCutoff = 1,
               verbose      = FALSE)



kk2 <- gseKEGG(geneList     = genesList,eps = 1e-50, minGSSize = 20,
               organism     = 'mmu',
               pvalueCutoff = 0.05,nPermSimple = 10000,
               verbose      = FALSE)
head(kk2,2)


mkk2 <- gseMKEGG(geneList = genesList,
                 organism = 'mmu',
                 pvalueCutoff = 1)
head(mkk2)

mkk2[order(mkk2$NES),]

kk2[order(kk2$rank, decreasing = T),]



library("pathview")
hsa04110 <- pathview(gene.data  = genesList,
                     pathway.id = "mmu04724",
                     species    = "mmu",
                     limit      = list(gene=max(abs(genesList)), cpd=1))



list(gene=max(abs(genesList)), cpd=1)

browseKEGG(kk2, 'mmu03010')


hsa04110 <- pathview(gene.data  = genesList,
                     pathway.id = "mmu03010",
                     species    = "mmu",
                     limit      = list(gene=max(abs(genesList)), cpd=1))


