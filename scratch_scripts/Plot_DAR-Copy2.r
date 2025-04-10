library(ggplot2)
library(data.table)
library(UpSetR)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

library(scales)

setwd("~/projects/combined_all/region_DARs_redone/")

mod_sex_colors = data.table(ATAC_Male="#2ba0f8", RNA_Female="#f76dd2", ATAC_Female = "#e62578")

rmeta = read.table("../Figures/rna_meta.csv")
rmeta$celltype_region = paste(rmeta$celltype_final, "--", rmeta$region,sep = "")
rmeta$celltype_region = gsub(" " , "_", rmeta$celltype_region )

meta = read.csv("../final_meta.csv")
meta$celltype_final = gsub("Lymphoid","Lymphoid NN", meta$celltype_final)
meta = meta[-which(meta$celltype_final=="doublet"),]
meta$celltype_final=paste(meta$celltype_final)
meta$celltype_final[which(meta$best_celltype_fixed=="IOL")] = "IOL NN"

meta$celltype_region = paste(meta$celltype_final, meta$region, sep = ":")
meta$celltype_region=gsub(" ", "_", meta$celltype_region)
meta$age = factor(meta$age , levels = c("2mo", "9mo", "18mo"))


meta$celltype_region_modality = paste(meta[,'celltype_region'], "ATAC", sep = "_")
rmeta$celltype_region_modality = paste(rmeta[,'celltype_region'], "RNA", sep = "_")

ameta = rbind(meta[,c('sample', 'celltype_region_modality','celltype_region', 'celltype_final', 'region')], rmeta[,c('sample', 'celltype_region_modality','celltype_region','celltype_final', 'region')])
ameta$celltype_region_modality = gsub("--", ":", ameta$celltype_region_modality)
ameta$celltype_region = gsub("--", ":", ameta$celltype_region)

ameta$Modality = "ATAC"
ameta$Modality[grep("RNA", ameta$celltype_region_modality)] = "RNA"


ameta$Sex = "Female"
ameta$Sex[grep("Male", ameta$sample)] ="Male"
ameta$Modality_Sex = paste(ameta$Modality, ameta$Sex) 

cl="Gaba"

setwd("~/projects/combined_all/region_DARs_redone/")
options(repr.plot.width=14, repr.plot.height=5)
files = list.files(".", cl)
files = files[grep(".csv", files)]
files = files[grep("diff", files)]

all = list()
for(f in files){
    curr = read.csv(f)
    id = gsub("diff_peaks_" , "", f)
    id = gsub("_2vs18.csv" , "", id)

    curr$id = id 
    all[[id]] = curr
    }


tab = do.call(rbind, all)
tab=as.data.frame(tab)
tab$sig = "no"
#tab$log2.fold_change. = -tab$log2.fold_change.
tab[which(tab$adjusted.p.value<0.01 & tab$log2.fold_change.< -0.25 ),"sig"] = "down"
tab[which(tab$adjusted.p.value<0.01 & tab$log2.fold_change.>0.25 ),"sig"] = "up"
tab$celltype = sapply(strsplit(as.character(tab$id), ":"), `[`, 1)
tab$region = sapply(strsplit(as.character(tab$id), ":"), `[`, 2)



setwd("~/projects/combined_all/female_RNA/DEG_results_latent_rep_mito/")
logfc = 0.25
files = list.files(".", cl)
files = files[grep(".csv", files)]
all = list()
for(f in files){
    curr = read.csv(f)
    id = gsub(".csv" , "", f)
    curr$id = id 
    all[[id]] = curr
}
rtab = do.call(rbind, all)
rtab=as.data.frame(rtab)
rtab$avg_log2FC = -rtab$avg_log2FC
rtab$sig = "no"
rtab[which(rtab$p_val_adj<0.01 & rtab$avg_log2FC<(-logfc) ),"sig"] = "down"
rtab[which(rtab$p_val_adj<0.01 & rtab$avg_log2FC>logfc ),"sig"] = "up"
#tab[which(tab$`pct.1`<0.05 |  tab$`pct.2`<0.05), "sig"] = "no"
rtab$Celltype= sapply(strsplit(as.character(rtab$id), "--"), `[`, 1)


head(rtab[order(rtab$avg_log2FC, decreasing = T),], 20)

rtab = rtab[-grep("Rpl", rtab$X),]

rtab$sig = "no"

rtab[which(rtab$p_val_adj<0.01 & rtab$avg_log2FC<(-logfc) ),"sig"] = "down"
rtab[which(rtab$p_val_adj<0.01 & rtab$avg_log2FC>logfc ),"sig"] = "up"

length(unique(rtab[which(rtab$sig=="up"), "X"]))

length(unique(rtab[which(rtab$sig=="down"), "X"]))

hist(table(rtab[which(rtab$sig=="down"), "X"]))

table(rtab[which(rtab$sig=="down"), "X"])[order(table(rtab[which(rtab$sig=="down"), "X"]),decreasing = T)]

hist(table(rtab[which(rtab$sig=="up"), "X"]))

options(repr.plot.width=10, repr.plot.height=11)

significant_points <- tab %>%
  filter(`adjusted.p.value` < 0.01 & abs(`log2.fold_change.`) > 0.25)

non_significant_points <- tab %>%
  filter(`adjusted.p.value` >= 0.01 | abs(`log2.fold_change.`) < 0.25) %>%
  group_by(celltype) %>%
  sample_n(min(n(), 1000))  # Sample max 1000 points per celltype

# Step 2: Combine the significant and sampled non-significant points
sampled_tab <- bind_rows(significant_points, non_significant_points)



significant_points <- rtab %>%
  filter(`p_val_adj` < 0.05 & abs(`avg_log2FC`) > 0.25)

non_significant_points <- rtab %>%
  filter(`p_val_adj` >= 0.05 | abs(`avg_log2FC`) < 0.25) %>%
  group_by(Celltype) %>%
  sample_n(min(n(), 1000))  # Sample max 1000 points per celltype

# Step 2: Combine the significant and sampled non-significant points
rsampled_tab <- bind_rows(significant_points, non_significant_points)

rsampled_tab$id = gsub("--", ":", rsampled_tab$id)


head(rsampled_tab)

g1 = ggplot(rsampled_tab[which(rsampled_tab$sig !="no"),], aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black", position = "dodge") +
  theme_minimal() + 
  labs(title = "Cell #", x = "", y = "Number of Cells") +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        legend.position = "right",plot.title = element_text(hjust = 0.5),
        strip.text.y = element_blank()) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                     labels = scales::label_number(scale_cut = cut_short_scale())) 
g1

keep = unique(rsampled_tab$id)[unique(rsampled_tab$id) %in% unique(sampled_tab$id)]


me =ameta[grep(cl,ameta$celltype_final),]
me$celltype_final = gsub("/", "-" , me$celltype_final)
me$celltype_region = gsub("/", "-" , me$celltype_region)
#keep = table(me$celltype_region)[which(table(me$celltype_region)>1000)]
#keep = names(keep)
me = me[which(me$celltype_region%in%keep),]
#me = me[which(me$celltype_final %in% c("Oligo NN", "OPC NN", "Microglia NN", "Astro-TE NN","Astro-NT NN" , "VLMC NN")),]
me$Modality_Sex = gsub(" ", "_", me$Modality_Sex)
ord =  names(table(me$celltype_final))[order(table(me$celltype_final), decreasing = T)]
me$celltype_final = factor(me$celltype_final, levels = ord)

sampled_tab$celltype = factor(sampled_tab$celltype, levels = gsub(" ", "_", ord))
sampled_tab = sampled_tab[which(sampled_tab$id %in%keep), ]
rsampled_tab$Celltype = factor(rsampled_tab$Celltype, levels = gsub(" ", "_", ord))
rsampled_tab = rsampled_tab[which(rsampled_tab$id %in%keep), ]


sampled_tab = sampled_tab[which(gsub("--", ":", sampled_tab$id) %in% unique(me$celltype_region)),]
rsampled_tab = rsampled_tab[which(gsub("--", ":", rsampled_tab$id) %in% unique(me$celltype_region)),]


# First plot (g1)
g1 = ggplot(me, aes(x = celltype_region, fill = Modality_Sex)) +
  geom_bar(stat = "count", color = "black", position = "dodge") +
  theme_minimal() + 
  scale_fill_manual(values = mod_sex_colors) +
  labs(title = "Cell #", x = "", y = "Number of Cells") +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        legend.position = "right",plot.title = element_text(hjust = 0.5),
        strip.text.y = element_blank()) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                     labels = scales::label_number(scale_cut = cut_short_scale())) +
  facet_grid(celltype_final ~ ., space = "free", scales = "free_y")

# Second plot (darlogfc)
darlogfc = ggplot(sampled_tab, 
  aes(x = id, y = `log2.fold_change.`, 
      color = sig)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "logFC Age-DARs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-10, 10)) +  # Set the y-axis limits here
  theme_minimal() +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        strip.text.y = element_blank(),plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.y = element_blank()) +  # Remove y-axis labels
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
  facet_grid(celltype ~ ., space = "free", scales = "free_y")



# Second plot (darlogfc)
# Select the top upregulated and downregulated gene per id (cell type)
top_degs <- rsampled_tab %>%
  filter(`p_val_adj` < 1e-10, (`avg_log2FC`) > 1) %>%
  group_by(id) %>%
  slice_max(order_by = avg_log2FC, n = 1, with_ties = FALSE) %>%  # Get top upregulated
  bind_rows(
    rsampled_tab %>%
      filter(`p_val_adj` < 1e-10, (`avg_log2FC`) < -1) %>%
      group_by(id) %>%
      slice_min(order_by = avg_log2FC, n = 1, with_ties = FALSE)  # Get top downregulated
  ) %>%
  ungroup() %>%
  mutate(nudge = ifelse(`avg_log2FC` > 0, 10, -10))  # Nudging +1 for up, -1 for down

# Plot
deglogfc <- ggplot(rsampled_tab, 
  aes(x = id, y = `avg_log2FC`, color = sig)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "logFC Age-DEGs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-10, 10)) +  
  theme_minimal() +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        strip.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.y = element_blank(), 
        axis.text.y = element_blank()
       ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
  facet_grid(Celltype ~ ., space = "free", scales = "free_y") +

  # Add gene labels for the top up/down regulated genes
  geom_text_repel(data = top_degs, 
                  aes(label = X, y =  nudge), 
                  size = 3, color = "black")
# Combine the two plots side by side
#combined_plot <- g1 + darlogfc + deglogfc +plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "right")

# Show the combined plot
#combined_plot



full_celltypes <- data.frame(id = keep)  # List of all cell types

# Step 1: Summarize the number of up and down DARs for each cell type
summary_tab <- tab %>%
  filter(`adjusted.p.value` < 0.01 & abs(`log2.fold_change.`)>0.5) %>%
  group_by(id, celltype) %>%
  summarise(
    num_down = sum(sig == "down"),  # Count the down DARs
    num_up = sum(sig == "up")     # Count the up DARs
  )

summary_tab_full <- full_celltypes %>%
  left_join(summary_tab, by = "id") %>%
  replace_na(list(num_down = 0, num_up = 0))  # Fill missing values with 0

summary_tab_full$celltype = sapply(strsplit(as.character(summary_tab_full$id), ":"), `[`, 1)
summary_tab_full$celltype = factor(summary_tab_full$celltype, levels = gsub(" ", "_", ord))

# Step 2: Create a mirror plot with the summarized data
mirror_plot <- ggplot(summary_tab_full, aes(x = id)) +
  geom_bar(aes(y = -num_down), stat = "identity", fill = "blue") +  # Down DARs (left side)
  geom_bar(aes(y = num_up), stat = "identity", fill = "red") +      # Up DARs (right side)
  coord_flip() +  # Flip the axes for horizontal bars
  labs(title = "# Age-DARs", x = "Cell Type", y = "Number of DARs") +
  theme_minimal() +
  scale_y_continuous(limits = c(-max(summary_tab_full$num_down,summary_tab_full$num_up), max(summary_tab_full$num_down,summary_tab_full$num_up)),breaks = scales::pretty_breaks(n = 4),
                     labels = scales::label_number(scale_cut = cut_short_scale())) +

  theme(text = element_text(color = "black", size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",    strip.text.y = element_blank(),
        axis.title.y = element_blank(), axis.text.y = element_blank()
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +   # Add a dotted line at y = 0
  facet_grid(celltype ~ ., space = "free", scales = "free_y")


full_celltypes <- data.frame(id = keep)  # List of all cell types

# Step 1: Summarize the number of up and down DARs for each cell type
summary_tab <- rtab %>%
  filter(`p_val_adj` < 0.01 & abs(`avg_log2FC`)>0.5) %>%
  group_by(id, Celltype) %>%
  summarise(
    num_down = sum(sig == "down"),  # Count the down DARs
    num_up = sum(sig == "up")     # Count the up DARs
  )

summary_tab$id = gsub("--", ":", summary_tab$id)

summary_tab_full <- full_celltypes %>%
  left_join(summary_tab, by = "id") %>%
  replace_na(list(num_down = 0, num_up = 0))  # Fill missing values with 0

summary_tab_full$celltype = sapply(strsplit(as.character(summary_tab_full$id), ":"), `[`, 1)
summary_tab_full$celltype = factor(summary_tab_full$celltype, levels = gsub(" ", "_", ord))

# Step 2: Create a mirror plot with the summarized data
rmirror_plot <- ggplot(summary_tab_full, aes(x = id)) +
  geom_bar(aes(y = -num_down), stat = "identity", fill = "blue") +  # Down DARs (left side)
  geom_bar(aes(y = num_up), stat = "identity", fill = "red") +      # Up DARs (right side)
  coord_flip() +  # Flip the axes for horizontal bars
  scale_y_continuous(limits = c(-max(summary_tab_full$num_down,summary_tab_full$num_up), max(summary_tab_full$num_down,summary_tab_full$num_up)),breaks = scales::pretty_breaks(n = 4),
                     labels = scales::label_number(scale_cut = cut_short_scale())) +
  labs(title = "# Age-DEGs", x = "Cell Type", y = "Number of DARs") +
  theme_minimal() +
  theme(text = element_text(color = "black", size = 13),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",    strip.text.y = element_blank(),    axis.title.y = element_blank(), axis.text.y = element_blank()
) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +   # Add a dotted line at y = 0
  facet_grid(celltype ~ ., space = "free", scales = "free_y")

options(repr.plot.width=15, repr.plot.height=13)

combined_plot <- g1 + mirror_plot + darlogfc + rmirror_plot+ deglogfc +plot_layout(ncol = 5, guides = "collect") & theme(legend.position = "right")
combined_plot

ggsave(combined_plot, file = paste("~/projects/combined_all/Figures/Gaba_DAR_DEG_logFC.pdf", sep = "")
       , height = 9, width = 13)

head(rsampled_tab)
rsampled_tab$region = sapply(strsplit(as.character(rsampled_tab$id), ":"), `[`, 2)


darlogfc = ggplot(sampled_tab, 
  aes(x = id, y = `log2.fold_change.`, 
      color = sig)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "logFC Age-DARs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-10, 10)) +  # Set the y-axis limits here
  theme_minimal() +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        strip.text.y = element_blank(),plot.title = element_text(hjust = 0.5),
        legend.position = "none",
      #  axis.title.y = element_blank(),  # Remove y-axis title
        #axis.text.y = element_blank()
       ) +  # Remove y-axis labels
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
  facet_grid(celltype ~ ., space = "free", scales = "free_y")



library(dplyr)
library(ggplot2)
library(ggrepel)

# Identify top upregulated and downregulated genes per id
cur = rsampled_tab[-grep("Rik", rsampled_tab$X),]
cur = cur[-grep("Gm", cur$X),]
cur = cur[-which(cur$p_val_adj>0.05),]
cur = cur[-which(abs(cur$avg_log2FC)<1.5),]

top_genes <- cur %>%
  group_by(id) %>%
  summarise(
    top_up = X[which.max(avg_log2FC)[1]],   # Gene with the highest logFC
    top_down = X[which.min(avg_log2FC)[1]]  # Gene with the lowest logFC
  ) %>%
  pivot_longer(cols = c("top_up", "top_down"), names_to = "regulation", values_to = "gene")


# Merge with original data to get coordinates for top genes
label_data <- rsampled_tab %>%
  semi_join(top_genes, by = c("id", "X" = "gene"))
#label_data = label_data[-which(label_data$p_val_adj>0.05),]
label_data$avg_log2FC = label_data$avg_log2FC*1.25
#label_data = label_data[-grep("Rik", label_data$X),]

# Plot with text labels for top up and downregulated genes
deglogfc <- ggplot(rsampled_tab, 
  aes(x = id, y = `avg_log2FC`, 
      color = sig)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "logFC Age-DEGs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-10, 10)) +  # Set the y-axis limits here
  theme_minimal() +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        strip.text.y = element_blank(), plot.title = element_text(hjust = 0.5),
        legend.position = "none",
         axis.title.y = element_blank(),  # Remove y-axis title
         axis.text.y = element_blank()
       ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
  facet_grid(Celltype ~ ., space = "free", scales = "free_y") +
  geom_text_repel(data = label_data, aes(label = X), 
                  size = 3, nudge_x = 0, direction = "y", max.overlaps = Inf)


deglogfc <- deglogfc + theme(legend.position = "none")

combined_plot_short <- darlogfc +deglogfc  +plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
combined_plot_short

options(repr.plot.width=10, repr.plot.height=11)

combined_plot_short

ggsave(combined_plot_short, file = paste("~/projects/combined_all/Figures/Figure5-Glut/short_DAR_DEG_logFC.pdf", sep = "")
       , height = 11, width = 10)

ggsave(combined_plot_short, file = paste("~/projects/combined_all/Figures/",cl,"_DAR_DEG_logFC.svg", sep = "")
       , height = 10, width = 14)

ggsave(combined_plot_short, file = paste("~/projects/combined_all/Figures/",cl,"_DAR_DEG_logFC.pdf", sep = "")
       , height = 10, width = 14)

ggsave(combined_plot, file = paste("~/projects/combined_all/Figures/",cl,"_DAR_DEG_logFC_mirror.svg", sep = "")
       , height = 10, width = 17)

ggsave(combined_plot, file = paste("~/projects/combined_all/Figures/",cl,"_DAR_DEG_logFC_mirror.pdf", sep = "")
       , height = 10, width = 17)

options(repr.plot.width=13, repr.plot.height=12)

library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)
 
rna_sig = rtab[which(rtab$p_val_adj<0.05), ]
rna_sig$brain_region = sapply(strsplit(as.character(rna_sig$id), "--"), `[`, 2)
colnames(rna_sig)[1] = "gene"
# Define function to calculate Jaccard similarity between two sets of genes
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Create a list of genes by "id"
genes_by_id <- rna_sig %>%
  group_by(id, brain_region) %>%
  summarise(genes = list(gene), .groups = 'drop')

# Generate a similarity matrix
similarity_matrix <- outer(
  genes_by_id$genes, 
  genes_by_id$genes, 
  Vectorize(jaccard_similarity)
)

# Convert the matrix to a data frame for easy plotting
rownames(similarity_matrix) <- genes_by_id$id
colnames(similarity_matrix) <- genes_by_id$id

# Create an annotation row for brain regions
annotation_row <- data.frame(
  brain_region = genes_by_id$brain_region
)
rownames(annotation_row) <- genes_by_id$id

# Plot using pheatmap
p = pheatmap(similarity_matrix,number_format = "%.1f" , 
         cluster_rows = TRUE,  # Hierarchical clustering of rows
         cluster_cols = TRUE,  # Hierarchical clustering of columns
         display_numbers = TRUE,annotation_row = annotation_row,  # Option to display similarity values
         main = paste("Jaccard Similarity of DEGs", cl))

options(repr.plot.width=13, repr.plot.height=12)

library(dplyr)
library(tidyr)
library(purrr)
library(pheatmap)
 
atac_sig = tab[which(tab$adjusted.p.value<0.05), ]
atac_sig$brain_region = sapply(strsplit(as.character(atac_sig$id), ":"), `[`, 2)
head(atac_sig)
colnames(atac_sig)[1] = "gene"
# Define function to calculate Jaccard similarity between two sets of genes
jaccard_similarity <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

# Create a list of genes by "id"
genes_by_id <- atac_sig %>%
  group_by(id, brain_region) %>%
  summarise(genes = list(gene), .groups = 'drop')

# Generate a similarity matrix
similarity_matrix <- outer(
  genes_by_id$genes, 
  genes_by_id$genes, 
  Vectorize(jaccard_similarity)
)

# Convert the matrix to a data frame for easy plotting
rownames(similarity_matrix) <- genes_by_id$id
colnames(similarity_matrix) <- genes_by_id$id

# Create an annotation row for brain regions
annotation_row <- data.frame(
  brain_region = genes_by_id$brain_region
)
rownames(annotation_row) <- genes_by_id$id

# Plot using pheatmap
p_dar = pheatmap(similarity_matrix,number_format = "%.1f" , 
         cluster_rows = TRUE,  # Hierarchical clustering of rows
         cluster_cols = TRUE,  # Hierarchical clustering of columns
         display_numbers = TRUE,annotation_row = annotation_row,  # Option to display similarity values
         main = paste("Jaccard Similarity of DARs", cl))

setwd("~/projects/combined_all/Figures/")

pdf(file = paste(cl , "_jaccard_DEG.pdf", sep = ""), height = 12, width =13)
p
dev.off()

pdf(file = paste(cl , "_jaccard_DAR.pdf", sep = ""), height = 12, width =13)
p_dar
dev.off()



##scratch 

ggsave(combined_plot, file = paste("~/projects/combined_all/Figures/",cl,"_DAR_DEG_logFC.svg", sep = "")
       , height = 9, width = 12)

ggsave(combined_plot, file = paste("~/projects/combined_all/Figures/",cl,"_DAR_DEG_logFC.pdf", sep = "")
       , height = 9, width = 12)

options(repr.plot.width=6, repr.plot.height=11)

dtab = tab[which(tab$adjusted.p.value<0.05),]
me =meta[grep(cl, meta$celltype_region),]
#ord= names(table(me$celltype_region))[order(table(me$celltype_region))]
#tab$id = factor(tab$id, levels = ord)
#dtab$id = factor(dtab$id, levels = ord)

#me$celltype_region = factor(me$celltype_region, levels = ord)

darp1 = ggplot(dtab, aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "DARs per Region",
       x = "Cell Type",
       y = "Number of DARs") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))
darp2 = ggplot(me, aes(x = celltype_region, fill= age)) +
  geom_bar(stat = "count", color = "black") +scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(title = "Cell # per Region",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))
darp3 = ggplot(tab[which(abs(tab$log2.fold_change.) > 0.15),], 
  aes(x = id, y = -`log2.fold_change.`, 
  color = ifelse(`adjusted.p.value` < 0.05, ifelse(`log2.fold_change.` < 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age DARs logFC",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") + ylim(c(-6,6))+
  theme_minimal() + coord_flip()+theme(text = element_text(color = "black", size = 13))
options(repr.plot.width=15, repr.plot.height=4)

#grid.arrange(darp2, darp1,darp3, ncol = 3)

RLP = dtab$feature.name[grep("RLP",dtab$id )]
AMY = dtab$feature.name[grep("AMY",dtab$id )]
CP = dtab$feature.name[grep(":CP",dtab$id )]
ENT = dtab$feature.name[grep("ENT",dtab$id )]
FC = dtab$feature.name[grep("FC",dtab$id )]
HCP = dtab$feature.name[grep("HCP",dtab$id )]
HCA = dtab$feature.name[grep("HCA",dtab$id )]
NAC = dtab$feature.name[grep("NAC",dtab$id )]

listInput <- list(RLP=RLP, AMY=AMY, CP=CP, ENT=ENT,FC=FC, HCP=HCP,HCA=HCA,NAC=NAC)
options(repr.plot.width=14, repr.plot.height=7)

darup = upset(fromList(listInput), order.by = "freq", nsets=8)
darup

numdars = ggplot(dtab, aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "DARs per Region",
       x = "Cell Type",
       y = "Number of DARs") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))+ facet_grid(celltype ~ ., space = "free",scales = "free_y")

numdars
