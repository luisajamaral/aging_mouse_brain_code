library(ggplot2)
library(data.table)
library(UpSetR)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

library(scales)

setwd("~/projects/combined_all/combined_DARs_redone/")

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

cl="_"

setwd("~/projects/combined_all/combined_DARs_redone/")
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
tab[which(tab$adjusted.p.value<0.05 & tab$log2.fold_change.< -0.25 ),"sig"] = "down"
tab[which(tab$adjusted.p.value<0.05 & tab$log2.fold_change.>0.25 ),"sig"] = "up"
tab$celltype = sapply(strsplit(as.character(tab$id), ":"), `[`, 1)
tab$region = sapply(strsplit(as.character(tab$id), ":"), `[`, 2)

setwd("~/projects/combined_all/female_RNA/DEG_results_latent_rep_mito_together///")
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
rtab[which(rtab$p_val_adj<0.05 & rtab$avg_log2FC<(-logfc) ),"sig"] = "down"
rtab[which(rtab$p_val_adj<0.05 & rtab$avg_log2FC>logfc ),"sig"] = "up"
#tab[which(tab$`pct.1`<0.05 |  tab$`pct.2`<0.05), "sig"] = "no"
rtab$Celltype= sapply(strsplit(as.character(rtab$id), "--"), `[`, 1)


significant_points <- rtab %>%
  filter(`p_val_adj` < 0.05)

non_significant_points <- rtab %>%
  filter(`p_val_adj` >= 0.05) %>%
  group_by(Celltype) %>%
  sample_n(min(n(), 1000))  # Sample max 1000 points per celltype

# Step 2: Combine the significant and sampled non-significant points
rsampled_tab <- bind_rows(significant_points, non_significant_points)

rsampled_tab$id = gsub("--", ":", rsampled_tab$id)


options(repr.plot.width=10, repr.plot.height=11)

significant_points <- tab %>%
  filter(`adjusted.p.value` < 0.05)

non_significant_points <- tab %>%
  filter(`adjusted.p.value` >= 0.05) %>%
  group_by(celltype) %>%
  sample_n(min(n(), 1000))  # Sample max 1000 points per celltype

# Step 2: Combine the significant and sampled non-significant points
sampled_tab <- bind_rows(significant_points, non_significant_points)





head(rsampled_tab)

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

keep



me =ameta
me$celltype_final = gsub("/", "-" , me$celltype_final)
me$celltype_final = gsub(" ", "_" , me$celltype_final)

me$celltype_region = gsub("/", "-" , me$celltype_region)
#keep = table(me$celltype_region)[which(table(me$celltype_region)>1000)]
#keep = names(keep)
me = me[which(me$celltype_final%in%keep),]
me = me[which(me$celltype_final%in%names(table(me$celltype_final)[which(table(me$celltype_final)>=5000)])),]
#me = me[which(me$celltype_final %in% c("Oligo NN", "OPC NN", "Microglia NN", "Astro-TE NN","Astro-NT NN" , "VLMC NN")),]
me$Modality_Sex = gsub(" ", "_", me$Modality_Sex)
ord =  names(table(me$celltype_final))[order(table(me$celltype_final), decreasing = T)]
me$celltype_final = factor(me$celltype_final, levels = ord)
head(me)



keep = unique(me$celltype_final)
keep

sampled_tab$celltype = factor(sampled_tab$celltype, levels = gsub(" ", "_", ord))
sampled_tab = sampled_tab[which(sampled_tab$id %in%keep), ]
rsampled_tab$Celltype = factor(rsampled_tab$Celltype, levels = gsub(" ", "_", ord))
rsampled_tab = rsampled_tab[which(rsampled_tab$id %in%keep), ]


#sampled_tab = sampled_tab[which(gsub("--", ":", sampled_tab$id) %in% unique(me$celltype_region)),]
#rsampled_tab = rsampled_tab[which(gsub("--", ":", rsampled_tab$id) %in% unique(me$celltype_region)),]


head(sampled_tab)

# First plot (g1)
g1 = ggplot(me, aes(x = celltype_final, fill = Modality_Sex)) +
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
deglogfc = ggplot(rsampled_tab, 
  aes(x = id, y = `avg_log2FC`, 
      color = sig)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "logFC Age-DEGs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-10, 10)) +  # Set the y-axis limits here
  theme_minimal() +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        strip.text.y = element_blank(),plot.title = element_text(hjust = 0.5),
        legend.position = "right",
        axis.title.y = element_blank(), axis.text.y = element_blank()
       ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
  facet_grid(Celltype ~ ., space = "free", scales = "free_y")

# Combine the two plots side by side
#combined_plot <- g1 + darlogfc + deglogfc +plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "right")

# Show the combined plot
#combined_plot



full_celltypes <- data.frame(id = keep)  # List of all cell types

# Step 1: Summarize the number of up and down DARs for each cell type
summary_tab <- tab %>%
  filter(`adjusted.p.value` < 0.05 & abs(`log2.fold_change.`)>0.25) %>%
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
  filter(`p_val_adj` < 0.05 & abs(`avg_log2FC`)>0.25) %>%
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

options(repr.plot.width=15, repr.plot.height=11)

combined_plot <- g1 + mirror_plot + darlogfc + rmirror_plot+ deglogfc +plot_layout(ncol = 5, guides = "collect") & theme(legend.position = "right")
combined_plot

head(tab)



full_celltypes <- data.frame(id = keep)  # List of all cell types

# Step 1: Summarize the number of up and down DARs for each cell type
summary_tab <- tab %>%
  filter(`adjusted.p.value` < 0.0001 & abs(`log2.fold_change.`)>0.25) %>%
  group_by(id, celltype,mark) %>%
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
  filter(`p_val_adj` < 0.05 & abs(`avg_log2FC`)>0.25) %>%
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

options(repr.plot.width=15, repr.plot.height=11)

combined_plot <- g1 + mirror_plot + darlogfc + rmirror_plot+ deglogfc +plot_layout(ncol = 5, guides = "collect") & theme(legend.position = "right")
combined_plot

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

library(org.Mm.eg.db)
library(viridis)
library(ggrepel)
library(AnnotationDbi)
library(biomaRt)

h3k9 = read.table("../../histone_volcano/H3K9me3_reprocessed_P0_forebrain.bed")


# Connect to Ensembl
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")  # Change to "hsapiens_gene_ensembl" for human
# Retrieve gene types for all genes in your dataset
  gene_annotations <- getBM(attributes = c("external_gene_name", "gene_biotype","transcript_length","transcript_gencode_basic", 
                                          "chromosome_name", "start_position", "end_position" , "percentage_gene_gc_content"
                                          ),
                          filters = "external_gene_name", values = unique(rsampled_tab$X),
                          mart = mart)
    

  # Filter for lncRNAs
  lncRNA_genes <- unique(gene_annotations$external_gene_name[gene_annotations$gene_biotype == "lncRNA"])
  protein_coding_genes <- unique(gene_annotations$external_gene_name[gene_annotations$gene_biotype == "protein_coding"])
  


tes = read.table("../../scTE/mm10.te.bed")
nrow(tes)
head(tes)

library(GenomicRanges)

# Convert h3k9 to a GRanges object
tes_gr <- GRanges(seqnames = tes$V1, id = tes$V4,
                   ranges = IRanges(start = tes$V2, end = tes$V3))

# Find overlaps
overlaps <- findOverlaps(gene_annotations_gr, tes_gr)

# Add mark column, default to NA
gene_annotations$te <- NA

# Assign "h3k9me3" to genes that overlap with any h3k9 region
gene_annotations$te[queryHits(overlaps)] <- "TE"

# View updated gene_annotations
head(gene_annotations)


library(GenomicRanges)

# Convert h3k9 to a GRanges object
h3k9_gr <- GRanges(seqnames = h3k9$V2,
                   ranges = IRanges(start = h3k9$V3, end = h3k9$V4))

# Convert gene_annotations to a GRanges object
gene_annotations_gr <- GRanges(seqnames = paste0("chr", gene_annotations$chromosome_name),
                               ranges = IRanges(start = gene_annotations$start_position, 
                                                end = gene_annotations$end_position))

# Find overlaps
overlaps <- findOverlaps(gene_annotations_gr, h3k9_gr)

# Add mark column, default to NA
gene_annotations$mark <- NA

# Assign "h3k9me3" to genes that overlap with any h3k9 region
gene_annotations$mark[queryHits(overlaps)] <- "h3k9me3"

# View updated gene_annotations
head(gene_annotations)


h3k9_genes <- gene_annotations$external_gene_name[gene_annotations$mark == "h3k9me3"]


rsampled_tab$sig[which(rsampled_tab$p_val_adj>1e-10)] = "no"
rsampled_tab$mark = "" 
rsampled_tab$mark[which(rsampled_tab$X %in% h3k9_genes)]= "H3K9me3"
rsampled_tab$sig_mark = paste(rsampled_tab$sig, rsampled_tab$mark)
rsampled_tab$sig_mark[which(rsampled_tab$sig == "no") ] = "no"

rsampled_tab$type = "Other" 
rsampled_tab$type[which(rsampled_tab$X %in% protein_coding_genes)]= "protein_coding_genes"
rsampled_tab$type[which(rsampled_tab$X %in% lncRNA_genes)]= "lncRNA_genes"
rsampled_tab$type[which(rsampled_tab$sig == "no")]= "not significant"

rsampled_tab$type = factor(rsampled_tab$type, levels = c("not significant" , "Other", "protein_coding_genes", "lncRNA_genes"))

library(dplyr)
library(ggplot2)
library(ggrepel)

# Identify top upregulated and downregulated genes per id
cur = rsampled_tab#[-grep("Rik", rsampled_tab$X),]
#cur = cur[-grep("Gm", cur$X),]
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



rsampled_tab = rsampled_tab[order(rsampled_tab$type),]

# Plot with text labels for top up and downregulated genes
deglogfc <- ggplot(rsampled_tab, 
  aes(x = id, y = `avg_log2FC`, 
      color = type)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.9) +
  scale_color_manual(values = c("not significant" =  "grey", "Other" =  "pink", "protein_coding_genes" =  "darkgrey", "lncRNA_genes" = "red")) +
  labs(title = "logFC Age-DEGs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-10, 10)) +  # Set the y-axis limits here
  theme_minimal() +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        strip.text.y = element_blank(), plot.title = element_text(hjust = 0.5),
       ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
  facet_grid(Celltype ~ ., space = "free", scales = "free_y") +
  geom_text_repel(data = label_data, aes(label = X), 
                  size = 3, nudge_x = 0, direction = "y", max.overlaps = Inf)

options(repr.plot.width=7, repr.plot.height=10)

deglogfc

library(dplyr)
library(ggplot2)
library(ggrepel)

rsampled_tab$sig = "no"
rsampled_tab$sig[which(rsampled_tab$p_val_adj<1e-5 & rsampled_tab$avg_log2FC>.5)] = "up"
rsampled_tab$sig[which(rsampled_tab$p_val_adj<1e-5 & rsampled_tab$avg_log2FC< -.5)] = "down"

rsampled_tab$mark = "lncRNA" 
rsampled_tab$mark[which(rsampled_tab$X %in% lncRNA_genes)]= "lncRNA"
rsampled_tab$mark[which(rsampled_tab$X %in% protein_coding_genes)]= "protein_coding"

#rsampled_tab$mark[which(rsampled_tab$sig == "no") ] = "not significant"
rsampled_tab$mark=factor(rsampled_tab$mark, levels = c("not significant" , "not lncRNA","protein_coding", "lncRNA"))
# Identify top upregulated and downregulated genes per id
cur = rsampled_tab#[-grep("Rik", rsampled_tab$X),]
#cur = cur[-grep("Gm", cur$X),]
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

#rsampled_tab = rsampled_tab[order(rsampled_tab$mark,decreasing = F),]
#rsampled_tab = rsampled_tab[sample(nrow(rsampled_tab)), ]

# Plot with text labels for top up and downregulated genes
deglogfc <- ggplot(rsampled_tab, 
  aes(x = id, y = `avg_log2FC`, 
      color = mark)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.6) +
  #scale_color_manual(values = c("not significant" =  "grey", "not H3K9me3" = "steelblue1", "H3K9me3" = "red")) +
  scale_color_manual(values = c("grey", "red")) +

  labs(title = "logFC Age-DEGs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-10, 10)) +  # Set the y-axis limits here
  theme_minimal() +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        strip.text.y = element_blank(), plot.title = element_text(hjust = 0.5),
       ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
 # facet_grid(Celltype ~ mark, space = "free", scales = "free_y") +
  geom_text_repel(data = label_data, aes(label = X), 
                  size = 3, nudge_x = 0, direction = "y", max.overlaps = Inf)

options(repr.plot.width=7, repr.plot.height=10)

deglogfc

library(dplyr)
library(ggplot2)
library(ggrepel)

rsampled_tab$sig = "no"
rsampled_tab$sig[which(rsampled_tab$p_val_adj<1e-5 & rsampled_tab$avg_log2FC>0.25)] = "up"
rsampled_tab$sig[which(rsampled_tab$p_val_adj<1e-5 & rsampled_tab$avg_log2FC< -0.25)] = "down"

rsampled_tab$mark = "not H3K9me3" 
rsampled_tab$mark[which(rsampled_tab$X %in% h3k9_genes)]= "H3K9me3"
#rsampled_tab$mark[which(rsampled_tab$sig == "no") ] = "not significant"
rsampled_tab$mark=factor(rsampled_tab$mark, levels = c("not significant" , "not H3K9me3", "H3K9me3"))
# Identify top upregulated and downregulated genes per id
cur = rsampled_tab#[-grep("Rik", rsampled_tab$X),]
#cur = cur[-grep("Gm", cur$X),]
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

#rsampled_tab = rsampled_tab[order(rsampled_tab$mark,decreasing = F),]
#rsampled_tab = rsampled_tab[sample(nrow(rsampled_tab)), ]

# Plot with text labels for top up and downregulated genes
deglogfc <- ggplot(rsampled_tab, 
  aes(x = id, y = `avg_log2FC`, 
      color = mark)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.6) +
  #scale_color_manual(values = c("not significant" =  "grey", "not H3K9me3" = "steelblue1", "H3K9me3" = "red")) +
  scale_color_manual(values = c("grey", "red")) +

  labs(title = "logFC Age-DEGs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-5, 5)) +  # Set the y-axis limits here
  theme_minimal() +
  coord_flip() +
  theme(text = element_text(color = "black", size = 13),
        strip.text.y = element_blank(), plot.title = element_text(hjust = 0.5),
       ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
#  facet_grid(Celltype ~ mark, space = "free", scales = "free_y") +
  geom_text_repel(data = label_data, aes(label = X), 
                  size = 3, nudge_x = 0, direction = "y", max.overlaps = Inf)

options(repr.plot.width=7, repr.plot.height=10)

deglogfc

head(sampled_tab)


tes_gr[subjectHits(overlaps)[1:5]]


tes_gr <- GRanges(seqnames = tes$V1, id = tes$V4,
                   ranges = IRanges(start = tes$V2, end = tes$V3))

peaks_gr <- GRanges(seqnames = sampled_tab$chromosome, id = sampled_tab$feature.name,
                   ranges = IRanges(start = sampled_tab$start, end = sampled_tab$end))


# Find overlaps
overlaps <- findOverlaps(peaks_gr, tes_gr)

# Add mark column, default to NA
sampled_tab$te <- NA

# Assign "h3k9me3" to genes that overlap with any h3k9 region
sampled_tab$te[queryHits(overlaps)] <- "TE"

# View updated gene_annotations
head(sampled_tab)




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
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
  facet_grid( ~ te , space = "free", scales = "free_y")


darlogfc

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
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") + 
  facet_grid( ~ mark , space = "free", scales = "free_y")


darlogfc

# Second plot (darlogfc)
darlogfc = ggplot(sampled_tab, 
  aes(x = id, y = `log2.fold_change.`, 
      color = mark)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("red", "grey")) +
  labs(title = "logFC Age-DARs", x = "Cell Type", y = "log2 Fold Change", color = "adj p < 0.05") +
  scale_y_continuous(limits = c(-10, 10)) +  # Set the y-axis limits here
  theme_minimal() +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black")# + 
  #facet_grid( ~ mark , space = "free", scales = "free_y")


darlogfc


