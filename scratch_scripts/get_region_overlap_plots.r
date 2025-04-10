library(ggplot2)
library(data.table)
library(UpSetR)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("~/projects/combined_all/region_DARs_redone/")

length(unique(rmeta$celltype_final))
nrow(rmeta)

age_color = read.csv("../Figures/color_scheme/Age.csv")
color_vector <- setNames(age_color$Color, age_color$Age)


meta = read.csv("../final_meta.csv")
meta$celltype_final = gsub("Lymphoid","Lymphoid NN", meta$celltype_final)
meta = meta[-which(meta$celltype_final=="doublet"),]
meta$celltype_final=paste(meta$celltype_final)
meta$celltype_final[which(meta$best_celltype_fixed=="IOL")] = "IOL NN"




meta$celltype_region = paste(meta$celltype_final, meta$region, sep = ":")
#meta$celltype_region=gsub("Astro-TE", "Astro", meta$celltype_region)
#meta$celltype_region=gsub("Astro-NT", "Astro", meta$celltype_region)
meta$celltype_region=gsub(" ", "_", meta$celltype_region)
meta$age = factor(meta$age , levels = c("2mo", "9mo", "18mo"))


meta$celltype_region_modality = paste(meta[,'celltype_region'], "ATAC", sep = "_")
rmeta$celltype_region_modality = paste(rmeta[,'celltype_region'], "RNA", sep = "_")

ameta = rbind(meta[,c('sample', 'celltype_region_modality','celltype_region', 'celltype_final', 'region')], rmeta[,c('sample', 'celltype_region_modality','celltype_region','celltype_final', 'region')])
ameta$celltype_region_modality = gsub("--", ":", ameta$celltype_region_modality)
ameta$celltype_region = gsub("--", ":", ameta$celltype_region)

ameta$Modality = "ATAC"
ameta$Modality[grep("RNA", ameta$celltype_region_modality)] = "RNA"

head(ameta)
tail(ameta)

#setwd("../../region_DARs_redone/")
options(repr.plot.width=14, repr.plot.height=5)
cl = "NN"
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
tab[which(tab$adjusted.p.value<0.05 & tab$log2.fold_change.<0 ),"sig"] = "up"
tab[which(tab$adjusted.p.value<0.05 & tab$log2.fold_change.>0 ),"sig"] = "down"


head(me)
me =ameta[grep(cl,ameta$celltype_final),]
me = me[which(me$celltype_final %in% c("Oligo NN", "OPC NN", "Microglia NN", "Astro-TE NN","Astro-NT NN" , "VLMC NN")),]

ggplot(me, aes(x = celltype_region, fill= Modality)) +
  geom_bar(stat = "count", color = "black", position = "dodge") +#scale_fill_manual(values = color_vector) +
  theme_minimal() +
  labs(title = "Cell # per Region",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))+ facet_grid(celltype_final ~ ., space = "free",scales = "free_y")

options(repr.plot.width=5, repr.plot.height=10)

ggplot(dtab, aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +  scale_fill_manual(values = c("blue", "red")) +
  labs(title = "DARs per Region",
       x = "Cell Type",
       y = "Number of DARs") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))+ facet_grid(celltype ~ ., space = "free",scales = "free_y")



options(repr.plot.width=5, repr.plot.height=10)

ggplot(tab[which(abs(tab$log2.fold_change.) > 0.75),], 
  aes(x = id, y = -`log2.fold_change.`, 
  color = ifelse(`adjusted.p.value` < 0.05, ifelse(`log2.fold_change.` < 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age DARs logFC",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") + ylim(c(-6,6))+
  theme_minimal() + coord_flip()+theme(text = element_text(color = "black", size = 13))+ facet_grid(celltype ~ ., space = "free",scales = "free_y")

options(repr.plot.width=5, repr.plot.height=10)

ggplot(tab[which(abs(tab$avg_log2FC) > 0.25),], 
  aes(x = id, y = -`avg_log2FC`, 
  color = ifelse(`p_val_adj` < 0.05, ifelse(`avg_log2FC` < 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Genes",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip()+theme(text = element_text(color = "black", size = 13))+ facet_grid(Celltype ~ ., space = "free",scales = "free_y")

head(tab)

tab$celltype = sapply(strsplit(as.character(tab$id), ":"), `[`, 1)
tab$region = sapply(strsplit(as.character(tab$id), ":"), `[`, 2)
tab = tab[which(tab$celltype %in% c("Oligo_NN", "OPC_NN", "Microglia_NN", "Astro-TE_NN","Astro-NT_NN" , "VLMC_NN")),]

dtab = tab[which(tab$adjusted.p.value<0.05),]
me =meta[grep(cl, meta$celltype_region),]
#ord= names(table(me$celltype_region))[order(table(me$celltype_region))]
#tab$id = factor(tab$id, levels = ord)
#dtab$id = factor(dtab$id, levels = ord)

me$celltype_region = factor(me$celltype_region, levels = ord)

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

grid.arrange(darp2, darp1,darp3, ncol = 3)

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

darup

# Set the directory containing the result files
result_dir <- "~/projects/combined_all/region_DARs_redone/"
# Load all result files
files <- list.files(result_dir, "knownResults.txt",full.names = F, recursive = T)
files = files[grep(cl,files)]
files

df = fread(files[1])
head(df)


# Identify the top 10 upregulated and top 10 downregulated pathways by lowest p-value for any cell type
top_upregulated <- combined_df %>%
  filter(`q-value ` < 0.05 & `P-value` <0.001 & Direction=="Up") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 30) %>%
  ungroup()

top_downregulated <- combined_df %>%
    filter(`q-value ` < 0.05 & `P-value` <0.001 & Direction=="Down") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 30) %>%
  ungroup()

# Combine the top upregulated and downregulated pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated)

unique_pathways <- unique(top_pathways$Motif)
unique_pathways




# Filter the combined data frame to include all NES and p-values for the selected pathways
plot_df <- combined_df %>%
  filter(Motif %in% unique_pathways)

# Create the dotplot
options(repr.plot.width=18, repr.plot.height=8)
plot_df$Description <- substr(plot_df$Description, 1, 50)

plot_df$significance_direction = plot_df$`Log P-value`
plot_df$significance_direction[which(plot_df$Direction=="Up")] = plot_df$`Log P-value`*-1
gsea = ggplot(plot_df, aes(x = reorder(Motif, `Log P-value`, decreasing = T), y = CellType,  size = -`Log P-value`,color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Significance x Direction") +
  scale_size_continuous(range = c(3, 10), name = "-log10(Adjusted p-value)") +
  theme_minimal() +
  xlab("Cell Type") +
  ylab("Pathway") +
  ggtitle(paste("Age - Motifs",cl)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        axis.text.y = element_text(size = 8),legend.position = "left")
gsea

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

# Set the directory containing the result files
result_dir <- "~/projects/combined_all/region_DARs_redone/"
# Load all result files
files <- list.files(result_dir, "knownResults.txt",full.names = F, recursive = T)
files = files[grep(cl,files)]
files

# Initialize an empty list to store data frames
all_results <- list()

# Loop through each file and read the data
for (file in files) {
  # Extract the cell type from the file name
  celltype <- gsub("_motifs/knownResults.txt", "", (file))
  celltype <- gsub("_peaks", "", (celltype))

  # Read the data
  df <- fread(file)
  if(nrow(df)<1){
      next
  }
  # Add cell type as a column
  df$CellType <- celltype
  df$Direction = "Down"
  if(length(grep("up", file))>0){
    df$Direction = "Up"
  }
  colnames(df)  = sapply(strsplit(as.character(colnames(df)), "[(]"), `[`, 1)

  # Store the results
  all_results[[celltype]] <- df
}

# Combine all data frames into one
combined_df <- bind_rows(all_results)
combined_df$Motif = sapply(strsplit(as.character(combined_df$`Motif Name`), "/"), `[`, 1)

# Identify the top 10 upregulated and top 10 downregulated pathways by lowest p-value for any cell type
top_upregulated <- combined_df %>%
  filter(`q-value ` < 0.0001 & `P-value` <0.00001 & Direction=="Up") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 10) %>%
  ungroup()

top_downregulated <- combined_df %>%
    filter(`q-value ` < 0.0001 & `P-value` <0.00001 & Direction=="Down") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 10) %>%
  ungroup()

# Combine the top upregulated and downregulated pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated)

unique_pathways <- unique(top_pathways$Motif)
plot_df <- combined_df %>%
  filter(Motif %in% unique_pathways)

# Perform hierarchical clustering on the motifs based on Log P-value
# First, create a matrix of Log P-values with motifs as rows and cell types as columns
motif_matrix <- plot_df %>%
  select(Motif, CellType, `Log P-value`) %>%
  spread(CellType, `Log P-value`, fill = 0)

# Convert the Motif column to row names (base R way)
motif_matrix <- as.data.frame(motif_matrix)
rownames(motif_matrix) <- motif_matrix$Motif
motif_matrix$Motif <- NULL  # Remove the Motif column after assigning it to rownames

# Compute the distance matrix and hierarchical clustering
motif_dist <- dist(motif_matrix)
motif_clust <- hclust(motif_dist)

# Reorder motifs based on the clustering
plot_df$Motif <- factor(plot_df$Motif, levels = rownames(motif_matrix)[motif_clust$order])

# Create the dotplot with clustered motifs
options(repr.plot.width=18, repr.plot.height=8)
plot_df$Description <- substr(plot_df$Description, 1, 50)

plot_df$significance_direction = plot_df$`Log P-value`
plot_df$significance_direction[which(plot_df$Direction=="Up")] = plot_df$`Log P-value`[which(plot_df$Direction=="Up")]*-1

plot_df$significance_direction[which(plot_df$significance_direction>100)] = 100
plot_df$significance_direction[which(plot_df$significance_direction< -100)] = -100


motifs= ggplot(plot_df, aes(x = Motif, y = CellType, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, c(-100,100), name = "Significance x Direction") +
  scale_size_continuous(range = c(3, 10), name = "-log10(Adjusted p-value)") +
  theme_minimal() +
  xlab("Motif (Clustered)") +
  ylab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        axis.text.y = element_text(size = 8), legend.position = "left") + facet_wrap(~Direction,nrow = 2, scales = "free_y")
motifs


rmeta = read.table("../Figures/rna_meta.csv")

rmeta$celltype_region = paste(rmeta$celltype_final, "--", rmeta$region,sep = "")
rmeta$celltype_region = gsub(" " , "_", rmeta$celltype_region )
head(rmeta)


files = list.files(".", "NN")
files = files[grep(".csv", files)]
files

setwd("~/projects/combined_all/female_RNA/DEG_results_latent_rep_mito/")
cl = "NN"
logfc = 0
files = list.files(".", cl)
files = files[grep(".csv", files)]
all = list()
for(f in files){
    curr = read.csv(f)
    id = gsub(".csv" , "", f)
    curr$id = id 
    all[[id]] = curr
    }
tab = do.call(rbind, all)
tab=as.data.frame(tab)
tab$sig = "no"
tab[which(tab$p_val_adj<0.05 & tab$avg_log2FC<(-logfc) ),"sig"] = "up"
tab[which(tab$p_val_adj<0.05 & tab$avg_log2FC>logfc ),"sig"] = "down"
#tab[which(tab$`pct.1`<0.05 |  tab$`pct.2`<0.05), "sig"] = "no"
tab$Celltype= sapply(strsplit(as.character(tab$id), "--"), `[`, 1)
tab = tab[which(tab$Celltype %in% c("Oligo_NN", "OPC_NN", "Microglia_NN", "Astro-TE_NN","Astro-NT_NN" , "VLMC_NN")),]

dtab = tab[which(tab$p_val_adj<0.05 & abs(tab$avg_log2FC)>logfc),]
me =rmeta[grep(cl, rmeta$celltype_region),]
bc = names(table(me$celltype_region)[which(table(me$celltype_region)>100)])
me = me[which(me$celltype_region%in%bc),]
#ord= names(table(me$celltype_region))[order(table(me$celltype_region))]
#tab$id = factor(tab$id, levels = ord)
#dtab$id = factor(dtab$id, levels = ord)
#me$celltype_region = factor(me$celltype_region, levels = ord)

p1 = ggplot(dtab, aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = c("blue", "red")) +
  theme_minimal() +  

  labs(title = "DEGs per Region",
       x = "Cell Type",
       y = "Number of DEGs") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))

me$age = factor(me$age , levels = c("2mo", "9mo", "18mo"))

p2 = ggplot(me, aes(x = celltype_region, fill= age)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +scale_fill_manual(values = color_vector) +
  labs(title = "Cell # per Region",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))
p3 = ggplot(tab[which(abs(tab$avg_log2FC) > 0.15),], 
  aes(x = id, y = -`avg_log2FC`, 
  color = ifelse(`p_val_adj` < 0.01, ifelse(`avg_log2FC` < 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Genes",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip()+theme(text = element_text(color = "black", size = 13))
options(repr.plot.width=15, repr.plot.height=4)

grid.arrange(p2, p1,p3, ncol = 3)

RLP = dtab$X[grep("RLP",dtab$id)]
AMY = dtab$X[grep("AMY",dtab$id)]
CP = dtab$X[grep("CP",dtab$id)]
ENT = dtab$X[grep("ENT",dtab$id)]
FC = dtab$X[grep("FC",dtab$id)]
HCP = dtab$X[grep("HCP",dtab$id)]
HCA = dtab$X[grep("HCA",dtab$id)]
NAC = dtab$X[grep("NAC",dtab$id)]
options(repr.plot.width=12, repr.plot.height=7)
listInput <- list(RLP=RLP, AMY=AMY, CP=CP, ENT=ENT,FC=FC, HCP=HCP,HCA=HCA,NAC=NAC)
upset(fromList(listInput), order.by = "freq", nsets=8, text.scale = 1.5)



options(repr.plot.width=15, repr.plot.height=6)

grid.arrange(p2, p1,p3, ncol = 3)




files

setwd("~/projects/combined_all/female_RNA/DEG_results_latent_rep_mito/")
# Set the directory containing the result files
result_dir <- "."

# Load all result files
files <- list.files(result_dir, "GSEA_table",full.names = F)

files=files[grep(cl, files)]


# Initialize an empty list to store data frames
all_results <- list()

# Loop through each file and read the data
for (file in files) {
  # Extract the cell type from the file name
  celltype <- gsub("_GSEA_table.txt", "", (file))
  
  # Read the data
  df <- read.table(file)
  if(nrow(df)<1){
      next
  }
  # Add cell type as a column
  df$CellType <- celltype
  
  # Store the results
  all_results[[celltype]] <- df
}

# Combine all data frames into one
combined_df <- bind_rows(all_results)

# Identify the top 10 upregulated and top 10 downregulated pathways by lowest p-value for any cell type
top_upregulated <- combined_df %>%
  filter(NES > 1 & p.adjust<0.0001 ) %>%
  group_by(CellType) %>%
  arrange(p.adjust) %>%
  slice_head(n = 30) %>%
  ungroup()

top_upregulated = rbind(top_upregulated, combined_df[which(combined_df$Description=="cytokine binding"),],
                       combined_df[which(combined_df$Description=="oxidative phosphorylation"),])

top_downregulated <- combined_df %>%
  filter(NES < -1& p.adjust<0.0001
        ) %>%
  group_by(CellType) %>%
  arrange(p.adjust) %>%
  slice_head(n = 30) %>%
  ungroup()

# Combine the top upregulated and downregulated pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated)

# Get the unique pathways from the top lists
unique_pathways <- unique(top_pathways$Description)

# Filter the combined data frame to include all NES and p-values for the selected pathways
plot_df <- combined_df %>%
  filter(Description %in% unique_pathways)

# Create the dotplot
options(repr.plot.width=12, repr.plot.height=8)
plot_df$Description <- substr(plot_df$Description, 1, 50)

gsea = ggplot(plot_df, aes(x = reorder(Description, NES, decreasing = T), y = CellType, size = -log10(p.adjust), color = NES)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
  scale_size_continuous(range = c(3, 10), name = "-log10(Adjusted p-value)") +
  theme_minimal() +
  xlab("Cell Type") +
  ylab("Pathway") +
  ggtitle(paste("Age - GSEA",cl)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        axis.text.y = element_text(size = 8),legend.position = "left")

options(repr.plot.width=12, repr.plot.height=4)


grid.arrange(p2, p1,p3, ncol = 3)
options(repr.plot.width=12, repr.plot.height=6)


upset(fromList(listInput), order.by = "freq", nsets=8, text.scale = 1)
options(repr.plot.width=12, repr.plot.height=6)
plot_df$CellType = factor(plot_df$CellType, levels = ord)
gsea
options(repr.plot.width=12, repr.plot.height=6)
grid.arrange(darp2, darp1,darp3, ncol = 3)
darup

options(repr.plot.width=18, repr.plot.height=6)
motifs



