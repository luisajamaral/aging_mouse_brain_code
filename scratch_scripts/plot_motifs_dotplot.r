library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
setwd("~/ps-renlab2/projects/combined_all/Figures")




rmeta = read.table("rna_meta.csv")
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

cl = ""

me =ameta[grep(cl,ameta$celltype_final),]
me$celltype_final = gsub("/", "-" , me$celltype_final)
me$celltype_region = gsub("/", "-" , me$celltype_region)
ord =  names(table(me$celltype_final))[order(table(me$celltype_final), decreasing = T)]
ord

input_motifs  =c("AP-1(bZIP)", "Fra1(bZIP)","CTCF(Zf)","BORIS(Zf)","RFX(HTH)", "Olig2(bHLH)", "NeuroG2(bHLH)", "Mef2b(MADS)" , "Tgif1(Homeobox)", "Sox10(HMG)", "NRF(NRF)", 
          "Sp5(Zf)","IRF1(IRF)","ETS(ETS)","Klf4(Zf)","NFY(CCAAT)","Eomes(T-box)","PU.1:IRF8(ETS:IRF)", "Lhx1(Homeobox)","Tbx5(T-box)",
          "Snail1(Zf)","Sp1(Zf)","Foxo3(Forkhead)", "Twist2(bHLH)","CREB5(bZIP)","Oct6(POU,Homeobox)","KLF14(Zf)","Sp1(Zf)","Sox8(HMG)","Sox6(HMG)","Sox3(HMG)","NeuroG2(bHLH)","YY1(Zf)","E2A(bHLH),near_PU.1",
           "Rfx1(HTH)", "X-box(HTH)", "Slug(Zf)", "BHLHA15(bHLH)", "Pbx3(Homeobox)", "Mef2c(MADS)","Mef2a(MADS)","Lhx2(Homeobox)", "KLF5(Zf)","Nrf2(bZIP)", "Nrf3(bZIP)","NRF3(bZIP)")


setwd("~/ps-renlab2/projects/combined_all/region_DARs_redone/motifs_genOn_GO_results/")

#setwd("~/ps-renlab2/projects/combined_all/region_DARs_redone/motifs_genOn_GO_results/")
setwd("~/ps-renlab2/projects/combined_all/combined_DARs_redone/motifs/")

# Set the directory containing the result files
#result_dir <- "~/ps-renlab2/projects/combined_all/region_DARs_redone/motifs_genOn_GO_results/"
result_dir <- "~/ps-renlab2/projects/combined_all/combined_DARs_redone/motifs/"

# Load all result files
files <- list.files(result_dir, "knownResults.txt",full.names = F, recursive = T)
files = files[grep(cl,files)]
#files = files[grep("NN",files)]
#files = files[-grep("IMN",files)]
#files = files[-grep("IOL",files)]


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
combined_df <- do.call(rbind, all_results)
combined_df$Motif = sapply(strsplit(as.character(combined_df$`Motif Name`), "/"), `[`, 1)
combined_df$`P-value`= as.numeric(combined_df$`P-value`)
combined_df$`q-value`= as.numeric(combined_df$`q-value`)

# Identify the top 10 upregulated and top 10 downregulated pathways by lowest p-value for any cell type
top_upregulated <- combined_df %>%
  filter(`q-value ` < 1e-10 & `P-value` < 1e-3  & Direction=="Up") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 58) %>%
  ungroup()

top_downregulated <- combined_df %>%
    filter(`q-value ` < 1e-10 & `P-value` < 1e-3 & Direction=="Down") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n =58) %>%
  ungroup()

# Combine the top upregulated and downregulated pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated)

unique_pathways <- unique(top_pathways$Motif)
unique_pathways = unique(c(unique_pathways, input_motifs))
#unique_pathways = input_motifs[which(input_motifs %in%unique(top_pathways$Motif))]
plot_df <- combined_df %>%
  filter(Motif %in% unique_pathways)

plot_df <- plot_df %>%
  distinct(Motif, CellType, .keep_all = TRUE)

plot_df$significance_direction = plot_df$`Log P-value`
plot_df$significance_direction[which(plot_df$Direction=="Up")] = plot_df$`Log P-value`[which(plot_df$Direction=="Up")]*-1

plot_df$significance_direction[which(plot_df$significance_direction>100)] = 100
plot_df$significance_direction[which(plot_df$significance_direction< -100)] = -100


# Perform hierarchical clustering on the motifs based on Log P-value
# First, create a matrix of Log P-values with motifs as rows and cell types as columns
motif_matrix <- plot_df %>%
  select(Motif, CellType, significance_direction) %>%
  spread(CellType, significance_direction, fill = 0)

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
options(repr.plot.width=8, repr.plot.height=5.5)

plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)
plot_df$ct = gsub("up_", "", plot_df$ct )
plot_df$ct = gsub("down_", "", plot_df$ct )



-log(0.05)


plot_df$CellType = gsub("up_", "", plot_df$CellType)
plot_df$CellType = gsub("down_", "", plot_df$CellType)

#plot_df = plot_df %>%
#  group_by(CellType) %>%
#  filter(any(`q-value ` < 0.1)) %>%  # Keep only CellTypes with q-value < 0.05
#  ungroup()



plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)

plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))

plot_df = plot_df[which(plot_df$`Log P-value`< -3),]


options(repr.plot.width=20, repr.plot.height=10
       )
plot_df$significance_direction[which(plot_df$significance_direction>10)] = 10
plot_df$significance_direction[which(plot_df$significance_direction< -10)] = -10
plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))

plot_df = plot_df[order(plot_df$`Log P-value`),]

motifs= ggplot(plot_df, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,  name = "Significance x Direction") +
  scale_size_continuous( name = "-log10(p-value)" ,range = c(1,5)) +
  theme_minimal() +
  ylab("Motif") +
  xlab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10), strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), legend.position = "right") + coord_flip() +facet_grid(ct~., scales = "free" ,space = "free")
motifs

#ggsave(motifs, height = 5.5, width = 8.5, file = "~/projects/combined_all/Figures/Figure4-Gaba/Gaba_Motifs.pdf")
#ggsave(motifs, height = 6.5, width = 10, file = "~/projects/combined_all/Figures/Figure3-NN/NN_Motifs.pdf")
#ggsave(motifs, height = 6.5, width = 10, file = "~/projects/combined_all/Figures/Figure3-NN/NN_Motifs.svg")


ggsave(motifs, height = 10.5, width = 10, file = "~/projects/combined_all/Figures/Figure5-Glut/Glut_Motifs_vertical.pdf")


options(repr.plot.width=9, repr.plot.height=9
       )


sub = plot_df[which(plot_df$Motif %in% input_motifs ), ]
sub$clade= sapply(strsplit(as.character(sub$CellType), ":"), `[`, 2)


motifs = ggplot(sub, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,  name = "Significance\n+ up in aging\n- down in aging") +
  scale_size_continuous( name = "-log10(p-value)" ,range = c(1,7)) +
  theme_minimal() +
  ylab("Motif") +
  xlab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10), strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), legend.position = "right") + coord_flip() +facet_grid(clade~., scales = "free" ,space = "free")
motifs
#ggsave(motifs, height = 5.5, width = 8.5, file = "~/projects/combined_all/Figures/Figure4-Gaba/Gaba_Motifs.pdf")
#ggsave(motifs, height = 6.5, width = 10, file = "~/projects/combined_all/Figures/Figure3-NN/NN_Motifs.pdf")
#ggsave(motifs, height = 6.5, width = 10, file = "~/projects/combined_all/Figures/Figure3-NN/NN_Motifs.svg")




ggsave(motifs, height = 9, width = 9, file = "~/projects/combined_all/Figures/Figure5-Glut/Glut_Motifs_small_orderclade.pdf")

ggsave(motifs, height = 9, width = 12, file = "~/projects/combined_all/Figures/Figure5-Glut/Glut_Motifs.pdf")

ggsave(motifs, height = 13, width = 18, file = "~/projects/combined_all/Figures/Figure3/NN_Motifs.pdf")

plot_df$CellType = gsub("up_", "", plot_df$CellType)
plot_df$CellType = gsub("down_", "", plot_df$CellType)
plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)

plot_df = plot_df %>%
  group_by(CellType) %>%
  filter(any(`q-value ` < 0.05)) %>%  # Keep only CellTypes with q-value < 0.05
  ungroup()
options(repr.plot.width=11, repr.plot.height=6)
motifs = ggplot(plot_df, aes(x = Motif, y = factor(CellType, levels = rev(unique(CellType))), fill = significance_direction)) +
  geom_tile() +
  #scale_fill_gradient(low = "blue", high = "red", name = "-log10(p-value)") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Significance x Direction") +
  theme_minimal() +
  ylab("Cell Type") +
  xlab("Motif") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        axis.text.y = element_text(size = 8), legend.position = "right") + 
  facet_grid(Direction~., scales = "free", space = "free",)

#ggsave(motifs, height = 6, width = 11, file = "~/projects/combined_all/Figures/Figure5/Gaba_Motifs.svg")
#ggsave(motifs, height = 6, width = 11, file = "~/projects/combined_all/Figures/Figure5/Gaba_Motifs.pdf")


options(repr.plot.width=15, repr.plot.height=9)

motifs


ggplot(plot_df, aes(x = Motif, y = factor(CellType, levels = rev(unique(CellType))), fill = significance_direction)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Significance x Direction") +
  theme_minimal() +
  ylab("Cell Type") +
  xlab("Motif") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(
    axis.text.x = element_text(angle = 65, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    strip.text.y = element_blank(),  # Hides repeated row facet labels
    strip.placement = "outside"      # Places facet labels outside the panel
  ) + 
  facet_grid(~Direction+ct , scales = "free", space = "free")

levels(plot_df$CellType)

options(repr.plot.width=6, repr.plot.height=14)

plot_df$significance_direction[which(plot_df$significance_direction>20)] = 20
plot_df$significance_direction[which(plot_df$significance_direction< -20)] = -20


motifs= ggplot(plot_df, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, c(-100,100), name = "Significance x Direction") +
  scale_size_continuous( name = "-log10(Adjusted p-value)") +
  theme_minimal() +
  xlab("Motif (Clustered)") +
  ylab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        axis.text.y = element_text(size = 8), legend.position = "left",text = element_text(family = "Arial")) + facet_wrap(~Direction,nrow = 1, scales = "free_x")
motifs


pdf("~/projects/combined_all/Figures/Figure3/Oligo_Motifs_big.pdf", height = 6, width = 19)
ggplot(plot_df, aes(x = Motif, y = CellType, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, c(-100,100), name = "Significance x Direction") +
  scale_size_continuous(range = c(3, 10), name = "-log10(Adjusted p-value)") +
  theme_minimal() +
  xlab("Motif (Clustered)") +
  ylab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1),
        axis.text.y = element_text(size = 8), legend.position = "left") + facet_wrap(~Direction,nrow = 2, scales = "free_y")
dev.off()

ggsave(motifs, height = 6, width = 19, file = "~/projects/combined_all/Figures/Figure3/Oligo_Motifs_big.svg")



library(svglite)


install.packages("svglite")

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(dendextend)

cl = "Oligo"

setwd("~/projects/combined_all/region_DARs_redone/")
result_dir <- "~/projects/combined_all/region_DARs_redone/"
files <- list.files(result_dir, "knownResults.txt", full.names = F, recursive = T)
files = files[grep(cl, files)]

all_results <- list()

for (file in files) {
  celltype <- gsub("_motifs/knownResults.txt", "", (file))
  celltype <- gsub("_peaks", "", (celltype))
  df <- fread(file)
  
  if(nrow(df) < 1) next
  df$CellType <- celltype
  df$Direction = ifelse(grepl("up", file), "Up", "Down")
  colnames(df) = sapply(strsplit(as.character(colnames(df)), "[(]"), `[`, 1)
  
  all_results[[celltype]] <- df
}

combined_df <- do.call(rbind, all_results)
combined_df$Motif = sapply(strsplit(as.character(combined_df$`Motif Name`), "/"), `[`, 1)

top_upregulated <- combined_df %>%
  filter(`q-value ` < 0.05 & `P-value` < 0.05 & Direction == "Up") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 30) %>%
  ungroup()

top_downregulated <- combined_df %>%
  filter(`q-value ` < 0.05 & `P-value` < 0.05 & Direction == "Down") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 30) %>%
  ungroup()

top_pathways <- bind_rows(top_upregulated, top_downregulated)
unique_pathways <- unique(top_pathways$Motif)

plot_df <- combined_df %>%
  filter(Motif %in% unique_pathways) %>%
  distinct(Motif, CellType, .keep_all = TRUE)

motif_matrix <- plot_df %>%
  select(Motif, CellType, `Log P-value`) %>%
  spread(CellType, `Log P-value`, fill = 0)

motif_matrix <- as.data.frame(motif_matrix)
rownames(motif_matrix) <- motif_matrix$Motif
motif_matrix$Motif <- NULL

# Perform hierarchical clustering
motif_dist <- dist(motif_matrix)
motif_clust <- hclust(motif_dist)

# Cut the tree to form clusters of similar motifs (set height threshold)
clusters <- cutree(motif_clust, h = 10)  # Adjust `h` for similarity threshold

# Add cluster IDs to plot_df
plot_df$Cluster <- clusters[match(plot_df$Motif, rownames(motif_matrix))]

# For each cluster, keep only the most significant motif (lowest p-value)
plot_df <- plot_df %>%
  group_by(Cluster) %>%
  slice(which.min(`P-value`)) %>%
  ungroup()

# Reorder motifs based on clustering order
plot_df$Motif <- factor(plot_df$Motif, levels = rownames(motif_matrix)[motif_clust$order])

most_significant_motifs <- plot_df %>%
  group_by(Cluster) %>%
  filter(`P-value` == min(`P-value`)) %>%
  ungroup() %>%
  pull(Motif)

# Now filter plot_df to retain all rows corresponding to the most significant motifs
plot_df_filtered <- plot_df %>%
  filter(Motif %in% most_significant_motifs)

# Reorder motifs based on the clustering order
plot_df_filtered$Motif <- factor(plot_df_filtered$Motif, levels = rownames(motif_matrix)[motif_clust$order])

# Create the dot plot
plot_df_filtered$significance_direction <- plot_df_filtered$`Log P-value`
plot_df_filtered$significance_direction[plot_df_filtered$Direction == "Up"] <- plot_df_filtered$`Log P-value`[plot_df_filtered$Direction == "Up"] * -1
plot_df_filtered$significance_direction <- pmin(pmax(plot_df_filtered$significance_direction, -100), 100)

motifs <- ggplot(plot_df_filtered, aes(x = Motif, y = CellType, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-100, 100), name = "Significance x Direction") +
  scale_size_continuous(range = c(3, 10), name = "-log10(Adjusted p-value)") +
  theme_minimal() +
  xlab("Motif (Clustered)") +
  ylab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size = 8), legend.position = "left") +
  facet_wrap(~Direction, nrow = 2, scales = "free_y")

motifs


plot_df$Motif <- factor(plot_df$Motif, levels = rownames(motif_matrix)[motif_clust$order])

# Create the dot plot
#plot_df$Description <- substr(plot_df$Description, 1, 50)
plot_df$significance_direction <- plot_df$`Log P-value`
plot_df$significance_direction[plot_df$Direction == "Up"] <- plot_df$`Log P-value`[plot_df$Direction == "Up"] * -1
plot_df$significance_direction <- pmin(pmax(plot_df$significance_direction, -100), 100)

motifs <- ggplot(plot_df, aes(x = Motif, y = CellType, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-100, 100), name = "Significance x Direction") +
  scale_size_continuous(range = c(3, 10), name = "-log10(Adjusted p-value)") +
  theme_minimal() +
  xlab("Motif (Clustered)") +
  ylab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 65, hjust = 1), axis.text.y = element_text(size = 8), legend.position = "left") +
  facet_wrap(~Direction, nrow = 2, scales = "free_y")

motifs


