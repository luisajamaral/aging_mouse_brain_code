library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
setwd("~/projects/combined_all/Figures")




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

cl = "All"

me =ameta#[grep(cl,ameta$celltype_final),]
me$celltype_final = gsub("/", "-" , me$celltype_final)
me$celltype_region = gsub("/", "-" , me$celltype_region)
ord =  names(table(me$celltype_final))[order(table(me$celltype_final), decreasing = T)]


input_motifs  =c("AP-1(bZIP)", "Fra1(bZIP)","CTCF(Zf)","BORIS(Zf)","RFX(HTH)", "Olig2(bHLH)", "NeuroG2(bHLH)", "Mef2b(MADS)" , "Tgif1(Homeobox)", "Sox10(HMG)", "NRF(NRF)", 
          "Sp5(Zf)","IRF1(IRF)","ETS(ETS)","Klf4(Zf)","NFY(CCAAT)","Eomes(T-box)","PU.1:IRF8(ETS:IRF)", "Lhx1(Homeobox)", "YY1(Zf)","Tbx5(T-box)",
          "Snail1(Zf)","Sp1(Zf)","Foxo3(Forkhead)", "Twist2(bHLH)","CREB5(bZIP)","Oct6(POU,Homeobox)","KLF14(Zf)","Sp1(Zf)","Sox3(HMG)","E2A(bHLH),near_PU.1",
          "NeuroG2(bHLH)", "Rfx1(HTH)", "X-box(HTH)", "Slug(Zf)", "BHLHA15(bHLH)", "Pbx3(Homeobox)", "Mef2c(MADS)","Mef2a(MADS)","Lhx2(Homeobox)", "KLF5(Zf)","Nrf2(bZIP)", "Nrf3(bZIP)","NRF3(bZIP)")


setwd("~/projects/combined_all/combined_DARs_redone/")
# Set the directory containing the result files
result_dir <- "~/projects/combined_all/combined_DARs_redone/"
# Load all result files
files <- list.files(result_dir, "knownResults.txt",full.names = F, recursive = T)
files = files#[grep(cl,files)]
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
  slice_head(n = 8) %>%
  ungroup()

top_downregulated <- combined_df %>%
    filter(`q-value ` < 1e-10 & `P-value` < 1e-3 & Direction=="Down") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n =8) %>%
  ungroup()

# Combine the top upregulated and downregulated pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated)

unique_pathways <- unique(top_pathways$Motif)
unique_pathways = unique(c(unique_pathways, input_motifs))
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




plot_df$CellType = gsub("up_", "", plot_df$CellType)
plot_df$CellType = gsub("down_", "", plot_df$CellType)

plot_df = plot_df %>%
  group_by(CellType) %>%
  filter(any(`q-value ` < 0.1)) %>%  # Keep only CellTypes with q-value < 0.05
  ungroup()



plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)

plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))

plot_df = plot_df[which(plot_df$`Log P-value`< -3),]


plot_df$clade= "NN"
plot_df$clade[grep("Glut", plot_df$CellType)] = "Glut"
plot_df$clade[grep("Gaba", plot_df$CellType)] = "Gaba"


options(repr.plot.width=17, repr.plot.height=9
       )
plot_df$significance_direction[which(plot_df$significance_direction>25)] = 25
plot_df$significance_direction[which(plot_df$significance_direction< -25)] = -25
plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))

plot_df = plot_df[order(plot_df$`Log P-value`),]
motifs= ggplot(plot_df, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,  name = "Significance") +
  scale_size_continuous( name = "-log10(p-value)" ,range = c(1,6)) +
  theme_minimal() +
  ylab("Motif (Clustered)") +
  xlab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10), strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), legend.position = "bottom") + coord_flip() + facet_grid(clade~., scales= "free", space = "free")
motifs
#ggsave(motifs, height = 5.5, width = 8.5, file = "~/projects/combined_all/Figures/Figure4-Gaba/Gaba_Motifs.pdf")
#ggsave(motifs, height = 6.5, width = 10, file = "~/projects/combined_all/Figures/Figure3-NN/NN_Motifs.pdf")
#ggsave(motifs, height = 6.5, width = 10, file = "~/projects/combined_all/Figures/Figure3-NN/NN_Motifs.svg")


unique(plot_df$Motif)



options(repr.plot.width=9, repr.plot.height=9
       )

motifs  =c("AP-1(bZIP)", "Fra1(bZIP)","CTCF(Zf)","BORIS(Zf)","RFX(HTH)", "Olig2(bHLH)", "NeuroG2(bHLH)", "Mef2b(MADS)" , "Tgif1(Homeobox)", "Sox10(HMG)", "NRF(NRF)", 
          "Sp5(Zf)","IRF1(IRF)","ETS(ETS)","Klf4(Zf)","NFY(CCAAT)","Eomes(T-box)","PU.1:IRF8(ETS:IRF)", "Lhx1(Homeobox)", "YY1(Zf)","Tbx5(T-box)",
          "Snail1(Zf)","Sp1(Zf)","Foxo3(Forkhead)", "Twist2(bHLH)","CREB5(bZIP)","Oct6(POU,Homeobox)","KLF14(Zf)","Sp1(Zf)","Sox3(HMG)","E2A(bHLH),near_PU.1",
          "NeuroG2(bHLH)", "Rfx1(HTH)", "X-box(HTH)", "Slug(Zf)", "BHLHA15(bHLH)", "Pbx3(Homeobox)", "Lhx2(Homeobox)", "KLF5(Zf)","Nrf2(bZIP)", "Nrf3(bZIP)","NRF3(bZIP)")
sub = plot_df[which(plot_df$Motif %in% motifs ), ]

motifs= ggplot(sub, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance_direction)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,  name = "Significance") +
  scale_size_continuous( name = "-log10(p-value)" ,range = c(1,6)) +
  theme_minimal() +
  ylab("Motif (Clustered)") +
  xlab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10), strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), legend.position = "right") + coord_flip() + facet_grid(clade~., scales= "free", space = "free")
motifs
#ggsave(motifs, height = 5.5, width = 8.5, file = "~/projects/combined_all/Figures/Figure4-Gaba/Gaba_Motifs.pdf")
#ggsave(motifs, height = 6.5, width = 10, file = "~/projects/combined_all/Figures/Figure3-NN/NN_Motifs.pdf")
#ggsave(motifs, height = 6.5, width = 10, file = "~/projects/combined_all/Figures/Figure3-NN/NN_Motifs.svg")


ggsave(motifs, height = 9, width = 9, file = "~/projects/combined_all/Figures/All_Motifs_combined_small.pdf")



