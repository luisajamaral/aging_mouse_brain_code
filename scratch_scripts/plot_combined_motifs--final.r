library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(RColorBrewer)
setwd("~/projects/combined_all/Figures")






cts = read.table("celltypes_figure6.txt")
head(cts)

ord = cts$x

input_motifs  =c("AP-1(bZIP)", "Fra1(bZIP)","CTCF(Zf)","BORIS(Zf)","RFX(HTH)", "Olig2(bHLH)", "NeuroG2(bHLH)", "Mef2b(MADS)" , "Tgif1(Homeobox)", "Sox10(HMG)", "NRF(NRF)", 
          "Sp5(Zf)","IRF1(IRF)","ETS(ETS)","Klf4(Zf)","NFY(CCAAT)","Eomes(T-box)","PU.1:IRF8(ETS:IRF)", "Lhx1(Homeobox)", "YY1(Zf)","Tbx5(T-box)",
          "Snail1(Zf)","Sp1(Zf)","Foxo3(Forkhead)", "Twist2(bHLH)","CREB5(bZIP)","Oct6(POU,Homeobox)","KLF14(Zf)","Sp1(Zf)","Sox3(HMG)","E2A(bHLH),near_PU.1",
          "NeuroG2(bHLH)", "Rfx1(HTH)", "X-box(HTH)", "Slug(Zf)", "BHLHA15(bHLH)", "Pbx3(Homeobox)", "Mef2c(MADS)","Mef2a(MADS)","Lhx2(Homeobox)", "KLF5(Zf)","Nrf2(bZIP)", "Nrf3(bZIP)","NRF3(bZIP)")


ord = c(ord , "IOL_NN")

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

#setwd("~/projects/combined_all/combined_DARs_redone/motifs")
setwd("~/projects/combined_all/Figures/Figure6-Het-TEs/motifs_h3k9/")

input_motifs  =c("AP-1(bZIP)", "Fos(bZIP)","JunB(bZIP)","Fra1(bZIP)","CTCF(Zf)","Foxo3(Forkhead)","IRF1(IRF)","Klf4(Zf)")
#input_motifs = c("Mef2b(MADS)",  "NeuroG2(bHLH)","Rfx1(HTH)", "Lhx2(Homeobox)", "Dlx2(Homeobox)",  
# "Olig2(bHLH)", "Sox10(HMG)", "Tbr1(T-box)", "Tgif1(Homeobox)", "Mef2a(MADS)")

# Set the directory containing the result files
#result_dir <- "~/projects/combined_all/combined_DARs_redone/motifs"
result_dir <- "~/projects/combined_all/Figures/Figure6-Het-TEs/motifs_h3k9/"


files <- list.files(result_dir, "knownResults.txt",full.names = F, recursive = T)
files = files#[grep(cl,files)]
files = files[grep("up", files)]


# Initialize an empty list to store data frames
all_results <- list()

# Loop through each file and read the data
for (file in files) {
  # Extract the cell type from the file name
  celltype <- gsub("_motifs/knownResults.txt", "", (file))
  celltype <- gsub("_peaks", "", (celltype))
  celltype <- gsub("H3K9me3_", "", (celltype))
  celltype <- gsub(".bed", "", (celltype))


  # Read the data
  df <- fread(file)
  if(nrow(df)<1){
      next
  }
  # Add cell type as a column
  df$Direction = "Up"
  colnames(df)  = sapply(strsplit(as.character(colnames(df)), "[(]"), `[`, 1)
  df$CellType = celltype
  # Store the results
  all_results[[celltype]] <- df
}

# Combine all data frames into one
combined_df <- do.call(rbind, all_results)
combined_df$Motif = sapply(strsplit(as.character(combined_df$`Motif Name`), "/"), `[`, 1)
combined_df$`P-value`= as.numeric(combined_df$`P-value`)
combined_df$`q-value`= as.numeric(combined_df$`q-value`)



#input_motifs  =c("AP-1(bZIP)", "Fos(bZIP)","JunB(bZIP)","Fra1(bZIP)","CTCF(Zf)","Foxo3(Forkhead)","IRF1(IRF)","Klf4(Zf)")
input_motifs  =c("AP-1(bZIP)", "Fos(bZIP)","JunB(bZIP)","Fra1(bZIP)")


# Identify the top 10 upregulated and top 10 downregulated pathways by lowest p-value for any cell type
top_upregulated <- combined_df %>%
  filter(`q-value ` < 1e-15 & `P-value` < 1e-15 ) %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 3) %>%
  ungroup()


# Combine the top upregulated and downregulated pathways
top_pathways <- bind_rows(top_upregulated)#, top_downregulated)

unique_pathways <- unique(top_pathways$Motif)
unique_pathways = unique(c(unique_pathways, input_motifs))
plot_df <- combined_df %>%
  filter(Motif %in% unique_pathways)

plot_df <- plot_df %>%
  distinct(Motif, CellType, .keep_all = TRUE)

plot_df$significance = -plot_df$`Log P-value`
plot_df$significance[which(plot_df$significance>200)] = 200

head(plot_df)
# Perform hierarchical clustering on the motifs based on Log P-value
# First, create a matrix of Log P-values with motifs as rows and cell types as columns
motif_matrix <- plot_df %>%
  select(Motif, CellType, significance) %>%
  spread(CellType, significance, fill = 0)

# Convert the Motif column to row names (base R way)
motif_matrix <- as.data.frame(motif_matrix)
rownames(motif_matrix) <- motif_matrix$Motif
motif_matrix$Motif <- NULL  # Remove the Motif column after assigning it to rownames

# Compute the distance matrix and hierarchical clustering
motif_dist <- dist(motif_matrix)
motif_clust <- hclust(motif_dist, method = "ward.D2")

# Reorder motifs based on the clustering
plot_df$Motif <- factor(plot_df$Motif, levels = rownames(motif_matrix)[motif_clust$order])

# Create the dotplot with clustered motifs
options(repr.plot.width=8, repr.plot.height=5.5)

plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)
plot_df$ct = gsub("up_", "", plot_df$ct )
plot_df$ct = gsub("down_", "", plot_df$ct )
plot_df$ct = gsub(".bed", "", plot_df$ct )
head(plot_df)

plot_df$CellType = gsub("up_", "", plot_df$CellType)

plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)

plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))

#plot_df = plot_df[which(plot_df$`Log P-value`< -1),]
plot_df$clade= "NN"
plot_df$clade[grep("Glut", plot_df$CellType)] = "Glut"
plot_df$clade[grep("Gaba", plot_df$CellType)] = "Gaba"
plot_df$clade[grep("Neur", plot_df$CellType)] = "Gaba"

plot_df$clade = factor(plot_df$clade, levels = c("NN","Gaba", "Glut"))


plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))
plot_df = plot_df[which(plot_df$ct %in% cts$x),]
plot_df$CellType = factor(plot_df$CellType, levels = rev(unique(factor(plot_df$CellType))))


options(repr.plot.width=5, repr.plot.height=8
       )
#plot_df$significance[which(plot_df$significance>25)] = 25


plot_df$`Log P-value`[which(plot_df$`Log P-value` < -50)] = -50

plot_df = plot_df[order(plot_df$`Log P-value`),]

plot_df$target = as.numeric(gsub("%", "" , plot_df$`% of Target Sequences with Motif`))


plot_df = plot_df[which(plot_df$Motif%in%input_motifs),]
plot_df$Motif = factor(plot_df$Motif, levels = unique(input_motifs))


options(repr.plot.width=6.5, repr.plot.height=8.5)
motifs <- ggplot(plot_df, aes(x = CellType, y = Motif, color = -`Log P-value`, size = target)) +
  geom_point() +
  scale_color_distiller(palette = "Reds", direction = 1, name = "-log10(p-value)") +  # Apply OrRd from RColorBrewer
  scale_size_continuous(name = "% of Target cCREs", range = c(2,6)) +
  theme_bw() +
  ylab("Motif") +
  xlab("Cell Type") +theme_bw()+
  ggtitle("Motifs Enriched in Age-Up cCREs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), text = element_text(size =15), 
        #strip.text.y = element_blank(),
        axis.text.y = element_text(size = 14), 
        legend.position = "right",panel.border = element_rect(colour="grey")) + 
  coord_flip() + 
 facet_grid(clade ~ ., scales = "free", space = "free")

motifs
pdf("../../../Figures/up_motifs_H3K9me3.pdf", height =8.5, width =6.5)
print(motifs)
dev.off()



setwd("~/projects/combined_all/combined_DARs_redone/motifs")
#setwd("~/projects/combined_all/Figures/Figure6-Het-TEs/motifs_h3k9/")

#input_motifs  =c("Mef2b", "Fra1(bZIP)","CTCF(Zf)","BORIS(Zf)","YY1(Zf)")
input_motifs = c("Mef2b(MADS)",  "NeuroG2(bHLH)","Rfx1(HTH)", "Lhx2(Homeobox)", "Dlx2(Homeobox)",  
 "Olig2(bHLH)", "Sox10(HMG)", "Tbr1(T-box)", "Tgif1(Homeobox)", "Mef2a(MADS)")

# Set the directory containing the result files
result_dir <- "~/projects/combined_all/combined_DARs_redone/motifs"
#result_dir <- "~/projects/combined_all/Figures/Figure6-Het-TEs/motifs_h3k9/"


files <- list.files(result_dir, "knownResults.txt",full.names = F, recursive = T)
files = files#[grep(cl,files)]
files = files[grep("down", files)]


# Initialize an empty list to store data frames
all_results <- list()

# Loop through each file and read the data
for (file in files) {
  # Extract the cell type from the file name
  celltype <- gsub("_motifs/knownResults.txt", "", (file))
  celltype <- gsub("_peaks", "", (celltype))
  celltype <- gsub("H3K9me3_", "", (celltype))
  celltype <- gsub(".bed", "", (celltype))


  # Read the data
  df <- fread(file)
  if(nrow(df)<1){
      next
  }
  # Add cell type as a column
  df$Direction = "Down"
  colnames(df)  = sapply(strsplit(as.character(colnames(df)), "[(]"), `[`, 1)
  df$CellType = celltype
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
  filter(`q-value ` < 1e-25 & `P-value` < 1e-25 ) %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 1) %>%
  ungroup()


# Combine the top upregulated and downregulated pathways
top_pathways <- bind_rows(top_upregulated)#, top_downregulated)

unique_pathways <- unique(top_pathways$Motif)
unique_pathways = unique(c(unique_pathways, input_motifs))
plot_df <- combined_df %>%
  filter(Motif %in% unique_pathways)

plot_df <- plot_df %>%
  distinct(Motif, CellType, .keep_all = TRUE)

plot_df$significance = -plot_df$`Log P-value`
plot_df$significance[which(plot_df$significance>200)] = 200

head(plot_df)
# Perform hierarchical clustering on the motifs based on Log P-value
# First, create a matrix of Log P-values with motifs as rows and cell types as columns
motif_matrix <- plot_df %>%
  select(Motif, CellType, significance) %>%
  spread(CellType, significance, fill = 0)

# Convert the Motif column to row names (base R way)
motif_matrix <- as.data.frame(motif_matrix)
rownames(motif_matrix) <- motif_matrix$Motif
motif_matrix$Motif <- NULL  # Remove the Motif column after assigning it to rownames

# Compute the distance matrix and hierarchical clustering
motif_dist <- dist(motif_matrix)
motif_clust <- hclust(motif_dist, method = "ward.D2")

# Reorder motifs based on the clustering
plot_df$Motif <- factor(plot_df$Motif, levels = rownames(motif_matrix)[motif_clust$order])

# Create the dotplot with clustered motifs
options(repr.plot.width=8, repr.plot.height=5.5)

plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)
plot_df$ct = gsub("up_", "", plot_df$ct )
plot_df$ct = gsub("down_", "", plot_df$ct )
plot_df$ct = gsub(".bed", "", plot_df$ct )
head(plot_df)

plot_df$CellType = gsub("down_", "", plot_df$CellType)

plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)

plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))

#plot_df = plot_df[which(plot_df$`Log P-value`< -1),]
plot_df$clade= "NN"
plot_df$clade[grep("Glut", plot_df$CellType)] = "Glut"
plot_df$clade[grep("Gaba", plot_df$CellType)] = "Gaba"
plot_df$clade[grep("Neur", plot_df$CellType)] = "Gaba"

plot_df$clade = factor(plot_df$clade, levels = c("NN","Gaba", "Glut"))


plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))
plot_df = plot_df[which(plot_df$ct %in% cts$x),]
plot_df$CellType = factor(plot_df$CellType, levels = rev(unique(factor(plot_df$CellType))))


options(repr.plot.width=5, repr.plot.height=8
       )
#plot_df$significance[which(plot_df$significance>25)] = 25


plot_df$`Log P-value`[which(plot_df$`Log P-value` < -50)] = -50

plot_df = plot_df[order(plot_df$`Log P-value`),]

plot_df$target = as.numeric(gsub("%", "" , plot_df$`% of Target Sequences with Motif`))


plot_df = plot_df[which(plot_df$Motif%in%input_motifs),]
#plot_df$Motif = factor(plot_df$Motif, levels = unique(input_motifs))


options(repr.plot.width=8, repr.plot.height=8.5)
motifs <- ggplot(plot_df, aes(x = CellType, y = Motif, color = -`Log P-value`, size = target)) +
  geom_point() +
  scale_color_distiller(palette = "Blues", direction = 1, name = "-log10(p-value)") +  # Apply OrRd from RColorBrewer
  scale_size_continuous(name = "% of Target cCREs", range = c(1,6)) +
  theme_bw() +
  ylab("Motif") +
  xlab("Cell Type") +theme_bw()+
  ggtitle("Motifs Enriched in Age-Down cCREs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), text = element_text(size =15), 
        #strip.text.y = element_blank(),
        axis.text.y = element_text(size = 14), 
        legend.position = "right",panel.border = element_rect(colour="grey")) + 
  coord_flip() + 
 facet_grid(clade ~ ., scales = "free", space = "free")

motifs

#pdf("../../Figures/down_motifs_2.pdf", height =8.5, width =8)
print(motifs)
#dev.off()



options(repr.plot.width=8, repr.plot.height=8.5)
motifs <- ggplot(plot_df, aes(x = CellType, y = Motif, color = -`Log P-value`, size = target)) +
  geom_point() +
  scale_color_distiller(palette = "Blues", direction = 1, name = "-log10(p-value)") +  # Apply OrRd from RColorBrewer
  scale_size_continuous(name = "% of Target cCREs", range = c(1,6)) +
  theme_bw() +
  ylab("Motif") +
  xlab("Cell Type") +theme_bw()+
  ggtitle("Motifs Enriched in Age-Down cCREs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16), text = element_text(size =15), 
        #strip.text.y = element_blank(),
        axis.text.y = element_text(size = 14), 
        legend.position = "right",panel.border = element_rect(colour="grey")) + 
  coord_flip() + 
 facet_grid(clade ~ ., scales = "free", space = "free")

motifs

pdf("../../Figures/down_motifs_2.pdf", height =8.5, width =8)
print(motifs)
dev.off()



plot_df = plot_df[which(plot_df$Motif%in%input_motifs),]




options(repr.plot.width=4, repr.plot.height=8.5)

motifs <- ggplot(plot_df, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance)) +
  geom_point() +
  scale_color_distiller(palette = "Blues", direction = 1, name = "Significance") +  # Apply OrRd from RColorBrewer
  scale_size_continuous(name = "-log10(p-value)", range = c(1,6)) +
  theme_bw() +
  ylab("Motif (Clustered)") +
  xlab("Cell Type") +
  ggtitle("Motifs Enriched in Age-Down DARs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), 
        legend.position = "top",panel.border = element_rect(colour="grey")) + 
  coord_flip() + 
 facet_grid(clade ~ ., scales = "free", space = "free")

motifs


pdf("../../Figures/down_motifs_toplegend.pdf", height =8.5, width =4)
print(motifs)
dev.off()

pdf("../../Figures/down_motifs.pdf", height =8.5, width =6)
print(motifs)
dev.off()

#setwd("~/projects/combined_all/Figures/Figure6-Het-TEs/motifs_h3k9/")
setwd("~/projects/combined_all/combined_DARs_redone/motifs")
input_motifs  =c("AP-1(bZIP)", "Fra1(bZIP)","CTCF(Zf)","BORIS(Zf)")#,"YY1(Zf)","Klf4(Zf)","Foxo3(Forkhead)","IRF1(IRF)")
                 
# Set the directory containing the result files
result_dir <- "~/projects/combined_all/combined_DARs_redone/motifs"
#result_dir <- "~/projects/combined_all/Figures/Figure6-Het-TEs/motifs_h3k9/"


files <- list.files(result_dir, "knownResults.txt",full.names = F, recursive = T)
files = files#[grep(cl,files)]
files = files[grep("up", files)]


# Initialize an empty list to store data frames
all_results <- list()

# Loop through each file and read the data
for (file in files) {
  # Extract the cell type from the file name
  celltype <- gsub("_motifs/knownResults.txt", "", (file))
  celltype <- gsub("_peaks", "", (celltype))
  celltype <- gsub("H3K9me3_", "", (celltype))
  celltype <- gsub(".bed", "", (celltype))

  # Read the data
  df <- fread(file)
  if(nrow(df)<1){
      next
  }
  # Add cell type as a column
  df$Direction = "Up"
  colnames(df)  = sapply(strsplit(as.character(colnames(df)), "[(]"), `[`, 1)
  df$CellType = celltype
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
  filter(`q-value ` < 1e-15 & `P-value` < 1e-2  & Direction=="Up") %>%
  group_by(CellType) %>%
  arrange(`P-value`) %>%
  slice_head(n = 1) %>%
  ungroup()


# Combine the top upregulated and downregulated pathways
top_pathways <- bind_rows(top_upregulated)#, top_downregulated)

unique_pathways <- unique(top_pathways$Motif)
unique_pathways = unique(c(unique_pathways, input_motifs))
plot_df <- combined_df %>%
  filter(Motif %in% unique_pathways)

plot_df <- plot_df %>%
  distinct(Motif, CellType, .keep_all = TRUE)

plot_df$significance = -plot_df$`Log P-value`

plot_df$significance[which(plot_df$significance>50)] = 50

head(plot_df)
# Perform hierarchical clustering on the motifs based on Log P-value
# First, create a matrix of Log P-values with motifs as rows and cell types as columns
motif_matrix <- plot_df %>%
  select(Motif, CellType, significance) %>%
  spread(CellType, significance, fill = 0)

# Convert the Motif column to row names (base R way)
motif_matrix <- as.data.frame(motif_matrix)
rownames(motif_matrix) <- motif_matrix$Motif
motif_matrix$Motif <- NULL  # Remove the Motif column after assigning it to rownames

# Compute the distance matrix and hierarchical clustering
motif_dist <- dist(motif_matrix)
motif_clust <- hclust(motif_dist, method = "ward.D2")

# Reorder motifs based on the clustering
plot_df$Motif <- factor(plot_df$Motif, levels = rownames(motif_matrix)[motif_clust$order])

# Create the dotplot with clustered motifs
options(repr.plot.width=8, repr.plot.height=5.5)

plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)
plot_df$ct = gsub("up_", "", plot_df$ct )
plot_df$ct = gsub("down_", "", plot_df$ct )
plot_df$ct = gsub(".bed", "", plot_df$ct )
head(plot_df)

plot_df$CellType = gsub("up_", "", plot_df$CellType)

plot_df$ct = sapply(strsplit(as.character(plot_df$CellType), ":"), `[`, 1)

plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))

#plot_df = plot_df[which(plot_df$`Log P-value`< -1),]
plot_df$clade= "NN"
plot_df$clade[grep("Glut", plot_df$CellType)] = "Glut"
plot_df$clade[grep("Gaba", plot_df$CellType)] = "Gaba"
plot_df$clade[grep("Neur", plot_df$CellType)] = "Gaba"

plot_df$clade = factor(plot_df$clade, levels = c("NN","Gaba", "Glut"))


plot_df$ct = factor(plot_df$ct, levels = gsub(" ", "_", ord))
plot_df = plot_df[which(plot_df$ct %in% cts$x),]
plot_df$CellType = factor(plot_df$CellType, levels = rev(unique(factor(plot_df$CellType))))


options(repr.plot.width=5, repr.plot.height=8
       )
#plot_df$significance[which(plot_df$significance>25)] = 25


plot_df$`Log P-value`[which(plot_df$`Log P-value` < -50)] = -50

plot_df = plot_df[order(plot_df$`Log P-value`),]

#plot_df = plot_df[which(plot_df$Motif%in%input_motifs),]
#plot_df$Motif = factor(plot_df$Motif, levels = unique(input_motifs))


options(repr.plot.width=5, repr.plot.height=8.5)
motifs <- ggplot(plot_df, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance)) +
  geom_point() +
  scale_color_distiller(palette = "Reds", direction = 1, name = "Significance") +  # Apply OrRd from RColorBrewer
  scale_size_continuous(name = "-log10(p-value)", range = c(1,6)) +
  theme_bw() +
  ylab("Motif") +
  xlab("Cell Type") +
  ggtitle("Motifs Enriched in Age-Up DARs in H3K9me3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), 
        legend.position = "right",panel.border = element_rect(colour="grey")) + 
  coord_flip() + 
 facet_grid(clade ~ ., scales = "free", space = "free")

motifs

#pdf("../../../Figures/up_motifs_in_h3k9me3.pdf", height =8.5, width =5)
#print(motifs)
#dev.off()


top_upregulated[order(top_upregulated$`P-value`),]

options(repr.plot.width=4.5, repr.plot.height=8.5)
motifs <- ggplot(plot_df, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance)) +
  geom_point() +
  scale_color_distiller(palette = "Reds", direction = 1, name = "Significance") +  # Apply OrRd from RColorBrewer
  scale_size_continuous(name = "-log10(p-value)", range = c(1,6)) +
  theme_bw() +
  ylab("Motif (Clustered)") +
  xlab("Cell Type") +
  ggtitle("Motifs Enriched in Age-Up DARs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
        strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), 
        legend.position = "top",panel.border = element_rect(colour="grey")) + 
  coord_flip() + 
 facet_grid(clade ~ ., scales = "free", space = "free")

motifs

pdf("../../Figures/up_motifs_toplegend1.pdf", height =8.5, width =3.5)
print(motifs)
dev.off()

# Ensure `% of Target Sequences with Motif` is numeric
plot_df$TargetPct <- as.numeric(gsub("%", "", plot_df$`% of Target Sequences with Motif`)) / 100

# Cap extreme values to avoid overly large/small dots
plot_df$TargetPct[plot_df$TargetPct > 0.5] <- 0.5  # Optional: Cap at 50% if necessary

options(repr.plot.width=5, repr.plot.height=8.5)

motifs <- ggplot(plot_df, aes(x = CellType, y = Motif, size = TargetPct, color = significance)) +
  geom_point() +
  scale_color_distiller(palette = "Reds", direction = 1, name = "Significance") +  
  scale_size_continuous(name = "% Target Sequences", breaks = seq(0, 0.5, 0.1), range = c(1, 6)) +
  theme_bw() +
  ylab("Motif (Clustered)") +
  xlab("Cell Type") +
  ggtitle("Motifs Enriched in Age-Up DARs") +  
  theme(
    legend.position = "right",  # Place legend at the top
    legend.box = "vertical",  # Stack size and fill legends vertically
    legend.spacing.y = unit(0.3, "cm"),  # Add space between legend elements
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    strip.text.y = element_blank(),
    axis.text.y = element_text(size = 10),
    panel.background = element_blank()
  ) + coord_flip() +
  facet_grid(clade ~ ., scales = "free", space = "free")

motifs

pdf("../../Figures/up_motifs.pdf", height = 8.5, width = 4)
print(motifs)
dev.off()


head(df)

options(repr.plot.width=5, repr.plot.height=9.5)

motifs <- ggplot(plot_df, aes(x = CellType, y = Motif, size = -`Log P-value`, color = significance)) +
  geom_point() +
  scale_color_distiller(palette = "Reds", direction = 1, name = "Significance") +  # Apply OrRd from RColorBrewer
  scale_size_continuous(name = "-log10(p-value)", range = c(1,6)) +
  theme_bw() +
  ylab("Motif") +
  xlab("Cell Type") +
  ggtitle("Motifs in Age-Up DARs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), 
        strip.text.y = element_blank(),
        axis.text.y = element_text(size = 14), text = element_text(size = 14),
        legend.position = "bottom",panel.border = element_rect(colour="grey")) + 
  coord_flip() + 
 facet_grid(clade ~ ., scales = "free", space = "free")
motifs

pdf("../../Figures/up_motifs.pdf", height =8.5, width =6.5)
print(motifs)
dev.off()

motifs_heatmap <- ggplot(plot_df, aes(x = CellType, y = Motif, fill = significance_direction)) +
  geom_tile() +  # Use tiles instead of points
  scale_fill_gradient2(low = "royalblue", mid = "white", high = "tomato2", midpoint = 0, name = "-log10(p-value)") +
  theme_bw() +
  ylab("Motif (Clustered)") +
  xlab("Cell Type") +
  ggtitle(paste("Age - Motifs", cl)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10), 
        strip.text.y = element_blank(),
        axis.text.y = element_text(size = 10), 
        legend.position = "right") +  coord_flip() +
  facet_grid(clade ~ ., scales = "free", space = "free")

motifs_heatmap


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
  scale_size_continuous( name = "-log10(p-value)" ,range = c(2,6)) +
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



