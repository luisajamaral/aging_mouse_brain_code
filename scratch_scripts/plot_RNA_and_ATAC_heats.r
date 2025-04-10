library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(reshape2)


meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")


meta = read.csv("~/ps-renlab2/projects/combined_all/female_RNA/meta_final.csv")


head(meta)

meta$celltype_final = gsub("Lymphoid", "Lymphoid NN", meta$celltype_final)



lkey=read.table("~/projects/combined_all/scTE/out.saf")
#lkey = lkey[-which(lkey$V5=="Retroposon"),]
ukey = lkey[-which(duplicated(lkey$V1)),]
rownames(ukey) = ukey$V1

bkey  = ukey[grep("_", rownames(ukey)), ]
rownames(bkey) = gsub("_", "-", rownames(bkey))
ukey = rbind(ukey,bkey)
ukey$V6 = gsub("[?]", "", ukey$V6)



files


setwd("~/projects/combined_all/female_RNA/SoloTE/subfamily_DESeq2/")

files = list.files(".", "only.txt", full.names = F)
cur = files[1]
tab = read.table(cur)
tab$gene = rownames(tab)
f = gsub("_DE_TEs_only.txt", "", files[1])

tab$celltype = f
#head(tab)
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_DE_TEs_only.txt", "", files[i])

    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}

tab$logFC = -tab$logFC


tab$gene = gsub("SoloTE-", "", tab$gene)
#tab$TE_type = ukey[paste(tab$gene), 6]


tab$type = "Gene"
tab$type[grep("SoloTE", rownames(tab))] = "SoloTE"
tab$direction = "Up in aging"
tab$direction[which(tab$logFC<0)] = "Down in aging"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("IMN", tab$celltype)]="NN"

#tab= tab[-grep("IT", tab$celltype),]

    
    
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","_", ord)


options(repr.plot.width=17, repr.plot.height=4)


d = tab
d = d[which(d$PValue<0.05),]
d = d[which(d$type!="Gene"),]
d$celltype = factor(d$celltype, levels = rev(ord))


g1 = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = position_stack(), color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "P < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))
print(g1)
    
#d = tab
##d = d[which(d$PValue<0.05),]
#d = d[which(d$type=="Gene"),]
#d$celltype = factor(d$celltype, levels = rev(ord))


#g = ggplot(d, aes(x = celltype, ,fill = direction)) +
#  geom_bar(stat = "count",position = "dodge", color = "black") +
#  theme_minimal() +
#  labs(title = "",
#       x = " Type",
#       y = "Number of DE") +
#  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "P < 0.05 genes"))+
#  theme(text = element_text(size = 18))
    
    #print(g)




options(repr.plot.width=18, repr.plot.height=8)

ggplot(d, aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", position = position_dodge(), color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "P < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))

head(tab[which(is.na(tab$TE_type)),])
head(tab)

options(repr.plot.width=7, repr.plot.height=5)

ggplot(tab, aes(x = logFC, color = clade)) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log2(fold_change) by type",
       x = "log2(fold_change)",
       y = "eCDF",
       color = "Type",
       linetype = "Major") +xlim(c(-2,2))+
  theme_minimal() 



head(tab)


rna = tab

setwd("~/projects/combined_all/scTE/celltype_DEseq2_results/")

files = list.files(".", "txt", full.names = F)
f = files[1]
tab = read.table(f)
tab$gene = rownames(tab)
f = gsub("_2vs18.txt", "", files[1])

tab$celltype = f
#head(tab)
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_2vs18.txt", "", files[i])

    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}

head(tab)

tab$logFC = -tab$logFC
tab$type = ukey[paste(tab$gene),"V6"]

tab$celltype = gsub("L2_3", "L2-3", tab$celltype)
tab$celltype = gsub("L2_3", "L2-3", tab$celltype)
tab$celltype = gsub("L5_6", "L5-6", tab$celltype)
tab$celltype = gsub("L6b_CT", "L6b-CT", tab$celltype)
tab$celltype = gsub("Lymphoid", "Lymphoid_NN", tab$celltype)


tab = tab[which(tab$celltype!="doublet"),]
unique(tab$celltype)[which(!unique(tab$celltype) %in% ord)]


tab$type = "Gene"
tab$type[grep("SoloTE", rownames(tab))] = "SoloTE"
#tab$direction = "Up in aging"
#tab$direction[which(tab$logFC>0)] = "Down in aging"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("IMN", tab$celltype)]="NN"
#tab$clade[grep("Inh-IMN", tab$celltype)]="Gaba"

#tab= tab[-grep("IT", tab$celltype),]

    
    
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","_", ord)


options(repr.plot.width=17, repr.plot.height=4)


d = tab
d = d[which(d$fdr<0.05),]
d$celltype = factor(d$celltype, levels = rev(ord))

options(repr.plot.width=17, repr.plot.height=7)





head(d)

options(repr.plot.width=18, repr.plot.height=8)

ggplot(d, aes(x = celltype, fill = dir)) +
  geom_bar(stat = "count", position = position_dodge(), color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "fdr < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))

head(tab)

head(rna)

sig= tab[which(tab$fdr<1 & abs(tab$logFC) > 0 ),]
sig_rna= rna[which(rna$fdr<1 & abs(rna$logFC) > 0),]


sig_genes = (unique(sig$gene))
length(sig_genes)

sig_genes_rna = (unique(sig_rna$gene))
length(sig_genes_rna)

sig_genes = sig_genes_rna[which(sig_genes_rna %in% sig_genes)]
length(sig_genes)

tab$type = "ATAC"
rna$type = "RNA"


rna = rna[,-c(11)]

head(rna)
head(tab)


rna$gene = ukey[paste(rna$gene ),1]
head(rna)

ukey[grep("MuRR" , rownames(ukey)),]

all = rbind(tab, rna)
#all$gene = gsub("-", "_", all$gene)

unique(tab$gene)[which(!unique(tab$gene) %in% unique(rna$gene))]

unique(rna$gene)[which(!unique(rna$gene) %in% unique(tab$gene))]

all$celltype_modality = paste(all$celltype, all$type)

head(all)

all = all[which(!is.na(all$gene)),]

all[which(all$gene =="L1MA5A"& all$celltype=="CA3_Glut"),]

library(dplyr)
atac_data <- all %>% filter(type == "ATAC")
rna_data <- all %>% filter(type == "RNA")

# Merge the dataframes on the gene and celltype columns
merged_data <- merge(atac_data, rna_data, by = c("gene", "celltype"), suffixes = c("_ATAC", "_RNA"))


merged_data$TE_type = ukey[paste(merged_data$gene),6]


merged_data = merged_data[which(merged_data$fdr_ATAC<0.01 | merged_data$PValue_RNA<0.01),]


nrow(merged_data)

head(merged_data)

cor(merged_data$logFC_RNA,merged_data$logFC_ATAC)

options(repr.plot.width=12, repr.plot.height=5)

ggplot(merged_data, aes(x = logFC_ATAC, y = logFC_RNA, color =significance_RNA )) +
  geom_point() + xlim(c(-1,1
                       ))+ylim(c(-3,3))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +  
    geom_text_repel(max.overlaps = 100, aes(label = ifelse((logFC_ATAC > 0.25 & logFC_RNA >1.25) |  
                                                           (logFC_ATAC < -0.25 & logFC_RNA < -1.25) & 
        (fdr_ATAC<0.05 &fdr_RNA<0.05) , gene, "")),
            hjust = 1.1, vjust = 1.1, check_overlap = TRUE) +

  labs(x = "logFC (ATAC)", y = "logFC (RNA)", title = "Scatter plot of logFC for ATAC and RNA (TEs)") +
  theme_minimal()+facet_wrap(~clade_RNA)

options(repr.plot.width=12, repr.plot.height=5)

ggplot(merged_data, aes(x = logFC_ATAC, y = logFC_RNA, color =significance_RNA )) +
  geom_point() + xlim(c(-1,1
                       ))+ylim(c(-3,3))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +  
    geom_text_repel(max.overlaps = 100, aes(label = ifelse((logFC_ATAC > 0.25 & logFC_RNA >1.25) |  
                                                           (logFC_ATAC < -0.25 & logFC_RNA < -1.25) & 
        (fdr_ATAC<0.05 &fdr_RNA<0.05) , gene, "")),
            hjust = 1.1, vjust = 1.1, check_overlap = TRUE) +

  labs(x = "logFC (ATAC)", y = "logFC (RNA)", title = "Scatter plot of logFC for ATAC and RNA (TEs)") +
  theme_minimal()+facet_wrap(~clade_RNA)

# Determine the threshold for labeling
label_threshold_ATAC <- .25
label_threshold_RNA <- 1.25

# Calculate fractions of points in each quadrant for each facet
calculate_fractions <- function(df) {
  n_total <- nrow(df)
  n_q1 <- nrow(df %>% filter(logFC_ATAC > 0 & logFC_RNA > 0))
  n_q2 <- nrow(df %>% filter(logFC_ATAC < 0 & logFC_RNA > 0))
  n_q3 <- nrow(df %>% filter(logFC_ATAC < 0 & logFC_RNA < 0))
  n_q4 <- nrow(df %>% filter(logFC_ATAC > 0 & logFC_RNA < 0))
  
  data.frame(
    clade_RNA = df$clade_RNA[1],
    frac_q1 = n_q1 / n_total,
    frac_q2 = n_q2 / n_total,
    frac_q3 = n_q3 / n_total,
    frac_q4 = n_q4 / n_total
  )
}

fractions <- merged_data %>%
  group_by(clade_RNA) %>%
  do(calculate_fractions(.))

# Create the scatter plot
p = ggplot(merged_data, aes(x = logFC_ATAC, y = logFC_RNA, color = TE_type)) +
  geom_point() +
  xlim(c(-1, 1)) + ylim(c(-3, 3)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  #geom_text_repel(max.overlaps = 100, aes(label = ifelse(gene %in% c("L1MA5A"), gene, "")),
  #                hjust = 1.1, vjust = 1.1, check_overlap = TRUE) +
  labs(x = "logFC (ATAC)", y = "logFC (RNA)", title = "Scatter plot of logFC for ATAC and RNA (TEs)") +
  theme_minimal() +
  facet_wrap(~clade_RNA) +
  geom_text(data = fractions, aes(x = 0.75, y = 2.95, label = paste("Q1: ", round(frac_q1, 2))), inherit.aes = FALSE, color = "black") +
  geom_text(data = fractions, aes(x = -0.75, y = 2.95, label = paste("Q2: ", round(frac_q2, 2))), inherit.aes = FALSE, color = "black") +
  geom_text(data = fractions, aes(x = -0.75, y = -2.95, label = paste("Q3: ", round(frac_q3, 2))), inherit.aes = FALSE, color = "black") +
  geom_text(data = fractions, aes(x = 0.75, y = -2.95, label = paste("Q4: ", round(frac_q4, 2))), inherit.aes = FALSE, color = "black")

p

ggsave(p,file="~/projects/combined_all/female_RNA/SoloTE/subfamily_DESeq2/logfc_rna_atac_p_0.05_both.pdf", width = 10, height = 4)

options(repr.plot.width=15, repr.plot.height=24)

ggplot(merged_data, aes(x = logFC_ATAC, y = logFC_RNA, color =TE_type )) +
  geom_point() + xlim(c(-1,1))+ylim(c(-3.2,3.2))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +  
    geom_text_repel(max.overlaps = 100, aes(label = ifelse(logFC_ATAC > 0.25
                                 & logFC_RNA >1.25 & (fdr_ATAC<0.05 |fdr_RNA<0.05) , gene, "")),
            hjust = 1.1, vjust = 1.1, check_overlap = TRUE) +

  labs(x = "logFC (ATAC)", y = "logFC (RNA)", title = "Scatter plot of logFC for ATAC and RNA") +
  theme_minimal()+facet_wrap(~celltype)

options(repr.plot.width=15, repr.plot.height=24)

ggplot(merged_data, aes(x = logFC_ATAC, y = logFC_RNA, color =clade )) +
  geom_point() + xlim(c(-1,1))+ylim(c(-3.2,3.2))+
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +  
    geom_text_repel(max.overlaps = 100, aes(label = ifelse(logFC_ATAC > 0.25
                                 & logFC_RNA >1.25 & (fdr_ATAC<0.05 |fdr_RNA<0.05) , gene, "")),
            hjust = 1.1, vjust = 1.1, check_overlap = TRUE) +

  labs(x = "logFC (ATAC)", y = "logFC (RNA)", title = "Scatter plot of logFC for ATAC and RNA") +
  theme_minimal()



cts = sapply(strsplit(as.character(colnames(heatmap_data)), " "), `[`, 1)
type = sapply(strsplit(as.character(colnames(heatmap_data)), " "), `[`, 2)


bad1 = cts[which(type == "ATAC")][which(!cts[which(type == "ATAC")]%in%cts[which(type == "RNA")])]
bad2 = cts[which(type == "RNA")][which(!cts[which(type == "RNA")]%in%cts[which(type == "ATAC")])]
bads = c(bad1,bad2)


heatmap_data = heatmap_data[,which(cts%in%bads)]
head(heatmap_data)

length(unique(all$gene)) 

#sig= tab[which(tab$fdr<.01 & abs(tab$logFC) > .3 ),]
#sig_rna= rna[which(tab$fdr<.01 & abs(tab$logFC) > .3 ),]

#sig_genes = (unique(sig$gene))
#sig_genes_rna = (unique(sig_rna$gene))

#sig_genes = sig_genes_rna[which(sig_genes_rna %in% sig_genes)]
all$modality = all$type

filtered_df <- all %>%
  filter(PValue < 0.01 & abs(logFC)> 0.05 )

# Get the top 5 genes for each celltype_modality by fdr
top_genes <- filtered_df %>%
  group_by(celltype_modality, ) %>%
  arrange(fdr) %>%
  slice_head(n = 100) %>%
  ungroup()

# Print the result
nrow(top_genes)


sig_genes =unique(top_genes$gene)
sig_genes
sig_genes = c(sig_genes , "L1MA5A")
options(repr.plot.width=15, repr.plot.height=12)
# Reshape the data for the heatmap
sig_tab = all[which(all$gene %in% c(sig_genes)),]
df = sig_tab
df = df[which(df$PValue<0.01),]
#df$gene = gsub("-", "_", df$gene)
heatmap_data <- dcast(df, gene ~ celltype_modality, value.var = "logFC")
row.names(heatmap_data) <- heatmap_data$gene
heatmap_data <- heatmap_data[,-1]

df$type = ukey[paste(df$gene), "V6"]

# Annotations
gene_annotation <- df[, c("gene", "type")]
gene_annotation <- unique(gene_annotation)
row.names(gene_annotation) <- gene_annotation$gene
gene_annotation <- gene_annotation[,-1, drop = FALSE]

celltype_annotation <- df[, c("celltype_modality", "clade")]
celltype_annotation <- unique(celltype_annotation)
row.names(celltype_annotation) <- celltype_annotation$celltype
celltype_annotation <- celltype_annotation[,-1, drop = FALSE]

# Colors for annotations
ann_colors <- list(
  type = brewer.pal(8, "Set2")[1:length(unique(gene_annotation$type))],
  clade = brewer.pal(8, "Set1")[1:length(unique(celltype_annotation$clade))]
)
names(ann_colors$type) <- unique(gene_annotation$type)
names(ann_colors$clade) <- c("Glut", "NN","Gaba")



# Define the color palette for the heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
heatmap_data[which(heatmap_data>1, arr.ind = T)] = 1
heatmap_data[which(heatmap_data< -1, arr.ind = T)] = -1

# Define breaks for the heatmap so that 0 is the midpoint
max_abs_val <- max(abs(heatmap_data), na.rm = TRUE)
breaks <- seq(-max_abs_val, max_abs_val, length.out = 101)

heatmap_data = heatmap_data[,c(grep("Glut",colnames(heatmap_data)), 
                               grep("Gaba",colnames(heatmap_data)),
                                    grep("IMN",colnames(heatmap_data)),
                                         grep("NN",colnames(heatmap_data)))]

heatmap_data = heatmap_data[,c(grep("ATAC",colnames(heatmap_data)), grep("RNA",colnames(heatmap_data)))]


head(heatmap_data)

ann_colors$clade = ann_colors$clade[c(1,3,2)]

L1MA5A = heatmap_data[grep("L1MA5A" , rownames(heatmap_data)), ] 



heatmap_data <- heatmap_data[rowMeans(is.na(heatmap_data)) <= .9, ]

heatmap_data = heatmap_data[which(apply(abs(heatmap_data), 1, mean, na.rm = TRUE)>.1),]
heatmap_data = heatmap_data[which(apply(abs(heatmap_data), 1, sum, na.rm = TRUE)>1),]
heatmap_data = heatmap_data[which(apply(abs(heatmap_data), 1, max, na.rm = TRUE)>.1),]
heatmap_data = rbind(heatmap_data, L1MA5A) 

heatmap_data = heatmap_data[,which(apply(heatmap_data, 2, max, na.rm = TRUE)>0.1)]
heatmap_data = heatmap_data[,which(apply(heatmap_data, 2, sum, na.rm = TRUE)>.1)]
heatmap_data <- heatmap_data[,colMeans(is.na(heatmap_data)) <= .9]


cts = sapply(strsplit(as.character(colnames(heatmap_data)), " "), `[`, 1)
type = sapply(strsplit(as.character(colnames(heatmap_data)), " "), `[`, 2)


bad1 = cts[which(type == "ATAC")][which(!cts[which(type == "ATAC")]%in%cts[which(type == "RNA")])]
bad2 = cts[which(type == "RNA")][which(!cts[which(type == "RNA")]%in%cts[which(type == "ATAC")])]
bads = c(bad1,bad2)
heatmap_data = heatmap_data[,which(!cts%in%bads)]


celltype_annotation$clade = factor(celltype_annotation$clade, levels = rev(c("Gaba","Glut", "NN")))

#heatmap_data = heatmap_data[,which(colMaxs(heatmap_data)>0)]

heatmap_data = heatmap_data[order(rownames(heatmap_data)),]
#heatmap_data = heatmap_data[,order(celltype_annotation$clade)]


options(repr.plot.width=10, repr.plot.height=15)


# Create the heatmap
pheatmap(
  as.matrix(heatmap_data),gaps_col = c(grep("RNA",colnames(heatmap_data))[1]-1),
  cluster_rows = F,
  cluster_cols = F,
  annotation_row = gene_annotation,
  annotation_col = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Age logFC - up in RNA"
)

options(repr.plot.width=12, repr.plot.height=12)

min_val <- -.75  # Adjust based on your data range for increased intensity
max_val <- .75

# Create custom breaks with a narrower range and more intervals for intensity
breaks <- seq(min_val, max_val, length.out = length(heatmap_colors) + 1)


p = pheatmap(
  t(as.matrix(heatmap_data)),gaps_row = c(grep("RNA",colnames(heatmap_data))[1]-1),
  cluster_rows = F,
  cluster_cols = F,
 # annotation_col = gene_annotation,
  annotation_row = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = T,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Age logFC"
)



 pheatmap(
  t(as.matrix(heatmap_data)),gaps_row = c(grep("RNA",colnames(heatmap_data))[1]-1),
  cluster_rows = F,
  cluster_cols = F,
 # annotation_col = gene_annotation,
  annotation_row = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = T,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "TE subfamily Age logFC", file = "~/projects/combined_all/Figures/Figure6-Het-TEs/TE_heatmap_2.pdf", height = 10, width = 10,
)


# Remove rows with too many NAs (e.g., threshold of 50% missing values)
heatmap_data <- heatmap_data[rowMeans(is.na(heatmap_data)) <= 0.5, ]

# Run pheatmap with row clustering enabled
pheatmap(
  t(as.matrix(heatmap_data)),
  gaps_row = c(grep("RNA", colnames(heatmap_data))[1] - 1),
  cluster_rows = TRUE,    # Enable row clustering
  cluster_cols = FALSE,
  annotation_col = gene_annotation,
  annotation_row = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Age logFC"
)


 pheatmap(
  t(as.matrix(heatmap_data)),gaps_row = c(grep("RNA",colnames(heatmap_data))[1]-1),
  cluster_rows = F,
  cluster_cols = F,
 # annotation_col = gene_annotation,
  annotation_row = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Age logFC"
)

pdf("~/projects/combined_all/Figures/Figure6-Het-TEs/TEs_heatmap1.pdf", height = 11, width = 9)
 pheatmap(
  t(as.matrix(heatmap_data)),gaps_row = c(grep("RNA",colnames(heatmap_data))[1]-1),
  cluster_rows = F,
  cluster_cols = F,
 # annotation_col = gene_annotation,
  annotation_row = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Age logFC"
)
dev.off()

while(1>0){dev.off()}

pheatmap(
  t(as.matrix(heatmap_data)),gaps_row = c(grep("RNA",colnames(heatmap_data))[1]-1),
  cluster_rows = F,
  cluster_cols = F,
 # annotation_col = gene_annotation,
  annotation_row = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,height = 11, width = 9,
  main = "Age logFC", filename = "~/projects/combined_all/Figures/Figure6-Het-TEs/TEs_heatmap.pdf"
)

head(heatmap_data_long)
unique(heatmap_data_long$Gene)

head(all)

head(rna_data)

library(ggplot2)
library(dplyr)

# Define the gene of interest and filter for it at the very start
gene <- "L1MA5A"
#sub <- all[which(all$gene ==gene),]# Filter for the specific gene first
sub <- all[which(all$PValue<0.05),]# Filter for the specific gene first

# Separate ATAC and RNA data
atac_data <- sub %>%
  filter(modality == "ATAC") %>%
  rename(logFC_ATAC = logFC)

rna_data <- sub %>%
  filter(modality == "RNA") %>%
  rename(logFC_RNA = logFC)

# Merge ATAC and RNA data by gene, celltype, and clade
merged_data <- inner_join(atac_data, rna_data, by = c("gene", "celltype", "clade"))

# Check if merged data has points to plot
if (nrow(merged_data) == 0) {
  print("No overlapping RNA and ATAC data points for this gene and celltype.")
} else {
  # Create scatter plot
  p <- ggplot(merged_data, aes(x = logFC_RNA, y = logFC_ATAC, color = clade)) +
    geom_point(aes(size = -log10(PValue.x)), alpha = 0.8) +  # Use PValue from either ATAC or RNA, as needed
    theme_minimal() +
    labs(
      title = paste("Scatter Plot of", gene, "Age Difference"),
      x = "log2 Fold Change RNA",
      y = "log2 Fold Change ATAC",
      color = "Clade",
      size = "-log10(PValue)"
    ) + xlim(c(-2,2))+ ylim(c(-.75,.75))+
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 14)
    )+facet_wrap(~

  # Print the plot
  print(p)
}


library(ggplot2)
library(dplyr)

# Ensure 'modality' distinguishes ATAC and RNA properly
gene <- "IAPLTR3-int"
sub <- all %>%
  filter(gene == gene) %>%
  pivot_wider(names_from = modality, values_from = logFC, names_prefix = "logFC_") %>%
  filter(!is.na(logFC_ATAC) & !is.na(logFC_RNA))  # Remove rows with NA in either logFC column

# If `sub` is empty, this indicates no overlap between RNA and ATAC data for the gene.
if (nrow(sub) == 0) {
  print("No overlapping RNA and ATAC data points for this gene.")
} else {
  # Create scatter plot
  p <- ggplot(sub, aes(x = logFC_RNA, y = logFC_ATAC, color = clade)) +
    geom_point(aes(size = -log10(PValue)), alpha = 0.8) +
    theme_minimal() +
    labs(
      title = paste("Scatter Plot of", gene, "Age Difference"),
      x = "log2 Fold Change RNA",
      y = "log2 Fold Change ATAC",
      color = "Clade",
      size = "-log10(PValue)"
    ) +
    scale_color_manual(values = c("clade1_color", "clade2_color", "clade3_color")) +  # Adjust colors for clades
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 14)
    )

  # Print the plot
  print(p)
}


options(repr.plot.width=14, repr.plot.height=10)
gene = "IAPLTR3-int"
sub = all[which(all$gene==gene),]

p <- ggplot(sub, aes(x = celltype, y = -log10(PValue),size= -log10(PValue), color = logFC)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Dot Plot of",gene, "age difference"),
    x = "Cell Type",
    y = "-log10(PValue)",
    color = "logFC"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_wrap(~type+clade, scales = "free") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

# Print the plot
print(p)

library(ggplot2)

# Subset the data for the gene of interest
gene <- "IAPLTR3-int"
sub <- all[all$gene == gene, ]

# Create scatter plot
p <- ggplot(sub, aes(x = logFC_RNA, y = logFC_ATAC, color = clade)) +
  geom_point(aes(size = -log10(PValue)), alpha = 0.8) +
  theme_minimal() +
  labs(
    title = paste("Scatter Plot of", gene, "age difference"),
    x = "log2 Fold Change RNA",
    y = "log2 Fold Change ATAC",
    color = "Clade",
    size = "-log10(PValue)"
  ) +
  scale_color_manual(values = c("clade1_color", "clade2_color", "clade3_color")) +  # Replace with your clade colors
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 14)
  )

# Print the plot
print(p)


# Load necessary library
library(ggplot2)
library(dplyr)

# Set plot dimensions
options(repr.plot.width = 14, repr.plot.height = 10)

# Filter the data for the specific gene
gene = "MuRRS-int"
sub = all[which(all$gene == gene),]

# Create a custom order for cell types within each clade
sub <- sub %>%
  group_by(clade) %>%
  mutate(celltype_order = factor(celltype, levels = unique(celltype))) %>%
  ungroup()

# Create the dot plot
p <- ggplot(sub, aes(x = celltype_order, y = -log10(PValue), size = -log10(PValue), color = logFC)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Dot Plot of", gene, "age difference"),
    x = "Cell Type",
    y = "-log10(PValue)",
    color = "logFC"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_wrap(~ type + clade, scales = "free") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

# Print the plot
print(p)


options(repr.plot.width=12, repr.plot.height=7)
# Filter the data for the specific gene
gene = "MuRRS-int"
sub = all[which(all$gene == gene),]

# Create a custom order for cell types within each clade
sub <- sub %>%
  group_by(clade) %>%
  mutate(celltype_order = factor(celltype, levels = unique(celltype))) %>%
  ungroup()

p <- ggplot(sub, aes(x = celltype_order, y = -log10(PValue), size = -log10(PValue), color = logFC)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Dot Plot of", gene, "age difference"),
    x = "Cell Type",
    y = "-log10(PValue)",
    color = "logFC"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_grid(type ~ clade, scales = "free", space = "free_x") +  # "free_x" allows different clade widths
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 

p

options(repr.plot.width=12, repr.plot.height=7)
# Filter the data for the specific gene
gene = "L1MA5A"
sub = all[which(all$gene == gene),]

# Create a custom order for cell types within each clade
sub <- sub %>%
  group_by(clade) %>%
  mutate(celltype_order = factor(celltype, levels = unique(celltype))) %>%
  ungroup()

p <- ggplot(sub, aes(x = celltype_order, y = -log10(PValue), size = -log10(PValue), color = logFC)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Dot Plot of", gene, "age difference"),
    x = "Cell Type",
    y = "-log10(PValue)",
    color = "logFC"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_grid(type ~ clade, scales = "free", space = "free_x") +  # "free_x" allows different clade widths
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 

p

options(repr.plot.width=7, repr.plot.height=10)

gene = "MuRRS-int"
sub = all[which(all$gene == gene),]

# Create a custom order for cell types within each clade
sub <- sub %>%
  group_by(clade) %>%
  mutate(celltype_order = factor(celltype, levels = unique(celltype))) %>%
  ungroup()
p <- ggplot(sub, aes(x = celltype_order, y = abs(logFC), size = -log10(PValue), color = logFC)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Dot Plot of", gene, "age difference"),
    x = "Cell Type",
    y = "LogFC (18mo/2mo)",
    color = "logFC"
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Adds horizontal dotted line at y=0
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) + coord_flip()+
  facet_grid(clade ~ type, scales = "free", space = "free_x") +  # "free_x" allows different clade widths
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,limits = c(-1, 1),oob = scales::squish) 

p

options(repr.plot.width=10, repr.plot.height=6)

gene = "IAPLTR3-int"
sub = all[which(all$gene == gene),]

# Create a custom order for cell types within each clade
sub <- sub %>%
  group_by(clade) %>%
  mutate(celltype_order = factor(celltype, levels = unique(celltype))) %>%
  ungroup()
p <- ggplot(sub, aes(x = celltype_order, y = logFC, size = -log10(PValue), color = logFC)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Dot Plot of", gene, "age difference"),
    x = "Cell Type",
    y = "LogFC (18mo/2mo)",
    color = "logFC"
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +  # Adds horizontal dotted line at y=0
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_grid(type ~ clade, scales = "free", space = "free_x") +  # "free_x" allows different clade widths
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,limits = c(-.75, .75),oob = scales::squish) 

p

setwd("~/projects/combined_all/Figures/Figure6-Het-TEs/")

ggsave(p, height = 6,width =10, file = "MuRRS-int.pdf")

ggsave(p, height = 6,width =10, file = "IAPLTR3-int.pdf")

# Step 1: Calculate y-axis limits for each type
library(dplyr)

# Find y-axis limits for each type (e.g., RNA and ATAC)
y_limits <- sub %>%
  group_by(type) %>%
  summarise(
    y_min = min(-log10(PValue), na.rm = TRUE),
    y_max = max(-log10(PValue), na.rm = TRUE)
  )

# Step 2: Create separate plots for each type with their specific limits
p <- ggplot(sub, aes(x = celltype_order, y = -log10(PValue), size = -log10(PValue), color = logFC)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Dot Plot of", gene, "age difference"),
    x = "Cell Type",
    y = "-log10(PValue)",
    color = "logFC"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_grid(type ~ clade, scales = "free_x", space = "free") +  # "free_x" allows different clade widths
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_y_continuous(
    limits = function(limits) {
      # Set y limits based on the "type"
      if (sub$type[1] == "RNA") {
        c(y_limits$y_min[y_limits$type == "RNA"], y_limits$y_max[y_limits$type == "RNA"])
      } else {
        c(y_limits$y_min[y_limits$type == "ATAC"], y_limits$y_max[y_limits$type == "ATAC"])
      }
    }
  )




p

# Load necessary library
library(ggplot2)
library(dplyr)

# Set plot dimensions
options(repr.plot.width = 14, repr.plot.height = 10)

# Filter the data for the specific gene
gene = "MuRRS-int"
sub = all[which(all$gene == gene),]

# Create a custom order for cell types within each clade
sub <- sub %>%
  group_by(clade) %>%
  mutate(celltype_order = factor(celltype, levels = unique(celltype))) %>%
  ungroup()

# Create the dot plot
p <- ggplot(sub, aes(x = celltype_order, y = -log10(PValue), size = -log10(PValue), color = logFC)) +
  geom_point() +
  theme_minimal() +
  labs(
    title = paste("Dot Plot of", gene, "age difference"),
    x = "Cell Type",
    y = "-log10(PValue)",
    color = "logFC"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  facet_wrap(~ type + clade, scales = "free_x") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

# Print the plot
print(p)



options(repr.plot.width=15, repr.plot.height=12)
# Reshape the data for the heatmap
sig_tab = rna[which(rna$gene %in% c(sig_genes)),]
df = sig_tab
df$gene = gsub("SoloTE-", "", df$gene)
#df$gene = gsub("-", "_", df$gene)
heatmap_data <- dcast(df, gene ~ celltype, value.var = "logFC")
row.names(heatmap_data) <- heatmap_data$gene
heatmap_data <- heatmap_data[,-1]

df$type = ukey[paste(df$gene), "V6"]

# Annotations
gene_annotation <- df[, c("gene", "type")]
gene_annotation <- unique(gene_annotation)
row.names(gene_annotation) <- gene_annotation$gene
gene_annotation <- gene_annotation[,-1, drop = FALSE]

celltype_annotation <- df[, c("celltype", "clade")]
celltype_annotation <- unique(celltype_annotation)
row.names(celltype_annotation) <- celltype_annotation$celltype
celltype_annotation <- celltype_annotation[,-1, drop = FALSE]

# Colors for annotations
ann_colors <- list(
  type = brewer.pal(8, "Set2")[1:length(unique(gene_annotation$type))],
  clade = brewer.pal(8, "Set1")[1:length(unique(celltype_annotation$clade))]
)
names(ann_colors$type) <- unique(gene_annotation$type)
names(ann_colors$clade) <- unique(celltype_annotation$clade)



# Define the color palette for the heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks for the heatmap so that 0 is the midpoint
max_abs_val <- max(abs(heatmap_data), na.rm = TRUE)
breaks <- seq(-max_abs_val, max_abs_val, length.out = 101)

# Create the heatmap
pheatmap(
  as.matrix(heatmap_data),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = gene_annotation,
  annotation_col = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Heatmap of logFC"
)


options(repr.plot.width=15, repr.plot.height=12)
# Reshape the data for the heatmap
sig_tab = rna[which(rna$gene %in% c(sig_genes)),]
df = sig_tab
df$gene = gsub("SoloTE-", "", df$gene)
#df$gene = gsub("-", "_", df$gene)
heatmap_data <- dcast(df, gene ~ celltype, value.var = "logFC")
row.names(heatmap_data) <- heatmap_data$gene
heatmap_data <- heatmap_data[,-1]

df$type = ukey[paste(df$gene), "V6"]

# Annotations
gene_annotation <- df[, c("gene", "type")]
gene_annotation <- unique(gene_annotation)
row.names(gene_annotation) <- gene_annotation$gene
gene_annotation <- gene_annotation[,-1, drop = FALSE]

celltype_annotation <- df[, c("celltype", "clade")]
celltype_annotation <- unique(celltype_annotation)
row.names(celltype_annotation) <- celltype_annotation$celltype
celltype_annotation <- celltype_annotation[,-1, drop = FALSE]

# Colors for annotations
ann_colors <- list(
  type = brewer.pal(8, "Set2")[1:length(unique(gene_annotation$type))],
  clade = brewer.pal(8, "Set1")[1:length(unique(celltype_annotation$clade))]
)
names(ann_colors$type) <- unique(gene_annotation$type)
names(ann_colors$clade) <- unique(celltype_annotation$clade)



# Define the color palette for the heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks for the heatmap so that 0 is the midpoint
max_abs_val <- max(abs(heatmap_data), na.rm = TRUE)
breaks <- seq(-max_abs_val, max_abs_val, length.out = 101)

# Create the heatmap
pheatmap(
  as.matrix(heatmap_data),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = gene_annotation,
  annotation_col = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Heatmap of logFC"
)
