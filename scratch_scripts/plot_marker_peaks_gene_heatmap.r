library(data.table)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(enrichR)
websiteLive <- getOption("enrichR.live")
dbs <- c( "GO_Biological_Process_2023")

meta = read.csv("../ABC/gene_meta.csv")
rownames(meta) = meta$`geneID.1`

marker_peaks = read.csv("../marker_peaks.001.csv", header = F)

marker_peaks$celltype = gsub(" ", "_",marker_peaks$V1)
nrow(marker_peaks)

files = list.files(".", "8wk")

marker_peaks$celltype = gsub("2-3", "23", marker_peaks$celltype)
marker_peaks$celltype = gsub("6b-", "6b", marker_peaks$celltype)
marker_peaks$celltype = gsub("D12", "D1", marker_peaks$celltype)

files

all_mp_genes = list()
all_mp_locs = list()
all_mp = list()
for(f in files){
    f2 = gsub("8wk", "18mo", f)
    t = read.csv(f)
    t$loc = paste(sapply(strsplit(as.character(t$X), "-"), `[`, 1), "-", sapply(strsplit(as.character(t$X), "-"), `[`, 2), sep ="")
    head(t)

    t$gene = sapply(strsplit(as.character(t$X), "-"), `[`, 3)
    t$gene_name = meta[paste(t$gene), "gene_name"]
    rownames(t) = t$X
    t = t[order(t$abc_score, decreasing = T),]
    #t = t[which(!duplicated(t$gene_name)),]
    t = t[which(!duplicated(t$loc)),]

    t = t[which(t$abc_score>0.05),]
    
    ct = gsub("filtered.", "", f)
    ct = gsub(".8wk.abc_score.csv", "", ct)
    
    t2 = read.csv(f2)
    t2$loc = paste(sapply(strsplit(as.character(t2$X), "-"), `[`, 1), "-", sapply(strsplit(as.character(t2$X), "-"), `[`, 2), sep ="")
    head(t2)

    t2$gene = sapply(strsplit(as.character(t2$X), "-"), `[`, 3)
    t2$gene_name = meta[paste(t2$gene), "gene_name"]
    rownames(t2) = t2$X
    t2 = t2[order(t2$abc_score, decreasing = T),]
    #t = t[which(!duplicated(t$gene_name)),]
    t2 = t2[which(!duplicated(t2$loc)),]

    t2 = t2[which(t2$abc_score>0.05),]
    
    
    if(!ct %in% unique(marker_peaks$celltype)) {
        cat(ct, " not found\n")
    }
    t = rbind(t,t2)
    t = t[order(t$abc_score, decreasing = T),]
    t = t[which(!duplicated(t$loc)),]

    
    mp = marker_peaks[which(marker_peaks$celltype==ct),]
    t[which(t$loc%in% mp$V2),]
    t$celltype = ct
    if (!ct %in% all_mp){
        all_mp_locs[[ct]] = unique(t[which(t$loc%in% mp$V2),"loc"])
        all_mp_genes[[ct]] = unique(t[which(t$loc%in% mp$V2),"gene_name"])
        all_mp[[ct]] = t[which(t$loc%in% mp$V2),c("celltype", "gene_name", "loc")]
    }
}

all = do.call(rbind,all_mp)
head(all)
all$gene_name

#write.table(collapsed_table, "connected_marker_genes.001.txt", sep = "\t")

obj = readRDS("../female_RNA/RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_age = paste(obj$celltype_final, obj$age )


Idents(obj) = "celltype_final"
av = AverageExpression(obj)

head(unique(collapsed_table$id))
head(colnames(av$RNA))

table(all$celltype)

length(unique(collapsed_table$element))

tab= av$RNA
colnames(tab)= gsub(" ", "_", colnames(tab))
colnames(tab)= gsub("/", "", colnames(tab))

tab = tab[all$gene_name,unique(all$celltype)]

tab = tab[all$gene_name,gsub("D1_G", "D12_G",unique(all$celltype))]

pheatmap(tab, main = nrow(tab),cluster_rows = F, cluster_cols = F,  scale = "row", show_rownames = F) 

library(viridis)

# Generate the heatmap with viridis color palette
ph = pheatmap(
  tab, 
  main = nrow(tab), 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,  
  scale = "row", 
  show_rownames = FALSE,
  color = viridis(100)  # Use viridis color palette with 100 color breaks
)

pdf("out_001_RNA_marker_heat.pdf")
print(ph)
dev.off()

pheatmap(log(tab+1), cluster_rows = F, cluster_cols = F,  scale = "none", show_rownames = F) 

library(zellkonverter)

h = readH5AD("../celltype_RPM.h5ad")

count_table <-h@assays@data$X
colnames(count_table)

tab= av$RNA
colnames(tab)= gsub(" ", "_", colnames(tab))
tab = tab[all$gene_name,unique(all$celltype)]

tab= count_table
colnames(tab)= gsub(" ", "_", colnames(tab))
colnames(tab) = gsub("2-3", "23", colnames(tab))
colnames(tab) = gsub("6b-", "6b", colnames(tab))
tab = tab[all$loc,unique(all$celltype)]

colnames(tab)[which(!colnames(tab)%in%unique(all$celltype))]

library(viridis)

# Generate the heatmap with viridis color palette
pheatmap(
  tab, 
  main = nrow(tab), 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,  
  scale = "row", 
  show_rownames = FALSE,
  color = viridis(100)  # Use viridis color palette with 100 color breaks
)

ph = pheatmap(
  tab, 
  main = nrow(tab), 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,  
  scale = "row", 
  show_rownames = FALSE,
  color = viridis(100)  # Use viridis color palette with 100 color breaks
)

pdf("out_001_ATAC_marker_heat.pdf")
print(ph)
dev.off()

while(3>1){dev.off()}

pheatmap(log(tab+1), cluster_rows = F, cluster_cols = F,  scale = "none", show_rownames = F) 


