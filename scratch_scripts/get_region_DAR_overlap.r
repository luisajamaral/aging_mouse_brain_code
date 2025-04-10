library(ggplot2)
library(data.table)
library(UpSetR)
library(gridExtra)
setwd("../../h5ads_final/")

install.packages("UpSetR")

meta = read.csv("../final_meta.csv")
meta$celltype_final = gsub("/","-", meta$celltype_final)
meta$celltype_final = gsub("Lymphoid","Lymphoid NN", meta$celltype_final)
meta = meta[-which(meta$celltype_final=="doublet"),]
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)


meta$celltype_region = paste(meta$celltype_final, meta$region, sep = ":")
meta$celltype_region=gsub("Astro-TE", "Astro", meta$celltype_region)
meta$celltype_region=gsub("Astro-NT", "Astro", meta$celltype_region)
meta$celltype_region=gsub(" ", "_", meta$celltype_region)
meta$age = factor(meta$age , levels = c("2mo", "9mo", "18mo"))

setwd("region_DARs/")

head(me)

unique(tab$id)

p1 = ggplot(dtab, aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "DARs per Region",
       x = "Cell Type",
       y = "Number of DARs") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))
p1

cl = "Astro-"
files = list.files(".", cl)
files

options(repr.plot.width=14, repr.plot.height=5)
cl = "Astro_NN"
files = list.files(".", cl)
all = list()
for(f in files){
    curr = read.csv(f)
    id = gsub("diff_peaks_" , "", f)
    id = gsub("_2vs18_iter.csv" , "", id)

    curr$id = id 
    all[[id]] = curr
    }
tab = do.call(rbind, all)
tab=as.data.frame(tab)
tab$sig = "no"
tab[which(tab$adjusted.p.value<0.05 & tab$log2.fold_change.<0 ),"sig"] = "up"
tab[which(tab$adjusted.p.value<0.05 & tab$log2.fold_change.>0 ),"sig"] = "down"


dtab = tab[which(tab$adjusted.p.value<0.05),]
me =meta[grep(cl, meta$celltype_region),]
ord= names(table(me$celltype_region))[order(table(me$celltype_region))]
tab$id = factor(tab$id, levels = ord)
dtab$id = factor(dtab$id, levels = ord)
me$celltype_region = factor(me$celltype_region, levels = ord)

p1 = ggplot(dtab, aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "DARs per Region",
       x = "Cell Type",
       y = "Number of DARs") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))
p2 = ggplot(me, aes(x = celltype_region, fill= age)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "Cell # per Region",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))
p3 = ggplot(tab[which(abs(tab$log2.fold_change.) > 0.15),], 
  aes(x = id, y = -`log2.fold_change.`, 
  color = ifelse(`adjusted.p.value` < 0.05, ifelse(`log2.fold_change.` < 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age DARs logFC",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip()+theme(text = element_text(color = "black", size = 13))

grid.arrange(p2, p1,p3, ncol = 3)

RLP = dtab$feature.name[which(dtab$id == paste(cl,':RLP', sep = ""))]
AMY = dtab$feature.name[which(dtab$id == paste(cl,':AMY', sep = ""))]
CP = dtab$feature.name[which(dtab$id == paste(cl,':CP', sep = ""))]
ENT = dtab$feature.name[which(dtab$id == paste(cl,':ENT', sep = ""))]
FC = dtab$feature.name[which(dtab$id == paste(cl,':FC', sep = ""))]
HCP = dtab$feature.name[which(dtab$id == paste(cl,':HCP', sep = ""))]
HCA = dtab$feature.name[which(dtab$id == paste(cl,':HCA', sep = ""))]
NAC = dtab$feature.name[which(dtab$id == paste(cl,':NAC', sep = ""))]


listInput <- list(RLP=RLP, AMY=AMY, CP=CP, ENT=ENT,FC=FC, HCP=HCP,HCA=HCA,NAC=NAC)
options(repr.plot.width=14, repr.plot.height=7)

upset(fromList(listInput), order.by = "freq", nsets=8,text.scale = 2)

options(repr.plot.width=14, repr.plot.height=5)
cl = "Astro-"
files = list.files(".", cl)
all = list()
for(f in files){
    curr = read.csv(f)
    id = gsub("diff_peaks_" , "", f)
    id = gsub("_2vs18_iter.csv" , "", id)

    curr$id = id 
    all[[id]] = curr
    }
tab = do.call(rbind, all)
tab=as.data.frame(tab)
tab$sig = "no"
tab[which(tab$adjusted.p.value<0.05 & tab$log2.fold_change.<0 ),"sig"] = "up"
tab[which(tab$adjusted.p.value<0.05 & tab$log2.fold_change.>0 ),"sig"] = "down"


dtab = tab[which(tab$adjusted.p.value<0.05),]
me =meta[grep(cl, meta$celltype_region),]
ord= names(table(me$celltype_region))[order(table(me$celltype_region))]
tab$id = factor(tab$id, levels = ord)
dtab$id = factor(dtab$id, levels = ord)
me$celltype_region = factor(me$celltype_region, levels = ord)

p1 = ggplot(dtab, aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "DARs per Region",
       x = "Cell Type",
       y = "Number of DARs") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))
p2 = ggplot(me, aes(x = celltype_region, fill= age)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "Cell # per Region",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))
p3 = ggplot(tab[which(abs(tab$log2.fold_change.) > 0.15),], 
  aes(x = id, y = -`log2.fold_change.`, 
  color = ifelse(`adjusted.p.value` < 0.05, ifelse(`log2.fold_change.` < 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age DARs logFC",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip()+theme(text = element_text(color = "black", size = 13))

grid.arrange(p2, p1,p3, ncol = 3)

RLP = dtab$feature.name[which(dtab$id == paste(cl,':RLP', sep = ""))]
AMY = dtab$feature.name[which(dtab$id == paste(cl,':AMY', sep = ""))]
CP = dtab$feature.name[which(dtab$id == paste(cl,':CP', sep = ""))]
ENT = dtab$feature.name[which(dtab$id == paste(cl,':ENT', sep = ""))]
FC = dtab$feature.name[which(dtab$id == paste(cl,':FC', sep = ""))]
HCP = dtab$feature.name[which(dtab$id == paste(cl,':HCP', sep = ""))]
HCA = dtab$feature.name[which(dtab$id == paste(cl,':HCA', sep = ""))]
NAC = dtab$feature.name[which(dtab$id == paste(cl,':NAC', sep = ""))]


listInput <- list(RLP=RLP, AMY=AMY, CP=CP, ENT=ENT,FC=FC, HCP=HCP,HCA=HCA,NAC=NAC)
options(repr.plot.width=14, repr.plot.height=7)

upset(fromList(listInput), order.by = "freq", nsets=8,text.scale = 2)

which(dtab$id == paste(cl,':CP', sep = ""))

listInput <- list( CP=CP, NAC=NAC)
options(repr.plot.width=14, repr.plot.height=7)

upset(fromList(listInput), order.by = "freq", nsets=8,text.scale = 2,show.numbers = TRUE)

NAC

setwd("region_DARs")

files = list.files(".", "diff_peaks")

for(f in files){
    # Read DAR file
    dar_data <- read.table(f, header=TRUE, sep=',', stringsAsFactors=FALSE)
    
    # Filter rows with significant changes
    significant_changes <- subset(dar_data, adjusted.p.value < 0.05)
    
    # Split into upregulated and downregulated based on log2.fold_change
    upregulated_changes <- subset(significant_changes, log2.fold_change. > 0)
    downregulated_changes <- subset(significant_changes, log2.fold_change. < 0)
    
    # Function to extract chromosome, start, and end from feature.name
    extract_coordinates <- function(feature_name) {
      parts <- unlist(strsplit(feature_name, "[:-]"))
      chr <- parts[1]
      start <- as.numeric(parts[2])
      end <- as.numeric(parts[3])
      return(c(chr, start, end))
    }
    
    # Apply the function to create BED files
    upregulated_bed <- t(sapply(upregulated_changes$feature.name, extract_coordinates))
    downregulated_bed <- t(sapply(downregulated_changes$feature.name, extract_coordinates))

    ct = gsub("diff_peaks_", "", f)
    ct = gsub("_iter.csv", "",ct)
    # Save BED files
    if(nrow(upregulated_bed)>1){
        write.table(upregulated_bed, paste(ct,'_up.bed', sep = ""), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
    if(nrow(downregulated_bed)>1){
        write.table(downregulated_bed, paste(ct, '_down.bed', sep=""), sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
}

##for i in *.bed ; do findMotifsGenome.pl $i mm10 ${i%.bed}_motifs_bg -size 250  -p 15 -nomotif -bg ../../all_peaks/${i%:*}_allpeaks.bed ; done
## remember the changes are reversed!

for i in Astro_NN*.bed ; do findMotifsGenome.pl $i mm10 ${i%.bed}_motifs_bg -size 250  -p 15 -nomotif -bg ../../all_peaks/Astro-TE_allpeaks.bed ; done
## remember the changes are reversed!

setwd("../female_RNA/DEG_results_region/") 

files = list.files(".", ".csv")

files

for(f in files){
    # Read DAR file
    dar_data <- read.table(f, header=TRUE, sep=',', stringsAsFactors=FALSE)
    
    # Filter rows with significant changes
    significant_changes <- subset(dar_data, p_val_adj < 0.05)
    
    # Split into upregulated and downregulated based on log2.fold_change
    upregulated_changes <- subset(significant_changes, avg_log2FC < 0)
    downregulated_changes <- subset(significant_changes, avg_log2FC > 0)
    
   
    
    # Apply the function to create BED files
    upregulated_bed <- upregulated_changes$X
    downregulated_bed <- downregulated_changes$X
    all = significant_changes$X
    ct = gsub(".csv", "", f)
    # Save BED files
    if(length(upregulated_bed)>1){
        write.table(upregulated_bed, paste(ct,'_up.txt', sep = ""), sep='\n', col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
    if(length(downregulated_bed)>1){
        write.table(downregulated_bed, paste(ct, '_down.txt', sep=""), sep='\n', col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
    write.table(all, paste(ct, '_all.txt', sep=""), sep='\n', col.names=FALSE, row.names=FALSE, quote=FALSE)

}

##for i in *down.txt; do findGO.pl $i mouse ${i%.txt} -bg ${i%down.txt}all.txt -cpu 5 ; done

library(Seurat)

obj = readRDS("../../female_RNA/RNA_final_SCT.RDS")

obj

meta = obj@meta.data

head(meta)
meta$celltype_region = paste(meta$celltype_final, meta$region, sep = "--")
meta$celltype_region = gsub(" ", "_", meta$celltype_region)

rmeta = meta

rmeta$celltype_region=gsub("Astro-TE", "Astro", rmeta$celltype_region)
rmeta$celltype_region=gsub("Astro-NT", "Astro", rmeta$celltype_region)



files = list.files(".", "STR")
files

setwd("~/projects/combined_all/female_RNA/DEG_results_region/")
cl = "Astro_NN"
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

dtab = tab[which(tab$p_val_adj<0.05 & abs(tab$avg_log2FC)>logfc),]
me =rmeta[grep(cl, rmeta$celltype_region),]
bc = names(table(me$celltype_region)[which(table(me$celltype_region)>100)])
me = me[which(me$celltype_region%in%bc),]
ord= names(table(me$celltype_region))[order(table(me$celltype_region))]
tab$id = factor(tab$id, levels = ord)
dtab$id = factor(dtab$id, levels = ord)
me$celltype_region = factor(me$celltype_region, levels = ord)

p1 = ggplot(dtab, aes(x = id, fill = sig)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "DEGs per Region",
       x = "Cell Type",
       y = "Number of DEGs") +
  coord_flip()+theme(text = element_text(color = "black", size = 13))

p2 = ggplot(me, aes(x = celltype_region, fill= age)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
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


head(dtab[which(dtab$id=="Microglia_NN--RLP"),])

obj

library(Seurat)

options(repr.plot.width=23, repr.plot.height=5)

VlnPlot(obj, group.by= "celltype_final", split.by="age", "Cst3", pt.size = 0,split.plot = TRUE)

head(obj$celltype_final)

library(Seurat)


filtered_seurat <- subset(obj, cells = which(obj$celltype_final == "Microglia NN"))
filtered_seurat
filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat, nfeatures = 5000)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(filtered_seurat))
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:40)  # You can adjust the number of dimensions
filtered_seurat <- FindClusters(filtered_seurat, resolution = 1)  # Adjust the resolution as needed
filtered_seurat <- RunUMAP(filtered_seurat, reduction = "pca", dims = 1:40)  # Adjust the number of dimensions as needed


filtered_seurat <- FindClusters(filtered_seurat, resolution = .5)  # Adjust the resolution as needed


options(repr.plot.width=6, repr.plot.height=6)
DimPlot(filtered_seurat)

options(repr.plot.width=12, repr.plot.height=6)

DimPlot(filtered_seurat, group.by = "age")

FeaturePlot(filtered_seurat, "Nsmce2",split.by = "age")

options(repr.plot.width=16, repr.plot.height=6)

DimPlot(filtered_seurat, group.by = "age",split.by = "region")
DimPlot(filtered_seurat, group.by = "age",split.by = "rep")

options(repr.plot.width=12, repr.plot.height=12)

DimPlot(filtered_seurat, split.by = "sample",ncol = 6, group.by = "age")


ggplot(filtered_seurat@meta.data, aes(x = sample, fill = seurat_clusters)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "DEGs per Region",
       x = "Cell Type",
       y = "Number of DEGs") +
  coord_flip()


library(dplyr)

sample_names <- unique(filtered_seurat$sample)
genes = rownames(filtered_seurat@assays$RNA$counts)
zero_matrix <- matrix(0, nrow = length(genes), ncol = length(sample_names))
rownames(zero_matrix) <- genes
colnames(zero_matrix) <- sample_names

for (sample in sample_names[2:length(sample_names)]) {
  subset_indices <- which(filtered_seurat$sample == sample)
  zero_matrix[, sample] <- rowSums(filtered_seurat@assays$RNA$counts[, subset_indices, drop = FALSE])
}

zero_matrix = zero_matrix[,which(!is.na(colnames(zero_matrix)))]

tpm = (zero_matrix/colSums(zero_matrix))*1e6

dim(tpm)

length(unique(dtab$X))

head(tpm[unique(dtab$X),])

library(pheatmap)



head(dtab)



options(repr.plot.width=8, repr.plot.height=6)
rlp = tpm[unique(dtab$X[which(dtab$id=="Microglia_NN--RLP")]),grep("^RLP", colnames(tpm))]
pheatmap(rlp, show_rownames = F, scale = "row")

unique(dtab$id)

filtered_seurat$age_celltype <- gsub(" ", "_", paste(filtered_seurat$age, filtered_seurat$celltype_final, sep = "_"))
filtered_seurat$age_celltype_region <- gsub(" ", "_", paste(filtered_seurat$age, filtered_seurat$celltype_final, filtered_seurat$region, sep = "_"))


unique(filtered_seurat$age_celltype_region)

Idents(filtered_seurat) = "age_celltype"

o = FindMarkers(filtered_seurat, `ident.1` = "18mo_Microglia_NN", `ident.2` = "2mo_Microglia_NN", test.use = "MAST")

nrow(rlp)

options(repr.plot.width=8, repr.plot.height=6)
rlp = tpm[rownames(o)[which(o$p_val_adj<0.01)],]
pheatmap(rlp, show_rownames = F, scale = "row")

o = FindMarkers(filtered_seurat, `ident.1` = "18mo_Microglia_NN_RLP", `ident.2` = "2mo_Microglia_NN_RLP", )

options(repr.plot.width=8, repr.plot.height=6)
rlp = tpm[rownames(o)[which(o$p_val_adj<0.01)],grep("^RLP", colnames(tpm))]
pheatmap(rlp, show_rownames = F, scale = "row")

head(o)

m =FindMarkers(filtered_seurat, `ident.1` = "18mo_Microglia_NN_RLP", `ident.2` = "2mo_Microglia_NN_RLP",test.use = "MAST" )

options(repr.plot.width=8, repr.plot.height=6)
rlp = tpm[rownames(m)[which(m$p_val_adj<0.01)],grep("^RLP", colnames(tpm))]
pheatmap(rlp, show_rownames = F, scale = "row")

head(m)

library(edgeR)

BiocManager::install("edgeR")


library(edgeR)

head(zero_matrix)

old = zero_matrix
zero_matrix = zero_matrix[,-grep("9mo", colnames(zero_matrix))]

reg = sapply(strsplit(as.character(colnames(zero_matrix)), "_"), "[[", 1)
age = sapply(strsplit(as.character(colnames(zero_matrix)), "_"), "[[", 2)
rep = sapply(strsplit(as.character(colnames(zero_matrix)), "_"), "[[", 3)



y <- DGEList(counts=zero_matrix,group=age)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)


age

reg

design <- model.matrix(~reg+age)
y <- estimateDisp(y,design)


head(design)

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit, contrast = c(1,0,0,0,0,0,0,0,-1))
topTags(qlf)

res = qlf$table[which(qlf$table$PValue<0.01),]

nrow(res)

options(repr.plot.width=8, repr.plot.height=6)
rlp = tpm[rownames(res),grep("^RLP", colnames(tpm))]
pheatmap(rlp, show_rownames = F, scale = "row")

options(repr.plot.width=8, repr.plot.height=6)
rlp = tpm[rownames(res),grep("^RLP", colnames(tpm))]
pheatmap(rlp, show_rownames = F, scale = "row")

obj

DimPlot(obj, label =T)

options(repr.plot.width=18, repr.plot.height=6)

DimPlot(obj, group.by = "celltype_final",label =T)

DimPlot(obj, group.by = "best_celltype_from_ATAC",label =T)

obj=readRDS("../female_RNA/RNA_final.RDS")

sub = subset(obj, cells = which(obj$seurat_clusters %in%c(7,52, 35,50)))
sub <-  SCTransform(sub)
    sub <- RunPCA(sub, features = VariableFeatures(sub))
    sub <- FindNeighbors(sub, dims = 1:45)  # You can adjust the number of dimensions
    sub <- FindClusters(sub, resolution = 2, )  # Adjust the resolution as needed
    sub <- RunUMAP(sub, reduction = "pca", dims = 1:45)  # Adjust the number of dimensions as needed
    sub
    sub@meta.data$sub_leiden = paste(cl, sub$seurat_clusters)
  
    print(DimPlot(sub, group.by = "seurat_clusters", reduction = "umap",label = T))
    print(DimPlot(sub, group.by = "best_celltype_from_ATAC",reduction = "umap", label = T))
    print(DimPlot(sub, group.by = "celltype_final", reduction = "umap",label = T))
     print(FeaturePlot(sub, reduction = "umap","log1p_total_counts", max.cutoff = 2000))
   # dev.off()

  

sub <- FindClusters(sub, resolution = .5 )  # Adjust the resolution as needed

print(DimPlot(sub, group.by = "seurat_clusters", reduction = "umap",label = T))


fm = FindMarkers(sub,`ident.1` = 10)

fm3 = FindMarkers(sub,`ident.1` = 3)

head(fm3)

head(fm)

    print(DimPlot(sub, group.by = "celltype_final",split.by = "age",reduction = "umap", label = T))


    print(DimPlot(sub, group.by = "best_celltype_fixed",reduction = "umap", label = T))


FeaturePlot(sub,"Clic6", split.by = "age")


