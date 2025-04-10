library(ggplot2)
library(data.table)
library(UpSetR)
library(pheatmap)
library(Seurat)
library(enrichR)


obj

obj = readRDS("RNA_final_SCT.RDS")

obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))

head(obj@meta.data$orig.ident)

obj@meta.data$region = sapply(strsplit(as.character(obj$orig.ident), "_"), `[`, 1)

options(repr.plot.width=20, repr.plot.height=6)

FeaturePlot(obj, "Lypd6b", split.by = "age")

obj$celltype_final[which(obj$seurat_clusters==46)

obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)

Idents(obj)= "celltype_region_age"

ag = AggregateExpression(obj,assays = "SCT", return.seurat = T, group.by = c("celltype_final", "region", "age", "orig.ident"))


ag$celltype.age <- paste(ag$celltype_final, ag$age, sep = "_")


Idents(ag) = "celltype_final"

VlnPlot(ag, features = c("Snca", "Mal"), idents = "Oligo NN", group.by = "age" , split.by = "region") 


aoligo = subset(ag, cells = which(ag$celltype_final=="Oligo NN"))

aoligo=RunPCA(aoligo, features = VariableFeatures(obj))

head(ag@meta.data)
ag$sample = sapply(strsplit(as.character(ag$orig.ident), "mo_"), `[`, 2)

options(repr.plot.width=20, repr.plot.height=5)

DimPlot(aoligo, group.by = "sample",label = T)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

dir.create("DESeq2_pseudobulk")
setwd("DESeq2_pseudobulk")

write.table(fm, file = "Oligo_combined.txt")

ag$celltype.age.region = paste(ag$celltype.age, ag$region)

r

all_genes = genes
genes_list = list()
Idents(ag) <- "celltype.age.region"
for(r in unique(ag$region)) {
fm <- FindMarkers(object = ag, 
                         ident.1 = paste("Oligo NN_2mo",r ), 
                         ident.2 = paste("Oligo NN_18mo",r ),
                         test.use = "DESeq2",min.cells.group=2 )
write.table(fm, file = paste("Oligo_", r, "_2vs18.txt", sep = ""))
head(fm, n = 15)
    all_genes = c(all_genes, rownames(fm)[which(fm$p_val_adj<0.05)])
    genes_list[[r]] = rownames(fm)[which(fm$p_val_adj<0.05)]
}

upset(fromList(genes_list), order.by = "freq", nsets = 20
     )


length(all_genes)

all_genes = (unique(all_genes))

fm[grep("", rownames(fm)),]

genes = rownames(fm)[which(fm$p_val_adj<0.05)]

genes = all_genes

DEGS = genes

obj$celltype_final[grep("Astro", obj$celltype_final)] = "Astro NN"

obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)

Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)


obj$celltype_sample = paste(obj$celltype_final, obj$orig.ident)
Idents(obj)= "celltype_sample"
av = AverageExpression(obj)


setwd("female_RNA/")

files = list.files("DEG_results_region/", ".csv", full.names = T)
cfiles = list.files("DEG_results_.01_.01/", ".csv", full.names = T)


cl = "Micro"
afiles = c(files,cfiles)
afiles = afiles[grep(cl, afiles)]
atabs= list()
list_genes = list()

for(f in afiles){
    # Read DAR file
    dar_data <- read.table(f, header=TRUE, sep=',', stringsAsFactors=FALSE)
    dar_data$cl= gsub(".csv","",f)
    dar_data$cl= gsub("DEG_results_.01_.01//","Combined_",dar_data$cl)
    dar_data$cl= gsub("DEG_results_region//","",dar_data$cl)

    # Filter rows with significant changes
    significant_changes <- subset(dar_data, p_val_adj < 0.05 & abs(avg_log2FC)>0.2 & (`pct.1`>0.05 | `pct.2`>0.05))
    atabs[[f]] = significant_changes
    #list_genes[[dar_data$cl[1]]] = rownames(significant_changes)
    list_genes[[dar_data$cl[1]]] = significant_changes$X

    
}
tabs = do.call(rbind, atabs)
tabs = tabs[which(tabs$p_val_adj<0.05& abs(tabs$avg_log2FC)>0.2 & (tabs$`pct.1`>0.05 | tabs$`pct.2`>0.05)),]
genes = unique(tabs$X)
length(genes)
upset(fromList(list_genes), order.by = "freq", nsets = 20
     )

#obj$celltype_sample = paste(obj$celltype_final, obj$orig.ident)
#Idents(obj) = "celltype_region_age"
ct = "Micro"
oligo = av$RNA[genes,grep(ct, colnames(av$RNA))]


options(repr.plot.width=12, repr.plot.height=5)

cts = obj@meta.data
cts$orig.ident = gsub("2mo", "02mo", cts$orig.ident)
cts$orig.ident = gsub("9mo", "09mo", cts$orig.ident)


cts$celltype_region_age = gsub("2mo", "02mo", cts$celltype_region_age)
cts$celltype_region_age = gsub("9mo", "09mo", cts$celltype_region_age)


cts = cts[grep(ct, cts$celltype_final),]
region_levels = names(table(cts$region)[order(table(cts$region), decreasing = T)])
cts$region = factor(cts$region, levels = region_levels)

ggplot(cts, aes(x = orig.ident, fill = age, color = region)) +
  geom_bar(stat = "count") +
  theme_minimal() +
  labs(title = "Cells per cluster",
       x = "Cluster",
       y = "\nCell #")+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+ facet_grid(~region, scales = "free")

upset(fromList(list_genes), order.by = "freq", nsets = 20)

ggplot(cts, aes(x = celltype_region_age, fill = age, color = region)) +
  geom_bar(stat = "count") +
  theme_minimal() +
  labs(title = "Cells number",
       x = "id",
       y = "\nCell #")+ theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=1))+ facet_grid(~region, scales = "free")


oligo = as.matrix(oligo)
#oligo = oligo[which(rowSums(oligo)>10),]
oligo = oligo[,c(grep("2mo", colnames(oligo)),grep("9mo", colnames(oligo)),grep("18mo", colnames(oligo)))]
region_levels[which(region_levels=="CP")]= " CP"
oligo = oligo[,c(grep(region_levels[1], colnames(oligo)),
  grep(region_levels[2], colnames(oligo)),
  grep(region_levels[3], colnames(oligo)),
  grep(region_levels[4], colnames(oligo)),
  grep(region_levels[5], colnames(oligo)), 
  grep(region_levels[6], colnames(oligo)),
  grep(region_levels[7], colnames(oligo)),
  grep(region_levels[8], colnames(oligo)))]

region <- colnames(oligo)
region[grep("AMY", region)] = "AMY"
region[grep("RLP", region)] = "RLP"
region[grep("HCA", region)] = "HCA"
region[grep("FC", region)] = "FC"
region[grep("NAC", region)] = "NAC"
region[grep(" CP", region)] = "CP"
region[grep("HCP", region)] = "HCP"

region[grep("ENT", region)] = "ENT"
age <- colnames(oligo)
age[grep("18mo", age)] = "18mo"
age[grep("2mo", age)] = "2mo"
age[grep("9mo", age)] = "9mo"

gaps =c(
grep(region_levels[1], colnames(oligo))[length(grep(region_levels[1], colnames(oligo)))],
grep(region_levels[2], colnames(oligo))[length(grep(region_levels[2], colnames(oligo)))],
grep(region_levels[3], colnames(oligo))[length(grep(region_levels[3], colnames(oligo)))],
grep(region_levels[4], colnames(oligo))[length(grep(region_levels[4], colnames(oligo)))],
grep(region_levels[5], colnames(oligo))[length(grep(region_levels[5], colnames(oligo)))],
grep(region_levels[6], colnames(oligo))[length(grep(region_levels[6], colnames(oligo)))],
grep(region_levels[7], colnames(oligo))[length(grep(region_levels[7], colnames(oligo)))]
)
annotation_col = data.frame(
    age = age,
    region = region
    
  )
annotation_col$age =factor(annotation_col$age, levels = c("2mo","9mo", "18mo"))

rownames(annotation_col) = colnames(oligo)
oligo = as.matrix(oligo)
options(repr.plot.width=12, repr.plot.height=5)

brks <- seq(-3,3,length.out=40) 


ph = pheatmap(oligo, kmeans_k = 6,show_rownames = F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps, cluster_rows = F, breaks = brks,main = length(genes))


#ph = pheatmap(oligo, kmeans_k = 6,show_rownames = F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
#        gaps_col = gaps, breaks = brks,main = length(genes))


pheatmap(oligo, show_rownames = F,cutree_rows = 6,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps, breaks = brks, main = length(genes))


pheatmap(oligo, show_rownames = F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps, breaks = brks)



ph

table(ph$kmeans$cluster)

cl4 = names(ph$kmeans$cluster[which(ph$kmeans$cluster==4)])

cat(cl4, sep = "','")


dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human")
if (websiteLive) {
    enriched <- enrichr(cl4, dbs)
}


if (websiteLive) head(enriched[["GO_Biological_Process_2023"]],50)


if (websiteLive) head(enriched[["KEGG_2021_Human"]],30)



amy = oligo[, grep("AMY", colnames(oligo))]
amy = t(scale(t(amy)))
rlp = oligo[, grep("RLP", colnames(oligo))]
rlp = t(scale(t(rlp)))
hca = oligo[, grep("HCA", colnames(oligo))]
hca = t(scale(t(hca)))
hcp = oligo[, grep("HCP", colnames(oligo))]
hcp = t(scale(t(hcp)))
fc = oligo[, grep("FC", colnames(oligo))]
fc = t(scale(t(fc)))
ent = oligo[, grep("ENT", colnames(oligo))]
ent = t(scale(t(ent)))
nac = oligo[, grep("NAC", colnames(oligo))]
nac = t(scale(t(nac)))
cp = oligo[, grep(" CP", colnames(oligo))]
cp = t(scale(t(cp)))
o = cbind(amy,rlp,nac,cp,fc,ent,hca,hcp)
#o = cbind(cp, amy,nac)

pheatmap(o, kmeans_k = 10,show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F)#,
        #gaps_col = gaps)#, gaps_col = c(8,14,20,27,33,41,48))

pheatmap(o, show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F)#,
       # gaps_col = gaps)#, , breaks = brksgaps_col = c(8,14,20,27,33,41,48))


pheatmap(o[-which(is.na(rowSums(o))),], kmeans_k = 6,show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps)#,
      

pheatmap(o[-which(is.na(rowSums(o))),],show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps)#,
      


