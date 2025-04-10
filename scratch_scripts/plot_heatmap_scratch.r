library(ggplot2)
library(data.table)
library(UpSetR)
library(pheatmap)
library(Seurat)
library(enrichR)

mer = readRDS("../../../merfish/filtered.RDS")

options(repr.plot.width=6, repr.plot.height=6)
Idents(mer) = "predicted_id_ext"
DimPlot(mer, label =T)+NoLegend()


rownames(mer)[which(rownames(mer)%in%genes)]

genes

options(repr.plot.width=12, repr.plot.height=4)


FeaturePlot(mer, "Atp10b", split.by = "age")
FeaturePlot(mer, "Rpa3", split.by = "age")
FeaturePlot(mer, "Rps23", split.by = "age")
FeaturePlot(mer, "Rps9", split.by = "age")


options(repr.plot.width=16, repr.plot.height=4)

VlnPlot(mer, "Atp10b", split.by = "age")


websiteLive <- getOption("enrichR.live")
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human", "SynGO_2022")

websiteLive <- getOption("enrichR.live")
dbs <- c( "GO_Biological_Process_2023")

obj = readRDS("../RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))


DimPlot(obj, group.by = "celltype_final", label = T)+NoLegend()


options(repr.plot.width=12, repr.plot.height=4)


FeaturePlot(obj, "Atp1b1", split.by = "age")
FeaturePlot(obj, "Nrf1", split.by = "age")
FeaturePlot(obj, "Ets2", split.by = "age")


setwd("../DEG_results_latent_rep_mito_together/")

files = list.files(".", ".csv", full.names = T)
#files = files[grep("Oligo", files)]
atabs= list()
list_genes = list()
list_genes_up = list()
list_genes_down = list()

for(f in files){
    # Read DAR file
    dar_data <- read.table(f, header=TRUE, sep=',', stringsAsFactors=FALSE)
    dar_data$cl= gsub(".csv","",f)
    dar_data$cl= gsub("./","",dar_data$cl)

    # Filter rows with significant changes
    significant_changes <- subset(dar_data, p_val_adj < 0.05 & abs(avg_log2FC)>0.1 & (`pct.1`>0.01 | `pct.2`>0.01))
    atabs[[f]] = significant_changes
    #list_genes[[dar_data$cl[1]]] = rownames(significant_changes)
    list_genes[[dar_data$cl[1]]] = significant_changes$X
    list_genes_up[[dar_data$cl[1]]] = significant_changes$X[which(significant_changes$avg_log2FC<(-.1))]
    list_genes_down[[dar_data$cl[1]]] = significant_changes$X[which(significant_changes$avg_log2FC>(.1))]

    
}
tabs = do.call(rbind, atabs)
tabs = tabs[which(tabs$p_val_adj<0.05& abs(tabs$avg_log2FC)>0.1 & (tabs$`pct.1`>0.01 | tabs$`pct.2`>0.01)),]
genes = unique(tabs$X)
length(genes)


grep("IOL",files)


head(tabs)

options(repr.plot.width=12, repr.plot.height=12)
upset(fromList(list_genes), order.by = "freq", nsets = 20)

head(tabs)

all_genes <- unlist(list_genes)
genes = names(table(all_genes)[which(table(all_genes)>=4)])

options(repr.plot.width=12, repr.plot.height=12)

upset(fromList(list_genes_up), order.by = "freq", nsets = 20)

upset(fromList(list_genes_down), order.by = "freq", nsets = 20)

obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)


obj$celltype_age = paste(obj$celltype_final, obj$age)
Idents(obj) = "celltype_age"
avca = AverageExpression(obj)


ogenes = read.table("../DEG_results_latent_rep_mito/oxphos_ribo_genes.txt")

ogenes = ogenes$V1

head(oligo)

obj$celltype_region_age_rep = paste(obj$celltype_final, obj$region, obj$age , obj$rep)
Idents(obj) = "celltype_region_age_rep"
av = AverageExpression(obj)


#obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
#Idents(obj) = "celltype_region_age"
ag = AggregateExpression(obj)



cs = colSums(ag$RNA[,grep(ct, colnames(ag$RNA))])
cs[which(cs>10000)]

11,133,532 5,260,362 9,771,326 34,577,097 33,707,208 35,893,337

11M 5M 9M 34M 33M 35M

2061 2134 3753 2492 6914 2911





obj = AddModuleScore(obj, list(ogenes), name = "oxphos_ribo")

obj$oxphos_ribo1[1:10]

u = unique(obj$orig.ident)
u = u[c(grep("2mo", u),grep("9mo", u),grep("18mo", u))]

options(repr.plot.width=20, repr.plot.height=6)
obj$orig.ident = factor(obj$orig.ident , levels = c(u))
Idents(obj) = "age"
VlnPlot(obj, "oxphos_ribo1", split.by = "orig.ident", group.by = "region", pt.size = 0)

head(mer@meta.data$)

options(repr.plot.width=6, repr.plot.height=6)

VlnPlot(mer , "nCount_originalexp",group.by = "batch", pt.size = 0)
VlnPlot(mer , "nCount_originalexp",group.by = "batch", pt.size = 0)

mean(mer$nCount_originalexp)


Idents(mer) = "predicted_id_ext"

avm = AverageExpression(mer)

Idents(obj) = "celltype_final"

av = AverageExpression(obj)

length(which(rownames(av$RNA) %in% rownames(avm$originalexp)))

length(which(colnames(av$RNA) %in% colnames(avm$originalexp)))

gs = rownames(av$SCT)[which(rownames(av$SCT) %in% rownames(avm$originalexp))]
cs = colnames(av$SCT)[which(colnames(av$SCT) %in% colnames(avm$originalexp))]


rna_500 = obj[gs,]
rna_500

length(gs)

rna_500 = subset(obj, features = gs)

VariableFeatures(rna_500) = gs

rna_500 <- SCTransform(rna_500, verbose = FALSE)


gene_meta = read.csv("../../ABC/gene_meta.csv")
rownames(gene_meta) = gene_meta$geneID
head(gene_meta)

gene_meta = gene_meta[-which(duplicated(gene_meta$gene_name)),] 
rownames(gene_meta) = gene_meta$gene_name
gene_meta$length =  gene_meta$end - gene_meta$start
head(gene_meta)

matched_lengths = gene_meta[paste(rownames(ag)),"length"]

# Load necessary libraries
library(biomaRt)

# Use Ensembl to get gene lengths
# Specify a different mirror if necessary
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", mirror = "uswest")

# Get gene lengths
genes <- rownames(counts)
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'transcript_length'),
                   filters = 'external_gene_name',
                   values = genes,
                   mart = ensembl)

# Average transcript length per gene if multiple transcripts exist
gene_lengths <- aggregate(transcript_length ~ external_gene_name, gene_info, mean)
colnames(gene_lengths) <- c("gene", "length")


rownames(gene_lengths) = gene_lengths$gene

Idents(mer) = "predicted_id_ext"

avm = AggregateExpression(mer,return.seurat = T, normalization.method = "CLR")
mcounts = avm@assays$originalexp$data
total_counts <- colSums(mcounts)
mrpkm = mcounts
#mrpkm <- sweep(mcounts, 2, total_counts, FUN = "/") * 1e6


avm$orig.ident$

Idents(obj) = "celltype_final"

ag = AggregateExpression(rna_500,  return.seurat = T, normalization.method = "CLR")

counts = ag@assays$RNA$data
matched_lengths = gene_lengths[paste(rownames(counts)),"length"]
counts = counts[-which(is.na(matched_lengths)),]
matched_lengths = gene_lengths[paste(rownames(counts)),"length"]
total_counts <- colSums(counts)
#rpkm <- sweep(counts, 1, matched_lengths, FUN = "/") * 1e3
rpkm = counts

rpkm = rpkm[which(rownames(rpkm)%in%rownames(mrpkm)),cs]

mrpkm = mrpkm[rownames(rpkm),colnames(rpkm)]


# Create data frames
df1 <- as.data.frame(mrpkm)
df2 <- as.data.frame(rpkm)

# Calculate the correlation between corresponding columns
correlation <- sapply(1:ncol(df1), function(i) cor(df1[, i], df2[, i]))

# Display the correlation
correlation[1:5]

median(correlation)


options(repr.plot.width=6, repr.plot.height=7)

cot <- cbind(celltype = colnames(df1), correlation = correlation)
cot <- as.data.frame(cot)
cot$correlation <- as.numeric(cot$correlation)
cot$celltype <- factor(cot$celltype, levels = cot$celltype[order(cot$correlation)])

# Calculate median correlation
median_cor <- median(cot$correlation)

# Create the bar plot
ggplot(cot, aes(x = celltype, y = correlation)) +
  geom_bar(stat = "identity") +
  #geom_hline(yintercept = median_cor, linetype = "dotted", color = "red") +  # Add dotted line at median
  geom_text(aes(label = sprintf("%.2f", correlation)), vjust = .5, hjust = -.05,size = 3) +  # Add rounded numbers above bars
  coord_flip() +
  xlab("Cell Type") +
  ylab("Correlation") +
  ggtitle("Correlation by Cell Type") +
  theme_classic()

head(df1)
hist(df2[,1])

clade = sapply(strsplit(rownames(cd1), " "), function(x) tail(x, 1))
clade[which(clade == "IMN")]= "NN"
clade[which(clade == "Gaba-Chol")]= "Gaba"
clade[which(clade == "Neur")]= "Gaba"

ann_row = data.frame(clade = clade)
rownames(ann_row) = rownames(cd1)

options(repr.plot.width=9, repr.plot.height=8)
cd1 = cor(df1, df2,  method = "pearson")

hc <- hclust(as.dist(1 - cd1))

# Determine the order of the variables
order <- hc$order

# Reorder the correlation matrix
cd1_ordered <- cd1[order, order]

# Plot the heatmap
pheatmap(cd1_ordered, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = ann_row)


options(repr.plot.width=9, repr.plot.height=8)
cd1 = cor(df1,  method = "pearson")

cd1_ordered <- cd1[order, order]

# Plot the heatmap
pheatmap(cd1_ordered, cluster_rows = F, cluster_cols = F,annotation_row = ann_row)


options(repr.plot.width=9, repr.plot.height=8)
cd1 = cor(df2,  method = "pearson")

cd1_ordered <- cd1[order, order]

# Plot the heatmap
pheatmap(cd1_ordered, cluster_rows = F, cluster_cols = F,annotation_row = ann_row)


Idents(obj) = obj$celltype_final
av = AverageExpression(obj)
rn = av$SCT
rn = rn[VariableFeatures(obj),]
df = as.data.frame(rn)
cor(df)
pheatmap(cor(df))

obj
meta = read.csv("../../../combined_all/final_meta.csv")
head(meta)

head(meta[which(meta$batch=="Female"),])

meta$umap_RNA_x = "NA"
meta$umap_RNA_y = "NA"


umap_coordinates <- Embeddings(obj, reduction = "umap")
head(umap_coordinates)

obj$cell_id = paste("Female:", obj$barcode, "-1",sep = "")


head(which(obj$cell_id%in%meta$cell_id))

umap_coordinates <- Embeddings(obj, reduction = "umap")

umap_coordinates = as.data.frame(umap_coordinates)
umap_coordinates$cell_id = obj$cell_id
#umap_coordinates = umap_coordinates[-which(is.na(umap_coordinates$cell_id)),]
rownames(umap_coordinates) = umap_coordinates$cell_id
nrow((umap_coordinates))

head(umap_coordinates)

meta$umap_RNA_x = umap_coordinates[paste(meta$cell_id), "umap_1"]
meta$umap_RNA_y = umap_coordinates[paste(meta$cell_id), "umap_2"]

head(meta)

write.table(cor(df) , file = "../../Figures/RNA_celltype_correlation_table.txt")

cor(df)["IOL NN",][order(cor(df)["IOL NN",], decreasing = T)]

dim(av$SCT)

rn = av$SCT[gs,cs]

ms = avm$SCT[gs,cs]

avm$SCT["Sirt2",]

ms[1:5,1:5]
rn[1:5,1:5]

ms[1:5,1:5]
rn[1:5,1:5]

ms["Sirt2",]


# Create data frames
df1 <- as.data.frame(ms)
df2 <- as.data.frame(rn)

# Calculate the correlation between corresponding columns
correlation <- sapply(1:ncol(df1), function(i) cor(df1[, i], df2[, i]))

# Display the correlation
correlation[1:5]

hist(correlation)

correlation

cd1 = cor(df1,  method = "pearson")
pheatmap(cd1)

cd1 = cor(df2,  method = "pearson")
pheatmap(cd1)

options(repr.plot.width=8, repr.plot.height=8)
cd1 = cor(df1, df2,  method = "pearson")

pheatmap(cd1, cluster_rows = F , cluster_cols = F)

all(rownames(df1) == rownames(df2))

hc <- hclust(as.dist(1 - cd1))

# Determine the order of the variables
order <- hc$order

# Reorder the correlation matrix
cd1_ordered <- cd1[order, order]

# Plot the heatmap
pheatmap(cd1_ordered, cluster_rows = FALSE, cluster_cols = FALSE)


cor(ms[,1],rn[,1])

mean(mer$nCount_originalexp)

obj$oxphos_ribo1[1:10]

options(repr.plot.width=7, repr.plot.height=5)
t = obj@meta.data[which(obj$celltype_final=="Oligo NN" & obj$region=="RLP"),]
t$orig.ident = factor(t$orig.ident , levels = c("RLP_2mo_1", "RLP_2mo_2", "RLP_9mo_1",  "RLP_9mo_2","RLP_9mo_3", "RLP_18mo_1", "RLP_18mo_2" ))
                      
p <- ggplot(t, aes(x = orig.ident, y = nCount_RNA, fill = orig.ident)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin Plot of nCount_RNA by orig.ident",
       x = "orig.ident",
       y = "nCount_RNA") +
  theme_minimal()
p

ggplot(t, aes(x = orig.ident, y = nCount_SCT, fill = orig.ident)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin Plot of nCount_SCT by orig.ident",
       x = "orig.ident",
       y = "nCount_SCT") +
  theme_minimal()

ggplot(t, aes(x = orig.ident, y = percent.ribo, fill = orig.ident)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin Plot of percent.ribo by orig.ident",
       x = "orig.ident",
       y = "percent.ribo") +
  theme_minimal()

ggplot(t, aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_violin(trim = FALSE) +
  labs(title = "Violin Plot of percent.ribo by orig.ident",
       x = "orig.ident",
       y = "percent.mt") +
  theme_minimal()

obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

Idents(obj) = "age_celltype_region"


head(Idents(obj))

DimPlot(obj)+NoLegend()

options(repr.plot.width=8, repr.plot.height=8)

DimPlot(mer, label=T, group.by = "predicted_id_ext")+NoLegend()
DimPlot(mer, label=T, group.by = "predicted_cca_ext")+NoLegend()

options(repr.plot.width=14, repr.plot.height=5)

FeaturePlot(mer, "Eif5a", split.by = "age")
VlnPlot(mer, "Eif5a", split.by = "age", group.by= "age")

mer
mer$batch_celltype<- gsub(" ", "_", paste(mer$batch, mer$predicted_id_ext, sep = "_"))

Idents(mer) = "batch_celltype"

avm = AverageExpression(mer)

head(mer@meta.data$batch)

avm$originalexp[which(rownames(avm$originalexp)%in% tfs),]

obj$age_celltype<- gsub(" ", "_", paste(obj$age, obj$celltype_final, sep = "_"))

Idents(obj) = "age_celltype"

de_resultsCA3_Glut <- FindMarkers(
            obj,
            ident.1 = "2mo_CA3_Glut",
            ident.2 = "18mo_CA3_Glut",
            test.use = "MAST",
            logfc.threshold = 0.5, 
            min.pct = 0.01, latent.vars = c('rep', 'percent.mt', "region")
          )

head(de_resultsCA3_Glut)

head(de_resultsM)

de_results["Elk3",]





library(RColorBrewer)
library(ggrepel)

sig = out[which(out$p_val_adj<1e-10 & out$family!= "" & abs(out$avg_log2FC)>1),]

head(sig)

options(repr.plot.width=7, repr.plot.height=5)

out = de_results
fam = "Elk"

out$family = ""
out$family[grep("Elk", out$X)]="Elk"
#out$family[which(out$X %in% fam)]=fam
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  fam,""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-1 & out$family!= "" & abs(out$avg_log2FC)>.1),]
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) 

options(repr.plot.width=7, repr.plot.height=5)

out = de_results
fam = "ETS"

out$family = ""
out$family[which(out$X%in% ets)]="ETS"
#out$family[which(out$X %in% fam)]=fam
#out$family[which(out$X %in% unclustered_protos)]="unclustered_protocadherins"
#out$family[which(out$X %in% major_cadherins)]="major_cadherins"

#out$family[grep("Atp", out$X)]="Atp"
#out$family[grep("Cdh", out$X)]="Cdh"

out$family = factor(out$family, levels = c(  fam,""))
out = out[order(out$family,decreasing = T),]
#out = out[order(out$p_val,decreasing = F),]

sig = out[which(out$p_val_adj<1e-1 & out$family!= "" & abs(out$avg_log2FC)>.1),]
#sig = out[which(out$X==" "),]




ggplot(out, aes(x = (out$`pct.1`+out$`pct.2`)/2, y = -avg_log2FC, color = family)) +
  geom_point(alpha=0.7) + scale_color_manual(values=c(brewer.pal(3, "Set1")[1],"lightgrey"))+
  theme_minimal() +geom_hline(yintercept=0,linetype="dashed")+ ylim(-10,10)+
  labs(x = "pct.cells", y = "Log Fold Change", color = "Gene Family") +  
geom_text_repel(data = sig, max.overlaps = 100,
                aes(label = sig$X, x = (sig$`pct.1`+sig$`pct.2`)/2, y = -avg_log2FC), color = "black", size = 3) 

head(de_results)

de_results$X = rownames(de_results)

19,771,326 - 234,577,097

cs = colSums(ag$RNA[,grep(ct, colnames(ag$RNA))])
cs[which(cs>10000)]

16393894 - 

 #dar_data <- read.table("../DEG_results_latent_rep_ribo/CEA-BST_Gaba--AMY.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
 

dar_data <- read.table("../DEG_results_latent_rep_ribo///O_NN--ENT.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
genes = dar_data$X[which(dar_data$p_val_adj<0.05)]

genes = tabs[which(abs(tabs$avg_log2FC)>0.25 & tabs$p_val_adj<0.001&tabs$cl=="STR_D12_Gaba"),"X"]
upgenes = tabs[which((tabs$avg_log2FC)>0.25 & tabs$p_val_adj<0.001&tabs$cl=="STR_D12_Gaba"),"X"]
downgenes = tabs[which((tabs$avg_log2FC)< -0.25 & tabs$p_val_adj<0.001&tabs$cl=="STR_D12_Gaba"),"X"]

ct = "D12"
oligo = av$RNA[genes,grep(ct, colnames(av$RNA))]
oligo = oligo[,c(grep("2mo",colnames(oligo)), grep("9mo",colnames(oligo)),grep("18mo",colnames(oligo)))]
age = rep("2mo", ncol(oligo))
age[grep("9mo", colnames(oligo))]  = "9mo"
age[grep("18mo", colnames(oligo))]  = "18mo"
age = factor(age, levels = c("2mo","9mo", "18mo"))

region = rep("HCA", ncol(oligo))
region[grep("HCP", colnames(oligo))]  = "HCP"
region[grep(" CP", colnames(oligo))]  = " CP"
region[grep("AMY", colnames(oligo))]  = "AMY"
region[grep("RLP", colnames(oligo))]  = "RLP"
region[grep("ENT", colnames(oligo))]  = "ENT"
region[grep("FC", colnames(oligo))]  = "FC"
region[grep("NAC", colnames(oligo))]  = "NAC"


#oligo = oligo[,-grep(paste(bads, collapse = '|'), colnames(oligo))]
annotation_col <- data.frame(age = age, region = region)
options(repr.plot.width=9, repr.plot.height=8)

rownames(annotation_col) = colnames(oligo)
ph =pheatmap(oligo[1:50,], cutree_rows = 2,show_rownames = T,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = "Oligo top Age-DEGs", cluster_cols = F, annotation_col=annotation_col)



ph = pheatmap(oligo, cutree_rows = 5,annotation_col=annotation_col,show_rownames = F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = length(genes), cluster_cols = F)


options(repr.plot.width=8, repr.plot.height=13)



merf = avm$originalexp[which(rownames(avm$originalexp)%in% genes),grep(ct, colnames(avm$originalexp))]
mgenes=  rownames(merf)
dir = rep("Up", length(mgenes))
dir[which(mgenes%in%upgenes)] = "Down"
annotation_row = data.frame(dir = dir)
rownames(annotation_row) = mgenes

pheatmap(merf, cutree_rows = 5,show_rownames = T,annotation_row=annotation_row,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = length(genes), cluster_cols = F)


pheatmap(merf, cutree_rows = 5,show_rownames = T,annotation_row=annotation_row,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = length(genes), cluster_cols = F)


colnames(avm$originalexp) = gsub(" ", "_",colnames(avm$originalexp)) 
colnames(avm$originalexp) = gsub("-", "_",colnames(avm$originalexp))
colnames(avm$originalexp) = gsub("/", "-",colnames(avm$originalexp)) 

colnames(avm$originalexp) = gsub("^g", "",colnames(avm$originalexp)) 

head(colnames(avm$originalexp))

genes <- tabs[which(abs(tabs$avg_log2FC) > 0.25 & tabs$p_val_adj < 0.01 & tabs$cl == ct), "X"]
genes

hist(mer$nCount_originalexp)
VlnPlot(mer$nCount_originalexp)

unique(tabs$cl)

colnames(av$RNA) = gsub(" ", "_" , colnames(av$RNA))

obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age )
Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)


#for (ct in unique(tabs$cl)) {
options(repr.plot.width=12, repr.plot.height=9)

for (ct in "Oligo_NN") {
    
  #de_resultsCA3_Glut$X= rownames(de_resultsCA3_Glut)
  #genes <- de_resultsCA3_Glut[which(abs(de_resultsCA3_Glut$avg_log2FC) > 0.25 & de_resultsCA3_Glut$p_val_adj < 0.01), "X"]
  #upgenes = de_resultsCA3_Glut[which((de_resultsCA3_Glut$avg_log2FC) > 0.25 & de_resultsCA3_Glut$p_val_adj < 0.01), "X"]
  #downgenes = de_resultsCA3_Glut[which((de_resultsCA3_Glut$avg_log2FC) < -0.25 & de_resultsCA3_Glut$p_val_adj < 0.01 ), "X"]
  
  genes <- tabs[which(abs(tabs$avg_log2FC) > 0.25 & tabs$p_val_adj < 0.001 & tabs$cl == ct), "X"]
  upgenes = tabs[which((tabs$avg_log2FC) > 0.25 & tabs$p_val_adj < 0.001 & tabs$cl == ct), "X"]
  downgenes = tabs[which((tabs$avg_log2FC) < -0.25 & tabs$p_val_adj < 0.001 & tabs$cl == ct), "X"]
  
  #genes = unique(c(genes, tfs))
  if (length(genes) < 10) {
    next
  }
  
  oligo <- av$RNA[genes, grep(ct, colnames(av$RNA))]
  oligo = oligo[which(rowSums(oligo)>0),]  

    
  ct = gsub("_NN", "", ct)
  ct = gsub("_Glut", "", ct)
  ct = gsub("_Gaba", "", ct)

  if (length(which(colnames(oligo) %in% bad)) > 0) {
    oligo <- oligo[, -which(colnames(oligo) %in% bad)]
  }
  
  if (ncol(oligo) < 3) {
    next
  }
  
  oligo <- oligo[, c(grep("2mo", colnames(oligo)), grep("9mo", colnames(oligo)), grep("18mo", colnames(oligo)))]
  
  region <- rep("HCA", ncol(oligo))
  region[grep("HCP", colnames(oligo))] <- "HCP"
  region[grep(" CP", colnames(oligo))] <- " CP"
  region[grep("AMY", colnames(oligo))] <- "AMY"
  region[grep("RLP", colnames(oligo))] <- "RLP"
  region[grep("ENT", colnames(oligo))] <- "ENT"
  region[grep("FC", colnames(oligo))] <- "FC"
  region[grep("NAC", colnames(oligo))] <- "NAC"
  
  br <- names(table(region)[which(table(region) < 3)])
  
  if (length(br) > 0) {
    oligo <- oligo[, -which(region %in% br)]
  }
  
  if (ncol(oligo) < 3) {
    next
  }
  
  age <- rep("2mo", ncol(oligo))
  age[grep("9mo", colnames(oligo))] <- "9mo"
  age[grep("18mo", colnames(oligo))] <- "18mo"
  age <- factor(age, levels = c("2mo", "9mo", "18mo"))
  
  region <- rep("HCA", ncol(oligo))
  region[grep("HCP", colnames(oligo))] <- "HCP"
  region[grep(" CP", colnames(oligo))] <- " CP"
  region[grep("AMY", colnames(oligo))] <- "AMY"
  region[grep("RLP", colnames(oligo))] <- "RLP"
  region[grep("ENT", colnames(oligo))] <- "ENT"
  region[grep("FC", colnames(oligo))] <- "FC"
  region[grep("NAC", colnames(oligo))] <- "NAC"
  
  annotation_col <- data.frame(age = age, region = region)
  rownames(annotation_col) <- colnames(oligo)
  
   #ph <- pheatmap(oligo[tfs[which(tfs%in%rownames(oligo))], ], cutree_rows = 2, show_rownames = TRUE, scale = "row",
   #                color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
   #                main = "top TFs", cluster_cols = FALSE, annotation_col = annotation_col)
   # pdf(paste0(ct, "_TFs.pdf"),height = 17,width = 7)
   # print(ph)
   # dev.off()
    
    
  if (nrow(oligo) > 50) {
    ph <- pheatmap(oligo[1:50, ], cutree_rows = 2, show_rownames = TRUE, scale = "row",
                   color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                   main = "top Age-DEGs", cluster_cols = FALSE, annotation_col = annotation_col)
    pdf(paste0(ct, "_top50AgeDEGs.pdf"),height = 7,width = 10)
    print(ph)
    dev.off()
  }
  
  show_rows <- ifelse(nrow(oligo) < 50, TRUE, FALSE)
  
  ph <- pheatmap(oligo, cutree_rows = 7, annotation_col = annotation_col, show_rownames = show_rows,
                 scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                 main = length(genes), cluster_cols = FALSE)
  
  clusters <- cutree(ph$tree_row, k = 7)
  
  cluster_genes <- vector("list", length = max(clusters))
  
  for (i in seq_along(clusters)) {
    cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
  }
  
  cluster_annotation <- data.frame(Cluster = as.character(clusters))
  rownames(cluster_annotation) <- rownames(oligo)
  
  ph_with_annotation <- pheatmap(oligo, cutree_rows = 7, show_rownames = show_rows, scale = "row",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes), ct), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation,annotation_col=annotation_col)
  
  pdf(paste0(ct, "_withAnnotation.pdf"), height = 7,width = 10)
  print(ph_with_annotation)
  dev.off()
  
  for (i in 1:7) {
    enriched <- enrichr(cluster_genes[i][[1]], dbs)
    write(cluster_genes[i][[1]], file = paste0(ct, "_Cluster", i, "_genes.txt"), sep = "\n")
    if (nrow(enriched$GO_Biological_Process_2023) > 0) {
      pdf(paste0(ct, "_Cluster", i, "_enrichment.pdf"), height = 4, width = 6)
      print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("Cluster", i)))
      dev.off()
    }
  }
  
  options(repr.plot.width=5, repr.plot.height=10)
    if(length(which(rownames(avm$originalexp)%in% genes))<2 | length(grep(ct, colnames(avm$originalexp)))<2){ 
        cat(ct, " no merf \n")
        next
    }
        merf = avm$originalexp[which(rownames(avm$originalexp)%in% genes),grep(ct, colnames(avm$originalexp))]
        mgenes = rownames(merf)
        dir = rep("TF", length(mgenes))
        dir[which(mgenes%in%upgenes)] = "Down"
        dir[which(mgenes%in%downgenes)] = "Up"

        annotation_row = data.frame(dir = dir)
        rownames(annotation_row) = mgenes

        ph = pheatmap(merf, cutree_rows = 2,show_rownames = T,annotation_row=annotation_row,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
                 cluster_rows = T,main = nrow(merf), cluster_cols = F)

      pdf(paste0(ct, "_MERFISH.pdf"), height = 9,width = 5)
      print(ph)
      dev.off()
  
    
}


 ph <- pheatmap(oligo, cutree_rows = 7, annotation_col = annotation_col, show_rownames = show_rows,
                 scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                 main = length(genes), cluster_cols = FALSE)
  
  clusters <- cutree(ph$tree_row, k = 7)
  
  cluster_genes <- vector("list", length = max(clusters))
  
  for (i in seq_along(clusters)) {
    cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
  }
  
  cluster_annotation <- data.frame(Cluster = as.character(clusters))
  rownames(cluster_annotation) <- rownames(oligo)
  
  ph_with_annotation <- pheatmap(oligo, cutree_rows = 7, show_rownames = show_rows, scale = "row",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes), ct), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation,annotation_col=annotation_col)




mer$age_celltype<- gsub(" ", "_", paste(mer$age, mer$predicted_id_ext, sep = "_"))

Idents(mer) = "age_celltype"
ct = "Oligo_NN"

de_results <- FindMarkers(
            mer,
            ident.1 = paste("2mo", ct , sep = "_"),
            ident.2 = paste("18mo", ct , sep = "_"),
            test.use = "MAST",
            logfc.threshold = 0.5, 
            min.pct = 0.05
          )

write.table(de_results, paste("MERFISH_DEG_",ct, ".txt", sep = ""), sep = "\t", quote = F)

de_results = de_results[which(de_results$p_val_adj<0.01),]
merf = avm$originalexp[which(rownames(avm$originalexp)%in% genes),grep(ct, colnames(avm$originalexp))]
merf = merf[which(rownames(merf) %in% rownames(de_results)),]
mgenes = rownames(merf)
        dir = rep("TF", length(mgenes))
        dir[which(mgenes%in%upgenes)] = "Down"
        dir[which(mgenes%in%downgenes)] = "Up"

        annotation_row = data.frame(dir = dir)
        rownames(annotation_row) = mgenes

        ph = pheatmap(merf, cutree_rows = 2,show_rownames = T,annotation_row=annotation_row,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
                 cluster_rows = T,main = nrow(merf), cluster_cols = F)



df_wide <- reshape2::dcast(df, CellType ~ Age, value.var = "Expression")

head(df_wide)


head(df_wide)

options(repr.plot.width=10, repr.plot.height=10)

DimPlot(obj,group.by= "celltype_final", label = T)+NoLegend()

options(repr.plot.width=10, repr.plot.height=10)

FeaturePlot(obj,"Sirt2", label = T)+NoLegend()



obj$celltype_final[which(obj$final_clusters_final == "IOL")] = "IOL NN"


obj$celltype_age<- gsub(" ", "_", paste(obj$celltype_final, obj$age, sep = "_"))

Idents(obj) = "celltype_age"

av = AverageExpression(obj)
sct = av$SCT[which(rowSums(av$SCT)>0.25),]
nrow(sct)

head(df)

options(repr.plot.width=5, repr.plot.height=10)
g = "Enpp6"
eif5a_data <- t(as.matrix(av$SCT[g, , drop=FALSE]))

# Create a data frame
df <- data.frame(
  CellType_Age = rownames(eif5a_data),
  Expression = as.numeric(eif5a_data)
)



df$Age <- "2mo"
df$Age[grep("9mo", df$CellType_Age)] <- "9mo"
df$Age[grep("18mo", df$CellType_Age)] <- "18mo"

df$CellType <- gsub("-18mo" , "" ,df$CellType_Age)
df$CellType <- gsub("-9mo" , "" ,df$CellType)
df$CellType <- gsub("-2mo" , "" ,df$CellType)


df$CellType_Age <- NULL

# Reshape data using melt and dcast
df_long <- melt(df, id.vars = c("CellType", "Age"), measure.vars = "Expression")
df_wide <- dcast(df_long, CellType ~ Age, value.var = "value")

# Set row names and plot heatmap
row.names(df_wide) <- df_wide$CellType
df_wide <- df_wide[,-1]
df_wide = df_wide[,c(2,3,1)]
df_wide = df_wide[which(rowSums(df_wide)>0),]
pheatmap(df_wide, cluster_cols = F, main = g,color=colorRampPalette(c("navy", "white", "red"))(40))
pheatmap(df_wide, cluster_cols = F, scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40))


obj

table(obj$celltype_final)

cnt = 0
s2_ex = s2
s2_ex[which(s2>cnt)] = TRUE
s2_ex[which(s2<=cnt)] = FALSE
table(s2_ex)
d = cbind(sample = paste(obj$orig.ident), celltype = paste(obj$celltype_final), age = paste(obj$age), s2_ex , s2)
d = as.data.frame(d)
d$age = factor(d$age, levels = c("2mo", "9mo", "18mo"))
#head(d)

library(dplyr)
d_fraction <- d %>%
  group_by(celltype) %>%
  summarize(fraction_s2_ex = mean(s2 > cnt)) %>%
  arrange(fraction_s2_ex)
d$celltype = factor(d$celltype, levels = rev(d_fraction$celltype))


s2_ex = s2
s2_ex[which(s2>cnt)] = TRUE
s2_ex[which(s2<=cnt)] = FALSE
table(s2_ex)

d_fraction <- d %>%
  group_by(celltype, age) %>%
  summarize(fraction_s2_ex_1 = mean(s2_ex ==1), .groups = 'drop')

# Reorder the celltype factor based on the fraction of cells with s2_ex equal to 1
d_fraction <- d_fraction %>%
  arrange(age, fraction_s2_ex_1)

d_fraction$celltype <- factor(d_fraction$celltype, levels = unique(d_fraction$celltype))

options(repr.plot.width=17, repr.plot.height=8)

ggplot(d_fraction, aes(x = celltype, y = fraction_s2_ex_1, fill = as.factor(age))) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = scales::percent(fraction_s2_ex_1, accuracy = 1), y = fraction_s2_ex_1 + 0.02),
            vjust = 0, size = 3) +
  theme_minimal() +
  labs(title = paste(g,"expression per Cell Type"),
       x = "Cell Type",
       y = paste("Fraction of cells with counts >" , cnt)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~age, ncol = 1) +
  scale_fill_manual(values = c("1" = "blue"), labels = c("1" = "Age Group"))

#g = "Sirt2"
#s2 = obj@assays$SCT$counts[g,]

cnt = 1
s2_ex = s2
s2_ex[which(s2>cnt)] = TRUE
s2_ex[which(s2<=cnt)] = FALSE
table(s2_ex)
d = cbind(sample = paste(obj$orig.ident), celltype = paste(obj$celltype_final), age = paste(obj$age),region=paste(obj$region),s2_ex , s2)
d = as.data.frame(d)
d$age = factor(d$age, levels = c("2mo", "9mo", "18mo"))
#head(d)

library(dplyr)
d_fraction <- d %>%
  group_by(celltype) %>%
  summarize(fraction_s2_ex = mean(s2 > cnt)) %>%
  arrange(fraction_s2_ex)
d$celltype = factor(d$celltype, levels = rev(d_fraction$celltype))


s2_ex = s2
s2_ex[which(s2>cnt)] = TRUE
s2_ex[which(s2<=cnt)] = FALSE
table(s2_ex)

d_fraction <- d %>%
  group_by(celltype, age,region) %>%
  summarize(fraction_s2_ex_1 = mean(s2_ex ==1), .groups = 'drop')

# Reorder the celltype factor based on the fraction of cells with s2_ex equal to 1
d_fraction <- d_fraction %>%
  arrange(age, fraction_s2_ex_1)

d_fraction$celltype <- factor(d_fraction$celltype, levels = unique(d_fraction$celltype))

options(repr.plot.width=27, repr.plot.height=12)

ggplot(d_fraction, aes(x = celltype, y = fraction_s2_ex_1, fill = as.factor(age))) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = scales::percent(fraction_s2_ex_1, accuracy = 1), y = fraction_s2_ex_1 + 0.02),
            vjust = 0, size = 3) +
  theme_minimal() +
  labs(title = paste(g,"expression per Cell Type"),
       x = "Cell Type",
       y = paste("Fraction of cells with counts >" , cnt)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(age~region, ncol = 8) +
  scale_fill_manual(values = c("1" = "blue"), labels = c("1" = "Age Group"))

length(which(obj$final_clusters_final == "IOL"))
all(which(obj$celltype_final == "IOL NN")==which(obj$final_clusters_final == "IOL"))


options(repr.plot.width=24, repr.plot.height=6)

VlnPlot(obj, "Sirt2", group.by="celltype_final", pt.size = 0, slot = "counts",sort = T)+NoLegend()


head(obj$celltype_age)

iol_markers = FindMarkers(obj, `ident.1`  = 'IOL_NN_18mo', `ident.2` = 'IOL_NN_2mo', latent.vars = "region", test.use = "MAST")

h

iol_markers["Sirt2",]

cor(obj@assays$RNA$scale.data["Sirt2",],obj@assays$RNA$scale.data["Lypd1",])



Idents(obj) = "celltype_final"
iol_nn_markers = FindMarkers(obj, `ident.1`  = 'IOL NN' )

sig_IOL_age = iol_markers[which(iol_markers$p_val_adj<0.5),]
sig_IOL = iol_nn_markers[which(iol_nn_markers$p_val_adj<0.1),]

head(sig_IOL_age)

head(sig_IOL)

sig_IOL_age[which(rownames(sig_IOL_age) %in% rownames(sig_IOL)),]

length(unique(tabs$X))

obj$celltype_final[which(obj$final_clusters_final=="IOL")] = "IOL NN"
table(obj$celltype_final)

head(obj$celltype_age)


obj$celltype_age<- gsub(" ", "--", paste(obj$celltype_final, obj$age, sep = "-"))

Idents(obj) = "celltype_age"

av = AverageExpression(obj,)


head(paste(names(table(obj$celltype_age)[which(table(obj$celltype_age)>150)])))
head(colnames(av$RNA))

avt = av$RNA[,paste(names(table(obj$celltype_age)[which(table(obj$celltype_age)>150)]))]

options(repr.plot.width=7, repr.plot.height=7)

genes = unique(tabs$X)
oligo = avt[genes,]
oligo = oligo[,c(grep("2mo",colnames(oligo)), grep("9mo",colnames(oligo)),grep("18mo",colnames(oligo)))]
age = rep("2mo", ncol(oligo))
age[grep("9mo", colnames(oligo))]  = "9mo"
age[grep("18mo", colnames(oligo))]  = "18mo"
age = factor(age, levels = c("2mo","9mo", "18mo"))

region = rep("HCA", ncol(oligo))
region[grep("HCP", colnames(oligo))]  = "HCP"
region[grep(" CP", colnames(oligo))]  = " CP"
region[grep("AMY", colnames(oligo))]  = "AMY"
region[grep("RLP", colnames(oligo))]  = "RLP"
region[grep("ENT", colnames(oligo))]  = "ENT"
region[grep("FC", colnames(oligo))]  = "FC"
region[grep("NAC", colnames(oligo))]  = "NAC"


#oligo = oligo[,-grep(paste(bads, collapse = '|'), colnames(oligo))]
annotation_col <- data.frame(age = age, region = region)

rownames(annotation_col) = colnames(oligo)
ph =pheatmap(oligo[1:50,], cutree_rows = 2,show_rownames = F, show_colnames=F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = "Oligo top Age-DEGs", cluster_cols = F, annotation_col=annotation_col)



ph =pheatmap(oligo, cutree_rows = 10,show_rownames = F, show_colnames=F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = "Oligo top Age-DEGs", cluster_cols = F, annotation_col=annotation_col)

clusters <- cutree(ph$tree_row, k = 10)

# Create a list to store genes for each cluster
cluster_genes <- vector("list", length = max(clusters))

# Assign genes to clusters
for (i in seq_along(clusters)) {
  cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
}

# Convert cluster assignment to data frame for annotation
cluster_annotation <- data.frame(Cluster = as.character(clusters))
rownames(cluster_annotation) = rownames(oligo)
# Label the pheatmap with cluster assignments
ph_with_annotation <- pheatmap(oligo, cutree_rows = 10, show_rownames = F, 
                               scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40),
                               cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = FALSE,
                               annotation_row = cluster_annotation, show_colnames=F)


#ph = pheatmap(oligo, cutree_rows = 5,annotation_col=annotation_col,show_rownames = F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
#         cluster_rows = T,main = length(genes), cluster_cols = F)



grep("Sirt2", cluster_genes)
cluster_genes[[8]]

clust8 = cluster_genes[[8]]
oligo8 = oligo[paste(clust8),]

clust4= cluster_genes[[4]]
oligo4 = oligo[paste(clust4),]

options(repr.plot.width=25, repr.plot.height=17)

ph <- pheatmap(oligo8, cutree_rows = 10, show_rownames = F, 
                               scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40),
                               cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = T,
                               annotation_row = cluster_annotation, show_colnames=T)


clusters <- cutree(ph$tree_row, k = 10)

# Create a list to store genes for each cluster
cluster_genes <- vector("list", length = max(clusters))

# Assign genes to clusters
for (i in seq_along(clusters)) {
  cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], rownames(oligo8)[i])
}

# Convert cluster assignment to data frame for annotation
cluster_annotation <- data.frame(Cluster = as.character(clusters))
rownames(cluster_annotation) = rownames(oligo8)
# Label the pheatmap with cluster assignments
ph_with_annotation <- pheatmap(oligo8, cutree_rows = 10, show_rownames = F, 
                               scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40),
                               cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = T,
                               annotation_row = cluster_annotation, show_colnames=F)


oligo88 = oligo8[cluster_genes[[8]],]

options(repr.plot.width=25, repr.plot.height=12)

ph <- pheatmap(oligo88, cutree_rows = 10, show_rownames = T, 
                               scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40),
                               cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = T,
                               annotation_row = cluster_annotation, show_colnames=T)


clusters <- cutree(ph$tree_row, k = 10)

# Create a list to store genes for each cluster
cluster_genes <- vector("list", length = max(clusters))

# Assign genes to clusters
for (i in seq_along(clusters)) {
  cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], rownames(oligo88)[i])
}

# Convert cluster assignment to data frame for annotation
cluster_annotation <- data.frame(Cluster = as.character(clusters))
rownames(cluster_annotation) = rownames(oligo88)
# Label the pheatmap with cluster assignments
ph_with_annotation <- pheatmap(oligo88, cutree_rows = 10, show_rownames = T, 
                               scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40),
                               cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = T,
                               annotation_row = cluster_annotation, show_colnames=F)


cat(rownames(oligo88), sep = "\n")

options(repr.plot.width=15, repr.plot.height=12)

ph <- pheatmap(oligo88[,c(grep("IOL", colnames(oligo88)),grep("OPC", colnames(oligo88)),grep("Oligo", colnames(oligo88)))], cutree_rows = 10, show_rownames = T, 
                               scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40),
                               cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = T,
                               annotation_row = cluster_annotation, show_colnames=T)


length(clust8)

grep("Sirt2", cluster_genes)
#length(cluster_genes[[8]])

ph =pheatmap(oligo, cutree_rows = 2,show_rownames = F, show_colnames=F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = "Oligo top Age-DEGs", cluster_cols = F, annotation_col=annotation_col)


head(iol_nn_markers)

iol_markers

iol_markers["Sirt2",]


mer$celltype_age<-  paste(gsub(" ", "-", mer$predicted_id_ext), mer$age, sep = " ")

Idents(mer) = "celltype_age"

avm = AverageExpression(mer)
sct = avm$SCT[which(rowSums(avm$SCT)>0.25),]
nrow(sct)



colnames(avm$originalexp ) = gsub("-", "_", colnames(avm$originalexp ))

head(colnames(avm$originalexp ))
head(df_long)
head(df_wide)


options(repr.plot.width=5, repr.plot.height=7)
g="Bcas1"
eif5a_data <- t(as.matrix(avm$originalexp[g, , drop=FALSE]))

# Create a data frame
df <- data.frame(
  CellType_Age = rownames(eif5a_data),
  Expression = as.numeric(eif5a_data)
)

# Split CellType_Age into Age and CellType
split_celltype_age <- function(celltype_age) {
  parts <- unlist(strsplit(celltype_age, " "))
  celltype <- parts[1]
  age <- paste(parts[-1], collapse = " ")
  return(c(age, celltype))
}

split_data <- t(sapply(df$CellType_Age, split_celltype_age))
df$Age <- split_data[, 1]
df$CellType <- split_data[, 2]
df$CellType_Age <- NULL

# Reshape data using melt and dcast
df_long <- melt(df, id.vars = c("CellType", "Age"), measure.vars = "Expression")
df_wide <- dcast(df_long, CellType ~ Age, value.var = "value")

# Set row names and plot heatmap
row.names(df_wide) <- df_wide$CellType
df_wide <- df_wide[,-1]
df_wide = df_wide[,c(2,3,1)]

pheatmap(df_wide, cluster_cols = F, main = g,color=colorRampPalette(c("navy", "white", "red"))(40))
pheatmap(df_wide, cluster_cols = F, scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40))







de_results = de_results[which(de_results$p_val_adj<0.001),]
merf = avm$originalexp[which(rownames(avm$originalexp)%in% genes),grep(ct, colnames(avm$originalexp))]
merf = merf[which(rownames(merf) %in% rownames(de_results)),]
mgenes = rownames(merf)
        dir = rep("TF", length(mgenes))
        dir[which(mgenes%in%upgenes)] = "Down"
        dir[which(mgenes%in%downgenes)] = "Up"

        annotation_row = data.frame(dir = dir)
        rownames(annotation_row) = mgenes

        ph = pheatmap(merf, cutree_rows = 2,show_rownames = T,annotation_row=annotation_row,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
                 cluster_rows = T,main = nrow(merf), cluster_cols = F)


options(repr.plot.width=5, repr.plot.height=6)

de_results = de_results[which(de_results$p_val_adj<0.0001 &abs(de_results$avg_log2FC)>.65 ),]
merf = avm$originalexp[which(rownames(avm$originalexp)%in% genes),grep(ct, colnames(avm$originalexp))]
merf = merf[which(rownames(merf) %in% rownames(de_results)),]
mgenes = rownames(merf)
        dir = rep("TF", length(mgenes))
        dir[which(mgenes%in%upgenes)] = "Down"
        dir[which(mgenes%in%downgenes)] = "Up"

        annotation_row = data.frame(dir = dir)
        rownames(annotation_row) = mgenes

        ph = pheatmap(merf, cutree_rows = 2,show_rownames = T,annotation_row=annotation_row,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
                 cluster_rows = T,main = nrow(merf), cluster_cols = F)


de_results["Tmem209",]

 ph = pheatmap(merf, cutree_rows = 2,show_rownames = T,annotation_row=annotation_row,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40),
                 cluster_rows = T,main = nrow(merf), cluster_cols = F)


colnames(merf)



merf

genes = tabs[which(abs(tabs$avg_log2FC)>0.25 & tabs$p_val_adj<0.01&tabs$cl=="Oligo_NN"),"X"]
ct = "Oligo"
oligo = av$RNA[genes,grep(ct, colnames(av$RNA))]
oligo = oligo[,c(grep("2mo",colnames(oligo)), grep("9mo",colnames(oligo)),grep("18mo",colnames(oligo)))]
age = rep("2mo", ncol(oligo))
age[grep("9mo", colnames(oligo))]  = "9mo"
age[grep("18mo", colnames(oligo))]  = "18mo"
age = factor(age, levels = c("2mo","9mo", "18mo"))

region = rep("HCA", ncol(oligo))
region[grep("HCP", colnames(oligo))]  = "HCP"
region[grep(" CP", colnames(oligo))]  = " CP"
region[grep("AMY", colnames(oligo))]  = "AMY"
region[grep("RLP", colnames(oligo))]  = "RLP"
region[grep("ENT", colnames(oligo))]  = "ENT"
region[grep("FC", colnames(oligo))]  = "FC"
region[grep("NAC", colnames(oligo))]  = "NAC"


#oligo = oligo[,-grep(paste(bads, collapse = '|'), colnames(oligo))]
annotation_col <- data.frame(age = age, region = region)

rownames(annotation_col) = colnames(oligo)
ph =pheatmap(oligo[1:50,], cutree_rows = 2,show_rownames = T,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = "Oligo top Age-DEGs", cluster_cols = F, annotation_col=annotation_col)



ph = pheatmap(oligo, cutree_rows = 5,annotation_col=annotation_col,show_rownames = F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = length(genes), cluster_cols = F)





setwd("heatmaps_wg/")





unique(tabs$cl)

bad = names(table(obj$celltype_region_age)[which(table(obj$celltype_region_age)<50)])
bad

obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)


colnames(av$SCT) = gsub(" ", "_", colnames(av$SCT))
colnames(av$SCT) = gsub("/", "-", colnames(av$SCT))

obj$celltype_region_age = gsub(" ", "_", obj$celltype_region_age)
obj$celltype_region_age = gsub("/", "-", obj$celltype_region_age)




kr = fread("Oligo_NN_Cluster1_genes_Motifs/knownResults.txt", sep = "\t")
kr$tf = sapply(strsplit(as.character(kr$`Motif Name`), "[(]"), `[`, 1)

capitalize <- function(x) {
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

# Apply the function to each element in the vector
kr$tf <- sapply(kr$tf, capitalize)

head(kr)




ets = kr[grep("ETS", kr$`Motif Name`),"tf"]
ets = unique(ets$tf)

tfs = unique(kr$tf)
tfs = tfs[which(tfs%in%rownames(av$SCT))]

cor(av$SCT["Ubb",] , av$SCT["Jun",])
cor(av$SCT["Ubb",] , av$SCT["Elk3",])
cor(av$SCT["Ubb",] , av$SCT["Cdh8",])
cor(av$SCT["Ubb",] , av$SCT["Atp1b1",])

cor(av$SCT["Eif5a",] , av$SCT["Elk3",])
cor(av$SCT["Eif5a",] , av$SCT["Jun",])
cor(av$SCT["Eif5a",] , av$SCT["Atp1b1",])
cor(av$SCT["Eif5a",] , av$SCT["Ubb",])
cor(av$SCT["Eif5a",] , av$SCT["Elk3",])


cor(sct["Rps8",] , sct["Tubb2a",])
cor(sct["Rps8",] , sct["Elk3",])
cor(sct["Rps8",] , sct["Jun",])
cor(sct["Rps8",] , sct["Ets1",])
cor(sct["Rps8",] , sct["Ets2",])

cor(sct["Rps8",] , sct["Prkn",])

Idents(obj) = "celltype_region_age_rep"
av = AverageExpression(obj)
sct = av$SCT[which(rowSums(av$SCT)>0.25),]
nrow(sct)


obj$oxphos_ribo

cor(obj@assays$SCT@data["Rplp1",], obj@assays$SCT@data["Jun",])

cor(obj@assays$SCT@data["Atp1b1",], obj@assays$SCT@data["Elk3",])

cor(obj@assays$SCT@data["Rplp1",], obj@assays$SCT@data["Rps8",])

cor(obj$oxphos_ribo1, obj@assays$SCT@data["Rps8",])

cor(obj$oxphos_ribo1, obj@assays$SCT@data["Jun",])
cor(obj$oxphos_ribo1, obj@assays$SCT@data["Elk3",])
cor(obj$oxphos_ribo1, obj@assays$SCT@data[tfs,])

tfs

length(which(rowSums(obj@assays$SCT@data)>5))

#correlations <- apply(sct, 1, function(gene_expression) cor(av$SCT["Rplp1",], gene_expression))

correlations <- apply(obj@assays$SCT@data[tfs,], 1, function(gene_expression) cor(obj$oxphos_ribo1, gene_expression))
            
                      
# Sort the correlations in descending order
sorted_correlations <- sort(correlations, decreasing = TRUE)

# Print the sorted correlations
head(sorted_correlations)

plot(sorted_correlations[1:20])

obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_region_age"
ag = AggregateExpression(obj,return.seurat=T)


ag_d12 = ag[,grep("D12", colnames(ag))]
ag_d12 = AddModuleScore(ag_d12, list(ogenes), name = "oxphos_ribo_de")


head(ag_d12@assays$RNA$data)

d12genes = ogenes[which(ogenes%in%rownames(de_results[which(de_results$avg_log2FC>0),]))]


ag_d12 = ag[,grep("Oligo", colnames(ag))]
ag_d12 = AddModuleScore(ag_d12, list(d12genes), name = "oxphos_ribo_ded12")


correlations <- apply(ag_d12@assays$RNA$data[ets,], 1, function(gene_expression) cor(ag_d12$oxphos_ribo_ded121, gene_expression))
            
                      
# Sort the correlations in descending order
sorted_correlations <- sort(correlations, decreasing = TRUE)

# Print the sorted correlations
head(sorted_correlations)

sorted_correlations[1:20]


plot(sorted_correlations)

plot(sorted_correlations)

#correlations <- apply(sct, 1, function(gene_expression) cor(av$SCT["Rplp1",], gene_expression))

correlations <- apply(av$SCT[tfs,], 1, function(gene_expression) cor(ag$oxphos_ribo_de1, gene_expression))
            
                      
# Sort the correlations in descending order
sorted_correlations <- sort(correlations, decreasing = TRUE)

# Print the sorted correlations
head(sorted_correlations)

#correlations <- apply(sct, 1, function(gene_expression) cor(av$SCT["Rplp1",], gene_expression))

correlations <- apply(obj@assays$SCT@data[which(rowSums(obj@assays$SCT@data)>10),], 1, function(gene_expression) cor(obj$oxphos_ribo1, gene_expression))
            
                      
# Sort the correlations in descending order
sorted_correlations <- sort(correlations, decreasing = TRUE)

# Print the sorted correlations
head(sorted_correlations)

head(sorted_correlations)

names(sorted_correlations)[1:200]


sorted_correlations["Ubb"]

system("pwd", intern=T)
unique(tabs$cl)

co = cor(av$SCT["Ubb",] , av$SCT[,])

Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)
sct = av$SCT[which(rowSums(av$SCT)>0.25),]
nrow(sct)

colnames(av$SCT) = gsub(" ", "_", colnames(av$SCT))
colnames(av$RNA) = gsub(" ", "_", colnames(av$RNA))

Idents(obj) = "celltype_region_age"


av$RNA[1:5,1:5]
av$RNA[1:5,1:5]

#for (ct in unique(tabs$cl)) {
options(repr.plot.width=12, repr.plot.height=6)

for (ct in "Microglia_NN") {

  genes <- tabs[which(abs(tabs$avg_log2FC) > 0.25 & tabs$p_val_adj < 0.01 & tabs$cl == ct), "X"]
  genes = rownames(de_resultsM[which(abs(de_resultsM$avg_log2FC) > 0.25 & de_resultsM$p_val_adj < 0.01),])
  genes = unique(c(genes))
  #  genes = unique(c(genes, tfs))
  if (length(genes) < 10) {
    next
  }
  
  oligo <- av$RNA[genes, grep(ct, colnames(av$RNA))]
  oligo = oligo[which(rowSums(oligo)>0),]  

  if (length(which(colnames(oligo) %in% bad)) > 0) {
    oligo <- oligo[, -which(colnames(oligo) %in% bad)]
  }
  
  if (ncol(oligo) < 3) {
    next
  }
  
  oligo <- oligo[, c(grep("2mo", colnames(oligo)), grep("9mo", colnames(oligo)), grep("18mo", colnames(oligo)))]
  
  region <- rep("HCA", ncol(oligo))
  region[grep("HCP", colnames(oligo))] <- "HCP"
  region[grep(" CP", colnames(oligo))] <- " CP"
  region[grep("AMY", colnames(oligo))] <- "AMY"
  region[grep("RLP", colnames(oligo))] <- "RLP"
  region[grep("ENT", colnames(oligo))] <- "ENT"
  region[grep("FC", colnames(oligo))] <- "FC"
  region[grep("NAC", colnames(oligo))] <- "NAC"
  
  br <- names(table(region)[which(table(region) < 3)])
  
  if (length(br) > 0) {
    oligo <- oligo[, -which(region %in% br)]
  }
  
  if (ncol(oligo) < 3) {
    next
  }
  
  age <- rep("2mo", ncol(oligo))
  age[grep("9mo", colnames(oligo))] <- "9mo"
  age[grep("18mo", colnames(oligo))] <- "18mo"
  age <- factor(age, levels = c("2mo", "9mo", "18mo"))
  
  region <- rep("HCA", ncol(oligo))
  region[grep("HCP", colnames(oligo))] <- "HCP"
  region[grep(" CP", colnames(oligo))] <- " CP"
  region[grep("AMY", colnames(oligo))] <- "AMY"
  region[grep("RLP", colnames(oligo))] <- "RLP"
  region[grep("ENT", colnames(oligo))] <- "ENT"
  region[grep("FC", colnames(oligo))] <- "FC"
  region[grep("NAC", colnames(oligo))] <- "NAC"
  
  annotation_col <- data.frame(age = age, region = region)
  rownames(annotation_col) <- colnames(oligo)
  
 #  ph <- pheatmap(oligo[tfs[which(tfs%in%rownames(oligo))], ], cutree_rows = 2, show_rownames = TRUE, scale = "row",
 #                  color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
 #                  main = "top TFs", cluster_cols = FALSE, annotation_col = annotation_col)
 #   pdf(paste0(ct, "_TFs.pdf"),height = 17,width = 7)
 #   print(ph)
 #   dev.off()
    
    
  if (nrow(oligo) > 50) {
    ph <- pheatmap(oligo[1:50, ], cutree_rows = 2, show_rownames = TRUE, scale = "row",
                   color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                   main = "top Age-DEGs", cluster_cols = FALSE, annotation_col = annotation_col)
    pdf(paste0(ct, "_top50AgeDEGs.pdf"),height = 7,width = 10)
    print(ph)
    dev.off()
  }
  
  show_rows <- ifelse(nrow(oligo) < 50, TRUE, FALSE)
  
  ph <- pheatmap(oligo, cutree_rows = 5, annotation_col = annotation_col, show_rownames = show_rows,
                 scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                 main = length(genes), cluster_cols = FALSE)
  
  clusters <- cutree(ph$tree_row, k = 4)
  
  cluster_genes <- vector("list", length = max(clusters))
  
  for (i in seq_along(clusters)) {
    cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
  }
  
  cluster_annotation <- data.frame(Cluster = as.character(clusters))
  rownames(cluster_annotation) <- rownames(oligo)
  
  ph_with_annotation <- pheatmap(oligo, cutree_rows = 4, show_rownames = show_rows, scale = "row",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes), ct), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation,annotation_col=annotation_col)
  
  pdf(paste0(ct, "_withAnnotation.pdf"), height = 7,width = 10)
  print(ph_with_annotation)
  dev.off()
  
  for (i in 1:4) {
    enriched <- enrichr(cluster_genes[i][[1]], dbs)
    write(cluster_genes[i][[1]], file = paste0(ct, "_Cluster", i, "_genes.txt"), sep = "\n")
    if (nrow(enriched$GO_Biological_Process_2023) > 0) {
      pdf(paste0(ct, "_Cluster", i, "_enrichment.pdf"), height = 4, width = 6)
      print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("Cluster", i)))
      dev.off()
    }
  }
}


#for (ct in unique(tabs$cl)) {
options(repr.plot.width=12, repr.plot.height=6)

for (ct in "Oligo_NN") {

  genes <- tabs[which(abs(tabs$avg_log2FC) > 0.25 & tabs$p_val_adj < 0.01 & tabs$cl == ct), "X"]
  genes = unique(c(genes, tfs))
  if (length(genes) < 10) {
    next
  }
  
  oligo <- av$SCT[genes, grep(ct, colnames(av$SCT))]
  oligo = oligo[which(rowSums(oligo)>0),]  

  if (length(which(colnames(oligo) %in% bad)) > 0) {
    oligo <- oligo[, -which(colnames(oligo) %in% bad)]
  }
  
  if (ncol(oligo) < 3) {
    next
  }
  
  oligo <- oligo[, c(grep("2mo", colnames(oligo)), grep("9mo", colnames(oligo)), grep("18mo", colnames(oligo)))]
  
  region <- rep("HCA", ncol(oligo))
  region[grep("HCP", colnames(oligo))] <- "HCP"
  region[grep(" CP", colnames(oligo))] <- " CP"
  region[grep("AMY", colnames(oligo))] <- "AMY"
  region[grep("RLP", colnames(oligo))] <- "RLP"
  region[grep("ENT", colnames(oligo))] <- "ENT"
  region[grep("FC", colnames(oligo))] <- "FC"
  region[grep("NAC", colnames(oligo))] <- "NAC"
  
  br <- names(table(region)[which(table(region) < 3)])
  
  if (length(br) > 0) {
    oligo <- oligo[, -which(region %in% br)]
  }
  
  if (ncol(oligo) < 3) {
    next
  }
  
  age <- rep("2mo", ncol(oligo))
  age[grep("9mo", colnames(oligo))] <- "9mo"
  age[grep("18mo", colnames(oligo))] <- "18mo"
  age <- factor(age, levels = c("2mo", "9mo", "18mo"))
  
  region <- rep("HCA", ncol(oligo))
  region[grep("HCP", colnames(oligo))] <- "HCP"
  region[grep(" CP", colnames(oligo))] <- " CP"
  region[grep("AMY", colnames(oligo))] <- "AMY"
  region[grep("RLP", colnames(oligo))] <- "RLP"
  region[grep("ENT", colnames(oligo))] <- "ENT"
  region[grep("FC", colnames(oligo))] <- "FC"
  region[grep("NAC", colnames(oligo))] <- "NAC"
  
  annotation_col <- data.frame(age = age, region = region)
  rownames(annotation_col) <- colnames(oligo)
  
   ph <- pheatmap(oligo[tfs[which(tfs%in%rownames(oligo))], ], cutree_rows = 2, show_rownames = TRUE, scale = "row",
                   color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                   main = "top TFs", cluster_cols = FALSE, annotation_col = annotation_col)
    pdf(paste0(ct, "_TFs.pdf"),height = 17,width = 7)
    print(ph)
    dev.off()
    
    
  if (nrow(oligo) > 50) {
    ph <- pheatmap(oligo[1:50, ], cutree_rows = 2, show_rownames = TRUE, scale = "row",
                   color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                   main = "top Age-DEGs", cluster_cols = FALSE, annotation_col = annotation_col)
    pdf(paste0(ct, "_top50AgeDEGs.pdf"),height = 7,width = 10)
    print(ph)
    dev.off()
  }
  
  show_rows <- ifelse(nrow(oligo) < 50, TRUE, FALSE)
  
  ph <- pheatmap(oligo, cutree_rows = 5, annotation_col = annotation_col, show_rownames = show_rows,
                 scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                 main = length(genes), cluster_cols = FALSE)
  
  clusters <- cutree(ph$tree_row, k = 6)
  
  cluster_genes <- vector("list", length = max(clusters))
  
  for (i in seq_along(clusters)) {
    cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
  }
  
  cluster_annotation <- data.frame(Cluster = as.character(clusters))
  rownames(cluster_annotation) <- rownames(oligo)
  
  ph_with_annotation <- pheatmap(oligo, cutree_rows = 6, show_rownames = show_rows, scale = "row",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes), ct), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation,annotation_col=annotation_col)
  
  pdf(paste0(ct, "_withAnnotation.pdf"), height = 7,width = 10)
  print(ph_with_annotation)
  dev.off()
  
  for (i in 1:6) {
    enriched <- enrichr(cluster_genes[i][[1]], dbs)
    write(cluster_genes[i][[1]], file = paste0(ct, "_Cluster", i, "_genes.txt"), sep = "\n")
    if (nrow(enriched$GO_Biological_Process_2023) > 0) {
      pdf(paste0(ct, "_Cluster", i, "_enrichment.pdf"), height = 4, width = 6)
      print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("Cluster", i)))
      dev.off()
    }
  }
}


ph_with_annotation <- pheatmap(oligo, cutree_rows = 3, show_rownames = show_rows, scale = "row",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes), ct), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation,annotation_col=annotation_col)
 

pdf(paste0(ct, "_Cluster", i, "_enrichment.pdf"), height = 4, width = 6)
      print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("Cluster", i)))
      dev.off()
      dev.off()

      print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("Cluster", i)))


pdf("plots.pdf")
for(ct in unique(tabs$cl)){
genes = tabs[which(abs(tabs$avg_log2FC)>0.25 & tabs$p_val_adj<0.01&tabs$cl==ct),"X"]
if(length(genes)<10){
    next
}
oligo = av$RNA[genes,grep(ct, colnames(av$RNA))]
if(length(which(colnames(oligo)%in%bad))>0){
    oligo=oligo[,-which(colnames(oligo)%in%bad)]
}
if(ncol(oligo)<3) {
    next
}
oligo = oligo[,c(grep("2mo",colnames(oligo)), grep("9mo",colnames(oligo)),grep("18mo",colnames(oligo)))]

region = rep("HCA", ncol(oligo))
region[grep("HCP", colnames(oligo))]  = "HCP"
region[grep(" CP", colnames(oligo))]  = " CP"
region[grep("AMY", colnames(oligo))]  = "AMY"
region[grep("RLP", colnames(oligo))]  = "RLP"
region[grep("ENT", colnames(oligo))]  = "ENT"
region[grep("FC", colnames(oligo))]  = "FC"
region[grep("NAC", colnames(oligo))]  = "NAC"
br = names(table(region)[which(table(region)<3)])    
if(length(br)>0){
    oligo = oligo[,-which(region%in%br)]
}
if(ncol(oligo)<3) {
    next
}
#pdf(paste(ct, ".pdf", sep = ""), height =7 ,width =7)
age = rep("2mo", ncol(oligo))
age[grep("9mo", colnames(oligo))]  = "9mo"
age[grep("18mo", colnames(oligo))]  = "18mo"
age = factor(age, levels = c("2mo","9mo", "18mo"))

region = rep("HCA", ncol(oligo))
region[grep("HCP", colnames(oligo))]  = "HCP"
region[grep(" CP", colnames(oligo))]  = " CP"
region[grep("AMY", colnames(oligo))]  = "AMY"
region[grep("RLP", colnames(oligo))]  = "RLP"
region[grep("ENT", colnames(oligo))]  = "ENT"
region[grep("FC", colnames(oligo))]  = "FC"
region[grep("NAC", colnames(oligo))]  = "NAC"
#oligo = oligo[,-grep(paste(bads, collapse = '|'), colnames(oligo))]
annotation_col <- data.frame(age = age, region = region)

rownames(annotation_col) = colnames(oligo)
if(nrow(oligo)>50){
ph =pheatmap(oligo[1:50,], cutree_rows = 2,show_rownames = T,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = paste( ct), cluster_cols = F, annotation_col=annotation_col)

print(ph)
}
show_rows = F
if(nrow(oligo)<50){
    show_rows = T
}
ph = pheatmap(oligo, cutree_rows = 5,annotation_col=annotation_col,show_rownames = show_rows,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main =paste(length(genes), ct), cluster_cols = F)


#options(repr.plot.width=7, repr.plot.height=7)

# Assuming 'genes' is your gene list and 'ph' is your pheatmap object

# Get cluster assignments
clusters <- cutree(ph$tree_row, k = 4)

# Create a list to store genes for each cluster
cluster_genes <- vector("list", length = max(clusters))

# Assign genes to clusters
for (i in seq_along(clusters)) {
  cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
}

# Convert cluster assignment to data frame for annotation
cluster_annotation <- data.frame(Cluster = as.character(clusters))
rownames(cluster_annotation) = rownames(oligo)
# Label the pheatmap with cluster assignments
ph_with_annotation <- pheatmap(oligo, cutree_rows = 4, show_rownames = show_rows, 
                               scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40),
                               cluster_rows = TRUE, main = paste(length(genes), ct), cluster_cols = FALSE,
                               annotation_row = cluster_annotation)
print(ph_with_annotation)
#options(repr.plot.width=6, repr.plot.height=4)
for(i in 1:4) {
    enriched <- enrichr(cluster_genes[i][[1]], dbs)
    if(nrow(enriched$GO_Biological_Process_2023)>0){
        print(plotEnrich(enriched$GO_Biological_Process_2023)+ggtitle(paste("Cluster",i)))
    }
}
#dev.off()
}
dev.off()

plotEnrich(enriched$GO_Biological_Process_2023)+ggtitle(paste("Cluster",i))

while(3>1){dev.off()}

options(repr.plot.width=6, repr.plot.height=4)
for(i in 1:4) {
    enriched <- enrichr(cluster_genes[i][[1]], dbs)
    print(plotEnrich(enriched$GO_Biological_Process_2023)+ggtitle(paste("Cluster",i)))
}

bads = gsub("_"," ",names(table(obj$orig.ident)[which(table(obj$orig.ident)<3000)]))
bads

## options(repr.plot.width=7, repr.plot.height=12)
#dar_data <- read.table("../DEG_results_region/DG_Glut--HCA.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
#genes = dar_data$X[which(dar_data$p_val_adj<1e-50)]

#genes = list_genes$`Oligo_NN--AMY`
ct = "DG Glut"
reg = " HCA"
oligo = av$SCT[genes,grep(ct, colnames(av$SCT))]
oligo = oligo[,grep(reg, colnames(oligo))]
oligo = oligo[,c(grep("2mo",colnames(oligo)), grep("9mo",colnames(oligo)),grep("18mo",colnames(oligo)))]
#oligo = oligo[,-grep(paste(bads, collapse = '|'), colnames(oligo))]
ph = pheatmap(oligo, cutree_rows = 5,show_rownames = T,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40),
         cluster_rows = T,main = length(genes), cluster_cols = F)






options(repr.plot.width=7, repr.plot.height=7)


# Assuming 'genes' is your gene list and 'ph' is your pheatmap object

# Get cluster assignments
clusters <- cutree(ph$tree_row, k = 4)

# Create a list to store genes for each cluster
cluster_genes <- vector("list", length = max(clusters))

# Assign genes to clusters
for (i in seq_along(clusters)) {
  cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
}

# Convert cluster assignment to data frame for annotation
cluster_annotation <- data.frame(Cluster = as.character(clusters))
rownames(cluster_annotation) = rownames(oligo)
# Label the pheatmap with cluster assignments
ph_with_annotation <- pheatmap(oligo, cutree_rows = 4, show_rownames = F, 
                               scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40),
                               cluster_rows = TRUE, main = paste(length(genes), reg), cluster_cols = FALSE,
                               annotation_row = cluster_annotation)

# Print genes for each cluster
#for (i in seq_along(cluster_genes)) {
##  cat("Genes in cluster", i, ":", cluster_genes[[i]], "\n")
#}


enriched <- enrichr(cluster_genes[1][[1]], dbs)

#head(enriched$GO_Biological_Process_2023,10)
plotEnrich(enriched$GO_Biological_Process_2023)

head(enriched$SynGO_2022,10)

obj$nCount_RNA

options(repr.plot.width=5, repr.plot.height=4)

Idents(obj)= "rep"
DefaultAssay(obj)="RNA"
VlnPlot(obj, features = "percent.mt", split.by = "age", pt.size = 0)
VlnPlot(obj, features = "percent.ribo", split.by = "age", pt.size = 0)
VlnPlot(obj, features = "nFeature_RNA", split.by = "age", pt.size = 0)

max(obj$percent.ribo)

hist(obj$percent.ribo, breaks = 100)

length(which(obj$percent.ribo>0))

table(obj$orig.ident[which(obj$percent.mt>1 & obj$percent.ribo>1)])

table(obj$orig.ident)


