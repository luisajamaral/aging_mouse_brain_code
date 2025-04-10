library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(edgeR)
library(enrichR)

meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")

setwd("~/projects/fmouse_multiome/SoloTE_out/")

all_directories <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
data_directories <- all_directories[grep("subfamily", all_directories)]

reg = "HCP"
data_directories = data_directories[grep(reg, data_directories)]

data_directories

# Initialize an empty Seurat object to store the merged data
setwd("~/projects/fmouse_multiome/SoloTE_out/")
directory = data_directories[1]
setwd(directory)
sample = strsplit(directory, "/", )[[1]][2]
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
merged_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100, project = sample)
merged_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100,project = sample)
merged_seuratobj$sample = gsub("_SoloTE_output","",merged_seuratobj$orig.ident)
merged_seuratobj$barcode = paste(merged_seuratobj$sample, sapply(strsplit(rownames(merged_seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
merged_seuratobj$barcode = gsub("02mo","2mo",merged_seuratobj$barcode)
merged_seuratobj$barcode = gsub("09mo","9mo",merged_seuratobj$barcode)
mat = match(merged_seuratobj$barcode, meta$barcode)
merged_seuratobj$celltype = meta$celltype_final[mat]
merged_seuratobj = merged_seuratobj[,which(!is.na(merged_seuratobj$celltype))]


# Loop through each directory and load the data
for (directory in data_directories[-1]) {
  setwd("~/projects/fmouse_multiome/SoloTE_out/")
  setwd(directory)
  sample = strsplit(directory, "/", )[[1]][2]
  cat(sample, "\n")
  solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
  if(dim(solote_matrix)[2]<25) {
      next
  }
  # Create a temporary Seurat object for the current directory
  temp_seuratobj <- CreateSeuratObject(count = solote_matrix, min.features = 100,project = sample)
  temp_seuratobj$sample = gsub("_SoloTE_output","",temp_seuratobj$orig.ident)
  temp_seuratobj$barcode = paste(temp_seuratobj$sample, sapply(strsplit(rownames(temp_seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
  temp_seuratobj$barcode = gsub("02mo","2mo",temp_seuratobj$barcode)
  temp_seuratobj$barcode = gsub("09mo","9mo",temp_seuratobj$barcode)
  mat = match(temp_seuratobj$barcode, meta$barcode)
  temp_seuratobj$celltype = meta$celltype_final[mat]
  temp_seuratobj = temp_seuratobj[,which(!is.na(temp_seuratobj$celltype))]

  # Merge the current Seurat object with the merged_seuratobj
  merged_seuratobj <- merge(x=merged_seuratobj, y=temp_seuratobj)
}


merged_seuratobj$sample = gsub("_SoloTE_output","",merged_seuratobj$orig.ident)
seuratobj = merged_seuratobj
head(merged_seuratobj@meta.data)

seuratobj=JoinLayers(seuratobj)

meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")
seuratobj$barcode = paste(seuratobj$sample, sapply(strsplit(rownames(seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
seuratobj$barcode = gsub("02mo","2mo",seuratobj$barcode)
seuratobj$barcode = gsub("09mo","9mo",seuratobj$barcode)
mat = match(seuratobj$barcode, meta$barcode)
seuratobj$celltype = meta$celltype_final[mat]

seuratobj = seuratobj[,which(!is.na(seuratobj$celltype))]

rs <- rowSums(seuratobj@assays$RNA$counts)

length(which(rs>4))

seuratobj = seuratobj[which(rs>=4),]

seuratobj

seuratobj <- NormalizeData(seuratobj)
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seuratobj)
seuratobj <- ScaleData(seuratobj, features = all.genes)
seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))
seuratobj <- FindNeighbors(seuratobj, dims = 1:45)
seuratobj <- FindClusters(seuratobj, resolution = 1)
seuratobj <- RunUMAP(seuratobj, dims = 1:45)
DimPlot(seuratobj,reduction="umap")



options(repr.plot.width=15, repr.plot.height=5)

DimPlot(seuratobj,reduction="umap", group.by = "celltype", label = T)

options(repr.plot.width=10, repr.plot.height=5)

DimPlot(seuratobj,reduction="umap", group.by = "age", label = T)

saveRDS(seuratobj, '~/projects/combined_all/female_RNA/SoloTE/HCP_SCT.RDS')

seuratobj = readRDS("~/projects/combined_all/female_RNA/SoloTE/HCA_SCT.RDS")

VariableFeatures(seuratobj) = rownames(seuratobj)

#meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")
seuratobj$barcode = paste(seuratobj$sample, sapply(strsplit(rownames(seuratobj@meta.data), "-"), "[[", 1), sep  = ":")
seuratobj$barcode = gsub("02mo","2mo",seuratobj$barcode)
seuratobj$barcode = gsub("09mo","9mo",seuratobj$barcode)
mat = match(seuratobj$barcode, meta$barcode)
seuratobj$celltype = meta$celltype_final[mat]

seuratobj

which(is.na(seuratobj$celltype))

#seuratobj = seuratobj[,which(!is.na(seuratobj$celltype))]

rownames(seuratobj)[grep("SoloTE",rownames(seuratobj))[1:10]]

table(seuratobj$sample )
seuratobj$age = sapply(strsplit(as.character(seuratobj$sample), "_"), "[[", 2)
table(seuratobj$age )

ps <- AggregateExpression(seuratobj, assays = "RNA", return.seurat = T, group.by = c("celltype", "sample"))
ps$age = sapply(strsplit(as.character(ps$sample), "-"), "[[", 2)
ps$celltype_age = paste(ps$celltype,ps$age, sep = "_")



counts = ps@assays$RNA$counts

dim(counts)

bin.cov = log10(Matrix::rowSums(counts)+1)
bin.cutoff = quantile(bin.cov[bin.cov > 1], 0.99)
idy = which(bin.cov <= bin.cutoff & bin.cov > 1)

hist(bin.cov)

unique(seuratobj$celltype)

setwd("~/projects/combined_all/female_RNA/SoloTE/")

for(ct in unique(seuratobj$celltype)){
cl = gsub("/", "-", ct)
cl = gsub(" ", "-", cl)
cl = paste("HCP_", cl,sep = "")
colnames(counts)[grep(ct, colnames(counts))]
counts_L23 = counts[,grep(ct, colnames(counts))]
groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)

options(repr.plot.width=6, repr.plot.height=6)

#hist(rowSums(counts_L23))
min_reads = 50000
colSums(counts_L23)
good = names(table(regs[which(colSums(counts_L23)>min_reads)]))[which(table(regs[which(colSums(counts_L23)>min_reads)])>5)]
good
if(length(good)<1){
    next
}
counts_L23 = counts_L23[,which(regs %in%good)]
bin.cov = log10(Matrix::rowSums(counts_L23)+1)
#hist(bin.cov)
bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.99)
idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
#length(idy)
counts_L23=counts_L23[idy,]
groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
groups = factor(groups, levels = c("18mo","09mo", "02mo"))
regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)


#if(length(unique(regs))>1){
#    design <- model.matrix(~groups+regs)
#} else {
 design <- model.matrix(~groups)
#}
head(design)
# Create a DGEList object for the subset
d <- DGEList(counts = counts_L23, group = groups)

# Filter genes based on counts in the current subset
cpm = cpm(d)


keep = (rowSums(cpm>1)>=2)
table(keep)
d = d[keep,]
# Normalize the filtered subset
d_normalized <- calcNormFactors(d)
d_disp <- estimateDisp(d_normalized)


#plotMDS(d, gene.selection = "common", labels = groups)
#plotMDS(d, gene.selection = "common", labels = regs)

# Perform ANOVA analysis
fit <- glmFit(d_disp, design)

fit_contrasts <- glmLRT(fit, coef = 2:3)
fit_contrasts$table$FDR <- p.adjust(fit_contrasts$table$PValue, method = "BH")  # Adjust p-values for FDR

fit_contrasts$table = fit_contrasts$table[order(fit_contrasts$table$PValue),]
head(fit_contrasts$table )
sig = fit_contrasts$table[which(fit_contrasts$table$PValue<0.05),]
nrow(sig)
head(sig[grep("Solo", rownames(sig)),])

write.table(sig, paste("../SoloTE/",cl, "_DEG_TE.txt",sep = ""), sep = "\t")

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human", "SynGO_2022")

rn = rownames(sig[which(sig$logFC.groups02mo>0),])
rn = sapply(strsplit(as.character(rn), ";"), tail, 1)

enriched_down <- enrichr(rn, dbs)


rn = rownames(sig[which(sig$logFC.groups02mo<0),])
rn = sapply(strsplit(as.character(rn), ";"), tail, 1)

enriched_up <- enrichr(rn, dbs)

head(enriched_up[["GO_Biological_Process_2023"]],5)
head(enriched_down[["GO_Biological_Process_2023"]],5)

up = enriched_up[["GO_Biological_Process_2023"]]
write.table(up[which(up$P.value<0.05),], paste("../SoloTE/",cl, "_GO_BP_up.txt",sep = ""), sep = "\t")


down = enriched_down[["GO_Biological_Process_2023"]]
write.table(down[which(down$P.value<0.05),], paste("../SoloTE/",cl, "_GO_BP_down.txt",sep = ""), sep = "\t")


up = enriched_up[["SynGO_2022"]]
write.table(up[which(up$P.value<0.05),], paste("../SoloTE/",cl, "_SynGO_up.txt",sep = ""), sep = "\t")


down = enriched_down[["SynGO_2022"]]
write.table(down[which(down$P.value<0.05),], paste("../SoloTE/",cl, "_SynGO_down.txt",sep = ""), sep = "\t")


th  = cpm[rownames(sig),]
type= rep("gene", nrow(th))
type[grep("Solo", rownames(th))] = "TE"
annotation_row = data.frame(
    type = type
  )
rownames(annotation_row) = rownames(th)

annotation_col= data.frame(
    age = groups,
    region = regs
  )
rownames(annotation_col) = colnames(th)

options(repr.plot.width=15, repr.plot.height=12)
g=c()
for(i in 2:length(groups)){
    g = c(g, groups[i]!=groups[i-1])
}

groups = factor(groups, levels = c("02mo","09mo", "18mo"))

tt = t(scale(t(th[which(type == "gene"),order(groups)])))
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(tt, n = 50)
if(nrow(tt)>2){
ph = pheatmap(tt, show_rownames = F,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
  breaks = mat_breaks)
pdf(paste(cl, "_gene_heat.pdf",sep = ""), height =10, width =12)
print(ph)
dev.off()
}
options(repr.plot.width=12, repr.plot.height=12)

tt = t(scale(t(th[which(type == "TE"),order(groups)])))
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(tt, n = 50)
sho=TRUE
if(nrow(tt)>50){
    sho = FALSE
}
if(nrow(tt)>2){
ph = pheatmap(tt, show_rownames = sho,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
  breaks= mat_breaks)
pdf(paste(cl, "_SoloTE_heat.pdf",sep = ""), height =10, width =12)
print(ph)
dev.off()
}
options(repr.plot.width=8, repr.plot.height=7)

markers = fit_contrasts$table
markers$gene = rownames(markers)
markers$type = ifelse(grepl("SoloTE", markers$gene), "SoloTE", "Gene")
volcano_plot <- ggplot(markers, aes(x = -logFC.groups02mo, y = -log10(PValue+1e-1000))) +
  geom_point(aes(color = type), size = 2) +
  scale_color_manual(values = c("Gene" = "black", "SoloTE" = "red")) +
  labs(title = paste(ct, "Volcano Plot\n 2mo vs 18mo"), x = "Average Log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()+
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue")  # Add horizontal dotted line at y = -log10(0.01)


top_significant_points <- markers[order(markers$PValue), ][1:25, ]

# Add text labels to the top 15 significant points
v1 = volcano_plot +
  geom_text_repel(data = top_significant_points, aes(x = -logFC.groups02mo, y = -log10(PValue), label = gene), 
            vjust = -0.5, hjust = 0.5, size = 3, color = "black")+
  geom_point(data = subset(markers, grepl("SoloTE", gene)),
             aes(x = -logFC.groups02mo, y = -log10(PValue + 1e-1000)), color = "red", size = 2)

volcano_plot <- ggplot(markers, aes(x = -logFC.groups09mo, y = -log10(PValue+1e-1000))) +
  geom_point(aes(color = type), size = 2) +
  scale_color_manual(values = c("Gene" = "black", "SoloTE" = "red")) +
  labs(title = paste(ct, "Volcano Plot\n 9mo vs 18mo"), x = "Average Log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()+
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue")  # Add horizontal dotted line at y = -log10(0.01)


top_significant_points <- markers[order(markers$PValue), ][1:25, ]

# Add text labels to the top 15 significant points
 v2 = volcano_plot +
  geom_text_repel(data = top_significant_points, aes(x = -logFC.groups09mo, y = -log10(PValue), label = gene), 
            vjust = -0.5, hjust = 0.5, size = 3, color = "black")+
  geom_point(data = subset(markers, grepl("SoloTE", gene)),
             aes(x = -logFC.groups09mo, y = -log10(PValue + 1e-1000)), color = "red", size = 2)
pdf(paste(cl, "_volcano.pdf",sep = ""), height =6, width =6)
print(v1)
print(v2)
dev.off()


dev.off()
}


rownames(fit_contrasts$table)

pdf(paste(cl, "_volcano.pdf",sep = ""), height =6, width =6)
print(v1)
print(v2)
dev.off()



v1

dev.off()

setwd("~/projects/combined_all/female_RNA/SoloTE")


write.table(sig, paste("../SoloTE/",cl, "_DEG_TE.txt",sep = ""), sep = "\t")

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human", "SynGO_2022")

rn = rownames(sig[which(sig$logFC.groups02mo>0),])
rn = sapply(strsplit(as.character(rn), ";"), tail, 1)

enriched_down <- enrichr(rn, dbs)


rn = rownames(sig[which(sig$logFC.groups02mo<0),])
rn = sapply(strsplit(as.character(rn), ";"), tail, 1)

enriched_up <- enrichr(rn, dbs)

head(enriched_up[["GO_Biological_Process_2023"]],5)
head(enriched_down[["GO_Biological_Process_2023"]],5)

up = enriched_up[["GO_Biological_Process_2023"]]
write.table(up[which(up$P.value<0.05),], paste("../SoloTE/",cl, "_GO_BP_up.txt",sep = ""), sep = "\t")


down = enriched_down[["GO_Biological_Process_2023"]]
write.table(down[which(down$P.value<0.05),], paste("../SoloTE/",cl, "_GO_BP_down.txt",sep = ""), sep = "\t")


up = enriched_up[["SynGO_2022"]]
write.table(up[which(up$P.value<0.05),], paste("../SoloTE/",cl, "_SynGO_up.txt",sep = ""), sep = "\t")


down = enriched_down[["SynGO_2022"]]
write.table(down[which(down$P.value<0.05),], paste("../SoloTE/",cl, "_SynGO_down.txt",sep = ""), sep = "\t")


th  = cpm[rownames(sig),]
type= rep("gene", nrow(th))
type[grep("Solo", rownames(th))] = "TE"
annotation_row = data.frame(
    type = type
  )
rownames(annotation_row) = rownames(th)

annotation_col= data.frame(
    age = groups,
    region = regs
  )
rownames(annotation_col) = colnames(th)

options(repr.plot.width=15, repr.plot.height=12)
g=c()
for(i in 2:length(groups)){
    g = c(g, groups[i]!=groups[i-1])
}

groups = factor(groups, levels = c("02mo","09mo", "18mo"))

tt = t(scale(t(th[which(type == "gene"),order(groups)])))
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(tt, n = 50)

ph = pheatmap(tt, show_rownames = F,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
  breaks = mat_breaks)
pdf(paste(cl, "_gene_heat.pdf",sep = ""), height =10, width =12)
print(ph)
dev.off()
options(repr.plot.width=12, repr.plot.height=12)

tt = t(scale(t(th[which(type == "TE"),order(groups)])))
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(tt, n = 50)

ph = pheatmap(tt, show_rownames = T,gaps_col = which(g),scale = "none", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,color   = colorRampPalette(c("navy", "white", "red"))(length(mat_breaks) - 1),
  breaks= mat_breaks)
pdf(paste(cl, "_SoloTE_heat.pdf",sep = ""), height =10, width =12)
print(ph)
dev.off()
options(repr.plot.width=8, repr.plot.height=7)

markers = fit_contrasts$table
markers$gene = rownames(markers)
markers$type = ifelse(grepl("SoloTE", markers$gene), "SoloTE", "Gene")
volcano_plot <- ggplot(markers, aes(x = -logFC.groups02mo, y = -log10(PValue+1e-1000))) +
  geom_point(aes(color = type), size = 2) +
  scale_color_manual(values = c("Gene" = "black", "SoloTE" = "red")) +
  labs(title = paste(ct, "Volcano Plot\n 2mo vs 18mo"), x = "Average Log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()+
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue")  # Add horizontal dotted line at y = -log10(0.01)


top_significant_points <- markers[order(markers$PValue), ][1:25, ]

# Add text labels to the top 15 significant points
v1 = volcano_plot +
  geom_text_repel(data = top_significant_points, aes(x = -logFC.groups02mo, y = -log10(PValue), label = gene), 
            vjust = -0.5, hjust = 0.5, size = 3, color = "black")+
  geom_point(data = subset(markers, grepl("SoloTE", gene)),
             aes(x = -logFC.groups02mo, y = -log10(PValue + 1e-1000)), color = "red", size = 2)
print(v)

volcano_plot <- ggplot(markers, aes(x = -logFC.groups09mo, y = -log10(PValue+1e-1000))) +
  geom_point(aes(color = type), size = 2) +
  scale_color_manual(values = c("Gene" = "black", "SoloTE" = "red")) +
  labs(title = paste(ct, "Volcano Plot\n 9mo vs 18mo"), x = "Average Log2 Fold Change", y = "-log10(p-value)") +
  theme_minimal()+
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "blue")  # Add horizontal dotted line at y = -log10(0.01)


top_significant_points <- markers[order(markers$PValue), ][1:25, ]

# Add text labels to the top 15 significant points
 v2 = volcano_plot +
  geom_text_repel(data = top_significant_points, aes(x = -logFC.groups09mo, y = -log10(PValue), label = gene), 
            vjust = -0.5, hjust = 0.5, size = 3, color = "black")+
  geom_point(data = subset(markers, grepl("SoloTE", gene)),
             aes(x = -logFC.groups09mo, y = -log10(PValue + 1e-1000)), color = "red", size = 2)
pdf(paste(cl, "_volcano.pdf",sep = ""), height =6, width =6)
print(v1)
print(v2)
dev.off()


dev.off()

seuratobj@meta.data$age= sapply(strsplit(as.character(seuratobj@meta.data$sample), "_"), "[[", 2)
DimPlot(seuratobj,reduction="umap", group.by = "age", label =T)

#Rename cells based on the expression of known marker genes for the 2C-like group ("Zscan4c", "Tcstv3")
Idents(seuratobj) = "age"
newmarkers <- FindAllMarkers(seuratobj,only.pos=TRUE)
newmarkers_signif <- newmarkers[which(newmarkers$p_val_adj<=0.05),]

allgenes = rownames(newmarkers_signif[which(newmarkers_signif$p_val_adj<0.05),])
allgenes = allgenes[-grep("Solo", allgenes)]
length(allgenes)

newmarkers_signif_te <- newmarkers_signif[grep("SoloTE",newmarkers_signif$gene),]
head(newmarkers_signif_te)

up18 = rownames(newmarkers_signif_te[which(newmarkers_signif_te$cluster=="18mo" & newmarkers_signif_te$p_val_adj<0.1),])
length(up18)

up2 = rownames(newmarkers_signif_te[which(newmarkers_signif_te$cluster=="02mo" & newmarkers_signif_te$p_val_adj<0.1),])
length(up2)

up9 = rownames(newmarkers_signif_te[which(newmarkers_signif_te$cluster=="09mo" & newmarkers_signif_te$p_val_adj<0.1),])
length(up9)

alldte = rownames(newmarkers_signif_te[which(newmarkers_signif_te$p_val_adj<0.05),])
length(alldte)

seuratobj$celltype_age = paste(seuratobj$celltype, seuratobj$age, sep = "_")


all = rownames(newmarkers_signif[which(newmarkers_signif$p_val_adj<0.05),])


seuratobj$celltype_sample = paste(seuratobj$celltype, seuratobj$sample, sep = "_")
ae = AverageExpression(seuratobj, group.by = "celltype_sample")
cm = ae$RNA
cmct = colnames(cm)
cmct = gsub("-02mo", "", cmct)
cmct = gsub("-09mo", "", cmct)
cmct = gsub("-18mo", "", cmct)
cmct = gsub("-HCP-2", "", cmct)
cmct = gsub("-HCP-1", "", cmct)
cmct = gsub("-HCP-3", "", cmct)

goodct = names(table(seuratobj$celltype)[which(table(seuratobj$celltype)>1000)])
cm = cm[,which(cmct%in% goodct)]


vf = FindVariableFeatures(seuratobj,nfeatures = 10000)

cmct = colnames(cm)
cmct = gsub("-02mo", "", cmct)
cmct = gsub("-09mo", "", cmct)
cmct = gsub("-18mo", "", cmct)
cmct = gsub("-NAC-2", "", cmct)
cmct = gsub("-NAC-1", "", cmct)
cmct = gsub("-NAC-3", "", cmct)
clade = sapply(strsplit(as.character(cmct), " "), tail, 1)
age=rep("02mo", ncol(cm))
age[grep("09", colnames(cm))]="09mo"
age[grep("18", colnames(cm))]="18mo"
annotation_col = data.frame(
    age = age,
    clade = clade
    
  )

rownames(annotation_col) = colnames(cm)
type= rep("gene", nrow(cm))
type[grep("Solo", rownames(cm))] = "TE"
annotation_row = data.frame(
    type = type
  )
rownames(annotation_row) = rownames(cm)

options(repr.plot.width=15, repr.plot.height=10)

pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "NAC"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col,annotation_row= annotation_row )


options(repr.plot.width=15, repr.plot.height=10)

pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "NAC"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col,annotation_row= annotation_row )


ps <- AggregateExpression(seuratobj, assays = "RNA", return.seurat = T, group.by = c("celltype", "sample"))
ps$age = sapply(strsplit(as.character(ps$sample), "-"), "[[", 2)
ps$celltype_age = paste(ps$celltype,ps$age, sep = "_")



Idents(ps) <- "age"

deseq_2_18 <- FindMarkers(object = ps,only.pos = FALSE, min.cells.group = 2,
                         ident.1 = "18mo", 
                         ident.2 = "02mo",
                         test.use = "DESeq2")
head(deseq_2_18, n = 15)

des = c(rownames(deseq_2_18[which(deseq_2_18$p_val_adj<0.05),]))
length(des)

cm = ae$SCT[which(rownames(ae$SCT)%in%des),]
cm=cm[,grep("STR D12 Gaba", colnames(cm))]
rownames(cm)[1:10]
type= rep("gene", nrow(cm))
type[grep("Solo", rownames(cm))] = "TE"
annotation_row = data.frame(
    type = type
  )
rownames(annotation_row) = rownames(cm)

pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "NAC"),annotation_row=annotation_row,
         color=colorRampPalette(c("navy", "white", "red"))(40))


dim(cm)

length(des)


