library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(enrichR)
setwd("~/projects/fmouse_multiome/SoloTE_out/")

all_directories <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
data_directories <- all_directories[grep("_family", all_directories)]

#data_directories = data_directories[grep("FC", data_directories)]
data_directories

meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")

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


subfamily = seuratobj

seuratobj = merged_seuratobj
which(is.na(seuratobj$celltype))


merged_seuratobj@meta.data %>% 
  ggplot(aes(color=sample, x=nFeature_RNA, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 100)

seuratobj=JoinLayers(seuratobj)

seuratobj <- NormalizeData(seuratobj)
seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seuratobj)
seuratobj <- ScaleData(seuratobj, features = all.genes)
#seuratobj <- SCTransform(seuratobj)
seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))
seuratobj <- FindNeighbors(seuratobj, dims = 1:45)
seuratobj <- FindClusters(seuratobj, resolution = 1)
seuratobj <- RunUMAP(seuratobj, dims = 1:45)
DimPlot(seuratobj,reduction="umap")

seuratobj

saveRDS(seuratobj, '~/projects/combined_all/female_RNA/SoloTE/all_family.RDS')

#seuratobj = readRDS("~/projects/combined_all/female_RNA/SoloTE/FC_SCT.RDS")

meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")

options(repr.plot.width=9, repr.plot.height=6)

DimPlot(seuratobj,reduction="umap", label = T)

options(repr.plot.width=20, repr.plot.height=9)

DimPlot(seuratobj,reduction="umap",group.by = "celltype",label = T)

options(repr.plot.width=20, repr.plot.height=9)

DimPlot(seuratobj,reduction="umap",group.by = "celltype",label = T)

DimPlot(seuratobj,reduction="umap",group.by = "celltype",label = T)

table(seuratobj@meta.data$sample)

#seuratobj@meta.data$celltype= sapply(strsplit(as.character(seuratobj@meta.data$orig.ident), "[.]"), "[[", 1)
#seuratobj@meta.data$sample= sapply(strsplit(as.character(seuratobj@meta.data$orig.ident), ":"), "[[", 2)
seuratobj@meta.data$region= sapply(strsplit(as.character(seuratobj@meta.data$sample), "_"), "[[", 1)
seuratobj@meta.data$age= sapply(strsplit(as.character(seuratobj@meta.data$sample), "_"), "[[", 2)
seuratobj@meta.data$rep= sapply(strsplit(as.character(seuratobj@meta.data$sample), "_"), "[[", 3)

options(repr.plot.width=8, repr.plot.height=6)

DimPlot(seuratobj,reduction="umap", group.by = "age")
DimPlot(seuratobj,reduction="umap", group.by = "region")
DimPlot(seuratobj,reduction="umap", group.by = "rep")



grep("L1MA5A", rownames(seuratobj))

rownames(seuratobj)[26670]

seuratobj$age_rep = paste(seuratobj$age, seuratobj$rep)

options(repr.plot.width=17, repr.plot.height=6)

FeaturePlot(seuratobj,"SoloTE-L1MA5A",reduction="umap", split.by = "age")

FeaturePlot(seuratobj,"SoloTE-RLTR18B",reduction="umap", split.by = "age",order = TRUE,pt.size = 3)


options(repr.plot.width=17, repr.plot.height=6)

FeaturePlot(seuratobj,"SoloTE-Lx2B",reduction="umap", split.by = "age",order = TRUE,pt.size = 1)



VlnPlot(seuratobj,"SoloTE-Lx2B",group.by = "region", split.by = "age")

FeaturePlot(seuratobj,"SoloTE-chr2-158385845-158385934-RLTR18B:ERVK:LTR-20.2-+",reduction="umap", split.by = "sample")
FeaturePlot(seuratobj,"SoloTE-chr10-62792618-62792761-B1-Mus1:Alu:SINE-7.6-+",reduction="umap", split.by = "sample")
FeaturePlot(seuratobj,"SoloTE-chr11-93378660-93378818-L1MA5A:L1:LINE-22.0--",reduction="umap", split.by = "sample")


options(repr.plot.width=20, repr.plot.height=5)

FeaturePlot(seuratobj,"Banf2",reduction="umap", split.by = "age")
FeaturePlot(seuratobj,"Meg3",reduction="umap", split.by = "age")
FeaturePlot(seuratobj,"Cmss1",reduction="umap", split.by = "age")



FeaturePlot(seuratobj,"SoloTE-L1M3f",reduction="umap", split.by = "age")


#FeaturePlot(seuratobj,"SoloTE-RLTR18B",reduction="umap", split.by = "age")
FeaturePlot(seuratobj,"SoloTE-L1",reduction="umap", split.by = "age")

VlnPlot(seuratobj,"SoloTE-L1",split.by = "age", group.by = "region")

#Rename cells based on the expression of known marker genes for the 2C-like group ("Zscan4c", "Tcstv3")
Idents(seuratobj) = "age"
newmarkers <- FindAllMarkers(seuratobj,only.pos=TRUE,test.use = "LR",min.pct = .1, logfc.threshold = .25)
newmarkers_signif <- newmarkers[which(newmarkers$p_val_adj<=0.05),]

head(newmarkers_signif,5)

allgenes = rownames(newmarkers_signif[which(newmarkers_signif$p_val_adj<0.05),])

allgenes = allgenes[-grep("Solo", allgenes)]
length(allgenes)

newmarkers_signif_te <- newmarkers_signif[grep("SoloTE",newmarkers_signif$gene),]
head(newmarkers_signif_te)

newmarkers_signif_te[order(newmarkers_signif_te$p_val_adj)[1:15],]

up18 = rownames(newmarkers_signif_te[which(newmarkers_signif_te$cluster=="18mo" & newmarkers_signif_te$p_val_adj<0.05),])
length(up18)

up2 = rownames(newmarkers_signif_te[which(newmarkers_signif_te$cluster=="02mo" & newmarkers_signif_te$p_val_adj<0.05),])
length(up2)

up9 = rownames(newmarkers_signif_te[which(newmarkers_signif_te$cluster=="09mo" & newmarkers_signif_te$p_val_adj<0.05),])
length(up9)

alldte = rownames(newmarkers_signif_te[which(newmarkers_signif_te$p_val_adj<0.05),])

seuratobj$celltype_age = paste(seuratobj$celltype, seuratobj$age, sep = "_")



seuratobj$celltype_age_region = paste(seuratobj$celltype, seuratobj$age, seuratobj$region,sep = "_")

ae = AverageExpression(seuratobj, group.by = "celltype_age_region")

nrow(ae$RNA)

all =  rownames(newmarkers_signif[which(newmarkers_signif$p_val_adj<0.01),])
length(all)



length(alldte)

cm = ae$RNA[which(rownames(ae$RNA)%in%alldte),]


cm = ae$RNA[grep("SoloTE",rownames(ae$RNA)),]


cmct = colnames(cm)
cmct = gsub("-02mo", "", cmct)
cmct = gsub("-09mo", "", cmct)
cmct = gsub("-18mo", "", cmct)
cmct

seuratobj$celltype_region =paste(seuratobj$celltype, seuratobj$region, sep = "-")

goodct = names(table(seuratobj$celltype)[which(table(seuratobj$celltype)>3000)])
goodct
cm = cm[,which(cmct%in% goodct)]



cmct = colnames(cm)
cmct = gsub("-02mo", "", cmct)
cmct = gsub("-09mo", "", cmct)
cmct = gsub("-18mo", "", cmct)

clade = sapply(strsplit(as.character(cmct), " "), tail, 1)
clade[which(clade == 'Neur')] = "Gaba"
age=rep("02mo", ncol(cm))
age[grep("09", colnames(cm))]="09mo"
age[grep("18", colnames(cm))]="18mo"
cm = cm[,order(clade,age)]
annotation_col = data.frame(
    age = age,
    clade = clade
    
  )

rownames(annotation_col) = colnames(cm)
dim(cm)

rownames(cm)[1:10]
type= rep("gene", nrow(cm))
type[grep("Solo", rownames(cm))] = "TE"
annotation_row = data.frame(
    type = type
  )
rownames(annotation_row) = rownames(cm)
dim(cm)

meta$region_type = "none"
meta$region_type[which(meta$region %in% c("FC","ENT"))]="Cortex"
meta$region_type[which(meta$region %in% c("HCA","HCP"))]="Hippocampus"

table(meta$region_type)
t = table(meta$celltype_final,meta$region_type)/rowSums(table(meta$celltype_final,meta$region_type))
t[1:10,]
cortical= rownames(t[which(t[,1]>0.5),])
hippo= rownames(t[which(t[,2]>0.5),])
other= rownames(t[which(t[,3]>0.5),])

#cortical = gsub("/","-", cortical)
#cortical = gsub(" ","-", cortical)

#hippo = gsub("/","-", hippo)
#hippo = gsub(" ","-", hippo)


celt= sapply(strsplit(colnames(cm), "-"), "[[", 1)

options(repr.plot.width=10, repr.plot.height=10)

pheatmap(cm[,which(clade == "Glut" & celt %in%cortical )], show_rownames = T,cluster_cols = F,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row, gaps_col = c(6,12))


options(repr.plot.width=7, repr.plot.height=10)

pheatmap(cm[,which(clade == "Glut" & celt %in%hippo )], show_rownames = T,cluster_cols = F,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row)


options(repr.plot.width=7, repr.plot.height=10)

pheatmap(cm[,which(clade == "Glut" & celt %in%other )], show_rownames = T,cluster_cols = F,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row)


options(repr.plot.width=7, repr.plot.height=10)

pheatmap(cm[,which(clade == "Gaba" & celt %in%other & celt!="RLP Neur")], show_rownames = T,cluster_cols = F,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row)


pheatmap(cm[,which(clade == "Gaba" & !(celt %in% other ))], show_rownames = T,cluster_cols = F,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row)


options(repr.plot.width=15, repr.plot.height=10)

pheatmap(cm[,], show_rownames = T,cluster_cols = F,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row)



pheatmap(cm[,which(clade == "NN")], show_rownames = T,gaps_col  = c(5,10),cluster_cols = F,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col, annotation_row= annotation_row)


options(repr.plot.width=15, repr.plot.height=10)

pheatmap(cm[,which(clade == "NN")], show_rownames = T,cluster_cols = F,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row)


options(repr.plot.width=15, repr.plot.height=10)

pheatmap(cm, show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "locus TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row)


options(repr.plot.width=10, repr.plot.height=5)

pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "locus TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col, annotation_row= annotation_row)


table(sapply(strsplit(as.character(up9), "-"), "[[", 5))[order(table(sapply(strsplit(as.character(up9), "-"), "[[", 5)), decreasing = T)[1:20]]

table(sapply(strsplit(as.character(up18), "-"), "[[", 5))[order(table(sapply(strsplit(as.character(up18), "-"), "[[", 5)), decreasing = T)[1:20]]

table(sapply(strsplit(as.character(up2), "-"), "[[", 5))[order(table(sapply(strsplit(as.character(up2), "-"), "[[", 5)), decreasing = T)[1:20]]

gm = ae$SCT[which(rownames(ae$SCT)%in%allgenes),]
cmct = colnames(gm)
cmct = gsub("-02mo", "", cmct)
cmct = gsub("-09mo", "", cmct)
cmct = gsub("-18mo", "", cmct)

gm = gm[,which(cmct%in% goodct)]

pheatmap(gm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "DEGs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col=annotation_col)


ae = AggregateExpression(seuratobj, group.by = "celltype_sample", normalization.method = "RC", scale.factor = 1e6)


seuratobj$celltype_sample = paste(seuratobj$celltype, seuratobj$sample, sep = "_")
ae = AggregateExpression(seuratobj, group.by = "celltype_sample", normalization.method = "CLR", scale.factor = 1e6)
cm = ae$RNA[which(rownames(ae$RNA)%in%alldte),]
cmct = colnames(cm)
cmct = gsub("-02mo", "", cmct)
cmct = gsub("-09mo", "", cmct)
cmct = gsub("-18mo", "", cmct)
cmct = gsub("-FC-2", "", cmct)
cmct = gsub("-FC-1", "", cmct)

goodct = names(table(seuratobj$celltype)[which(table(seuratobj$celltype)>1000)])
cm = cm[,which(cmct%in% goodct)]
cm = cm[,grep("Glut", colnames(cm))]
pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "locus TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))


seuratobj

seuratobj$celltype_sample = paste(seuratobj$celltype, seuratobj$sample, sep = "_")
ae = AverageExpression(seuratobj, group.by = "celltype_sample")
cm = ae$RNA[which(rownames(ae$RNA)%in%all[1:2000]),]
cmct = colnames(cm)
cmct = gsub("-02mo", "", cmct)
cmct = gsub("-09mo", "", cmct)
cmct = gsub("-18mo", "", cmct)
cmct = gsub("-FC-2", "", cmct)
cmct = gsub("-FC-1", "", cmct)

goodct = names(table(seuratobj$celltype)[which(table(seuratobj$celltype)>1000)])
cm = cm[,which(cmct%in% goodct)]
cm = cm[,grep("Glut", colnames(cm))]
pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "locus TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))


dim(cm)

pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "locus TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))


options(repr.plot.width=7, repr.plot.height=5)

ae = AverageExpression(seuratobj, group.by = "sample")
cm = ae$SCT[which(rownames(ae$SCT)%in%alldte),]

pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "locus TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))




Idents(seuratobj) = "celltype_age"


unique(Idents(seuratobj) )



#Rename cells based on the expression of known marker genes for the 2C-like group ("Zscan4c", "Tcstv3")
Idents(seuratobj) = "celltype_age"
newmarkers <- FindMarkers(seuratobj,`ident.1` = 'L2/3 IT CTX Glut_18mo', `ident.2` = 'L2/3 IT CTX Glut_02mo' )
newmarkers_signif <- newmarkers[which(newmarkers$p_val_adj<=0.05),]

head(newmarkers_signif)

newmarkers_signif_te <- newmarkers[grep("SoloTE",rownames(newmarkers)),]
alldte = rownames(newmarkers_signif_te[which(newmarkers_signif_te$p_val_adj<0.05),])

allgenes =  rownames(newmarkers_signif[which(newmarkers_signif$p_val_adj<0.05),])
allgenes = allgenes[-grep("Solo", allgenes)]
allgenes

table(seuratobj$sample)

options(repr.plot.width=5, repr.plot.height=5)

ae = AverageExpression(seuratobj, group.by = "celltype_sample")
cm = ae$SCT[which(rownames(ae$SCT)%in%allgenes),]
cm=cm[,grep("L2/3 IT CTX Glut", colnames(cm))]
pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "locus TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))


options(repr.plot.width=5, repr.plot.height=5)

ae = AverageExpression(seuratobj, group.by = "celltype_sample")
cm = ae$SCT[which(rownames(ae$SCT)%in%alldte),]
cm=cm[,grep("L2/3 IT CTX Glut", colnames(cm))]
pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "locus TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))


unique(seuratobj$celltype)

ps <- AggregateExpression(seuratobj, assays = "RNA", return.seurat = T, group.by = c("celltype", "sample"))
ps$age = sapply(strsplit(as.character(ps$sample), "-"), "[[", 2)
ps$celltype_age = paste(ps$celltype,ps$age, sep = "_")



ps

counts = ps@assays$RNA$counts

bin.cov = log10(Matrix::rowSums(counts)+1)
bin.cutoff = quantile(bin.cov[bin.cov > 1], 0.99)
idy = which(bin.cov <= bin.cutoff & bin.cov > 1)

length(idy)

hist(bin.cov)

unique(seuratobj$celltype)

setwd("~/projects/combined_all/female_RNA/SoloTE/")

for(ct in unique(seuratobj$celltype)){
#ct = "CA3"
    
cl = gsub("/", "-", ct)
    cl = gsub(" ", "-", cl)

colnames(counts)[grep(ct, colnames(counts))]
counts_L23 = counts[,grep(ct, colnames(counts))]
groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)

options(repr.plot.width=6, repr.plot.height=6)

#hist(rowSums(counts_L23))
min_reads = 100000
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

# Perform DEG analysis
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
}

print("IN")

files = list.files(".", "TE.txt")
files
cur = files[1]
tab = fread(cur)
head(tab)

cur = files[1]
tab = read.table(cur)
head(tab)

table(seuratobj$region[which(seuratobj$celltype=="Astro-NT NN")])
table(seuratobj$region[which(seuratobj$celltype=="Astro-TE NN")])
table(seuratobj$region[which(seuratobj$celltype=="CA3 Glut")])

for(i in 2:length(files)){
    cur = read.table(files[i])
    cur$celltype = 
    tab = rbind(tab, cur)
        
}

options(repr.plot.width=4, repr.plot.height=2)

d = sig[which(sig$PValue<0.05),]
d$type = "gene"
d$type[grep("Solo", rownames(d))] = "TE"
d$dir = "Down"
d$dir[which(d$logFC.groups02mo<0)]= "Up"
ggplot(d, aes(x = type, fill = dir)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip()

options(repr.plot.width=10, repr.plot.height=10)

pheatmap(th[which(type == "TE"),order(groups)], show_rownames = T,scale = "row", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,
        color=colorRampPalette(c("navy", "white", "red"))(40))

pheatmap(th[which(type == "gene"),order(groups)], show_rownames = F,scale = "row", cluster_cols = F,annotation_row= annotation_row,annotation_col=annotation_col,
        color=colorRampPalette(c("navy", "white", "red"))(40))

ps$clade = sapply(strsplit(as.character(ps$celltype), "[ ]"), tail, 1)
ps$clade_age = paste(ps$clade, ps$age)

Idents(ps) <- "clade_age"

deseq_2_18 <- FindMarkers(object = ps, min.cells.group = 2,
                         ident.1 = "Glut 18mo", 
                         ident.2 = "Glut 02mo",
                         test.use = "DESeq2")
head(deseq_2_18, n = 15)

Idents(ps) <- "clade_age"

deseq_2_18_Gaba <- FindMarkers(object = ps, min.cells.group = 2,
                         ident.1 = "Gaba 18mo", 
                         ident.2 = "Gaba 02mo",
                         test.use = "DESeq2")
head(deseq_2_18_Gaba, n = 15)

d = deseq_2_18_Gaba[which(deseq_2_18_Gaba$p_val<0.05),]
d$type = "gene"
d$type[grep("Solo", rownames(d))] = "TE"
d$dir = "Down"
d$dir[which(d$avg_log2FC>0)]= "Up"
ggplot(d, aes(x = type, fill = dir)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip()

d = deseq_2_18_NN[which(deseq_2_18_NN$p_val<0.05),]
d$type = "gene"
d$type[grep("Solo", rownames(d))] = "TE"
d$dir = "Down"
d$dir[which(d$avg_log2FC>0)]= "Up"
ggplot(d, aes(x = type, fill = dir)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip()

d = deseq_2_18[which(deseq_2_18$p_val<0.05),]
d$type = "gene"
d$type[grep("Solo", rownames(d))] = "TE"
d$dir = "Down"
d$dir[which(d$avg_log2FC>0)]= "Up"
ggplot(d, aes(x = type, fill = dir)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip()


deseq_9_18 <- FindMarkers(object = ps, min.cells.group = 2,
                         ident.1 = "Glut 18mo", 
                         ident.2 = "Glut 09mo",
                         test.use = "DESeq2")
head(deseq_9_18, n = 15)

BiocManager::install("MAST")


Idents(seuratobj) = "celltype_age"
deseq_2_18 <- FindMarkers(object = seuratobj, min.pct = .05,logfc.threshold = .2,min.diff.pct = .01,min.cells.group = 10,
                         ident.1 = "L5 IT CTX Glut_02mo", 
                         ident.2 = "L5 IT CTX Glut_18mo",
                         test.use = "LR")
head(deseq_2_18, n = 15)

length(rownames(deseq_2_18[which(deseq_2_18$p_val_adj<0.0000000000000000001),]))

rownames(deseq_2_18)[(which(rownames(deseq_2_18[which(deseq_2_18$p_val_adj<0.01),])%in%rownames(cpm)))]


th  = cpm[rownames(deseq_2_18)[(which(rownames(deseq_2_18[which(deseq_2_18$p_val_adj<0.000001),])%in%rownames(cpm)))],]
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



th  = cpm[grep("Solo", rownames(cpm)),]
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


dim(th)

seuratobj$clade = sapply(strsplit(as.character(seuratobj$celltype), "[ ]"), tail, 1)
seuratobj$clade_age = paste(seuratobj$clade, seuratobj$age)
seuratobj$clade_sample = paste(seuratobj$clade, seuratobj$sample)
ae = AverageExpression(seuratobj, group.by = "clade_sample")


options(repr.plot.width=15, repr.plot.height=5)

ae = AverageExpression(seuratobj, group.by = "clade_sample")
cm = ae$RNA[which(rownames(ae$RNA)%in%des[grep("Solo", des)]),]
cm=cm[,grep("Glut", colnames(cm))]
pheatmap(cm, show_rownames = T,cutree_rows = 4,scale = "row",cluster_cols = T,main = paste(nrow(cm), "subfamily TEs FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))


cm = ae$RNA[which(rownames(ae$RNA)%in%des[-grep("Solo", des)]),]
cm=cm[,grep("Glut", colnames(cm))]
pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "genes FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))


cm = ae$SCT[which(rownames(ae$SCT)%in%des),]
cm=cm[,grep("L2/3 IT CTX Glut", colnames(cm))]
pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), "genes FC"),
         color=colorRampPalette(c("navy", "white", "red"))(40))


rownames(cm)[1:10]
type= rep("gene", nrow(cm))
type[grep("Solo", rownames(cm))] = "TE"
annotation_row = data.frame(
    type = type
  )
rownames(annotation_row) = rownames(cm)


seuratobj

cm = ae$SCT
cm=cm[,grep("L2/3 IT CTX Glut", colnames(cm))]
pheatmap(cm, show_rownames = F,cutree_rows = 4,scale = "row",main = paste(nrow(cm), " FC"),annotation_row=annotation_row,
         color=colorRampPalette(c("navy", "white", "red"))(40))



