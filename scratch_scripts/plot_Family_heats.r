library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(enrichR)
setwd("~/projects/combined_all/ATAC_Object_for_Signac/")

seuratobj = readRDS( '~/projects/combined_all/female_RNA/SoloTE/all_family.RDS')

seuratobj@meta.data$region= sapply(strsplit(as.character(seuratobj@meta.data$sample), "_"), "[[", 1)
seuratobj@meta.data$age= sapply(strsplit(as.character(seuratobj@meta.data$sample), "_"), "[[", 2)
seuratobj@meta.data$rep= sapply(strsplit(as.character(seuratobj@meta.data$sample), "_"), "[[", 3)

seuratobj$celltype_age_region_rep = paste(seuratobj$celltype, seuratobj$age, seuratobj$region,seuratobj$rep,sep = "-")

ae = AverageExpression(seuratobj, group.by = "celltype_age_region")

cm = ae$RNA[grep("SoloTE",rownames(ae$RNA)),]
cmct = colnames(cm)

age=rep("02mo", ncol(cm))
age[grep("09", colnames(cm))]="09mo"
age[grep("18", colnames(cm))]="18mo"
clade = rep(NA, ncol(cm))
clade[grep("Glut", colnames(cm))]= "Glut"
clade[grep("NN", colnames(cm))]= "NN"
clade[grep("Gaba", colnames(cm))]= "Gaba"
clade[grep("IMN", colnames(cm))]= "IMN"

annotation_col = data.frame(
    age = age,
    clade = clade
    
  )

rownames(annotation_col) = colnames(cm)
dim(cm)

gn = names(table(seuratobj$celltype_age_region)[which(table(seuratobj$celltype_age_region)>200)])


##pdf("all_family_region_heatmaps.pdf")
options(repr.plot.width=12, repr.plot.height=11)

for (r in unique(seuratobj$region)){
    fc = cm[,grep(paste("-",r,sep = ""), colnames(cm))]
    fc = fc[,which(colnames(fc)%in%gn)]

    ct= sapply(strsplit(colnames(fc), "-"), "[[", 1)
    fc =fc[,which(ct%in%names(table(ct)[which(table(ct)>2)]))]
    fc=as.data.frame(fc)
    #fc[which(fc>400, arr.ind = T)]= 400
    p = pheatmap(fc, show_rownames = T,cluster_cols = T,scale = "none",main = paste(nrow(cm), "family TEs", r),
             color=colorRampPalette(c("navy", "white", "red"))(99),annotation_col=annotation_col)
    p = pheatmap(fc, show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs", r),
             color=colorRampPalette(c("navy", "white", "red"))(99),annotation_col=annotation_col)

    print(p)
    
}
#dev.off()

while(3>1){dev.off()}

fc = cm[,grep("-ENT", colnames(cm))]
fc = fc[,which(colnames(fc)%in%gn)]

ct= sapply(strsplit(colnames(fc), "-"), "[[", 1)
fc =fc[,which(ct%in%names(table(ct)[which(table(ct)==3)]))]
clade = 

pheatmap(fc, show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs ENT"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col)


fc = cm[,grep("-FC", colnames(cm))]
fc = fc[,which(colnames(fc)%in%gn)]

ct= sapply(strsplit(colnames(fc), "-"), "[[", 1)
fc =fc[,which(ct%in%names(table(ct)[which(table(ct)==3)]))]
pheatmap(fc[,grep("Glut", colnames(fc))], show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col)


fc = cm[,grep("-HCP", colnames(cm))]
fc = fc[,which(colnames(fc)%in%gn)]

ct= sapply(strsplit(colnames(fc), "-"), "[[", 1)
fc =fc[,which(ct%in%names(table(ct)[which(table(ct)==3)]))]
pheatmap(fc[,grep("Glut", colnames(fc))], show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col)


fc = cm[,grep("-HCA", colnames(cm))]
fc = fc[,which(colnames(fc)%in%gn)]
pheatmap(fc[,grep("Glut", colnames(fc))], show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col)


fc = cm[,grep("-HCP", colnames(cm))]
fc = fc[,which(colnames(fc)%in%gn)]
pheatmap(fc[,grep("Glut", colnames(fc))], show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col)


fc = cm[,grep("-AMY", colnames(cm))]
fc = fc[,which(colnames(fc)%in%gn)]
pheatmap(fc[,grep("NN", colnames(fc))], show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col)


fc = cm[,grep("-NAC", colnames(cm))]
fc = fc[,which(colnames(fc)%in%gn)]
pheatmap(fc[,grep("NN", colnames(fc))], show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col)


options(repr.plot.width=10, repr.plot.height=10)


pheatmap(cm[,grep("L2", colnames(cm))], show_rownames = T,cluster_cols = T,scale = "row",main = paste(nrow(cm), "family TEs"),
         color=colorRampPalette(c("navy", "white", "red"))(40),annotation_col=annotation_col)


annotation_col


