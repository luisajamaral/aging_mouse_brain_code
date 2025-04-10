library(data.table)
library(ggplot2)
library(UpSetR)
library(pheatmap)
library(Seurat)
library(enrichR)
websiteLive <- getOption("enrichR.live")
dbs <- c( "GO_Biological_Process_2023")

obj = readRDS("../../female_RNA/RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))


obj$celltype_age = paste(obj$celltype_final, obj$age )
Idents(obj) = "celltype_age"
av = AverageExpression(obj)

dg_RNA = av$SCT[,grep("DG Glut", colnames(av$SCT))]

rpm = fread("../../h5ads_final/celltype_age_RPM_files/DG.Glut_RPM.txt")
rpm = as.data.frame(rpm)
rownames(rpm) = rpm$V1

deg = read.csv("../../female_RNA/DEG_results_latent_rep_mito_together/DG_Glut.csv")

dar = read.csv("../../h5ads_final/combined_diff/diff_csvs/diff_peaks_DG_Glut_2vs18_iter.csv")

dar$log2.fold_change. = - dar$log2.fold_change.
rownames(dar) = dar$feature.name

rownames(deg) = deg$X
deg$avg_log2FC = -deg$avg_log2FC
head(deg)

head(deg)

deg$chr = gmeta[paste(deg$X), "chrom"]
deg$start = gmeta[paste(deg$X), "start"]
deg$end = gmeta[paste(deg$X), "end"]
deg$length = deg$end-deg$start
head(deg)

options(repr.plot.width=5, repr.plot.height=4)
ggplot(deg, aes( x=dir, y=log(length))) +
geom_boxplot( alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
#scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
#scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="DAR", y = "Density")+

theme_classic()

meta = read.csv("../gene_meta.csv")
rownames(meta) = meta$geneID

head(meta)
gmeta= meta[which(!duplicated(meta$gene_name)),]
rownames(gmeta) = gmeta$gene_name
head(gmeta)




dg_2mo = fread("DG_Glut.8wk.abc_score.csv", sep = ",")
dg_2mo = dg_2mo[which(dg_2mo$abc_score>.04),]

dg_2mo = fread("DG_Glut.18mo.abc_score.csv", sep = ",")
dg_2mo = dg_2mo[which(dg_2mo$abc_score>.04),]

nrow(dg_2mo)

rownames(dg_2mo) = dg_2mo$V1


merged_df = dg_2mo
merged_df$gene = sapply(strsplit(as.character(merged_df$V1), "-"), `[`, 3)
merged_df$loc = paste(sapply(strsplit(as.character(merged_df$V1), "-"), `[`, 1), "-", sapply(strsplit(as.character(merged_df$V1), "-"), `[`, 2), sep ="")
merged_df$gene_name = meta[paste(merged_df$gene), "gene_name"]
rownames(merged_df) = merged_df$V1
merged_df= merged_df[order(merged_df$abc_score, decreasing = T),]


hist(merged_df$activity)
head(merged_df[order(merged_df$activity),])

hist(table(merged_df$gene_name), main = "peaks per gene")

mean(table(merged_df$gene_name))

nrow(dg_2mo)

hist(dg_2mo$abc_score, breaks = 100)


deg$dir = "NS"
deg$dir[which( deg$avg_log2FC>0.15)[1:1000]] = "Up"
deg$dir[which( deg$avg_log2FC< -0.15)[1:1000]] = "Down"
table(deg$dir)
merged_df$deg_logfc = deg[paste(merged_df$gene_name), "avg_log2FC"]
merged_df$deg_sig = deg[paste(merged_df$gene_name), "dir"]


dar$dir = "NS"
dar$dir[which(dar$adjusted.p.value<0.001 & dar$log2.fold_change.>0.5)] = "Up"
dar$dir[which(dar$adjusted.p.value<0.001 & dar$log2.fold_change.< -0.5)] = "Down"
table(dar$dir)

merged_df$dar_logfc = dar[paste(merged_df$loc), "log2.fold_change."]
merged_df$dar_sig = dar[paste(merged_df$loc), "dir"]




tab = merged_df[which(merged_df$dar_sig%in%c("Up", "Down", "NS") & merged_df$deg_sig%in%c("Up", "Down", "NS")),]

nrow(tab)

loc_tab = merged_df[which(!duplicated(merged_df$loc)),]
loc_tab = loc_tab[which(loc_tab$dar_sig%in%c("Up", "Down", "NS")),]
loc_tab = loc_tab[which(!is.na(loc_tab$deg_sig)),]



options(repr.plot.width=5, repr.plot.height=4)
ggplot(tab, aes( x=dar_sig, fill=deg_sig)) +
geom_bar(position = "fill", alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="DAR", y = "Density")+

theme_classic()

options(repr.plot.width=5, repr.plot.height=4)
ggplot(tab, aes( x=deg_sig, fill=dar_sig)) +
geom_bar(position = "fill", alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="DAR", y = "Density")+

theme_classic()

loc_tab = merged_df[which(!duplicated(merged_df$gene_name)),]
loc_tab = loc_tab[which(loc_tab$dar_sig%in%c("Up", "Down", "NS")),]
loc_tab = loc_tab[which(!is.na(loc_tab$deg_sig)),]
options(repr.plot.width=5, repr.plot.height=4)
ggplot(loc_tab, aes( x=dar_sig, fill=deg_sig)) +
geom_bar(position = "fill", alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="DAR", y = "Density")+

theme_classic()

options(repr.plot.width=5, repr.plot.height=4)
ggplot(loc_tab, aes( x=dar_sig, y=deg_logfc)) +
geom_violin( alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
#scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
#scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="DAR", y = "Density")+

theme_classic()

head(rpm)



ugs = loc_tab$gene_name[which(loc_tab$dar_sig%in%c("Up"))]
head(dg_RNA[ugs,])
up_RNA= dg_RNA[ugs,]
up_RNA = up_RNA[which(rowSums(up_RNA)>0),]

up_RNA

options(repr.plot.width=5, repr.plot.height=8)
ugs = tab$gene_name[which(tab$dar_sig%in%c("Down"))]
ugs = unique(ugs)

head(dg_RNA[ugs,])
up_RNA= dg_RNA[ugs,]
up_RNA = up_RNA[which(rowSums(up_RNA)>0),]
pheatmap(up_RNA, scale= "row", cluster_cols = F, show_rownames = F,cutree_rows = 5)

ugs = tab$gene_name[which(tab$dar_sig%in%c("Up"))]
ugs = unique(ugs)
head(dg_RNA[ugs,])
up_RNA= dg_RNA[ugs,]
up_RNA = up_RNA[which(rowSums(up_RNA)>0),]
ph =pheatmap(up_RNA, scale= "row", cluster_cols = F, show_rownames = F,cutree_rows = 5)

options(repr.plot.width=5, repr.plot.height=8)

ugs = tab$gene_name[which(tab$dar_sig%in%c("Down"))]
ugs = unique(ugs)
head(dg_RNA[ugs,])
up_RNA= dg_RNA[ugs,]
up_RNA = up_RNA[which(rowSums(up_RNA)>0.1),]
ph =pheatmap(up_RNA, scale= "row", cluster_cols = F, show_rownames = F,cutree_rows = 3)
genes = rownames(up_RNA)
ph <- pheatmap(up_RNA, cutree_rows = 3, 
                 scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                 main = length(genes), cluster_cols = FALSE)
  
  clusters <- cutree(ph$tree_row, k = 3)
  
  cluster_genes <- vector("list", length = max(clusters))
  
  for (i in seq_along(clusters)) {
    cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
  }
  
  cluster_annotation <- data.frame(Cluster = as.character(clusters))
  rownames(cluster_annotation) <- rownames(up_RNA)
  
  ph_with_annotation <- pheatmap(up_RNA, cutree_rows = 3, scale = "row",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation)
pheatmap(up_RNA, cutree_rows = 3, scale = "none",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation)

options(repr.plot.width=12, repr.plot.height=6)

for (i in 1:3) {
    enriched <- enrichr(cluster_genes[i][[1]], dbs)
    #write(cluster_genes[i][[1]], file = paste0(ct, "_Cluster", i, "_genes.txt"), sep = "\n")
    if (nrow(enriched$GO_Biological_Process_2023) > 0) {
      #pdf(paste0(ct, "_Cluster", i, "_enrichment.pdf"), height = 4, width = 6)
      print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("Cluster", i)))
      #dev.off()
    }
  }



rpm = rpm[,-c(1)]
head(rpm)

head(ugs)
head(up_RNA)

options(repr.plot.width=5, repr.plot.height=8)

ugs = tab$loc[which( tab$deg_sig%in%c("Down"))]
ugs = unique(ugs)
head(rpm[ugs,])
up_RNA= rpm[ugs,]
up_RNA = up_RNA[which(rowSums(up_RNA)>0.1),]
ph =pheatmap(up_RNA, scale= "row", cluster_cols = F, show_rownames = F,cutree_rows = 3)
genes = rownames(up_RNA)
ph <- pheatmap(up_RNA, cutree_rows = 3, 
                 scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                 main = length(genes), cluster_cols = FALSE, show_rownames = F)
  
  clusters <- cutree(ph$tree_row, k = 3)
  
  cluster_genes <- vector("list", length = max(clusters))
  
  for (i in seq_along(clusters)) {
    cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
  }
  
  cluster_annotation <- data.frame(Cluster = as.character(clusters))
  rownames(cluster_annotation) <- rownames(up_RNA)
  
  ph_with_annotation <- pheatmap(up_RNA, cutree_rows = 3, scale = "row",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation, show_rownames=F)
pheatmap(up_RNA, cutree_rows = 3, scale = "none",
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = FALSE,
                                 annotation_row = cluster_annotation, show_rownames=F)




cluster_genes[[1]][1:5]
chromosomes <- c()
starts <- c()
ends <- c()

for (region in cluster_genes[[1]]) {
  parts <- unlist(strsplit(region, "[:-]"))
  chromosomes <- c(chromosomes, parts[1])
  starts <- c(starts, as.numeric(parts[2]))
  ends <- c(ends, as.numeric(parts[3]))
}

# Create a data frame in BED format (0-based start, 1-based end)
bed_df <- data.frame(
  chromosome = chromosomes,
  start = starts - 1, # BED format uses 0-based start positions
  end = ends
)

head(bed_df)

nrow(bed_df)
write.table(bed_df, file = "DG_peaks_linked_down_genes_down.bed", quote=F, row.names=F, col.names=F, sep = "\t")

meta
