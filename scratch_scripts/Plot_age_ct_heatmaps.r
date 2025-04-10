library(data.table)
library(ggplot2)
library(Seurat)
library(pheatmap)
library(enrichR)
websiteLive <- getOption("enrichR.live")
dbs <- c( "GO_Biological_Process_2023")

library(zellkonverter)
h = readH5AD("../h5ads_final/celltype_batch_age_PMAT_RPM1.h5ad")

count_table <-h@assays@data$X
colnames(count_table)

obj = readRDS("../female_RNA/RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_age = paste(obj$celltype_final, obj$age )
obj$celltype_age_region = paste(obj$celltype_final, obj$age ,obj$region)
Idents(obj) = "celltype_age_region"
av = AverageExpression(obj)

ct = "DG Glut"

rpm_ct = count_table[,grep(ct, colnames(count_table))]


Oligo_RNA = av$RNA[,grep(ct, colnames(av$RNA))]

good_groups = names(table(Idents(obj))[which(table(Idents(obj))>300)])
Oligo_RNA= Oligo_RNA[,which(colnames(Oligo_RNA)%in%good_groups)]
colnames(Oligo_RNA)



#dar = read.csv(paste("../h5ads_final/combined_diff/diff_csvs/diff_peaks_",gsub(" ", "_", ct),"_2vs18_iter.csv", sep = ""))

dar = read.csv(paste("../h5ads_final/female_diff/diff_csvs/diff_peaks_",gsub(" ", "_", ct),":Female_2vs18_iter.csv", sep = ""))


dar$log2.fold_change. = - dar$log2.fold_change.
rownames(dar) = dar$feature.name

nrow(dar[which(dar$adjusted.p.value<0.05),])

rpm = fread(paste("../h5ads_final/celltype_age_RPM_files/",gsub(" ", ".", ct),"_RPM.txt", sep =""))
rpm = as.data.frame(rpm)
rownames(rpm) = rpm$V1
rpm = rpm[,-c(1)]

deg = read.csv(paste("../female_RNA/DEG_results_latent_rep_mito_together/",gsub(" ", "_", ct),".csv", sep = ""))



rownames(deg) = deg$X
deg$avg_log2FC = -deg$avg_log2FC


head(dar)

meta = read.csv("../ABC/gene_meta.csv")
rownames(meta) = meta$`geneID.1`

o2 = read.csv(paste("filtered.",gsub(" ", "_" , ct),".8wk.abc_score.csv",sep =""))
o18 = read.csv(paste("filtered.",gsub(" ", "_" , ct),".18mo.abc_score.csv",sep =""))

rownames(o2) = o2$X
rownames(o18) = o18$X

merged_df <-merge(o2, o18, by = "X", suffixes = c("_2mo", "_18mo"))
merged_df$gene = sapply(strsplit(as.character(merged_df$X), "-"), `[`, 3)
merged_df$loc = paste(sapply(strsplit(as.character(merged_df$X), "-"), `[`, 1), "-", sapply(strsplit(as.character(merged_df$X), "-"), `[`, 2), sep ="")
merged_df$gene_name = meta[paste(merged_df$gene), "gene_name"]
rownames(merged_df) = merged_df$X
merged_df = merged_df[order(merged_df$abc_score_2mo, decreasing = T),]
merged_df = merged_df[-which(duplicated(merged_df$loc)),]
#merged_df = merged_df[-which(duplicated(merged_df$gene_name)),]

merged_df$deg_logfc = deg[paste(merged_df$gene_name), "avg_log2FC"]

down = deg[which( deg$avg_log2FC< -.25), ]
down = down[1:min(nrow(down),1000),"X"]
up = deg[which( deg$avg_log2FC>.25), ]
up = up[1:min(nrow(up),1000),"X"]


merged_df$DE = "no"
merged_df[which(merged_df$gene_name%in% down), "DE"] = "down"
merged_df[which(merged_df$gene_name%in% up), "DE"] = "up"

merged_df$dar_logfc = dar[paste(merged_df$loc), "log2.fold_change."]

down = dar[which( dar$log2.fold_change.< -.15), ]
down = down[1:min(nrow(down),2000),"feature.name"]
up = dar[which( dar$log2.fold_change.>.15), ]
up = up[1:min(nrow(up),2000),"feature.name"]

merged_df$DAR = "no"
merged_df[which(merged_df$loc%in% down), "DAR"] = "down"
merged_df[which(merged_df$loc%in% up), "DAR"] = "up"


merged_df$abc_diff = merged_df$abc_score_18mo-merged_df$abc_score_2mo
merged_df$rpm_2mo = rpm[paste(merged_df$loc), paste(ct,":2mo", sep = "")]
merged_df$rpm_9mo = rpm[paste(merged_df$loc), paste(ct,":9mo", sep = "")]
merged_df$rpm_18mo = rpm[paste(merged_df$loc), paste(ct,":18mo", sep = "")]
merged_df$rpm_diff = merged_df$rpm_18mo-merged_df$rpm_2mo 
merged_df$log_rpm_diff = sign(merged_df$rpm_diff) * log2(abs(merged_df$rpm_diff)+1e-10)

head(merged_df)
table(merged_df$DAR)
table(merged_df$DE)


nrow(merged_df[which(merged_df$abc_score_2mo>.02),])/length(unique(merged_df[which(merged_df$abc_score_2mo>.02),]$gene_name))

length(unique(merged_df$loc))

hist(table(merged_df$gene_name), main = "peaks per gene")

nrow(merged_df)

tab = merged_df[which(merged_df$DAR%in%c("up", "down", "no") & merged_df$DE%in%c("up", "down", "no")),]

tab = tab[which(tab$abc_score_2mo>.1),]

options(repr.plot.width=5, repr.plot.height=4)
ggplot(tab, aes( x=DAR, fill=DE)) +
geom_bar(position = "fill", alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="DAR", y = "Density")+

theme_classic()

options(repr.plot.width=5, repr.plot.height=4)
ggplot(tab, aes( x=DE, fill=DAR)) +
geom_bar(position = "fill", alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="DEG", y = "Density")+

theme_classic()

plot_heatmap <- function(features, mat, enrich, num_groups = 3) {
    options(repr.plot.width=8, repr.plot.height=8)
    features = unique(features)
    tab = mat[features,]
    tab = tab[which(rowSums(tab)>0),]
    age= rep("2mo", ncol(tab))
    age[grep("9mo", colnames(tab))] = "9mo"
    age[grep("18mo", colnames(tab))] = "18mo"

    col_annot = data.frame(age = age)
    rownames(col_annot) = colnames(tab)
    
    ph <- pheatmap(tab, cutree_rows = num_groups, show_rownames = F,
                 scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                 main = length(features), cluster_cols = T, annotation_col = col_annot)
  
  clusters <- cutree(ph$tree_row, k = num_groups)
  
  cluster_genes <- vector("list", length = max(clusters))
  
  for (i in seq_along(clusters)) {
    cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], features[i])
  }
  
  cluster_annotation <- data.frame(Cluster = as.character(clusters))
  rownames(cluster_annotation) <- rownames(tab)
  options(repr.plot.width=8, repr.plot.height=8)
  ph_with_annotation <- pheatmap(tab, cutree_rows = num_groups, scale = "row",show_rownames = F,
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(features)), cluster_cols = T,
                                 annotation_row = cluster_annotation, annotation_col = col_annot)

    #options(repr.plot.width=8, repr.plot.height=4)
    if(enrich){
    for (i in 1:num_groups) {
        enriched <- enrichr(cluster_genes[i][[1]], dbs)
        #write(cluster_genes[i][[1]], file = paste0(ct, "_Cluster", i, "_genes.txt"), sep = "\n")
        if (nrow(enriched$GO_Biological_Process_2023) > 0) {
          #pdf(paste0(ct, "_Cluster", i, "_enrichment.pdf"), height = 4, width = 6)
          print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("Cluster", i)))
          #dev.off()
        }
      }
    enriched <- enrichr(rownames(tab), dbs)
    print(plotEnrich(enriched$GO_Biological_Process_2023) + ggtitle(paste("All")))
    }
} 
    
    

genes = tab$gene_name[which(tab$DAR%in%c("up"))]
mat = Oligo_RNA
plot_heatmap(genes,mat, TRUE,1)

genes = tab$gene_name[which(tab$DAR%in%c("down"))]
mat = Oligo_RNA
plot_heatmap(genes,mat, TRUE,1)

genes = tab$gene_name[which(tab$DAR%in%c("up") & tab$DE%in%c("up"))]
mat = Oligo_RNA
plot_heatmap(genes,mat, TRUE,1)

genes = tab$gene_name[which(tab$DAR%in%c("down") & tab$DE%in%c("down"))]
mat = Oligo_RNA
plot_heatmap(genes,mat, TRUE,1)

features = tab$loc[which(tab$DE%in%c("up"))]
head(features)
mat = rpm_ct
plot_heatmap(features,mat, F)

features = tab$loc[which(tab$DE%in%c("down"))]
mat = rpm_ct
plot_heatmap(features,mat, F)

dir = "up"
features = tab$loc[which(tab$DE%in%c(dir))]
head(features)

convert_to_bed <- function(feature) {
  # Split the feature string by colon and hyphen
  parts <- unlist(strsplit(feature, '[:-]'))
  data.frame(
    chr = parts[1],
    start = as.numeric(parts[2]),
    end = as.numeric(parts[3]),
    stringsAsFactors = FALSE
  )
}

# Apply the function to all features and combine into a single data frame
bed_df <- do.call(rbind, lapply(features, convert_to_bed))
head(bed_df)
# Save as a BED file
write.table(bed_df, file = paste(gsub(" ", "_", ct),dir,"linked_peaks.bed", sep = "_"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

nrow(bed_df)

features = tab$loc[which(tab$DE%in%c("up"))]
head(features)
mat = rpm_ct
plot_heatmap(features,mat, F)

features = tab$loc[which(tab$DE%in%c("down"))]
mat = rpm_ct
plot_heatmap(features,mat, F)

head(rpm)
head(ugs)

options(repr.plot.width=8, repr.plot.height=8)

ugs = merged_df$loc[which(merged_df$DE%in%c("up") & merged_df$abc_score_2mo>0.04)]
ugs = unique(ugs)
head(rpm[ugs,])
up_RNA= as.data.frame(rpm[ugs,])
up_RNA = up_RNA[which(rowSums(up_RNA)>0.1),]
genes = rownames(up_RNA)
ph <- pheatmap(up_RNA, cutree_rows = 3, show_rownames = F,
                 scale = "row", color = colorRampPalette(c("navy", "white", "red"))(40), cluster_rows = TRUE,
                 main = length(genes), cluster_cols = T)
  
  clusters <- cutree(ph$tree_row, k = 3)
  
  cluster_genes <- vector("list", length = max(clusters))
  
  for (i in seq_along(clusters)) {
    cluster_genes[[clusters[i]]] <- c(cluster_genes[[clusters[i]]], genes[i])
  }
  
  cluster_annotation <- data.frame(Cluster = as.character(clusters))
  rownames(cluster_annotation) <- rownames(up_RNA)
  
  ph_with_annotation <- pheatmap(up_RNA, cutree_rows = 3, scale = "row",show_rownames = F,
                                 color = colorRampPalette(c("navy", "white", "red"))(40),
                                 cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = T,
                                 annotation_row = cluster_annotation)
#pheatmap(up_RNA, cutree_rows = 3, scale = "none",
#                                 color = colorRampPalette(c("navy", "white", "red"))(40),
#                                 cluster_rows = TRUE, main = paste(length(genes)), cluster_cols = T,
#                                 annotation_row = cluster_annotation)

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


up_RNA


