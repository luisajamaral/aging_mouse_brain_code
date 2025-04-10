library(data.table)

rpm = fread("../h5ads_final/celltype_age_RPM_files/Oligo.NN_RPM.txt")
rpm = as.data.frame(rpm)
rownames(rpm) = rpm$V1

deg = read.csv("../female_RNA/DEG_results_latent_rep_mito_together/Oligo_NN.csv")

meta = read.csv("gene_meta.csv")

head(meta)
rownames(meta) = meta$geneID

require(biomaRt)
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
genes = meta$geneID
annot <- getBM(
  attributes = c(
    'mgi_symbol',
    'external_gene_name',
    'ensembl_gene_id',
    'entrezgene_id',
    'gene_biotype','transcript_length'),
  filters = 'ensembl_gene_id',
  values = genes,
  mart = ensembl)

annot <- merge(
  x = as.data.frame(genes),
  y =  annot,
  by.y = 'ensembl_gene_id',
  all.x = T,
  by.x = 'genes')

o2 = read.csv("Oligo_NN/Oligo_NN.8wk.abc_score.csv")
o9 = read.csv("Oligo_NN/Oligo_NN.9mo.abc_score.csv")
o18 = read.csv("Oligo_NN/Oligo_NN.18mo.abc_score.csv")


rownames(o2) = o2$X
rownames(o9) = o9$X
rownames(o18) = o18$X


merged_df <-merge(o2, o9, by = "X", suffixes = c("_2mo", "_9mo"))

merged_df <-merge(merged_df, o18, by = "X",suffixes = c("", "_18mo"))
head(merged_df)

hist(o2$abc_score, xlim = c(0,5), breaks = 60000)

length(which(o2$abc_score>1))/nrow(o2)

rownames(deg) = deg$X
head(deg)

deg$avg_log2FC = -deg$avg_log2FC

merged_df$deg_logfc = deg[paste(merged_df$gene_name), "avg_log2FC"]

head(merged_df)

merged_df$gene = sapply(strsplit(as.character(merged_df$X), "-"), `[`, 3)
merged_df$loc = paste(sapply(strsplit(as.character(merged_df$X), "-"), `[`, 1), "-", sapply(strsplit(as.character(merged_df$X), "-"), `[`, 2), sep ="")
merged_df$gene_name = meta[paste(merged_df$gene), "gene_name"]
rownames(merged_df) = merged_df$X
down = deg[which(deg$p_val_adj<0.001 & deg$avg_log2FC< -.25), "X"]
up = deg[which(deg$p_val_adj<0.001 & deg$avg_log2FC> .25), "X"]
merged_df$DE = "no"
merged_df[which(merged_df$gene_name%in% down), "DE"] = "down"
merged_df[which(merged_df$gene_name%in% up), "DE"] = "up"

length(which( up %in% merged_df$gene_name))

merged_df$abc_diff = merged_df$abc_score-merged_df$abc_score_2mo



head(merged_df[order(merged_df$abc_diff),])

merged_df$rpm_2mo = rpm[paste(merged_df$loc), "Oligo NN:2mo"]
merged_df$rpm_9mo = rpm[paste(merged_df$loc), "Oligo NN:9mo"]
merged_df$rpm_18mo = rpm[paste(merged_df$loc), "Oligo NN:18mo"]
merged_df$rpm_diff = merged_df$rpm_18mo-merged_df$rpm_2mo 
merged_df$log_rpm_diff = sign(merged_df$rpm_diff) * log2(abs(merged_df$rpm_diff)+1e-10)
head(merged_df)

head(merged_df)

nrow(merged_df[which(merged_df$abc_score>.02),])/length(unique(merged_df[which(merged_df$abc_score>.02),]$gene_name))

nrow(o2[which(o2$abc_score>.02),])/length(unique(o2[which(o2$abc_score>.02),]$gene_name))

nrow(merged_df[which(merged_df$abc_score>.02),])/length(unique(merged_df[which(merged_df$abc_score>.02),]$gene_name))

length(unique(merged_df$loc))

nrow(merged_df)

hist(table(merged_df$gene_name), main = "peaks per gene")

merged_df$changes_abc = "NO"
merged_df$changes_abc[which(merged_df$log_abc_diff>5)] = "UP"
merged_df$changes_abc[which(merged_df$log_abc_diff< -5)] = "DOWN"


merged_df$changes_RNA = "NO"
merged_df$changes_RNA[which(merged_df$deg_logfc>.5)] = "UP"
merged_df$changes_RNA[which(merged_df$deg_logfc< -.5)] = "DOWN"


library(ggplot2)

head(log2(merged_df$abc_diff))

merged_df$log_abc_diff = sign(merged_df$abc_diff)*log2(abs(merged_df$abc_diff)+1e-10)

options(repr.plot.width=5, repr.plot.height=4)

ggplot(merged_df, aes( x=changes_abc, fill=changes_RNA)) +
geom_bar(position = "fill", alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="changes_abc", y = "Density")+

theme_classic()


ggplot(merged_df, aes(y=DE, x=changes_abc, fill=changes_abc)) +
geom_violin(aes(), alpha=0.5)+
#geom_density(alpha=0.6)+ xlim(c(-50,50))+ 
scale_color_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
scale_fill_manual(values=c( "#E69F00","#999999", "#56B4E9"))+
labs(title="abc diff",x="diff", y = "Density")+
geom_hline(aes(yintercept=0),
           linetype="dashed")+
theme_classic()

rpm$diff_18_9=rpm$`Oligo NN:18mo`-rpm$`Oligo NN:9mo`
rpm$diff_9_2=rpm$`Oligo NN:9mo`-rpm$`Oligo NN:2mo`
rpm$diff=rpm$`Oligo NN:18mo`-rpm$`Oligo NN:2mo`

options(repr.plot.width=5, repr.plot.height=4)

hist(rpm$diff, breaks = 150, xlim = c(-15,15), ylim = c(0,70000))
hist(rpm$diff_18_9, breaks = 100, xlim = c(-15,15), ylim = c(0,70000))
hist(rpm$diff_9_2, breaks = 100, xlim = c(-15,15), ylim = c(0,70000))


urpm = drpm

options(repr.plot.width=5, repr.plot.height=4)

hist(urpm$diff, breaks = 50, xlim = c(-15,15))
hist(drpm$diff, breaks = 50, xlim = c(-15,15))
hist(rpm$diff, breaks = 150, xlim = c(-15,15), ylim = c(0,70000))
hist(rpm$diff, breaks = 150, xlim = c(-15,15), ylim = c(0,20000))


o2 = read.csv("Oligo_NN/Oligo_NN.18mo.abc_score.csv")


head(rpm)

o2$gene = sapply(strsplit(as.character(o2$X), "-"), `[`, 3)
o2$loc = paste(sapply(strsplit(as.character(o2$X), "-"), `[`, 1), "-", sapply(strsplit(as.character(o2$X), "-"), `[`, 2), sep ="")
o2$gene_name = meta[paste(o2$gene), "gene_name"]
rownames(o2) = o2$X
head(o2)
length(unique(o2$gene_name))


options(repr.plot.width=5, repr.plot.height=7)

down = deg[which(deg$p_val_adj<0.00001 & deg$avg_log2FC< -.25), "X"]
d2 = o2[which(o2$gene_name%in% down), ]
d2 = d2[order(d2$abc_score, decreasing = T),]
d2 = d2[-which(duplicated(d2$gene_name)),]
d2 = d2[which(d2$abc_score>5),]
drpm = rpm[paste(d2$loc),]
drpm$gene = d2$gene_name
rownames(drpm) = drpm$gene
head(drpm)
pheatmap(drpm[,2:4], show_rownames = F, scale = "row")


down = deg[which(deg$p_val_adj<0.00001 & deg$avg_log2FC>.25), "X"]
d2 = o2[which(o2$gene_name%in% down), ]



up = deg[which(deg$p_val_adj<0.01 & deg$avg_log2FC< -.25), "X"]

d18 = o18[which(o18$gene_name%in% down), ]

u18 = o18[which(o18$gene_name%in% up), ]

tail(d2)

options(repr.plot.width=5, repr.plot.height=7)

down = deg[which(deg$p_val_adj<0.00001 & deg$avg_log2FC> .25), "X"]
d2 = o2[which(o2$gene_name%in% down), ]
d2 = d2[order(d2$abc_score, decreasing = T),]
d2 = d2[-which(duplicated(d2$gene_name)),]
d2 = d2[which(d2$abc_score>5),]
drpm = rpm[paste(d2$loc),]
drpm$gene = d2$gene_name
rownames(drpm) = drpm$gene
head(drpm)
pheatmap(drpm[,2:4], show_rownames = F, scale = "row")

hist(rpm$`Oligo NN:18mo`)
hist(rpm$`Oligo NN:2mo`)

u2 = o2[which(o2$gene_name%in% up), ]
u2 = u2[order(u2$abc_score, decreasing = T),]
u2 = u2[-which(duplicated(u2$gene_name)),]
head(u2)



u18 = u18[order(u18$abc_score, decreasing = T),]
u18 = u18[-duplicated(u18$gene_name),]
head(u18)

nrow(o18)

nrow(o2)

urpm = rpm[paste(u2$loc),]



head(urpm)

head(drpm)

nrow(drpm)

tail(urpm)

library(pheatmap)

pheatmap(urpm[,2:4], show_rownames = F, scale = "row")

pheatmap(drpm[,2:4], show_rownames = F, scale = "row")


