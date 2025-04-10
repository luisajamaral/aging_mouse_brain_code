library(Seurat)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(ggrepel)

obj = readRDS("all_subfamily.RDS")

obj@meta.data$region= sapply(strsplit(as.character(obj@meta.data$sample), "_"), "[[", 1)
obj@meta.data$age= sapply(strsplit(as.character(obj@meta.data$sample), "_"), "[[", 2)
obj@meta.data$rep= sapply(strsplit(as.character(obj@meta.data$sample), "_"), "[[", 3)

DimPlot(obj
    
        ,reduction="umap")

meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")


options(repr.plot.width=20, repr.plot.height=9)

DimPlot(obj,reduction="umap",group.by = "celltype",label = T)

head(obj@meta.data)


FeaturePlot(obj,"SoloTE-MuRRS-int",reduction="umap", split.by = "age",order = TRUE,pt.size = 3)


FeaturePlot(obj,"SoloTE-L1MA5A",reduction="umap", split.by = "age")

options(repr.plot.width=12, repr.plot.height=5)
VlnPlot(obj,"SoloTE-L1MA5A", pt.size = 0,group.by = "region", split.by = "age",y.max = 2)

VlnPlot(obj,"SoloTE-MuRRS-int",group.by = "region", split.by = "age", pt.size = 0, y.max = .0001)

head(ps@meta.data)


ps <- AggregateExpression(obj, assays = "RNA", return.seurat = T, group.by = c("celltype", "sample"))
ps$age = sapply(strsplit(as.character(ps$sample), "-"), "[[", 2)
ps$celltype_age = paste(ps$celltype,ps$age, sep = "_")


counts = ps@assays$RNA$counts

ct = 'L2/3 IT CTX Glut'
#ct = 'Oligo NN'
#ct = 'OPC NN'

cl = gsub("/", "-", ct)
cl = gsub(" ", "-", cl)

colnames(counts)[grep(ct, colnames(counts))]
counts_L23 = counts[,grep(ct, colnames(counts))]
counts_L23 = counts_L23[,-grep("9mo", colnames(counts_L23))]
groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)


groups

head(grep("Solo",rownames(counts_L23)))

length(rowSums((counts_L23) >= 10) >= 2 )
length(!grep("Solo",rownames(counts_L23)))

head( rowSums((counts_L23) >= 10) >= 2 & !grep("Solo",rownames(counts_L23)))

for (ct in unique(obj$celltype)) {
   # ct = 'IT AON-TT-DP Glut'
    cat(ct, "\n")
    cl = gsub("/", "-", ct)
    cl = gsub(" ", "-", cl)
    
    colnames(counts)[grep(ct, colnames(counts))]
    counts_L23 = counts[,grep(ct, colnames(counts))]
   # counts_L23 = counts_L23[grep("Solo" ,rownames(counts_L23)),]
    counts_L23 = counts_L23[,-grep("9mo", colnames(counts_L23))]
    groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
    groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
    regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
    regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)
    
    min_reads = 10000
    colSums(counts_L23)
    good = names(table(regs[which(colSums(counts_L23) > min_reads)]))[which(table(regs[which(colSums(counts_L23) > min_reads)]) > 3)]
    good
    if(length(good) < 1) {
        next
    }
    counts_L23 = counts_L23[,which(regs %in% good)]
    
    groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
    groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
    groups = factor(groups, levels = c("18mo", "02mo"))
    regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
    regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)
    
    y <- DGEList(counts = counts_L23, group = groups)
    
    keep <- rowSums((counts_L23) >= 10) >= 2
    length(keep)
    y <- y[keep,, keep.lib.sizes = FALSE]
    keep = grep("Solo",rownames(y$counts))
    y <- y[-keep,, keep.lib.sizes = F]

    y <- calcNormFactors(y)
    if(length(unique(regs)) > 1) {
        design <- model.matrix(~groups + regs)
    } else {
        design <- model.matrix(~groups)
    }
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = "groups02mo")
    
    fdr = p.adjust(lrt$table$PValue, method = "BH")
    out = cbind(lrt$table, fdr)
    out = out[order(out$fdr),]
    
    # out$logFC = -out$logFC
    out$significance <- ifelse(out$fdr < 0.05, "significant", "not significant")
    out$dir = "ns"
    out$dir[which(out$logFC > 0 & out$fdr < 0.05)] = "down"
    out$dir[which(out$logFC < 0 & out$fdr < 0.05)] = "up"
    out = out[order(out$fdr),]
    table(out$dir)
    
    cl = gsub(" ", "_", ct)
    cl = gsub("/", "-", cl)
    
    write.table(out, file = paste(cl, "DE_genes_only.txt", sep = "_"), sep = "\t")
    out_TE = out
    top10 <- out_TE[order(out_TE$PValue), ][1:10, ]
    top10$TE = rownames(top10)
    
    p = ggplot(out_TE, aes(x = -logFC, y = -log10(fdr), color = dir)) +
        geom_point() +
        scale_color_manual(values = c("up" = "red", "down" = "blue")) +
        theme_minimal() +
        labs(title = paste(cl, "Genes"),
             x = "log2 Fold Change",
             y = "-log10(P-Value)") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        # geom_vline(xintercept = c(-.25, .25), linetype = "dashed", color = "black") +
        theme(legend.position = "top") + # xlim(c(-.5,.5))+
        geom_text_repel(data = top10, aes(label = TE), vjust = 1.5, hjust = 1.5)
    p
    ggsave(p, file = paste(cl, "volcano_genes_only.pdf", sep = "_"))
}


for (ct in unique(obj$celltype)) {
   # ct = 'IT AON-TT-DP Glut'
    cat(ct, "\n")
    cl = gsub("/", "-", ct)
    cl = gsub(" ", "-", cl)
    
    colnames(counts)[grep(ct, colnames(counts))]
    counts_L23 = counts[,grep(ct, colnames(counts))]
   # counts_L23 = counts_L23[grep("Solo" ,rownames(counts_L23)),]
    counts_L23 = counts_L23[,-grep("9mo", colnames(counts_L23))]
    groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
    groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
    regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
    regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)
    
    min_reads = 10000
    colSums(counts_L23)
    good = names(table(regs[which(colSums(counts_L23) > min_reads)]))[which(table(regs[which(colSums(counts_L23) > min_reads)]) > 3)]
    good
    if(length(good) < 1) {
        next
    }
    counts_L23 = counts_L23[,which(regs %in% good)]
    
    groups = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
    groups = sapply(strsplit(as.character(groups), "-"), "[[", 2)
    groups = factor(groups, levels = c("18mo", "02mo"))
    regs = sapply(strsplit(as.character(colnames(counts_L23)), "_"), "[[", 2)
    regs = sapply(strsplit(as.character(regs), "-"), "[[", 1)
    
    y <- DGEList(counts = counts_L23, group = groups)
    
    keep <- rowSums((counts_L23) >= 10) >= 2
    length(keep)
    y <- y[keep,, keep.lib.sizes = FALSE]
    keep = grep("Solo",rownames(y$counts))
    y <- y[keep,, keep.lib.sizes = TRUE]

    y <- calcNormFactors(y)
    if(length(unique(regs)) > 1) {
        design <- model.matrix(~groups + regs)
    } else {
        design <- model.matrix(~groups)
    }
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef = "groups02mo")
    
    fdr = p.adjust(lrt$table$PValue, method = "BH")
    out = cbind(lrt$table, fdr)
    out = out[order(out$fdr),]
    
    # out$logFC = -out$logFC
    out$significance <- ifelse(out$fdr < 0.05, "significant", "not significant")
    out$dir = "ns"
    out$dir[which(out$logFC > 0 & out$fdr < 0.05)] = "down"
    out$dir[which(out$logFC < 0 & out$fdr < 0.05)] = "up"
    out = out[order(out$fdr),]
    table(out$dir)
    
    cl = gsub(" ", "_", ct)
    cl = gsub("/", "-", cl)
    
    write.table(out, file = paste(cl, "DE_TEs_only.txt", sep = "_"), sep = "\t")
    out_TE = out[grep("Solo", rownames(out)),]
    top10 <- out_TE[order(out_TE$PValue), ][1:10, ]
    top10$TE = rownames(top10)
    
    p = ggplot(out_TE, aes(x = -logFC, y = -log10(fdr), color = dir)) +
        geom_point() +
        scale_color_manual(values = c("up" = "red", "down" = "blue")) +
        theme_minimal() +
        labs(title = paste(cl, "TEs"),
             x = "log2 Fold Change",
             y = "-log10(P-Value)") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        # geom_vline(xintercept = c(-.25, .25), linetype = "dashed", color = "black") +
        theme(legend.position = "top") + # xlim(c(-.5,.5))+
        geom_text_repel(data = top10, aes(label = TE), vjust = 1.5, hjust = 1.5)
    p
    ggsave(p, file = paste(cl, "volcano_TEs_only.pdf", sep = "_"))
}


unique(obj$celltype)

setwd("subfamily_DESeq2/")

files = list.files(".", "only.txt", full.names = F)
cur = files[1]
tab = read.table(cur)
tab$gene = rownames(tab)
f = gsub("_DE_TEs_only.txt", "", files[1])

tab$celltype = f
#head(tab)
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_DE_TEs_only.txt", "", files[i])

    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}

tab$logFC = -tab$logFC

head(tab[which(tab$celltype == "Oligo_NN"),])

unique(tab$celltype)[which(!unique(tab$celltype) %in% ord)]


tab$type = "Gene"
tab$type[grep("SoloTE", rownames(tab))] = "SoloTE"
tab$direction = "Up in aging"
tab$direction[which(tab$logFC<0)] = "Down in aging"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("IMN", tab$celltype)]="NN"

#tab= tab[-grep("IT", tab$celltype),]

    
    
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","_", ord)


options(repr.plot.width=17, repr.plot.height=4)


d = tab
d = d[which(d$PValue<0.01),]
d = d[which(d$type!="Gene"),]
d$celltype = factor(d$celltype, levels = rev(ord))


g1 = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = position_stack(), color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "P < 0.05 All Subfamily TEs"))+
  theme(text = element_text(size = 18))
print(g1)
    
d = tab
d = d[which(d$PValue<0.05),]
d = d[which(d$type=="Gene"),]
d$celltype = factor(d$celltype, levels = rev(ord))


g = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "stack", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free",ncol = 3)+ ggtitle(paste( "P < 0.05 genes"))+
  theme(text = element_text(size = 18))
    
    #print(g)




options(repr.plot.width=18, repr.plot.height=9)

g1


options(repr.plot.width=16, repr.plot.height=6)
tab$celltype = factor(tab$celltype, levels = rev(ord))

ggplot(tab[which( tab$type=="SoloTE"),], 
  aes(x = celltype, y = logFC, 
  color = ifelse(`PValue` < 0.01, ifelse(`logFC` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential TE Accessibility",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "p < 0.05") +
  theme_minimal() + coord_flip() + facet_wrap(~clade, scales = "free", ncol = 3)

options(repr.plot.width=12, repr.plot.height=5)
tab$celltype = factor(tab$celltype, levels = rev(ord))

ggplot(tab[which( tab$type=="SoloTE"),], 
  aes(x = celltype, y = logFC, 
  color = ifelse(`PValue` < 0.01, ifelse(`logFC` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential TE Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "p < 0.01") +
  theme_minimal() + coord_flip() + facet_wrap(~clade, scales = "free", ncol = 3)

out = tab[which(tab$celltype=="L2-3_IT_ENT_Glut"),]
  out_TE = out[grep("Solo", rownames(out)),]
    top10 <- out_TE[order(out_TE$PValue), ][1:10, ]
    top10$TE = rownames(top10)
    
    p = ggplot(out_TE, aes(x = logFC, y = -log10(PValue), color = dir)) +
        geom_point() +
        scale_color_manual(values = c("up" = "red", "down" = "blue")) +
        theme_minimal() +
        labs(title = paste(cl, "TEs"),
             x = "log2 Fold Change",
             y = "-log10(P-Value)") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        # geom_vline(xintercept = c(-.25, .25), linetype = "dashed", color = "black") +
        theme(legend.position = "top") + # xlim(c(-.5,.5))+
        geom_text_repel(data = top10, aes(label = TE), vjust = 1.5, hjust = 1.5)
    p


ggplot(tab[which(abs(tab$logFC)>0.1 & tab$type=="Gene"),], 
  aes(x = celltype, y = logFC, 
  color = ifelse(`fdr` < 0.05, ifelse(`logFC` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential TE Accessibility",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip() + facet_wrap(~clade, scales = "free", ncol = 3)



up_tab = table(tab[which(tab$dir=="up"& tab$type == "SoloTE"),"gene"])
up_tab[order(up_tab,decreasing = T)]

head(tab)

tab$dir_pval= "ns"
tab$dir_pval[which(tab$fdr<0.1 &tab$logFC>0)] = "up"
tab$dir_pval[which(tab$fdr<0.1 &tab$logFC<0)] = "down"

lkey=read.table("~/projects/combined_all/scTE/out.saf")
#lkey = lkey[-which(lkey$V5=="Retroposon"),]
ukey = lkey[-which(duplicated(lkey$V1)),]
rownames(ukey) = ukey$V1
LINES = unique(paste(lkey[which(lkey$V6%in%c("LINE", "LINE?")),1]))
LTR = unique(paste(lkey[which(lkey$V6%in%c("LTR", "LTR?")),1]))
SINE = unique(paste(lkey[which(lkey$V6%in%c("SINE","SINE?")),1]))
Retroposon = unique(paste(lkey[which(lkey$V6%in%c("Retroposon")),1]))


ukey[grep("B1", ukey[,"V1"]),]
head(ukey)

ukey

bkey  = ukey[grep("_", rownames(ukey)), ]
rownames(bkey) = gsub("_", "-", rownames(bkey))
ukey = rbind(ukey,bkey)

ukey$V6 = gsub("[?]", "", ukey$V6)

sig= tab[which(tab$fdr<.05 & abs(tab$logFC) > 0.25 & tab$clade=="NN"),]


sig_genes = (unique(sig$gene))

sig_genes

head(tab)

sig_tab = tab[which(tab$gene %in% c(sig_genes)),]

library(pheatmap)
library(RColorBrewer)


options(repr.plot.width=15, repr.plot.height=12)
# Reshape the data for the heatmap
df = sig_tab
df$gene = gsub("SoloTE-", "", df$gene)
#df$gene = gsub("-", "_", df$gene)
heatmap_data <- dcast(df, gene ~ celltype, value.var = "logFC")
row.names(heatmap_data) <- heatmap_data$gene
heatmap_data <- heatmap_data[,-1]

df$type = ukey[paste(df$gene), "V6"]

# Annotations
gene_annotation <- df[, c("gene", "type")]
gene_annotation <- unique(gene_annotation)
row.names(gene_annotation) <- gene_annotation$gene
gene_annotation <- gene_annotation[,-1, drop = FALSE]

celltype_annotation <- df[, c("celltype", "clade")]
celltype_annotation <- unique(celltype_annotation)
row.names(celltype_annotation) <- celltype_annotation$celltype
celltype_annotation <- celltype_annotation[,-1, drop = FALSE]

# Colors for annotations
ann_colors <- list(
  type = brewer.pal(8, "Set2")[1:length(unique(gene_annotation$type))],
  clade = brewer.pal(8, "Set1")[1:length(unique(celltype_annotation$clade))]
)
names(ann_colors$type) <- unique(gene_annotation$type)
names(ann_colors$clade) <- unique(celltype_annotation$clade)



# Define the color palette for the heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Define breaks for the heatmap so that 0 is the midpoint
max_abs_val <- max(abs(heatmap_data), na.rm = TRUE)
breaks <- seq(-max_abs_val, max_abs_val, length.out = 101)

# Create the heatmap
pheatmap(
  as.matrix(heatmap_data),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = gene_annotation,
  annotation_col = celltype_annotation,
  annotation_colors = ann_colors,
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = heatmap_colors,
  breaks = breaks,
  main = "Heatmap of logFC for each Gene and Cell Type Combination"
)

options(repr.plot.width=16, repr.plot.height=19)

library(reshape2)
# Reshape the data for the heatmap
heatmap_data <- dcast(sig_tab, gene ~ celltype, value.var = "logFC")

# Convert the data to long format for ggplot2
heatmap_data_long <- melt(heatmap_data, id.vars = "gene")

# Create the heatmap
ggplot(heatmap_data_long, aes(x = variable, y = gene, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(min(heatmap_data_long$value, na.rm = TRUE), max(heatmap_data_long$value, na.rm = TRUE)),
                       name = "logFC") +
  theme_minimal() +
  labs(x = "Cell Type", y = "Gene", title = "Heatmap of logFC for each Gene and Cell Type Combination") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

table(tab[which(tab$gene == "SoloTE-L1MA5A"),"direction" ] )

table(tab[which(tab$gene == "SoloTE-"),"direction" ] )

tab[which(tab$gene == "SoloTE-MuRRS-int"), ]


