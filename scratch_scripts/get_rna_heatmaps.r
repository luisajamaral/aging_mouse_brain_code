library(ggplot2)
library(data.table)
library(UpSetR)
library(pheatmap)
library(Seurat)
library(enrichR)

obj = readRDS("RNA_final_SCT.RDS")
obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)

length(unique(enrichment_results$log_p_value))

files = list.files("DEG_results_region/", ".csv", full.names = T)
cfiles = list.files("DEG_results_.01_.01/", ".csv", full.names = T)
cl = "Glut"
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
    significant_changes <- subset(dar_data, p_val_adj < 0.05 & abs(avg_log2FC)>0.15 & (`pct.1`>0.05 | `pct.2`>0.05))
    atabs[[f]] = significant_changes
    #list_genes[[dar_data$cl[1]]] = rownames(significant_changes)
    list_genes[[dar_data$cl[1]]] = significant_changes$X

    
}
tabs = do.call(rbind, atabs)
tabs = tabs[which(tabs$p_val_adj<0.05& abs(tabs$avg_log2FC)>0.15 & (tabs$`pct.1`>0.01 | tabs$`pct.2`>0.01)),]
genes = unique(tabs$X)
length(genes)
upset(fromList(list_genes), order.by = "freq", nsets = 20)



#obj$celltype_sample = paste(obj$celltype_final, obj$orig.ident)
#Idents(obj) = "celltype_region_age"
ct = "Glut"
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


ph = pheatmap(oligo, kmeans_k = 6,show_rownames = T,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps, cluster_rows = F, breaks = brks,main = length(genes))


#ph = pheatmap(oligo, kmeans_k = 6,show_rownames = F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
#        gaps_col = gaps, breaks = brks,main = length(genes))


pheatmap(oligo, show_rownames = F,cutree_rows = 6,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps, breaks = brks, main = length(genes))


pheatmap(oligo, show_rownames = F,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps, breaks = brks)






head(annotation_row)

hp=pheatmap(oligo, show_rownames = F,cutree_rows = 4,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps, breaks = brks, main = length(genes))


cluster = cutree(hp$tree_row, k = 4)
annotation_row = data.frame(
    cluster = paste(cluster)
  )
rownames(annotation_row) = names(cluster)

hp=pheatmap(oligo, show_rownames = F,cutree_rows = 4,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), ,annotation_row = annotation_row, annotation_col = annotation_col, cluster_cols = F,
        gaps_col = gaps, breaks = brks, main = length(genes))


pheatmap(oligo, show_rownames = F,cutree_rows = 4,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), ,annotation_row = annotation_row, annotation_col = annotation_col, cluster_cols = T,
        gaps_col = gaps, breaks = brks, main = length(genes))


cluster["Nrf1"]

websiteLive <- getOption("enrichR.live")


websiteLive

#cl4 = names(ph$kmeans$cluster[which(ph$kmeans$cluster==3)])
cl4 = names(cluster[which(cluster==4)])
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human", "SynGO_2022")
enriched <- enrichr(cl4, dbs)



if (websiteLive) head(enriched[["SynGO_2022"]],10)


if (websiteLive) head(enriched[["KEGG_2021_Human"]],10)


if (websiteLive) head(enriched[["GO_Biological_Process_2023"]],30)


library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507


listAttributes(ensembl)

gene.data <- getBM(attributes=c('ensembl_gene_id','mgi_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0006119', mart = ensembl)

unique(gene.data$mgi_symbol)

length(unique(gene.data$mgi_symbol))

synvess = fread("ox_phos.txt", sep = "\t")
head(synvess)

gs = unique(synvess$`MGI Gene/Marker ID`)

length(gs)

obj = AddModuleScore(obj, features =list(gs), name = "Ox Phos")

obj = AddModuleScore(obj, features =list(unique(gene.data$mgi_symbol)), name = "Ox Phos")

obj$`Ox Phos1`

options(repr.plot.width=12, repr.plot.height=5)

Idents(obj)= "celltype_final"
plot <- FeaturePlot(obj,
                features = "Nrf1", split.by = "age",label = F, repel = TRUE,order = TRUE)# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot

plot <- VlnPlot(obj,
                features = "Nrf1", ,group.by = "region",split.by = "age")# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot

obj$major = "NN"
obj$major[grep("Gaba",obj$celltype_final)] = "Gaba"
obj$major[grep("Glut",obj$celltype_final)] = "Glut"
obj$major[grep("Neur",obj$celltype_final)] = "RLP"


options(repr.plot.width=8, repr.plot.height=5)


plot <- VlnPlot(obj,
                features = "Nrf1", ,group.by = "major",split.by = "age")# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot

options(repr.plot.width=12, repr.plot.height=5)

plot <- VlnPlot(obj[,which(obj$major=="NN")],
                features = "Nrf1", split.by = "age")# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot

options(repr.plot.width=12, repr.plot.height=5)

plot <- VlnPlot(obj[,which(obj$major=="Glut")],
                features = "Nrf1", split.by = "age",)# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot

options(repr.plot.width=12, repr.plot.height=5)

plot <- VlnPlot(obj[,which(obj$major=="Gaba")],
                features = "Nrf1", split.by = "age")# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot


plot <- VlnPlot(obj[,which(obj$celltype_final=="Oligo NN")],
                features = "Nrf1", group.by  = "region",split.by = "age")# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot

unique(obj@meta.data$celltype_final)

ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="Oligo NN"),], aes(x=region, fill=age ,y=`Ox Phos1`)) + 
  geom_boxplot(notch=TRUE)+ggtitle("Oligo Ox Phos")
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


d12= obj@meta.data[which(obj@meta.data$celltype_final=="DG Glut"),]
d12 = d12[which(d12$region %in%c("HCA","HCP")),]
ggplot(d12, aes(x=region, fill=age ,y=`Ox Phos1`)) + 
  geom_boxplot(notch=TRUE)+ggtitle("DG Glut Ox Phos")


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="STR D12 Gaba"),], aes(x=region, fill=age ,y=`Ox Phos1`)) + 
  geom_boxplot(notch=TRUE)+ggtitle("STR D12 Gaba Ox Phos")
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="DG Glut"),], aes(x=region, fill=age ,y=`Ox Phos1`)) + 
  geom_boxplot(notch=TRUE)+ggtitle("Oligo Ox Phos")
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="Oligo NN"),], aes(x=region, fill=age,color = paste(rep) ,y=`Ox Phos1`)) + 
  geom_boxplot(notch=TRUE)+ggtitle("Oligo Ox Phos")
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="Oligo NN"),], aes(x=region, fill=age ,y=`Ox Phos1`)) + 
  geom_boxplot(notch=TRUE)+ggtitle("Oligo Ox Phos")
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="STR D12 Gaba"),], aes(x=region, fill=age,color = paste(rep) ,y=`Ox Phos1`)) + 
  geom_boxplot(notch=TRUE)+ggtitle("STR D12 Gaba Ox Phos1")
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="DG Glut"),], aes(x=region, fill=age,color = paste(rep) ,y=`Cellular respiration1`)) + 
  geom_boxplot(notch=TRUE)+ggtitle("DG Glut Cellular respiration1")
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="Oligo NN"),], aes(x=region, fill=age ,y=`Cytoplasmic Translation1`)) + 
  geom_boxplot(notch=TRUE)
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="Oligo NN"),], aes(x=region, fill=paste(rep), color = age ,y=`Cytoplasmic Translation1`)) + 
  geom_boxplot(notch=TRUE)
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="OPC NN"),], aes(x=region, fill=age ,y=`Double-Strand Break Repair1`)) + 
  geom_boxplot(notch=TRUE)
  

ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="Microglia NN"),], aes(x=region, fill=age ,y=`Double-Strand Break Repair1`)) + 
  geom_boxplot(notch=TRUE)
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


ggplot(obj@meta.data[which(obj@meta.data$celltype_final=="Microglia NN"),], aes(x=region, color = age,fill= paste(rep),y=`DNA Damage Response1`)) + 
  geom_boxplot(notch=TRUE)
  
  # Add median lines # geom_signif(comparisons = list(c("2mo", "18mo")), map_signif_level = TRUE, textsize = 5, vjust = -0.5)


module_scores <- obj@meta.data$`Double-Strand Break Repair1`

# Extract data for each group
group_data <- split(module_scores, list(obj$region, obj$age))

# Calculate median and quartiles for each group
medians <- sapply(group_data, function(x) median(x))
q1 <- sapply(group_data, function(x) quantile(x, 0.25))
q3 <- sapply(group_data, function(x) quantile(x, 0.75))

# Perform statistical test (e.g., Wilcoxon rank-sum test)
p_value <- wilcox.test(module_scores[obj$region == "RLP" & obj$age == "2mo"],
                       module_scores[obj$region == "RLP" & obj$age == "18mo"])$p.value

p_value

options(repr.plot.width=12, repr.plot.height=5)

plot <- VlnPlot(obj[,which(obj$celltype_final=="Microglia NN")],
                features = "Double-Strand Break Repair1", group.by = "region",split.by = "age", pt.size = 0)# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot

options(repr.plot.width=12, repr.plot.height=5)

plot <- VlnPlot(obj[,which(obj$celltype_final=="Microglia NN")],
                features = "DNA Damage Response1", group.by = "region",split.by = "age", pt.size = 0)# +
               # scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot


