library(ggplot2)


meta = read.csv("../final_meta.csv")
meta$celltype_final = gsub("/","-", meta$celltype_final)
meta$celltype_final = gsub("Lymphoid","Lymphoid NN", meta$celltype_final)
meta = meta[-which(meta$celltype_final=="doublet"),]
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)


head(meta[which(meta$batch=="Female"),])

d12 = read.table("../female_RNA/subclustering3/D12MSN_split_meta.txt")

rownames(d12) = paste("Female:", d12$barcode, "-1",sep = "")

head(d12)

mat = match(meta$cell_id,rownames(d12))

head(mat)

length(mat)

meta$label = meta$celltype_final


meta$label= d12$lab[mat]

table(meta$label)

head(meta)

meta$label[which(is.na(meta$label))]= paste(meta$celltype_final[which(is.na(meta$label))])

head(meta$label)

dsub = read.csv("subcluster/STR D12 Gaba_subcluster_meta.csv")

head(dsub)

mat = match(dsub$cell_id,rownames(d12))

dsub$label= d12$lab[mat]

nrow(dsub)

ggplot(dsub[which(!is.na(dsub$label)),],aes(x=umap_x,y=umap_y,color = label)) + geom_point(alpha = 0.9)


ggplot(dsub[which(!is.na(dsub$label)),],aes(x=umap_x,y=umap_y,color = label)) + geom_point(alpha = 0.9)+theme_classic()


ggplot(dsub[which(!is.na(dsub$label)),],aes(x=umap_x,y=umap_y,color = age)) + geom_point(alpha = 0.7)+theme_classic()


ggplot(dsub[which(!is.na(dsub$label)),],aes(x=umap_x,y=umap_y,color = region)) + geom_point(alpha = 0.7)+theme_classic()



