library(ggplot2)

meta = read.csv("../ABC/gene_meta.csv")
rownames(meta) = meta$geneID

head(meta)
gmeta= meta[which(!duplicated(meta$gene_name)),]
rownames(gmeta) = gmeta$gene_name
head(gmeta)

files = list.files("AMY_region_edgeR//", "DEG.txt")

files

f = files[11]
cat(f)
ct = gsub("_DEG.txt", "", f )
deg = read.table(paste("AMY_region_edgeR/", f, sep = ""))
head(deg)  

nrow(deg)

deg$dir = "NS"
deg$dir[which( deg$logFC>0.15 & deg$PValue<0.05)] = "Down"
    deg$dir[which( deg$logFC< -0.15 & deg$PValue<0.05 )] = "Up"
    table(deg$dir)
    deg$chr = gmeta[rownames(deg), "chrom"]
    deg$start = gmeta[rownames(deg), "start"]
    deg$end = gmeta[rownames(deg), "end"]
    deg$length = deg$end-deg$start
    head(deg)


options(repr.plot.width=8, repr.plot.height=4.5)

ggplot(deg, aes(x = length, color = dir
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of length",
       x = "length",
       y = "eCDF",
       color = "Type") +#xlim(c(-15,15))+
  theme_minimal()


options(repr.plot.width=8, repr.plot.height=4.5)

ggplot(deg, aes(x = logCPM, color = dir
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of logCPM",
       x = "logCPM",
       y = "eCDF",
       color = "Type") +#xlim(c(-15,15))+
  theme_minimal()


files = list.files("DEG_results_latent_rep_mito_together/", ".csv")

head(deg)

all = list()
for(f in files){ 
    ct = gsub(".csv", "", f )
    deg = read.csv(paste("DEG_results_latent_rep_mito_together/", f, sep = ""))
    
    deg$dir = "NS"
    deg$dir[which( deg$avg_log2FC>0.15 & deg$p_val_adj<0.01)] = "Down"
    deg$dir[which( deg$avg_log2FC< -0.15 & deg$p_val_adj<0.01 )] = "Up"
    table(deg$dir)
    deg$chr = gmeta[paste(deg$X), "chrom"]
    deg$start = gmeta[paste(deg$X), "start"]
    deg$end = gmeta[paste(deg$X), "end"]
    deg$length = deg$end-deg$start
    head(deg)
    deg$ct = ct
    options(repr.plot.width=5, repr.plot.height=4)
    g = ggplot(deg, aes( x=dir, y=log(length))) +
    geom_boxplot( alpha=0.5)+
    labs(title=ct,x="dir", y = "log length")+
    theme_classic()
    print(g)
    all[[ct]] = deg
}

diff = do.call(rbind, all)


head(diff)
diff$major = rep("Gaba", nrow(diff))
diff$major[grep("Glut" , diff$ct)] = "Glut"
diff$major[grep("NN" , diff$ct)] = "NN"
diff$major[grep("IMN" , diff$ct)] = "IMN"

diff = diff[which(!is.na(diff$length)),]

head(diff[which(diff$ct=="Oligo_NN"),])


diff$dir2 = "NS"
#deg$dir[which( deg$avg_log2FC>0.15)[1:1000]] = "Up"
#deg$dir[which( deg$avg_log2FC< -0.15)[1:1000]] = "Down"

diff$dir2[which( diff$avg_log2FC>1 & diff$p_val_adj<0.001)] = "Down"
diff$dir2[which( diff$avg_log2FC< -.5 & diff$p_val_adj<0.001 )] = "Up"
table(diff$dir2)


il = diff[grep("Il", diff$X),]
il[which(il$dir !="NS"),]

Sirt = diff[grep("Sirt", diff$X),]
Sirt[which(Sirt$dir !="NS"),]

diff$pct = (diff$`pct.1`+diff$`pct.2`)/2

diff_seurat = diff

options(repr.plot.width=8, repr.plot.height=4.5)

ggplot(diff[which(diff$ct == "STR_D12_Gaba"),], aes(x = pct, color = dir2
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of pct",
       x = "pct",
       y = "eCDF",
       color = "Type") +#xlim(c(-15,15))+
  theme_minimal()


options(repr.plot.width=8, repr.plot.height=4.5)

ggplot(diff[which(diff$ct == "STR_D12_Gaba"),], aes(x = log(length), color = dir2
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log(gene_length)",
       x = "log(gene_length)",
       y = "eCDF",
       color = "Type") +#xlim(c(-15,15))+
  theme_minimal()


options(repr.plot.width=8, repr.plot.height=4.5)

ggplot(diff, aes(x = pct, color = dir2
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of pct",
       x = "pct",
       y = "eCDF",
       color = "Type") +#xlim(c(-15,15))+
  theme_minimal()+ facet_wrap(~major)


options(repr.plot.width=8, repr.plot.height=4.5)

ggplot(diff, aes(x = log(length), color = dir2
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log(gene_length)",
       x = "log(gene_length)",
       y = "eCDF",
       color = "Type") +#xlim(c(-15,15))+
  theme_minimal()+ facet_wrap(~major)


head(diff[which(diff$ct == "STR_D12_Gaba"),])


options(repr.plot.width=8, repr.plot.height=4.5)

ggplot(diff[which(diff$ct == "DG_Glut"),], aes(x = pct, color = dir
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of pct",
       x = "pct",
       y = "eCDF",
       color = "Type") +#xlim(c(-15,15))+
  theme_minimal()


options(repr.plot.width=8, repr.plot.height=4.5)

ggplot(diff, aes(x = log(length), color = dir
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log(gene_length)",
       x = "log(gene_length)",
       y = "eCDF",
       color = "Type") +#xlim(c(-15,15))+
  theme_minimal()+ facet_wrap(~major)


deg = read.csv("DEG_results_latent_rep_mito_together/Oligo_NN.csv")


deg$dir = "NS"
#deg$dir[which( deg$avg_log2FC>0.15)[1:1000]] = "Up"
#deg$dir[which( deg$avg_log2FC< -0.15)[1:1000]] = "Down"

deg$dir[which( deg$avg_log2FC>0.15 & deg$p_val_adj<0.01)] = "Down"
deg$dir[which( deg$avg_log2FC< -0.15 & deg$p_val_adj<0.01 )] = "Up"
table(deg$dir)


deg$chr = gmeta[paste(deg$X), "chrom"]
deg$start = gmeta[paste(deg$X), "start"]
deg$end = gmeta[paste(deg$X), "end"]
deg$length = deg$end-deg$start
head(deg)




