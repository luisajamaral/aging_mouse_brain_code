library(ggplot2)

library(ggplot2)
meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")


setwd("~/projects/combined_all/female_RNA/SoloTE/FC_locus/")

reg = "NAC_"
setwd(paste("../",reg,"locus/", sep = ""))
files = list.files(".", "TE.txt", full.names = F)
cur = files[1]
tab = read.table(cur)
tab$gene = rownames(tab)
f = gsub("_DEG_TE.txt", "", files[1])
f = gsub(reg, "", f)

tab$celltype = f
#head(tab)
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_DEG_TE.txt", "", files[i])
    cur$celltype = gsub(reg, "", cur$celltype)

    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}
nrow(tab)
length(which(abs(tab$PValue)<0.01))

tab$type = "Gene"
tab$type[grep("SoloTE", rownames(tab))] = "SoloTE"
tab$direction = "Up in aging"
tab$direction[which(tab$logFC>0)] = "Down in aging"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("Ex-IMN", tab$celltype)]="Glut"
tab$clade[grep("Inh-IMN", tab$celltype)]="Gaba"

tab= tab[-grep("IT", tab$celltype),]

ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","-", ord)
reg = gsub("_","", reg)


options(repr.plot.width=17, repr.plot.height=4)


d = tab
d = d[which(d$PValue<0.05),]
d = d[which(d$type!="Gene"),]
d$celltype = factor(d$celltype, levels = rev(ord))


ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.05 All Locus TEs"))+
  theme(text = element_text(size = 18))

d = tab
d = d[which(d$PValue<0.05),]
d = d[which(d$type=="Gene"),]
d$celltype = factor(d$celltype, levels = rev(ord))


ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.05 genes"))+
  theme(text = element_text(size = 18))


d = tab
d = d[which(d$PValue<0.05),]
d = d[grep("SoloTE",d$gene),]
d$celltype = factor(d$celltype, levels = rev(ord))
te_type = rep("other",nrow(d))
te_type[grep("LINE", d$gene)] = "LINE"
te_type[grep("SINE", d$gene)] = "SINE"

d$te_type = te_type
d = d[grep("other",d$te_type),]

ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.05 LINE elements"))+
  theme(text = element_text(size = 18))


d = tab
d = d[which(d$PValue<0.05),]
d = d[grep("SoloTE",d$gene),]
d$celltype = factor(d$celltype, levels = rev(ord))
te_type = rep("other",nrow(d))
te_type[grep("LINE", d$gene)] = "LINE"
te_type[grep("SINE", d$gene)] = "SINE"
d$te_type = te_type
d = d[grep("SINE",d$te_type),]

ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.05 SINE elements"))+
  theme(text = element_text(size = 18))



options(repr.plot.width=6, repr.plot.height=4)

d = tab
d = d[which(d$PValue<0.05),]
#d = d[grep("SoloTE",d$gene),]
d$celltype = factor(d$celltype, levels = rev(ord))
te_type = rep("gene",nrow(d))
te_type[grep("SoloTE", d$gene)] = "other TE"
te_type[grep("LINE", d$gene)] = "LINE"
te_type[grep("SINE", d$gene)] = "SINE"
te_type = factor(te_type, levels = c("LINE","SINE", "other TE", "gene"))
d$te_type= te_type
table(d$direction, d$te_type)
counts = data.frame(table(d$direction, d$te_type, d$clade))
colnames(counts) = c("Category", "Type", "Clade","Count")

counts <- counts %>%
  group_by(Type) %>%
  mutate(Percentage = Count / sum(Count) * 100)

counts

# Create plot
ggplot(counts, aes(x = Type, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Overlap Percentages by Category and Type",
       y = "Percentage",
       x = "TE type") +
  scale_fill_manual(values = c("#FF9999", "#66CCFF", "#99FF99")) +
  theme_minimal()+
  theme(text = element_text(size = 16)) + ggtitle(reg)

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC>0 & tab$type == "SoloTE"),]
down_tab$TE=sapply(strsplit(as.character(down_tab$gene), "-"), "[[", 5)
df = as.data.frame(t(table(down_tab$TE)[order(table(down_tab$TE), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("TE","Freq")
head(df,10)
hist(table(down_tab$TE),breaks = 10, main = "Number of celltypes vs TEs down DE",col = "tomato")

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC<0 & tab$type == "SoloTE"),]
down_tab$TE=paste( sapply(strsplit(as.character(down_tab$gene), "-"), "[[", 5))
df = as.data.frame(t(table(down_tab$TE)[order(table(down_tab$TE), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("TE","Freq")
head(df,20)
hist(table(down_tab$TE),breaks = 10, main = "Number of celltypes vs TEs down DE",col = "tomato")

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC>0 & tab$type == "SoloTE"),]
down_tab$TE=paste( sapply(strsplit(as.character(down_tab$gene), "-"), "[[", 5))
df = as.data.frame(t(table(down_tab$TE)[order(table(down_tab$TE), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("TE","Freq")
head(df,20)
hist(table(down_tab$TE),breaks = 10, main = "Number of celltypes vs TEs down DE",col = "tomato")

head(tab)

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC>0 & tab$type == "SoloTE"),]
down_tab$TE=paste(sapply(strsplit(as.character(down_tab$gene), "-"), "[[", 5), sapply(strsplit(as.character(down_tab$gene), "-"), "[[", 6))
df = as.data.frame(t(table(down_tab$TE)[order(table(down_tab$TE), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("TE","Freq")
head(df,20)
hist(table(down_tab$TE),breaks = 10, main = "Number of celltypes vs TEs down DE",col = "tomato")

head(down_tab$gene)

unique(sapply(strsplit(as.character(down_tab$gene), "-"), "[[", 5))


down_tab = tab[which(tab$logFC>0 ),]
df = as.data.frame(t(table(down_tab$gene)[order(table(down_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("Gene","Freq")
df$Gene = sapply(strsplit(as.character(df$Gene), "[;]"), "[[", 1)
head(df,15)
hist(table(down_tab$gene),breaks = 10, main = "Number of celltypes vs genes down DE",col = "tomato")

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC<0 ),]
df = as.data.frame(t(table(down_tab$gene)[order(table(down_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("Gene","Freq")
df$Gene = sapply(strsplit(as.character(df$Gene), "[;]"), "[[", 1)
head(df,15)
hist(table(down_tab$gene),breaks = 10, main = "Number of celltypes vs genes up DE",col = "tomato")

s

setwd("~/projects/combined_all/female_RNA/FC_region_edgeR/")

reg = "AMY_"
setwd(paste("../",reg,"region_edgeR/", sep = ""))
files = list.files(".", "DEG.txt", full.names = F)
cur = files[1]
tab = read.table(cur)
tab$gene = rownames(tab)
f = gsub("_DEG.txt", "", files[1])
f = gsub(reg, "", f)

tab$celltype = f
head(tab)
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_DEG.txt", "", files[i])
    cur$celltype = gsub(reg, "", cur$celltype)

    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}
nrow(tab)
length(which(abs(tab$PValue)<0.01))

tab$type = "Gene"
tab$direction = "Up in aging"
tab$direction[which(tab$logFC>0)] = "Down in aging"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("Ex-IMN", tab$celltype)]="Glut"
tab$clade[grep("Inh-IMN", tab$celltype)]="Gaba"

ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","-", ord)
reg = gsub("_","", reg)
options(repr.plot.width=15, repr.plot.height=6)

d = tab
d = d[which(d$PValue<0.025),]
d$celltype = factor(d$celltype, levels = rev(ord))


ggplot(d, aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.01"))+
  theme(text = element_text(size = 18))


d = tab
d = d[which(d$PValue<0.05),]
d$celltype = factor(d$celltype, levels = rev(ord))

ggplot(d, aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.05"))+
  theme(text = element_text(size = 18))

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC>0 ),]
df = as.data.frame(t(table(down_tab$gene)[order(table(down_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("Gene","Freq")
df$Gene = sapply(strsplit(as.character(df$Gene), "[;]"), "[[", 1)
head(df,15)
hist(table(down_tab$gene),breaks = 10, main = "Number of celltypes vs genes down DE",col = "tomato")

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC<0 ),]
df = as.data.frame(t(table(down_tab$gene)[order(table(down_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("Gene","Freq")
df$Gene = sapply(strsplit(as.character(df$Gene), "[;]"), "[[", 1)
head(df,15)
hist(table(down_tab$gene),breaks = 10, main = "Number of celltypes vs genes up DE",col = "tomato")

cat(df[which(df$Freq>2),"Gene"], sep = "\n")

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC>0 ),]
df = as.data.frame(t(table(down_tab$gene)[order(table(down_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("Gene","Freq")
df$Gene = sapply(strsplit(as.character(df$Gene), "[;]"), "[[", 1)
head(df,15)
hist(table(down_tab$gene),breaks = 10, main = "Number of celltypes vs genes down DE",col = "tomato")

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC<0 ),]
df = as.data.frame(t(table(down_tab$gene)[order(table(down_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("Gene","Freq")
df$Gene = sapply(strsplit(as.character(df$Gene), "[;]"), "[[", 1)
head(df,15)
hist(table(down_tab$gene),breaks = 10, main = "Number of celltypes vs genes up DE",col = "tomato")


