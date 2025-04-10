library(ggplot2)

#setwd("subfamily_analysis/")
files = list.files(".", "TE.txt", full.names = F)
cur = files[1]
tab = read.table(cur)
tab$gene = rownames(tab)
head(tab)
f = gsub("_DEG_TE.txt", "", files[1])
tab$celltype = f

for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_DEG_TE.txt", "", files[i])
    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}

table(sign(tab$logFC.groups02mo) == sign(tab$logFC.groups09mo))
table(abs(tab$logFC.groups02mo)>0.5)
#tab = tab[which(sign(tab$logFC.groups02mo) == sign(tab$logFC.groups09mo) & abs(tab$logFC.groups02mo)>0.25),]
tab = tab[which(abs(tab$logFC.groups02mo)>0.5),]
tab = tab[which(abs(tab$PValue)<0.01),]

table(abs(tab$logFC.groups02mo)>0.25)

tab$type = "Gene"
tab$type[grep("SoloTE", rownames(tab))] = "SoloTE"
tab$direction = "Up in aging"
tab$direction[which(tab$logFC.groups02mo>0)] = "Down in aging"
tab$clade="Glut"
tab$clade[grep("NN", tab$celltype)]="NN"
tab$clade[grep("Gaba", tab$celltype)]="Gaba"
tab$clade[grep("Neur", tab$celltype)]="Gaba"
tab$clade[grep("Ex-IMN", tab$celltype)]="Glut"
tab$clade[grep("Inh-IMN", tab$celltype)]="Gaba"


meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")

ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","-", ord)



7.308934e-01+ 5.920344e-02+ 2.099031e-01

meta$region_type = "none"
meta$region_type[which(meta$region %in% c("FC","ENT"))]="Cortex"
meta$region_type[which(meta$region %in% c("HCA","HCP"))]="Hippocampus"

table(meta$region_type)
t = table(meta$celltype_final,meta$region_type)/rowSums(table(meta$celltype_final,meta$region_type))
t
cortical= rownames(t[which(t[,1]>0.75),])
hippo= rownames(t[which(t[,2]>0.75),])


cortical = gsub("/","-", cortical)
cortical = gsub(" ","-", cortical)

hippo = gsub("/","-", hippo)
hippo = gsub(" ","-", hippo)


options(repr.plot.width=10, repr.plot.height=14)

d = tab
d$celltype = factor(d$celltype, levels = rev(ord))
d$region_type = "NA"
d$region_type[which(d$celltype%in%cortical)]="Cortex"
d$region_type[which(d$celltype%in%hippo)]="Hippocampus"

ggplot(d, aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade+type,scales = "free",ncol = 2)

head(tab)

options(repr.plot.width=7, repr.plot.height=5)

d = tab
d$celltype = factor(d$celltype, levels = ord)
d$region_type = "NA"
d$region_type[which(d$celltype%in%cortical)]="Cortex"
d$region_type[which(d$celltype%in%hippo)]="Hippocampus"

ggplot(d[which(d$region_type=="Cortex"),], aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~type+clade, scales = "free")

options(repr.plot.width=7, repr.plot.height=4)

d = tab
d$celltype = factor(d$celltype, levels = ord)
d$region_type = "NA"
d$region_type[which(d$celltype%in%cortical)]="Cortex"
d$region_type[which(d$celltype%in%hippo)]="Hippocampus"

ggplot(d[which(d$region_type=="Hippocampus"),], aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~type+clade, scales = "free")

options(repr.plot.width=7, repr.plot.height=3)

d = tab
d$celltype = factor(d$celltype, levels = ord)
d$region_type = "NA"
d$region_type[which(d$celltype%in%cortical)]="Cortex"
d$region_type[which(d$celltype%in%hippo)]="Hippocampus"

ggplot(d[which(d$clade=="NN"),], aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~type, scales = "free")

length(unique(tab$gene))



length(up_tab$gene)
length(unique(up_tab$gene))


options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC.groups02mo>0 & tab$type == "Gene"),]
df = as.data.frame(t(table(down_tab$gene)[order(table(down_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("Gene","Freq")
df$Gene = sapply(strsplit(as.character(df$Gene), "[;]"), "[[", 1)
head(df,25)
hist(table(down_tab$gene),breaks = 30, main = "Number of celltypes vs genes down DE",col = "tomato")

options(repr.plot.width=4, repr.plot.height=4)

down_tab = tab[which(tab$logFC.groups02mo>0 & tab$type == "SoloTE"),]
df = as.data.frame(t(table(down_tab$gene)[order(table(down_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("TE","Freq")
head(df,25)
hist(table(down_tab$gene),breaks = 10, main = "Number of celltypes vs TEs down DE",col = "tomato")

options(repr.plot.width=4, repr.plot.height=4)

up_tab = tab[which(tab$logFC.groups02mo<0 & tab$type == "Gene"),]
df = as.data.frame(t(table(up_tab$gene)[order(table(up_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("Gene","Freq")
df$Gene = sapply(strsplit(as.character(df$Gene), "[;]"), "[[", 1)
head(df,25)
hist(table(up_tab$gene),breaks = 30,main = "Number of celltypes vs genes up DE" ,col = "tomato")

options(repr.plot.width=4, repr.plot.height=4)

up_tab = tab[which(tab$logFC.groups02mo<0 & tab$type == "SoloTE"),]
df = as.data.frame(t(table(up_tab$gene)[order(table(up_tab$gene), decreasing = T)]))
df = df[,2:3]
colnames(df)=c("TE","Freq")
head(df,25)
hist(table(up_tab$gene),breaks = 10,main = "Number of celltypes vs TE up DE" ,col = "tomato")

tab[which(tab$gene=="SoloTE-Lx2B"),]

up_tab = tab[which(tab$logFC.groups02mo<0),]
table(up_tab$gene)[order(table(up_tab$gene), decreasing = T)[1:15]]

hist(table(up_tab$gene))





unique(tab$celltype)

dg = tab[which(tab$celltype=="DG-Glut"),]
dg$Significant <- ifelse(dg$PValue < 0.01, "Yes", "No")

top15 <- dg[order(dg$PValue), ][1:15, ]

# Create an MA plot with color and labels for the top 15 genes
ggplot(dg, aes(x = logCPM, y = logFC.groups02mo, color = Significant, label = gene)) +
  geom_point() +
  geom_text_repel(data = top15, aes(label = gene), nudge_y = 0.2) +
  labs(title = "MA Plot with Significance (Top 15 Genes)",
       x = "logCPM",
       y = "logFC(groups02mo)") +
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) +
  theme_minimal()

options(repr.plot.width=10, repr.plot.height=7)

d = tab
d$celltype = factor(d$celltype, levels = ord)
d$region_type = "NA"
d$region_type[which(d$celltype%in%cortical)]="Cortex"
d$region_type[which(d$celltype%in%hippo)]="Hippocampus"

ggplot(d[which(d$clade=="Glut"),], aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~type, scales = "free")

options(repr.plot.width=10, repr.plot.height=6)

d = tab
d$celltype = factor(d$celltype, levels = ord)
d$region_type = "NA"
d$region_type[which(d$celltype%in%cortical)]="Cortex"
d$region_type[which(d$celltype%in%hippo)]="Hippocampus"

ggplot(d[which(d$clade=="Gaba"),], aes(x = celltype, fill = direction)) +
  geom_bar(stat = "count", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~type, scales = "free")


