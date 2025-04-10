library(ggplot2)
meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")
setwd("~/projects/combined_all/female_RNA/SoloTE/HCP_subfamily/")


#for(reg in c("FC_", "ENT_", "AMY_", "RLP_", "NAC_","HCA_","HCP_", "CP_")) { 

for(reg in c( "NAC_","HCA_","HCP_", "CP_")) { 
    
    #reg = "FC_"
setwd(paste("../",reg,"subfamily/", sep = ""))
files = list.files(".", "TE.txt", full.names = F)
cur = files[1]
tab = read.table(cur)
tab$gene = rownames(tab)
f = gsub("_subfamily_DEG_TE.txt", "", files[1])
f = gsub(reg, "", f)

tab$celltype = f
#head(tab)
for(i in 2:length(files)){
    cur = read.table(files[i])

    cur$celltype = gsub("_subfamily_DEG_TE.txt", "", files[i])
    cur$celltype = gsub(reg, "", cur$celltype)

    cur$gene = rownames(cur)
    tab = rbind(tab, cur)
        
}
#nrow(tab)
#length(which(abs(tab$PValue)<0.05))
#head(tab)
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

#tab= tab[-grep("IT", tab$celltype),]

    
    
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","-", ord)
reg = gsub("_","", reg)


options(repr.plot.width=17, repr.plot.height=4)


d = tab
d = d[which(d$PValue<0.05),]
d = d[which(d$type!="Gene"),]
d$celltype = factor(d$celltype, levels = rev(ord))


g1 = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = position_stack(), color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~clade,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.05 All Subfamily TEs"))+
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
  coord_flip() + facet_wrap(~clade,scales = "free_y",ncol = 3)+ ggtitle(paste(reg, "P < 0.05 genes"))+
  theme(text = element_text(size = 18))
    
    #print(g)
}



(unique(tab$gene[grep("Solo", tab$gene)]))

d$gene

tab[grep("SoloTE-MuRRS-int", tab$gene),]




