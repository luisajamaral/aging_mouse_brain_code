library(ggplot2)

meta = read.csv("~/projects/combined_all/female_RNA/meta_final.csv")
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final),decreasing = T)])
ord = gsub("/","-", ord)
ord = gsub(" ","_", ord)

setwd("../DEG_results_latent_rep_mito/")
files = list.files(".", ".csv", full.names = F)

curr = read.csv(files[1])

curr$celltype = strsplit(files[1] , "--")[[1]][1]
curr$gene = curr$X
curr$region = strsplit(files[1] , "--")[[1]][2]
curr$region = gsub(".csv","", curr$region)

head(curr)

tab = curr 
for(i in 2:length(files)){
    cur = read.csv(files[i])

    cur$celltype = strsplit(files[i] , "--")[[1]][1]
    cur$gene = cur$X
    cur$region = strsplit(files[i] , "--")[[1]][2]
    cur$region = gsub(".csv","", cur$region)

    tab = rbind(tab, cur)
        
}

head(tab)
nrow(tab)

tab = tab[which(tab$p_val_adj<0.05),]

nrow(tab)


tab = tab[which(abs(tab$avg_log2FC)>0.1),]
tab = tab[which((tab$`pct.1`+tab$`pct.2`)/2 > 0.01),]


head(tab[which(tab$celltype=="Oligo_NN"),])

tab$direction = "Down"
tab$direction[which(tab$avg_log2FC<0)] = "Up"

ord

table(tab$celltype)

table(tab$celltype)

options(repr.plot.width=25, repr.plot.height=12)

d = tab
d$celltype = factor(d$celltype, levels = rev(ord))
d= d[which(!is.na(d$celltype)),]
d = d[-which(d$p_val_adj>0.01),]
g = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~region,ncol = 8)+ ggtitle(paste("MAST FDR < 0.01 All DEGs"))+
  theme(text = element_text(size = 18))
    print(g)


unique(tab$celltype)[which(!unique(tab$celltype) %in% ord)]

options(repr.plot.width=25, repr.plot.height=12)

d = tab
d$celltype = factor(d$celltype, levels = rev(ord))
d= d[which(!is.na(d$celltype)),]
d = d[-which(d$p_val_adj>0.01),]
g = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~region,ncol = 8)+ ggtitle(paste("MAST FDR < 0.01 All DEGs"))+
  theme(text = element_text(size = 18))
    print(g)


head(d$celltype)

options(repr.plot.width=25, repr.plot.height=12)
#mito ribo rep
d = tab
d$celltype = factor(d$celltype, levels = rev(ord))
d= d[which(!is.na(d$celltype)),]
#d = d[-which(d$p_val_adj>0.01),]
g = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~region,ncol = 8)+ ggtitle(paste("MAST FDR < 0.01 All DEGs"))+
  theme(text = element_text(size = 18))
    print(g)


options(repr.plot.width=25, repr.plot.height=12)

d = tab
d$celltype = factor(d$celltype, levels = rev(ord))
d= d[which(!is.na(d$celltype)),]
#d = d[-which(d$p_val_adj>0.01),]
g = ggplot(d, aes(x = celltype, ,fill = direction)) +
  geom_bar(stat = "count",position = "dodge", color = "black") +
  theme_minimal() +ylim(c(0,1000))+
  labs(title = "",
       x = " Type",
       y = "Number of DE") +
  coord_flip() + facet_wrap(~region,ncol = 8)+ ggtitle(paste("MAST FDR < 0.01 All DEGs"))+
  theme(text = element_text(size = 18))
    print(g)



