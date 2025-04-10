library(ggplot2)
library(Seurat)
library(reshape2) 


ameta = read.csv("/home/lamaral/ps-renlab2/projects/combined_all/final_meta.csv")

head(ameta)

ameta$celltype_final[grep("Astro-", ameta$celltype_final)] = "Astro NN"
ameta$celltype_final[grep("IOL", ameta$best_celltype_fixed)] = "IOL NN"
ameta$celltype_final[grep("RGL", ameta$best_celltype_fixed)] = "RGL NN"


ameta$age_rep = paste(ameta$age, ameta$rep,ameta$region, ameta$batch)

tab = table(ameta$celltype_final,ameta$age_rep)
#head(tab)
tab = sweep(tab,2,colSums(tab),'/')
#head(tab)
melted = melt(tab)
#head(melted)

melted$rep = sapply(strsplit(as.character(melted$Var2), " "), "[[", 2)
melted$age = sapply(strsplit(as.character(melted$Var2), " "), "[[", 1)
melted$region = sapply(strsplit(as.character(melted$Var2), " "), "[[", 3)
melted$age = factor(melted$age, levels = c("2mo", "9mo", "18mo"))
melted$sex = sapply(strsplit(as.character(melted$Var2), " "), "[[", 4)


options(repr.plot.width=10, repr.plot.height=6)

cl = "IOL NN"

pdf("cell_fraction.pdf")
for (cl in unique(melted$Var1)) {

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl) + facet_wrap(~sex, nrow = 2)

print(g2)
    }

dev.off()

options(repr.plot.width=10, repr.plot.height=6)

cl = "OB-STR-CTX Inh IMN"

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl) + facet_wrap(~sex, nrow = 2)

g2

options(repr.plot.width=10, repr.plot.height=6)

cl = "Astro NN"

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl) + facet_wrap(~sex, nrow = 2)

g2

obj = readRDS("RNA_final.RDS")

meta = obj@meta.data

melted$rep = sapply(strsplit(as.character(melted$Var2), " "), "[[", 2)
melted$age = sapply(strsplit(as.character(melted$Var2), " "), "[[", 1)
melted$region = sapply(strsplit(as.character(melted$Var2), " "), "[[", 3)
melted$age = factor(melted$age, levels = c("2mo", "9mo", "18mo"))


options(repr.plot.width=10, repr.plot.height=6)

cl = "IOL NN"

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl)

g2

meta$age_rep = paste(meta$age, meta$rep,meta$region)

tab = table(meta$celltype_final,meta$age_rep)
head(tab)
tab = sweep(tab,2,colSums(tab),'/')
head(tab)
melted = melt(tab)
head(melted)

head(melted)

melted$rep = sapply(strsplit(as.character(melted$Var2), " "), "[[", 2)
melted$age = sapply(strsplit(as.character(melted$Var2), " "), "[[", 1)
melted$region = sapply(strsplit(as.character(melted$Var2), " "), "[[", 3)

melted$age = factor(melted$age, levels = c("2mo", "9mo", "18mo"))

options(repr.plot.width=10, repr.plot.height=6)

cl = "Oligo NN"

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl)

g2

cl = "DG Glut"

pdf("cell_fraction.pdf", width = 7, height = )
for (cl in unique(melted$Var1)) {

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl)

print(g2)
    }
dev.off()



cl = "Astro NN"

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl)

g2


cl = "PAG Pou4f1 Ebf2 Glut"

g2 = ggplot(melted[which(melted$Var1 == cl ),]) + 
  geom_col(aes(factor(region),value,fill=age, color =rep ),position="dodge") + 
  xlab("Cell Types") + ylab("Fraction") + 
  guides(fill=guide_legend(title="Age Group"), color = guide_legend(title="Rep"))  +
  theme_bw() +ggtitle(cl)

g2


unique(melted$Var1)



