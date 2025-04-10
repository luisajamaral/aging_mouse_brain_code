library(data.table)
library(ggplot2)

setwd("../combined_DARs_redone/")

meta = read.csv("../final_meta.csv")
meta$celltype_final = gsub("/","-", meta$celltype_final)
meta$celltype_final = gsub("Lymphoid","Lymphoid NN", meta$celltype_final)
meta = meta[-which(meta$celltype_final=="doublet"),]
meta$celltype_final=paste(meta$celltype_final)
meta$celltype_final[which(meta$best_celltype_fixed=="IOL")] = "IOL NN"
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final))])


tab = fread("background_peaks_L2-3_IT_CTX_Glut_2vs18_ann.txt")
tab$ann = sapply(strsplit(as.character(tab$Annotation), " "), `[`, 1)
tab$type = "enhancer"
tab$type[which(tab$ann=="promoter-TSS")] = "promoter"

up = fread("up_peaks_L2-3_IT_CTX_Glut.ann.txt")
up$ann = sapply(strsplit(as.character(up$Annotation), " "), `[`, 1)
up$type = "enhancer"
up$type[which(up$ann=="promoter-TSS")] = "promoter"

down = fread("down_peaks_L2-3_IT_CTX_Glut.ann.txt")
down$ann = sapply(strsplit(as.character(down$Annotation), " "), `[`, 1)
down$type = "enhancer"
down$type[which(down$ann=="promoter-TSS")] = "promoter"

back_pt  = fread(paste("background_peaks_",ct,"_2vs18.csv", sep = ""))
head(back_pt)
cd$`feature name`[which(!cd$`feature name` %in% back_pt$background_peaks)]


setwd("../combined_DARs_redone/")
files = list.files(".", "diff_peaks")

all_df  = list()
for(f in files){
    cd = fread(f)
    cd$`log2(fold_change)` = cd$`log2(fold_change)`
    ct = gsub("diff_peaks_", "", f)
    ct = gsub("_2vs18.csv", "", ct)
    anno  = fread(paste("background_peaks_",ct,"_2vs18_ann.txt", sep = ""))
    anno$ann = sapply(strsplit(as.character(anno$Annotation), " "), `[`, 1)
    anno=as.data.frame(anno)
    rownames(anno) = paste(anno$Chr,":",anno$Start-1,"-",anno$End, sep = "")
    back_pt  = fread(paste("background_peaks_",ct,"_2vs18.csv", sep = ""))
    back_pt = as.data.frame(back_pt)
    rownames(back_pt) = back_pt[,1]
    cd$pct_2mo = back_pt[paste(cd$`feature name`), 'pct_cells_2mo']
    cd$pct_18mo = back_pt[paste(cd$`feature name`), 'pct_cells_18mo']
    cd$annotation = anno[paste(cd$`feature name`),"ann"]
    cat(ct)

    #cat(cd$`feature name`[which(! paste(cd$`feature name`) %in% rownames(anno))])
    cd$celltype = ct
    #head(cd)
    all_df[[ct]] = cd
}
    

diff = do.call(rbind, all_df)


head(diff)

diff$celltype = gsub("Lymphoid", "Lymphoid_NN" , diff$celltype)
diff$celltype = factor(diff$celltype , levels = gsub(" " , "_", ord))


diff$type = "non-TSS"
diff$type[which(diff$annotation=="promoter-TSS")]  = "promoter-TSS"

diff$major = sapply(strsplit(as.character(diff$celltype), "_"), function(x) tail(x, 1))
diff$major[which(diff$major%in%"Neur")] = "Gaba"
diff$major[which(diff$major%in%"Gaba-Chol")] = "Gaba"
diff$major[which(diff$major%in%"Glut-Sero")] = "Glut"
diff$major[which(diff$major%in%"Dopa")] = "Glut"
diff$major[which(diff$major%in%"IMN")] = "NN"


table(diff$major)

options(repr.plot.width=11, repr.plot.height=5)

ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.25 & diff$type == "non-TSS"),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Non-TSS Age DARs",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

options(repr.plot.width=11, repr.plot.height=5)

ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.2 & diff$type == "promoter-TSS"),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "TSS Age-DARs",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

options(repr.plot.width=12, repr.plot.height=6)

library(dplyr)
peak_counts <- diff[which(diff$type == "promoter-TSS"),] %>%
  group_by(celltype) %>%
  summarise(up = sum(`log2(fold_change)` > .25 & `adjusted p-value`<1),
            down = sum(`log2(fold_change)` < -.25 & `adjusted p-value`<1))

# Melt the data for easy plotting
melted_peak_counts <- melt(peak_counts, id.vars = "celltype")


melted_peak_counts$major = sapply(strsplit(as.character(melted_peak_counts$celltype), "_"), function(x) tail(x, 1))
melted_peak_counts$major[which(melted_peak_counts$major%in%"Neur")] = "Gaba"
melted_peak_counts$major[which(melted_peak_counts$major%in%"Gaba-Chol")] = "Gaba"
melted_peak_counts$major[which(melted_peak_counts$major%in%"Glut-Sero")] = "Glut"
melted_peak_counts$major[which(melted_peak_counts$major%in%"Dopa")] = "Glut"
melted_peak_counts$major[which(melted_peak_counts$major%in%"IMN")] = "NN"

table(melted_peak_counts$major)



# Create the bar plot
ggplot(melted_peak_counts, aes(x = celltype, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "promoter-TSS DARs",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

head(diff[which(diff$celltype=="L2-3_IT_CTX_Glut" & diff$annotation=="promoter-TSS"),])

options(repr.plot.width=12, repr.plot.height=6)

library(dplyr)
peak_counts <- diff[which(diff$type == "non-TSS"),] %>%
  group_by(celltype) %>%
  summarise(up = sum(`log2(fold_change)` > .25 & `adjusted p-value`<1 ),
            down = sum(`log2(fold_change)` < -.25 & `adjusted p-value`<1))

# Melt the data for easy plotting
melted_peak_counts <- melt(peak_counts, id.vars = "celltype")


melted_peak_counts$major = sapply(strsplit(as.character(melted_peak_counts$celltype), "_"), function(x) tail(x, 1))
melted_peak_counts$major[which(melted_peak_counts$major%in%"Neur")] = "Gaba"
melted_peak_counts$major[which(melted_peak_counts$major%in%"Gaba-Chol")] = "Gaba"
melted_peak_counts$major[which(melted_peak_counts$major%in%"Glut-Sero")] = "Glut"
melted_peak_counts$major[which(melted_peak_counts$major%in%"Dopa")] = "Glut"
melted_peak_counts$major[which(melted_peak_counts$major%in%"IMN")] = "NN"

table(melted_peak_counts$major)



# Create the bar plot
ggplot(melted_peak_counts, aes(x = celltype, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  labs(title = "non-TSS DARs",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

library(stats)

#op <- par(mfrow = c(3, 1), mgp = c(1.5, 0.8, 0), mar =  .1+c(3,3,2,1))

F10 <- ecdf(diff$`log2(fold_change)`[which(diff$celltype=="DG_Glut")])
F0 <- ecdf(diff$`log2(fold_change)`[which(diff$celltype=="Oligo_NN")])

summary(F10)
summary(F0)


head(diff)


head(diff)

diff$log_pct_FC = log2((diff$pct_18mo+0.001)/(diff$pct_2mo+0.001))

head(diff)

options(repr.plot.width=14, repr.plot.height=4.5)

ggplot(diff, aes(x = log_pct_FC, color = type
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log2(pct_fold_change)",
       x = "log2(pct_fold_change)",
       y = "eCDF",
       color = "Type") +xlim(c(-15,15))+
  theme_minimal()+ facet_wrap(~major)


diff$dir = "ns"
diff$dir[which(diff$`adjusted p-value`<0.05 & diff$`log2(fold_change)`>0)] = "up"
diff$dir[which(diff$`adjusted p-value`<0.05 & diff$`log2(fold_change)`<0)] = "down"


options(repr.plot.width=14, repr.plot.height=4.5)

ggplot(diff, aes(x = (pct_2mo+pct_18mo / 2), color = dir, linetype = type
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of %cells accessible at cCRE",
       x = "% of cells accessible at cCRE",
       y = "eCDF",
       color = "Change in aging",linetype = "cCRE type") +
  theme_minimal()+ facet_wrap(~major)


options(repr.plot.width=14, repr.plot.height=4.5)

ggplot(diff, aes(x = (pct_2mo+pct_18mo / 2), color = dir, linetype = type
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log2(pct_fold_change)",
       x = "% of cells accessible at cCRE",
       y = "eCDF",
       color = "Change in aging",linetype = "cCRE type") +
  theme_minimal()+ facet_wrap(~major)


options(repr.plot.width=9, repr.plot.height=4.5)

ggplot(diff, aes(x = major,y = (pct_2mo+pct_18mo / 2), color = dir, fill = type
                )) +
  geom_boxplot() +
  labs(title = "vln of % accessibility",
       y = "% of cells accessible at cCRE",
       x = "",
       color = "Change in aging",linetype = "cCRE type") +
      scale_color_manual(values = c("blue", "grey", "red")) +

  theme_minimal()


options(repr.plot.width=14, repr.plot.height=4.5)

ggplot(diff[which(diff$log_pct_FC>0),], aes(x = log_pct_FC, color = type
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log_pct_FC - Age-UP DARs",
       x = "log_pct_FC",
       y = "eCDF",
       color = "Type") +
  theme_minimal()+ facet_wrap(~major)


options(repr.plot.width=14, repr.plot.height=4.5)

ggplot(diff[which(diff$`log2(fold_change)`>0.1),], aes(x = `log2(fold_change)`, color = type
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log2(fold_change) - Age-UP DARs",
       x = "log2(fold_change)",
       y = "eCDF",
       color = "Type") +
  theme_minimal()+ facet_wrap(~major)


ggplot(diff[which(diff$`log2(fold_change)`<0),], aes(x = (`log2(fold_change)`), color = type
                )) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log2(fold_change) - Age-Down DARs",
       x = "log2(fold_change)",
       y = "eCDF",
       color = "Type") +
  theme_minimal()+ facet_wrap(~major)


ggplot(diff[which(diff$major=="Gaba"),], aes(x = `log2(fold_change)`, color = type)) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log2(fold_change) - Gaba",
       x = "log2(fold_change)",
       y = "eCDF",
       color = "Type") +
  theme_minimal()


ggplot(diff[which(diff$major=="Glut"),], aes(x = `log2(fold_change)`, color = type)) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log2(fold_change) - Glut",
       x = "log2(fold_change)",
       y = "eCDF",
       color = "Type") +
  theme_minimal()


ggplot(diff[which(diff$major=="NN"),], aes(x = `log2(fold_change)`, color = type)) +
  stat_ecdf(geom = "step") +
  labs(title = "eCDF of log2(fold_change) - NN",
       x = "log2(fold_change)",
       y = "eCDF",
       color = "Type") +
  theme_minimal()


head(diff)


library(tidyr)
library(dplyr)
data_bed <- diff %>%
  separate(`feature name`, into = c("chrom", "start_end"), sep = ":") %>%
  separate(start_end, into = c("start", "end"), sep = "-")




head(data_bed)

for (ct in unique(diff$celltype)) {
    up_TSS = diff[which(diff$celltype==ct & diff$`adjusted p-value`< 0.05 &diff$`log2(fold_change)`>0.25 & diff$type == "promoter-TSS"),]
    up_nonTSS = diff[which(diff$celltype==ct & diff$`adjusted p-value`< 0.05 & diff$`log2(fold_change)`> 0.25 & diff$type == "non-TSS"),]
    down_TSS = diff[which(diff$celltype==ct & diff$`adjusted p-value`< 0.05 & diff$`log2(fold_change)`< -0.25 & diff$type == "promoter-TSS"),]
    down_nonTSS = diff[which(diff$celltype==ct& diff$`adjusted p-value`< 0.05 & diff$`log2(fold_change)`< -0.25 & diff$type == "non-TSS"),]

    data_bed <- up_TSS %>%
      separate(`feature name`, into = c("chrom", "start_end"), sep = ":") %>%
      separate(start_end, into = c("start", "end"), sep = "-")
    write.table(data_bed[,c(1:3)], row.names = F, col.names = F, quote = F, sep = "\t", 
                file = paste("../split_promoters_enhancers/combined_diff_promoters_adjp_0.05/", ct,"--up.bed", sep = "")
               )
    data_bed <- down_TSS %>%
      separate(`feature name`, into = c("chrom", "start_end"), sep = ":") %>%
      separate(start_end, into = c("start", "end"), sep = "-")
    write.table(data_bed[,c(1:3)], row.names = F, col.names = F, quote = F, sep = "\t", 
                file = paste("../split_promoters_enhancers/combined_diff_promoters_adjp_0.05/", ct,"--down.bed", sep = "")
               )

    data_bed <- up_nonTSS %>%
      separate(`feature name`, into = c("chrom", "start_end"), sep = ":") %>%
      separate(start_end, into = c("start", "end"), sep = "-")
    write.table(data_bed[,c(1:3)], row.names = F, col.names = F, quote = F, sep = "\t", 
                file = paste("../split_promoters_enhancers/combined_diff_enhancers_adjp_0.05/", ct,"--up.bed", sep = "")
               )

    data_bed <- down_nonTSS %>%
      separate(`feature name`, into = c("chrom", "start_end"), sep = ":") %>%
      separate(start_end, into = c("start", "end"), sep = "-")
    write.table(data_bed[,c(1:3)], row.names = F, col.names = F, quote = F, sep = "\t", 
                file = paste("../split_promoters_enhancers/combined_diff_enhancers_adjp_0.05/", ct,"--down.bed", sep = "")
               )
    }
    

files[1]

files = list.files(".", "background_peaks_")
files = files[grep("_2vs18_ann.txt", files)]
for(f in files){
    ct = gsub("_2vs18_ann.txt", "", f)
    ct = gsub("background_peaks_", "", ct)

    anno = fread(f)
    anno$ann = sapply(strsplit(as.character(anno$Annotation), " "), `[`, 1)
    anno=as.data.frame(anno)
    anno$Start = anno$Start-1
    back_tss=anno[which(anno$ann == "promoter-TSS"), c(2,3,4)]
    write.table(back_tss, row.names = F, col.names = F, quote = F, sep = "\t", 
                file = paste("../split_promoters_enhancers/background_promoters/", ct,".bed", sep = "")
               )

    back_tss=anno[which(anno$ann != "promoter-TSS"), c(2,3,4)]
    write.table(back_tss, row.names = F, col.names = F, quote = F, sep = "\t", 
                file = paste("../split_promoters_enhancers/background_enhancers/", ct,".bed", sep = "")
               )
            
}

head(anno)




anno  = fread(paste("background_peaks_",ct,"_2vs18_ann.txt", sep = ""))
    anno$ann = sapply(strsplit(as.character(anno$Annotation), " "), `[`, 1)
    anno=as.data.frame(anno)
    rownames(anno) = paste(anno$Chr,":",anno$Start-1,"-",anno$End, sep = "")
    cd$annotation = anno[paste(cd$`feature name`),"ann"]
    cat(ct)

   
