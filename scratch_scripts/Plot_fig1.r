library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(ggrepel)

setwd("../h5ads_final/")

setwd("../../h5ads_final/")

options(repr.plot.width=4, repr.plot.height=10)


meta = read.csv("../final_meta.csv")
meta$celltype_final = gsub("/","-", meta$celltype_final)
meta$celltype_final = gsub("Lymphoid","Lymphoid NN", meta$celltype_final)
meta = meta[-which(meta$celltype_final=="doublet"),]
meta$celltype_final=paste(meta$celltype_final)
meta$celltype_final[which(meta$best_celltype_fixed=="IOL")] = "IOL NN"
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)

# Create a vertical bar plot
ggplot(meta, aes(x = celltype_final)) +
  geom_bar(stat = "count", color = "black", fill = "tomato") +
  theme_minimal() +
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()

meta$major = sapply(strsplit(as.character(meta$celltype_final), " "), function(x) tail(x, 1))
meta$major[which(meta$major%in%"Neur")] = "Gaba"
meta$major[which(meta$major%in%"Gaba-Chol")] = "Gaba"
meta$major[which(meta$major%in%"Glut-Sero")] = "Glut"
meta$major[which(meta$major%in%"Dopa")] = "Glut"
meta$major[which(meta$major%in%"IMN")] = "NN"



table(meta$major)

ord = names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)
ggplot(meta, aes(x = celltype_final)) +
  geom_bar(stat = "count", color = "black", fill = "tomato") +
  theme_minimal() +
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()

#key[which( key$CellType %in%meta$celltype_final ),]
setwd("~/projects/combined_all/Figures/Figure3/")
key_file = "../color_scheme/updated_celltype_palette.csv"
key = read.csv(key_file)
head(key)
color_vector <- setNames(key$Color, key$CellType)
unique(meta$celltype_final)[which(! unique(meta$celltype_final) %in% key$CellType  )]

options(repr.plot.width=15, repr.plot.height=7)

ggplot(meta, aes(x = celltype_final)) +
  geom_bar(stat = "count", color = "black", fill = "tomato") +
  theme_minimal() +
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +coord_flip()+
 facet_wrap(~major, scales="free") + scale_fill_manual(color_vector)

library(scales)
options(repr.plot.width=7, repr.plot.height=6)


m = meta[-grep("NN", meta$celltype_final),]
m = m[-grep("IMN", m$celltype_final),]

g <- ggplot(m, aes(x = celltype_final, fill = celltype_final)) +
  geom_bar(stat = "count", color = "black") +
  scale_fill_manual(values = color_vector) +
  theme_minimal() + 
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip() + 
  geom_text(stat='count', aes(label=label_number(scale_cut = cut_si(""), accuracy = 1)(after_stat(count))), 
            position=position_stack(vjust=1), hjust=-0.1) + 
  scale_y_continuous(breaks = c(0, 50000, 100000, 150000, 200000), 
                     labels = label_number(scale_cut = cut_si(""), accuracy = 1),
                     limits = c(0, 225000),  # Adjusts the upper limit to 110% of the max value
                     expand = expansion(mult = c(0, 0.05))) +  # Adds some space at the top
  facet_wrap(~major, scales = "free") +
  theme(legend.position = "none")  # Remove the legend entirely

# Display the plot
print(g)
svg("../Figures/Figure3/Neuron_cell_count.svg", height = 6, width = 7)
print(g)
dev.off()


pdf("celltype_count_bar.pdf", height = 10, width = 4)
ggplot(meta, aes(x = celltype_final, fill = )) +
  geom_bar(stat = "count", color = "black", fill = "tomato") +
  theme_minimal() +
  labs(title = "Cells per Cell Type",
       x = "Cell Type",
       y = "Number of Cells") +
  coord_flip()
dev.off()

options(repr.plot.width=5, repr.plot.height=10)
ggplot(meta, aes(x = celltype_final, fill = region)) +
  geom_bar(stat = "count", color = "black",position = "fill" ) +
  theme_minimal() +
  labs(title = "Cells per region",
       x = "Cell Type",
       y = "Fraction of Cells") +
  coord_flip()

pdf("celltype_region_bar.pdf", height = 10, width = 5)
ggplot(meta, aes(x = celltype_final, fill = region)) +
  geom_bar(stat = "count", color = "black",position = "fill" ) +
  theme_minimal() +
  labs(title = "Cells per region",
       x = "Cell Type",
       y = "Fraction of Cells") +
  coord_flip()
dev.off()

peaks = fread("peaks_0.01.csv" )
peaks = as.data.frame(peaks)
rownames(peaks) = peaks$Peaks

sums <- colSums(peaks[,-c(1)])
df <- data.frame(feature = c(names(sums),'CB Granule Glut' ), value = c(sums,0))

# Order the data frame by column sums
#df = rbind(df, c('CB Granule Glut', 0 ))

df <- df[order(df$value, decreasing = TRUE),]
df$feature = gsub("Lymphoid", "Lymphoid NN", df$feature)
df$feature = factor(df$feature,levels = ord)

# Create a vertical bar plot
pdf("celltype_peaks_bar.pdf", height = 10, width = 4)

ggplot(df, aes(x = feature, y = value)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Peaks called",
       x = "Cell type",
       y = "Number of Peaks")+
  coord_flip()
dev.off()

options(repr.plot.width=5, repr.plot.height=10)

ggplot(df, aes(x = feature, y = value)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "Peaks called",
       x = "Cell type",
       y = "Number of Peaks")+
  coord_flip()

#write background bed file
#for(i in 2:ncol(peaks)) {
#    rn = peaks[which(peaks[,i] == TRUE),]
#    rn = separate(rn, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
#    ct = gsub(" ", "_", colnames(peaks)[i])
#    write.table(rn[,c(1:3)] , file = paste(ct, "_allpeaks.bed", sep = ""),row.names = F, sep = "\t" , quote = F) 
#}

head(meta)

options(repr.plot.width=5, repr.plot.height=10)
colors =c(brewer.pal(3, "RdPu")[2], brewer.pal(3, "PuBu")[2])
pdf("celltype_sex_bar.pdf", height = 10, width = 4.5)
#colors <- brewer.pal(3, "Purples")
ggplot(meta, aes(x = celltype_final, fill = batch)) +
  geom_bar(stat = "count", position = "fill", color = "black") +
  theme_minimal() +
  scale_fill_manual(values = colors) + 
  labs(title = "Cells per Sex",
       x = "Cell Type",
       y = "% of Cells") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 100))
dev.off()




options(repr.plot.width=7, repr.plot.height=12)

meta$age = factor(meta$age , levels = c("2mo", "9mo", "18mo"))
# Create a vertical bar plot
pdf("celltype_age_bar.pdf", height = 10.5, width = 6)

colors <- brewer.pal(3, "Blues")

ggplot(meta, aes(x = celltype_final, fill = age)) +
  geom_bar(stat = "count", position = "fill", color = "black") +
  theme_minimal() +
  labs(
       x = "Cell type",
       y = "% of cells") +
  scale_fill_manual(values = colors) + 
  coord_flip() +
  facet_wrap(~ batch) +
  scale_y_continuous(labels = scales::percent_format(scale = 100))
dev.off()

display.brewer.all()

colors <- brewer.pal(3, "Purples")

ggplot(meta, aes(x = celltype_final, fill = age)) +
  geom_bar(stat = "count", position = "fill", color = "black") +
  theme_minimal() +
  labs(title = "Cells per Age",
       x = "Cell type",
       y = "% of cells") +
  scale_fill_manual(values = colors) + 
  coord_flip() +
  facet_wrap(~ batch) +
  scale_y_continuous(labels = scales::percent_format(scale = 100))

files


setwd("../../region_DARs_redone//")
files = list.files(".", "diff_peaks")

all_df  = list()
for(f in files){
    cd = fread(f)
    cd$`log2(fold_change)` = cd$`log2(fold_change)`
    #up = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`>0),"feature name"]
    #down = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`<0),"feature name"]

    #up <- separate(up, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
    #down <- separate(down, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
        
    #fwrite(down, file = gsub("_2vs18_iter.csv" , "_down.bed" ,f),sep = "\t", col.names = F, row.names = F)
    #fwrite(up, file = gsub("_2vs18_iter.csv" , "_up.bed" ,f),sep = "\t", col.names = F, row.names = F)

    ct = gsub("diff_peaks_", "", f)
    ct = gsub("_2vs18.csv", "", ct)
    bed = cd[which(cd$`adjusted p-value`<0.05),"feature name"]
    cd$celltype = ct
    all_df[[ct]] = cd
}
    

diff = do.call(rbind, all_df)
head(diff)

head(diff)


diff$region = sapply(strsplit(as.character(diff$celltype), ":"), function(x) tail(x, 1))
                     
diff$celltype = sapply(strsplit(as.character(diff$celltype), ":"),`[`, 1)
diff$celltype = gsub("Lymphoid", "Lymphoid_NN" , diff$celltype)
diff$celltype = factor(diff$celltype , levels = gsub(" " , "_", ord))


diff$major = sapply(strsplit(as.character(diff$celltype), "_"), function(x) tail(x, 1))

diff$major[which(diff$major%in%"Neur")] = "Gaba"
diff$major[which(diff$major%in%"Gaba-Chol")] = "Gaba"
diff$major[which(diff$major%in%"Glut-Sero")] = "Glut"
diff$major[which(diff$major%in%"Dopa")] = "Glut"
diff$major[which(diff$major%in%"IMN")] = "IMN"
head(diff)



options(repr.plot.width=11, repr.plot.height=5)

for(reg in unique(diff$region)){
  g= ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.25 & diff$region ==reg),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = paste("Age DARs", reg),
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip() + facet_grid(~region, scales = "free")
  print(g)
}

print(g)


options(repr.plot.width=14, repr.plot.height=7)

ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.25),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Accessible Regions",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +ylim(c(-8,8))+
  theme_minimal() + coord_flip() + facet_wrap(~region, scales = "free",nrow=2)

nn = diff[which(!diff$major%in%  c("NN", "IMN")),]
table(nn$celltype)


options(repr.plot.width=12, repr.plot.height=5)

g = ggplot(nn[which(abs(nn$`log2(fold_change)`) > 0.5),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Accessible Regions",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +ylim(c(-8,8))+
  theme_minimal() + coord_flip() + facet_wrap(~region, scales = "free",nrow=2)



pdf("~/projects/combined_all/Figures/Figure3/neuron_DAR_plot.pdf", height =8, width = 13)
print(g)
dev.off()

options(repr.plot.width=14, repr.plot.height=7)

ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.25),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Accessible Regions",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +ylim(c(-8,8))+
  theme_minimal() + coord_flip() + facet_wrap(~region, scales = "free",nrow=2)


setwd("../combined_DARs_redone/")
files = list.files(".", "diff_peaks")

all_df  = list()
for(f in files){
    cd = fread(f)
    cd$`log2(fold_change)` = cd$`log2(fold_change)`
    #up = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`>0),"feature name"]
    #down = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`<0),"feature name"]

    #up <- separate(up, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
    #down <- separate(down, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
        
    #fwrite(down, file = gsub("_2vs18_iter.csv" , "_down.bed" ,f),sep = "\t", col.names = F, row.names = F)
    #fwrite(up, file = gsub("_2vs18_iter.csv" , "_up.bed" ,f),sep = "\t", col.names = F, row.names = F)

    ct = gsub("diff_peaks_", "", f)
    ct = gsub("_2vs18.csv", "", ct)
    bed = cd[which(cd$`adjusted p-value`<0.05),"feature name"]
    cd$celltype = ct
    all_df[[ct]] = cd
}
    

diff = do.call(rbind, all_df)
diff$celltype = gsub(":Male", "" , diff$celltype)

diff$celltype = gsub("Lymphoid", "Lymphoid_NN" , diff$celltype)
diff$celltype = factor(diff$celltype , levels = gsub(" " , "_", ord))
options(repr.plot.width=5.5, repr.plot.height=10)

options(repr.plot.width=11, repr.plot.height=5)

diff$major = sapply(strsplit(as.character(diff$celltype), "_"), function(x) tail(x, 1))
diff$major[which(diff$major%in%"Neur")] = "Gaba"
diff$major[which(diff$major%in%"Gaba-Chol")] = "Gaba"
diff$major[which(diff$major%in%"Glut-Sero")] = "Glut"
diff$major[which(diff$major%in%"Dopa")] = "Glut"
diff$major[which(diff$major%in%"IMN")] = "NN"


table(diff$major)
                    

ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.25),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Accessible Regions",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)


head(diff)

options(repr.plot.width=12, repr.plot.height=6)

library(dplyr)
peak_counts <- diff %>%
  group_by(celltype) %>%
  summarise(up = sum(`log2(fold_change)` > .15 & `adjusted p-value`<0.05 ),
            down = sum(`log2(fold_change)` < -.15 & `adjusted p-value`<0.05))

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
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Total Number of Differential Peaks",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

options(repr.plot.width=6
        , repr.plot.height=10)
ggplot(melted_peak_counts, aes(x = celltype, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Total Number of Differential Peaks",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() 


setwd("../../male_diff/diff_csvs/")
files = list.files(".", "2vs18")

all_df  = list()
for(f in files){
    cd = fread(f)
    cd$`log2(fold_change)` = -cd$`log2(fold_change)`
    #up = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`>0),"feature name"]
    #down = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`<0),"feature name"]

    #up <- separate(up, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
    #down <- separate(down, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
        
    #fwrite(down, file = gsub("_2vs18_iter.csv" , "_down.bed" ,f),sep = "\t", col.names = F, row.names = F)
    #fwrite(up, file = gsub("_2vs18_iter.csv" , "_up.bed" ,f),sep = "\t", col.names = F, row.names = F)

    ct = gsub("diff_peaks_", "", f)
    ct = gsub("_2vs18_iter.csv", "", ct)
    bed = cd[which(cd$`adjusted p-value`<0.05),"feature name"]
    cd$celltype = ct
    all_df[[ct]] = cd
}
    

diff = do.call(rbind, all_df)

head(diff)

diff$celltype = gsub(":Male", "" , diff$celltype)

diff$celltype = gsub("Lymphoid", "Lymphoid_NN" , diff$celltype)
diff$celltype = factor(diff$celltype , levels = gsub(" " , "_", ord))

head(diff)

options(repr.plot.width=5.5, repr.plot.height=10)

ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.15),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Accessible Regions",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip()

options(repr.plot.width=11, repr.plot.height=5)

diff$major = sapply(strsplit(as.character(diff$celltype), "_"), function(x) tail(x, 1))
diff$major[which(diff$major%in%"Neur")] = "Gaba"
diff$major[which(diff$major%in%"Gaba-Chol")] = "Gaba"
diff$major[which(diff$major%in%"Glut-Sero")] = "Glut"
diff$major[which(diff$major%in%"Dopa")] = "Glut"
diff$major[which(diff$major%in%"IMN")] = "NN"


table(diff$major)
                    

ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.25),], 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Accessible Regions",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

options(repr.plot.width=12, repr.plot.height=6)

library(dplyr)
peak_counts <- diff %>%
  group_by(celltype) %>%
  summarise(up = sum(`log2(fold_change)` > .15 & `adjusted p-value`<0.05 ),
            down = sum(`log2(fold_change)` < -.15 & `adjusted p-value`<0.05))

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
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Total Number of Differential Peaks",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

setwd("../../male_diff/diff_csvs/")
files = list.files(".", "2vs18")

all_df  = list()
for(f in files){
    cat(f)
    cd = fread(f)
    cd$`log2(fold_change)` = -cd$`log2(fold_change)`
   # up = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`>0),"feature name"]
   # down = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`<0),"feature name"]

   # up <- separate(up, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
   # down <- separate(down, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
        
   # fwrite(down, file = gsub("_2vs18_iter.csv" , "_down.bed" ,f),sep = "\t", col.names = F, row.names = F)
   # fwrite(up, file = gsub("_2vs18_iter.csv" , "_up.bed" ,f),sep = "\t", col.names = F, row.names = F)

    ct = gsub("diff_peaks_", "", f)
    ct = gsub("_2vs18_iter.csv", "", ct)
    bed = cd[which(cd$`adjusted p-value`<0.05),"feature name"]
    cd$celltype = ct
    all_df[[ct]] = cd
}
male_diff = do.call(rbind, all_df)

setwd("../../female_diff/diff_csvs/")
files = list.files(".", "2vs18")

all_df  = list()
for(f in files){
    cat(f)
    cd = fread(f)
    cd$`log2(fold_change)` = -cd$`log2(fold_change)`
   # up = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`>0),"feature name"]
   # down = cd[which(cd$`adjusted p-value`<0.05 & cd$`log2(fold_change)`<0),"feature name"]

   # up <- separate(up, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
   # down <- separate(down, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
        
   # fwrite(down, file = gsub("_2vs18_iter.csv" , "_down.bed" ,f),sep = "\t", col.names = F, row.names = F)
   # fwrite(up, file = gsub("_2vs18_iter.csv" , "_up.bed" ,f),sep = "\t", col.names = F, row.names = F)

    ct = gsub("diff_peaks_", "", f)
    ct = gsub("_2vs18_iter.csv", "", ct)
    bed = cd[which(cd$`adjusted p-value`<0.05),"feature name"]
    cd$celltype = ct
    all_df[[ct]] = cd
}
female_diff = do.call(rbind, all_df)

female_diff$ct = gsub(":Female", "", female_diff$celltype)
female_diff$ct_dir = paste(female_diff$ct,"Up", sep = "_")
female_diff$ct_dir[which(female_diff$`log2(fold_change)`<0)] = paste(female_diff$ct,"Down", sep = "_")
female_diff$ct_dir_loc = paste(female_diff$ct_dir, female_diff$`feature name`, sep = "_")

male_diff$ct = gsub(":Male", "", male_diff$celltype)
male_diff$ct_dir = paste(male_diff$ct,"Up", sep = "_")
male_diff$ct_dir[which(male_diff$`log2(fold_change)`<0)] = paste(male_diff$ct,"Down", sep = "_")
male_diff$ct_dir_loc = paste(male_diff$ct_dir, male_diff$`feature name`, sep = "_")

female_diff = female_diff[which(female_diff$`adjusted p-value`<0.05),]
male_diff = male_diff[which(male_diff$`adjusted p-value`<0.05),]
#diff = diff[which(diff$`adjusted p-value`<0.05),]

diff = diff[,-c(7)]

diff$ct_dir = paste(diff$celltype,"Up", sep = "_")
diff$ct_dir[which(diff$`log2(fold_change)`<0)] = paste(diff$celltype,"Down", sep = "_")
diff$ct_dir_loc = paste(diff$ct_dir, diff$`feature name`, sep = "_")

diff$sig = "Sig in combined"
diff$sig[which(diff$ct_dir_loc%in%female_diff$ct_dir_loc)] = "Sig in combined + Female"
diff$sig[which(diff$ct_dir_loc%in%male_diff$ct_dir_loc)] = "Sig in combined + Male"
diff$sig[which(diff$ct_dir_loc%in%female_diff$ct_dir_loc & diff$ct_dir_loc%in%male_diff$ct_dir_loc )] = "Sig in combined + Male + Female"


options(repr.plot.width=9, repr.plot.height=12)

ggplot(diff[which(abs(diff$`log2(fold_change)`) > 0.15),], aes(x = celltype, y = `log2(fold_change)`, color = sig)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 3) +
  #scale_color_manual(values = c("grey", "red")) +
  labs(title = "Differential Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_minimal() + coord_flip()



female_diff$sex = "Female"
male_diff$sex = "Male"

bdiff = rbind(female_diff, male_diff)

peak_counts <- bdiff %>%
  group_by(ct, sex) %>%
  summarise(up = sum(`log2(fold_change)` > .15 & `adjusted p-value`<0.05 ),
            down = sum(`log2(fold_change)` < -.15 & `adjusted p-value`<0.05))

head(peak_counts)
melted_peak_counts <- melt(peak_counts)
head(melted_peak_counts)

options(repr.plot.width=8, repr.plot.height=10)

bdiff$ct = factor(bdiff$ct , levels = gsub(" " , "_", ord))
peak_counts <- bdiff %>%
  group_by(ct, sex) %>%
  summarise(up = sum(`log2(fold_change)` > .15 & `adjusted p-value`<0.05 ),
            down = sum(`log2(fold_change)` < -.15 & `adjusted p-value`<0.05))

# Melt the data for easy plotting
melted_peak_counts <- melt(peak_counts)


# Create the bar plot
ggplot(melted_peak_counts, aes(x = ct, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Total Number of Differential Peaks",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() + facet_wrap(~sex, ncol = 2)



setwd("../../female_RNA/DEG_results_latent_rep_mito_together/")


setwd("../../female_RNA/DEG_results_latent_rep_mito_together_wilcox_.01pct_0logfc//")


setwd("../../female_RNA/DEG_results_latent_orig.indent//")


#write background bed file
#for(i in 2:ncol(peaks)) {
#    rn = peaks[which(peaks[,i] == TRUE),]
#    rn = separate(rn, 1, into = c("chromosome", "start", "end"), sep = "[:-]")
#    ct = gsub(" ", "_", colnames(peaks)[i])
#    write.table(rn[,c(1:3)] , file = paste(ct, "_allpeaks.bed", sep = ""),row.names = F, sep = "\t" , quote = F) 
#}

#write.table(file = paste("up/", ct, ".txt", sep = ""), cd[which(cd$avg_log2FC>0& cd$p_val_adj<0.05),"V1"], quote = F, row.names = F, col.names = T)

files = list.files(".", ".csv")

all_df  = list()
for(f in files){
    cd = fread(f)
    cd$avg_log2FC = -cd$avg_log2FC
    ct = gsub(".csv", "", f)
    cd$celltype = ct
    all_df[[ct]] = cd
    #write.table(file = paste("up/", ct, ".txt", sep = ""), cd[which(cd$avg_log2FC>0& cd$p_val_adj<0.05),"V1"], quote = F, row.names = F, col.names = T)
    #write.table(file = paste("down/", ct, ".txt", sep = ""), cd[which(cd$avg_log2FC<0& cd$p_val_adj<0.05),"V1"], quote = F, row.names = F, col.names = T)
    #write.table(file = paste("bg/", ct, ".txt", sep = ""), cd[,"V1"], quote = F, row.names = F, col.names = T)
}
diff = do.call(rbind, all_df)

diff=diff[-which(is.na(diff$celltype)),]

unique(diff$celltype)[which(!unique(diff$celltype) %in%  gsub(" " , "_", ord))]

#diff = diff[-which(diff$celltype=="doublet"),]
diff$celltype = gsub("Lymphoid", "Lymphoid_NN",diff$celltype) 
diff$celltype = factor(diff$celltype , levels = gsub(" " , "_", ord))

diff$avg_pct = (diff$`pct.1`+diff$`pct.2`)/2

diff$label = diff$V1
diff$label[which(abs(diff$avg_log2FC)<2 | diff$p_val_adj>0.05)] = ""
head(diff)


options(repr.plot.width=12, repr.plot.height=6)


diff$major = sapply(strsplit(as.character(diff$celltype), "_"), function(x) tail(x, 1))
diff$major[which(diff$major%in%"Neur")] = "Gaba"
diff$major[which(diff$major%in%"Gaba-Chol")] = "Gaba"
diff$major[which(diff$major%in%"Glut-Sero")] = "Glut"
diff$major[which(diff$major%in%"Dopa")] = "Glut"
diff$major[which(diff$major%in%"IMN")] = "NN"
diff$major = factor(diff$major, levels = c("Glut", "Gaba", "NN"))
                    
                    

ggplot(diff[which(abs(diff$avg_log2FC) > 0.15 & diff$avg_pct>0.01),], 
  aes(x = celltype, y = avg_log2FC,
  color = ifelse(p_val_adj < 0.05, ifelse(avg_log2FC > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_classic() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

pdf("DEG_logfc_plot.pdf", width = 12, height = 6)
ggplot(diff[which(abs(diff$avg_log2FC) > 0.15 & diff$avg_pct>0.1),], 
  aes(x = celltype, y = avg_log2FC,
  color = ifelse(p_val_adj < 0.05, ifelse(avg_log2FC > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_classic() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

dev.off()

selected_points <- diff[which(diff$label != ""), ]

# Create the ggplot without geom_text_repel
p <- ggplot(diff[which(abs(diff$avg_log2FC) > 0.15 & diff$avg_pct > 0.1), ], 
            aes(x = celltype, y = avg_log2FC,
                color = ifelse(p_val_adj < 0.05, ifelse(avg_log2FC > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_classic() + ylim(c(-8, 8)) +
  coord_flip() + 
  facet_wrap(~major, scales = "free", ncol = 3)

# Check if the base plot renders correctly
print(p)

# Attempt to add labels for selected points using geom_text_repel with minimal parameters
#p + geom_text_repel(data = selected_points, aes(label = label))

dev.off()

selected_points <- diff[which(diff$label !=""), ]

# Create the ggplot with selected points
ggplot(diff[which(abs(diff$avg_log2FC) > 0.15 & diff$avg_pct>0.1),], aes(x = celltype, y = avg_log2FC,
                 color = ifelse(p_val_adj < 0.05, ifelse(avg_log2FC > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_classic() + ylim(c(-10,10))+
  coord_flip() + 
  facet_wrap(~major, scales = "free", ncol = 3) +

  # Add labels for selected points
  geom_text(data = selected_points, aes(label = label), 
            hjust = -0.2, vjust = -0.2, size = 3, color = "black", fontface = "bold")

selected_points <- diff[which(diff$label !=""), ]

ggplot(diff[which(abs(diff$avg_log2FC) > 0.15 & diff$avg_pct > 0.1), ], 
       aes(x = celltype, y = avg_log2FC,
           color = ifelse(p_val_adj < 0.05, ifelse(avg_log2FC > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Expression",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +
  theme_classic() + 
  ylim(c(-10, 10)) +
  coord_flip() + 
  facet_wrap(~major, scales = "free", ncol = 3) +

  # Add labels for selected points
  geom_text(data = selected_points, aes(label = label), 
            nudge_x = 0.5, nudge_y = 0.5, size = 2, color = "black", fontface = "bold")

head(diff)
genes  = diff$V1

require(biomaRt)
ensembl <- useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
#genes = meta$geneID
annot <- getBM(
  attributes = c(
    'mgi_symbol',
    'external_gene_name',
    'ensembl_gene_id',
    'entrezgene_id',
    'gene_biotype','transcript_length'),
  mart = ensembl)
head(annot)
#annot <- merge(
#  x = as.data.frame(genes),
#  y =  annot,
#  by.y = 'external_gene_name',
#  all.x = T,
#  by.x = 'genes')

head(as.data.frame(genes))

head(annot)

listFilters(ensembl)

mgi

mgi <- useMart('ENSEMBL_MART_MOUSE')

# Retrieve gene information
annot <- getBM(
  attributes = c(
    'mgi_symbol',
    'entrezgene_id',
    'gene_biotype',
    'transcript_length'
  ),
  filters = 'mgi_symbol',
  values = genes,
  mart = mgi
)

# Merge the data with 'diff'
annot <- merge(
  x = data.frame(genes),
  y = annot,
  by.x = 'genes',
  by.y = 'mgi_symbol',
  all.x = TRUE
)

listMarts()

listDatasets(ensembl)$dataset

unique(annot$transcript_length)

peak_counts <- diff %>%
  group_by(celltype) %>%
  summarise(up = sum(avg_log2FC > .25 & p_val_adj<0.01 & avg_pct >0.001),
            down = sum(avg_log2FC < -.25 & p_val_adj<0.01  & avg_pct >0.001))

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
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Total Number of Differential Genes",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

peak_counts <- diff %>%
  group_by(celltype) %>%
  summarise(up = sum(avg_log2FC > .25 & p_val_adj<0.01 & avg_pct >0.1),
            down = sum(avg_log2FC < -.25 & p_val_adj<0.01  & avg_pct >0.1))

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
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Total Number of Differential Genes",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)

peak_counts <- diff %>%
  group_by(celltype) %>%
  summarise(up = sum(avg_log2FC > 0 & p_val_adj<0.05 ),
            down = sum(avg_log2FC < -0 & p_val_adj<0.05))

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
  geom_bar(stat = "identity", position = "stack", color = "black") +
  labs(title = "Total Number of Differential Genes",
       x = "Cell Type",
       y = "Count",
       fill = "Direction") +
  theme_minimal() + coord_flip() + facet_wrap(~major, scales = "free", ncol = 3)


d = diff[which(diff$V1%in%target_genes),]


library(ggrepel)

unique_celltypes

options(repr.plot.width=12, repr.plot.height=6)

df = diff
df = as.data.frame(df)

target_genes <- c("Xist", "Tsix")

pdf("Xist_volcanos.pdf")
# Loop through each unique cell type
unique_celltypes <- unique(df$celltype)
for (celltype in unique_celltypes) {
  
  # Subset data for the specific cell type
  subset_df <- df[which(df$celltype == celltype),]
  
  # Create a volcano plot with label repulsion
  volcano_plot <- ggplot(subset_df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = ifelse(V1 %in% target_genes, "chartreuse1", 
                                  ifelse(p_val_adj < 0.05, ifelse(avg_log2FC > 0, "red", "blue"), "black")))) +
    geom_text_repel(
      data = subset_df %>% filter(V1 %in% target_genes),
      aes(label = V1),
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    scale_color_identity() +
    labs(title = paste("Volcano Plot for", celltype),
         x = "Average log2 Fold Change",
         y = "-log10(Adjusted p-value)") +
    theme_minimal()
  
  # Print or save the plot as needed
  print(volcano_plot)
}
dev.off()

dev.off()

df = diff
df = as.data.frame(df)

head(df)

pdf("volcanos.pdf")

unique_celltypes <- unique(df$celltype)
for (celltype in unique_celltypes) {
  
  # Subset data for the specific cell type
  subset_df <- df[which(df$celltype == celltype),]
  
  # Sort data by adjusted p-values
  subset_df <- subset_df %>% arrange(p_val_adj)
  
  # Select the top 20 genes for labeling
  top_genes <- head(subset_df$V1, 20)
  
  # Create a volcano plot with label repulsion
  volcano_plot <- ggplot(subset_df, aes(x = avg_log2FC, y = -log10(p_val_adj+1e-200))) +
    geom_point(aes(color = ifelse(p_val_adj < 0.05, ifelse(avg_log2FC > 0, "red", "blue"), "black"))) +
    geom_text_repel(max.overlaps = 35,
      data = subset_df %>% filter(V1 %in% top_genes),
      aes(label = V1),
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    scale_color_identity() +
    labs(title = paste("Volcano Plot for", celltype),
         x = "Average log2 Fold Change",
         y = "-log10(Adjusted p-value)") +
    theme_minimal()
  
  # Print or save the plot as needed
# pdf(paste(celltype, "_volcano.pdf", sep = ""))
  print(volcano_plot)
 # dev.off()
}
dev.off()

df

pdf("MAplots.pdf")

unique_celltypes <- unique(df$celltype)
for (celltype in unique_celltypes) {
  
  # Subset data for the specific cell type
  subset_df <- df[which(df$celltype == celltype),]
  
  labeled_genes <- head(subset_df$V1[subset_df$p_val_adj < 0.05], 50)
  
  # Create an MA plot
  ma_plot <- ggplot(subset_df, aes(x = (pct.1 + pct.2) / 2, y = avg_log2FC)) +
     geom_point(
      aes(color = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.15,
                         ifelse(avg_log2FC > 0, "red", "blue"), "black"))
        ) +    geom_text_repel(
      data = subset_df %>% filter(V1 %in% labeled_genes),
      aes(label = V1),
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    scale_color_identity() +
    labs(title = paste("MA Plot for", celltype),
         x = "Average Abundance (A)",
         y = "Log2 Fold Change (M)") +
    theme_minimal()
  
  # Print or save the plot as needed
  print(ma_plot)
}
dev.off()

while(3>1){
dev.off()}

rmeta = read.csv("../meta_celltypes.csv")

big_diff = diff

downdiff[which(downdiff$V1=="Xist"),"celltype"]

updiff[which(updiff$V1=="Xist"),"celltype"]

options(repr.plot.width=6, repr.plot.height=6)

downdiff = diff[which(diff$avg_log2FC<(-.15) & diff$p_val_adj<0.05),]
down_tab = table(downdiff$V1)[order(table(downdiff$V1) , decreasing = T)]
hist(down_tab, xlab = "Number of celltypes a gene is DE in")
head(table(downdiff$V1)[order(table(downdiff$V1) , decreasing = T)],50)

down_tab[1:100]

df = as.data.frame(head(table(downdiff$V1)[order(table(downdiff$V1) , decreasing = T)],20))
colnames(df) = c("Gene","Freq")
df

updiff = diff[which(diff$avg_log2FC>.1& diff$p_val_adj<0.05),]
head(table(updiff$V1)[order(table(updiff$V1) , decreasing = T)],50)
up_tab = table(updiff$V1)[order(table(updiff$V1) , decreasing = T)]
hist(table(updiff$V1)[order(table(updiff$V1) , decreasing = T)])

options(repr.plot.width=5, repr.plot.height=6)

updiff = diff[which(diff$avg_log2FC>.1& diff$p_val_adj<0.05),]
updiff = updiff[which(updiff$major=="NN"),]
tab = (table(updiff$V1,updiff$celltype)[order(rowSums(table(updiff$V1,updiff$celltype)), decreasing =T),])
mtab = melt(tab[1:10,])
mtab = mtab[which(mtab$value>0),]


mtab$clade = sapply(strsplit(as.character(mtab$Var2), "_"), tail, n = 1)
mtab$clade[which(mtab$clade %in%c("Dopa", "Glut-Sero"))] = "Glut"
mtab$clade[which(mtab$clade %in%c("Gaba-Chol"))] = "Gaba"

mtab$clade[which(mtab$clade %in%c("Neur"))] = "Gaba"
ggplot(mtab, aes(x = Var1, fill=Var2)) +
  geom_bar() + theme_minimal()+coord_flip()+labs(x = "Gene", y = "Occurances DE up in aging", color = "Celltype")


options(repr.plot.width=15, repr.plot.height=6)

updiff = diff[which(diff$avg_log2FC>.1& diff$p_val_adj<0.05),]
updiff = updiff[which(updiff$major=="Glut"),]
tab = (table(updiff$V1,updiff$celltype)[order(rowSums(table(updiff$V1,updiff$celltype)), decreasing =T),])
mtab = melt(tab[1:10,])
mtab = mtab[which(mtab$value>0),]


mtab$clade = sapply(strsplit(as.character(mtab$Var2), "_"), tail, n = 1)
mtab$clade[which(mtab$clade %in%c("Dopa", "Glut-Sero"))] = "Glut"
mtab$clade[which(mtab$clade %in%c("Gaba-Chol"))] = "Gaba"

mtab$clade[which(mtab$clade %in%c("Neur"))] = "Gaba"
ggplot(mtab, aes(x = Var1, fill=Var2)) +
  geom_bar() + theme_minimal()+coord_flip()+labs(x = "Gene", y = "Occurances DE up in aging", color = "Celltype")


options(repr.plot.width=15, repr.plot.height=6)

updiff = diff[which(diff$avg_log2FC<(-.1)& diff$p_val_adj<0.05),]
updiff = updiff[which(updiff$major=="NN"),]
tab = (table(updiff$V1,updiff$celltype)[order(rowSums(table(updiff$V1,updiff$celltype)), decreasing =T),])
mtab = melt(tab[1:10,])
mtab = mtab[which(mtab$value>0),]


mtab$clade = sapply(strsplit(as.character(mtab$Var2), "_"), tail, n = 1)
mtab$clade[which(mtab$clade %in%c("Dopa", "Glut-Sero"))] = "Glut"
mtab$clade[which(mtab$clade %in%c("Gaba-Chol"))] = "Gaba"

mtab$clade[which(mtab$clade %in%c("Neur"))] = "Gaba"
ggplot(mtab, aes(x = Var1, fill=Var2)) +
  geom_bar() + theme_minimal()+coord_flip()+labs(x = "Gene", y = "Occurances DE up in aging", color = "Celltype")


updiff = diff[which(diff$avg_log2FC>.1& diff$p_val_adj<0.05),]
#updiff = updiff[which(updiff$major=="NN"),]
tab = (table(updiff$V1,updiff$celltype)[order(rowSums(table(updiff$V1,updiff$celltype)), decreasing =T),])
mtab = melt(tab[1:40,])
mtab = mtab[which(mtab$value>0),]
head(mtab)


mtab$clade = sapply(strsplit(as.character(mtab$Var2), "_"), tail, n = 1)
mtab$clade[which(mtab$clade %in%c("Dopa", "Glut-Sero"))] = "Glut"
mtab$clade[which(mtab$clade %in%c("Gaba-Chol"))] = "Gaba"

mtab$clade[which(mtab$clade %in%c("Neur"))] = "Gaba"
#pdf("Up_in_aging_genes_NN.pdf", height= 6, width = 4)
ggplot(mtab, aes(x = Var1, fill=clade)) +
  geom_bar() + theme_minimal()+coord_flip()+labs(x = "Gene", y = "Occurances DE up in aging", color = "Clade")
#dev.off()

levels(diff$celltype)

options(repr.plot.width=11, repr.plot.height=7)

colors <- (colorRampPalette(brewer.pal(9, "PuBu")[4:9])(100))

gene="Snca"
cdh8 = diff[which(diff$V1==gene & diff$p_val_adj<0.05  ),]
ggplot(cdh8, aes(x = celltype, y = avg_log2FC, fill = -log10(p_val_adj+1e-100))) +
  geom_bar(stat = "identity") +  
  scale_fill_gradientn(colors = colors) + theme_minimal()+
# You can change this to geom_bar or any other geom based on your preference
  labs(title = paste(gene,"Age LogFC (18mo/2mo)"), x = "Region", y = "Avg LogFC (18mo/2mo)", fill = "-log10(adj p-value)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))  # Rotate x-axis labels for better readability



tab = (table(updiff$V1,updiff$celltype)[order(rowSums(table(updiff$V1,updiff$celltype)), decreasing =T),])
mtab = melt(tab[1:40,])
mtab = mtab[which(mtab$value>0),]
head(mtab)


mtab$clade = sapply(strsplit(as.character(mtab$Var2), "_"), tail, n = 1)
mtab$clade[which(mtab$clade %in%c("Dopa", "Glut-Sero"))] = "Glut"
mtab$clade[which(mtab$clade %in%c("Gaba-Chol"))] = "Gaba"

mtab$clade[which(mtab$clade %in%c("Neur"))] = "Gaba"
pdf("Up_in_aging_genes.pdf", height= 6, width = 4)
ggplot(mtab, aes(x = Var1, fill=clade)) +
  geom_bar() + theme_minimal()+coord_flip()+labs(x = "Gene", y = "Occurances DE up in aging", color = "Clade")
dev.off()

options(repr.plot.width=4, repr.plot.height=6)

tab = (table(downdiff$V1,downdiff$celltype)[order(rowSums(table(downdiff$V1,downdiff$celltype)), decreasing =T),])
mtab = melt(tab[1:40,])
mtab = mtab[which(mtab$value>0),]
head(mtab)


mtab$clade = sapply(strsplit(as.character(mtab$Var2), "_"), tail, n = 1)
mtab$clade[which(mtab$clade %in%c("Dopa", "Glut-Sero"))] = "Glut"
mtab$clade[which(mtab$clade %in%c("Gaba-Chol"))] = "Gaba"

mtab$clade[which(mtab$clade %in%c("Neur"))] = "Gaba"
pdf("Down_in_aging_genes.pdf", height= 6, width = 4)
ggplot(mtab, aes(x = Var1, fill=clade)) +
  geom_bar() + theme_minimal()+coord_flip()+labs(x = "Gene", y = "Occurances DE down in aging", color = "Clade")
dev.off()

table(updiff$V1)[order(table(updiff$V1) , decreasing = T)][grep("Pcdh", names(table(updiff$V1)[order(table(updiff$V1) , decreasing = T)]))]


nrow(diff)
head(table(diff$V1)[order(table(diff$V1) , decreasing = T)])
hist(table(diff$V1)[order(table(diff$V1) , decreasing = T)])

mer_genes = read.table("../../../merfish/all_mergenes.txt")

mer_genes = mer_genes$V1

length(which(names(down_tab)%in%mer_genes))

length(which(names(up_tab)%in%mer_genes))

up_tab[which(names(up_tab)%in%mer_genes)]

down_tab[which(names(down_tab)%in%mer_genes)]

diff$V1[which(diff$celltype==i)]

dir ="down"
if(dir == "up" ) {
    diff = big_diff[which(big_diff$avg_log2FC>0 & big_diff$p_val_adj<0.05),]
} else {
    diff = big_diff[which(big_diff$avg_log2FC<0 & big_diff$p_val_adj<0.05),]
}
mat = matrix(0,nrow=length(unique(diff$celltype)),ncol=length(unique(diff$celltype)))
x = 1
for (i in unique(diff$celltype)){
    y=1
    set1 = diff$V1[which(diff$celltype==i)]
    # coors1 = set1$V1
    for (j in unique(diff$celltype)){
     # print(paste(i,j))
      set2 = diff$V1[which(diff$celltype==j)]
      #  coors2 = set2$V1
      mat[x,y] = length(which(set1 %in% set2))/length(union(set1, set2))
      y = y+1
    }
    x = x+1
}
  
melted = melt(mat)
  
  
mat[which(mat==1, arr.ind = T)] = 0.25
mat[which(mat>0.25, arr.ind = T)] = 0.25
rownames(mat) = unique(diff$celltype)
colnames(mat) = unique(diff$celltype)  



#annotation_col = data.frame(
#    Clade = clades,
#    Sex = sex
#  )
#rownames(annotation_col) = colnames(mat) 

clade = sapply(strsplit(as.character(colnames(mat)), "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"



annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

library(pheatmap)

options(repr.plot.width=12.5, repr.plot.height=10)

pheatmap(mat, main = paste( "Jaccard of differential genes up"),show_colnames = F , annotation_row = annotation_col)


options(repr.plot.width=12.5, repr.plot.height=10)

pheatmap(mat, main = paste( "Jaccard of differential genes down"),show_colnames = F , annotation_row = annotation_col)



