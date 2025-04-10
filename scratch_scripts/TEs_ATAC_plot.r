library(data.table)
library(edgeR)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

lkey=read.table("out.saf")
#lkey = lkey[-which(lkey$V5=="Retroposon"),]
ukey = lkey[-which(duplicated(lkey$V1)),]
rownames(ukey) = ukey$V1
LINES = unique(paste(lkey[which(lkey$V6%in%c("LINE", "LINE?")),1]))
LTR = unique(paste(lkey[which(lkey$V6%in%c("LTR", "LTR?")),1]))
SINE = unique(paste(lkey[which(lkey$V6%in%c("SINE","SINE?")),1]))
Retroposon = unique(paste(lkey[which(lkey$V6%in%c("Retroposon")),1]))


meta = fread("../final_meta.csv")
head(meta)

sum_by_group <- meta %>%
  group_by(sample, celltype_final) %>%
  summarise(sum_of_values = sum(n_fragment))
sum_by_group = as.data.frame(sum_by_group)
sum_by_group$celltype_final=gsub(" ", "_", sum_by_group$celltype_final)
sum_by_group$celltype_final=gsub("/", "_", sum_by_group$celltype_final)
sum_by_group$sample=gsub("8wk", "2mo", sum_by_group$sample)

rownames(sum_by_group) = paste(sum_by_group$celltype_final, sum_by_group$sample, sep = '.' )
head(sum_by_group)

head(sum_by_group)

tab = fread("allTE.feature.counts.txt")
tab = as.data.frame(tab)
tab = tab[,-which(duplicated(colnames(tab)))]



dim(tab)

tab = tab[,-c(2:5)]

colnames(tab) = gsub("../female_split_bams/", "", colnames(tab))
colnames(tab) = gsub("../male_split_bams/", "", colnames(tab))
colnames(tab) = gsub(".bam", "", colnames(tab))

write.table(tab, "TE_counts_by_sample.txt", sep = "\t", quote = F)

rownames(tab) = tab$Geneid
btab = as.data.frame(tab[,-c(1:2)])
btab = btab[-grep("Astroepen", colnames(btab))]
rownames(btab) = tab$Geneid
sample = colnames(btab)

df =btab



norm_df = df %>%
    mutate(across(everything(), ~ . / sum(.)))

L1ma5a = norm_df["MuRRS-int",]

head(paste(L1ma5a[1,]))

age= rep("2mo",length(L1ma5a))
age[grep("9mo", names(L1ma5a))] = "9mo"
age[grep("18mo", names(L1ma5a))] = "18mo"

clade= rep("NN",length(L1ma5a))
clade[grep("Glut", names(L1ma5a))] = "Glut"
clade[grep("Gaba", names(L1ma5a))] = "Gaba"
clade[grep("Neur", names(L1ma5a))] = "Gaba"

region= rep("AMY",length(L1ma5a))
region[grep("HCA", names(L1ma5a))] = "HCA"
region[grep("NAC", names(L1ma5a))] = "NAC"
region[grep("RLP", names(L1ma5a))] = "RLP"
region[grep("ENT", names(L1ma5a))] = "ENT"
region[grep("FC", names(L1ma5a))] = "FC"
region[grep("CP", names(L1ma5a))] = "CP"
region[grep("HCP", names(L1ma5a))] = "HCP"


sex= rep("Male",length(L1ma5a))
sex[grep("Female", names(L1ma5a))] = "Female"


cur = cbind(names(L1ma5a),norm_acc= paste(L1ma5a[1,]), age,clade,region, sex)
head(cur)

cur = as.data.frame(cur)


cur$norm_acc = as.numeric(as.character(cur$norm_acc))

cur$age= factor(cur$age, levels = c("2mo", "9mo", "18mo"))

library(ggplot2)
options(repr.plot.width=10
        , repr.plot.height=10)
p <- ggplot(cur, aes(x=clade, y=(norm_acc), fill = age)) + 
  geom_violin(draw_quantiles = c(.25,.5,.75)) + facet_wrap(~region) #+ylim(c(0,0.001))
p

p <- ggplot(cur, aes(x=clade, y=norm_acc, fill=age, color=sex)) + 
  geom_violin(draw_quantiles=c(0.25, 0.5, 0.75)) + 
  facet_wrap(~region) +
  stat_summary(fun=median, geom="point", shape=20, size=2, color="black", fill="black") + 
  theme_minimal() +
  labs(title="Violin Plot of Normalized Accessibility by Clade and Region",
       x="Clade",
       y="Normalized Accessibility")

print(p)

options(repr.plot.width=5
        , repr.plot.height=5)
ggplot(cur[which(cur$region=="FC"),], aes(x=clade, y=(norm_acc), fill = age)) + 
  geom_boxplot(draw_quantiles = c(.5))  

colnames(df)

 col_names <- colnames(df)
  cell_types <- unique(sub("\\..*", "", col_names))
  ages <- c("2mo", "9mo", "18mo")
  
  # Initialize a new data frame to store the combined results with the same number of rows as df
  combined_data <- df[0, , drop = FALSE] # Initialize with zero rows but keep column structure
  combined_data <- data.frame(matrix(ncol = 0, nrow = nrow(df))) # Create an empty data frame with the same number of rows
  
  for (cell_type in cell_types) {
    for (age in ages) {
      pattern <- paste0(cell_type, ".*", age)
      matching_cols <- col_names[grepl(pattern, col_names)]
      if (length(matching_cols) > 0) {
        new_col_name <- paste(cell_type, age, sep = "_")
        combined_data[[new_col_name]] <- rowSums(df[, matching_cols, drop = FALSE], na.rm = TRUE)
      }
    }
  }
  raw = combined_data
  # Normalize each value by dividing it by the sum of its column
  combined_data <- combined_data %>%
    mutate(across(everything(), ~ . / sum(.)))
  



rownames(combined_data)=rownames(df)

sums = (sum_by_group[paste(sample),])

sample[which(!sample%in%rownames(sum_by_group))]
#rownames(sum_by_group)[which(!rownames(sum_by_group) %in%sample)]

head(sums$sum_of_values)

head(combined_data)

#nbtab = combined_data
#for(i in 1:ncol(btab)){
#    nbtab[,i] = btab[,i]/sums$sum_of_values[i]*1000000
    
#}
#head(nbtab)
#btab = nbtab

btab = combined_data
hist(rowSums(btab))

btab=btab[-which(rowSums(btab)<.1),]

clade = rep("IMN", ncol(btab))
clade[grep("_NN.", colnames(btab))] = "NN"
clade[grep("_Glut.", colnames(btab))] = "Glut"
clade[grep("_Gaba.", colnames(btab))] = "Gaba"

age = rep("2mo", ncol(btab))
age[grep("9mo", colnames(btab))] = "9mo"
age[grep("18mo", colnames(btab))] = "18mo"




head(cpm)

library(MatrixGenerics)

hist(colSums(raw))
gc = names(colSums(raw)[which(colSums(raw)>40000000)])
gc = gsub("_2mo", "", gc)
gc = gsub("_9mo", "", gc)
gc = gsub("_18mo", "", gc)
gc = unique(gc)
gc = gc[-which(gc=="doublet")]

cpm = btab

cts = colnames(cpm)
cts = gsub("_2mo", "", cts)
cts = gsub("_9mo", "", cts)
cts = gsub("_18mo", "", cts)

cpm = cpm[,which(cts %in%gc)]
cpm=cpm[which(as.vector(rowMaxs(as.matrix(cpm)))>.0025),]

type = ukey[paste(rownames(cpm)),6]

cpm = cpm[which(type%in%c("SINE",  "LINE", "LTR")),]
type = ukey[paste(rownames(cpm)),6]

anno_col<-data.frame( Type=type)
rownames(anno_col) = rownames(cpm)



clade = rep("Gaba", ncol(cpm))
clade[grep("_NN.", colnames(cpm))] = "NN"
clade[grep("_Glut.", colnames(cpm))] = "Glut"
clade[grep("_Gaba.", colnames(cpm))] = "Gaba"

age = rep("2mo", ncol(cpm))
age[grep("9mo", colnames(cpm))] = "9mo"
age[grep("18mo", colnames(cpm))] = "18mo"

anno<-data.frame(Age = age, Clade=clade)
rownames(anno) = colnames(cpm)

dim(cpm)
cpm = as.data.frame(cpm)

#cpm=cpm[,-which(colSums(cpm)<25)]

options(repr.plot.width=15, repr.plot.height=10)

pheatmap((t(cpm)), scale = "column", show_colnames = T, cluster_rows = F, gaps_row = seq(3,48, 3),
         annotation_row = anno, annotation_col = anno_col )




head(cpm,20)

anno

rownames(cpm) == rownames(anno_col)

#cpm = cpm(btab)
#cpm  = cpm[-which(rowSums(cpm) ==0),]
#cpm[which(cpm>5000, arr.ind = T)] = 5000
type = ukey[paste(rownames(cpm)),6]

cpm = cpm[which(type%in%c("SINE",  "LINE", "LTR")),]
type = ukey[paste(rownames(cpm)),6]

anno<-data.frame(row.names=rownames(cpm), Group=type)
colnames(anno) = "Type"
# define the colours
acc = brewer.pal(3,"Accent")
set1 = brewer.pal(3,"Set1")


clade = rep("IMN", ncol(btab))
clade[grep("_NN.", colnames(btab))] = "NN"
clade[grep("_Glut.", colnames(btab))] = "Glut"
clade[grep("_Gaba.", colnames(btab))] = "Gaba"

age = rep("2mo", ncol(btab))
age[grep("9mo", colnames(btab))] = "9mo"
age[grep("18mo", colnames(btab))] = "18mo"


sex = rep("Male", ncol(btab))
sex[grep("Female", colnames(btab))] = "Female"


#head(cpm)
cpm = btab
cpm = cpm[,c(grep("DG_", colnames(cpm)), grep("Oligo", colnames(cpm)))]
cpm = cpm[,grep("Female", colnames(cpm))]
cpm = cpm[,grep(":HCA", colnames(cpm))]

#cpm = cpm[,grep("ENT", colnames(cpm))]
cpm=cpm[-which(rowSums(cpm)<1000),]
dim(cpm)
#cpm=cpm[,-which(colSums(cpm)<25)]

options(repr.plot.width=12, repr.plot.height=12)

pheatmap(t(cpm), scale = "column", show_colnames = F,annotation_col = anno )



which(rowSums(cpm)<500)




colnames(clade) = "Clade"
annoCol<-list(Type=c(SINE="#009E73", LINE="#CC79A7", LTR="#0072B2"), Clade=c(ExN = acc[1], InN = acc[2], Glia = acc[3]))

pdf("FC_scTE_heat_te_allcl.pdf", height = 10, width = 7)
pheatmap(t(cpm), scale = "column", show_colnames = F, annotation_row  = clade ,
         annotation_col = anno,annotation_colors = annoCol)

pheatmap(t(cpm), scale = "column", show_colnames = F, annotation_row  = clade ,
         annotation_col = anno,annotation_colors = annoCol,color=colorRampPalette(c("navy", "white", "red"))(50))
dev.off()

mcpm = list()
for (i in unique(cl)) {
  mini = cpm[,which(cl == i)]
  rm = rowMeans(mini)
  mcpm[[i]] = rm
}
mecpm = do.call(cbind, mcpm)
clade = key[colnames(mecpm),4]
clade <-data.frame(row.names=colnames(mecpm), Group=clade)
colnames(clade) = "Clade"

pdf("FC_scTE_heat_te_allcl_tog.pdf", height = 7, width = 8)
pheatmap(t(mecpm), scale = "column", show_colnames = F, annotation_row  = clade , labels_row = key[colnames(mecpm),3],
         annotation_col = anno,annotation_colors = annoCol,color=colorRampPalette(c("navy", "white", "red"))(50))
pheatmap(t(mecpm), scale = "column", show_colnames = F, annotation_row  = clade , labels_row = key[colnames(mecpm),3],
         annotation_col = anno,annotation_colors = annoCol)

dev.off()



