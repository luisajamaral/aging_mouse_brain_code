library(data.table)
library(dplyr)
library(ggplot2)
setwd("../../h5ads_final/all_peaks/")

setwd("../../combined_DARs_redone/")



files = list.files(".","ann.txt")
files = files[grep("up", files)]
files[1]

f = 'down_peaks_Oligo_NN.ann.txt'
ann = fread(f)
    ct = gsub(".ann.txt" , "", f)
    ct = gsub("up_peaks_" , "", ct)

    ann$celltype = ct
    ann$ann = sapply(strsplit(as.character(ann$Annotation), "[(]"), `[`, 1)
    ann$ann = gsub(" ", "", ann$ann)
    res = cbind(names(table(ann$ann)), table(ann$ann),ct)
head(ann)

files = list.files(".","ann.txt")
files = files[grep("up", files)]
uall = list()
for(f in files){
    ann = fread(f)
    ct = gsub(".ann.txt" , "", f)
    ct = gsub("up_peaks_" , "", ct)

    ann$celltype = ct
    ann$ann = sapply(strsplit(as.character(ann$Annotation), "[(]"), `[`, 1)
    ann$ann = gsub(" ", "", ann$ann)
    res = cbind(names(table(ann$ann)), table(ann$ann),ct)
    uall[[ct]]=res
}

files = list.files(".","ann.txt")
files = files[grep("down", files)]
dall = list()
for(f in files){
    ann = fread(f)
    ct = gsub(".ann.txt" , "", f)
    ct = gsub("down_peaks_" , "", ct)

    ann$celltype = ct
    ann$ann = sapply(strsplit(as.character(ann$Annotation), "[(]"), `[`, 1)
    ann$ann = gsub(" ", "", ann$ann)
    res = cbind(names(table(ann$ann)), table(ann$ann),ct)
    dall[[ct]]=res
}

udata = do.call(rbind, uall)
udata = as.data.frame(udata)
colnames(udata) = c("Annotation", "count", "celltype")
head(udata)

ddata = do.call(rbind, dall)
ddata = as.data.frame(ddata)
colnames(ddata) = c("Annotation", "count", "celltype")
head(ddata)

ddata$direction = "Down"
udata$direction = "Up"


adata= rbind(ddata, udata)

head(adata)

chromHMM_colors <- c(
  "3'UTR" = "#ff7f00",         # Orange
  "miRNA" = "#377eb8",        # Blue
  "ncRNA" = "#4daf4a",        # Green
  "TTS" = "#984ea3",          # Purple
  "pseudo" = "#e41a1c",       # Red
  "exon" = "#f781bf",         # Pink
  "intron" = "#a65628",       # Brown
  "Intergenic" = "#999999",   # Gray
  "promoter-TSS" = "#ffff33",     # Yellow
  "5'UTR" = "#ff7f00",         # Orange (same as 3UTR, or you can choose another distinct color)
  "snoRNA" = "#984ea3",       # Purple (same as TTS, or you can choose another distinct color)
  "rRNA" = "#4daf4a"          # Green (same as ncRNA, or you can choose another distinct color)
)

adata$count = as.numeric(adata$count)
head(adata$count[which(adata$direction=="Down")])
adata$count[which(adata$direction=="Down")] = -(adata$count[which(adata$direction=="Down")])


head(adata[which(adata$celltype=="Oligo_NN"),])

meta = read.csv("../final_meta.csv")
meta$celltype_final = gsub("/","-", meta$celltype_final)
meta$celltype_final = gsub("Lymphoid","Lymphoid NN", meta$celltype_final)
meta = meta[-which(meta$celltype_final=="doublet"),]
meta$celltype_final=paste(meta$celltype_final)
meta$celltype_final[which(meta$best_celltype_fixed=="IOL")] = "IOL NN"
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)


head(adata)

options(repr.plot.width = 16, repr.plot.height = 7)
adata$celltype = gsub("Lymphoid","Lymphoid_NN", adata$celltype)

adata$major = sapply(strsplit(as.character(adata$celltype), "_"), function(x) tail(x, 1))
adata$major[which(adata$major%in%"Neur")] = "Gaba"
adata$major[which(adata$major%in%"Gaba-Chol")] = "Gaba"
adata$major[which(adata$major%in%"Glut-Sero")] = "Glut"
adata$major[which(adata$major%in%"Dopa")] = "Glut"
adata$major[which(adata$major%in%"IMN")] = "NN"

#adata$celltype = factor(adata$celltype, levels = ord)
ggplot(adata, aes(x = celltype, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "stack") +
  #scale_fill_manual(values = chromHMM_colors) +  # Apply the chromHMM color palette
  theme_minimal() +
  labs(
    title = "Age-DARs by Cell Type",
    x = "Cell Type",
    y = "num peaks",
    fill = "Annotation"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) + 
  coord_flip() + 
  facet_wrap(~ major, scales = "free")

options(repr.plot.width=14, repr.plot.height=7)

ggplot(adata, 
  aes(x = celltype, y = `log2(fold_change)`, 
  color = ifelse(`adjusted p-value` < 0.05, ifelse(`log2(fold_change)` > 0, "up", "down"), "not sig"))) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("blue", "grey", "red")) +
  labs(title = "Age Differential Accessible Regions",
       x = "Cell Type",
       y = "log2 Fold Change",
       color = "adj p < 0.05") +ylim(c(-8,8))+
  theme_minimal() + coord_flip() + facet_wrap(~region, scales = "free",nrow=2)

setwd("../h5ads_final/all_peaks/")
files = list.files(".","annStats")

files


meta$celltype_final = gsub(" ", "_", meta$celltype_final)
ord = names(table(meta$celltype_final)[order(table(meta$celltype_final))])
meta$celltype_final <- factor(meta$celltype_final, levels = ord)




f = "DTN-LDT-IPN_Otp_Pax3_Gaba_allpeaks_annStats"
ann = fread(f)
   # ann = ann[1:(grep("Ann", ann$Annotation)-1),]
    ct = gsub("_allpeaks_annStats" , "", f)
    ann$celltype = ct

head(ann)

all = list()
for(f in files){
    ann = fread(f)
    ann = ann[1:(grep("Ann", ann$Annotation)-1),]
    ct = gsub("_allpeaks_annStats" , "", f)
    ann$celltype = ct
    all[[ct]]=ann
}

data = do.call(rbind, all)
data <- data %>%
  mutate(`Number of peaks` = as.numeric(`Number of peaks`))

data$celltype = gsub("Lymphoid","Lymphoid_NN", data$celltype)

data$major = sapply(strsplit(as.character(data$celltype), "_"), function(x) tail(x, 1))
data$major[which(data$major%in%"Neur")] = "Gaba"
data$major[which(data$major%in%"Gaba-Chol")] = "Gaba"
data$major[which(data$major%in%"Glut-Sero")] = "Glut"
data$major[which(data$major%in%"Dopa")] = "Glut"
data$major[which(data$major%in%"IMN")] = "NN"



data[which(is.na(data$celltype)),]

data = data[which(!is.na(data$`Number of peaks`)),]

chromHMM_colors <- c(
  "3UTR" = "#ff7f00",         # Orange
  "miRNA" = "#377eb8",        # Blue
  "ncRNA" = "#4daf4a",        # Green
  "TTS" = "#984ea3",          # Purple
  "pseudo" = "#e41a1c",       # Red
  "Exon" = "#f781bf",         # Pink
  "Intron" = "#a65628",       # Brown
  "Intergenic" = "#999999",   # Gray
  "Promoter" = "#ffff33",     # Yellow
  "5UTR" = "#ff7f00",         # Orange (same as 3UTR, or you can choose another distinct color)
  "snoRNA" = "#984ea3",       # Purple (same as TTS, or you can choose another distinct color)
  "rRNA" = "#4daf4a"          # Green (same as ncRNA, or you can choose another distinct color)
)

data$celltype = factor(data$celltype, levels = ord)

data$type = "Non-TSS"
data$type[which(data$Annotation=="Promoter")] = "TSS"



options(repr.plot.width = 20, repr.plot.height = 7)

ggplot(data, aes(x = celltype, y = `Number of peaks`, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  #scale_fill_manual(values = chromHMM_colors) +  # Apply the chromHMM color palette
  theme_minimal() +
  labs(
    title = "Annotation Distribution by Cell Type",
    x = "Cell Type",
    y = "MACS3 called Peaks",
    fill = "Annotation"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) + 
  coord_flip() + 
  facet_wrap(~ major, scales = "free")

g = ggplot(data, aes(x = celltype, y = `Number of peaks`, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = chromHMM_colors) +  # Apply the chromHMM color palette
  theme_minimal() +
  labs(
    title = "Called Peaks by Cell Type",
    x = "Cell Type",
    y = "MACS3 called Peaks",
    fill = "Annotation"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) + 
  coord_flip() + 
  facet_wrap(~ major, scales = "free")

ggsave(g, file = "~/projects/combined_all/Figures/Figure1-overview/peak_annotation_bar.pdf", height = 10, width = 20)

options(repr.plot.width = 5, repr.plot.height = 5)

ggplot(data, aes(x = "", y = `Number of peaks`, fill = Annotation)) +
  geom_bar(stat = "identity", width = 1, position = "fill") +
  scale_fill_manual(values = chromHMM_colors) +  # Apply the chromHMM color palette
  theme_minimal() +
  labs(
    title = "Annotation Distribution by Cell Type",
    x = "",
    y = "",
    fill = "Annotation"
  ) +
  theme(
    #axis.title = element_text(size = 14),
    #axis.text = element_text(size = 12),
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text = element_blank()    # Remove axis text
  ) +
  coord_polar(theta = "y") +  # Convert to pie chart using polar coordinates
  geom_text(aes(label = paste0(round((..count../sum(..count..))*100, 1), "%")),
            position = position_stack(vjust = 0.5), stat = "count")  # Add labels
# +  # Convert to pie chart using polar coordinates
  #facet_wrap(~ major)  # One pie chart per cell type

library(dplyr)
library(ggplot2)

# Calculate the percentage for each slice and position for the labels
data <- data %>%
  group_by(celltype, Annotation) %>%
  summarise(`Number of peaks` = sum(`Number of peaks`), .groups = 'drop') %>%
  mutate(percentage = `Number of peaks` / sum(`Number of peaks`) * 100) %>%
  arrange(desc(Annotation)) %>%
  mutate(label_pos = cumsum(percentage) - 0.5 * percentage)

# Plot the pie chart
options(repr.plot.width = 5, repr.plot.height = 5)

ggplot(data, aes(x = "", y = `Number of peaks`, fill = Annotation)) +
  geom_bar(stat = "identity", width = 1, position = "fill") +
  scale_fill_manual(values = chromHMM_colors) +  # Apply the chromHMM color palette
  theme_minimal() +
  labs(
    title = "Annotation Distribution for all peaks",
    x = "",
    y = "",
    fill = "Annotation"
  ) +
  theme(
    panel.grid = element_blank(),  # Remove gridlines
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.text = element_blank()    # Remove axis text
  ) +
  coord_polar(theta = "y") +  # Convert to pie chart using polar coordinates
  geom_text(aes(label = paste0(round(percentage, 1), "%"), y = label_pos), 
            position = position_stack(vjust = 0.5), size = 3)  # Add labels


options(repr.plot.width = 20, repr.plot.height = 7)

ggplot(data, aes(x = celltype, y = `Number of peaks`, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = chromHMM_colors) +  # Apply the chromHMM color palette
  theme_minimal() +
  labs(
    title = "Annotation Distribution by Cell Type",
    x = "Cell Type",
    y = "Proportion of Peaks",
    fill = "Annotation"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) + 
  coord_flip() + 
  facet_wrap(~ major, scales = "free")

options(repr.plot.width = 6, repr.plot.height = 12)

ggplot(data, aes(x = major, y = `Number of peaks`, fill = Annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = chromHMM_colors) +  # Apply the chromHMM color palette
  theme_minimal() +
  labs(
    title = "Called peaks",
    x = "Cell Type",
    y = "Proportion of Peaks",
    fill = "Annotation"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) + 
  coord_flip() 

options(repr.plot.width=20, repr.plot.height=7)
ggplot(data, aes(x = celltype, y = `Number of peaks`, fill = Annotation)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(
    title = "Annotation Distribution by Cell Type",
    x = "Cell Type",
    y = "Proportion of Peaks",
    fill = "Annotation"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  ) + coord_flip()+ 
  facet_wrap(~ major, scales = "free") 

unique(data$Annotation)


sum(ann[,2])


