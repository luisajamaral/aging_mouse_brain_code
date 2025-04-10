library(ggplot2)
library(data.table)
library(UpSetR)
library(pheatmap)
library(Seurat)
library(zellkonverter)

setwd("../../h5ads_final/")

h = readH5AD("celltype_major_batch_age_region_PMAT_RPM.h5ad")


h = readH5AD("../h5ads_final/celltype_batch_age_region_PMAT_RPM.h5ad")


count_table <-h@assays@data$X


setwd("../combined_DARs_redone/")
files = list.files(".", "diff_peaks", full.names = T)
#rfiles = list.files("../region_DARs_redone/", "diff_peaks", recursive = T, full.names=T)
#files
#rfiles

f = files[grep("Oligo", files)]
dar_data = read.table(f, header=TRUE, sep=',', stringsAsFactors=FALSE)

nrow(dar_data)

significant_changes <- subset(dar_data, adjusted.p.value < 0.01 & abs(log2.fold_change.)>0.5)
significant_changes = significant_changes[-grep("chrY", significant_changes$feature.name),]


nrow(significant_changes)

oligo = count_table[paste(significant_changes[,"feature.name"]),grep("DG Glut", colnames(count_table))]


dim(oligo)

oligo = oligo[,c(grep("HC", colnames(oligo)))]
head(oligo)

o = oligo
region = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 4)
age = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 3)
sex = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 2)
sex = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 2)
celltype = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 1)

#age = sapply(strsplit(as.character(age), "_"), "[[", 1)
annotation_col = data.frame(
    age = age,
    region = region,
    sex = sex,
    celltype = celltype
    
  )
annotation_col$age =factor(annotation_col$age, levels = c("2mo","9mo", "18mo"))
rownames(annotation_col) = colnames(o)

ph = pheatmap(oligo[1:4000,],
              main = paste(nrow(oligo), "DARs"),
              show_rownames = F,scale = "row",
             # annotation_col=annotation_col,
              color=colorRampPalette(c("navy", "white", "red"))(40),  
              cluster_rows = T,cluster_cols = F)#,


which(rowSums(oligo[,grep("Female", colnames(oligo))])==0)

pheatmap(oligo[,grep("Male", colnames(oligo))],
              main = paste(nrow(oligo), "DARs"),
              show_rownames = F,scale = "row",
              #annotation_col=annotation_col,
              color=colorRampPalette(c("navy", "white", "red"))(40),  
              cluster_rows = T,cluster_cols = F)#,


color_age = read.csv("../Figures/color_scheme/Age.csv")
color_region = read.csv("../Figures/color_scheme/MajorRegion-id.csv")
color_sex = read.csv("../Figures/color_scheme/Sex.csv")

# Convert to named color vectors
age_colors <- setNames(color_age$Color, color_age$Age)
region_colors <- setNames(color_region$Color, color_region$Region)
sex_colors <- setNames(color_sex$Color, color_sex$Sex)

# Combine color annotations into a list
annotation_colors <- list(
  age = age_colors,
  region = region_colors,
  sex = sex_colors
)

annotation_col = annotation_col[,c(1:3)] 

BuRd <- color("BuRd")
plot_scheme(BuRd(9), colours = TRUE, size = 0.9)


PRGn <- color("PRGn")
plot_scheme(PRGn(9), colours = TRUE, size = 0.9)

significant_changes$direction = "up"
significant_changes$direction[which(significant_changes$log2.fold_change.<0)] ="down" 
head(significant_changes)

annotation_col$

options(repr.plot.width=9, repr.plot.height=7)

# Split the data into male and female samples
female_samples <- annotation_col$sex == "Female"
male_samples <- annotation_col$sex == "Male"

# Scale each subset separately (row-wise scaling)
scaled_female <- t(scale(t(oligo[, female_samples])))  # Scale females
scaled_male <- t(scale(t(oligo[, male_samples])))      # Scale males

# Combine scaled data into a single matrix
scaled_oligo <- cbind(scaled_female, scaled_male)

# Ensure columns are in the original order
scaled_oligo <- scaled_oligo[, order(c(which(female_samples), which(male_samples)))]
colnames(scaled_oligo) = gsub("_", " ", colnames(scaled_oligo)) 
rownames(annotation_col) = gsub("_", " ", rownames(annotation_col)) 


# Plot the heatmap
p = pheatmap(cutree_rows = 2,clustering_method = "ward.D2",
  scaled_oligo[,],
  main = 
  paste("Oligodendrocyte Age-DARs\n", 
        table(significant_changes$direction)[[2]], " up | ",
        table(significant_changes$direction)[[1]]," down", sep = "") ,

  show_rownames = FALSE,
  scale = "none",  # Data is already scaled
  annotation_col = annotation_col,
  breaks = seq(-3,3,length.out=100), color = BuRd(100),
  cluster_rows = TRUE,gaps_col = c(24),
  cluster_cols = FALSE,  annotation_colors = annotation_colors  # Add the color schemes

)

pdf("../Figures/Figure3-oligos/Oligo_DAR_heatmap.pdf", width = 9)
p
dev.off()

while (1){dev.off()}

head(significant_changes)

# Split the data into up and down peaks
up_peaks <- significant_changes[significant_changes$direction == "up", ]
down_peaks <- significant_changes[significant_changes$direction == "down", ]

# Function to parse feature.name into BED columns
parse_feature_name <- function(data) {
  bed <- data.frame(
    chr = sub(":(.*)", "", data$feature.name),                  # Extract chromosome
    start = as.numeric(sub(".*:(\\d+)-.*", "\\1", data$feature.name)),  # Extract start
    end = as.numeric(sub(".*-(\\d+)", "\\1", data$feature.name))        # Extract end
  )
  return(bed)
}

# Create BED data for up and down peaks
up_bed <- parse_feature_name(up_peaks)
down_bed <- parse_feature_name(down_peaks)

# Write the BED files
write.table(up_bed, file = "../Figures/Figure3-oligos/up_peaks.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(down_bed, file = "../Figures/Figure3-oligos/down_peaks.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


