library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringdist)
library(stats)
library(MatrixGenerics)
library(pheatmap)
library(khroma)
BuRd <- color("BuRd")


database = "motifs"


setwd(paste("../../../combined_all/combined_DAR_yng_vs_old_progenitors/",sep = ""))

pval = 1e-5
files = list.files(".", pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)
files = paste(files, sep = "")
files = files[which(file.exists(files))]
length(files)
ind_files = files[grep("IOL",files)]
files

files = files[-grep("0.1", files)]
#files= files[grep("OB-STR-CTX_Inh_IMN",files)]

database = "motifs"
pval = 0.01
gkey = c()
f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
curr = curr[which(as.numeric(curr$`P-value`) < pval)[1:5],]

gkey = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  #  show(head(curr))
  curr = curr[which(as.numeric(curr$`P-value`) < pval),]
  curr = curr[1:min(5, nrow(curr)),]
  #cat(f, nrow(curr),"\n")

  gkey = c(gkey , paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1)))

}
gkey = gkey[which(!duplicated(gkey))]
cat(length(gkey))

gkey = gkey[which(gkey != 'NA')]
gkey
# getting significance for each term + file
mat = matrix(0, nrow=length(gkey),ncol=length(files))
smat = matrix('', nrow=length(gkey),ncol=length(files))

x = 1
for(i in 1:length(gkey)) {
  y = 1
  for(f in files) {
    curr = fread(f, sep = "\t", fill=T)
    curr = as.data.frame(curr)
    names = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
    which(duplicated(names))
    curr = curr[which(!duplicated(names)),]
    rownames(curr) = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
    if(paste(gkey[i]) %in% paste(rownames(curr))) {
      # mat[x,y]  = -log10(curr[paste(gkey[i]), "P.Value"]+1e-100)
      mat[x,y]  = -curr[paste(gkey[i]),"Log P-value"]
      if(curr[paste(gkey[i]),"q-value (Benjamini)"]<0.05){
        smat[x,y] = '*'
      }
    } else {
      mat[x,y]  = NA
    }
    y= y+1
  }
  x=x+1
}


rownames(mat) = gkey
colnames(mat) = gsub("_motifs/knownResults.txt", "",files)
colnames(mat) = gsub("./motifs_genOn_GO_results/", "", colnames(mat))
colnames(mat) = gsub("peaks_", "", colnames(mat))
colnames(mat) = gsub("../combined_DARs_redone//motifs/", "combined_", colnames(mat))





snmat = smat
nmat = mat 
#snmat = smat[which(rowMaxs(mat)>5),which(colMaxs(mat)>5)]
#nmat = mat[which(rowMaxs(mat)>5),which(colMaxs(mat)>5)]
#nmat[which(nmat>100, arr.ind = T)] = 100

dir = rep("down" , ncol(nmat))
dir[grep("up", as.character(colnames(nmat)))] = "up"

region = sapply(strsplit(colnames(nmat), ":"), `[`, 2)

ct = sapply(strsplit(colnames(nmat), ":"), `[`, 1)
ct = sapply(strsplit(ct, "_"), `[`, 2)


nannotation_col = data.frame(
    Direction = dir, Celltype= ct, Region = region
  )
rownames(nannotation_col) = colnames(nmat) 

snmat = snmat[,order(dir, ct,region)]
nmat = nmat[,order(dir,ct,region)]

#nmat[which(nmat>100, arr.ind = T)] = 100

nannotation_col = nannotation_col[,-c(2)]
color_age = read.csv("../Figures/color_scheme/Age.csv")
color_region = read.csv("../Figures/color_scheme/MajorRegion-id.csv")
color_sex = read.csv("../Figures/color_scheme/Sex.csv")
color_region = rbind(color_region , c("Combined", "white"))
# Convert to named color vectors
age_colors <- setNames(color_age$Color, color_age$Age)
region_colors <- setNames(color_region$Color, color_region$Region)


color_sex$Sex = c("up", "down")
color_sex$Color = c("red", "blue")
colnames(color_sex)[1] = "Direction"
sex_colors <- setNames(color_sex$Color, color_sex$Sex)
sex_colors <- setNames(color_sex$Color, color_sex$Direction)

# Combine color annotations into a list
annotation_colors <- list(
  Direction = sex_colors,
  Region = region_colors
)

nannotation_col = nannotation_col[,c(2,1)]

rownames(nannotation_col) = gsub("_", " ", rownames(nannotation_col))
rownames(nannotation_col) = gsub(":", " ", rownames(nannotation_col))


colnames(nmat) = gsub("_", " ", colnames(nmat))
colnames(nmat) = gsub(":", " ", colnames(nmat))

colnames(nmat) = gsub("./", "", colnames(nmat))

snmat = snmat[,order(colnames(nmat))]
nmat = nmat[,order(colnames(nmat))]

snmat = snmat[,c(grep("down", colnames(nmat)), grep("up", colnames(nmat)))]
nmat = nmat[,c(grep("down", colnames(nmat)), grep("up", colnames(nmat)))]


options(repr.plot.width=4.5, repr.plot.height=5.5)
nmat[which(nmat>30, arr.ind = T)] = 30
colnames(nmat) = gsub("motif", "" ,colnames(nmat))

heatmap <- pheatmap(nmat[nrow(nmat):1,],
                    gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    #annotation_col = nannotation_col[,2, drop = F], 
                    cluster_cols = F, color = BuRd(20)[2:20],
                    main = "OB-STR-CTX Inh IMN\nAge-DAR motif enrichment",              
                    clustering_method = 'ward.D2',#gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    #annotation_colors = annotation_colors[c(2,1)],
                    display_numbers = snmat[nrow(nmat):1,], 
                    cluster_rows = T)

pdf("../Figures/Figure2-IMN-IOL/all_motifs.pdf", height = 5.5, width = 4.5)
heatmap
dev.off()

options(repr.plot.width=4, repr.plot.height=3)

# Load necessary library
library(ggplot2)

# Create a data frame
data <- data.frame(
  Category = c("% of Age-Up DARs with Motif", "% of all peaks with Motif"),
  Percentage = c(30.96, 18.56)
)

# Create the bar plot
p <- ggplot(data, aes(x = Category, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5) +
  labs(title = "Barplot of Given Percentages", y = "Percentage") +
  theme_minimal()

# Save the plot as a PDF
ggsave("../Figures/Figure2-IMN-IOL/barplot.pdf", plot = p, width = 6, height = 4)
p

options(repr.plot.width=2.5, repr.plot.height=3)
#nmat[which(nmat>100, arr.ind = T)] = 100
colnames(nmat) = gsub("OB-STR-CTX Inh IMN", "" ,colnames(nmat))
colnames(nmat) = gsub("motif", "" ,colnames(nmat))

heatmap <- pheatmap(nmat[nrow(nmat):1,],cutree_rows = 2, 
                    gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    #annotation_col = nannotation_col[,2, drop = F], 
                    cluster_cols = F, color = BuRd(20)[2:20],
                    main = "OB-STR-CTX Inh IMN\nAge-DAR motif enrichment",              
                    clustering_method = 'ward.D2',#gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    #annotation_colors = annotation_colors[c(2,1)],
                    display_numbers = snmat[nrow(nmat):1,], 
                    cluster_rows = F)

pdf("../Figures/Figure2-IMN-IOL/OB-STR_motifs.pdf", height = 3, width = 2.5)
heatmap
dev.off()

options(repr.plot.width=6, repr.plot.height=3.5)
#nmat[which(nmat>25, arr.ind = T)] = 25
nmat[which(nmat>100, arr.ind = T)] = 100

colnames(nmat) = gsub("OB-STR-CTX Inh IMN", "" ,colnames(nmat))

heatmap <- pheatmap(t(nmat[,]),cutree_cols = 2, 
                    gaps_col = c(grep("up", colnames(nmat))[1]-1),
                     
                    cluster_cols = T,
                    main = "OB-STR-CTX Inh IMN\nAge-DAR motif enrichment",              
                    clustering_method = 'ward.D2',#gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    
                    display_numbers = t(snmat[,]), 
                    cluster_rows = F)



library(khroma)
library(pheatmap)

# Define the color palette
BuRd_palette <- BuRd(100)  # Generate 100 colors

# Define breaks to force midpoint at 5
breaks_list <- seq(0, max(nmat), length.out = 100)  # Equal intervals
midpoint_index <- which.min(abs(breaks_list - 5))  # Find the closest break to 5
breaks_list[midpoint_index] <- 5  # Force 5 to be the middle value

# Generate heatmap
heatmap <- pheatmap(t(nmat), 
                    cutree_cols = 2, 
                    gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    cluster_cols = TRUE, 
                    color = BuRd_palette,
                    breaks = breaks_list,  # Apply custom breaks
                    main = "OB-STR-CTX Inh IMN\nAge-DAR motif enrichment",              
                    clustering_method = 'ward.D2',
                    display_numbers = t(snmat), 
                    cluster_rows = FALSE)


svg("../Figures/Figure2-IMN-IOL/motifs_all.svg",width=5.5, height =8)
heatmap
dev.off()

annotation_colors[2]

nannotation_col[,2, drop = F]


