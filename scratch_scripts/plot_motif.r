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


setwd(paste("../region_DARs_redone/",sep = ""))

files = list.files(".", pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)
files = paste(files, sep = "")
files = files[which(file.exists(files))]
length(files)
ind_files = files[grep("DG_Glut",files)] 


files = list.files("../combined_DARs_redone/", pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)
files = paste(files, sep = "")
files = files[which(file.exists(files))]
length(files)
files = files[grep("DG_Glut",files)] 
comb_files = files

files = c(ind_files , comb_files)

#files=files[-grep("within_pcdh",files)]
files

database = "motifs"
pval = 1e-10
gkey = c()
f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)

curr = curr[which(as.numeric(curr$`P-value`) < pval),]
head(curr)
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
  head(curr)

  curr = curr[1:min(15, nrow(curr)),]
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
      if(curr[paste(gkey[i]),"q-value (Benjamini)"]<0.01){
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
snmat = smat[which(rowMaxs(mat)>5),which(colMaxs(mat)>5)]
nmat = mat[which(rowMaxs(mat)>5),which(colMaxs(mat)>5)]
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

head(nmat)

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

nannotation_col["combined_down_DG_Glut",1] = "Combined"
nannotation_col["combined_up_DG_Glut",1] = "Combined"


rownames(nannotation_col) = gsub("_", " ", rownames(nannotation_col))
rownames(nannotation_col) = gsub(":", " ", rownames(nannotation_col))


colnames(nmat) = gsub("_", " ", colnames(nmat))
colnames(nmat) = gsub(":", " ", colnames(nmat))

snmat = snmat[,order(colnames(nmat))]
nmat = nmat[,order(colnames(nmat))]

colnames(nmat)[1:2] =c("down DG Glut Combined","up DG Glut Combined" )

snmat = snmat[,c(grep("down", colnames(nmat)), grep("up", colnames(nmat)))]
nmat = nmat[,c(grep("down", colnames(nmat)), grep("up", colnames(nmat)))]


rownames(nannotation_col)[nrow(nannotation_col):(nrow(nannotation_col)-1)] = c("up DG Glut Combined","down DG Glut Combined" )
nannotation_col[nrow(nannotation_col):(nrow(nannotation_col)-1),1] = c("Combined","Combined" )

annotation_colors

options(repr.plot.width=5, repr.plot.height=8.5)
nmat[which(nmat>50, arr.ind = T)] = 50

heatmap <- pheatmap(nmat[nrow(nmat):1,],cutree_rows = 5, gaps_row = c(5), gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    annotation_col = nannotation_col, 
                    cluster_cols = F, color = BuRd(100)[1:100],
                    main = "DG Glut Age-DAR motif enrichment",              
                    clustering_method = 'ward.D2',#gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    annotation_colors = annotation_colors,
                    display_numbers = snmat[nrow(nmat):1,], 
                    cluster_rows = T)



svg("../Figures/Figure5-DG/motifs.svg", height  =8.5, width = 5)
print(heatmap)
dev.off()

library(pheatmap)
nmat[which(nmat>100, arr.ind = T)] = 100

# Generate BuRd color palette with 100 colors
BuRd_palette <- colour("BuRd")(100)  # Generates a smooth BuRd gradient

# Define custom breaks to shift the scale
max_val <- max(nmat, na.rm = TRUE)  # Get the max value in the data

# Shift white downward and ensure red starts at 20
breaks <- c(seq(0, 20, length.out = 40),  # More resolution below 10 (shift white down)
            seq(20.1, 20, length.out = 30),  # White-to-red transition happens here
            seq(20.1, max_val, length.out = 30))  # Red dominates after 20

# Plot heatmap with the BuRd color scale and adjusted breaks
heatmap <- pheatmap(nmat[nrow(nmat):1,], 
                    cutree_rows = 5, 
                    gaps_row = c(5), 
                    annotation_col = nannotation_col, 
                    cluster_cols = FALSE, 
                    color = BuRd_palette,  # Use BuRd palette
                    breaks = breaks,  # Apply custom breaks
                    main = "DG Glut Age-DAR motif enrichment",              
                    clustering_method = 'ward.D2',
                    annotation_colors = annotation_colors,
                    display_numbers = snmat[nrow(nmat):1,], 
                    cluster_rows = TRUE)





