library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

database = "motifs"
setwd(paste("female_diff/motif_results/",sep = ""))
files = list.files(".", pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)
files = paste(files, sep = "")
files = files[which(file.exists(files))]
f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$`P-value` < 1e-25),]

gkey = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  #  show(head(curr))
  curr = curr[which(curr$`P-value` < 1e-15),]
  #cat(f, nrow(curr),"\n")

  gkey = c(gkey , paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1)))
 # gkey = c(gkey ,paste(curr[,1]))
  # cat((gkey))
}
gkey = gkey[which(!duplicated(gkey))]
cat(length(gkey))
fgkey = gkey

head(curr)

# getting significance for each term + file
mat = matrix(0, nrow=length(gkey),ncol=length(files))
smat = matrix('', nrow=length(gkey),ncol=length(files))

x = 1
for(i in 1:length(gkey)) {
  y = 1
  for(f in files) {
    if (length(grep("Female", f)) == 1) {
        f = paste("../../female_diff/motif_results/", f, sep = "")
    }

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
      #  cat(mat[x,y], " ")
    } else {
      mat[x,y]  = NA
    }
    y= y+1
  }
  x=x+1
}



rownames(mat) = gkey
colnames(mat) = gsub("_motifs_bg/knownResults.txt", "",files)
colnames(mat) = gsub("./", "", colnames(mat))
#colnames(mat) = gsub(paste(dir, "_homer_bg/", sep = ""), "", colnames(mat))
#colnames(mat) = gsub(database, "", colnames(mat))
#colnames(mat) = gsub("_", "", colnames(mat))
ct = sapply(strsplit(as.character(colnames(mat)), ":"), "[[", 1)


clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))

annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

nmat = mat[which(rowMaxs(mat)>17.5),which(colMaxs(mat)>7.5)]
nsmat = smat[which(rowMaxs(mat)>17.5),which(colMaxs(mat)>7.5)]
               
nmat[which(nmat>30, arr.ind = T)] = 30
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
dir = rep("down" , length(ct))
dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    Clade = clade,
    Dir = dir 
  )
rownames(nannotation_col) = colnames(nmat) 

options(repr.plot.width=21, repr.plot.height=28)

#database = gsub(".txt", "", database)
database = "motifs"
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database, display_numbers = nsmat)
heatmap

options(repr.plot.width=21, repr.plot.height=28)

umat = mat[,grep("up",colnames(mat))]
usmat = smat[,grep("up",colnames(mat))]

nmat = umat[which(rowMaxs(umat)>20),which(colMaxs(umat)>10)]
nsmat = usmat[which(rowMaxs(umat)>20),which(colMaxs(umat)>10)]
               
nmat[which(nmat>50, arr.ind = T)] = 50
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
dir = rep("down" , length(ct))
dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    Clade = clade,
    Dir = dir 
  )
rownames(nannotation_col) = colnames(nmat) 
               
options(repr.plot.width=10, repr.plot.height=25)

heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database, display_numbers = nsmat)
heatmap

options(repr.plot.width=18, repr.plot.height=22)

umat = mat[,grep("up",colnames(mat))]
usmat = smat[,grep("up",colnames(mat))]

nmat = umat[which(rowMaxs(umat)>20),which(colMaxs(umat)>10)]
nsmat = usmat[which(rowMaxs(umat)>20),which(colMaxs(umat)>10)]
               
nmat[which(nmat>50, arr.ind = T)] = 50
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
dir = rep("down" , length(ct))
dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    Clade = clade,
    Dir = dir 
  )
rownames(nannotation_col) = colnames(nmat) 
               
options(repr.plot.width=8, repr.plot.height=22)

heatmap <- pheatmap(nmat, annotation_col = nannotation_col, display_numbers = nsmat, main = "motifs female up")
heatmap

rownames(mat)

options(repr.plot.width=8
        , repr.plot.height=7)

database = "motifs"
nmat = mat
#cs = c(grep("Oligo", colnames(nmat)), grep("Micro", colnames(nmat)), grep("D12", colnames(nmat)), grep("L2-3", colnames(nmat)))
#nmat = nmat[, cs]
nmat = nmat[, grep("up", colnames(nmat))]
nmat = nmat[, grep(":", colnames(nmat))]

#nmat = mat[,grep("_NN",colnames(mat))]
#nmat = nmat[,grep("_NN",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>30.5), ]
nmat = nmat[,which(colMaxs(nmat)>30.5)]
fir = substr(rownames(nmat), 1, 5)
nmat = nmat[-which(duplicated(fir)),]

type = sapply(strsplit(rownames(nmat), "[(]"), `[`, 2)
nmat = nmat[-which(duplicated(type)),]


nmat[which(nmat>50, arr.ind = T)] = 50
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
#dir = rep("down" , length(ct))
#dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = "Motifs enriched in age-down")
heatmap

fir

type

options(repr.plot.width=8
        , repr.plot.height=10)

database = "motifs"
nmat = mat
#cs = c(grep("Oligo", colnames(nmat)), grep("Micro", colnames(nmat)), grep("D12", colnames(nmat)), grep("L2-3", colnames(nmat)))
#nmat = nmat[, cs]
nmat = nmat[, grep("up", colnames(nmat))]
snmat = smat[, grep("up", colnames(mat))]
#nmat = mat[,grep("_NN",colnames(mat))]
#nmat = nmat[,grep("_NN",colnames(nmat))]

snmat = snmat[which(rowMaxs(nmat)>15),which(colMaxs(nmat)>10)]

nmat = nmat[which(rowMaxs(nmat)>15),which(colMaxs(nmat)>10)]

fir = substr(rownames(nmat), 1, 5)
nmat = nmat[-which(duplicated(fir)),]
snmat = snmat[-which(duplicated(fir)),]

type = sapply(strsplit(rownames(nmat), "[(]"), `[`, 2)

rm = as.data.frame(rowMaxs(nmat))
rm$motif = rownames(rm)
rm = rm[order(rm[,1], decreasing = T),]
rm$type = sapply(strsplit(rownames(rm), "[(]"), `[`, 2)

library(dplyr)

filtered_df <- rm %>%
  group_by(type) %>%
  filter(row_number() <= 3 | n() < 4) %>%
  ungroup()

# View the filtered data frame
print(filtered_df$motif)

snmat = snmat[which(rownames(nmat)%in% filtered_df$motif),]
nmat = nmat[which(rownames(nmat)%in% filtered_df$motif),]

nmat[which(nmat>65, arr.ind = T)] = 65
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
#dir = rep("down" , length(ct))
#dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = "Motifs enriched in age-down", display_numbers = snmat)
heatmap



rm


