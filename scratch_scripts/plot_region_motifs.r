library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

setwd(paste("beds",sep = ""))


files = list.files(".", pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)
dir = "down"
files = files[grep("Astro-NT", files)]
files = files[grep(dir, files)]

files = files[which(file.exists(files))]
files



database = "motifs"
files = list.files(".", pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)
dir = "up"
files = files[grep("Oligo_NN", files)]
files = files[grep(dir, files)]
pval = 1e-10

files = files[which(file.exists(files))]
f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$`P-value` < pval),]

gkey = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  #  show(head(curr))
  curr = curr[which(curr$`P-value` < pval),]
  #cat(f, nrow(curr),"\n")

  gkey = c(gkey , paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1)))
 # gkey = c(gkey ,paste(curr[,1]))
  # cat((gkey))
}
gkey = gkey[which(!duplicated(gkey))]
cat(length(gkey))


# getting significance for each term + file
mat = matrix(0, nrow=length(gkey),ncol=length(files))
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

nmat = mat[which(rowMaxs(mat)>3), ]
nmat = nmat[,which(colMaxs(nmat)>3)]
nmat[which(nmat>30, arr.ind = T)] = 30

options(repr.plot.width=5, repr.plot.height=12
       )

nmat = mat[which(rowMaxs(mat)>15), ]
#nmat = nmat[,which(colMaxs(nmat)>5)]
nmat[which(nmat>1000, arr.ind = T)] = 1000

heatmap <- pheatmap(nmat, main = dir)
heatmap



options(repr.plot.width=10, repr.plot.height=20)

database = "motifs"

nmat = mat[,grep("_NN",colnames(mat))]
#nmat = nmat[,grep("_NN",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>17.5), ]
nmat = nmat[,which(colMaxs(nmat)>7.5)]
nmat[which(nmat>50, arr.ind = T)] = 50
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
dir = rep("down" , length(ct))
dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    dir = dir,
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=10, repr.plot.height=20)

database = "motifs"

nmat = mat[,grep("_Glut",colnames(mat))]
nmat = nmat[,grep("Female",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>20.5), ]
nmat = nmat[,which(colMaxs(nmat)>20.5)]
nmat[which(nmat>100, arr.ind = T)] = 100
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
dir = rep("down" , length(ct))
dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    dir = dir,
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=10, repr.plot.height=18)

database = "motifs"

nmat = mat[,grep("_NN",colnames(mat))]
nmat = nmat[,grep("Female",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>25.5), ]
nmat = nmat[,which(colMaxs(nmat)>20.5)]
nmat[which(nmat>100, arr.ind = T)] = 100
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
dir = rep("down" , length(ct))
dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    dir = dir,
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=10, repr.plot.height=18)

database = "motifs"

nmat = mat[,grep("_Glut",colnames(mat))]
nmat = nmat[,grep("Male",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>25.5), ]
nmat = nmat[,which(colMaxs(nmat)>20.5)]
nmat[which(nmat>100, arr.ind = T)] = 100
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
dir = rep("down" , length(ct))
dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    dir = dir,
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=10, repr.plot.height=18)

database = "motifs"

nmat = mat[,grep("_Gaba",colnames(mat))]
nmat = nmat[,grep("Male",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>25.5), ]
nmat = nmat[,which(colMaxs(nmat)>20.5)]
nmat[which(nmat>100, arr.ind = T)] = 100
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "",ct)
ct = gsub("_up", "",ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
dir = rep("down" , length(ct))
dir[grep("up", as.character(colnames(nmat)))] = "up"
nannotation_col = data.frame(
    dir = dir,
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=10, repr.plot.height=16)

database = "motifs"

nmat = mat[,grep("down",colnames(mat))]
nmat = nmat[,grep("_NN",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>17.5), ]
nmat = nmat[,which(colMaxs(nmat)>7.5)]
nmat[which(nmat>50, arr.ind = T)] = 50
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "", ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
nannotation_col = data.frame(
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=10, repr.plot.height=13)

database = "motifs"

nmat = mat[,grep("down",colnames(mat))]
nmat = nmat[,grep("_Glut",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>18.5), ]
nmat = nmat[,which(colMaxs(nmat)>18.5)]
nmat[which(nmat>100, arr.ind = T)] = 100
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_down", "", ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
nannotation_col = data.frame(
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=12, repr.plot.height=18)

database = "motifs"

nmat = mat[,grep("up",colnames(mat))]
nmat = nmat[,grep("_Glut",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>25.5), ]
nmat = nmat[,which(colMaxs(nmat)>7.5)]
nmat[which(nmat>100, arr.ind = T)] = 100
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_up", "", ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
nannotation_col = data.frame(
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=10, repr.plot.height=13)

database = "motifs"

nmat = mat[,grep("down",colnames(mat))]
nmat = nmat[,grep("_Gaba",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>20.5), ]
nmat = nmat[,which(colMaxs(nmat)>7.5)]
nmat[which(nmat>100, arr.ind = T)] = 100
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap

options(repr.plot.width=12, repr.plot.height=18)

database = "motifs"

nmat = mat[,grep("up",colnames(mat))]
nmat = nmat[,grep("_Gaba",colnames(nmat))]

nmat = nmat[which(rowMaxs(nmat)>25.5), ]
nmat = nmat[,which(colMaxs(nmat)>7.5)]
nmat[which(nmat>50, arr.ind = T)] = 50
ct = sapply(strsplit(as.character(colnames(nmat)), ":"), "[[", 1)
ct = gsub("_up", "", ct)
clade = sapply(strsplit(ct, "_"), function(x) tail(x, 1))
nannotation_col = data.frame(
    celltype = ct
  )
rownames(nannotation_col) = colnames(nmat) 
# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)
heatmap


