library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

t = "Female.CHN.DMG"
setwd("../Female.CHN.DMG/")

#setwd("~/ps-renlab2/projects/combined_all/female_RNA/DEG_results_.01_.01/")
database = "GO_BP"
dir = "Hypo"
files = list.files(pattern = database, path = ".")
files = files[grep(dir, files)]
files = files[grep(".txt", files)]

files = files[which(file.exists(files))]
pval = .05
f = files[1]
curr = read.table(f, sep = "_",fill = T)
curr = as.data.frame(curr)



curr = curr[which(curr$Adjusted.P.value < pval),]

gkey = paste(curr[,1])
gkey = c()
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = read.table(f, sep = "_",fill = T)
  curr = as.data.frame(curr)
  #  show(head(curr))
  curr = curr[which(curr$Adjusted.P.value < pval),]
  gkey = c(gkey ,paste(curr[,1]))
  # cat((gkey))
}

gkey = gkey[which(!duplicated(gkey))]
cat(length(gkey))




# getting significance for each term + file
mat = matrix(0, nrow=length(gkey),ncol=length(files))
smat = matrix('', nrow=length(gkey),ncol=length(files))

x = 1
for(i in 1:length(gkey)) {
  y = 1
  for(f in files) {
    curr = read.table(f, sep = "_", fill=T)
    curr = as.data.frame(curr)
    names = curr[,2]
    which(duplicated(names))
    curr = curr[which(!duplicated(names)),]
    curr = curr[which(!is.na(curr$Term)),]

    rownames(curr) = curr$Term
    if(paste(gkey[i]) %in% paste(rownames(curr))) {
       mat[x,y]  = -log10(curr[paste(gkey[i]), "Adjusted.P.value"]+1e-100)
      if(curr[paste(gkey[i]),"Adjusted.P.value"]<0.05) {
            smat[x,y]= '*'
        }
      #  cat(mat[x,y], " ")
    } else {
      mat[x,y]  = 0
    }
    y= y+1
  }
  x=x+1
}

options(repr.plot.width=14, repr.plot.height=11)
rownames(mat) = gkey
cn = gsub(database, "", files)
cn = gsub(".txt", "", cn)
cn = gsub("/", "", cn)


colnames(mat) = cn

clade = sapply(strsplit(cn, "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"


pheatmap(mat)

files

options(repr.plot.width=11, repr.plot.height=7
       )


nmat = mat
nsmat = smat[which(rowMaxs(nmat)>2),which(colMaxs(nmat)>1) ]

nmat = nmat[which(rowMaxs(nmat)>2), which(colMaxs(nmat)>1)]
nmat[which(nmat>50, arr.ind = T)] = 50
               
clade = sapply(strsplit(colnames(nmat), "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"
               
major = colnames(nmat)
major[grep("Gaba", major)] = "Gaba"             
major[grep("Glut", major)] = "Glut"             
major[grep("NN", major)] = "NN"             
major[grep("IMN", major)] = "IMN"  


nannotation_col = data.frame(
    Direction = clade, Clade= major
  )

colnames(nmat) = gsub(t, '', colnames(nmat))
colnames(nmat) = gsub("_$", '', colnames(nmat))

rownames(nannotation_col) = colnames(nmat) 

heatmap <- pheatmap(nmat, main = paste(database, dir, t), annotation_col = nannotation_col,display_numbers = nsmat)


pdf(paste(dir, database,"heatmap.pdf", sep = "_"), width = 12, height = 8)
print(heatmap)
dev.off()




