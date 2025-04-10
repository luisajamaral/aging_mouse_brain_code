library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)









reg = "NAC"

setwd(paste("~/projects/combined_all/female_RNA/",reg,"_region_edgeR/", sep = ""))

database = "GO_BP"
dir = "up"
files = list.files(pattern = "", path = ".")
files = files[grep(database, files)]
files = files[grep(dir, files)]
files
pval = .0001

f = files[1]
curr = read.table(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
head(curr)
curr = curr[which(curr$P.value < pval),]
curr= curr[which(as.numeric(sapply(strsplit(as.character(curr$Overlap), "/"), "[[", 1))>10),]
gkey = paste(curr[,c(1)])

gkey


for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  cat(f)
  curr = read.table(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  #  show(head(curr))
    
  #print(plotEnrich(curr, showTerms = 20, numChar = 35, y = "Count", orderBy = "P.value", title = paste(f, reg,"GO_BP")))
    

  curr = curr[which(curr$P.value < pval),]
  curr= curr[which(as.numeric(sapply(strsplit(as.character(curr$Overlap), "/"), "[[", 1))>10),]
  
  #  curr = curr[which(curr$Genes.in.Term < 400),]
  # gkey = c(gkey , paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1)))
  gkey = c(gkey ,paste(curr[,1]))
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
    curr = read.table(f, sep = "\t", fill=T)
    curr = as.data.frame(curr)
    names = curr$Term
    which(duplicated(names))
    curr = curr[which(!duplicated(names)),]
    curr = curr[which(!is.na(curr$Term)),]

    rownames(curr) = curr$Term
    if(paste(gkey[i]) %in% paste(rownames(curr))) {
      # mat[x,y]  = -log10(curr[paste(gkey[i]), "P.Value"]+1e-100)
      mat[x,y]  = -log10(curr[paste(gkey[i]),"P.value"]+0.0000000001)
      #  cat(mat[x,y], " ")
    } else {
      mat[x,y]  = 0
    }
    y= y+1
  }
  x=x+1
}


rownames(mat) = gkey
cn = gsub(database, "", files)
cn = gsub(".txt", "", cn)
cn = gsub("/", "", cn)


colnames(mat) = cn
height = 9
if(nrow(mat)<15){
    height = 5
}

options(repr.plot.width=10, repr.plot.height=height)


clade = rep("Gaba",length(cn))
clade[grep("NN",cn)] = "NN"
clade[grep("Glut",cn)] = "Glut"
clade[grep("Dopa",cn)] = "Glut"
clade[grep("IMN",cn)] = "Glut"

annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

heatmap <- pheatmap(mat, main = paste(reg, dir,"in aging\n", database), annotation_col = annotation_col)


options(repr.plot.width=12, repr.plot.height=12)
heatmap



clade

library(enrichR)

enriched_up[[1]]

options(repr.plot.width=7, repr.plot.height=5)

plotEnrich(curr, showTerms = 20, numChar = 35, y = "Count", orderBy = "P.value", title = paste(reg,"GO_BP"))


reg="HCP"


