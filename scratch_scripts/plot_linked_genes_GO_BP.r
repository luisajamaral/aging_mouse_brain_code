library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)


setwd("linked_peak_beds/")

files = list.files(path = ".", pattern = "BP.txt", recursive = T, full.names = T)

files[1]

files = files[grep("down", files)]
file_info <- file.info(files)
files = files[which(file_info$size>0)]

pval=0.05



f = files[1]
curr = read.table(f)
curr = as.data.frame(curr)
curr = curr[which(curr$`Adjusted.P.value`< pval),]
curr= curr[which(as.numeric(sapply(strsplit(as.character(curr$`Overlap`), "[/]"), `[`, 1))>10),]
head(curr)
gkey = curr$Name

for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = read.table(f)
  curr = as.data.frame(curr)
  curr = curr[which(curr$`Adjusted.P.value`< pval),]
  curr= curr[which(as.numeric(sapply(strsplit(as.character(curr$`Overlap`), "/"), `[`, 1))>10),]
  gkey = c(gkey ,paste(curr[,1]))
  # cat((gkey))
}

gkey = gkey[which(!duplicated(gkey))]
cat(length(gkey))
gkey

# getting significance for each term + file
mat = matrix(0, nrow=length(gkey),ncol=length(files))
smat = matrix('', nrow=length(gkey),ncol=length(files))

x = 1
for(i in 1:length(gkey)) {
  y = 1
  for(f in files) {
    curr = read.table(f)
    curr = as.data.frame(curr)

    rownames(curr) = curr$Term
    if(paste(gkey[i]) %in% paste(rownames(curr))) {
      # mat[x,y]  = -log10(curr[paste(gkey[i]), "P.Value"]+1e-100)
      mat[x,y]  = -log(curr[paste(gkey[i]),"P.value"])
      if(curr[paste(gkey[i]),"Adjusted.P.value"] < 0.05){
          smat[x,y] = '*'
        }
      #  cat(mat[x,y], " ")
    } else {
      mat[x,y]  = 0
    }
    y= y+1
  }
  x=x+1
}

options(repr.plot.width=10, repr.plot.height=18)

rownames(mat) = gkey
#cn = gsub(":Male_up_genome_ontology/repeats.genomeOntology.txt", "", files)
cn = gsub("_down_linked_GO_BP.txt", "", files)
cn = gsub("./", "", cn)
colnames(mat) = cn


clade = sapply(strsplit(colnames(mat), "_"), function(x) tail(x, 1))
clade[grep("Gaba", clade )] = "Gaba"
clade[grep("Dopa", clade )] = "Glut"

clade[grep("Neur", clade )] = NA
               
               
               
annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

options(repr.plot.width=15, repr.plot.height=15)

pheatmap(mat,annotation_col=annotation_col,colorRampPalette(c("blue", "white", "red"))(100),display_numbers = smat)

pheatmap(mat,annotation_col=annotation_col,display_numbers = smat)

options(repr.plot.width=15, repr.plot.height=13 )


minn = 10
rmin = 10
nsmat = smat[which(rowMaxs(mat)>rmin) , which(colMaxs(mat)>minn )]
nmat = mat[which(rowMaxs(mat)>rmin), which(colMaxs(mat)>minn)]
nmat[which(nmat>25, arr.ind = T)] = 25
#nmat[which(nmat<(-70), arr.ind = T)] = -70

heatmap <- pheatmap(nmat, display_numbers = nsmat,annotation_col=annotation_col)


nmat


