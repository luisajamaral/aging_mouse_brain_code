library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

files = list.files(".", pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)
files = paste(files, sep = "")
files = files[which(file.exists(files))]
length(files)

files = files[grep("adj", files)]
length(files)

files = files[grep("promoters", files)]
length(files)

pval = 1e-150

f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$`q-value (Benjamini)` < pval),]

gkey = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  #  show(head(curr))
  curr = curr[which(curr$`q-value (Benjamini)` < pval),]
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
     # mat[x,y]  = -log10(curr[paste(gkey[i]), "q-value (Benjamini)"]+1e-1000)
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
colnames(mat) = gsub("./combined_diff_", "", colnames(mat))
colnames(mat) = gsub("adjp_0.05/", "", colnames(mat))
colnames(mat) = gsub("promoters_", "", colnames(mat))
colnames(mat) = gsub("enhancers_", "", colnames(mat))

colnames(mat) = gsub("peaks_", "", colnames(mat))

snmat = smat[which(rowMaxs(mat)>150),which(colMaxs(mat)>25)]
nmat = mat[which(rowMaxs(mat)>150),which(colMaxs(mat)>25)]
#nmat[which(nmat>100, arr.ind = T)] = 100

dir = rep("down" , ncol(nmat))
dir[grep("up", as.character(colnames(nmat)))] = "up"

clade = rep("Gaba" , ncol(nmat))
clade[grep("Glut", as.character(colnames(nmat)))] = "Glut"
clade[grep("NN", as.character(colnames(nmat)))] = "NN"
clade[grep("IMN", as.character(colnames(nmat)))] = "IMN"

clade = factor(clade, levels = c("Gaba", "Glut", "NN", "IMN"))

#type = rep("Promoter" , ncol(nmat))
#type[grep("enhancer", as.character(colnames(nmat)))] = "Enhancer"


nannotation_col = data.frame(
    Dir = dir, Clade = clade#, Type = type
  )
rownames(nannotation_col) = colnames(nmat) 

#snmat = snmat[,order(dir,clade,colnames(nmat))]
#nmat = nmat[,order(dir,clade,colnames(nmat) )]

snmat = snmat[,order(dir,clade)]
nmat = nmat[,order(dir,clade)]

#snmat = snmat[,order(dir, clade)]
#nmat = nmat[,order(clade )]





options(repr.plot.width=8, repr.plot.height=8)
#snmat = smat[which(rowMaxs(mat)>75),which(colMaxs(mat)>15)]
#nmat = mat[which(rowMaxs(mat)>75),which(colMaxs(mat)>15)]

nmat[which(nmat>250, arr.ind = T)] = 250


heatmap <- pheatmap(nmat[,],  clustering_method = "ward.D2", gaps_col =  c(grep("up", colnames(nmat))[1]-1), 
                    cluster_cols = F, 
                     annotation_col = nannotation_col,
                    display_numbers = snmat, 
                    cluster_rows = T, main = "promoter-TSS age DAR top motifs")

#pdf("combined_motifs_nop.pdf", width=15, height = 18)
#print(heatmap)
#dev.off()






