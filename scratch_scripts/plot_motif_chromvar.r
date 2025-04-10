library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

setwd("motif_chromvar/")

files = list.files(".", ".txt")

files

pval = 0.05
f = files[1]
curr = read.table(f, sep = " ")
curr = as.data.frame(curr)
curr = curr[which(curr$p_val_adj < pval),]
head(curr)

database = "motifs"
#dir = "up"
#files = files[grep("Oligo_NN", files)]
#files = files[grep(dir, files)]
pval = 1e-50

files = files[which(file.exists(files))]
f = files[1]
curr = read.table(f, sep = " ")
curr = as.data.frame(curr)
curr = curr[which(curr$p_val_adj < pval),]

gkey = curr$motif
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = read.table(f, sep = " ")
  curr = as.data.frame(curr)
  #  show(head(curr))
  curr = curr[which(curr$p_val_adj < pval),]
  #cat(f, nrow(curr),"\n")

  gkey = c(gkey , curr[,"motif"])
 # gkey = c(gkey ,paste(curr[,1]))
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
    curr = read.table(f, sep = " ")
    curr = as.data.frame(curr)
    rownames(curr) = curr$motif
    if(paste(gkey[i]) %in% paste(rownames(curr))) {
      # mat[x,y]  = -log10(curr[paste(gkey[i]), "P.Value"]+1e-100)
      mat[x,y]  = -log(curr[paste(gkey[i]),"p_val_adj"])*sign(curr[paste(gkey[i]),'avg_diff'])
      if(curr[paste(gkey[i]),"p_val_adj"]<0.01){
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
colnames(mat) = gsub("_18v2_motifs.txt", "",files)
colnames(mat) = gsub("./", "", colnames(mat))

mat[which(is.na(mat),arr.ind=T)] = 0 

clade = sapply(strsplit(colnames(mat), " "), function(x) tail(x, 1))
clade[grep("Gaba", clade )] = "Gaba"
clade[grep("Neur", clade )] = NA
               
               
               
annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 


options(repr.plot.width=15, repr.plot.height=12.5)

heatmap <- pheatmap(mat,display_numbers = smat,annotation_col=annotation_col,colorRampPalette(c("blue", "white", "red"))(100))


options(repr.plot.width=10, repr.plot.height=13
       )


minn = 15
rmin = 35
nsmat = smat[which(rowMaxs(mat)>rmin | rowMins(mat)<(-rmin)), which(colMaxs(mat)>minn | colMins(mat)<(-minn))]
nmat = mat[which(rowMaxs(mat)>rmin | rowMins(mat)<(-rmin)), which(colMaxs(mat)>minn | colMins(mat)<(-minn))]
nmat[which(nmat>70, arr.ind = T)] = 70
nmat[which(nmat<(-70), arr.ind = T)] = -70

heatmap <- pheatmap(nmat, display_numbers = nsmat,annotation_col=annotation_col,colorRampPalette(c("blue", "white", "red"))(100))


max(nmat)

mat_breaks

quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(nmat, n = 200)

heatmap <- pheatmap(nmat,breaks= mat_breaks,display_numbers = nsmat,colorRampPalette(c("blue", "white", "red"))(100))





