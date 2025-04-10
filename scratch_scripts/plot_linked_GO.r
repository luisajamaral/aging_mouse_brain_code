library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

setwd("results_logfc.15_adjpval.1/")

file.size(f)

setwd(".")
database = "reactome.txt"
dir = "DE_DAR"
#setwd(paste( dir,"/GO",sep = ""))
files = list.files(pattern = dir, path = ".")
files = paste(files, paste("/",database, sep =  ""), sep = "")
files = files[which(file.exists(files))]
files = files[grep("DE_DAR",files)]
pval = 1e-3



f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$Enrichment < pval),]
gkey = paste(curr[,c(2)])
all_files = files[1]
for(f in files[2:length(files)]) {
  if(!file.exists(f) | file.size(f)==0){
    cat(f, " NO FILE\n")
    next()
  }
  all_files = c(all_files, f)
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  #  show(head(curr))
  curr = curr[which(curr$Enrichment < pval),]
  #  curr = curr[which(curr$Genes.in.Term < 400),]
  # gkey = c(gkey , paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1)))
  gkey = c(gkey ,paste(curr[,2]))
  # cat((gkey))
}

gkey = gkey[which(!duplicated(gkey))]
cat(length(gkey))



gkey

  mat = matrix(0, nrow=length(gkey),ncol=length(all_files))
    x = 1
    for(i in 1:length(gkey)) {
      y = 1
      for(f in all_files) {
        curr = fread(f, sep = "\t", fill=T)
        curr = as.data.frame(curr)
        names = curr$Term
        which(duplicated(names))
        curr = curr[which(!duplicated(names)),]
        curr = curr[which(!is.na(curr$Term)),]
        curr = curr[which(curr$`Target Genes in Term`>2),]

        rownames(curr) = curr$Term
        if(paste(gkey[i]) %in% paste(rownames(curr))) {
          # mat[x,y]  = -log10(curr[paste(gkey[i]), "P.Value"]+1e-100)
          mat[x,y]  = -curr[paste(gkey[i]),"logP"]
          #  cat(mat[x,y], " ")
        } else {
          mat[x,y]  = NA
        }
        y= y+1
      }
      x=x+1
    }

 options(repr.plot.width=14
         , repr.plot.height=10)
    rownames(mat) = gkey
    cn = gsub(database, "", files)
    cn = gsub(".txt", "", cn)
    cn = gsub("/", "", cn)
    cn = gsub("_DE_", "",cn)
    cn = gsub("DAR", "",cn)


    colnames(mat) = cn
    mat[which(is.na(mat), arr.ind = T)] = 0
    nmat = mat
    nmat = nmat[which(rowMaxs(nmat)>3
                     ), ]
    nmat = nmat[order(rowMaxs(nmat), decreasing = T),]
   # if(nrow(nmat)>45) {
   #     nmat= nmat[1:45,]
   # }
    nmat = nmat[,which(colMaxs(nmat)>1.5)]
    nmat[which(nmat>35, arr.ind = T)] = 35

    clade = colnames(nmat)
    clade[grep("Neur", clade)] = "Gaba"
    clade[grep("Gaba", clade)] = "Gaba"
    clade[grep("Glut", clade)] = "Glut"
    clade[grep("NN", clade)] = "NN"
    clade[grep("IMN", clade)] = "IMN"


    dir = colnames(nmat)
    dir[grep("up", dir)] = "up"
    dir[grep("down", dir)] = "down"

    nannotation_col = data.frame(
        Clade = clade, Dir = dir
      )
    rownames(nannotation_col) = colnames(nmat) 

    nmat = nmat[,order(dir, clade)]

    display_matrix <- ifelse(nmat > 3, '*', '')

    heatmap <- pheatmap(nmat,cluster_cols = F,gaps_col = c(grep("up", colnames(nmat))[1]-1),
                        clustering_method = "average",,display_numbers = display_matrix, main = paste(database, "DAR linked genes"), annotation_col = nannotation_col)
    print(heatmap)


 options(repr.plot.width=15
         , repr.plot.height=10)
    rownames(mat) = gkey
    cn = gsub(database, "", files)
    cn = gsub(".txt", "", cn)
    cn = gsub("/", "", cn)
    cn = gsub("_DE_", "",cn)
    cn = gsub("GO", "",cn)


    colnames(mat) = cn
    mat[which(is.na(mat), arr.ind = T)] = 0
    nmat = mat
    nmat = nmat[which(rowMaxs(nmat)>3
                     ), ]
    nmat = nmat[order(rowMaxs(nmat), decreasing = T),]
   # if(nrow(nmat)>45) {
   #     nmat= nmat[1:45,]
   # }
    nmat = nmat[,which(colMaxs(nmat)>3)]
    nmat[which(nmat>35, arr.ind = T)] = 35

    clade = colnames(nmat)
    clade[grep("Neur", clade)] = "Gaba"
    clade[grep("Gaba", clade)] = "Gaba"
    clade[grep("Glut", clade)] = "Glut"
    clade[grep("NN", clade)] = "NN"
    clade[grep("IMN", clade)] = "IMN"


    dir = colnames(nmat)
    dir[grep("up", dir)] = "up"
    dir[grep("down", dir)] = "down"

    nannotation_col = data.frame(
        Clade = clade, Dir = dir
      )
    rownames(nannotation_col) = colnames(nmat) 

    nmat = nmat[,order(dir, clade)]

    display_matrix <- ifelse(nmat > 3, '*', '')

    heatmap <- pheatmap(nmat,cluster_cols = F,gaps_col = c(grep("up", colnames(nmat))[1]-1),
                        clustering_method = "average",,display_numbers = display_matrix, main = paste(database), annotation_col = nannotation_col)
    print(heatmap)





