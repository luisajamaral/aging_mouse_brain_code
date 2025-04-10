library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

database = "reactome.txt"
setwd(paste("../combined_DARs_redone/",sep = ""))


files = list.files(".", pattern = database, full.names = TRUE, recursive = TRUE)
files = paste(files, sep = "")
files = files[which(file.exists(files))]
length(files)

pval = 1e-5

f = files[1]
curr = fread(f, sep = "\t", fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$Enrichment< pval),]
curr= curr[which(curr$`Genes in Term`>5),]
head(curr)
gkey = curr$Term
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  curr = curr[which(curr$Enrichment< pval),]
  curr= curr[which(curr$`Genes in Term`>5),]
  curr = curr[1:min(20, nrow(curr)),]
  cat(f , length(paste(curr[,2])), "\n")
  gkey = c(gkey ,paste(curr[,2]))
  # cat((gkey))
}

gkey = gkey[which(!duplicated(gkey))]
cat(length(gkey))



# getting significance for each term + file
smat = matrix('', nrow=length(gkey),ncol=length(files))
mat = matrix(0, nrow=length(gkey),ncol=length(files))
x = 1
for(i in 1:length(gkey)) {
  y = 1
  for(f in files) {
    curr = fread(f, sep = "\t", fill=T)
    curr = as.data.frame(curr)
    names = curr$Term
    which(duplicated(names))
    curr = curr[which(!duplicated(names)),]
    curr = curr[which(!is.na(curr$Term)),]

    rownames(curr) = curr$Term
    if(paste(gkey[i]) %in% paste(rownames(curr))) {
      # mat[x,y]  = -log10(curr[paste(gkey[i]), "P.Value"]+1e-100)
      mat[x,y]  = -curr[paste(gkey[i]),"logP"]
      #  cat(mat[x,y], " ")
    if(curr[paste(gkey[i]),"Enrichment"]<0.01){
        smat[x,y] = '*'
      }
    } else {
      mat[x,y]  = 0
    }
    y= y+1
  }
  x=x+1
}


rownames(mat) = gkey
colnames(mat) = gsub("_go", "",files)
colnames(mat) = gsub("/", "",colnames(mat))

colnames(mat) = gsub(database, "",colnames(mat))

colnames(mat) = gsub("./", "", colnames(mat))
colnames(mat) = gsub("peaks_", "", colnames(mat))
colnames(mat) = gsub("[.]", "", colnames(mat))


snmat = smat[which(rowMaxs(mat)>15),which(colMaxs(mat)>15)]
nmat = mat[which(rowMaxs(mat)>15
                ),which(colMaxs(mat)>15)]
#nmat[which(nmat>100, arr.ind = T)] = 100

dir = rep("down" , ncol(nmat))
dir[grep("up", as.character(colnames(nmat)))] = "up"


clade= sapply(sapply(strsplit(colnames(nmat), ":"), `[`, 1), function(x) {
  parts <- strsplit(x, "_")[[1]]  # Split the string by space
  last_element <- parts[length(parts)]  # Get the last element
  return(last_element)
})
clade[which(clade == "Neur")] = "Gaba"
nannotation_col = data.frame(
    Dir = dir,  Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
               
nmat = nmat[,order(dir,  clade)]
snmat = snmat[,order(dir,  clade)]

options(repr.plot.width=12, repr.plot.height=12)

breaks <- c(seq(min(nmat), 5, length.out = 100), seq(5.01, max(nmat), length.out = 100))
colors <- c(colorRampPalette(c("gray", "yellow"))(100), colorRampPalette(c("yellow1", "red2"))(100))

heatmap <- pheatmap(nmat, 
                    annotation_col = nannotation_col, 
                    cluster_cols = F, 
                    clustering_method = 'ward.D2',
                    main = database, 
                    display_numbers = snmat, 
                    cluster_rows = T)
pdf("Reactome_DAR_GO.pdf", width =14, height = 12)
print(heatmap)
dev.off()


