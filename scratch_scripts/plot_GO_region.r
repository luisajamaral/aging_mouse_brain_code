library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

setwd("~/projects/combined_all/female_RNA/DEG_results_region//")
database = "biological_process.txt"
dir = "up"
ct = "ENT"
#setwd(paste( dir,"/GO",sep = ""))
files = list.files(pattern = dir, path = ".")
files = paste(files, paste("/",database, sep =  ""), sep = "")
files = files[which(file.exists(files))]
files=files[grep(ct, files)]
files2 = list.files(pattern = ct, path = paste("../DEG_results_.01_.01/",dir,"/GO/", sep = ""), full.names = T)
files2 = paste(files2, paste("/",database, sep =  ""), sep = "")
files = c(files, files2)
pval = 1e-3
f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$Enrichment < pval),]
gkey = paste(curr[,c(2)])
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
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



files = list.files(pattern = dir, path = ".")


#setwd("~/ps-renlab2/projects/combined_all/female_RNA/DEG_results_region//")
for(ct in c("--ENT", "--HCP", "--AMY", "--NAC", "--RLP", "--HCA", "--FC", "--CP")) {
    database = "biological_process.txt"
    dir = "down"
   # ct = "--ENT"
    #setwd(paste( dir,"/GO",sep = ""))
    files = list.files(pattern = dir, path = ".")
    files = paste(files, paste("/",database, sep =  ""), sep = "")
    files = files[which(file.exists(files))]
    files=files[grep(ct, files)]
    files2 = list.files(pattern = ct, path = paste("../DEG_results_.01_.01/",dir,"/GO/", sep = ""), full.names = T)
    files2 = paste(files2, paste("/",database, sep =  ""), sep = "")
    files = c(files, files2)
    pval = 1e-2
    f = files[1]
    curr = fread(f, sep = "\t",fill = T)
    curr = as.data.frame(curr)
    curr = curr[which(curr$Enrichment < pval),]
    gkey = paste(curr[,c(2)])
    for(f in files[2:length(files)]) {
      if(!file.exists(f)){
       # cat(f, " NO FILE\n")
        next()
      }
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
    if(length(gkey)<2){
        next
    }

    files=files[-length(files)]


    # getting significance for each term + file
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
        } else {
          mat[x,y]  = NA
        }
        y= y+1
      }
      x=x+1
    }

    options(repr.plot.width=15, repr.plot.height=10)
    rownames(mat) = gkey
    cn = gsub(database, "", files)
    cn = gsub(".txt", "", cn)
    cn = gsub("/", "", cn)
    cn = gsub("..DEG_results_.01_.01upGO", "combined_",cn)
    cn = gsub("..DEG_results_.01_.01downGO", "combined_",cn)


    colnames(mat) = cn
    mat[which(is.na(mat), arr.ind = T)] = 0
    nmat = mat
    nmat = nmat[which(rowMaxs(nmat)>2), ]
    nmat = nmat[order(rowMaxs(nmat), decreasing = T),]
    if(nrow(nmat)>45) {
        nmat= nmat[1:45,]
    }
    #nmat = nmat[,which(colMaxs(nmat)>7)]
    nmat[which(nmat>35, arr.ind = T)] = 35

    clade = colnames(nmat)
    clade[grep("Neur", clade)] = "Gaba"
    clade[grep("Gaba", clade)] = "Gaba"
    clade[grep("Glut", clade)] = "Glut"
    clade[grep("NN", clade)] = "NN"
    clade[grep("IMN", clade)] = "IMN"

    nannotation_col = data.frame(
        Clade = clade
      )
    rownames(nannotation_col) = colnames(nmat) 

    heatmap <- pheatmap(nmat, main = paste(gsub("--", "",ct), dir, database), annotation_col = nannotation_col)
    print(heatmap)
}

#setwd("~/ps-renlab2/projects/combined_all/female_RNA/DEG_results_region//")
for(ct in c("--HCP", "--AMY", "--NAC", "--RLP", "--HCA", "--FC", "--ENT", "--CP")) {
    database = "kegg.txt"
    dir = "up"
    #ct = "--HCP"
    #setwd(paste( dir,"/GO",sep = ""))
    files = list.files(pattern = dir, path = ".")
    files = paste(files, paste("/",database, sep =  ""), sep = "")
    files = files[which(file.exists(files))]
    files=files[grep(ct, files)]
    files2 = list.files(pattern = ct, path = paste("../DEG_results_.01_.01/",dir,"/GO/", sep = ""), full.names = T)
    files2 = paste(files2, paste("/",database, sep =  ""), sep = "")
    files = c(files, files2)
    pval = 1e-4
    f = files[1]
    curr = fread(f, sep = "\t",fill = T)
    curr = as.data.frame(curr)
    curr = curr[which(curr$Enrichment < pval),]
    gkey = paste(curr[,c(2)])
    for(f in files[2:length(files)]) {
      if(!file.exists(f)){
       # cat(f, " NO FILE\n")
        next()
      }
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
    #cat(length(gkey))
    if(length(gkey)<2){
        next
    }

    files=files[-length(files)]


    # getting significance for each term + file
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
        } else {
          mat[x,y]  = NA
        }
        y= y+1
      }
      x=x+1
    }

    options(repr.plot.width=8, repr.plot.height=5)
    rownames(mat) = gkey
    cn = gsub(database, "", files)
    cn = gsub(".txt", "", cn)
    cn = gsub("/", "", cn)
    cn = gsub("..DEG_results_.01_.01upGO", "combined_",cn)
    cn = gsub("..DEG_results_.01_.01downGO", "combined_",cn)


    colnames(mat) = cn
    mat[which(is.na(mat), arr.ind = T)] = 0
    nmat = mat
    nmat = nmat[which(rowMaxs(nmat)>4), ]
    #nmat = nmat[,which(colMaxs(nmat)>7)]
    nmat[which(nmat>15, arr.ind = T)] = 15

    clade = colnames(nmat)
    clade[grep("Neur", clade)] = "Gaba"
    clade[grep("Gaba", clade)] = "Gaba"
    clade[grep("Glut", clade)] = "Glut"
    clade[grep("NN", clade)] = "NN"
    clade[grep("IMN", clade)] = "IMN"

    nannotation_col = data.frame(
        Clade = clade
      )
    rownames(nannotation_col) = colnames(nmat) 

    heatmap <- pheatmap(nmat, main = paste(gsub("--", "",ct), dir, database), annotation_col = nannotation_col)
    print(heatmap)
}





options(repr.plot.width=14, repr.plot.height=8)
heatmap

options(repr.plot.width=7, repr.plot.height=5)

nmat = mat
nmat = nmat[,-grep("combined", colnames(nmat))]
nmat = nmat[which(rowMaxs(nmat)>2), ]
#nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(nmat>15, arr.ind = T)] = 15
     
heatmap <- pheatmap(nmat, main = database)


options(repr.plot.width=6, repr.plot.height=10)

nmat = mat
#
#nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(is.na(nmat), arr.ind = T)] = 0
nmat = nmat[which(rowMaxs(nmat[,grep("--", colnames(nmat))])>10), ]
nmat[which(nmat>50, arr.ind = T)] = 50
     
heatmap <- pheatmap(nmat, main = database)


options(repr.plot.width=8, repr.plot.height=17)

nmat = mat
#
#nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(is.na(nmat), arr.ind = T)] = 0
nmat = nmat[which(rowMaxs(nmat)>5.3), ]
nmat[which(nmat>25, arr.ind = T)] = 25
     
heatmap <- pheatmap(nmat, main = database)


options(repr.plot.width=12, repr.plot.height=10)

nmat = mat
nmat = nmat[which(rowMaxs(nmat)>9), ]
#nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(nmat>15, arr.ind = T)] = 15
               
clade = sapply(strsplit(colnames(nmat), "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"

nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
heatmap <- pheatmap(nmat, main = database, annotation_col = nannotation_col)


heatmap

options(repr.plot.width=14, repr.plot.height=10)


nmat = mat
nmat = nmat[which(rowMaxs(nmat)>13), ]
nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(nmat>25, arr.ind = T)] = 25
               
clade = sapply(strsplit(colnames(nmat), "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"

nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
heatmap <- pheatmap(nmat, main = database, annotation_col = nannotation_col)


options(repr.plot.width=14, repr.plot.height=10)


nmat = mat
nmat = nmat[which(rowMaxs(nmat)>13), ]
nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(nmat>25, arr.ind = T)] = 25
               
clade = sapply(strsplit(colnames(nmat), "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"

nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
heatmap <- pheatmap(nmat, main = database, annotation_col = nannotation_col)



options(repr.plot.width=14, repr.plot.height=10)
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

annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

nmat = mat
nmat = nmat[which(rowMaxs(nmat)>35), ]
nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(nmat>55, arr.ind = T)] = 55
               
clade = sapply(strsplit(colnames(nmat), "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"

nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
heatmap <- pheatmap(nmat, main = database, annotation_col = nannotation_col)


setwd("~/ps-renlab2/projects/combined_all/female_RNA/DEG_results_.01_.01/")
pval = 1e-15
database = "biological_process.txt"
dir = "down"
setwd(paste( dir,"/GO",sep = ""))
files = list.files(pattern = "", path = ".")
files = paste(files, paste("/",database, sep =  ""), sep = "")
files = files[which(file.exists(files))]
files = files[-grep("doublet", files)]

f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$Enrichment < pval),]
gkey = paste(curr[,c(2)])
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
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

# getting significance for each term + file
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
    } else {
      mat[x,y]  = NA
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

annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

nmat = mat
nmat = nmat[which(rowMaxs(nmat)>9), ]
nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(nmat>50, arr.ind = T)] = 50
               
clade = sapply(strsplit(colnames(nmat), "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"

nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
heatmap <- pheatmap(nmat, main = database, annotation_col = nannotation_col)



options(repr.plot.width=18, repr.plot.height=12)
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

annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

nmat = mat
nmat = nmat[which(rowMaxs(nmat)>9), ]
nmat = nmat[,which(colMaxs(nmat)>7)]
nmat[which(nmat>50, arr.ind = T)] = 50
               
clade = sapply(strsplit(colnames(nmat), "_"), function(x) tail(x, 1))
clade[which(clade%in%"Neur")] = "Gaba"
clade[which(clade%in%"Gaba-Chol")] = "Gaba"
clade[which(clade%in%"Glut-Sero")] = "Glut"
clade[which(clade%in%"Dopa")] = "Glut"
clade[which(clade%in%"IMN")] = "NN"

nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
heatmap <- pheatmap(nmat, main = database, annotation_col = nannotation_col)


options(repr.plot.width=20, repr.plot.height=13)
dir = "down"
rownames(mat) = gkey
cn = gsub(paste("_2mo_vs_18mo.edger_",dir,"_motifs/", sep = ""), "", files)
cn = gsub(database, "", cn)
cn = gsub(".txt", "", cn)

colnames(mat) = cn

clade = key[colnames(mat), "L2.1"]
annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

nmat = mat[which(rowMaxs(mat)>15.5), ]
#nmat = nmat[,which(colMaxs(nmat)>5.5)]
nmat[which(nmat>20, arr.ind = T)] = 20
clade = key[colnames(nmat), "L2.1"]
nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
heatmap <- pheatmap(nmat, main = database)

options(repr.plot.width=20, repr.plot.height=13)

rownames(mat) = gkey
cn = gsub(paste("_2mo_vs_18mo.edger_",dir,"_motifs/", sep = ""), "", files)
cn = gsub(database, "", cn)
cn = gsub(".txt", "", cn)

colnames(mat) = cn

clade = key[colnames(mat), "L2.1"]
annotation_col = data.frame(
    Clade = clade
  )
rownames(annotation_col) = colnames(mat) 

nmat = mat[which(rowMaxs(mat)>20.5), ]
nmat = nmat[,which(colMaxs(nmat)>7.5)]
nmat[which(nmat>20, arr.ind = T)] = 20
clade = key[colnames(nmat), "L2.1"]
nannotation_col = data.frame(
    Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 

heatmap <- pheatmap(nmat, main = database)


database = gsub(".txt", "", database)

# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)


database = gsub(".txt", "", database)

# Create the heatmap plot
heatmap <- pheatmap(nmat, annotation_col = nannotation_col, main = database)

# Save the plot as a PNG image
ggsave(filename = paste(database, "_", dir, "_", "heat.png", sep = ""), plot = heatmap, width = 14, height = 11)


pheatmap(mat, annotation_col = annotation_col, main = paste(database, dir))
pheatmap(nmat, annotation_col = nannotation_col, main = paste(database, dir))





