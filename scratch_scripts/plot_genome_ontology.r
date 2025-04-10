library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

setwd("genome_ontology/")

database = "basic"
dir = "down"
files = list.files(pattern = ".", path = ".", recursive = T, full.names = T)
files = files[grep(database, files)]
files = files[grep(dir, files)]
files=files[-grep('2000', files)]

pval = .01

file_info <- file.info(files)
file_info$size
files = files[which(file_info$size>0)]


f = files[1]
curr = fread(f, sep = "\t", fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$`P-value(vs. Control)`< pval),]
curr= curr[which(curr$`Overlap(#peaks)`>2),]
head(curr)
gkey = curr$Name

gkey

for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  curr = curr[which(curr$`P-value(vs. Control)`< pval),]
  curr= curr[which(curr$`Overlap(#peaks)`>2),]
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
    curr = fread(f, sep = "\t", fill=T)
    curr = as.data.frame(curr)
    names = curr$Name
    which(duplicated(names))
    curr = curr[which(!duplicated(names)),]
    curr = curr[which(!is.na(curr$Name)),]

    rownames(curr) = curr$Name
    if(paste(gkey[i]) %in% paste(rownames(curr))) {
      # mat[x,y]  = -log10(curr[paste(gkey[i]), "P.Value"]+1e-100)
      mat[x,y]  = -curr[paste(gkey[i]),"Log P-value(+ underrepresented, vs. Control)"]
      if(curr[paste(gkey[i]),'P-value(vs. Control)']<0.01){
         smat[x,y]='*' 
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
cn = gsub("_motifs_bg/basic.genomeOntology.txt", "", files)

cn = gsub("diff_peaks_", "", cn)

cn = gsub("./", "", cn)
colnames(mat) = cn

options(repr.plot.width=15, repr.plot.height=15)

pheatmap(mat)

options(repr.plot.width=8, repr.plot.height=8)
#png(paste(database, "_heatmaps.png", sep = ""),width = 7,height = 7)

for(reg in c("ENT", "FC", "AMY", "RLP", "NAC", "HCA", "HCP", ":CP")){

    nmat = mat[,grep(reg, colnames(mat))]
    nsmat = smat[,grep(reg, colnames(mat))]
    nsmat = nsmat[which(rowMaxs(nmat)>5),which(colMaxs(nmat)>3)]
    nmat = nmat[which(rowMaxs(nmat)>5),which(colMaxs(nmat)>3)]
    nmat[which(nmat>30, arr.ind = T)] = 30
    nmat[which(nmat<(-30), arr.ind = T)] = -30


    region = sapply(strsplit(colnames(nmat), ":"), function(x) tail(x, 1))
    clade = rep("NN",ncol(nmat))
    clade[grep("Glut", colnames(nmat))] = "Glut"
    clade[grep("Gaba", colnames(nmat))] = "Gaba"
    clade[grep("IMN", colnames(nmat))] = "IMN"
    clade[grep("Neur", colnames(nmat))] = "Gaba"

    sex = rep("Female", ncol(nmat))
    sex[grep("Male", colnames(nmat))]="Male"               

    nannotation_col = data.frame(
        Clade = clade,
        Sex = sex
      )
    rownames(nannotation_col) = colnames(nmat) 

    p= pheatmap(nmat, annotation_col = nannotation_col, breaks = sequence(61,from = -30, by = 1),
             display_numbers = nsmat,main = paste(reg,"down in aging\n", database,"enrichment (homer)"),colorRampPalette(c("blue", "white", "red"))(61))
    print(p)
}
#dev.off()
                

options(repr.plot.width=8, repr.plot.height=8)
#png(paste(database, "_heatmaps.png", sep = ""),width = 7,height = 7)

for(reg in c("ENT", "FC", "AMY", "RLP", "NAC", "HCA", "HCP", ":CP")){

    nmat = mat[,grep(reg, colnames(mat))]
    nsmat = smat[,grep(reg, colnames(mat))]
    nsmat = nsmat[which(rowMaxs(nmat)>5),which(colMaxs(nmat)>3)]
    nmat = nmat[which(rowMaxs(nmat)>5),which(colMaxs(nmat)>3)]
    nmat[which(nmat>30, arr.ind = T)] = 30
    nmat[which(nmat<(-30), arr.ind = T)] = -30


    region = sapply(strsplit(colnames(nmat), ":"), function(x) tail(x, 1))
    clade = rep("NN",ncol(nmat))
    clade[grep("Glut", colnames(nmat))] = "Glut"
    clade[grep("Gaba", colnames(nmat))] = "Gaba"
    clade[grep("IMN", colnames(nmat))] = "IMN"
    clade[grep("Neur", colnames(nmat))] = "Gaba"

    sex = rep("Female", ncol(nmat))
    sex[grep("Male", colnames(nmat))]="Male"               

    nannotation_col = data.frame(
        Clade = clade,
        Sex = sex
      )
    rownames(nannotation_col) = colnames(nmat) 

    p= pheatmap(nmat, annotation_col = nannotation_col, breaks = sequence(61,from = -30, by = 1),
             display_numbers = nsmat,main = paste(reg,"up in aging\n", database,"enrichment (homer)"),colorRampPalette(c("blue", "white", "red"))(61))
    print(p)
}
#dev.off()
                

while(3>1){dev.off()}


options(repr.plot.width=10, repr.plot.height=10)

reg = "FC"

nmat = mat[,grep(reg, colnames(mat))]
nsmat = smat[,grep(reg, colnames(mat))]
nsmat = nsmat[which(rowMaxs(nmat)>5),which(colMaxs(nmat)>5)]
nmat = nmat[which(rowMaxs(nmat)>5),which(colMaxs(nmat)>5)]
nmat[which(nmat>50, arr.ind = T)] = 50
nmat[which(nmat<(-50), arr.ind = T)] = -50


region = sapply(strsplit(colnames(nmat), ":"), function(x) tail(x, 1))
clade = rep("NN",ncol(nmat))
clade[grep("Glut", colnames(nmat))] = "Glut"
clade[grep("Gaba", colnames(nmat))] = "Gaba"
clade[grep("IMN", colnames(nmat))] = "IMN"
clade[grep("Neur", colnames(nmat))] = "Gaba"

sex = rep("Female", ncol(nmat))
sex[grep("Male", colnames(nmat))]="Male"               

nannotation_col = data.frame(
    Clade = clade,
    Sex = sex
  )
rownames(nannotation_col) = colnames(nmat) 

pheatmap(nmat, annotation_col = nannotation_col, breaks = sequence(101,from = -50, by = 1),
         display_numbers = nsmat,main = paste(reg,"up in aging\n genomic enrichment (homer)"),colorRampPalette(c("blue", "white", "red"))(101))


sequence(-50,50,100)




