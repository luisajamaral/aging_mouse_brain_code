library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

database = "repeats.genomeOntology.txt"
setwd(paste("../combined_DARs_redone/",sep = ""))


files = list.files(".", pattern = database, full.names = TRUE, recursive = TRUE)
files = paste(files, sep = "")
files = files[which(file.exists(files))]
files = files[grep("bg", files)]

length(files)

pval = 1e-10

f = files[1]
curr = fread(f, sep = "\t", fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$`P-value`< pval),]
curr= curr[which(curr$`Overlap(#peaks)`>5),]
head(curr)
gkey = curr$Term
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  curr = curr[which(curr$`P-value`< pval),]
  curr= curr[which(curr$`Overlap(#peaks)`>5),]
  curr = curr[1:min(20, nrow(curr)),]
  cat(f , length(paste(curr[,1])), "\n")
  gkey = c(gkey ,paste(curr[,1]))
  # cat((gkey))
}

gkey = gkey[which(!duplicated(gkey))]
gkey = gkey[which(gkey!="NA")]
cat(length(gkey))

gkey

head(curr)


# getting significance for each term + file
smat = matrix('', nrow=length(gkey),ncol=length(files))
mat = matrix(0, nrow=length(gkey),ncol=length(files))
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
      #  cat(mat[x,y], " ")
    if(curr[paste(gkey[i]),"P-value(vs. Control)"]<0.001 & curr[paste(gkey[i]),"Log P-value(+ underrepresented, vs. Control)"]<0){
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
colnames(mat) = gsub("_genOn_bg", "",files)
colnames(mat) = gsub("/", "",colnames(mat))

colnames(mat) = gsub(database, "",colnames(mat))

colnames(mat) = gsub("./", "", colnames(mat))
colnames(mat) = gsub("peaks_", "", colnames(mat))
colnames(mat) = gsub("[.]", "", colnames(mat))


snmat = smat[which(rowMaxs(mat)>4),which(colMaxs(mat)>4)]
nmat = mat[which(rowMaxs(mat)>4 ),which(colMaxs(mat)>4)]
nmat[which(nmat>100, arr.ind = T)] = 100
nmat[which(nmat< 0, arr.ind = T)] = 0

dir = rep("down" , ncol(nmat))
dir[grep("up", as.character(colnames(nmat)))] = "up"


clade= sapply(sapply(strsplit(colnames(nmat), ":"), `[`, 1), function(x) {
  parts <- strsplit(x, "_")[[1]]  # Split the string by space
  last_element <- parts[length(parts)]  # Get the last element
  return(last_element)
})
clade[which(clade == "Neur")] = "Gaba"
clade[which(clade == "Gaba-Chol")] = "Gaba"
clade[which(clade == "Glut-Sero")] = "Glut"

nannotation_col = data.frame(
    Dir = dir,  Clade = clade
  )
rownames(nannotation_col) = colnames(nmat) 
               
nmat = nmat[,order(dir,  clade)]
snmat = snmat[,order(dir,  clade)]

options(repr.plot.width=19, repr.plot.height=17)


heatmap <- pheatmap(nmat, 
                    annotation_col = nannotation_col, 
                    cluster_cols = F, 
                    clustering_method = 'ward.D2',
                    main = database, display_numbers=snmat,
                    
                    cluster_rows = T)
#pdf("Reactome_DAR_GO.pdf", width =14, height = 12)
#print(heatmap)
#dev.off()

files

pval = 1e-2
cols = c('Name','PeakFile/Annotation','Overlap(#peaks)','Overlap(bp)','Expected Overlap(bp, gsize=2.00e+09)','Log Ratio Enrichment','Log P-value(+ underrepresented)','P-value','Control Log P-value(+ underrepresented)','Control P-value','Log Ratio Enrichment(vs. Control)','Log P-value(+ underrepresented, vs. Control)','P-value(vs. Control)','file')
f = files[1]
curr = fread(f, sep = "\t", fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$`P-value`< pval),]
curr= curr[which(curr$`Overlap(#peaks)`>5),]
big = curr
big$file = f
big = big[,cols]

gkey = curr$Term
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  curr = curr[which(curr$`P-value`< pval),]
  curr= curr[which(curr$`Overlap(#peaks)`>5),]
  if(nrow(curr)>0){
  curr$file = f 
  curr= curr[,cols]
  big = rbind(big,curr)
      }
}

head(big,20)


big$file = gsub("_genOn_bg", "",big$file)
big$file = gsub("/", "",big$file)

big$file = gsub(database, "",big$file)

big$file = gsub("./", "", big$file)
big$file = gsub("peaks_", "", big$file)
big$file = gsub("[.]", "", big$file)


big$dir = rep("down" , nrow(big))
big$dir[grep("up", as.character(big$file))] = "up"


big$clade= sapply(big$file, function(x) {
  parts <- strsplit(x, "_")[[1]]  # Split the string by space
  last_element <- parts[length(parts)]  # Get the last element
  return(last_element)
})
big$clade[which(big$clade == "Neur")] = "Gaba"
big$clade[which(big$clade == "Gaba-Chol")] = "Gaba"
big$clade[which(big$clade == "Glut-Sero")] = "Glut"

big$ct = gsub("down_", "",big$file)
big$ct = gsub("up_", "",big$ct)
big$ct = factor(big$ct, levels = unique(big$ct[order(big$clade)]))


options(repr.plot.width=7, repr.plot.height=7)

#big$file = factor(big$file, levels = unique(big$file[order(big$clade)]))

df = big[grep("H3K9me3_re", big$Name),]
ggplot(data=df, aes(x=`ct`, y=-`Log P-value(+ underrepresented, vs. Control)`, fill = `clade`)) + ggtitle("H3K9me3_broadpeaks_P0_forebrain")+
  geom_bar(stat="identity", width=0.5) + facet_wrap(~dir)+coord_flip() +theme_classic() + ylab("-Log P-value(+ overrepresented, vs. Control)")

unique(big$Name)

options(repr.plot.width=7, repr.plot.height=7)

ele = 'LTR'
df = big[grep(ele, big$Name),]
df = big[which(big$Name ==ele),]

ggplot(data=df, aes(x=`ct`, y=-`Log P-value(+ underrepresented, vs. Control)`, fill = `clade`)) + ggtitle(ele)+
  geom_bar(stat="identity", width=0.5) + facet_wrap(~dir)+coord_flip() +theme_classic() + ylab("-Log P-value(+ overrepresented, vs. Control)")

options(repr.plot.width=10, repr.plot.height=10)

df = cbind(dir, clade,H3K9me3_overlap=mat["CTCF_P0_forebrain.bed",])
df = as.data.frame(df)

df$H3K9me3_overlap=as.numeric(as.character(df$H3K9me3_overlap))
ggplot(data=df, aes(x=rownames(df), y=H3K9me3_overlap, fill = clade)) +
  geom_bar(stat="identity", width=0.5) + facet_wrap(~dir, scales="free_y")+coord_flip()

options(repr.plot.width=12, repr.plot.height=12)


heatmap <- pheatmap(nmat, 
                    annotation_col = nannotation_col, 
                    cluster_cols = F, 
                    clustering_method = 'ward.D2',
                    main = database, 
                    display_numbers = snmat, 
                    cluster_rows = T)
#pdf("Reactome_DAR_GO.pdf", width =14, height = 12)
#print(heatmap)
#dev.off()

max(nmat)



