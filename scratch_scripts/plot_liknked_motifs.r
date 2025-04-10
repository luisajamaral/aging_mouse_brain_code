library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)

head(curr$`# of Target Sequences with Motif(of 1)`)

setwd("results_logfc.15_adjpval.1/")

database = "motifs"
files = list.files(".", pattern = "knownResults.txt", full.names = TRUE, recursive = TRUE)
dir = "motifs"
#files = files[grep("Oligo_NN", files)]
files = files[grep(dir, files)]
files = files[grep("bg", files)]
files = files[-grep("logfc.15", files)]

pval = 1e-4

files = files[which(file.exists(files))]
f = files[1]
curr = fread(f, sep = "\t",fill = T)
curr = as.data.frame(curr)
curr = curr[which(curr$`P-value` < pval),]

gkey = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
for(f in files[2:length(files)]) {
  if(!file.exists(f)){
    cat(f, " NO FILE\n")
    next()
  }
  curr = fread(f, sep = "\t",fill = T)
  curr = as.data.frame(curr)
  #  show(head(curr))
  curr = curr[which(curr$`P-value` < pval),]
  #cat(f, nrow(curr),"\n")

  gkey = c(gkey , paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1)))
 # gkey = c(gkey ,paste(curr[,1]))
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
    curr = fread(f, sep = "\t", fill=T)
    curr = as.data.frame(curr)
    names = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
    which(duplicated(names))
    curr = curr[which(!duplicated(names)),]
    rownames(curr) = paste(sapply(strsplit(as.character(curr[,c(1)]), "/"), "[[", 1))
    if(paste(gkey[i]) %in% paste(rownames(curr))) {
      # mat[x,y]  = -log10(curr[paste(gkey[i]), "P.Value"]+1e-100)
      mat[x,y]  = -curr[paste(gkey[i]),"Log P-value"]
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
colnames(mat) = gsub("/knownResults.txt", "",files)
colnames(mat) = gsub("./", "", colnames(mat))
colnames(mat) = gsub("_locs_motifs", "", colnames(mat))
colnames(mat) = gsub("_bg", "", colnames(mat))

snmat = smat
nmat = mat
#nmat[which(nmat>100, arr.ind = T)] = 100

dir = rep("down" , ncol(nmat))
dir[grep("up", as.character(colnames(nmat)))] = "up"


clade= sapply(sapply(strsplit(colnames(nmat), ":"), `[`, 1), function(x) {
  parts <- strsplit(x, "_")[[1]]  # Split the string by space
  last_element <- parts[length(parts)-1]  # Get the last element
  return(last_element)
})

nannotation_col = data.frame(
     Clade = clade,  Dir = dir
  )
rownames(nannotation_col) = colnames(nmat) 
               

snmat = snmat[,order(dir,clade,colnames(nmat))]
nmat = nmat[,order(dir,clade,colnames(nmat) )]




options(repr.plot.width=10, repr.plot.height=10)
snmat = snmat[which(rowMaxs(nmat)>3),which(colMaxs(nmat)>4)]
nmat = nmat[which(rowMaxs(nmat)>3),which(colMaxs(nmat)>4)]
nmat[which(nmat>30, arr.ind = T)] = 30

heatmap <- pheatmap(nmat, 
                    annotation_col = nannotation_col, 
                    cluster_cols = T, 
                    main = "DEG linked DAR Motifs",  clustering_method = 'ward.D2',
                    display_numbers = snmat, 
                    cluster_rows = T)

options(repr.plot.width=7, repr.plot.height=10
       )

nmat = mat[which(rowMaxs(mat)>4), ]
nmat = nmat[,which(colMaxs(nmat)>3)]
nmat[which(nmat>30, arr.ind = T)] = 30

heatmap <- pheatmap(nmat, main = dir)
heatmap

options(repr.plot.width=7, repr.plot.height=10
       )

nmat = mat[which(rowMaxs(mat)>4), ]
nmat = nmat[,which(colMaxs(nmat)>3)]
nmat[which(nmat>30, arr.ind = T)] = 30

heatmap <- pheatmap(nmat, main = dir)
heatmap

cols = c('Motif Name','Consensus','P-value','Log P-value','q-value (Benjamini)','% of Target Sequences with Motif','% of Background Sequences with Motif','file')
pval = 1.1
f = files[2]
curr = fread(f, sep = "\t", fill = T)
curr = as.data.frame(curr)
#curr = curr[which(curr$`P-value`< pval),]
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
 # curr = curr[which(curr$`P-value`< pval),]
  if(nrow(curr)>0){
  curr$file = f 
  curr= curr[,cols]
  big = rbind(big,curr)
      }
}

head(big,5)

big = big[order(big$`Log P-value`),]
big = big[-which(duplicated(paste(big$file, big$Name))),]
head(big)



big$file = gsub("/knownResults.txt", "",big$file)
big$file = gsub("./", "", big$file)
big$file = gsub("_locs_motifs", "", big$file)
big$file = gsub("_bg", "",big$file)




big$dir = rep("down" , nrow(big))
big$dir[grep("up", as.character(big$file))] = "up"



big$clade= sapply(sapply(strsplit(big$file, ":"), `[`, 1), function(x) {
  parts <- strsplit(x, "_")[[1]]  # Split the string by space
  last_element <- parts[length(parts)-1]  # Get the last element
  return(last_element)
})
big$clade[which(big$clade == "Neur")] = "Gaba"
big$clade[which(big$clade == "Gaba-Chol")] = "Gaba"
big$clade[which(big$clade == "Glut-Sero")] = "Glut"

big$ct = gsub("_down", "",sapply(strsplit(big$file, ":"), `[`, 1))
big$ct = gsub("_down", "",big$file)
big$ct = gsub("up_", "",big$ct)
big$ct = factor(big$ct, levels = unique(big$ct[order(big$clade, big$dir)]))
big$Name = paste(sapply(strsplit(as.character(big$`Motif Name`), "/"), "[[", 1))

unique(big$`Motif Name`)

library(dplyr)
ele = 'AP-1(bZIP)'
#df = big[grep(ele, big$Name),]
df = big[which(big$Name ==ele),]


# Adjusting y-values to not exceed the ylim maximum and adding significance stars
df <- df %>%
  mutate(
    `Log P-value capped` = pmin(-`Log P-value`, 30),
    stars = case_when(
      `q-value (Benjamini)` < 1e-100 ~ "***",
      `q-value (Benjamini)` < 1e-10 ~ "**",
      `q-value (Benjamini)` < 0.01 ~ "*",
      TRUE ~ ""
    )
  )

# Creating the plot
ggplot(data=df, aes(x=`ct`, y=`Log P-value capped`, fill=`clade`)) + 
  ggtitle(paste(ele, "Motif")) + 
  ylim(c(0,32)) +
  geom_bar(stat="identity", width=0.5, position="dodge") + 
  geom_text(aes(label=stars), 
            position=position_dodge(width=0.5), 
            vjust=1.2, hjust=0.5,
            angle=90, 
            size=5) + 
  facet_wrap(~dir) +
  coord_flip() + 
  theme_classic() + 
  ylab("-Log P-value")


options(repr.plot.width=9, repr.plot.height=10
       )

nmat = mat[which(rowMaxs(mat)>4), ]
nmat = nmat[,which(colMaxs(nmat)>3)]
nmat[which(nmat>20, arr.ind = T)] = 20

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

    rownames(nannotation_col) = colnames(nmat) 

nmat = nmat[,order(dir, clade)]
display_matrix <- ifelse(nmat > 4.6, '*', '')
# Plot the heatmap with display_numbers parameter
library(pheatmap)
heatmap <- pheatmap(nmat, cluster_cols = F,gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    display_numbers = display_matrix, main = "DEG-Linked peaks Motifs" , annotation_col = nannotation_col)
heatmap

options(repr.plot.width=8, repr.plot.height=12
       )

nmat = mat[which(rowMaxs(mat)>4), ]
nmat = nmat[,which(colMaxs(nmat)>3)]
nmat[which(nmat>20, arr.ind = T)] = 20

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

    rownames(nannotation_col) = colnames(nmat) 

nmat = nmat[,order(dir, clade)]
display_matrix <- ifelse(nmat > 4.6, '*', '')
# Plot the heatmap with display_numbers parameter
library(pheatmap)
heatmap <- pheatmap(nmat, cluster_cols = F,gaps_col = c(grep("up", colnames(nmat))[1]-1),
                    display_numbers = display_matrix, main = database, annotation_col = nannotation_col)
heatmap

nmat






