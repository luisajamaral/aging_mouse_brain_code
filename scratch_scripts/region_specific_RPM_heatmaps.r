library(ggplot2)
library(data.table)
library(UpSetR)
library(pheatmap)
library(Seurat)

library(zellkonverter)

setwd("../h5ads_final/")

h = readH5AD("celltype_major_batch_age_region_PMAT_RPM.h5ad")
#h = readH5AD("celltype_batch_age_region_rep_PMAT_RPM.h5ad")

h = readH5AD("celltype_age_PMAT_RPM.h5ad")

rm(h)

count_table <-h@assays@data$X


colnames(count_table)

setwd("../h5ads_final/combined_diff/diff_csvs//")
files = list.files(".", "diff_peaks", full.names = T)
#rfiles = list.files("../region_DARs_redone/", "diff_peaks", recursive = T, full.names=T)
files
#rfiles

files[1]

options(repr.plot.width=11, repr.plot.height=8)
afiles = files
atabs= list()
list_genes = list()

tab = c("celltype", "dir", "comparison","fraction")


for(f in afiles){
    cl = gsub("_2vs18_iter.csv","",f)
    cl = gsub("./diff_peaks_","",cl)

    # Read DAR file
    dar_data <- read.table(f, header=TRUE, sep=',', stringsAsFactors=FALSE)
    dar_data$cl= cl
  
    # Filter rows with significant changes
    significant_changes <- subset(dar_data, adjusted.p.value < 0.01 & abs(log2.fold_change.)>0.25)
    if(nrow(significant_changes)<1){
        next
    }
    #if(nrow(significant_changes)>8000){
    #    significant_changes = significant_changes[1:8000,]
    #}
    ct = gsub("_", " ", cl)
    if(length(grep("^./diff" , f))>0) {
        options(repr.plot.width=4, repr.plot.height=6)

      #  oligo = count_table[paste(significant_changes$feature.name),grep(ct, colnames(count_table))]
      #  pheatmap(oligo[, which(colSums(oligo)>0)], scale = "row", show_rownames = F)
        
        
        
        oligo = count_table[paste(significant_changes[which(significant_changes$log2.fold_change.> 0.25),"feature.name"]),grep(ct, colnames(count_table))]
        if(length(significant_changes[which(significant_changes$log2.fold_change.> 0.25),"feature.name"])>10){
            diff29 = oligo[,2]-oligo[,1]
            diff218 = oligo[,3]-oligo[,1]
            diff918 = oligo[,3]-oligo[,2]
            
            cat(cl , "up 2v9 concordant",table(sign(diff29))['-1']/table(sign(diff218))[1], "\n")
            cat(cl , "up 9v18 concordant",table(sign(diff918))['-1']/table(sign(diff218))[1], "\n")
            tab = rbind(tab, c(cl, "up", "2v9", table(sign(diff29))['-1']/table(sign(diff218))[1]))
            tab = rbind(tab, c(cl, "up", "9v18", table(sign(diff918))['-1']/table(sign(diff218))[1]))

        }
        oligo = count_table[paste(significant_changes[which(significant_changes$log2.fold_change.<  -0.25),"feature.name"]),grep(ct, colnames(count_table))]
        if(length(significant_changes[which(significant_changes$log2.fold_change.<  -0.25),"feature.name"])>10){
            diff29 = oligo[,2]-oligo[,1]
            diff218 = oligo[,3]-oligo[,1]
            diff918 = oligo[,3]-oligo[,2]

            cat(cl ,"down 2v9 concordant", table(sign(diff29))['1']/table(sign(diff218))[1], "\n")
            cat(cl , "down 9v18 concordant",table(sign(diff918))['1']/table(sign(diff218))[1], "\n")
            tab = rbind(tab, c(cl, "down", "2v9", table(sign(diff29))['1']/table(sign(diff218))[1]))
            tab = rbind(tab, c(cl, "down", "9v18", table(sign(diff918))['1']/table(sign(diff218))[1]))

        }
        }
    
}


tab = as.data.frame(tab)
colnames(tab) =  c("celltype", "dir", "comparison","fraction")
tab = tab[-1,]
head(tab)


library(dplyr)
result <- tab %>%
  filter(dir == "up" & comparison == "2v9") %>%
  summarise(mean_fraction = mean(fraction))

# Print the result
print(result[1])

tab %>%
  filter(dir == "up" & comparison == "2v9") %>%
  summarise(mean_fraction = mean(fraction))

tab %>%
  filter(dir == "up" & comparison == "9v18") %>%
  summarise(mean_fraction = mean(fraction))

tab %>%
  filter(dir == "down" & comparison == "2v9") %>%
  summarise(mean_fraction = mean(fraction))

tab %>%
  filter(dir == "down" & comparison == "9v18") %>%
  summarise(mean_fraction = mean(fraction))

options(repr.plot.width=11, repr.plot.height=10)


tab$fraction = as.numeric(as.character(tab$fraction))
head(tab)

ggplot(tab, aes(x=celltype, y=fraction, fill=comparison)) + ggtitle("fraction of concordant changes") +

  geom_bar(stat="identity", position = 'dodge')+theme_minimal()+coord_flip()+facet_wrap(~dir)

options(repr.plot.width=5, repr.plot.height=5)

p <- ggplot(tab, aes(x=comparison, y=fraction,fill=dir)) + ylim(c(0,1))+
geom_boxplot()+ theme_classic()
p

options(repr.plot.width=5, repr.plot.height=5)

p <- ggplot(tab, aes(x=comparison, y=fraction,fill=dir)) + ylim(c(0,1))+
geom_boxplot()+ theme_classic()
p

options(repr.plot.width=11, repr.plot.height=8)

tab = as.data.frame(tab)
colnames(tab) =  c("celltype", "dir", "comparison","fraction")
tab = tab[-1,]
tab$fraction = as.numeric(as.character(tab$fraction))
head(tab)

ggplot(tab, aes(x=celltype, y=fraction, fill=comparison)) + ggtitle("fraction of concordant changes") 
  geom_bar(stat="identity", position = 'dodge')+theme_minimal()+coord_flip()+facet_wrap(~dir)

oligo


significant_changes

options(repr.plot.width=11, repr.plot.height=8)
afiles = c(rfiles,files)

cl="Oligo_NN"
afiles = afiles[grep(cl, afiles)]
atabs= list()
list_genes = list()

for(f in afiles){
    # Read DAR file
    dar_data <- read.table(f, header=TRUE, sep=',', stringsAsFactors=FALSE)
    dar_data$cl= gsub("_2vs18.csv","",f)
    dar_data$cl= gsub("../region_DARs_redone//","",dar_data$cl)
    dar_data$cl= gsub("diff_peaks_","",dar_data$cl)

    # Filter rows with significant changes
    significant_changes <- subset(dar_data, adjusted.p.value < 0.01 & abs(log2.fold_change.)>0.25)
    if(nrow(significant_changes)>8000){
        significant_changes = significant_changes[1:8000,]
    }
    atabs[[f]] = significant_changes
    #list_genes[[dar_data$cl[1]]] = rownames(significant_changes)
    list_genes[[dar_data$cl[1]]] = significant_changes$feature.name
    ct = gsub("_", " ", cl)
    if(length(grep("^./diff" , f))>0) {
        options(repr.plot.width=4, repr.plot.height=6)

        oligo = count_table[paste(significant_changes$feature.name),grep(ct, colnames(count_table))]
        pheatmap(oligo[, which(colSums(oligo)>0)], scale = "row", show_rownames = F)
        
        
        
        oligo = count_table[paste(significant_changes[which(significant_changes$log2.fold_change.> 0.25),"feature.name"]),grep(ct, colnames(count_table))]

        diff29 = oligo[,2]-oligo[,1]
        diff218 = oligo[,3]-oligo[,1]
        diff918 = oligo[,3]-oligo[,2]

        cat("up 2v9 concordant",table(sign(diff29))[2]/table(sign(diff218))[1], "\n")
        cat("up 9v18 concordant",table(sign(diff918))[2]/table(sign(diff218))[1], "\n")
        oligo = count_table[paste(significant_changes[which(significant_changes$log2.fold_change.<  -0.25),"feature.name"]),grep(ct, colnames(count_table))]

        diff29 = oligo[,2]-oligo[,1]
        diff218 = oligo[,3]-oligo[,1]
        diff918 = oligo[,3]-oligo[,2]

        cat("down 2v9 concordant", table(sign(diff29))[1]/table(sign(diff218))[1], "\n")
        cat("down 9v18 concordant",table(sign(diff918))[1]/table(sign(diff218))[1], "\n")
        }
    
}
tabs = do.call(rbind, atabs)
locs = unique(tabs$feature.name)
length(locs)
upset(fromList(list_genes), order.by = "freq", nsets = 30)



table(sign(diff29))[2]/table(sign(diff218))[1]

tabs = tabs[-which(duplicated(tabs$feature.name)),]
tabs = tabs[order(tabs$log2.fold_change.),]
head(tabs)

oligo = count_table[paste(significant_changes[which(significant_changes$log2.fold_change.> 0.25),"feature.name"]),grep(ct, colnames(count_table))]

diff29 = oligo[,2]-oligo[,1]
diff218 = oligo[,3]-oligo[,1]
diff918 = oligo[,3]-oligo[,2]

table(sign(diff29))[2]/table(sign(diff218))[1]
table(sign(diff918))[2]/table(sign(diff218))[1]

oligo = count_table[paste(significant_changes[which(significant_changes$log2.fold_change.<  -0.25),"feature.name"]),grep(ct, colnames(count_table))]

diff29 = oligo[,2]-oligo[,1]
diff218 = oligo[,3]-oligo[,1]
diff918 = oligo[,3]-oligo[,2]

table(sign(diff29))[1]/table(sign(diff218))[1]
table(sign(diff918))[1]/table(sign(diff218))[1]

table(sign(diff29))[2]/table(sign(diff218))[1]
table(sign(diff918))[2]/table(sign(diff218))[1]

age = sapply(strsplit(as.character(colnames(oligo)), ":"), "[[", 2)
annotation_col = data.frame(
    age = age
    
  )
rownames(annotation_col) = colnames(oligo)

age

pheatmap(oligo[, which(colSums(oligo)>4)], scale = "row", show_rownames = F)



setwd("../region_DARs/")
#files = list.files(".", ".csv", full.names = T)
cfiles = list.files("../combined_diff/diff_csvs", ".csv", full.names = T)
ffiles = list.files("../female_diff/diff_csvs", ".csv", full.names = T)
mfiles = list.files("../male_diff/diff_csvs", ".csv", full.names = T)




cl = "D12"
#afiles = c(files,cfiles,ffiles,mfiles)
afiles = c(cfiles,ffiles,mfiles)

afiles = afiles[grep(cl, afiles)]
afiles 


options(repr.plot.width=11, repr.plot.height=8)

atabs= list()
list_genes = list()

for(f in afiles){
    # Read DAR file
    dar_data <- read.table(f, header=TRUE, sep=',', stringsAsFactors=FALSE)
    dar_data$cl= gsub("_iter.csv","",f)
    dar_data$cl= gsub("../combined_diff/diff_csvs/diff_peaks_","Combined_",dar_data$cl)
    dar_data$cl= gsub("../female_diff/diff_csvs/diff_peaks_","Female_",dar_data$cl)
    dar_data$cl= gsub("../male_diff/diff_csvs/diff_peaks_","Male_",dar_data$cl)

    dar_data$cl= gsub("./diff_peaks_","",dar_data$cl)

    # Filter rows with significant changes
    significant_changes <- subset(dar_data, adjusted.p.value < 0.01 & abs(log2.fold_change.)>0.25)
    if(nrow(significant_changes)>8000){
        significant_changes = significant_changes[1:8000,]
    }
    atabs[[f]] = significant_changes
    #list_genes[[dar_data$cl[1]]] = rownames(significant_changes)
    list_genes[[dar_data$cl[1]]] = significant_changes$feature.name
    
}
tabs = do.call(rbind, atabs)
locs = unique(tabs$feature.name)
length(locs)
upset(fromList(list_genes), order.by = "freq", nsets = 30)



tabs = tabs[-which(duplicated(tabs$feature.name)),]
tabs$log2.fold_change. = -tabs$log2.fold_change.

tabs = tabs[order(tabs$log2.fold_change.),]
head(tabs)


locs = tabs$feature.name
length(locs)

ct = gsub("_", " ", cl)
#oligo = count_table[paste(locs),grep(ct, colnames(count_table))]
oligo = count_table[paste(locs),]

head(oligo)

head(oligo)

region = sapply(strsplit(as.character(colnames(oligo)), "_"), "[[", 4)
age = sapply(strsplit(as.character(colnames(oligo)), "_"), "[[", 3)
sex = sapply(strsplit(as.character(colnames(oligo)), "_"), "[[", 2)
#age = sapply(strsplit(as.character(age), "_"), "[[", 1)
annotation_col = data.frame(
    age = age,
    region = region,
    sex = sex
    
  )
annotation_col$age =factor(annotation_col$age, levels = c("2mo","9mo", "18mo"))
rownames(annotation_col) = colnames(oligo)

head(oligo)

outlist=list()
for(i in c("Oligo", "Gaba","Glut","Micro", "Astro_", "OPC")){
    cur = oligo[, grep(i, colnames(oligo))]
    fcur = cur[, grep("Female", colnames(cur))]
    fcur = t(scale(t(fcur)))
    mcur = cur[, grep("Male", colnames(cur))]
    mcur = t(scale(t(mcur)))
    outlist[[paste(i,"Female")]]=fcur
    outlist[[paste(i,"Male")]]=mcur
    }
o = do.call(cbind, outlist)
o <- apply(o, 2, function(x) replace(x, is.na(x), min(x, na.rm = TRUE)))
head(o)

#oligo = oligo[,c(grep("Oligo", colnames(oligo)),grep("Gaba", colnames(oligo)), grep("Glut", colnames(oligo)), grep("Micro", colnames(oligo)), grep("Astro_", colnames(oligo)), grep("OPC", colnames(oligo)))]

#fem = oligo[, grep("Female", colnames(oligo))]
#fem = t(scale(t(fem)))
#mal = oligo[, grep("Male", colnames(oligo))]
#mal = t(scale(t(mal)))
#o = cbind(fem,mal)
#o <- apply(o, 2, function(x) replace(x, is.na(x), min(x, na.rm = TRUE)))


region = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 4)
age = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 3)
sex = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 2)
sex = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 2)
celltype = sapply(strsplit(as.character(colnames(o)), "_"), "[[", 1)

#age = sapply(strsplit(as.character(age), "_"), "[[", 1)
annotation_col = data.frame(
    age = age,
    region = region,
    sex = sex,
    celltype = celltype
    
  )
annotation_col$age =factor(annotation_col$age, levels = c("2mo","9mo", "18mo"))
rownames(annotation_col) = colnames(o)

options(repr.plot.width=6, repr.plot.height=6)

calculate_wss <- function(data, k_max) {
  wss <- numeric(k_max)
  for (i in 1:k_max) {
    kmeans_result <- kmeans(data, centers=i)
    wss[i] <- kmeans_result$tot.withinss
  }
  return(wss)
}
# Define the maximum number of clusters to consider
k_max <- 10
# Compute WSS for each value of k
wss_values <- calculate_wss(o, k_max)
# Plot the elbow plot
plot(1:k_max, wss_values, type="b", pch=19, frame=FALSE,
     xlab="Number of clusters (k)", ylab="Total within-cluster sum of squares (WSS)",
     main="Elbow Method to Find Optimal K")
# Add text to indicate the optimal k
elbow_k <- 3
#text(elbow_k, wss_values[elbow_k], paste("Optimal k =", elbow_k), pos=3, col="red", cex=1.2)






brks <- seq(-3,3,length.out=40) 


options(repr.plot.width=36, repr.plot.height=6)

ph = pheatmap(o[,order(  annotation_col$celltype,annotation_col$sex,annotation_col$age ,annotation_col$region)],kmeans_k = 12, main = paste(nrow(o), "DARs"), show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks)#,


options(repr.plot.width=36, repr.plot.height=6)

ph = pheatmap(o[,order(  annotation_col$celltype,annotation_col$sex,annotation_col$age ,annotation_col$region)], main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   cluster_rows = T,cluster_cols = T,breaks = brks)#,


options(repr.plot.width=16, repr.plot.height=6)

c = o[,grep("Gaba", colnames(o))]
canno = annotation_col[grep("Gaba", rownames(annotation_col)),]
p2 = pheatmap(c[,order( canno$sex,canno$region,canno$age)], main = paste(nrow(o), "DARs"), cutree_rows = 2,show_rownames = F,gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3),scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   cluster_rows = T,cluster_cols = F,breaks = brks)#,


options(repr.plot.width=16, repr.plot.height=6)

c = o[,grep("Micro", colnames(o))]
canno = annotation_col[grep("Micro", rownames(annotation_col)),]
p2 = pheatmap(c[,order( canno$sex,canno$region,canno$age)], main = paste(nrow(o), "DARs"), cutree_rows = 2,show_rownames = F,gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3),scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   cluster_rows = T,cluster_cols = F,breaks = brks)#,


c = o[,grep("Micro", colnames(o))]
canno = annotation_col[grep("Micro", rownames(annotation_col)),]
p2 = pheatmap(c[,order( canno$sex,canno$region,canno$age)], main = paste(nrow(o), "DARs"), cutree_rows = 2,show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   cluster_rows = T,cluster_cols = T,breaks = brks)#,


options(repr.plot.width=8, repr.plot.height=6)
ct = "Micro"
c = o[,grep(ct, colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = ct, show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = T,breaks = brks)#,)#,


options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep(ct, colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = T,breaks = brks)#,)#,



options(repr.plot.width=8, repr.plot.height=6)
ct = "Micro"
c = o[,grep(ct, colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = ct, show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,


options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep(ct, colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



options(repr.plot.width=36, repr.plot.height=6)

ph = pheatmap(o[p2$tree_row$order,order(  annotation_col$celltype,annotation_col$sex,annotation_col$age ,annotation_col$region)], main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   cluster_rows = T,cluster_cols = T,breaks = brks)#,


options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep("Oligo", colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep("Oligo", rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = "Oligo", show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,


options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep("Oligo", colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep("Oligo", rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep("Oligo", colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep("Oligo", rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = "Oligo", show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,


options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep("Oligo", colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep("Oligo", rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



nrow(c)

options(repr.plot.width=8, repr.plot.height=6)
ct = "Micro"
c = o[,grep(ct, colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = paste(ct), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,


options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep(ct, colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)],  show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



options(repr.plot.width=8, repr.plot.height=6)
ct = "Gaba"
c = o[,grep(ct, colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = ct, show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,


options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep(ct, colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)],  show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



options(repr.plot.width=8, repr.plot.height=6)
ct = "Glut"
c = o[,grep(ct, colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = ct, show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,


options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep(ct, colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep(ct, rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep("Oligo", colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep("Oligo", rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep("Astro", colnames(o))]
c= c[,grep("Male", colnames(c))]
canno = annotation_col[grep("Astro", rownames(annotation_col)),]
canno = canno[grep("Male", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep("Gaba", colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep("Gaba", rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[p2$tree_row$order,order( canno$sex,canno$region,canno$age)], main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,



options(repr.plot.width=8, repr.plot.height=6)

c = o[,grep("Glut", colnames(o))]
c= c[,grep("Female", colnames(c))]
canno = annotation_col[grep("Glut", rownames(annotation_col)),]
canno = canno[grep("Female", rownames(canno)),]

pheatmap(c[ph$tree_row$order,order( canno$sex,canno$region,canno$age)], main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   
         cluster_rows = F,cluster_cols = F,breaks = brks, gaps_col = seq(from = 3, by = 3, length.out = ncol(c)/3))#,)#,




c = o[,grep("Gaba", colnames(o))]
pheatmap(c[ph$tree_row$order,], main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),   cluster_rows = F,cluster_cols = T,breaks = brks)#,


options(repr.plot.width=36, repr.plot.height=6)

ph = pheatmap(o[,order(  annotation_col$celltype,annotation_col$sex,annotation_col$age ,annotation_col$region)],kmeans_k = 12, main = paste(nrow(o), "DARs"), show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks)#,


options(repr.plot.width=16, repr.plot.height=6)

ol = o[,grep("Oligo",colnames(o))]
ph = pheatmap(ol, main = paste(nrow(o), "DARs"), show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_rows = F,cluster_cols = F,breaks = brks)#,


ph = pheatmap(ol, main = paste(nrow(o), "DARs"), show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_rows = T,cluster_cols = F,breaks = brks)#,


options(repr.plot.width=16, repr.plot.height=6)
ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],kmeans_k = 12, main = paste(nrow(o), "DARs"), show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,

ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],  cutree_rows = 2, main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,

ph = pheatmap(o[,],  cutree_rows = 2, main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#T

ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],  cutree_rows = 2, main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,

ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],  cutree_rows = 8, main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,

options(repr.plot.width=16, repr.plot.height=6)
ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],kmeans_k = 12, main = paste(nrow(o), "DARs"), show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,

ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],  cutree_rows = 8, main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,


#o = cbind(cp, amy,nac)
options(repr.plot.width=26, repr.plot.height=6)

ph = pheatmap(o[1:2000,], main = paste("top 2000 DARs"), kmeans_k = 12,show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col =c(24))#,
ph = pheatmap(o, main = paste(" DARs"), kmeans_k = 12,show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col =c(24))#,


ph = pheatmap(o[,], main = paste("top DARs"), kmeans_k = 12,show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = T,breaks = brks,gaps_col =c(24))#,


options(repr.plot.width=13, repr.plot.height=5)

ph = pheatmap(o, main = paste(nrow(o), "DARs"), kmeans_k = 5,show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col =c(24))#,


seq(from = 3, by = 3, length.out = ncol(o)/3)

ncol(o)/3

options(repr.plot.width=13, repr.plot.height=6)

pheatmap(o, main = paste(nrow(o), "DARs"),cutree_rows = 5,annotation_col=annotation_col,show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col =c(24))#,


pheatmap(o[,grep("Female", colnames(o))], main = paste(nrow(o), "DARs"),cutree_rows = 5,annotation_col=annotation_col,show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = T,breaks = brks)#,


pheatmap(o[,grep("Male", colnames(o))], main = paste(nrow(o), "DARs"),cutree_rows = 5,annotation_col=annotation_col,show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = T,breaks = brks)#,


ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],kmeans_k = 12, main = paste(nrow(o), "DARs"), show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,


ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)], main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,


options(repr.plot.width=13, repr.plot.height=6)
ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],kmeans_k = 12, main = paste(nrow(o), "DARs"), show_rownames = T,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,

ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)],  cutree_rows = 8, main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,


ph = pheatmap(o[,order( annotation_col$region ,annotation_col$sex)], cutree_rows = 8,main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,


ph = pheatmap(o[1:5000,order( annotation_col$region ,annotation_col$sex)],cutree_rows = 6, main = paste(nrow(o), "DARs"), show_rownames = F,scale = "none",annotation_col=annotation_col,color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col = seq(from = 3, by = 3, length.out = ncol(o)/3))#,


options(repr.plot.width=13, repr.plot.height=6)

pheatmap(o[1:3000,], main = paste(nrow(o), "DARs"),cutree_rows = 6,annotation_col=annotation_col,show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col =c(24))#,


pheatmap(o[1:2000,], main = paste(nrow(o), "DARs"),cutree_rows = 4,annotation_col=annotation_col,show_rownames = F,scale = "none",color=colorRampPalette(c("navy", "white", "red"))(40),  cluster_cols = F,breaks = brks,gaps_col =c(3))#,



