library(ggplot2)

for(ct in c("STR_D12_Gaba")){
    f1 = read.table(paste("male_diff/diff_csvs/diff_peaks_",ct,":Male_2vs18_iter.csv",sep = ""), sep = ",",header =T)
    f2 = read.table(paste("female_diff/diff_csvs/diff_peaks_",ct,":Female_2vs18_iter.csv",sep = ""), sep = ",",header =T)
    sig = c(f1$feature.name[which(f1$adjusted.p.value<0.01)], f2$feature.name[which(f2$adjusted.p.value<0.01)])
    length(sig)
    sig = unique(sig)
    length(sig)
    rownames(f1)= f1$feature.name
    rownames(f2)= f2$feature.name
    f1s = f1[sig,]
    f2s = f2[sig,]
    mat = cbind(sign(-f1s$log2.fold_change.)*-log(f1s$adjusted.p.value+1e-100), sign(-f2s$log2.fold_change.)*-log(f2s$adjusted.p.value+1e-100), rownames(f1s),sapply(strsplit(as.character(rownames(f1s)), ":"), `[`, 1))
    mat = as.data.frame(mat)
    mat$type = "autosomal"
    mat$type[which(mat$V4 == "chrX")] = "chrX"
    mat$type[which(mat$V4 == "chrY")] = "pseudoautosomal"

    mat$V1=as.numeric(mat$V1)
    mat$V2=as.numeric(mat$V2)
    g = ggplot(mat, aes(x=V1, y=V2, color=type)) + ggtitle(ct)+
    geom_point(size=1)+
    geom_point(data = subset(mat, type != "autosomal"), size = 1)
    print(g)
                
}

f1 = read.table(paste("male_diff/diff_csvs/diff_peaks_",ct,":Male_2vs18_iter.csv",sep = ""), sep = ",",header =T)
    f2 = read.table(paste("female_diff/diff_csvs/diff_peaks_",ct,":Female_2vs18_iter.csv",sep = ""), sep = ",",header =T)
    sig = c(f1$feature.name[which(f1$adjusted.p.value<0.01)], f2$feature.name[which(f2$adjusted.p.value<0.01)])
    length(sig)
    sig = unique(sig)
    length(sig)
    rownames(f1)= f1$feature.name
    rownames(f2)= f2$feature.name
    f1s = f1[sig,]
    f2s = f2[sig,]
    mat = cbind(-f1s$log2.fold_change., -f2s$log2.fold_change., rownames(f1s),sapply(strsplit(as.character(rownames(f1s)), ":"), `[`, 1))
    mat = as.data.frame(mat)
    mat$type = "autosomal"
    mat$type[which(mat$V4 == "chrX")] = "chrX"
    mat$type[which(mat$V4 == "chrY")] = "pseudoautosomal"

    mat$V1=as.numeric(mat$V1)
    mat$V2=as.numeric(mat$V2)
    g = ggplot(mat, aes(x=V1, y=V2, color=type)) + ggtitle(ct)+
    geom_point(size=1)+
    geom_point(data = subset(mat, type != "autosomal"), size = 1)
    print(g)
 

cor(mat$V1, mat$V2
   )

options(repr.plot.width=5, repr.plot.height=4)

ggplot(mat, aes(x=V2, y=V1, color=type)) + 
  ggtitle(ct) +
  geom_point(size=1) +
  geom_point(data = subset(mat, type != "autosomal"), size = 1) + theme_bw()+ 
  scale_color_manual(values = c("grey", "#c28aa5","blue")) +xlab("Female sig*dir")+ylab("Male sig*dir")


for(ct in c("Oligo_NN", "DG_Glut","L6_IT_CTX_Glut","L5_IT_CTX_Glut","L2-3_IT_PPP_Glut","CEA-BST_Gaba", "OPC_NN","STR_D12_Gaba","L2-3_IT_CTX_Glut","L2-3_IT_ENT_Glut")){
    f1 = read.table(paste("male_diff/diff_csvs/diff_peaks_",ct,":Male_2vs18_iter.csv",sep = ""), sep = ",",header =T)
    f2 = read.table(paste("female_diff/diff_csvs/diff_peaks_",ct,":Female_2vs18_iter.csv",sep = ""), sep = ",",header =T)
    sig = c(f1$feature.name[which(f1$adjusted.p.value<0.01)], f2$feature.name[which(f2$adjusted.p.value<0.01)])
    length(sig)
    sig = unique(sig)
    length(sig)
    rownames(f1)= f1$feature.name
    rownames(f2)= f2$feature.name
    f1s = f1[sig,]
    f2s = f2[sig,]
    mat = cbind(sign(-f1s$log2.fold_change.)*-log(f1s$adjusted.p.value+1e-100), sign(-f2s$log2.fold_change.)*-log(f2s$adjusted.p.value+1e-100), rownames(f1s),sapply(strsplit(as.character(rownames(f1s)), ":"), `[`, 1))
    mat = as.data.frame(mat)
    mat$type = "autosomal"
    mat$type[which(mat$V4 == "chrX")] = "chrX"
    mat$type[which(mat$V4 == "chrY")] = "chrY"

    mat$V1=as.numeric(mat$V1)
    mat$V2=as.numeric(mat$V2)
    g = ggplot(mat, aes(x=V1, y=V2, color=type)) + ggtitle(ct)+
    geom_point(size=1)+
    geom_point(data = subset(mat, type != "autosomal"), size = 1)
    print(g)
                
}

f1 = read.table("male_diff/diff_csvs/diff_peaks_STR_D12_Gaba:Male_2vs18_iter.csv", sep = ",",header =T)

f2 = read.table("female_diff/diff_csvs/diff_peaks_STR_D12_Gaba:Female_2vs18_iter.csv", sep = ",",header = T)

head(f2)

sig = c(f1$feature.name[which(f1$adjusted.p.value<0.01)], f2$feature.name[which(f2$adjusted.p.value<0.01)])

length(sig)

sig = unique(sig)
length(sig)

rownames(f1)= f1$feature.name

rownames(f2)= f2$feature.name

f1s = f1[sig,]
f2s = f2[sig,]


head(f1s)
head(f2s)

head(mat[which(mat$V4=="chrY"),])

mat = cbind(sign(-f1s$log2.fold_change.)*-log(f1s$adjusted.p.value+1e-100), sign(-f2s$log2.fold_change.)*-log(f2s$adjusted.p.value+1e-100), rownames(f1s),sapply(strsplit(as.character(rownames(f1s)), ":"), `[`, 1)
)

mat = as.data.frame(mat)

mat$type = "autosomal"
mat$type[which(mat$V4 == "chrX")] = "chrX"
mat$type[which(mat$V4 == "chrY")] = "chrY"
mat$V1=as.numeric(mat$V1)
mat$V2=as.numeric(mat$V2)



table(mat$type)
mat$type <- factor(mat$type, levels = (c("autosomal", "chrX", "chrY")))


ggplot(mat, aes(x=V1, y=V2, color=type)) + 
    geom_point(size=1)+
    geom_point(data = subset(mat, type != "autosomal"), size = 1)

ggplot(mat, aes(x=V1, y=V2, color=type)) + 
    geom_point(size=1)


