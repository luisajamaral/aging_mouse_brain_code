library(Seurat)
library(ggplot2)

obj = readRDS("RNA_final_SCT.RDS")

Idents(obj) = "celltype_final"

DefaultAssay(obj)="SCT"

obj

DimPlot(obj, label = T)+NoLegend()

FeaturePlot(obj, "percent.ribo")

data = obj
data[rownames(data) != "Malat1",] -> data.nomalat

apply(
  data.nomalat@assays$RNA@counts,
  2,
  max
) -> data.nomalat$largest_count

apply(
  data.nomalat@assays$RNA@counts,
  2,
  which.max
) -> data.nomalat$largest_index

rownames(data.nomalat)[data.nomalat$largest_index] -> data.nomalat$largest_gene

100 * data.nomalat$largest_count / data.nomalat$nCount_RNA -> data.nomalat$percent.Largest.Gene

data.nomalat$largest_gene -> data$largest_gene
data.nomalat$percent.Largest.Gene -> data$percent.Largest.Gene

rm(data.nomalat)

library(Matrix)
gc()
# Exclude MALAT1 (without making a full copy)
data.nomalat <- data
data.nomalat <- subset(data.nomalat, features = rownames(data.nomalat) != "Malat1")

# Use sparse matrix methods to find max values efficiently
largest_counts <- apply(data.nomalat@assays$RNA@counts, 2, max)
largest_indices <- apply(data.nomalat@assays$RNA@counts, 2, which.max)

# Get gene names corresponding to the max values
largest_genes <- rownames(data.nomalat)[largest_indices]

# Compute percent of total RNA contributed by the largest gene
percent_largest_gene <- 100 * largest_counts / data.nomalat$nCount_RNA

# Assign to metadata without creating a new object in memory
data$largest_gene <- largest_genes
data$percent.Largest.Gene <- percent_largest_gene

# Clean up temporary variables to free memory
rm(data.nomalat, largest_counts, largest_indices, largest_genes, percent_largest_gene)
gc()  # Force garbage collection


head(data.nomalat)

VlnPlot(data, group.by = "age",pt.size = 0,features=c("nCount_RNA","percent.mt", "percent.ribo","percent.Largest.Gene"))


as_tibble(
  data[[]],
  rownames="Cell.Barcode"
) -> qc.metrics

head(qc.metrics)

qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Example of plotting QC metrics") +
  geom_hline(yintercept = 750) +
  geom_hline(yintercept = 2000) 

qc.metrics %>%
  arrange(percent.ribo) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.ribo)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
  ggtitle("Example of plotting QC metrics") +
  geom_hline(yintercept = 750) +
  geom_hline(yintercept = 2000) 

qc.metrics %>%
  mutate(complexity=log10(nFeature_RNA) / log10(nCount_RNA))  -> qc.metrics

lm(log10(qc.metrics$nFeature_RNA)~log10(qc.metrics$nCount_RNA)) -> complexity.lm

complexity.lm

qc.metrics %>%
  mutate(
    complexity_diff = log10(nFeature_RNA) - ((log10(qc.metrics$nCount_RNA)*complexity.lm$coefficients[2])+complexity.lm$coefficients[1])
  ) -> qc.metrics

qc.metrics %>%
  ggplot(aes(x=complexity_diff)) +
  geom_density(fill="yellow")

min(c(max(qc.metrics$complexity_diff),0-min(qc.metrics$complexity_diff))) -> complexity_scale

qc.metrics %>%
  mutate(complexity_diff=replace(complexity_diff,complexity_diff< -0.1,-0.1)) %>%
  ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity_diff)) +
  geom_point(size=0.5) +
  geom_abline(slope=complexity.lm$coefficients[2], intercept = complexity.lm$coefficients[1]) +
  scale_colour_gradient2(low="blue2",mid="grey",high="red2")

ggplot(qc.metrics, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=complexity_diff)) +
  geom_point(size=0.5, alpha=0.6) +
  geom_smooth(method="loess", se=FALSE, color="black") +  
  scale_colour_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_minimal() +
  coord_fixed(ratio=1) +
  labs(title="Gene Complexity vs UMI Counts",
       x="log10(UMI Counts)",
       y="log10(Unique Features)",
       colour="Complexity Score")


qc.metrics %>%
  ggplot(aes(percent.ribo)) + 
  geom_histogram(binwidth = 0.5, fill="yellow", colour="black") +
  ggtitle("Distribution of Percentage Mitochondrion") +
  geom_vline(xintercept = 10)

qc.metrics %>%
  group_by(largest_gene) %>%
  count() %>%
  arrange(desc(n)) -> largest_gene_list

largest_gene_list

ggplot(mapping = aes(obj@assays$SCT@data["Gapdh",])) + 
  geom_histogram(binwidth = 0.05, fill="yellow", colour="black") + 
  ggtitle("GAPDH expression")

library(dplyr)
library(tidyr)
library(tidyverse)


as.tibble(
  obj@assays$SCT$data[,1:100]
) %>%
  pivot_longer(
    cols=everything(),
    names_to="cell",
    values_to="expression"
  ) %>%
  ggplot(aes(x=expression, group=cell)) +
  geom_density() +
  coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))

NormalizeData(data, normalization.method = "CLR", margin = 2) -> data


obj$celltype_final[which(obj$best_celltype_fixed == "IOL")] = "IOL NN"

obj$age_celltype <- gsub(" ", "_", paste(obj$age, obj$celltype_final, sep = "_"))


obj$age_celltype=gsub("Astro-TE","Astro",obj$age_celltype)
obj$age_celltype=gsub("Astro-NT","Astro",obj$age_celltype)


obj$age_celltype_region <- gsub(" ", "_", paste(obj$age_celltype, obj$region, sep = "_"))


#obj$age_celltype <- gsub(" ", "_", paste(obj$age, obj$celltype_final, sep = "_"))

obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))






#dir.create("DEG_results_region", showWarnings = FALSE)
setwd("DEG_results_region_latent_orig.ident")

dir.create("DEG_results_region_MAST", showWarnings = FALSE)
setwd("DEG_results_region_MAST")

unique(obj$celltype_final)

obj$celltype_final[which(obj$best_celltype_fixed=="IOL")] = "IOL NN"


obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) = "age_celltype_region"
for(r in unique(obj$region)){
# Loop through each unique cell type
    #for (celltype in unique(obj$celltype_final)) {
    for (celltype in unique(obj$celltype_final)) {
      # Define unique identity combining age and cell type
      file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"--",r, ".csv", sep = "")

      if(file.exists(file_name)){
          cat(file_name, "file exists")
          next
     }
      ident1 <- gsub(" ", "_",paste("2mo", celltype,r, sep = "_"))
      ident2 <- gsub(" ", "_",paste("18mo", celltype, r,sep = "_"))
      if(length(which(obj$age_celltype_region ==ident1))>20 & length(which(obj$age_celltype_region ==ident2))>20){
          cat("Running", file_name)
          de_results <- FindMarkers(
            obj,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            logfc.threshold = 0.01, 
            min.pct = 0.01, latent.vars = c('percent.mt')
          )
          head(de_results)
          # Save the results to a CSV file for each cell type
          write.csv(de_results, file = file_name, row.names = TRUE)
        }
    }
}


unique(obj$region)


obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) = "age_celltype_region"
for(r in unique(obj$region)){
# Loop through each unique cell type
    #for (celltype in unique(obj$celltype_final)) {
    for (celltype in unique(obj$celltype_final)) {
      # Define unique identity combining age and cell type
      file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"--",r, "_9vs18.csv", sep = "")

      if(file.exists(file_name)){
          cat(file_name, "file exists")
          next
     }
      ident1 <- gsub(" ", "_",paste("9mo", celltype,r, sep = "_"))
      ident2 <- gsub(" ", "_",paste("18mo", celltype, r,sep = "_"))
      if(length(which(obj$age_celltype_region ==ident1))>20 & length(which(obj$age_celltype_region ==ident2))>20){
          cat("Running", file_name)
          de_results <- FindMarkers(
            obj,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            logfc.threshold = 0.01, 
            min.pct = 0.01, latent.vars = c('percent.mt')
          )
          head(de_results)
          # Save the results to a CSV file for each cell type
          write.csv(de_results, file = file_name, row.names = TRUE)
        }
    }
}


num_genes_to_exclude <- 10  # Adjust the number as needed
top_expressed_genes <- rownames(
  head(
    x = obj@assays$RNA@counts[order(rowMeans(obj@assays$RNA@counts, na.rm = TRUE), decreasing = TRUE), ],
    n = num_genes_to_exclude
  )
)
top_expressed_genes

obj <- subset(obj, features = rownames(obj)[!(rownames(obj) %in% top_expressed_genes[1:3])])


obj <- subset(obj, features = rownames(obj)[!(rownames(obj) %in% c("Malat1", "Lsamp", "Kcnip4"))])


obj

r = "HCP"
celltype = "Oligo_NN"

which(obj$age_celltype_region ==ident1)

table(obj$age_celltype_region)

DefaultAssay(obj)="RNA"

celltype = "Oligo NN"
r = "ENT"
file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"--",r, "_logfc.25.csv", sep = "")

      if(file.exists(file_name)){
          cat(file_name, "file exists")
          next
     }
      ident1 <- gsub(" ", "_",paste("2mo", celltype,r, sep = "_"))
      ident2 <- gsub(" ", "_",paste("18mo", celltype, r,sep = "_"))
      if(length(which(obj$age_celltype_region ==ident1))>50 & length(which(obj$age_celltype_region ==ident2))>50){
          cat("Running", file_name)
          de_results <- FindMarkers(
            obj,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            logfc.threshold = 0.25, 
            min.pct = 0.01, latent.vars = 'orig.ident'
          )
          head(de_results)
          # Save the results to a CSV file for each cell type
          write.csv(de_results, file = file_name, row.names = TRUE)
          }

unique(obj$region)

write.table(obj@meta.data, file = "../meta.txt")

table(obj$class)
obj$age_clade_region=  gsub(" ", "_", paste(obj$age, obj$class, obj$region,sep = "_"))
table(obj$age_clade_region)

ident2

setwd("../DEG_results_regions_category/")


Idents(obj) = "age_clade_region"
for(r in unique(obj$region)){
# Loop through each unique cell type
    for (celltype in unique(obj$class)) {
    #for (celltype in "Astro_NN") {
      # Define unique identity combining age and cell type
      file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"--",r, ".csv", sep = "")

      if(file.exists(file_name)){
          cat(file_name, "file exists")
          next
     }
      ident1 <- gsub(" ", "_",paste("2mo", celltype,r, sep = "_"))
      ident2 <- gsub(" ", "_",paste("18mo", celltype, r,sep = "_"))
      if(length(which(obj$age_clade_region ==ident1))>50 & length(which(obj$age_clade_region ==ident2))>50){
          cat("Running", file_name)
          de_results <- FindMarkers(
            obj,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            logfc.threshold = 0.0, 
            min.pct = 0.001, latent.vars = c('percent.mt')
          )
          head(de_results)
          # Save the results to a CSV file for each cell type
          write.csv(de_results, file = file_name, row.names = TRUE)
        }
    }
}


setwd("DEG_results_latent_rep_mito/")

obj$celltype_final[which(obj$final_clusters_final=="IOL")] = "IOL NN"

length(which(obj$celltype_final=="IOL NN"))

unique(obj$celltype_final)

setwd("../DEG_results_latent_rep_mito_together/")

setwd("../progenitor_DEGs/")
obj$age_celltype <- gsub(" ", "_", paste(obj$age, obj$celltype_final, sep = "_"))
obj$age_celltype[which(obj$age_celltype=="18mo_IOL_NN")] = "Old_IOL_NN"
obj$age_celltype[which(obj$age_celltype=="9mo_IOL_NN")] = "Old_IOL_NN"

obj$age_celltype[which(obj$age_celltype=="18mo_OB-STR-CTX_Inh_IMN")] = "Old_OB-STR-CTX_Inh_IMN"
obj$age_celltype[which(obj$age_celltype=="9mo_OB-STR-CTX_Inh_IMN")] = "Old_OB-STR-CTX_Inh_IMN"

obj$age_celltype[which(obj$age_celltype=="18mo_DG-PIR_Ex_IMN")] = "Old_DG-PIR_Ex_IMN"
obj$age_celltype[which(obj$age_celltype=="9mo_DG-PIR_Ex_IMN")] = "Old_DG-PIR_Ex_IMN"

Idents(obj) = "age_celltype"
# Loop through each unique cell type
#for (celltype in unique(obj$celltype_final)) {
for (celltype in c("DG-PIR_Ex_IMN", "OB-STR-CTX_Inh_IMN" ,"IOL_NN" )) {
  # Define unique identity combining age and cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"2vsOld_default", ".csv", sep = "")

  if(file.exists(file_name)){
      cat(file_name, "file exists")
      next
 }
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("Old", celltype, sep = "_"))
  if(length(which(obj$age_celltype ==ident1))>20 & length(which(obj$age_celltype ==ident2))>20){
      cat("Running", file_name)
      de_results <- FindMarkers(
        obj,
        ident.1 = ident1,
        ident.2 = ident2,
        logfc.threshold = 0.05, 
        min.pct = 0.01, latent.vars = c('rep','percent.mt')
      )
      head(de_results)
      # Save the results to a CSV file for each cell type
      write.csv(de_results, file = file_name, row.names = TRUE)
    }
}


ident2

setwd("../progenitor_DEGs/")
obj$age_celltype <- gsub(" ", "_", paste(obj$age, obj$celltype_final, sep = "_"))
obj$age_celltype[which(obj$age_celltype=="18mo_IOL_NN")] = "Old_IOL_NN"
obj$age_celltype[which(obj$age_celltype=="9mo_IOL_NN")] = "Old_IOL_NN"

obj$age_celltype[which(obj$age_celltype=="18mo_OB-STR-CTX_Inh_IMN")] = "Old_OB-STR-CTX_Inh_IMN"
obj$age_celltype[which(obj$age_celltype=="9mo_OB-STR-CTX_Inh_IMN")] = "Old_OB-STR-CTX_Inh_IMN"

obj$age_celltype[which(obj$age_celltype=="18mo_DG-PIR_Ex_IMN")] = "Old_DG-PIR_Ex_IMN"
obj$age_celltype[which(obj$age_celltype=="9mo_DG-PIR_Ex_IMN")] = "Old_DG-PIR_Ex_IMN"

Idents(obj) = "age_celltype"
# Loop through each unique cell type
#for (celltype in unique(obj$celltype_final)) {
for (celltype in c("DG-PIR_Ex_IMN", "OB-STR-CTX_Inh_IMN" ,"IOL_NN" )) {
  # Define unique identity combining age and cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"2vsOld", ".csv", sep = "")

  if(file.exists(file_name)){
      cat(file_name, "file exists")
      next
 }
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("Old", celltype, sep = "_"))
  if(length(which(obj$age_celltype ==ident1))>20 & length(which(obj$age_celltype ==ident2))>20){
      cat("Running", file_name)
      de_results <- FindMarkers(
        obj,
        ident.1 = ident1,
        ident.2 = ident2,
        test.use = "MAST",
        logfc.threshold = 0.05, 
        min.pct = 0.01, latent.vars = c('rep','percent.mt')
      )
      head(de_results)
      # Save the results to a CSV file for each cell type
      write.csv(de_results, file = file_name, row.names = TRUE)
    }
}


obj$age_celltype <- gsub(" ", "_", paste(obj$age, obj$celltype_final, sep = "_"))


Idents(obj) = "age_celltype"
# Loop through each unique cell type
    #for (celltype in unique(obj$celltype_final)) {
    for (celltype in c("IOL_NN" )) {
      # Define unique identity combining age and cell type
      file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"2vsOld", ".csv", sep = "")

      if(file.exists(file_name)){
          cat(file_name, "file exists")
          next
     }
      ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
      ident2 <- gsub(" ", "_",paste("Old", celltype, sep = "_"))
      if(length(which(obj$age_celltype ==ident1))>20 & length(which(obj$age_celltype ==ident2))>20){
          cat("Running", file_name)
          de_results <- FindMarkers(
            obj,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            logfc.threshold = 0.05, 
            min.pct = 0.01, latent.vars = c('rep','percent.mt')
          )
          head(de_results)
          # Save the results to a CSV file for each cell type
          write.csv(de_results, file = file_name, row.names = TRUE)
        }
    }








obj$age_celltype_region <- gsub(" ", "_", paste(obj$age, obj$celltype_final, obj$region, sep = "_"))

head(obj$age_celltype_region)
Idents(obj) = "age_celltype_region"
for(r in unique(obj$region)){
# Loop through each unique cell type
    #for (celltype in unique(obj$celltype_final)) {
    for (celltype in "IOL_NN") {
      # Define unique identity combining age and cell type
      file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"--",r, ".csv", sep = "")

      if(file.exists(file_name)){
          cat(file_name, "file exists")
          next
     }
      ident1 <- gsub(" ", "_",paste("2mo", celltype,r, sep = "_"))
      ident2 <- gsub(" ", "_",paste("18mo", celltype, r,sep = "_"))
      if(length(which(obj$age_celltype_region ==ident1))>20 & length(which(obj$age_celltype_region ==ident2))>20){
          cat("Running", file_name)
          de_results <- FindMarkers(
            obj,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            logfc.threshold = 0.05, 
            min.pct = 0.01, latent.vars = c('rep','percent.mt')
          )
          head(de_results)
          # Save the results to a CSV file for each cell type
          write.csv(de_results, file = file_name, row.names = TRUE)
        }
    }
}


ident1

file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"--",r, ".csv", sep = "")
write.csv(de_results, file = file_name, row.names = TRUE)


head(de_results)

setwd("../DEG_results_latent_rep_mito_minpct005")

unique(obj$celltype_final)

obj$age_celltype <- gsub(" ", "_", paste(obj$age, obj$celltype_final, sep = "_"))



Idents(obj) = "age_celltype"
# Loop through each unique cell type
for (celltype in unique(obj$celltype_final)) {
#for (celltype in c("Astro-NT NN","DG Glut")) {
  # Define unique identity combining age and cell type
    
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("18mo", celltype, sep = "_"))

  de_results <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    logfc.threshold = 0.1, 
    min.pct = 0.01, latent.vars = c('rep','percent.mt')
  )
  head(de_results)
  # Save the results to a CSV file for each cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"csv", sep = ".")
  write.csv(de_results, file = file_name, row.names = TRUE)
}


setwd("../DEG_results_latent_rep_mito_together_wilcox_.001pct_0logfc")


Idents(obj) = "age_celltype"
# Loop through each unique cell type
for (celltype in unique(obj$celltype_final)) {
#for (celltype in c("Astro-NT NN","DG Glut")) {
  # Define unique identity combining age and cell type
    
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("18mo", celltype, sep = "_"))

  de_results <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "wilcox",
    logfc.threshold = 0, 
    min.pct = 0.001, latent.vars = c('rep','percent.mt')
  )
  head(de_results)
  # Save the results to a CSV file for each cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"csv", sep = ".")
  write.csv(de_results, file = file_name, row.names = TRUE)
}


curr = read.csv("Oligo_NN.csv", sep = ",")
head(curr$X)

nrow(curr)
nrow(curr[which(curr$`pct.2`>0.005 | curr$`pct.1` >0.005 ),])

write.table(paste(curr$X), file = "Oligo_genes.txt", sep = "\n", quote = F, row.names = F, col.names = F)

curr$avg_pct = (curr$pct.2+curr$pct.1)/2

curr[order(curr$avg_pct)[1:10],]

files = list.files(".", ".csv")
files

files = list.files(".", ".csv")

for( f in files) {
    curr = read.csv(f, sep = ",")
    write.table(paste(curr$X), file = gsub(".csv", "_genes.txt", f), sep = "\n", quote = F, row.names = F, col.names = F)
    }

head(de_results)

obj

write.table(obj@meta.data

meta = obj@meta.data

head(meta)

d12  = read.table("../../female_RNA/subclustering3/D12MSN_split_meta.txt")

head(d12)

mat = match(meta$barcode,d12$barcode)

head(mat)

obj$celltype_age = paste(obj$celltype_final, obj$age, sep = "_")
Idents(obj) = "celltype_age"
av = AverageExpression(obj)

save(av , file = "../female_RNA/celltype_age_average_expression.RData")

obj


