library(Seurat)


d12 = readRDS("2_sub.RDS")

d12 = FindClusters(d12, resolution = 4)

DimPlot(d12, label =T)

DimPlot(d12, label =T)

options(repr.plot.width=12, repr.plot.height=5)

VlnPlot(d12, "Drd1")
VlnPlot(d12, "Drd2")

options(repr.plot.width=7, repr.plot.height=8)

d1 = c(1,2,3,4,10,11,12,14,17,19,23,26,27,28,31,32,34,35,37,38,41,42,43,44,45,48,49,52,53,55,56,57,59,60,61,62)
d12$lab = "STR D2 Gaba"
d12$lab[which(d12$seurat_clusters%in%d1)] = "STR D1 Gaba"
d12$lab[which(d12$seurat_clusters%in%c(40,28,58,50,15))]="doublet"
DimPlot(d12, label =T, group.by = "lab")

head(d12@meta.data)

write.table(d12@meta.data, file = "D12MSN_split_meta.txt", sep ="\t")



d12$age_celltype <- gsub(" ", "_", paste(d12$age, d12$lab, sep = "_"))
d12$age_celltype_region <- gsub(" ", "_", paste(d12$age, d12$lab, d12$region, sep = "_"))


obj = d12
obj <- subset(obj, features = rownames(obj)[!(rownames(obj) %in% c("Malat1", "Lsamp", "Kcnip4"))])


setwd("../DEG_results_region")


Idents(obj) = "age_celltype_region"
for(r in unique(obj$region)){
# Loop through each unique cell type
    for (celltype in unique(obj$lab)[1:2]) {
      # Define unique identity combining age and cell type
      file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)),"--",r, ".csv", sep = "")

      if(file.exists(file_name)){
          next
     }
      ident1 <- gsub(" ", "_",paste("2mo", celltype,r, sep = "_"))
      ident2 <- gsub(" ", "_",paste("18mo", celltype, r,sep = "_"))
      if(length(which(obj$age_celltype_region ==ident1))>200 & length(which(obj$age_celltype_region ==ident2))>200){
          cat(ident1)
          de_results <- FindMarkers(
            obj,
            ident.1 = ident1,
            ident.2 = ident2,
            test.use = "MAST",
            logfc.threshold = 0.01, 
            min.pct = 0.01, 
          )
          head(de_results)
          # Save the results to a CSV file for each cell type
          write.csv(de_results, file = file_name, row.names = TRUE)
        }
    }
}


setwd("../DEG_results_.01_.01")


Idents(obj) = "age_celltype"
# Loop through each unique cell type
for (celltype in unique(obj$lab)) {
  # Define unique identity combining age and cell type
    
  ident1 <- gsub(" ", "_",paste("2mo", celltype, sep = "_"))
  ident2 <- gsub(" ", "_",paste("18mo", celltype, sep = "_"))

  de_results <- FindMarkers(
    obj,
    ident.1 = ident1,
    ident.2 = ident2,
    test.use = "MAST",
    logfc.threshold = 0.01, 
    min.pct = 0.01, 
  )
  head(de_results)
  # Save the results to a CSV file for each cell type
  file_name <- paste( gsub("/","-",gsub(" ", "_",celltype)), "csv", sep = ".")
  write.csv(de_results, file = file_name, row.names = TRUE)
}



