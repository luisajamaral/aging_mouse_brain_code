library(Seurat)

## Rachel celltypes:
meta_rach = read.csv("/home/lamaral/projects/brain_aging_mouse/analysis/integrtaion_adata/final_metadata.csv")

 unique(meta_rach$sample_name)

atac = read.table("../meta_with_subclassv3.txt")

atac = read.table("meta_from_RNA.txt")

head(atac)

atac = read.table("../meta")

#atac = atac[which(atac$batch== "Male"),]

atac$rep = paste("rep", atac$replicate, sep = "")

atac$agen=atac$age
atac$agen[which(atac$age == "2mo")] = "8wk"

atac$samplen = paste(atac$agen,  gsub("-", "+",sapply(strsplit(as.character(atac$sample_name), "_"), `[`, 2)), atac$rep)

which(unique(atac$samplen) %in% unique(meta_rach$sample_name))

atac$rsamp_barcode = paste(atac$samplen , ":", sapply(strsplit(as.character(atac$cell_id), ":"), `[`, 3), sep = "")

atac$rsamp_barcode[which(atac$batch=="Female")] = "Fem"

meta_rach$samp_barcode = paste(meta_rach$sample_name,":" ,meta_rach$barcode, sep = "")

head(meta_rach)


mat = match(atac$rsamp_barcode, meta_rach$samp_barcode, nomatch = NA)
atac$rachel_celltypes = meta_rach[mat,"new_celltypes"]

atac$Before = meta_rach[mat,"Before"]
atac$co_cluster_l2 = meta_rach[mat,"co_cluster_l2"]
atac$co_cluster_l1 = meta_rach[mat,"co_cluster_l1"]
atac$DissectionRegion = meta_rach[mat,"DissectionRegion"]

write.table(atac, file = "../meta_with_rachel_anno.txt", sep = "\t")

length(unique(atac$L2Annot.rough))






