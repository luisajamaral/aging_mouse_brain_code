library(Seurat)

atac_meta = read.csv("../meta_bestcelltype.csv")

head(atac_meta)

atac_meta = atac_meta[which(atac_meta$age == "2mo"),]
atac_meta = atac_meta[which(atac_meta$batch == "Male"),]


nrow(atac_meta)

samps = unique(atac_meta$sample_name)

rep = substr(samps, nchar(samps), nchar(samps))

reg = sapply(strsplit(as.character(samps), "_"), `[`, 2)

rep

akey = cbind(samps, reg, rep)

akey = as.data.frame(akey)
rownames(akey) = akey$samps

akey$samp_id = paste(akey$reg, "_",akey$rep, sep = "")

which(akey$samp_id %in% key$samp_id)

head(akey)



meta = "/projects/ps-renlab2/szu/projects/CEMBA2/supple.02.annotation.all.in.one.meta/mba.whole.cell.meta.v9.7.rds"
m = readRDS(meta)

parts <- sapply(strsplit(unique(m$sample), "_"), function(x) x[2])
count_vector <- ave(parts, parts, FUN = seq_along)

key = cbind(unique(m$sample), parts, count_vector)

key = as.data.frame(key)

key$samp_id = paste(key$parts , "_",key$count_vector, sep = "")

rownames(key) = unique(m$sample)

key

m$samp_id = key[paste(m$sample), "samp_id"]

head(m)

m$sampid_barcode = paste(m$samp_id , ":", m$barcode, sep = "")

rownames(m) = m$sampid_barcode

colnames(m)

atac = read.csv("../meta_bestcelltype.csv")
atac$samp_id = akey[paste(atac$sample_name), "samp_id"]
sampid_barcode = paste(atac$samp_id, ":", sapply(strsplit(as.character(atac$cell_id), ":"), `[`, 3), sep = "")
atac$sampid_barcode = sampid_barcode

atac_meta$samp_id = akey[paste(atac_meta$sample_name), "samp_id"]
sampid_barcode = paste(atac_meta$samp_id, ":", sapply(strsplit(as.character(atac_meta$cell_id), ":"), `[`, 3), sep = "")
atac_meta$sampid_barcode = sampid_barcode
rownames(atac_meta) = atac_meta$sampid_barcode
head(atac_meta)



mat = match(atac$sampid_barcode, rownames(m), nomatch = NA)
atac$L2Annot.rough = m[mat,"L2Annot.rough"]

write.table(atac, file = "../meta_with_subclassv3.txt", sep = "\t")

atac_meta$subclass_label_v3 = m[mat,"subclass_label_v3"]

head(atac)



atac$subclass_label_v3 = m[paste(atac$sampid_barcode),"subclass_label_v3"]

length(unique(atac$L2Annot.rough))

length(unique(atac$subclass_label_v3))

table(atac$L2Annot.rough)


