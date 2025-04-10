library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(qlcMatrix)
library(JASPAR2020)
library(data.table)
library(TFBSTools)

o = ReadMtx("peak_mat_129276_by_881231_high_quality.mtx", "barcodes.tsv","features.tsv",cell.column = 1,
  feature.column = 1,cell.sep = "\n",feature.sep = "\n")

obj  = readRDS("subsampled_high_quality_RNA_seurat.RDS")

ATAC= CreateAssayObject(o)
obj[["peaks"]]=ATAC

# Function to parse locations and create data.table
parse_location <- function(location) {
  parts <- unlist(strsplit(location, "[:-]"))
  chr <- parts[1]
  start <- as.numeric(parts[2])
  end <- as.numeric(parts[3])
  return(data.table(chr = chr, start = start, end = end))
}

# Apply function to each element of the locations vector
dt_list <- lapply(locations, parse_location)

# Combine into a single data.table
locations_dt <- rbindlist(dt_list)

# Convert data.table to GRanges object
granges_obj <- GRanges(
  seqnames = locations_dt$chr,
  ranges = IRanges(start = locations_dt$start, end = locations_dt$end)
)

# View the Granges object
print(granges_obj)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

peaks = CreateChromatinAssay(obj[["peaks"]]@counts,ranges = granges_obj, annotation = annotation )

obj[["peaks"]] = peaks

obj

DefaultAssay(obj) = "peaks"

obj

obj <- RegionStats(obj, genome = BSgenome.Mmusculus.UCSC.mm10)


lp = LinkPeaks(obj, peak.assay = "peaks", expression.assay = "RNA", expression.slot = "data",peak.slot = "counts")

lp@assays$peaks$counts[1:5,1:5]

lp@assays$RNA$counts[1:5,1:5]

links = Links(lp)

fwrite(as.data.frame(links), "Links_Granges.txt")

links

save(links, file = "links.RData")

l = as.data.frame(links)
l[order(l$pvalue)[1:30],]

obj = readRDS("lp_all.RDS")

load("links.RData")

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

obj

DefaultAssay(obj) <- "peaks"
obj <- FindTopFeatures(obj, min.cutoff = 5)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

options(repr.plot.width=26, repr.plot.height=12)

DimPlot(obj,group.by = "celltype_final", label = T)

options(repr.plot.width=8, repr.plot.height=7)

DimPlot(obj,group.by = "age", label = F)

obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(obj, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()


obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(obj) <- 'chromvar'


mn = as.data.frame(obj@assays$peaks@motifs@motif.names)
mn = (t(mn))

mn = as.data.frame(mn)
mn$id = rownames(mn)
#rownames(mn) = mn[,1]

mn[grep("MA1619.1",mn$id),]

mn[grep("JUN",mn$V1),]

obj$age_rep = paste(obj$age, obj$rep)

DefaultAssay(obj)="RNA"

options(repr.plot.width=15, repr.plot.height=10)

FeaturePlot(
  object = obj,
  features = c("Rfx1","Rfx2"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, split.by = "age"
)

obj$celltype_age = paste(obj$celltype_final, ":", obj$age)
Idents(obj) = "celltype_age"
Idents(obj)[1:10]

setwd("motif_chromvar/")

which(unique(obj$celltype_final)==ct)

ct = "Oligo NN"
differential.activity <- FindMarkers(
  object = obj,
  ident.1 = paste(ct,' : 2mo',sep = ""),
  ident.2 = paste(ct,' : 18mo',sep = ""),
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
differential.activity$motif = mn[paste(rownames(differential.activity)),1]

head(differential.activity)

Idents(obj) = "celltype_age"

for (ct in unique(obj$celltype_final)[71:length(unique(obj$celltype_final))]) {


differential.activity <- FindMarkers(
  object = obj,
  ident.1 = paste(ct,' : 18mo',sep = ""),
  ident.2 = paste(ct,' : 2mo',sep = ""),
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
differential.activity$motif = mn[paste(rownames(differential.activity)),1]
ct = gsub("/","-",ct)
write.table(differential.activity, file = paste(ct, "_18v2_motifs.txt", sep = ""))    
    
options(repr.plot.width=14, repr.plot.height=5)
pdf(paste(ct, "_18v2_motifs.pdf", sep = ""), height =5, width = 13)

p1 = MotifPlot(
  object = obj,
  motifs = head(rownames(differential.activity[which(differential.activity$avg_diff>0),])),
  assay = 'peaks'
)

p2= MotifPlot(
  object = obj,
  motifs = head(rownames(differential.activity[which(differential.activity$avg_diff<0),])),
  assay = 'peaks'
)


p3 = FeaturePlot(
  object = obj,
  features = head(rownames(differential.activity[which(differential.activity$avg_diff>0),]),1),
  min.cutoff = 'q10',
  max.cutoff = 'q80',
  pt.size = 0.1, split.by = "age"
)

p4 = FeaturePlot(
  object = obj,
  features = head(rownames(differential.activity[which(differential.activity$avg_diff<0),]),1),
  min.cutoff = 'q10',
  max.cutoff = 'q80',
  pt.size = 0.1, split.by = "age"
)
    
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
}

differential.activity$motif = mn[paste(rownames(differential.activity)),1]
differential.activity[1:10,]

MotifPlot(
  object = obj,
  motifs = head(rownames(differential.activity[which(differential.activity$avg_diff>0),])),
  assay = 'peaks'
)

MotifPlot(
  object = obj,
  motifs = head(rownames(differential.activity[which(differential.activity$avg_diff<0),])),
  assay = 'peaks'
)

Idents(obj)= "celltype_final"

VlnPlot(
  object = obj,
  features = "Rfx1",
  pt.size = 0, split.by = "age", group.by = "region"
)

saveRDS(obj, file = "obj.RDS")

obj@meta.data$age = factor(obj@meta.data$age, levels = c("2mo", "9mo", "18mo"))
obj@meta.data$region = factor(obj@meta.data$region, levels = c("AMY","RLP","HCA","HCP","FC","ENT","CP","NAC"))
obj$celltype_region_age = paste(obj$celltype_final, obj$region, obj$age)
Idents(obj) = "celltype_region_age"
av = AverageExpression(obj)

Idents(obj) = "celltype_age"


load("links.RData")


DARs = read.csv("../h5ads_final/combined_diff/diff_csvs/diff_peaks_Oligo_NN_2vs18_iter.csv")

head(DARs)

degs = read.csv("../female_RNA/DEG_results_.01_.01/Oligo_NN.csv")

head(degs)

degs = degs[order(degs$avg_log2FC),]
genes = degs$X[which(degs$p_val_adj<0.01)]

head(genes)

Idents(obj)="celltype_final"
fm=FindAllMarkers(obj,only.pos = T, logfc.threshold = .5, assay = "RNA")

head(fm)


fm %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>%
    ungroup() -> top10

DefaultAssay(obj)="RNA"
DoHeatmap(obj, features = top10$gene) + NoLegend()


av = AverageExpression(obj)

ct = "Oligo NN"

oligo = av$RNA[genes,grep(ct, colnames(av$RNA))]


linktable = as.data.frame(links)
linked_peaks=linktable[which(linktable$gene %in%genes),]
linked_peaks$gene = factor(linked_peaks$gene, levels = genes)
linked_peaks = linked_peaks[order(linked_peaks$gene, decreasing = F),]
head(linked_peaks)

linked_peaks = linked_peaks[,"peak"]



DARs = read.csv("../h5ads_final/combined_diff/diff_csvs/diff_peaks_Oligo_NN_2vs18_iter.csv")
diffs = DARs[which(DARs$p.value<0.05),]
diffs$peak = gsub(":","-", diffs$feature.name)


oligo_atac = av$peaks[linked_peaks,grep(ct, colnames(av$peaks))]
oligo_atac = oligo_atac[which(rownames(oligo_atac)%in%diffs$peak),]

#oligo_atac = oligo_atac[-which(rowSums(oligo_atac)<1),]
dim(oligo_atac)

age <- colnames(oligo)
age[grep("18mo", age)] = "18mo"
age[grep("2mo", age)] = "2mo"
age[grep("9mo", age)] = "9mo"
region <- colnames(oligo)
region[grep("AMY", region)] = "AMY"
region[grep("RLP", region)] = "RLP"
region[grep("HCA", region)] = "HCA"
region[grep("FC", region)] = "FC"
region[grep("NAC", region)] = "NAC"
region[grep(" CP", region)] = "CP"
region[grep("HCP", region)] = "HCP"


annotation_col = data.frame(
    age = age,
    region = region
    
  )
annotation_col$age =factor(annotation_col$age, levels = c("2mo","9mo", "18mo"))
rownames(annotation_col) = colnames(oligo)


pheatmap(oligo[,order(annotation_col$age)], show_rownames = F,cutree_rows = 6,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), cluster_cols =F,cluster_rows =F,
         main = length(genes), annotation_col=annotation_col)



pheatmap(oligo_atac[,order(annotation_col$age)], show_rownames = F,cutree_rows = 6,scale = "row",color=colorRampPalette(c("navy", "white", "red"))(40), cluster_cols =F,cluster_rows =F,
         main = nrow(oligo_atac), annotation_col=annotation_col)



min(rowSums(oligo_atac))
oligo_atac = oligo_atac[-which(rowSums(oligo_atac)==0),]

nrow(oligo_atac)

which(rowSums(oligo_atac)==0)

d12 = subset(obj, cells = which(obj$celltype_final=="Oligo NN"))

iol = subset(d12, cells = which(d12$seurat_clusters==46))

iol

iol <- FindMultiModalNeighbors(
  object = iol,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
iol <- RunUMAP(
  object = iol,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(iol, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()


Idents(iol)= "age"
DefaultAssay(iol) ="chromvar"
differential.activity <- FindMarkers(
  object = iol,
  ident.1 = paste('18mo',sep = ""),
  ident.2 = paste('2mo',sep = ""),
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

differential.activity$motif = mn[paste(rownames(differential.activity)),1]
head(differential.activity)

DimPlot(d12, label = TRUE, repel = TRUE, reduction = "umap") 


options(repr.plot.width=7, repr.plot.height=5)

DimPlot(d12, label = F, repel = TRUE, reduction = "umap", group.by="age")
DimPlot(d12, label = T, repel = TRUE, reduction = "umap", group.by="seurat_clusters")


DefaultAssay(d12)= "RNA"
FeaturePlot(
  object = d12,
  features = c("Mal", "Oligo1"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)

table(d12$age)

options(repr.plot.width=15, repr.plot.height=13)

FeaturePlot(
    
  object = d12,
  features = c("MA1564.1","Ctcf","MA0139.1"),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1, split.by = "age"
    
)



VlnPlot(
  object = d12,
  features = "MA1128.1",
  pt.size = 0.01, split.by = "rep", group.by = "age"
)


