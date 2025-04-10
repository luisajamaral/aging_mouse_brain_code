library(GenomicRanges)
library(ggplot2)
library(data.table)
library(UpSetR)
library(pheatmap)
library(Seurat)
library(enrichR)

dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2021_Human", "SynGO_2022")


load("links.RData")

split_and_convert <- function(strings) {
  # Split strings by "-"
  split_strings <- strsplit(strings, "-")
  
  # Extract components
  chr <- sapply(split_strings, `[`, 1)
  start <- sapply(split_strings, `[`, 2)
  end <- sapply(split_strings, `[`, 3)
  
  # Combine components into a data frame
  df <- data.frame(chr, start, end)
  
  return(df)
}


cts = gsub("diff_peaks_","",list.files("../h5ads_final/combined_diff/diff_csvs/",".csv"))
cts = gsub("_2vs18_iter.csv", "", cts)
cts

file.exists(paste("../female_RNA/DEG_results_.01_.01/",ct,".csv",sep = ""))

i

for(i in 43:length(cts)){
#ct = "DG_Glut"
ct = cts[i]
DARs = read.csv(paste("../h5ads_final/combined_diff/diff_csvs/diff_peaks_",ct,"_2vs18_iter.csv", sep = ""))
if(!file.exists(paste("../female_RNA/DEG_results_.01_.01/",ct,".csv",sep = ""))){
    next
}
degs = read.csv(paste("../female_RNA/DEG_results_.01_.01/",ct,".csv",sep = ""))


diffs = DARs
diffs$peak = gsub(":","-", diffs$feature.name)

diffs_up = diffs[which(diffs$log2.fold_change.<0.1),]
diffs_down = diffs[which(diffs$log2.fold_change.>0.1),]


dlinks = links[which(links$peak %in%diffs_down$peak ),]
ulinks = links[which(links$peak %in%diffs_up$peak ),]

cat("num down linked",length(dlinks), "\n")
cat("num up linked",length(ulinks), "\n")



nrow(degs)

rownames(degs) = degs$X
head(degs)

degs = degs[which(degs$p_val<0.05) ,]
data1 <- -degs[unique(ulinks$gene), "avg_log2FC"]
data2 <- -degs[unique(dlinks$gene), "avg_log2FC"]

# Create the histogram for data1 with transparent color
hist(data1, breaks = 40, col = rgb(1, 0, 0, 0.5), main = "Histogram of avg_log2FC", xlab = "avg_log2FC", ylab = "Frequency", ylim = c(0, max(c(hist(data1, breaks = 40, plot = FALSE)$counts, hist(data2, breaks = 40, plot = FALSE)$counts))))

# Add histogram for data2 with a different transparent color
hist(data2, breaks = 40, col = rgb(0, 0, 1, 0.5), add = TRUE)

# Add a legend
legend("topright", legend = c("ulinks", "dlinks"), fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))

up_genes = degs$X[which(degs$p_val<0.05 & degs$avg_log2FC<(0))]
down_genes = degs$X[which(degs$p_val<0.05 & degs$avg_log2FC>(0))]
all_genes = degs$X

length(up_genes)
length(down_genes)

enriched <- enrichr(up_genes[which(up_genes%in%unique(ulinks$gene))], dbs)
head(enriched[["GO_Biological_Process_2023"]],5)
write.table(enriched[["GO_Biological_Process_2023"]], file=paste(ct,"_up_linked_GO_BP.txt",sep = ""))
enriched <- enrichr(down_genes[which(down_genes%in%unique(dlinks$gene))], dbs)
head(enriched[["GO_Biological_Process_2023"]],5)
write.table(enriched[["GO_Biological_Process_2023"]], file=paste(ct,"_down_linked_GO_BP.txt",sep = ""))


overlapping_down= dlinks[which((dlinks$gene)%in%down_genes)]$peak
overlapping_up= ulinks[which((ulinks$gene)%in%up_genes)]$peak
overlapping_all= links[which((links$gene)%in%all_genes)]$peak
length(overlapping_all)

# Convert the vector of strings to a data frame
result_up <- split_and_convert(overlapping_up)
result_down <- split_and_convert(overlapping_down)
result_all <- split_and_convert(overlapping_all)

nrow(result_up)
nrow(result_down)
nrow(result_all)

write.table(result_up, paste(ct,"_overlap_linked_peaks_genes_up.bed", sep = ""), sep = "\t", row.names = F, col.names = F, quote=F)
write.table(result_down, paste(ct,"_overlap_linked_peaks_genes_down.bed", sep = ""), sep = "\t", row.names = F, col.names = F, quote=F)
write.table(result_all, paste(ct,"_overlap_linked_peaks_genes_all.bed", sep = ""), sep = "\t", row.names = F, col.names = F, quote=F)
    }

findMotifsGenome.pl Oligo_overlap_linked_peaks_genes_up.bed mm10  Oligo_overlap_linked_peaks_genes_up_bg_motifs -size 200  -p 15 -bg Oligo_overlap_linked_peaks_genes_all.bed
findMotifsGenome.pl D12_overlap_linked_peaks_genes_up.bed mm10  D12_overlap_linked_peaks_genes_up_bg_motifs -size 200  -p 15 -bg D12_overlap_linked_peaks_genes_all.bed
findMotifsGenome.pl L5_IT_CTX_Glut_overlap_linked_peaks_genes_down.bed mm10  L5_IT_CTX_Glut_overlap_linked_peaks_genes_down_bg_motifs -size 200  -p 15 -bg L5_IT_CTX_Glut_overlap_linked_peaks_genes_all.bed
findMotifsGenome.pl L5_IT_CTX_Glut_overlap_linked_nonDARpeaks_genes_up.bed mm10  L5_IT_CTX_Glut_overlap_linked_nonDARpeaks_genes_up.bed_bg_motifs -size 200  -p 15 -bg L5_IT_CTX_Glut_overlap_linked_peaks_genes_all.bed


findMotifsGenome.pl DG_Glut_overlap_linked_peaks_genes_down.bed mm10  DG_Glut_overlap_linked_peaks_genes_down_bg_motifs -size 200  -p 15 -bg DG_Glut_overlap_linked_peaks_genes_all.bed
findMotifsGenome.pl DG_Glut_overlap_linked_peaks_genes_up.bed mm10  DG_Glut_overlap_linked_peaks_genes_up_bg_motifs -size 200  -p 15 -bg DG_Glut_overlap_linked_peaks_genes_all.bed


findMotifsGenome.pl Microglia_NN_overlap_linked_peaks_genes_down.bed mm10  Microglia_NN_overlap_linked_peaks_genes_down_bg_motifs -size 200  -p 15 -bg Microglia_NN_overlap_linked_peaks_genes_all.bed


#findMotifsGenome.pl Oligo_overlap_linked_peaks_genes_down.bed mm10  Oligo_overlap_linked_peaks_genes_down_bg_motifs -size 200  -p 15 -bg Oligo_overlap_linked_peaks_genes_all.bed
#findMotifsGenome.pl Oligo_overlap_linked_peaks_genes_up.bed mm10  Oligo_overlap_linked_peaks_genes_up.bed_motifs -size 200  -p 15
