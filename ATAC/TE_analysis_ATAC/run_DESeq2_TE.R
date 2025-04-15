# DESeq2 TE Analysis at Region and Cell Type Levels

library(DESeq2)
library(tidyverse)
library(data.table)
library(pheatmap)

# ==== 1. Load count matrix and metadata ====
count_data <- fread("allTE.feature.counts.txt", skip=1)
count_data <- as.data.frame(count_data)

# Format rownames using the first column (gene/feature names)
rownames(count_data) <- count_data$Geneid
count_data <- count_data[, -(1:6)]  # Remove annotation columns

# Clean column names (remove BAM suffix)
colnames(count_data) <- gsub(".bam", "", colnames(count_data))
colnames(count_data) <- gsub("../female_split_bams/", "", colnames(count_data))

# ==== 2. Load sample metadata (region, age, and cell type) ====
# You should modify this metadata file accordingly
meta <- read.csv("sample_metadata.csv")
meta$sample <- gsub(".bam", "", meta$sample)

# Match and reorder metadata to count columns
meta <- meta[match(colnames(count_data), meta$sample), ]
stopifnot(all(meta$sample == colnames(count_data)))

# ==== 3. Run DESeq2: Region-level ====
meta$group <- factor(paste(meta$region, meta$age, sep="_"))
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = meta, design = ~ region + age)
dds <- DESeq(dds)
res <- results(dds, contrast = c("age", "18mo", "2mo"))

# Save results
write.csv(as.data.frame(res), "TE_deseq_region_18vs2.csv")

# ==== 4. Run DESeq2: Celltype-level with region as batch ====
meta$group_ct <- factor(paste(meta$celltype, meta$age, sep="_"))

# Iterate by cell type
unique_cts <- unique(meta$celltype)
for (ct in unique_cts) {
  ct_meta <- meta[meta$celltype == ct, ]
  ct_counts <- count_data[, ct_meta$sample]
  
  dds_ct <- DESeqDataSetFromMatrix(countData = ct_counts, colData = ct_meta, design = ~ region + age)
  dds_ct <- DESeq(dds_ct)
  res_ct <- results(dds_ct, contrast = c("age", "18mo", "2mo"))
  
  write.csv(as.data.frame(res_ct), file = paste0("TE_deseq_celltype_", gsub(" ", "_", ct), "_18vs2.csv"))
}
