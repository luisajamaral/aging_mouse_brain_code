# DESeq2 Analysis of Transposable Elements (TEs) by Cell Type using edgeR
# Description: Differential analysis of TE subfamilies across age groups (2mo vs 18mo) using counts from SoloTE.

library(edgeR)
library(Matrix)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)

# === Load Data ===
setwd("~/projects/combined_all/female_RNA/SoloTE")
seuratobj <- readRDS("all_subfamily.RDS")  # Replace with subfamily RDS if needed

# Extract counts and metadata
counts <- GetAssayData(seuratobj, slot = "counts")
meta <- seuratobj@meta.data

# === Differential TE Analysis by Cell Type ===
for (ct in unique(meta$celltype)) {
  cat(ct, "\n")
  cl <- gsub("[ /]", "-", ct)

  # Subset counts for cell type, remove 9mo samples
  ct_cols <- grep(ct, colnames(counts))
  counts_ct <- counts[, ct_cols]
  counts_ct <- counts_ct[, -grep("9mo", colnames(counts_ct))]

  # Extract group (age) and region from column names
  groups <- sapply(strsplit(colnames(counts_ct), "_"), function(x) strsplit(x[2], "-")[[1]][2])
  regs <- sapply(strsplit(colnames(counts_ct), "_"), function(x) strsplit(x[2], "-")[[1]][1])

  # Filter for samples with enough reads and region replicates
  min_reads <- 10000
  valid_regions <- names(which(table(regs[colSums(counts_ct) > min_reads]) > 3))
  if (length(valid_regions) < 1) next
  keep_cols <- which(regs %in% valid_regions)
  counts_ct <- counts_ct[, keep_cols]

  # Redefine groups and regions after filtering
  groups <- factor(sapply(strsplit(colnames(counts_ct), "_"), function(x) strsplit(x[2], "-")[[1]][2]),
                   levels = c("18mo", "02mo"))
  regs <- sapply(strsplit(colnames(counts_ct), "_"), function(x) strsplit(x[2], "-")[[1]][1])

  # Create edgeR object
  y <- DGEList(counts = counts_ct, group = groups)
  keep <- rowSums(counts_ct >= 10) >= 2
  y <- y[keep,, keep.lib.sizes = FALSE]
  y <- y[grep("Solo", rownames(y$counts)),, keep.lib.sizes = TRUE]

  y <- calcNormFactors(y)
  design <- if (length(unique(regs)) > 1) model.matrix(~groups + regs) else model.matrix(~groups)

  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = "groups02mo")

  # Save results
  fdr <- p.adjust(lrt$table$PValue, method = "BH")
  out <- cbind(lrt$table, fdr)
  out <- out[order(out$fdr), ]
  out$significance <- ifelse(out$fdr < 0.05, "significant", "not significant")
  out$dir <- "ns"
  out$dir[out$logFC > 0 & out$fdr < 0.05] <- "down"
  out$dir[out$logFC < 0 & out$fdr < 0.05] <- "up"

  write.table(out, file = paste0(cl, "_DE_TEs_only.txt"), sep = "\t")

  # Volcano plot
  out_TE <- out[grep("Solo", rownames(out)), ]
  top10 <- head(out_TE[order(out_TE$PValue), ], 10)
  top10$TE <- rownames(top10)

  p <- ggplot(out_TE, aes(x = -logFC, y = -log10(fdr), color = dir)) +
    geom_point() +
    scale_color_manual(values = c("up" = "red", "down" = "blue")) +
    theme_minimal() +
    labs(title = paste(cl, "TEs"), x = "log2 Fold Change", y = "-log10(FDR)") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme(legend.position = "top") +
    geom_text_repel(data = top10, aes(label = TE), vjust = 1.5, hjust = 1.5)

  ggsave(p, file = paste0(cl, "_volcano_TEs_only.pdf"))
}
