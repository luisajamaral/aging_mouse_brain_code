# GO and GSEA Enrichment Analysis using ClusterProfiler

# Required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(data.table)

# Set input and output directories
gene_list_dir <- "./dars_gene_lists/"
ranked_list_dir <- "./ranked_gene_lists/"  # for GSEA, optional
output_dir <- "./go_gsea_results/"
dir.create(output_dir, showWarnings = FALSE)

# Load expressed genes background per cell type (gene universe)
universe_file <- "./background_genes_by_celltype.csv"  # CSV with celltype, gene columns
universe_df <- fread(universe_file)

# Helper function for GO enrichment
run_go_enrichment <- function(gene_file, celltype, direction) {
  genes <- scan(gene_file, what = "character")
  universe <- universe_df[universe_df$celltype == celltype, gene]

  ego <- enrichGO(
    gene = genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    universe = universe
  )

  out_file <- file.path(output_dir, paste0(celltype, "_", direction, "_GO.csv"))
  write.csv(as.data.frame(ego), out_file, row.names = FALSE)
}

# Helper function for GSEA
run_gsea <- function(ranked_file, celltype, direction) {
  df <- fread(ranked_file)
  ranked_genes <- df$logFC
  names(ranked_genes) <- df$gene
  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  universe <- universe_df[universe_df$celltype == celltype, gene]
  ranked_genes <- ranked_genes[names(ranked_genes) %in% universe]

  gsea_result <- gseGO(
    geneList = ranked_genes,
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    verbose = FALSE
  )

  out_file <- file.path(output_dir, paste0(celltype, "_", direction, "_GSEA.csv"))
  write.csv(as.data.frame(gsea_result), out_file, row.names = FALSE)
}

# Loop through gene lists and run both GO and GSEA if ranked list is available
files <- list.files(gene_list_dir, pattern = "*.txt", full.names = TRUE)
for (gene_file in files) {
  fname <- basename(gene_file)
  parts <- strsplit(fname, "["]|_|[.]", perl = TRUE)[[1]]
  celltype <- parts[1]
  direction <- parts[2]

  cat("Running GO for:", fname, "\n")
  run_go_enrichment(gene_file, celltype, direction)

  ranked_file <- file.path(ranked_list_dir, paste0(celltype, "_", direction, "_ranked.csv"))
  if (file.exists(ranked_file)) {
    cat("Running GSEA for:", ranked_file, "\n")
    run_gsea(ranked_file, celltype, direction)
  }
}
