# Load libraries
library(Matrix)
library(data.table)
library(future)
library(doFuture)
library(foreach)
library(lme4)
library(ggplot2)
library(Seurat)
library(dplyr)
library(lme4)
library(Matrix)
if (!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel")
}
library(doParallel)
library(foreach)
library(parallel)

obj = readRDS("../female_RNA/RNA_final_SCT.RDS")

DefaultAssay(obj) = "RNA"

tracksheet= read.csv("tracksheet_with_paths_oct23.csv")
rownames(tracksheet) = tracksheet$Sample_name
obj$track_Sample= tracksheet[gsub("0", "",paste(obj$orig.ident)),"Sample"]
obj$track_Sample= tracksheet[gsub("0", "",paste(obj$orig.ident)),"Sample"]
obj$is_Frozen = "no"
obj$is_Frozen[grep("Frozen", obj$track_Sample)] = "yes"
DimPlot(obj,reduction="umap",group.by = "is_Frozen",label = T)

head(obj@meta.data[,c('orig.ident','nCount_RNA','nFeature_RNA','percent.mt','percent.ribo',
                     'class','celltype_final','barcode','region','age','rep','is_Frozen')])

colnames(obj@meta.data)

# Load required libraries
library(doParallel)
library(foreach)
library(parallel)
library(dplyr)
library(lme4)

# Functions
subset_data <- function(obj, region, celltype) {
  cat(paste0("Subsetting data for region: ", region, ", cell type: ", celltype, "\n"))
  
  # Subset metadata for the specific region and cell type
  metadata <- obj@meta.data %>%
    filter(region == region & celltype_final == celltype)
  
  # Subset expression matrix to match the metadata
  expression_matrix <- obj@assays$RNA@data[, rownames(metadata)]
  
  return(list(metadata = metadata, expression_matrix = expression_matrix))
}

filter_genes <- function(expr_matrix, min_cells = 10, min_fraction = 0.1) {
  cat("Filtering genes...\n")
  
  keep <- rowSums(expr_matrix > 0) >= min_cells & 
          rowMeans(expr_matrix > 0) >= min_fraction
  return(expr_matrix[keep, ])
}

fit_gene_model <- function(expr_matrix, metadata, fixed_effect, random_effects, latent_vars = NULL) {
  cat("Fitting mixed-effects models...\n")
  results <- foreach(gene = rownames(expr_matrix), .combine = rbind, .errorhandling = "pass") %dopar% {
    tryCatch({
      # Prepare the data
      data <- metadata
      data$gene_expression <- expr_matrix[gene, ]
      
      # Construct formula dynamically
      formula <- paste("gene_expression ~", fixed_effect, 
                       if (!is.null(latent_vars)) paste(latent_vars, collapse = " + "), 
                       paste("(1 |", random_effects, ")", collapse = " + "))
      
      # Fit model
      model <- lmer(as.formula(formula), data = data, REML = FALSE)
      
      # Extract p-value
      coef_summary <- summary(model)
      p_value <- coef_summary$coefficients[fixed_effect, "Pr(>|t|)"]
      
      # Return results
      data.frame(gene = gene, p_value = p_value)
    }, error = function(e) {
      # Return NA in case of error
      data.frame(gene = gene, p_value = NA)
    })
  }
  
  return(results)
}

# Setup parallel processing
num_cores <- 14
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Main Analysis
regions <- unique(obj@meta.data$region)
cell_types <- unique(obj@meta.data$celltype_final)

# Output directory
output_dir <- "differential_expression_results"
if (!dir.exists(output_dir)) dir.create(output_dir)

for (region in regions) {
  for (cell_type in cell_types) {
    cat(paste0("Processing region: ", region, ", cell type: ", cell_type, "\n"))
    
    # Subset data
    data <- subset_data(obj, region, cell_type)
    metadata <- data$metadata
    expr_matrix <- filter_genes(data$expression_matrix, min_cells = 10, min_fraction = 0.1)
    
    if (nrow(expr_matrix) == 0) {
      cat(paste0("No genes passed filtering for region: ", region, ", cell type: ", cell_type, "\n"))
      next
    }
    
    # Fit mixed-effects models
    results <- fit_gene_model(
      expr_matrix = expr_matrix,
      metadata = metadata,
      fixed_effect = "age",
      random_effects = "rep",  # Assuming `rep` is your sample replicate
      latent_vars = c("percent.mt", "percent.ribo")
    )
    
    # Save results to file
    output_file <- file.path(output_dir, paste0(region, "_", cell_type, "_DEGs.csv"))
    write.csv(results, output_file, row.names = FALSE)
    
    cat(paste0("Saved results for region: ", region, ", cell type: ", cell_type, 
               " -> ", nrow(results), " genes processed.\n"))
  }
}

# Stop parallel processing
stopCluster(cl)

cat("Analysis completed. Results saved in: ", output_dir, "\n")


library(doParallel)
library(Matrix)
library(lme4)
library(foreach)

# Initialize parallel backend
num_cores <- 14
cl <- makeCluster(num_cores, outfile = "worker_log.txt")  # Log to file
registerDoParallel(cl)

# Minimize metadata
metadata <- metadata[, c("age", "rep", "percent.mt", "percent.ribo", "is_Frozen")]

# Convert expression matrix to sparse format
expr_matrix <- as(expr_matrix, "dgCMatrix")

# Define parameters
fixed_effect <- "age"
random_effects <- "rep"
latent_vars <- c("percent.mt", "percent.ribo", "is_Frozen")

# Fitting models with parallel processing
cat("Fitting mixed-effects models...\n")
results <- foreach(
  gene = rownames(expr_matrix),
  .combine = rbind,
  .packages = c("lme4"),
  .errorhandling = "pass"
) %dopar% {
  tryCatch({
    # Prepare gene-specific data
    data <- metadata
    data$gene_expression <- expr_matrix[gene, ]
    
    # Construct formula dynamically
    formula <- as.formula(paste(
      "gene_expression ~", fixed_effect, 
      if (!is.null(latent_vars)) paste("+", paste(latent_vars, collapse = " + ")), 
      "+ (1 |", random_effects, ")"
    ))
    
    # Fit model
    model <- lmer(formula, data = data, REML = FALSE)
    
    # Extract fixed effect p-value
    coef_summary <- summary(model)
    p_value <- coef_summary$coefficients[fixed_effect, "Pr(>|t|)"]
    
    return(data.frame(gene = gene, p_value = p_value))
  }, error = function(e) {
    # Handle errors gracefully
    return(data.frame(gene = gene, p_value = NA))
  })
}

# Stop cluster
stopCluster(cl)

# Review results
cat("Model fitting complete.\n")
saveRDS(results, file = "differential_expression_results.rds")
head(results)


region = "HCP"
cell_type = "Oligo NN"

# Subset data
data <- subset_data(obj, region, cell_type)
data
metadata <- data$metadata
head(metadata)
expr_matrix <- filter_genes(data$expression_matrix, min_cells = 10, min_fraction = 0.1)

if (nrow(expr_matrix) == 0) {
  cat(paste0("No genes passed filtering for region: ", region, ", cell type: ", cell_type, "\n"))
  next
}

# Fit mixed-effects models
results <- fit_gene_model(
  expr_matrix = expr_matrix,
  metadata = metadata,
  fixed_effect = "age",
  random_effects = "rep",  # Assuming `rep` is your sample replicate
  latent_vars = c("percent.mt", "percent.ribo", "is_Frozen")
)

region = "HCP"
cell_type = "Oligo NN"

# Subset data
data <- subset_data(obj, region, cell_type)
data
metadata <- data$metadata
head(metadata)
expr_matrix <- filter_genes(data$expression_matrix, min_cells = 10, min_fraction = 0.1)
 
expr_matrix = expr_matrix
  metadata = metadata
  fixed_effect = "age"
  random_effects = "rep"  # Assuming `rep` is your sample replicate
  latent_vars = c("percent.mt", "percent.ribo", "is_Frozen")
cat("Fitting mixed-effects models...\n")
  results <- foreach(gene = rownames(expr_matrix), .combine = rbind, .errorhandling = "pass") %dopar% {
    tryCatch({
      # Prepare the data
      data <- metadata
      data$gene_expression <- expr_matrix[gene, ]
      
      # Construct formula dynamically
      formula <- paste("gene_expression ~", fixed_effect, 
                       if (!is.null(latent_vars)) paste(latent_vars, collapse = " + "), 
                       paste("(1 |", random_effects, ")", collapse = " + "))
      
      # Fit model
      model <- lmer(as.formula(formula), data = data, REML = FALSE)
      
      # Extract p-value
      coef_summary <- summary(model)
      p_value <- coef_summary$coefficients[fixed_effect, "Pr(>|t|)"]
      
      # Return results
      head(data.frame(gene = gene, p_value = p_value))
        
    }, error = function(e) {
      # Return NA in case of error
      data.frame(gene = gene, p_value = NA)
    })
  }

head(significant_results)

subset_data <- function(obj, region, celltype) {
  # Subset metadata for the specific region and cell type
  metadata <- obj@meta.data %>%
    filter(region == region & celltype_final == celltype)
  
  # Subset expression matrix to match the metadata
  expression_matrix <- obj@assays$RNA@data[, rownames(metadata)]
  
  return(list(metadata = metadata, expression_matrix = expression_matrix))
}


filter_genes <- function(expr_matrix, min_cells = 10, min_fraction = 0.1) {
  keep <- rowSums(expr_matrix > 0) >= min_cells & 
          rowMeans(expr_matrix > 0) >= min_fraction
  return(expr_matrix[keep, ])
}

fit_gene_model <- function(expr_matrix, metadata, fixed_effect, random_effects, latent_vars = NULL) {
  results <- list()
  
  for (gene in rownames(expr_matrix)) {
    tryCatch({
      # Prepare the data
      data <- metadata
      data$gene_expression <- expr_matrix[gene, ]
      
      # Construct formula dynamically
      formula <- paste("gene_expression ~", fixed_effect, 
                       paste(latent_vars, collapse = " + "), 
                       paste("(1 |", random_effects, ")", collapse = " + "))
      
      # Fit model
      model <- lmer(as.formula(formula), data = data, REML = FALSE)
      
      # Extract p-value
      coef_summary <- summary(model)
      p_value <- coef_summary$coefficients[fixed_effect, "Pr(>|t|)"]
      
      # Store results
      results[[gene]] <- list(coef = coef_summary$coefficients, p_value = p_value)
    }, error = function(e) {
      # Handle errors gracefully
      results[[gene]] <- list(error = as.character(e))
    })
  }
  
  return(results)
}




# Define regions and cell types
regions <- unique(obj@meta.data$region)
cell_types <- unique(obj@meta.data$celltype_final)

# Iterate through regions and cell types
for (region in regions) {
  for (cell_type in cell_types) {
    # Subset data
    data <- subset_data(obj, region, cell_type)
    metadata <- data$metadata
    expr_matrix <- filter_genes(data$expression_matrix, min_cells = 10, min_fraction = 0.1)
    
    if (nrow(expr_matrix) == 0) {
      cat(paste("No genes passed filtering for", cell_type, "in", region, "\n"))
      next
    }
    
    # Fit mixed-effects models
    results <- fit_gene_model(
      expr_matrix = expr_matrix,
      metadata = metadata,
      fixed_effect = "age",
      random_effects = "rep",  # Assuming `rep` is your sample replicate
      latent_vars = c("percent.mt", "percent.ribo")
    )
    
    # Save or process results
    significant_genes <- do.call(rbind, lapply(results, function(res) {
      if (!is.null(res$p_value) && res$p_value < 0.05) {
        return(data.frame(gene = names(res), p_value = res$p_value))
      }
    }))
    
    cat(paste("Processed:", cell_type, "in", region, "-> Significant genes:", nrow(significant_genes), "\n"))
  }
}


region = "HCP"
cell_type = "Oligo NN"
data <- subset_data(obj, region, cell_type)


metadata <- data$metadata
expr_matrix <- filter_genes(data$expression_matrix, min_cells = 10, min_fraction = 0.1)
    

# Convert expression matrix to a sparse matrix
expr_matrix <- as(expr_matrix, "dgCMatrix")


# Ensure metadata contains necessary columns
metadata$age <- factor(metadata$age)  # Convert age to a factor
metadata$rep <- factor(metadata$rep)  # Convert replicates to a factor


results <- fit_gene_model(
      expr_matrix = expr_matrix,
      metadata = metadata,
      fixed_effect = "age",
      random_effects = "rep",  # Assuming `rep` is your sample replicate
      latent_vars = c("percent.mt", "percent.ribo", "is_Frozen")
    )

 significant_genes <- do.call(rbind, lapply(results, function(res) {
      if (!is.null(res$p_value) && res$p_value < 0.05) {
        return(data.frame(gene = names(res), p_value = res$p_value))
      }
    }))
    

