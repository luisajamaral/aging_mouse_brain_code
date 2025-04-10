library(ggplot2)
library(pheatmap)
library(edgeR)


data = read.table("male_H3K9me3_reprocessed_P0_forebrain_counts.txt", header =T)

colnames(data) = gsub("...male_split_bams.bams.", "", colnames(data))
colnames(data) = gsub(".Male", "", colnames(data))
colnames(data) = gsub(".bam", "", colnames(data))
data = data[,-grep("merged",colnames(data))]
head(data)

rownames(data) = paste(data$Chr, ":",data$Start, "-",data$End, sep="")

all_counts = data

mapping_key <- data.frame(
  Codes = c("2A+3A", "8E+9H+8J+9J", "11E+11F+12E", "12D+13B", "3F+4E", 
            "5E+6E", "7H+8H+9G", "13D+14C"),
  Region = c("FC", "HCA", "HCP", "ENT", "NAC", "CP", "AMY", "RLP"),
  stringsAsFactors = FALSE
)


# Split Codes into individual entries for easy matching
mapping_key <- mapping_key %>%
  mutate(Codes = strsplit(Codes, "\\+")) %>%
  tidyr::unnest(Codes)

# Assign regions to columns
assign_region <- function(colname, mapping_key) {
  # Check if the column name contains any code from the mapping key
  match <- mapping_key$Codes[sapply(mapping_key$Codes, function(code) grepl(code, colname))]
  if (length(match) > 0) {
    # Return the corresponding region
    return(mapping_key$Region[which(mapping_key$Codes == match[1])])
  } else {
    return(NA)  # If no match found, return NA
  }
}

regions <- sapply(colnames(all_counts), assign_region, mapping_key = mapping_key)


# Function to extract cell type
extract_celltype <- function(colname) {
  # Use regex to capture everything before ".2mo", ".9mo", or ".18mo"
  celltype <- sub("\\.(2mo|9mo|18mo).*", "", colname)
  return(celltype)
}

# Load required libraries

for(reg in c("HCA", "HCP","AMY","ENT",  "NAC", "RLP", "CP","FC")){
dir.create(paste("../male_",reg, sep = ""), showWarnings = FALSE)
setwd(paste("../male_",reg, sep = ""))

# Load required libraries
library(edgeR)
library(ggplot2)

counts = (all_counts[,grep(reg, regions)])

# Define ages explicitly for matching
ages <- c("2mo", "9mo", "18mo")

# Function to parse column names
parse_colname <- function(colname) {
  # Extract cell type (everything before ".FC")
  celltype <-  extract_celltype(colname)
  
  # Extract age using predefined list
  age <- ages[sapply(ages, function(x) grepl(x, colname))]
  
  # Extract replicate (last part after splitting by "_")
  rep <- tail(strsplit(colname, "_")[[1]], 1)
  
  list(CellType = celltype, Age = age, Replicate = rep)
}

# Parse all column names
metadata <- do.call(rbind, lapply(colnames(counts), parse_colname))
metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
rownames(metadata) <- colnames(counts)

# Make age syntactically valid in metadata
metadata$Age <- make.names(metadata$Age)  # Converts "2mo" to "X2mo", "18mo" to "X18mo", etc.

metadata$Replicate = gsub("rep", "", metadata$Replicate)
                     
# Print metadata for validation
print(head(metadata))
print(tail(metadata))

# Get unique cell types
unique_cell_types <- unique(metadata$CellType)

# Loop through each cell type
for (cell_type in unique_cell_types) {
    # Subset counts for the current cell type
    cell_indices <- which(metadata$CellType == cell_type)
    if(length(cell_indices) < 6) {
        next
    }
    cell_counts <- counts[, cell_indices]
    
    # Extract group (age) information for this cell type
    group <- factor(metadata$Age[cell_indices], levels = unique(metadata$Age))
    if(any(table(group)<2)){
        next
    }
    # Create DGEList
    dge <- DGEList(counts = cell_counts, group = group)
    
    # Filter low counts
    keep <- filterByExpr(dge)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    # Normalize and estimate dispersions
    dge <- calcNormFactors(dge)
    
    # Create design matrix
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    
    # Estimate dispersion
    dge <- estimateDisp(dge, design)
    
    # Fit a GLM model
    fit <- glmQLFit(dge, design)
    
    # Pairwise comparisons
    comparisons <- list(
        `2mo_vs_9mo` = c("X9mo", "X2mo"),
        `9mo_vs_18mo` = c("X18mo", "X9mo"),
        `2mo_vs_18mo` = c("X18mo", "X2mo")
    )
    
    for (comparison in names(comparisons)) {
        groups <- comparisons[[comparison]]
        if (all(groups %in% levels(group))) {
            contrast <- makeContrasts(
                contrasts = paste(groups[1], "-", groups[2]), 
                levels = design
            )
            qlf <- glmQLFTest(fit, contrast = contrast)
            
            # Save results
            output_file <- paste0(cell_type, "_", comparison, "_differential_peaks.txt")
            results <- topTags(qlf, n = Inf)$table
            write.table(results, output_file, sep = "\t", quote = FALSE)
            
            # Generate volcano plot
            results$Significance <- "Neutral"
            results$Significance[results$logFC > 1 & results$PValue < 0.05] <- "Upregulated"
            results$Significance[results$logFC < -1 & results$PValue < 0.05] <- "Downregulated"
            
            max_abs_fc <- max(abs(results$logFC), na.rm = TRUE)
            x_range <- c(-max_abs_fc, max_abs_fc)
            
            volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = Significance)) +
                geom_point(alpha = 0.7) +
                scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Neutral" = "gray")) +
                xlim(x_range) +
                theme_minimal() +
                labs(
                    title = paste(cell_type, comparison),
                    x = "log2 Fold Change",
                    y = "-log10(P-value)"
                ) +
                geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
                geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")
            
            ggsave(
                filename = paste0(cell_type, "_", comparison, "_volcano.pdf"),
                plot = volcano_plot,
                height = 5, width = 5
            )
        }
    }
}



#plot merged plots:




# Define comparisons and their file patterns
comparisons <- c("2mo_vs_18mo", "9mo_vs_18mo", "2mo_vs_9mo")
comparison_files <- list(
    `2mo_vs_18mo` = list.files(pattern = "_2mo_vs_18mo_differential_peaks.txt$"),
    `9mo_vs_18mo` = list.files(pattern = "_9mo_vs_18mo_differential_peaks.txt$"),
    `2mo_vs_9mo` = list.files(pattern = "_2mo_vs_9mo_differential_peaks.txt$")
)

# Function to parse row names and extract domain size
parse_domain_size <- function(rowname) {
  # Split the row name by ":" and "-" to extract positions
  parts <- unlist(strsplit(rowname, "[:-]"))
  if (length(parts) >= 3) {
    # Ensure numeric conversion of positions
    start_pos <- as.numeric(parts[2])
    end_pos <- as.numeric(parts[3])
    if (!is.na(start_pos) && !is.na(end_pos)) {
      return(abs(end_pos - start_pos))  # Domain size is the length of the genomic region
    }
  }
  return(NA)  # Return NA if parsing fails
}

# Function to add clade information based on cell type
assign_clade <- function(cell_type) {
  if (grepl("Glut", cell_type)) {
    "Glut"
  } else if (grepl("Gaba", cell_type)) {
    "Gaba"
  } else {
    "NN"
  }
}

# Loop through comparisons and generate volcano plots
for (comparison in comparisons) {
    # Load and combine all files for the comparison
    files <- comparison_files[[comparison]]
    combined_data <- data.frame()
    
    for (file in files) {
        # Extract cell type from filename
        cell_type <- sub(paste0("_", comparison, "_differential_peaks.txt"), "", file)
        
        # Load the data
        data <- read.table(file, header = TRUE)
        data$CellType <- cell_type
        data$Clade <- assign_clade(cell_type)
        data$DomainSize <- sapply(rownames(data), parse_domain_size)
        combined_data <- rbind(combined_data, data)
    }
    
    # Filter out rows with NA values
    combined_data <- combined_data %>%
        filter(!is.na(logFC) & !is.na(PValue) & !is.na(DomainSize))
    
    # Add significance information
    combined_data$Significance <- "Neutral"
    combined_data$Significance[combined_data$logFC > 1 & combined_data$PValue < 0.05] <- "Upregulated"
    combined_data$Significance[combined_data$logFC < -1 & combined_data$PValue < 0.05] <- "Downregulated"
    
    # Determine x-axis range based on maximum absolute logFC
    max_abs_fc <- max(abs(combined_data$logFC), na.rm = TRUE)
    x_range <- c(-max_abs_fc, max_abs_fc)
    
    # Create a volcano plot
    volcano_plot <- ggplot(combined_data, aes(x = logFC, y = -log10(FDR), color = Clade, size = DomainSize)) +
        geom_point(alpha = 0.7) +

        scale_color_manual(values = c("Glut" = "firebrick1", "Gaba" = "green2", "NN" = "blue3")) +
        xlim(x_range) +
        theme_minimal() +
        labs(
            title = paste("Volcano Plot for", comparison),
            x = "log2 Fold Change",
            y = "-log10(FDR)"
        ) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black")
    
    # Save the plot
    ggsave(filename = paste0("volcano_plot_", comparison, ".pdf"), plot = volcano_plot, height = 5, width = 6.5)
}
                     }

                     





                     

groups


