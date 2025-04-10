library(pheatmap)
library(data.table)
library(matrixStats)
library(IRdisplay)
library(ggplot2)
library(RColorBrewer)

setwd("~/projects/combined_all/female_RNA/DEG_results_latent_isfrozen_sample_removed_lowq_samples/")

setwd("EnrichR")

files_down

# Load necessary libraries
library(ggplot2)
library(RColorBrewer)
library(scales)
library(dplyr)
library(data.table)
library(tidyr)

# Parameters
database <- "_GO_BP"  # Specify the database
pval <- 0.01  # Adjusted p-value threshold for significance
min_genes_in_pathway <- 5  # Minimum genes required in the pathway for inclusion

# Step 1: Define function to process UP and DOWN files
process_files <- function(files, direction) {
  data_list <- list()  # List to collect data frames for each file
  
  for (f in files) {
    # Read each file and convert to data frame
    curr <- fread(f, sep = ",", fill = TRUE)
    curr <- as.data.frame(curr)
    
    # Filter by significance threshold
    curr <- curr[curr$`Adjusted P-value` < pval, ]
    
    # Extract Celltype name from filename
    Celltype <- gsub(paste0(database, "|UP|DOWN|\\.csv|__"), "", f)
    
    # If there are significant pathways, append Celltype and Direction labels
    if (nrow(curr) > 0) {
      curr$Celltype <- Celltype
      curr$Direction <- direction
      data_list[[length(data_list) + 1]] <- curr
    }
  }
  
  # Combine all data frames into a single data frame
  return(do.call(rbind, data_list))
}

# Step 2: Load and label UP and DOWN files
files_up <- list.files(pattern = paste0("UP.*", database, ".*\\.csv"), path = ".")
files_down <- list.files(pattern = paste0("DOWN.*", database, ".*\\.csv"), path = ".")

# Process UP and DOWN files and combine them into a single dataset
up_data <- process_files(files_up, "UP")
down_data <- process_files(files_down, "DOWN")
combined_data <- rbind(up_data, down_data)

# Step 3: Annotate `clade` and `region` based on `Celltype`
combined_data$clade <- ""
combined_data$clade[grep("Glut", combined_data$Celltype)] <- "Glut"
combined_data$clade[grep("Dopa", combined_data$Celltype)] <- "Glut"
combined_data$clade[grep("Gaba", combined_data$Celltype)] <- "Gaba"
combined_data$clade[grep("Neur", combined_data$Celltype)] <- "Gaba"
combined_data$clade[grep("NN", combined_data$Celltype)] <- "NN"
combined_data$clade[grep("IMN", combined_data$Celltype)] <- "NN"

combined_data$clade <- factor(combined_data$clade, levels = c("Glut", "Gaba", "NN"))
combined_data$region <- sapply(strsplit(as.character(combined_data$Celltype), "--"), `[`, 2)

# Step 4: Calculate gene count and create Color_metric for plotting
combined_data$Gene_count <- sapply(strsplit(as.character(combined_data$Genes), ";"), length)
combined_data$Direction_numeric <- ifelse(combined_data$Direction == "UP", 1, -1)
combined_data$Color_metric <- -log10(combined_data$`Adjusted P-value`) * combined_data$Direction_numeric

# Cap Color_metric within the range [-5, 5] for consistent coloring
combined_data$Color_metric <- pmin(pmax(combined_data$Color_metric, -5), 5)

# Function to filter out terms with high gene overlap
filter_terms_by_gene_overlap <- function(data, overlap_threshold = 0.8) {
  # Sort terms by significance (Adjusted P-value), ascending
  sorted_terms <- data %>%
    arrange(`Adjusted P-value`) %>%
    distinct(Term, Genes, `Adjusted P-value`)
  
  # List to keep track of terms to keep
  terms_to_keep <- list()
  
  # Loop over each term in order of significance
  for (i in seq_len(nrow(sorted_terms))) {
    current_term <- sorted_terms$Term[i]
    current_genes <- strsplit(sorted_terms$Genes[i], ";")[[1]]
    
    # Check overlap with previously retained terms
    retain <- TRUE
    for (kept in terms_to_keep) {
      kept_genes <- strsplit(kept$Genes, ";")[[1]]
      overlap <- length(intersect(current_genes, kept_genes)) / length(current_genes)
      
      if (overlap > overlap_threshold) {
        retain <- FALSE
        break
      }
    }
    
    # If not redundant, add the term to terms_to_keep
    if (retain) {
      terms_to_keep[[length(terms_to_keep) + 1]] <- sorted_terms[i, ]
    }
  }
  
  # Return the final list of terms to keep
  return(data %>% filter(Term %in% sapply(terms_to_keep, function(x) x$Term)))
}

# Step 5: Define function to create plots by clade with clustering
plot_by_clade <- function(clade_name) {
    # Filter data for the specified clade
    clade_data <- combined_data %>% filter(clade == clade_name)
    
    # Determine significant terms for the current clade
    significant_terms <- clade_data %>%
      filter(`Adjusted P-value` < pval & Gene_count >= min_genes_in_pathway) %>%
      pull(Term) %>%
      unique()
    
    # Filter the data to include only significant terms for the current clade,
    # but retain all rows for these terms across cell types in the clade
    clade_filtered <- clade_data %>%
      filter(Term %in% significant_terms) %>%
      mutate(Significant = ifelse(`Adjusted P-value` < pval, TRUE, FALSE))
    
    # Apply gene overlap filtering to remove redundant terms
    clade_filtered <- filter_terms_by_gene_overlap(clade_filtered)
    
    # Ensure uniqueness by taking the most significant Color_metric for each Term + Celltype + Region
    clade_filtered <- clade_filtered %>%
      group_by(Term, Celltype, region) %>%
      slice_max(order_by = abs(Color_metric), n = 1) %>%
      ungroup()
    
    # Prepare data for clustering
    clustering_data <- clade_filtered %>%
      select(Term, Celltype, Color_metric) %>%
      spread(key = Celltype, value = Color_metric, fill = 0)
    
    # Hierarchical clustering on terms
    dist_matrix <- dist(clustering_data %>% select(-Term))
    term_order <- hclust(dist_matrix)$order
    ordered_terms <- clustering_data$Term[term_order]
    clade_filtered$Term <- factor(clade_filtered$Term, levels = ordered_terms)
    
    # Define color gradient for `Color_metric`
    color_palette <- scale_color_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        limits = c(-5, 5), oob = squish, name = "Significance * Direction"
    )
    
    # Create the dot plot with faceting by region in one row
    p <- ggplot(clade_filtered, aes(x = Celltype, y = Term, size = `Odds Ratio`, color = Color_metric)) +
        # Plot nonsignificant points with lower alpha and smaller size for visibility
        geom_point(data = subset(clade_filtered, Significant == FALSE), 
                   alpha = 0.4, size = 2, shape = 1) +
        # Plot significant points with full color and size
        geom_point(data = subset(clade_filtered, Significant == TRUE), 
                   alpha = 0.7) +
        color_palette +
        scale_size_continuous(name = "Odds Ratio", range = c(1, 6), breaks = c(1, 3, 5, 10)) +
        theme_minimal() +
        facet_grid(~ region,  scales = "free_x", space = "free") +
        labs(title = paste("Pathway Enrichment for", clade_name, "by Region -", database), 
             x = "Cell Type", 
             y = "Pathway") +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            strip.text = element_text(size = 10),
            plot.title = element_text(size = 14, face = "bold"),
            legend.position = "right"
        )
    
    # Display the plot
    print(p)
    
    # Save each plot as a PDF with the clade name
   # ggsave(paste0("Enrichment_Dotplot_", clade_name, "_", database, ".pdf"), plot = p, height = 8, width = 16)
}

# Generate plots for each clade
plot_by_clade("Glut")
plot_by_clade("Gaba")
plot_by_clade("NN")


options(repr.plot.width = 16, repr.plot.height = 6)

plot_by_clade("Glut")
plot_by_clade("Gaba")
plot_by_clade("NN")

head(df)

library(data.table)
library(ggplot2)
library(RColorBrewer)

options(repr.plot.width = 20, repr.plot.height = 20)

database <- "_Reactome_2022"
pval <- 0.05

# Load and label UP and DOWN files
files_up <- list.files(pattern = paste0("UP.*", database, ".*\\.csv"), path = ".")
files_down <- list.files(pattern = paste0("DOWN.*", database, ".*\\.csv"), path = ".")

# Function to process files and add direction label
process_files <- function(files, direction) {
  gkey <- c()
  data_list <- list()
  
  for (f in files) {
    curr <- fread(f, sep = ",", fill = TRUE)
    curr <- as.data.frame(curr)
    curr <- curr[which(curr$`Adjusted P-value` < pval), ]
    
    # Extract Celltype name from filename
    Celltype <- gsub(paste0(database, "|UP|DOWN|\\.csv|__"), "", f)
    
    # Append data with directional label
    if (nrow(curr) > 0) {
      curr$Celltype <- Celltype
      curr$Direction <- direction
      data_list[[length(data_list) + 1]] <- curr
      gkey <- c(gkey, unique(curr$Term))
    }
  }
  
  gkey <- unique(gkey)
  return(list(data = do.call(rbind, data_list), gkey = gkey))
}

# Process UP and DOWN files
up_data <- process_files(files_up, "UP")
down_data <- process_files(files_down, "DOWN")

# Combine data, keeping unique terms across both directions
combined_data <- rbind(up_data$data, down_data$data)
gkey <- unique(c(up_data$gkey, down_data$gkey))

# Keep only the top pathway per Celltype for both UP and DOWN directions
combined_data <- combined_data[combined_data$Term %in% gkey, ]
combined_data <- combined_data[order(-combined_data$`Adjusted P-value`), ]

# Generate plot data
df <- expand.grid(Celltype = unique(combined_data$Celltype), Term = gkey)
df <- merge(df, combined_data, by = c("Celltype", "Term"), all.x = TRUE)

# Set colors for UP (red) and DOWN (blue) with point layering logic
colors_up <- "red"
colors_down <- "blue"

# Plotting UP and DOWN data with layering
p <- ggplot(df, aes(x = Celltype, y = Term)) +
  geom_point(data = subset(df, Direction == "DOWN"), aes(size = -log10(`Adjusted P-value`), color = colors_down), alpha = 0.7) +
  geom_point(data = subset(df, Direction == "UP"), aes(size = -log10(`Adjusted P-value`), color = colors_up), alpha = 0.7) +
  scale_color_identity() +
  scale_size_continuous(range = c(0, 5), breaks = c(0, 5, 20, 25, 50, 75, 100)) +
  theme_minimal() +
  ggtitle(paste("Combined UP and DOWN Enrichment for", database)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print the plot
print(p)

library(ggplot2)
library(RColorBrewer)
library(data.table)
library(scales)

options(repr.plot.width = 20, repr.plot.height = 8)


# Annotate the clade and region columns based on Celltype
df$clade <- ""
df$clade[grep("Glut", df$Celltype)] <- "Glut"
df$clade[grep("Dopa", df$Celltype)] <- "Glut"
df$clade[grep("Gaba", df$Celltype)] <- "Gaba"
df$clade[grep("Neur", df$Celltype)] <- "Gaba"
df$clade[grep("NN", df$Celltype)] <- "NN"
df$clade[grep("IMN", df$Celltype)] <- "NN"

df$clade <- factor(df$clade, levels = c("Glut", "Gaba", "NN"))

# Extract region information
df$region <- sapply(strsplit(as.character(df$Celltype), "--"), `[`, 2)

# Filter pathways with significant Adjusted P-value, sufficient gene overlap, and presence in at least one cell type
df <- df[df$`Adjusted P-value` < 0.1 & !is.na(df$Genes), ]
df$Gene_count <- sapply(strsplit(as.character(df$Genes), ";"), length)
df <- df[df$Gene_count > 5, ]

# Keep pathways significant in at least one cell type
significant_terms <- unique(df$Term)
df <- df[df$Term %in% significant_terms, ]

# Create a combined metric for coloring: Significance * Direction
df$Direction_numeric <- ifelse(df$Direction == "UP", 1, -1)
df$Color_metric <- -log10(df$`Adjusted P-value`) * df$Direction_numeric  # Positive for UP, Negative for DOWN

# Cap the Color_metric at -5 and 5 for the color scale
df$Color_metric <- pmin(pmax(df$Color_metric, -5), 5)

# Define color gradient with capped range
color_palette <- scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-5, 5), oob = squish
)

# Function to generate plots for each clade and save as PDF
plot_clade <- function(clade_name) {
    # Filter data for the specified clade
    clade_df <- subset(df, clade == clade_name)
    
    # Create the plot with faceting by region
    p <- ggplot(clade_df, aes(x = Celltype, y = Term, size = `Combined Score`, color = Color_metric)) +
        geom_point(alpha = 0.7) +
        color_palette +
        scale_size_continuous(name = "Odds Ratio", range = c(2, 6), breaks = c(0, 5, 10, 20, 50, 100)) +
        theme_minimal() +
        facet_wrap(~ region, scales = "free_x",nrow =1) +
        ggtitle(paste("Enrichment for", clade_name, "by Brain Region")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    # Display the plot
    print(p)
    
    # Save each plot as a PDF
  #  ggsave(paste0(clade_name, "_Enrichment_by_Region.pdf"), plot = p, height = 9, width = 14)
}

# Generate plots for each clade
plot_clade("Glut")
plot_clade("Gaba")
plot_clade("NN")

# Optionally save as PDF
#ggsave(paste0("Combined_", database, "_Dotplot.pdf"), plot = p, height = 9, width = 14)




library(ggplot2)
library(RColorBrewer)
library(scales)
options(repr.plot.width = 20, repr.plot.height = 29)

# Ensure Color_metric is capped within the range [-5, 5]
df$Color_metric <- pmin(pmax(df$Color_metric, -5), 5)

# Define color gradient for capped Color_metric
color_palette <- scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-5, 5), oob = squish, name = "Significance * Direction"
)

# Create the dot plot with faceting by region
p <- ggplot(df, aes(x = Celltype, y = Term, size = `Odds Ratio`, color = Color_metric)) +
    geom_point(alpha = 0.7) +
    color_palette +
    scale_size_continuous(name = "Odds Ratio", range = c(1, 6), breaks = c(1, 3, 5, 10)) +
    theme_minimal() +
    facet_wrap(~ region, scales = "free_x") +
    labs(title = "Pathway Enrichment by Cell Type and Region", 
         x = "Cell Type", 
         y = "Pathway") +
    theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "right"
    )

# Display the plot
print(p)

# Optionally save the plot as a PDF
#ggsave("Enrichment_Dotplot_by_Region.pdf", plot = p, height = 10, width = 14)


df$clade = ""
df$clade[grep("Glut", df$Celltype)] = "Glut"
df$clade[grep("Dopa", df$Celltype)] = "Glut"

df$clade[grep("Gaba", df$Celltype)] = "Gaba"
df$clade[grep("Neur", df$Celltype)] = "Gaba"
df$clade[grep("NN", df$Celltype)] = "NN"
df$clade[grep("IMN", df$Celltype)] = "NN"

df$clade = factor(df$clade , levels = c("Glut", "Gaba", "NN"))


df$region = sapply(strsplit(as.character(df$Celltype), "--"), `[`, 2)


options(repr.plot.width=30, repr.plot.height=8)
colors <- c(brewer.pal(9, "YlOrRd")[2:9])
breaks <- c(0,1,2,3,4,5,10, 25, 50,75)
dir = "Down"

ggplot(df, aes(x = Celltype, y = Term, size = Size, color = Significance)) +
    geom_point() + 
     scale_color_gradientn(colors = colors,breaks=breaks,values = scales::rescale(breaks), 
                          limits = c(0, max(nmat))
                          ) +
    theme_minimal() +
    ggtitle(paste(dir, database)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+labs(size = "# genes overlap", color = "-log10(adj p-value)")+ facet_wrap(~region,scales = "free_x", nrow =1)


head(df)


